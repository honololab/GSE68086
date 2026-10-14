# 06b -- Non-parametric test of independence between batch and the data.
#
# Step 1 (this file so far): who is in each cell of the batch x condition table.
#   - Map each sample to its patient, batch and condition.
#   - Count samples AND patients in each batch x condition cell.
#   - Check whether any patient has samples in more than one batch.
#
# Why patients matter: if one patient gave samples in two batches, those samples are
# not independent, and a test that treats them as separate units would count the same
# person twice. If every patient has one sample, samples and patients are the same
# units and the test can work on samples.
#
# Step 4: design balance check (secondary). Fisher's exact test on the batch x
#   condition table. Not the main test: it is reported because an imbalance changes
#   how the main result is read.
# Steps 5-9: main tests, H0: batch independent of X | condition, by permutation of
#   batch inside each condition. 6 PERMANOVA (centres and spread), 7 HSIC (overall),
#   8 MMD (one pair of batches at a time).
#
# The subset .rds written by 02 keeps only sample, group and batch, so the patient ID
# is taken from the sample sheet written by 01 and joined by sample name.
#
# Base R only (SummarizedExperiment is used only to read the vst values from 04).
# Run from the project root:  Rscript "functions/06b_bivariate_batch_tests.R"

# ================== CHOICES -- edit here ==================
IN     <- "data/derived/subset_lung_gbm_hc.rds"   # written by 02
SHEET  <- "GSE68086_sample_sheet.csv"             # written by 01, has the patient column
ALPHA     <- 0.05                                 # significance level
FISHER_WS <- 2e8                                  # memory for the exact Fisher calculation
B_MC      <- 1e5                                  # Monte Carlo draws, only if exact Fisher fails
SEED      <- 1                                    
IN_VSD    <- "data/derived/vsd_blind.rds"         # written by 04: vst values used in the main tests
B_PERM    <- 1999                                 # shuffles per main test (smallest p = 1/2000)
SAVE   <- TRUE                                    # write the tables as .csv
OUTDIR <- "data/derived/06b_batch_tests"

# ---- 0. load ------------------------------------------------------------------
stopifnot(
  "subset not found -- run 02 first"        = file.exists(IN),
  "sample sheet not found -- run 01 first"  = file.exists(SHEET)
)
d     <- readRDS(IN) # reads .rds file
ss    <- d$ss # takes sample table from d (list) holds only the samples (Lung, HC, GBM) and columns sample, group, batch.
sheet <- read.csv(SHEET, stringsAsFactors = FALSE) # reads sample list into a data frame. holds the counts.
stopifnot(
  "counts columns and sample sheet rows are not aligned" =
    identical(colnames(d$counts), rownames(ss))
)

# ---- 1. map each sample to its patient and batch -----------------------------------
# match() keeps the row order of ss, so the map stays aligned with the counts columns.
idx <- match(ss$sample, sheet$sample) # matches both tables without changing the order.
stopifnot("some subset samples are missing from the sample sheet" = !anyNA(idx))

sample_map <- data.frame( 
  sample    = ss$sample,
  patient   = sheet$patient[idx],
  condition = ss$group,
  batch     = ss$batch,
  stringsAsFactors = FALSE
)
stopifnot(
  "some samples have no patient ID" = !anyNA(sample_map$patient) && all(sample_map$patient != ""),
  # The batch in the sheet and in the subset must agree, otherwise one of them is stale.
  "batch in the sample sheet does not match the subset" =
    identical(as.character(sample_map$batch), sheet$batch[idx])
)

# ---- 2. batch x condition table: samples and patients per cell ------------------------
# A patient is counted in a cell if they have at least one sample in it.
n_samples  <- table(batch = sample_map$batch, condition = sample_map$condition) # counts samples in each batch x condition pair.
n_patients <- tapply(sample_map$patient,
                     list(batch = sample_map$batch, condition = sample_map$condition),
                     function(p) length(unique(p))) 
# tapply() splits x=sample_map$patients into groups and run f=length(unique(p)). 
# Here it splits patient IDs by batch × condition and counts unique patients in each cell
n_patients[is.na(n_patients)] <- 0L   # empty cells come back as NA from tapply

# Add row and column totals. Patient totals count each patient once, so they can be
# smaller than the sum of the cells if a patient appears in more than one cell.
add_totals <- function(tab, total_row, total_col, total_all) {
  tab <- cbind(tab, Total = total_row)
  rbind(tab, Total = c(total_col, total_all))
}
uniq <- function(p) length(unique(p))
samples_tot  <- add_totals(unclass(n_samples),
                           rowSums(n_samples), colSums(n_samples), sum(n_samples))
patients_tot <- add_totals(n_patients,
                           tapply(sample_map$patient, sample_map$batch, uniq),
                           tapply(sample_map$patient, sample_map$condition, uniq),
                           uniq(sample_map$patient))

# One table for reading: "samples / patients" in each cell.
both <- matrix(paste(samples_tot, "/", patients_tot),
               nrow = nrow(samples_tot), dimnames = dimnames(samples_tot))

cat("===== Batch x condition: samples / patients in each cell =====\n")
print(noquote(both))

# ---- 3. are patients split across batches? ----------------------------------------
# For each patient: how many samples, and in how many different batches.
# per_patient <- data.frame(
#   patient   = names(table(sample_map$patient)),
#   n_samples = as.integer(table(sample_map$patient)),
#   n_batches = as.integer(tapply(sample_map$batch, sample_map$patient,
#                                 function(b) length(unique(b)))),
#   batches   = as.character(tapply(as.character(sample_map$batch), sample_map$patient,
#                                   function(b) paste(sort(unique(b)), collapse = ", "))),
#   stringsAsFactors = FALSE
# )
# split_pat <- per_patient[per_patient$n_batches > 1, ]
# multi_pat <- per_patient[per_patient$n_samples > 1, ]

# cat("\n===== Patients =====\n")
# cat(sprintf("Samples: %d   Patients: %d\n", nrow(sample_map), nrow(per_patient)))
# cat(sprintf("Patients with more than one sample:        %d\n", nrow(multi_pat)))
# cat(sprintf("Patients with samples in more than one batch: %d\n", nrow(split_pat)))
# if (nrow(split_pat) > 0) {
#   print(split_pat, row.names = FALSE)
# } else {
#   cat("No patient is split across batches.\n")
# }
# if (nrow(multi_pat) == 0) {
#   cat("Every patient has exactly one sample: samples and patients are the same units.\n")
# }

# ---- 4. design balance check: Fisher's exact test on batch x condition --------------
# Secondary check, NOT the main test. It asks whether the condition mix is the same in
# every batch. It does not look at the counts at all.
# Why report it: if the mix differs between batches, batch and condition partly move
# together. Then a difference between batches in the main test could come from the
# condition, and a difference between conditions could come from the batch -- so the
# main result has to be read with that in mind.
#
# Uses the sample counts. Here samples = patients (step 3), so both tables give the
# same answer. Fisher works on the table itself, with no large-sample approximation,
# so small cells (Lung in Batch03 = 7) are not a problem.
tab <- n_samples
fis <- tryCatch(
  fisher.test(tab, workspace = FISHER_WS),
  error = function(e) {   # table too large for the exact calculation
    set.seed(SEED)
    fisher.test(tab, simulate.p.value = TRUE, B = B_MC) # FISHER TEST
  }
)

# Where the imbalance is: adjusted residuals compare each cell with the count expected
# if the mix were the same in every batch. |value| > 2 = that cell is far from it.
resid <- suppressWarnings(chisq.test(tab, correct = FALSE))$stdres # CHI-SQUARED TEST for adjusted residuals

cat("\n===== Design balance check (secondary): Fisher's exact test =====\n")
cat("H0: the condition mix is the same in every batch.\n")
cat("Condition mix within each batch (row %):\n")
print(round(100 * prop.table(tab, 1)))
cat(sprintf("\nFisher: p = %.3g  (%s)\n", fis$p.value, fis$method)) 
  # Fisher: p = 0.138  (Fisher's Exact Test for Count Data)
cat("Adjusted residuals (|value| > 2 = cell far from a balanced design):\n")
print(round(resid, 1))
if (fis$p.value < ALPHA) {
  cat("Batch and condition are NOT balanced: read the main test with this in mind.\n")
} else {
  cat("No evidence of imbalance between batch and condition.\n")
}
 # No evidence of imbalance between batch and condition.

balance <- data.frame(
  test      = "Fisher exact, batch x condition (samples)",
  method    = fis$method,
  n         = sum(tab),
  p_value   = fis$p.value,
  max_resid = max(abs(resid))
)

# ===================================================================================
# MAIN TESTS: is batch independent of the expression data, given the condition?
#   H0:  batch independent of X | condition
# Three tests, all on the same data and with the same kind of permutation:
#   6. PERMANOVA      -- do the batch centres differ? and their spread (dispersion)?
#   7. HSIC           -- any dependence at all between batch and X (overall test)
#   8. MMD            -- one pair of batches at a time: do the two differ?
# ===================================================================================

# ---- 5. data and permutation scheme for the main tests --------------------------
# X = vst values from 04: the 6,994 genes kept by filterByExpr, with TMM size factors,
# on a log-like scale. Same values as the PCA in 05, so the tests describe what the
# PCA plots show. Without the log-like scale, distances between samples would be
# decided by the few most abundant genes.
stopifnot("vsd_blind.rds not found -- run 04 first" = file.exists(IN_VSD))
vsd <- readRDS(IN_VSD)
X   <- t(SummarizedExperiment::assay(vsd))   # load samples (rows) x genes (columns) *needs transpose*
stopifnot("samples in vsd are not the same, or not in the same order, as the subset" =
            identical(rownames(X), sample_map$sample))

cond  <- droplevels(factor(sample_map$condition)) # remove any condition without samples
batch <- droplevels(factor(sample_map$batch)) # remove any batch without samples
n     <- nrow(X) # number of samples

# Conditional permutation: shuffle the batch labels (test variable) INSIDE each condition (defined group 
# by the variable i am conditioning on. 
# Tested: batch. Its link to X is the one the shuffle has to break.
# Conditioned on: condition. It stays fixed and defines the groups the shuffle happens in.
# The batch x condition table stays exactly as it is, 
# so only the link between batch and X (for a given condition) is broken. 
# This is what makes all three tests conditional.
# null hypothesis is batch ⊥ X | condition

# Why shuffling batch inside each condition matches that null: take the 40 Lung samples, 
# of which 19 are in Batch02, 7 in Batch03 and 14 in Batch04. If batch has no effect 
# once you know the condition, then which Lung samples ended up in Batch02 has nothing to do 
# with their expression.
shuffle_within <- function(b, cond) {
  for (c in levels(cond)) {
    i <- which(cond == c)
    b[i] <- b[i][sample.int(length(i))] # random reordering of samples in which(cond ==c)
  }
  b
}
# Permutation p-value: share of shuffles with a statistic at least as large as the
# observed one. The +1 counts the observed data as one of the shuffles, so p is never 0.
perm_p <- function(obs, null) (1 + sum(null >= obs)) / (length(null) + 1)

# ---- 6. PERMANOVA: centres (location) and spread (dispersion) --------------------
# Location: pseudo-F for batch AFTER condition. With Euclidean distance, PERMANOVA is
# the same as comparing sums of squares in gene space, so it is computed from X directly:
#   RSS0 = squared distance of each sample to its condition mean        (~ condition)
#   RSS1 = squared distance of each sample to its condition-and-batch fit (~ condition + batch)
#   F    = [(RSS0 - RSS1) / df_batch] / [RSS1 / df_residual]
# A large F means knowing the batch moves the samples' centre, beyond the condition.
pseudo_F <- function(Y, b, cond) {
  d0 <- model.matrix(~ cond) # design for a model with condition only
  d1 <- model.matrix(~ cond + b) # design for a model with condition + batch
  rss0 <- sum(qr.resid(qr(d0), Y)^2) # squared residual for condition only
  rss1 <- sum(qr.resid(qr(d1), Y)^2) # squared residual for condition + batch
  c(F  = ((rss0 - rss1) / (ncol(d1) - ncol(d0))) / (rss1 / (nrow(Y) - ncol(d1))),
    R2 = (rss0 - rss1) / sum(scale(Y, scale = FALSE)^2))   # share of total variation due to batch
}

# Dispersion: distance of each sample to the centre of its condition x batch cell.
# If one batch is more spread out, its samples sit further from their centres. Tested
# with the same pseudo-F on these distances. This matters because PERMANOVA location
# can also turn significant when only the spread differs.
dist_to_centre <- function(Y, b, cond) {
  cell    <- interaction(cond, b, drop = TRUE)
  centres <- rowsum(Y, cell)
  centres <- centres / as.vector(table(cell)[rownames(centres)])
  sqrt(rowSums((Y - centres[as.character(cell), ])^2))
}

set.seed(SEED)
loc_obs  <- pseudo_F(X, batch, cond)
disp_obs <- pseudo_F(matrix(dist_to_centre(X, batch, cond)), batch, cond)
loc_null <- disp_null <- numeric(B_PERM)
for (k in seq_len(B_PERM)) {
  bk <- shuffle_within(batch, cond)
  loc_null[k]  <- pseudo_F(X, bk, cond)["F"]
  # the cells change with the shuffled labels, so the distances are recomputed
  disp_null[k] <- pseudo_F(matrix(dist_to_centre(X, bk, cond)), bk, cond)["F"]
}

# ---- 7. HSIC: overall dependence between batch and X, given condition -------------
# HSIC measures any kind of dependence between two variables through kernels
# (similarity matrices). It is 0 when they are independent.
#   K = Gaussian kernel on X: similarity of two samples' expression profiles (takes 0 when two variables are truly independent)
#   L = 1 if two samples share a batch, 0 otherwise
# Conditional version (condition has only 3 values): compute HSIC inside each
# condition and add them up, weighted by the condition's size. Batch effects that
# exist only inside the conditions are what counts.
#
### Kernel width (sigma): median distance between samples, over all samples. One sigma
# for all conditions, and the same kernel is reused for MMD in step 8.
D2    <- as.matrix(dist(X))^2 # gives EUCLIDIAN distance between all pairs of samples 
sigma <- median(sqrt(D2[upper.tri(D2)])) # median of those distances
K     <- exp(-D2 / (2 * sigma^2)) # Gaussian kernel: similarity of two samples: 1 for identical profiles, close to 0 for very different

idx_c <- split(seq_len(n), cond) # groups sample by condition
# K centred inside each condition (H K H, with H = I - 1/m). Done once, outside the loop.
Kc <- lapply(idx_c, function(i) { # runs a function on each group and returns a list
  H <- diag(length(i)) - 1 / length(i)
  H %*% K[i, i] %*% H
})
hsic_stat <- function(b) {
  sum(sapply(names(idx_c), function(c) {
    i <- idx_c[[c]]
    L <- outer(b[i], b[i], "==") * 1
    (length(i) / n) * sum(Kc[[c]] * L) / length(i)^2
  }))
}

hsic_obs  <- hsic_stat(batch)
hsic_null <- replicate(B_PERM, hsic_stat(shuffle_within(batch, cond)))

# ---- 8. MMD: one pair of batches at a time -----------------------------------------
# MMD (maximum mean discrepancy) compares two groups of samples through the same
# kernel K: mean similarity inside batch a + inside batch b - 2 x between a and b.
# It is 0 when the two batches have the same distribution.
# Conditional version: MMD between the two batches inside each condition, then a
# weighted sum. Shuffles move labels only between the two batches, inside each condition.
mmd2 <- function(Ks, is_a) {
  mean(Ks[is_a, is_a]) + mean(Ks[!is_a, !is_a]) - 2 * mean(Ks[is_a, !is_a])
}
mmd_stat <- function(b, cond_s, K_s, a) {
  w <- table(cond_s)
  sum(sapply(levels(cond_s), function(c) {
    i <- which(cond_s == c)
    (w[[c]] / sum(w)) * mmd2(K_s[i, i], b[i] == a)
  }))
}

pairs <- combn(levels(batch), 2)
mmd_rows <- list()
for (j in seq_len(ncol(pairs))) {
  s      <- batch %in% pairs[, j]          # samples in the two batches of this pair
  b_s    <- droplevels(batch[s])
  cond_s <- droplevels(cond[s])
  K_s    <- K[s, s]
  obs    <- mmd_stat(b_s, cond_s, K_s, pairs[1, j])
  null   <- replicate(B_PERM, mmd_stat(shuffle_within(b_s, cond_s), cond_s, K_s, pairs[1, j]))
  mmd_rows[[j]] <- data.frame(test = "MMD", comparison = paste(pairs[, j], collapse = " vs "),
                              statistic = obs, effect = NA, p_value = perm_p(obs, null))
}
mmd_res <- do.call(rbind, mmd_rows)
mmd_res$p_holm <- p.adjust(mmd_res$p_value, method = "holm")   # 3 pairs -> correct for 3 tests

# ---- 9. report the main tests --------------------------------------------------------
main_tests <- rbind(
  data.frame(test = "PERMANOVA location", comparison = "all batches",
             statistic = loc_obs[["F"]], effect = loc_obs[["R2"]],
             p_value = perm_p(loc_obs[["F"]], loc_null), p_holm = NA),
  data.frame(test = "PERMANOVA dispersion", comparison = "all batches",
             statistic = disp_obs[["F"]], effect = NA,
             p_value = perm_p(disp_obs[["F"]], disp_null), p_holm = NA),
  data.frame(test = "HSIC", comparison = "all batches",
             statistic = hsic_obs, effect = NA,
             p_value = perm_p(hsic_obs, hsic_null), p_holm = NA),
  mmd_res
)

cat("\n===== Main tests: batch independent of X | condition =====\n")
cat(sprintf("X = vst, %d genes x %d samples | %d shuffles of batch inside each condition\n",
            ncol(X), n, B_PERM))
cat("statistic: pseudo-F (PERMANOVA), HSIC, or MMD^2 | effect: R2 of batch (PERMANOVA location)\n")
cat("p_holm: only for the 3 MMD pairs\n")
print(transform(main_tests, statistic = signif(statistic, 3), effect = signif(effect, 3),
                p_value = signif(p_value, 3), p_holm = signif(p_holm, 3)), row.names = FALSE)




cat(sprintf("Smallest possible p with %d shuffles: %.4f\n", B_PERM, 1 / (B_PERM + 1)))
cat("Balance check (step 4) for reading these:", sprintf("Fisher p = %.3g\n", fis$p.value))

# ---- 10. save ---------------------------------------------------------------------
if (SAVE) {
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(sample_map,   file.path(OUTDIR, "sample_patient_batch_map.csv"), row.names = FALSE)
  write.csv(samples_tot,  file.path(OUTDIR, "batch_x_condition_samples.csv"))
  write.csv(patients_tot, file.path(OUTDIR, "batch_x_condition_patients.csv"))
  # write.csv(per_patient,  file.path(OUTDIR, "per_patient_batches.csv"), row.names = FALSE)  # step 3 is off
  write.csv(balance,      file.path(OUTDIR, "design_balance_fisher.csv"), row.names = FALSE)
  write.csv(main_tests,   file.path(OUTDIR, "main_tests.csv"), row.names = FALSE)
  cat("\nTables saved in", OUTDIR, "\n")
}
