# 04b -- One batch at a time: does a sample feature follow the batch or the condition?
#
# ---- the question ------------------------------------------------------------
# Figure 2 of the report suggests two things: sequencing depth follows the BATCH, and
# the share of genes with zero reads follows the CONDITION (HC / GBM / Lung). In those
# plots both sources of variation are mixed together, so looking at all samples at
# once cannot tell which one explains a difference. This script separates them by
# STRATIFYING: it fixes one variable and compares only inside it.
#
# Notation (same as section 2.3 of the report):
#   Y = a feature measured on each sample (depth or zero share)
#   C = condition (group)      B = batch
#
# ---- the three runs, in order --------------------------------------------------
#   A. Inside each batch, compare the conditions.      H0: Y independent of C | B = b
#      All samples in one batch were processed together, so a difference between
#      HC, GBM and Lung there cannot come from the batch.
#      Question answered: "does the CONDITION change Y?"
#
#   B. Inside each condition, compare every batch with the reference batch (Batch03).
#      H0: Y has the same distribution in batch b and in Batch03, within condition c.
#      Question answered: "does the BATCH change Y?"  (the mirror of A)
#
#   C. Interaction check, all samples in one model. Does the batch shift every
#      condition by the same amount? If not, adding batch as a plain covariate in the
#      later models does not fully remove it.
#
# ---- how to read the output ----------------------------------------------------
#   A feature that is TECHNICAL (follows the batch):   A mostly not significant,
#                                                      B significant.
#   A feature that is BIOLOGICAL (follows the group):  A significant, with the SAME
#                                                      ordering in the three batches,
#                                                      B mostly not significant.
#   Look at the direction and the medians, not only at the p-values: each batch has
#   only 37 to 49 samples (Lung in Batch03 has 7), so a non-significant test does not
#   prove there is no difference.
#
# ---- why zeros get a second, depth-adjusted test --------------------------------
#   With fewer reads, more genes get zero reads just by chance. If depth follows the
#   batch and zeros follow depth, part of the zero pattern could be technical. So run
#   A also asks: do the conditions differ in zeros AT EQUAL DEPTH?
#
# Base R only. Needs the .rds written by 02.
# Run from the project root:  Rscript "functions/04b_within_batch_tests.R"

# ================== CHOICES -- edit here ==================
IN         <- "data/derived/subset_lung_gbm_hc.rds"   # written by 02
REF_BATCH  <- "Batch03"                  # every batch is compared with this one in run B
FEATURES   <- c("lib_size", "zero_share")  # sample features to test (built in step 1)
P_ADJUST   <- "holm"                     # multiple-testing correction, see ?p.adjust
ALPHA      <- 0.05                       # only used to flag rows with "*" in the output
MIN_N      <- 3                          # smallest group x batch cell allowed
SAVE       <- TRUE                       # write the result tables as .csv
OUTDIR     <- "data/derived/04b_within_batch"

# ---- 0. load ---------------------------------------------------------------------
# d$counts: genes (rows) x samples (columns). d$ss: one row per sample, with group and
# batch, in the same order as the columns of counts (02 checked this before saving).
stopifnot("subset not found -- run 02 first" = file.exists(IN))
d      <- readRDS(IN)
counts <- d$counts
ss     <- d$ss
stopifnot(
  "counts columns and sample sheet rows are not aligned" =
    identical(colnames(counts), rownames(ss)),
  all(c("group", "batch") %in% names(ss)),
  "REF_BATCH is not one of the batches in the subset" = REF_BATCH %in% levels(ss$batch)
)

# ---- 1. build the sample features ----------------------------------------------------
# Same definitions as in 03, so the tests describe exactly what Figure 2 shows.
ss$lib_size   <- colSums(counts)        # depth: total reads in the sample
ss$n_zero     <- colSums(counts == 0)   # genes with no reads
ss$n_nonzero  <- colSums(counts > 0)    # genes with at least one read
ss$zero_share <- colMeans(counts == 0)  # n_zero as a share of all 57,736 genes

stopifnot("unknown name in FEATURES" = all(FEATURES %in% names(ss)))

# Every group x batch cell needs samples, or runs A and B have nothing to compare.
cells <- table(group = ss$group, batch = ss$batch)
print(cells)
stopifnot("a group x batch cell is smaller than MIN_N" = all(cells >= MIN_N))

# The batch factor with REF_BATCH first. R reads the first level of a factor as the
# baseline, so in run C every batch coefficient reads "<batch> vs Batch03".
ss$batch_ref <- relevel(ss$batch, ref = REF_BATCH)

# ---- 2. helpers ---------------------------------------------------------------------

# Ordering of the group medians, e.g. "HC < GBM < Lung". This is the DIRECTION of a
# difference, the thing to compare between batches in run A.
direction <- function(y, g) {
  m <- tapply(y, g, median)
  paste(names(sort(m)), collapse = " < ")
}

# Group medians as one readable string, e.g. "HC 2.52 | GBM 2.06 | Lung 2.64".
# Depth is shown in millions of reads, zero share in percent.
show_medians <- function(y, g, feature) {
  m <- tapply(y, g, median)
  m <- if (feature == "lib_size") m / 1e6 else 100 * m
  paste(sprintf("%s %.2f", names(m), m), collapse = " | ")
}

# Kruskal-Wallis: are the three groups the same? It compares RANKS, not raw values,
# so it needs no normality and one extreme sample cannot dominate it. Taking the log
# of depth would give the same p-value, because the log does not change the ranks.
kw_p <- function(y, g) kruskal.test(x = y, g = g)$p.value

# Zeros at equal depth. The response is, for each sample, "n_zero out of all genes",
# i.e. a proportion, so the natural model is a binomial GLM. It is "quasi" because
# genes are not independent coin flips: the real spread is larger than a binomial
# allows, and quasibinomial estimates that extra spread from the data instead of
# assuming it. log(lib_size) goes in first, so the group term is tested on what
# depth leaves unexplained. The F-test compares the model with and without group.
zeros_at_equal_depth <- function(dd) {
  without_group <- glm(cbind(n_zero, n_nonzero) ~ log(lib_size),
                       family = quasibinomial, data = dd)
  with_group    <- glm(cbind(n_zero, n_nonzero) ~ log(lib_size) + group,
                       family = quasibinomial, data = dd)
  p <- anova(without_group, with_group, test = "F")$`Pr(>F)`[2]
  # Odds ratio of a zero vs HC, at equal depth: > 1 means more zeros than HC.
  or <- exp(coef(with_group)[grep("^group", names(coef(with_group)))])
  names(or) <- sub("^group", "", names(or))
  list(p = p, or = paste(sprintf("%s %.2f", names(or), or), collapse = " | "))
}

# Mark rows under ALPHA after correction, just to make the tables easier to scan.
flag <- function(p) ifelse(p < ALPHA, "*", "")

# ---- 3. run A: inside each batch, compare the conditions ------------------------------
# One family of tests = every batch x every feature (3 x 2 = 6). The correction is
# applied over the whole family, because all six answer the same question.
rows <- list()
for (b in levels(ss$batch)) {
  dd <- ss[ss$batch == b, ]
  for (f in FEATURES) {
    rows[[length(rows) + 1]] <- data.frame(
      batch     = b,
      feature   = f,
      n         = nrow(dd),
      medians   = show_medians(dd[[f]], dd$group, f),
      direction = direction(dd[[f]], dd$group),
      p_raw     = kw_p(dd[[f]], dd$group)
    )
  }
}
runA <- do.call(rbind, rows)
runA$p_adj <- p.adjust(runA$p_raw, method = P_ADJUST)
runA$sig   <- flag(runA$p_adj)

# Zeros adjusted for depth: a separate family of 3 tests (one per batch). It is a
# second look at the zero share, not a new feature, so it is corrected on its own.
rows <- list()
for (b in levels(ss$batch)) {
  z <- zeros_at_equal_depth(ss[ss$batch == b, ])
  rows[[length(rows) + 1]] <- data.frame(batch = b, odds_ratio_vs_HC = z$or, p_raw = z$p)
}
runA_zeros <- do.call(rbind, rows)
runA_zeros$p_adj <- p.adjust(runA_zeros$p_raw, method = P_ADJUST)
runA_zeros$sig   <- flag(runA_zeros$p_adj)

# ---- 4. run B: inside each condition, each batch vs REF_BATCH -----------------------
# Wilcoxon rank-sum (Mann-Whitney): two groups, compared by ranks, for the same
# reasons as Kruskal-Wallis above. Using one reference batch means 2 comparisons per
# condition (Batch02 vs Batch03, Batch04 vs Batch03) instead of every possible pair.
# One family = 3 conditions x 2 batches x 2 features = 12 tests.
# If R warns "cannot compute exact p-value with ties", two samples share a value and
# R switched to the normal approximation; the result is still valid.
other_batches <- setdiff(levels(ss$batch), REF_BATCH)
rows <- list()
for (g in levels(ss$group)) {
  for (b in other_batches) {
    for (f in FEATURES) {
      y_b   <- ss[[f]][ss$group == g & ss$batch == b]
      y_ref <- ss[[f]][ss$group == g & ss$batch == REF_BATCH]
      rows[[length(rows) + 1]] <- data.frame(
        group      = g,
        comparison = paste(b, "vs", REF_BATCH),
        feature    = f,
        n_b        = length(y_b),
        n_ref      = length(y_ref),
        # > 1 means batch b has a higher median than the reference batch
        median_ratio = round(median(y_b) / median(y_ref), 2),
        p_raw      = wilcox.test(y_b, y_ref)$p.value
      )
    }
  }
}
runB <- do.call(rbind, rows)
runB$p_adj <- p.adjust(runB$p_raw, method = P_ADJUST)
runB$sig   <- flag(runB$p_adj)

# ---- 5. run C: does the batch shift every condition by the same amount? --------------
# "Additive" model: group + batch  -> the batch moves all groups by the same amount.
# "Interaction" model: group * batch -> each group can move differently in each batch.
# The test compares the two. A small p-value means the additive model is not enough.
#
# Depth: linear model on log(depth). The log turns "Batch04 has half the reads" into
# a constant shift, which is what "same amount" means for depth. lm assumes roughly
# normal residuals on that scale; check with plot(depth_int) if in doubt.
depth_add <- lm(log(lib_size) ~ group + batch_ref, data = ss)
depth_int <- lm(log(lib_size) ~ group * batch_ref, data = ss)
p_int_depth <- anova(depth_add, depth_int)$`Pr(>F)`[2]

# Zeros: the same quasibinomial model as above, with depth kept in both models.
zero_add <- glm(cbind(n_zero, n_nonzero) ~ log(lib_size) + group + batch_ref,
                family = quasibinomial, data = ss)
zero_int <- glm(cbind(n_zero, n_nonzero) ~ log(lib_size) + group * batch_ref,
                family = quasibinomial, data = ss)
p_int_zero <- anova(zero_add, zero_int, test = "F")$`Pr(>F)`[2]

runC <- data.frame(
  feature = c("lib_size (log)", "zero_share (adjusted for depth)"),
  p_interaction = c(p_int_depth, p_int_zero)
)
runC$sig <- flag(runC$p_interaction)

# ---- 6. print ------------------------------------------------------------------------
fmt <- function(x) { x[] <- lapply(x, function(v) if (is.numeric(v)) signif(v, 3) else v); x }

cat("\n========== A. Inside each batch: do the conditions differ? ==========\n")
cat("Medians: depth in millions of reads, zero share in %.",
    "Compare 'direction' between batches.\n")
print(fmt(runA), row.names = FALSE)

cat("\n---------- A (zeros). Same question at equal depth ----------\n")
cat("odds_ratio_vs_HC > 1: more zeros than HC with the same depth.\n")
print(fmt(runA_zeros), row.names = FALSE)

cat("\n========== B. Inside each condition: does each batch differ from",
    REF_BATCH, "? ==========\n")
print(fmt(runB), row.names = FALSE)

cat("\n========== C. Interaction group x batch (reference:", REF_BATCH, ") ==========\n")
print(fmt(runC), row.names = FALSE)
cat("Coefficients of the depth model with interaction (batch vs", REF_BATCH, "):\n")
print(round(coef(summary(depth_int)), 3))

cat(sprintf("\n'*' = %s-adjusted p < %.2f. Families: A = %d tests, A (zeros) = %d, B = %d.\n",
            P_ADJUST, ALPHA, nrow(runA), nrow(runA_zeros), nrow(runB)))

# ---- 7. save -------------------------------------------------------------------------
if (SAVE) {
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(runA,       file.path(OUTDIR, "A_within_batch.csv"),       row.names = FALSE)
  write.csv(runA_zeros, file.path(OUTDIR, "A_zeros_equal_depth.csv"),  row.names = FALSE)
  write.csv(runB,       file.path(OUTDIR, "B_vs_reference_batch.csv"), row.names = FALSE)
  write.csv(runC,       file.path(OUTDIR, "C_interaction.csv"),        row.names = FALSE)
  cat("Tables saved in", OUTDIR, "\n")
}
