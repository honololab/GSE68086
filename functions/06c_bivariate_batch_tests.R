# 06c -- More bivariate tests of batch against the data (follows 06b).
#
# Planned tests (not written yet):
#   - Chi-squared on the batch x condition table (alongside the Fisher test in 06b).
#   - Mann-Whitney: two batches at a time, inside each condition, on each gene or PC.
#   - Kruskal-Wallis: all batches at once, inside each condition, on each gene or PC.
#
# Run from the project root:  Rscript "functions/06c_bivariate_batch_tests.R"

library(tidyverse)

# ================== CHOICES -- edit here ==================
IN     <- "data/derived/subset_lung_gbm_hc.rds"   # written by 02
SHEET  <- "GSE68086_sample_sheet.csv"             # written by 01, has the patient column
IN_VSD <- "data/derived/vsd_blind.rds"            # written by 04: vst values to test
ALPHA  <- 0.05                                    # significance level
SEED   <- 1

# chi-squared
CHISQ_CORRECT <- FALSE                            # Yates correction (only used for 2x2 tables)
CHISQ_SIM     <- FALSE                            # TRUE = Monte Carlo p, if expected counts are small
B_MC          <- 1e5                              # Monte Carlo draws, only if CHISQ_SIM = TRUE

# Mann-Whitney / Kruskal-Wallis
UNIT     <- "pc"                                  # "pc" = test each PC, "gene" = test each gene
N_PC     <- 10                                    # number of PCs to test, if UNIT = "pc"
P_ADJUST <- "BH"                                  # correction across PCs/genes ("holm", "BH", ...)

SAVE   <- TRUE                                    # write the tables as .csv
OUTDIR <- "data/derived/06c_batch_tests"

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

class(ss)
head(ss)
class(sheet)
head(sheet)

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

class(sample_map) # data.frame
head(sample_map)
dim(sample_map) # sample patient condition   batch
sample_map[1, ]

# ---- 2. batch x condition table: samples and patients per cell ------------------------
# A patient is counted in a cell if they have at least one sample in it.
n_samples  <- table(batch = sample_map$batch, condition = sample_map$condition) # counts samples in each batch x condition pair.
n_samples
n_patients <- tapply(sample_map$patient,
                     list(batch = sample_map$batch, condition = sample_map$condition),
                     function(p) length(unique(p)))
class(n_patients)
n_patients
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
samples_tot
patients_tot
# One table for reading: "samples / patients" in each cell.
both <- matrix(paste(samples_tot, "/", patients_tot),
               nrow = nrow(samples_tot), dimnames = dimnames(samples_tot))

cat("===== Batch x condition: samples / patients in each cell =====\n")
print(noquote(both))
# each patient is only one sample.

# ---- 3. join sample_map and vsd (filtered and normalised genes)
vsd <- readRDS(IN_VSD)
m <- SummarizedExperiment::assay(vsd) # extract numerical matrix data from an object in Bioconductor
colnames(m) # samples
stopifnot("vst columns and sample_map are not aligned" =
  identical(colnames(m), sample_map$sample))
sample_map$value <- m["ENSG00000000419", ] # value = one gene's vst values, one per sample. see whether a gene i care about survived the filter

### TESTING

## Kruskal Wallis
# H₀ (Hipótesis Nula): No existen diferencias significativas en la distribución.
# H₁ (Hipótesis Alternativa): Al menos uno de los grupos presenta una distribución distinta.
# Nivel de significancia: α = 0.05 Criterio de decisión: Si p-value < 0.05, rechazo H₀.
hc <- subset(sample_map, condition == "HC")
lung <- subset(sample_map, condition == "Lung")
gbm <- subset(sample_map, condition == "GBM")
kruskal.test(value ~ batch, data = hc)
kruskal.test(value ~ batch, data = lung)
kruskal.test(value ~ batch, data = gbm)

## Mann Whitney
# H₀ (Hipótesis Nula): No existen diferencias significativas en la distribución.
# H₁ (Hipótesis Alternativa): Existen diferencias significativas en la distribución.
# Nivel de significancia que elijo: α = 0.05 Criterio de decisión: Si p-value < 0.05, rechazo H₀.
pair <- subset(sample_map, condition == "HC" & batch %in% c("Batch02", "Batch03"))
wilcox.test(value ~batch, data = pair)

# function that runs it for all conditions (hc, lung, gbm) and batches (02, 03, 04).
# Mann-Whitney for every pair of batches, inside each condition.
# Returns one row per test: condition, the two batches, their sizes, W and p-value.
mw_by_condition <- function(sample_map) {
  out <- list()
  for (cond in sort(unique(sample_map$condition))) {
    sub   <- sample_map[sample_map$condition == cond, ]
    pairs <- combn(sort(unique(as.character(sub$batch))), 2)   # one column per pair
    for (k in seq_len(ncol(pairs))) {
      a <- sub$value[sub$batch == pairs[1, k]]
      b <- sub$value[sub$batch == pairs[2, k]]
      t <- wilcox.test(a, b)                                   # unpaired = Mann-Whitney
      out[[length(out) + 1]] <- data.frame(
        condition = cond, batch_a = pairs[1, k], batch_b = pairs[2, k],
        n_a = length(a), n_b = length(b), W = unname(t$statistic), p = t$p.value)
    }
  }
  do.call(rbind, out)
}

mw_by_condition(sample_map)

