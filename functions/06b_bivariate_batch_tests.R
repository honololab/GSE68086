# 06b -- One batch at a time: do the conditions differ inside the same batch?
#
# Section 2.3 of the report ("Estrategia: un lote a la vez").
#
# All samples in one batch were processed together. So if HC, GBM and Lung differ
# INSIDE the same batch, that difference cannot come from the batch.
#
# Notation (same as the report):
#   Y = a feature measured on each sample (depth or zero share)
#   C = condition (HC / GBM / Lung)      B = batch
#
# For each batch b and each feature Y the null hypothesis is
#   H0(b): Y independent of C | B = b
# i.e. inside batch b, Y has the same distribution in the three conditions.
# Tested with Kruskal-Wallis. 3 batches x 2 features = 6 tests, Holm-corrected.
#
# One batch gives a local answer. Repeating it in the three batches shows if the
# pattern is CONSISTENT:
#   Same order of medians in all batches -> evidence of an effect of the condition.
#   Difference in one batch only         -> chance, low power, or batch x condition.
#   Order reversed between batches       -> the all-samples pattern may come from how
#                                           the samples were split across batches.
# Batch03 is the reference: it is printed first and the order of the other batches
# is compared with it.
#
# Small groups (Lung in Batch03 has 7 samples) can give a large p-value even when
# there is a real difference, so look at the medians and the order, not only at p.
#
# Base R only. Needs the .rds written by 02.
# Run from the project root:  Rscript "functions/06b_bivariate_batch_tests.R"

# ================== CHOICES -- edit here ==================
IN        <- "data/derived/subset_lung_gbm_hc.rds"   # written by 02
REF_BATCH <- "Batch03"                  # the batch the others are compared with
FEATURES  <- c("lib_size", "zero_share")  # sample features to test (built in step 1)
P_ADJUST  <- "holm"                     # multiple-testing correction, see ?p.adjust
ALPHA     <- 0.05                       # only used to mark rows with "*"
SAVE      <- TRUE                       # write the table as .csv
OUTDIR    <- "data/derived/06b_bivariate_batch_test"

# ---- 0. load ------------------------------------------------------------------
# d$counts: genes (rows) x samples (columns). d$ss: one row per sample, in the same
# order as the columns of counts (02 checked this before saving).
stopifnot("subset not found -- run 02 first" = file.exists(IN))
d      <- readRDS(IN)
counts <- d$counts
ss     <- d$ss
stopifnot(
  "counts columns and sample sheet rows are not aligned" =
    identical(colnames(counts), rownames(ss)),
  "REF_BATCH is not one of the batches in the subset" = REF_BATCH %in% levels(ss$batch)
)

# ---- 1. sample features --------------------------------------------------------
# Same definitions as in 03, so the tests describe what Figure 2 shows.
ss$lib_size   <- colSums(counts)        # depth: total reads in the sample
ss$zero_share <- colMeans(counts == 0)  # share of genes with no reads
stopifnot("unknown name in FEATURES" = all(FEATURES %in% names(ss)))

# How many samples in each condition x batch.
print(table(group = ss$group, batch = ss$batch))

# Put the reference batch first.
batches <- c(REF_BATCH, setdiff(levels(ss$batch), REF_BATCH))

# Median in readable units: depth in millions of reads, zero share in percent.
readable <- function(x, feature) if (feature == "lib_size") x / 1e6 else 100 * x

# ---- 2. inside each batch: are the three conditions the same? --------------------
# Kruskal-Wallis compares ranks, not raw values: it needs no normal distribution and
# one extreme sample cannot dominate it.
rows <- list()
for (f in FEATURES) {
  for (b in batches) {
    y <- ss[[f]][ss$batch == b]
    g <- ss$group[ss$batch == b]
    m <- tapply(y, g, median)
    rows[[length(rows) + 1]] <- data.frame(
      feature = f,
      batch   = b,
      n       = length(y),
      medians = paste(sprintf("%s %.2f", names(m), readable(m, f)), collapse = " | "),
      order   = paste(names(sort(m)), collapse = " < "),
      p_raw   = kruskal.test(y, g)$p.value
    )
  }
}
res       <- do.call(rbind, rows)
res$p_adj <- p.adjust(res$p_raw, method = P_ADJUST)
res$sig   <- ifelse(res$p_adj < ALPHA, "*", "")

# Is the order of the medians the same as in the reference batch?
ref_order       <- setNames(res$order[res$batch == REF_BATCH], res$feature[res$batch == REF_BATCH])
res$same_as_ref <- res$order == ref_order[res$feature]

# ---- 3. print --------------------------------------------------------------------
res_print <- res
res_print$p_raw <- signif(res$p_raw, 3)
res_print$p_adj <- signif(res$p_adj, 3)

cat("\n===== Inside each batch: do the conditions differ? (reference:", REF_BATCH, ") =====\n")
cat("Medians: depth in millions of reads, zero share in %.\n")
print(res_print, row.names = FALSE)
cat(sprintf("'*' = %s-adjusted p < %.2f, over all %d tests.\n",
            P_ADJUST, ALPHA, nrow(res)))

# ---- 4. save ---------------------------------------------------------------------
if (SAVE) {
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  write.csv(res, file.path(OUTDIR, "within_batch.csv"), row.names = FALSE)
  cat("Table saved in", OUTDIR, "\n")
}
