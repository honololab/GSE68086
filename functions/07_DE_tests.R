# 07 -- Differential expression: fit the model, extract the contrasts, shrink, save
#
# Reads the filtered, normalised object 04 saved (6,994 genes x 133 samples, design
# ~ batch + group) and does the inference. Steps, in this order:
#   1. DESeq()            fits the negative-binomial GLM for every gene
#   2. resultsNames()     the five coefficients (betas) the model estimated
#   3. results()          one table per contrast: Lung vs HC, GBM vs HC, Lung vs GBM
#   4. summary()          how many genes move, and how many were set aside
#   5. lfcShrink()        a second, more honest log fold change for ranking and plots
#   6. tables to disk     one CSV per contrast, ordered by padj
#   7. MA plots           one per contrast, raw vs shrunken side by side
#
# Run from the project root:  Rscript "functions/07_DE_tests.R"

library(DESeq2)
library(apeglm)      # shrinkage method for named coefficients
library(ashr)        # shrinkage method for contrasts that are not a single coefficient

IN     <- "data/derived/dds_filtered.rds"     # written by 04
OUTDIR <- "data/derived/de"                   # one CSV per contrast
FIGDIR <- "data/derived/figures"
ALPHA  <- 0.05                                # padj threshold used everywhere below

stopifnot("dds not found -- run 04 first" = file.exists(IN))
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(IN)
stopifnot(nrow(dds) < 10000)                  # a filtered object, not the 57,736-gene one
print(dds)

# ---- 1. fit -----------------------------------------------------------------------
# DESeq() does three things in sequence: size factors (already there from 04, kept),
# a dispersion per gene (shrunk towards a trend over all genes), and the GLM fit with
# a Wald test on each coefficient. Takes a minute or two on 6,994 genes.
dds <- DESeq(dds)

# ---- 2. the coefficients ----------------------------------------------------------
# One beta per column of the design matrix. HC and Batch02 are the reference levels
# (set in 02), so the names read "<level> vs <reference>".
print(resultsNames(dds))
# [1] "Intercept" "batch_Batch03_vs_Batch02" "batch_Batch04_vs_Batch02"
# [4] "group_GBM_vs_HC" "group_Lung_vs_HC"

# ---- 3. one table per contrast ----------------------------------------------------
# Always name the contrast. results(dds) with no arguments returns the LAST coefficient,
# which is easy to misread. alpha only changes what summary() counts and how the
# independent filtering is tuned; the padj column itself is the same.
#
# Lung vs HC and GBM vs HC are single coefficients -> name=.
# Lung vs GBM is the difference of two coefficients -> contrast=c(variable, num, denom).
res <- list(
  Lung_vs_HC  = results(dds, name = "group_Lung_vs_HC",  alpha = ALPHA),
  GBM_vs_HC   = results(dds, name = "group_GBM_vs_HC",   alpha = ALPHA),
  Lung_vs_GBM = results(dds, contrast = c("group", "Lung", "GBM"), alpha = ALPHA)
)

# Columns, for the record:
#   baseMean        mean normalised count over all 133 samples (how expressed the gene is)
#   log2FoldChange  the beta of the contrast: 1 = twice as much, -1 = half
#   lfcSE           its standard error
#   stat            Wald statistic = log2FoldChange / lfcSE
#   pvalue          from the Wald test
#   padj            Benjamini-Hochberg adjusted p-value; this is the one to filter on
#   NA in padj      gene set aside: too few counts for the test, or one extreme sample (Cook's)

# ---- 4. how much moves ------------------------------------------------------------
for (k in names(res)) {
  cat(sprintf("\n=== %s ===\n", k))
  summary(res[[k]])
  cat(sprintf("padj < %.2f: %d genes (%d up, %d down)\n", ALPHA,
              sum(res[[k]]$padj < ALPHA, na.rm = TRUE),
              sum(res[[k]]$padj < ALPHA & res[[k]]$log2FoldChange > 0, na.rm = TRUE),
              sum(res[[k]]$padj < ALPHA & res[[k]]$log2FoldChange < 0, na.rm = TRUE)))
}
# Expect the two cancer-vs-HC lists to be long: the exploratory section showed HC
# libraries differ in composition, not only in depth. Lung vs GBM is the cleaner contrast.

# ---- 5. shrunken log fold changes -------------------------------------------------
# The raw log2FoldChange of a low-count gene is noisy: 2 vs 8 reads is "four-fold" but
# could easily have been 4 vs 5. lfcShrink() pulls each estimate towards zero in
# proportion to its uncertainty, so abundant genes with solid moderate changes rank
# above sparse genes with huge, unreliable ones. It does NOT change pvalue or padj --
# those still come from the unshrunken test. Use padj to decide, shrunken LFC to rank.
#
# apeglm works on a named coefficient; for a contrast between two non-reference levels
# it cannot be used, so Lung vs GBM uses ashr instead.
shr <- list(
  Lung_vs_HC  = lfcShrink(dds, coef = "group_Lung_vs_HC", type = "apeglm"),
  GBM_vs_HC   = lfcShrink(dds, coef = "group_GBM_vs_HC",  type = "apeglm"),
  Lung_vs_GBM = lfcShrink(dds, contrast = c("group", "Lung", "GBM"), type = "ashr")
)

# ---- 6. tables to disk ------------------------------------------------------------
# One row per gene: raw and shrunken LFC side by side, ordered by padj. The Ensembl id
# is the row name; symbols can be added later with org.Hs.eg.db (00c shows how).
for (k in names(res)) {
  tab <- data.frame(
    gene        = rownames(res[[k]]),
    baseMean    = res[[k]]$baseMean,
    lfc_raw     = res[[k]]$log2FoldChange,
    lfc_shrunk  = shr[[k]]$log2FoldChange,
    lfcSE       = res[[k]]$lfcSE,
    pvalue      = res[[k]]$pvalue,
    padj        = res[[k]]$padj
  )
  tab <- tab[order(tab$padj, -abs(tab$lfc_shrunk)), ]
  out <- file.path(OUTDIR, paste0("de_", k, ".csv"))
  write.csv(tab, out, row.names = FALSE)
  cat(sprintf("written %s (%d genes, %d with padj < %.2f)\n",
              out, nrow(tab), sum(tab$padj < ALPHA, na.rm = TRUE), ALPHA))
  cat("top 10 by padj, then |shrunken LFC|:\n")
  print(head(tab[!is.na(tab$padj), c("gene", "baseMean", "lfc_raw", "lfc_shrunk", "padj")], 10),
        row.names = FALSE, digits = 3)
}

# ---- 7. MA plots: raw vs shrunken -------------------------------------------------
# x = baseMean (log scale), y = log2FoldChange, blue = padj < alpha. On the left the
# low-count genes fan out at the far left; on the right the shrinkage pulls that fan
# back towards zero while leaving the well-supported genes where they were.
for (k in names(res)) {
  png(file.path(FIGDIR, paste0("12_MA_", k, ".png")), width = 1800, height = 800, res = 150)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  plotMA(res[[k]], ylim = c(-4, 4), alpha = ALPHA, main = paste(k, "-- raw LFC"))
  plotMA(shr[[k]], ylim = c(-4, 4), alpha = ALPHA, main = paste(k, "-- shrunken LFC"))
  dev.off()
}
cat(sprintf("\nMA plots written to %s\n", FIGDIR))

saveRDS(dds, "data/derived/dds_fitted.rds")   # the fitted object, in case a later script needs it
