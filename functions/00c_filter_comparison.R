# 00c. Gene-filter comparison on the Lung / GBM / HC, Batch02-04 subset (read-only companion to 00b)
#
# Purpose: decide, with numbers, which low-count filter 04 should use and why, before any DE model.
# It reproduces the subset of 02_make subset.R in memory, applies several candidate rules to the
# same matrix, and reports what each rule keeps and *which group* the disagreements come from.
# Nothing is written to disk. Run from the project root (like the other scripts).
#
# Rules compared
#   filterByExpr(group)   edgeR default: CPM >= 10/median-lib-size*1e6 in >= min-group-size samples
#                         (with the large.n / min.prop tolerance), plus total count >= 15.
#   filterByExpr(design)  same, but min group size comes from the ~ batch + group hat matrix.
#   manual >=10 in >=38   the raw-count rule that 04 used before switching to filterByExpr.
#   CPM>1 in >=38         the rule in subset/make_subset.py.
#   >=5 in >=10%          a "paper-like" rule (Best 2015 excluded transcripts with < 5 reads).
suppressPackageStartupMessages({ library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi) })

counts <- as.matrix(read.csv("GSE68086_TEP_data_matrix.csv", row.names = 1, check.names = FALSE))
ss     <- read.csv("GSE68086_sample_sheet.csv", row.names = 1, check.names = FALSE)   # written by 01
stopifnot(setequal(rownames(ss), colnames(counts)))
counts <- counts[, rownames(ss)]

# ---- 1. subset exactly as 02 -----------------------------------------------------------------
keep_s <- ss$group %in% c("Lung", "GBM", "HC") & ss$batch %in% c("Batch02", "Batch03", "Batch04")
ss  <- droplevels(ss[keep_s, ]); cnt <- counts[, rownames(ss)]
ss$group <- relevel(factor(ss$group), ref = "HC"); ss$batch <- factor(ss$batch)
lib <- colSums(cnt)
cat("=== subset as in 02 ===\n"); print(table(ss$group, ss$batch))
cat("median library size by group:\n"); print(tapply(lib, ss$group, median))
cat("genes with any read in the subset:", sum(rowSums(cnt) > 0), "of", nrow(cnt), "\n")

# ---- 2. the author-excluded GBM sample is still inside this subset ----------------------------
# GEO !Sample_description flags VU398-GBM (and VU258-CRC) as "did not yield sufficient ... reads".
# 01/02 never drop them, so VU398-GBM travels into the Lung/GBM/HC subset.
cat("\n=== VU398-GBM (flagged by the authors) ===\n")
cat("in subset:", "VU398-GBM" %in% colnames(cnt),
    "| lib size:", lib["VU398-GBM"], "(rank", rank(lib)["VU398-GBM"], "of", length(lib), ")",
    "| genes detected:", sum(cnt[, "VU398-GBM"] > 0), "vs GBM median", median(colSums(cnt[, ss$group == "GBM"] > 0)), "\n")

# ---- 3. candidate filters ----------------------------------------------------------------------
y      <- DGEList(cnt, samples = ss, group = ss$group)
design <- model.matrix(~ batch + group, data = ss)
n_min  <- min(table(ss$group))

F <- list(
  "filterByExpr(group)"   = filterByExpr(y, group = ss$group),
  "filterByExpr(design)"  = filterByExpr(y, design = design),
  "manual >=10 in >=38"   = rowSums(cnt >= 10) >= n_min,
  "CPM>1 in >=38 (py)"    = rowSums(cpm(y) > 1) >= n_min,
  ">=5 in >=10% samples"  = rowSums(cnt >= 5) >= ceiling(0.10 * ncol(cnt))
)
cat("\n=== genes kept per rule ===\n"); print(sapply(F, sum))

med_lib <- median(y$samples$lib.size)
cpm_cut <- 10 / med_lib * 1e6
cat("\nfilterByExpr implied CPM cutoff:", round(cpm_cut, 2), "CPM",
    "= ~", round(cpm_cut * median(lib[ss$group == "HC"]) / 1e6, 1), "reads in a median HC library,",
    "~", round(cpm_cut * median(lib[ss$group == "Lung"]) / 1e6, 1), "reads in a median Lung library\n")
cat("filterByExpr(group)  min sample size: n_min =", n_min, "-> tolerance-adjusted", 10 + (n_min - 10) * 0.7, "\n")
cat("filterByExpr(design) min sample size: 1/max(hat) =", round(1 / max(hat(design)), 1),
    " (a batch x group cell, not a group -> more permissive)\n")

# ---- 4. where do filterByExpr(group) and the raw-count rule disagree? --------------------------
a <- F[["filterByExpr(group)"]]; b <- F[["manual >=10 in >=38"]]
cat("\n=== filterByExpr(group) vs manual >=10 in >=38 ===\n"); print(table(filterByExpr = a, manual = b))

logcpm   <- cpm(y, log = TRUE, prior.count = 2)
grp_mean <- sapply(levels(ss$group), function(g) rowMeans(logcpm[, ss$group == g, drop = FALSE]))
det10    <- sapply(levels(ss$group), function(g) rowSums(cnt[, ss$group == g, drop = FALSE] >= 10))
only_fbe <- a & !b
cat("\ngenes kept only by filterByExpr:", sum(only_fbe), "\n")
cat("  mean logCPM by group           :"); print(round(colMeans(grp_mean[only_fbe, , drop = FALSE]), 2))
cat("  median # samples with >=10 reads:"); print(apply(det10[only_fbe, , drop = FALSE], 2, median))
# -> these are genes expressed in the shallow HC libraries; a raw-count threshold drops them
#    simply because HC was sequenced ~half as deep, i.e. the raw rule is biased against the controls.

# ---- 5. composition: how different are HC and cancer libraries? -------------------------------
sym  <- suppressMessages(mapIds(org.Hs.eg.db, keys = rownames(cnt), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first"))
top10 <- head(sort(rowSums(cnt), decreasing = TRUE), 10)
cat("\n=== top-10 genes by total reads ===\n")
print(data.frame(symbol = sym[names(top10)], pct_reads = round(100 * top10 / sum(cnt), 1)))

ribo <- !is.na(sym) & grepl("^RP[SL]\\d", sym)
cat("\nribosomal-protein share of reads, median by group:\n")
print(round(tapply(colSums(cnt[ribo, ]) / lib, ss$group, median), 3))

# mitochondrial genes are intron-less, so intron-spanning counting should leave them at ~0
mt_ids <- c(`MT-ND1` = "ENSG00000198888", `MT-CO1` = "ENSG00000198804", `MT-ND4` = "ENSG00000198886",
            `MT-CYB` = "ENSG00000198727", `MT-ATP6` = "ENSG00000198899")
cat("MT genes present in matrix:", sum(mt_ids %in% rownames(cnt)), "of", length(mt_ids),
    "| total reads across them:", sum(cnt[intersect(mt_ids, rownames(cnt)), ]), "\n")

y2 <- normLibSizes(y[a, , keep.lib.sizes = FALSE])
nf <- y2$samples$norm.factors
cat("\nTMM normalisation factors after filterByExpr, by group:\n")
print(round(do.call(rbind, tapply(nf, ss$group, function(v) c(median = median(v), min = min(v), max = max(v)))), 3))
# -> HC factors sit well above 1 and Lung below 1: normalisation is absorbing a group-wide
#    composition difference (see the ribosomal share). Expect long cancer-vs-HC DE lists.

cat("\nunannotated Ensembl-75 ids (no symbol in org.Hs.eg.db):", sum(is.na(sym)), "of", nrow(cnt),
    "| among filterByExpr-kept:", sum(is.na(sym[a])), "of", sum(a), "\n")

# ---- 6. sex is not in the metadata but is (noisily) inferable from Y-linked genes -------------
ymark <- c(RPS4Y1 = "ENSG00000129824", DDX3Y = "ENSG00000067048", UTY = "ENSG00000183878")
yreads <- colSums(cnt[ymark, ])
cat("\n=== Y-linked reads (RPS4Y1 + DDX3Y + UTY) ===\n")
cat("XIST (ENSG00000229807) detected in", sum(cnt["ENSG00000229807", ] > 0), "samples -> too sparse to use\n")
print(table(group = ss$group, y_reads_gt0 = yreads > 0))

# ---- 7. reference point: the paper's population (283 samples, all groups) ---------------------
c283 <- counts[, setdiff(colnames(counts), c("VU258-CRC", "VU398-GBM"))]
s283 <- ss0 <- read.csv("GSE68086_sample_sheet.csv", row.names = 1, check.names = FALSE)[colnames(c283), ]
cat("\n=== 283-sample cohort (paper analysed 5,003 genes) ===\n")
cat("filterByExpr(group, 7 groups):", sum(filterByExpr(DGEList(c283), group = s283$group)), "\n")
cat(">=5 reads in >=50% of samples:", sum(rowSums(c283 >= 5) >= ceiling(0.5 * ncol(c283))), "\n")
cat(">=5 reads in >=10% of samples:", sum(rowSums(c283 >= 5) >= ceiling(0.1 * ncol(c283))), "\n")
