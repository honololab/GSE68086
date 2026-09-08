# GSE68086 — Beginner-friendly differential expression with DESeq2 and PCA (R)

**Goal.** Find which genes are expressed differently in two cancers (Lung, GBM) compared with healthy controls (HC), using DESeq2, and use PCA to see and check the structure of the data along the way.

**How to use this.** Run the code blocks in order in one R script. Each block says what it does, what to plot, and how to read the plot. Everything runs on a laptop in a few minutes.

---

## 0. Setup

```r
install.packages(c("BiocManager", "ggplot2", "pheatmap", "ggrepel"))
BiocManager::install(c("DESeq2", "apeglm", "ashr", "org.Hs.eg.db", "limma"))

library(DESeq2); library(ggplot2); library(pheatmap); library(ggrepel)
library(org.Hs.eg.db); library(limma)
```

Files needed in your working directory:
- `GSE68086_TEP_data_matrix.csv` — raw counts, 57,736 genes × 285 samples (from the Kaggle zip)
- `GSE68086_sample_sheet.csv` — one row per sample, already parsed (provided with this document). Its `source` column matches the count-matrix column names.

Quick facts about the data (so nothing surprises you later): 1 sample per patient; shallow sequencing (median ~2 M reads); 85 % of the matrix is zeros; a handful of platelet genes hold ~28 % of all reads; samples were processed in 6 batches and batch is partly mixed with disease.

---

## 1. Load counts and sample information

**What this does.** Reads both files, aligns them (same samples, same order), and keeps only the columns we need.

```r
counts <- as.matrix(read.csv("GSE68086_TEP_data_matrix.csv", row.names = 1, check.names = FALSE))
ss     <- read.csv("GSE68086_sample_sheet.csv", check.names = FALSE)

ss <- data.frame(sample = ss$source,
                 group  = ss$`cancer type`,
                 batch  = ss$batch,
                 row.names = ss$source)

stopifnot(setequal(rownames(ss), colnames(counts)))   # must be TRUE
counts <- counts[, rownames(ss)]                        # same order
dim(counts); table(ss$group)
```

---

## 2. Choose the subset

**What this does.** Keeps Lung, GBM and HC, and only batches 02–04, where all three groups are present. This lets us adjust for batch instead of having it confused with disease.

```r
keep_s <- ss$group %in% c("Lung", "GBM", "HC") & ss$batch %in% c("Batch02", "Batch03", "Batch04")
ss     <- droplevels(ss[keep_s, ])
counts <- counts[, rownames(ss)]

ss$group <- relevel(factor(ss$group), ref = "HC")      # HC = reference level
ss$batch <- factor(ss$batch)
table(ss$group, ss$batch)
```
Expected: GBM 38, HC 55, Lung 40 (133 samples).

---

## 3. First look at the raw data (plots)

**What this does.** Before any modelling, check how deep each sample was sequenced and how many genes it detects.

```r
ss$lib_size  <- colSums(counts)          # total reads per sample
ss$n_detected <- colSums(counts > 0)     # genes with at least 1 read

ggplot(ss, aes(group, lib_size / 1e6, fill = group)) + geom_boxplot() +
  labs(y = "Library size (millions of reads)", title = "Sequencing depth per group")

ggplot(ss, aes(group, n_detected, fill = group)) + geom_boxplot() +
  labs(y = "Genes detected (count > 0)", title = "Library complexity per group")

# Which genes dominate? (composition)
top10 <- head(sort(rowSums(counts), decreasing = TRUE), 10)
round(100 * sum(top10) / sum(counts), 1)     # % of all reads in 10 genes
barplot(colSums(counts[names(top10), ]) / colSums(counts), col = as.integer(ss$group),
        border = NA, main = "Fraction of reads in the top-10 genes, per sample")
```
*How to read it.* HC samples are shallower (median ~1 M vs ~2 M) but detect **more** genes. That is a real difference in library composition between HC and cancer samples — keep it in mind when interpreting HC contrasts. The top-10 fraction should be roughly similar across groups; if not, normalization (Step 5) matters even more.

---

## 4. Build the DESeq2 object and remove low-count genes

**What this does.** Puts counts + sample sheet + model into one object. The model `~ batch + group` means "explain each gene by batch and by disease group". Then drops genes with too few reads to test.

```r
dds <- DESeqDataSetFromMatrix(countData = counts, colData = ss, design = ~ batch + group)

# keep genes with >= 10 reads in at least 38 samples (38 = size of the smallest group)
keep_g <- rowSums(counts(dds) >= 10) >= 38
table(keep_g)
dds <- dds[keep_g, ]
nrow(dds)          # expect roughly 8,000–10,000 genes
```
*Why filter.* Genes with almost no reads cannot show a reliable difference; keeping them only makes the multiple-testing correction harsher and the plots noisier.

---

## 5. Normalization

**What this does.** DESeq2 computes one *size factor* per sample (median-of-ratios method) so that samples with more reads do not look like they have "more expression" everywhere.

```r
dds <- estimateSizeFactors(dds)
ss$size_factor <- sizeFactors(dds)

ggplot(ss, aes(lib_size / 1e6, size_factor, colour = group)) + geom_point() +
  labs(x = "Library size (M reads)", y = "DESeq2 size factor")

# Before vs after: distribution of log counts per sample
par(mfrow = c(1, 2))
boxplot(log2(counts(dds) + 1), las = 2, col = as.integer(ss$group), main = "Raw log2 counts", outline = FALSE, xaxt = "n")
boxplot(log2(counts(dds, normalized = TRUE) + 1), las = 2, col = as.integer(ss$group), main = "Normalized log2 counts", outline = FALSE, xaxt = "n")
par(mfrow = c(1, 1))
```
*How to read it.* Size factors should track library size closely (roughly a straight line). After normalization the per-sample boxes should line up at a similar level. If HC samples still sit visibly apart, that is the composition difference from Step 3 showing through.

---

## 6. Variance-stabilizing transform (needed for PCA and heatmaps)

**What this does.** Raw counts have variance that grows with the mean, which would let a few high-count genes dominate any PCA. `vst()` puts all genes on a comparable scale. Use it for **plots only** — the statistical test in Step 8 still uses raw counts.

```r
vsd <- vst(dds, blind = TRUE)     # blind = does not use the design; good for exploration
head(assay(vsd)[, 1:4])
```

---

## 7. PCA — seeing the structure

### 7a. Built-in PCA (quick)
```r
plotPCA(vsd, intgroup = "group", ntop = 500)
plotPCA(vsd, intgroup = "batch", ntop = 500)
```
*How to read it.* Each point is a sample. If groups separate along PC1, disease is the main source of variation; if batches separate, batch is. Here you should see HC on one side, Lung on the other, GBM in between, and some batch structure on PC2.

### 7b. PCA by hand (more control)
```r
top_var <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), 1000)   # 1000 most variable genes
pca     <- prcomp(t(assay(vsd)[top_var, ]), scale. = FALSE)
pct     <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

# scree plot: how much variance each PC carries
barplot(pct[1:10], names.arg = paste0("PC", 1:10), ylab = "% variance", main = "Scree plot")

pc <- data.frame(pca$x[, 1:4], ss)
ggplot(pc, aes(PC1, PC2, colour = group, shape = batch)) + geom_point(size = 3) +
  labs(x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)"))
ggplot(pc, aes(PC1, PC3, colour = group, shape = batch)) + geom_point(size = 3)
ggplot(pc, aes(PC1, PC2, colour = log10(lib_size))) + geom_point(size = 3) +
  labs(title = "Is PC1 just sequencing depth?")
```
*How to read it.* Scree: PC1 should carry ~35 %, then a sharp drop. Depth colouring: if PC1 changes smoothly with library size, depth is leaking into the main axis (expected mildly here, correlation ≈ −0.3).

### 7c. Which genes drive PC1?
```r
load1 <- sort(pca$rotation[, "PC1"])
head(load1, 15); tail(load1, 15)      # most negative / most positive loadings
```
Save these IDs — in Step 9 we check whether they are also differentially expressed.

### 7d. PCA after removing batch (for visualization only)
```r
mat_bc <- removeBatchEffect(assay(vsd), batch = ss$batch, design = model.matrix(~ group, ss))
pca_bc <- prcomp(t(mat_bc[top_var, ]))
pct_bc <- round(100 * pca_bc$sdev^2 / sum(pca_bc$sdev^2), 1)
ggplot(data.frame(pca_bc$x[, 1:2], ss), aes(PC1, PC2, colour = group, shape = batch)) +
  geom_point(size = 3) + labs(title = "PCA after batch removal (display only)",
  x = paste0("PC1 (", pct_bc[1], "%)"), y = paste0("PC2 (", pct_bc[2], "%)"))
```
*Never* feed batch-corrected values into DESeq2; the model already handles batch.

### 7e. Sample-to-sample distances
```r
d <- dist(t(assay(vsd)))
pheatmap(as.matrix(d), annotation_col = ss[, c("group", "batch")], show_rownames = FALSE,
         show_colnames = FALSE, main = "Sample distances (vst)")
```
*How to read it.* Blocks of similar colour = groups of similar samples. A single row/column that is far from everything = candidate outlier. Check whether blocks follow `group` or `batch`.

---

## 8. Differential expression with DESeq2

**What this does.** One call fits the model for every gene: estimates dispersion (how noisy each gene is), shares that information across genes, and tests each coefficient.

```r
dds <- DESeq(dds)
resultsNames(dds)
# "Intercept" "batch_Batch03_vs_Batch02" "batch_Batch04_vs_Batch02" "group_GBM_vs_HC" "group_Lung_vs_HC"

plotDispEsts(dds)    # black = per-gene, red = trend, blue = final (shrunk) estimates
```
*How to read it.* Dispersion should decrease with mean expression and the blue points should hug the red trend. Very high scattered black points at low counts are normal for shallow data.

### 8a. The three contrasts
```r
alpha <- 0.05
res_lung <- results(dds, name = "group_Lung_vs_HC", alpha = alpha)
res_gbm  <- results(dds, name = "group_GBM_vs_HC",  alpha = alpha)
res_lg   <- results(dds, contrast = c("group", "Lung", "GBM"), alpha = alpha)

summary(res_lung); summary(res_gbm); summary(res_lg)
```
`summary()` tells you how many genes go up/down at your FDR threshold, and how many were removed by outlier or low-count filtering.

### 8b. Shrink the fold changes (for ranking and plotting)
Raw log2 fold changes for low-count genes are noisy and often huge. Shrinkage pulls unreliable ones toward zero without changing the p-values.
```r
shr_lung <- lfcShrink(dds, coef = "group_Lung_vs_HC", type = "apeglm")
shr_gbm  <- lfcShrink(dds, coef = "group_GBM_vs_HC",  type = "apeglm")
shr_lg   <- lfcShrink(dds, contrast = c("group", "Lung", "GBM"), res = res_lg, type = "ashr")
```

### 8c. Add gene symbols
```r
add_symbols <- function(res) {
  res$symbol <- mapIds(org.Hs.eg.db, keys = rownames(res), column = "SYMBOL",
                       keytype = "ENSEMBL", multiVals = "first")
  res
}
shr_lung <- add_symbols(shr_lung); shr_gbm <- add_symbols(shr_gbm); shr_lg <- add_symbols(shr_lg)
```

---

## 9. Plotting the DE results

### 9a. p-value histogram (the first thing to check)
```r
par(mfrow = c(1, 3))
hist(res_lung$pvalue, breaks = 50, main = "Lung vs HC", xlab = "p-value")
hist(res_gbm$pvalue,  breaks = 50, main = "GBM vs HC",  xlab = "p-value")
hist(res_lg$pvalue,   breaks = 50, main = "Lung vs GBM", xlab = "p-value")
par(mfrow = c(1, 1))
```
*How to read it.* A healthy histogram is flat with a spike near 0 (the spike = real signal). A U-shape or a hump in the middle means something is wrong with the model (e.g., a missing covariate).

### 9b. MA plots
```r
par(mfrow = c(1, 3))
plotMA(shr_lung, ylim = c(-4, 4), main = "Lung vs HC")
plotMA(shr_gbm,  ylim = c(-4, 4), main = "GBM vs HC")
plotMA(shr_lg,   ylim = c(-4, 4), main = "Lung vs GBM")
par(mfrow = c(1, 1))
```
*How to read it.* x = average expression, y = log2 fold change, blue = significant. Significant genes should appear across the expression range, not only at the far right (that would point to a composition artefact).

### 9c. Volcano plots
```r
volcano <- function(res, title, n_label = 15) {
  df <- as.data.frame(res); df <- df[!is.na(df$padj), ]
  df$sig <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1, "DE", "not DE")
  lab <- head(df[order(df$padj), ], n_label)
  ggplot(df, aes(log2FoldChange, -log10(padj), colour = sig)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_text_repel(data = lab, aes(label = symbol), size = 3, colour = "black") +
    scale_colour_manual(values = c("DE" = "firebrick", "not DE" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = 2) + geom_hline(yintercept = -log10(0.05), linetype = 2) +
    labs(title = title, x = "log2 fold change (shrunken)", y = "-log10 adjusted p")
}
volcano(shr_lung, "Lung vs HC"); volcano(shr_gbm, "GBM vs HC"); volcano(shr_lg, "Lung vs GBM")
```

### 9d. Heatmap of the top genes
```r
top_genes <- head(rownames(shr_lung)[order(shr_lung$padj)], 40)
mat <- mat_bc[top_genes, ]                       # batch-removed vst values from Step 7d
rownames(mat) <- shr_lung[top_genes, "symbol"]
pheatmap(mat, scale = "row", annotation_col = ss[, c("group", "batch")],
         show_colnames = FALSE, main = "Top 40 genes, Lung vs HC (row-scaled)")
```
*How to read it.* Rows = genes, columns = samples, colour = expression relative to that gene's average. A good result shows a clear colour split that lines up with the `group` bar, not with the `batch` bar.

### 9e. Individual genes
```r
for (g in head(top_genes, 4))
  plotCounts(dds, gene = g, intgroup = "group", main = shr_lung[g, "symbol"])
```
*How to read it.* Normalized counts per sample. You want the groups to be separated *and* the points within a group to be reasonably tight; one or two extreme dots driving the difference is a warning sign.

### 9f. Do PCA and DE agree?
```r
pc1_genes <- names(c(head(load1, 100), tail(load1, 100)))            # top PC1 drivers
de_lung   <- rownames(shr_lung)[which(shr_lung$padj < 0.05)]
length(intersect(pc1_genes, de_lung)) / length(pc1_genes)            # fraction of PC1 drivers that are DE

# PCA using only DE genes
pca_de <- prcomp(t(mat_bc[de_lung, ]))
ggplot(data.frame(pca_de$x[, 1:2], ss), aes(PC1, PC2, colour = group)) + geom_point(size = 3) +
  labs(title = "PCA on Lung-vs-HC DE genes only")
```
*How to read it.* If most PC1 drivers are DE genes, the dominant axis of variation *is* the disease signal. The DE-only PCA should separate Lung from HC cleanly (that is partly circular, so treat it as a sanity check, not evidence).

---

## 10. Export the results

```r
write_res <- function(res, file) {
  df <- as.data.frame(res)[order(res$padj), ]
  df <- cbind(gene_id = rownames(df), df)
  write.csv(df, file, row.names = FALSE)
}
write_res(shr_lung, "DE_Lung_vs_HC.csv")
write_res(shr_gbm,  "DE_GBM_vs_HC.csv")
write_res(shr_lg,   "DE_Lung_vs_GBM.csv")

# genes significant in both cancers vs HC, same direction
both <- intersect(rownames(shr_lung)[which(shr_lung$padj < 0.05)],
                  rownames(shr_gbm)[which(shr_gbm$padj  < 0.05)])
same_dir <- both[sign(shr_lung[both, "log2FoldChange"]) == sign(shr_gbm[both, "log2FoldChange"])]
length(both); length(same_dir)
```

**What to report.** For each contrast: number of DE genes (padj < 0.05), top 20 by shrunken fold change with symbol, the p-value histogram, one volcano, one heatmap, and the PCA coloured by group. Plus one paragraph on the caveats below.

---

## 11. Caveats to keep in mind

- **HC differs in library complexity**, not only depth. Some "cancer vs HC" genes may reflect sample handling rather than biology. The **Lung vs GBM** contrast does not have this problem and is the cleaner disease-specific list.
- **Batch** is included additively in the model; it cannot be fully separated from group with these sample sizes, so a few batch-driven genes may remain.
- **Shallow counts** → use shrunken fold changes for ranking, never raw ones.
- **PCA is descriptive.** It shows structure; DESeq2 tests it. Use PCA to check assumptions and to sanity-check the DE list, not as evidence on its own.

## 12. Natural next steps (when the above runs cleanly)

1. A global 3-group test: `DESeq(dds, test = "LRT", reduced = ~ batch)` — "does group matter at all for this gene?"
2. Infer sex from `RPS4Y1` (ENSG00000129824) on the vst matrix, add it to `ss`, refit with `~ batch + sex + group`, compare lists.
3. Gene-set enrichment on the ranked shrunken LFC (`fgsea` or `clusterProfiler`).
4. Repeat on CRC instead of GBM to see how stable the HC-related genes are.
5. Later: the edgeR + comparison pipeline from the previous document, run on the same subset.
