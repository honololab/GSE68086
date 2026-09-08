# GSE68086 — Exploratory differential-expression pipeline (edgeR / DESeq2)

**Scope.** Methodological sandbox, unrelated to the hantavirus cohort. Goal: find which genes are most associated with two chosen cancers vs. healthy controls, using edgeR and DESeq2 side by side, without going deep into methodology yet.

---

## 0. What the data looks like (from a first inspection)

| Item | Observation |
|---|---|
| Matrix | 57,736 Ensembl genes × 285 samples, raw integer HTSeq counts, no missing values |
| Design | 1 sample per patient, single time point → **cross-sectional**, no longitudinal structure |
| Classes | Lung 60 · HC (healthy) 55 · CRC 42 · GBM 40 · Breast 39 · Pancreas 35 · Hepatobiliary 14 |
| Metadata | `tissue`, `cell type`, `patient id`, `cancer type`, `batch` (6 batches), `mutational subclass` (wt / KRAS / EGFR / HER2+ / PIK3CA / MET / TN) |
| Depth | Median 2.0 M reads/sample, range 0.22–7.3 M → **shallow**, ~10× spread |
| Sparsity | 85 % zeros; ~28.5 k genes ever detected; **~8.9 k genes with ≥10 counts in ≥20 samples** |
| Composition | Top 10 genes carry **28.5 %** of all counts (platelet transcripts — IDs point to PPBP, TMSB4X, B2M, ACTB, F13A1, CLU, SPARC, FTL; confirm with annotation) |
| Depth by class | HC median 1.04 M vs Lung 2.04 M vs GBM 1.73 M, yet HC detects the most genes (9.7 k vs 7.7 k in Lung) → HC libraries differ in complexity, not just size |
| Batch | Strongly tied to class: HC only in Batch02–04; Batch05/06 contain only Breast and Lung |
| Quick PCA (log2 CPM, Lung/GBM/HC) | PC1 (34 %) separates HC from Lung with GBM in between; PC2 and PC3 are **batch-driven** (R² ≈ 0.6). Restricting to Batch02–04: PC1 37 %, class-driven; batch influence drops (PC2 R² 0.42). No outlier samples |
| Sex | Not recorded; RPS4Y1 detected in ~half of samples → sex can be inferred and used as covariate |

Two practical traps found in the files:

1. **Column names ≠ patient id.** Count-matrix columns (e.g. `3-Breast-Her2-ampl`, `MGH-NSCLC-L23-TR524`) match `Sample_source_name_ch1` in the series matrix, *not* `Sample_title` or `patient id` (only 58/285 match those). Map through `source_name`.
2. **Shifted characteristics fields.** In the series matrix the `characteristics_ch1` lines are not aligned by position across samples. Parse each cell by its `key: value` prefix (GEOquery's `*:ch1` columns do this correctly).

---

## 1. Choice of groups

**Default: Lung (NSCLC) vs GBM vs HC, restricted to Batch02–04.**

| | Batch02 | Batch03 | Batch04 | Total |
|---|---|---|---|---|
| GBM | 12 | 15 | 11 | 38 |
| HC | 16 | 15 | 24 | 55 |
| Lung | 19 | 7 | 14 | 40 |
| | | | | **133** |

Why: (i) two of the largest groups; (ii) biologically distant tumours (lung vs. brain) so a contrast between them is meaningful; (iii) all three classes are present in each of these batches, so **batch can be adjusted rather than being confounded**. Lung samples in Batch05/06 (20 of 60) are dropped because no HC/GBM exist there — they would only inform the batch coefficient, not the group contrast.

Alternative: CRC (42) instead of GBM; same batch restriction gives CRC 38.

Contrasts of interest: `Lung − HC`, `GBM − HC`, `Lung − GBM`, plus a global 3-group test.

---

## 2. Pipeline

Each step: **purpose → what to do → what to look at**. Code is illustrative, not final.

### Step 1 · Ingest and build a clean sample sheet
**Purpose.** A count matrix and a sample sheet whose rows correspond exactly, with the covariates we will use (`group`, `batch`, `subclass`).

```r
library(GEOquery); library(edgeR)
counts <- as.matrix(read.csv("GSE68086_TEP_data_matrix.csv", row.names = 1, check.names = FALSE))
gse    <- getGEO("GSE68086", GSEMatrix = TRUE)[[1]]
pd     <- pData(gse)
ss <- data.frame(sample = pd$source_name_ch1,          # == colnames(counts)
                 group  = pd$`cancer type:ch1`,
                 batch  = pd$`batch:ch1`,
                 subclass = pd$`mutational subclass:ch1`)
stopifnot(setequal(ss$sample, colnames(counts)))
counts <- counts[, ss$sample]
```
**Look at.** `table(ss$group, ss$batch)`; any duplicated patient ids; that the top genes by total count are platelet genes (sanity check that the matrix is what it claims).

### Step 2 · Subset and define the design
**Purpose.** Fix the analysis population and the model before looking at results (avoids drifting choices later).

```r
keep_s <- ss$group %in% c("Lung","GBM","HC") & ss$batch %in% c("Batch02","Batch03","Batch04")
ss <- droplevels(ss[keep_s, ]); counts <- counts[, ss$sample]
ss$group <- relevel(factor(ss$group), ref = "HC")
design <- model.matrix(~ batch + group, data = ss)
```
**Look at.** The final crosstab (should match the table above). Decide now whether `batch` is additive only (yes for this exploration).

### Step 3 · QC and exploratory description
**Purpose.** Understand what drives variation *before* modelling: depth, complexity, batch, class, outliers.

```r
y <- DGEList(counts, samples = ss, group = ss$group)
lib   <- colSums(counts); det <- colSums(counts > 0)
boxplot(log10(lib) ~ ss$group); boxplot(det ~ ss$group)        # depth & complexity by class
logcpm <- cpm(y, log = TRUE, prior.count = 2)
plotMDS(logcpm, col = as.integer(ss$group), pch = as.integer(factor(ss$batch)))
cor_s <- cor(logcpm, method = "spearman"); pheatmap::pheatmap(cor_s, annotation_col = ss[, c("group","batch")])
```
**Look at.**
- Do HC samples cluster apart already on MDS dim 1 (expected from the PCA above)? How much of dims 1–2 is batch?
- Sample–sample correlation heatmap: any sample that correlates poorly with everything (outlier)? Any HC/cancer block structure that follows batch rather than class?
- Share of counts in the top 10–20 genes per sample: if it varies systematically by class, composition normalization matters (TMM / size factors will handle it, but note it).
- Sex inference: `logcpm["ENSG00000129824", ]` (RPS4Y1) bimodal? If yes, add `sex` to the sample sheet.

### Step 4 · Gene annotation
**Purpose.** Ensembl IDs → symbols and biotypes, so results are readable and so we can decide whether to keep only protein-coding genes.

```r
library(EnsDb.Hsapiens.v86)   # or biomaRt
ann <- AnnotationDbi::select(EnsDb.Hsapiens.v86, keys = rownames(counts),
                             columns = c("SYMBOL","GENEBIOTYPE"), keytype = "GENEID")
```
**Look at.** Biotype breakdown of the ~8.9 k expressed genes. Platelet RNA has many non-coding / mitochondrial / ribosomal transcripts; keep everything for now but flag biotype in the result tables.

### Step 5 · Filtering and normalization
**Purpose.** Remove genes with too few counts to test (they only inflate the multiple-testing burden and destabilize dispersion) and correct for depth and composition.

```r
# edgeR
keep_g <- filterByExpr(y, group = ss$group)     # expect roughly 8–10 k genes
y <- y[keep_g, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")         # normLibSizes() in edgeR ≥ 4
# DESeq2
library(DESeq2)
dds <- DESeqDataSetFromMatrix(counts[keep_g, ], colData = ss, design = ~ batch + group)
dds <- estimateSizeFactors(dds)
```
**Look at.** `summary(keep_g)`; distribution of TMM factors and DESeq2 size factors by class — if HC systematically get factors far from 1, that confirms the composition/complexity difference seen in QC. `plotMDS` again after normalization.

### Step 6 · Differential expression with edgeR (quasi-likelihood)
**Purpose.** First set of DE genes with a method that is robust for moderate n and handles the batch covariate in a GLM.

```r
y  <- estimateDisp(y, design, robust = TRUE); plotBCV(y)
fit <- glmQLFit(y, design, robust = TRUE);     plotQLDisp(fit)
con <- makeContrasts(Lung_vs_HC = groupLung, GBM_vs_HC = groupGBM,
                     Lung_vs_GBM = groupLung - groupGBM, levels = design)
res_e <- lapply(colnames(con), function(k) topTags(glmQLFTest(fit, contrast = con[, k]), n = Inf)$table)
names(res_e) <- colnames(con)
anova_e <- topTags(glmQLFTest(fit, coef = grep("^group", colnames(design))), n = Inf)$table  # global test
```
**Look at.** Common BCV (platelet data tends to be noisy, > 0.4 would not be surprising); number of genes at FDR < 0.05 per contrast; `plotMD` per contrast; whether DE genes are dominated by very-high-count platelet genes (composition artefact) or spread across the expression range.

### Step 7 · Differential expression with DESeq2
**Purpose.** Same questions with an independent estimation framework (median-of-ratios normalization, empirical-Bayes shrinkage of dispersion and LFC). Agreement between the two strengthens the gene list; disagreement is informative.

```r
dds <- DESeq(dds)
res_d <- list(
  Lung_vs_HC  = lfcShrink(dds, coef = "group_Lung_vs_HC", type = "apeglm"),
  GBM_vs_HC   = lfcShrink(dds, coef = "group_GBM_vs_HC",  type = "apeglm"),
  Lung_vs_GBM = results(dds, contrast = c("group","Lung","GBM")))
anova_d <- results(DESeq(dds, test = "LRT", reduced = ~ batch))     # global test
```
**Look at.** `plotDispEsts(dds)`; `plotMA` per contrast; count of DE genes at padj < 0.05; how many genes were flagged as outliers by Cook's distance (many flags → look at those samples).

### Step 8 · Compare the two tools
**Purpose.** This is the methodological core of the exercise: how much of the "signal" depends on the tool.

```r
cmp <- merge(res_e$Lung_vs_HC["logFC"], as.data.frame(res_d$Lung_vs_HC)["log2FoldChange"], by = 0)
cor(cmp$logFC, cmp$log2FoldChange, method = "spearman")
de_e <- rownames(res_e$Lung_vs_HC)[res_e$Lung_vs_HC$FDR < 0.05]
de_d <- rownames(res_d$Lung_vs_HC)[which(res_d$Lung_vs_HC$padj < 0.05)]
length(intersect(de_e, de_d)); length(setdiff(de_e, de_d)); length(setdiff(de_d, de_e))
```
**Look at.** LFC rank correlation (expect high, > 0.9); overlap of significant sets; genes significant in one tool only — are they low-count, high-dispersion, or driven by a few samples? Repeat for each contrast.

### Step 9 · Interpret and report
**Purpose.** Turn tables into a short, defensible description of "genes most related to each disease".

For each contrast, report the genes that are significant **in both tools**, ranked by shrunken LFC, with symbol, biotype, mean expression, LFC (both tools), FDR (both tools). Add:
- Volcano plot per contrast; heatmap of the top ~50 consensus genes (batch-corrected `logcpm` via `limma::removeBatchEffect`, **for display only**, never for testing).
- A light enrichment pass (`clusterProfiler::enrichGO` or `fgsea` on the ranked LFC) to see whether the consensus genes point to coherent biology (e.g., platelet activation, immune/interferon signatures, splicing) or look like noise.
- A specific note on the **Lung vs GBM** contrast: genes here are not confounded by the HC processing difference (see caveats), so they are the more trustworthy "disease-specific" candidates.

### Step 10 · Caveats to carry into the write-up
- **HC is not a clean control.** HC libraries are shallower but more complex, and HC only exists in three batches. Some "cancer vs HC" genes may reflect sample handling rather than tumour biology. The Lung vs GBM contrast partly sidesteps this.
- **Composition.** A handful of platelet genes dominate; normalization handles this but check that DE lists are not just the complement of those genes.
- **Shallow, sparse counts.** Low-count genes will have unstable LFC; prefer shrunken LFC when ranking.
- **Batch.** Modelled additively; interactions (batch × group) are not identifiable with these sizes and are not attempted.
- **No replication of patients.** One sample per patient, so no pseudoreplication issue — but also no within-patient information.

---

## 3. Suggestions and further steps

1. **Add limma-voom as a third method.** Cheap to run once the `DGEList` exists, and it is the linear-model machinery closest to what a later penalized / classification pipeline would use.
2. **Infer sex from RPS4Y1/XIST and include it** in the design; check whether it changes the top lists.
3. **Sub-questions inside Lung.** `mutational subclass` splits the 40 Lung samples into wt 14 · KRAS 12 · EGFR/MET 13 (23 KRAS in the full Lung set) — a second exploratory contrast (e.g., KRAS vs wt) once the main pipeline runs; sizes are small, so treat it as hypothesis-generating.
4. **Stability of the gene list.** Repeat Steps 5–8 on bootstrap or subsampled cohorts (e.g., 80 % of patients, 50 times) and report how often each gene stays significant. This is a small-n habit worth building now.
5. **Bridge to the hantavirus project.**
   - Collapse genes to **pathway / gene-set scores** (`GSVA`, ssGSEA) and re-run the same contrasts — tests the dimensionality-reduction step of the house pipeline.
   - Fit a **penalized multinomial model** (glmnet) on the 133 samples with **DE-based filtering done inside the CV folds**, not on the full data — a clean check against double dipping.
   - Try a **kernel independence test (HSIC)** between expression and group as a filter-free screening step, and compare its top genes with the edgeR/DESeq2 consensus.
6. **Hold-out check.** The original study used training/validation splits; setting aside ~30 % of patients before Step 5 and checking that the consensus genes separate groups in the held-out set is a low-cost credibility check.
7. **Later, if useful:** protein-coding-only analysis; explicit treatment of mitochondrial/ribosomal fraction as a QC covariate; comparison with the gene lists published for this dataset (PMID 26525104).

---

## Appendix · Setup

```r
install.packages(c("BiocManager","pheatmap"))
BiocManager::install(c("edgeR","DESeq2","limma","GEOquery","apeglm",
                       "EnsDb.Hsapiens.v86","clusterProfiler","fgsea","GSVA"))
```
Data: Kaggle mirror of GEO GSE68086 (`GSE68086_TEP_data_matrix.csv` + `GSE68086_series_matrix.txt`), ~70 MB unzipped. Downloadable without login via `https://www.kaggle.com/api/v1/datasets/download/samiraalipour/gene-expression-omnibus-geo-dataset-gse68086`.
