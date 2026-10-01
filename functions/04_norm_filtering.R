library(DESeq2)
library(edgeR)
# 04 -- Build the DESeq2 object, remove low-count genes, normalise with TMM, save
#
# Takes the 133-sample subset that 02 saved and prepares it for everything that follows.
# Four steps, in this order:
#   1. counts + sample sheet + model formula (~ batch + group) go into one object
#   2. genes with too few reads are removed (edgeR::filterByExpr)
#   3. one size factor per sample is estimated with TMM (edgeR::calcNormFactors)
#   4. counts are put on a log-like scale (vst), for plots only
# Then the two objects are saved: 05 reads them for the PCA, 07 for the DE tests.
# No differential-expression test runs here.

IN      <- "data/derived/subset_lung_gbm_hc.rds"     # written by 02
OUT_DDS <- "data/derived/dds_filtered.rds"           # 07 starts from this
OUT_VSD <- "data/derived/vsd_blind.rds"              # 05 starts from this

stopifnot("subset not found -- run 02 first" = file.exists(IN))

# ---- 0. load the subset from 02. gives counts and ss.
d <- readRDS(IN)
counts <- d$counts
ss <- d$ss

# 02 checked these before saving. Checking again means a stale .rds fails here.
stopifnot(identical(colnames(counts), rownames(ss)))   # samples in the same order
stopifnot(levels(ss$group)[1] == "HC")                 # HC is the reference level
print(dim(counts))                                     # 57736 x 133
print(table(ss$group, ss$batch))

# ---- 1. one object: counts + sample sheet + model

# ~ batch + group: "for every gene, explain its counts by which batch the sample was
# processed in and which group it belongs to". Batch is in the formula because depth
# follows the batch (03, figure 1) and the controls sit mostly in the shallow batch.
# HC is the first level of group, so every group coefficient reads "<group> vs HC".
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = ss,
                              design    = ~ batch + group)

# A DESeqDataSet bundles three things that until now were floating separately:
# the count matrix, the sample sheet, and the model formula. From here on they travel together,
# which is what prevents the labels from detaching from the columns.

# ---- 2. remove genes with too few reads to test

# filterByExpr keeps a gene if it reaches a cutoff in enough samples. The cutoff is
# in CPM (about 10 reads in the median library, ~5.9 CPM here), so a shallow HC
# sample and a deep Lung sample are asked the same question relative to their depth.
# "Enough samples" comes from the smallest group (38), reduced by edgeR's tolerance
# to about 30. The gene must also reach 15 reads in total.
n_any <- sum(rowSums(counts) > 0) # genes with any read at all
keep_g <- filterByExpr(counts(dds), group = dds$group)

cat(sprintf("\ngenes: %d in the matrix, %d with at least one read, %d kept by filterByExpr\n",
            nrow(counts), n_any, sum(keep_g)))
    # genes: 57736 in the matrix, 26177 with at least one read, 6994 kept by filterByExpr

kept_share <- colSums(counts[keep_g, ]) / colSums(counts)
cat("share of each sample's reads kept by the filter, median by group:\n")
print(round(tapply(kept_share, ss$group, median), 4))

dds <- dds[keep_g, ]

# counts(dds) extracts the raw count matrix (integers, pre-normalization).
# dds$group pulls the group column from colData(dds), since $ on a
# DESeqDataSet forwards to colData. It's shorthand for colData(dds)$group.
# The group argument is what makes the filter aware of your design:
# without it, the function would require a gene to be expressed in a
# fixed fraction of all samples, which would discard genes expressed only
# in the smallest condition.

# ---- 3. size factors with TMM: one number per sample

# DESeq2 divides each sample's counts by its size factor so that a sample with more
# reads does not look like it has more expression everywhere. TMM computes that
# number in two parts: the library size (total reads on the kept genes) times a
# correction for composition (norm.factors), so a few very abundant genes in one
# sample do not make all its other genes look lower. The product is the effective
# library size. Dividing by its geometric mean puts it on the scale DESeq2 uses
# (size factors around 1). Same filter + TMM as 04c, so the factors match.
y <- DGEList(counts(dds), group = dds$group)    # dds is already filtered, so lib.size
y <- calcNormFactors(y, method = "TMM")         # is the total on the kept genes
eff <- y$samples$lib.size * y$samples$norm.factors
sizeFactors(dds) <- eff / exp(mean(log(eff)))

cat("\nTMM factors, median by group:\n")
print(round(tapply(y$samples$norm.factors, ss$group, median), 3))
sf <- sizeFactors(dds)
cat("size factors, median by group:\n")
print(round(tapply(sf, ss$group, median), 3))

# ---- 4. vst: a log-like scale for plots

# Raw counts spread more the higher a gene's mean, so a PCA on them would be a PCA of
# the few most abundant genes. vst() removes that dependence and applies the size
# factors set above (TMM) on the way. blind = TRUE: it does not look at the design,
# which is right when the question is what structure the data has on its own. For
# plots only -- 07 tests the raw counts in dds, never these values.
vsd <- vst(dds, blind = TRUE)

# ---- 5. save

saveRDS(dds, OUT_DDS)
saveRDS(vsd, OUT_VSD)
cat(sprintf("\nwritten %s (%d genes x %d samples) and %s\n",
            OUT_DDS, nrow(dds), ncol(dds), OUT_VSD))

# ---- 6. filter summary

cat("\n=== resumen de filtros ===\n")
cat(sprintf("muestras: 285 en GEO -> %d en el subconjunto Lung/GBM/HC, Batch02-04\n", ncol(counts)))
cat(sprintf("genes:    %d en la matriz -> %d con alguna lectura -> %d tras filterByExpr\n",
            nrow(counts), n_any, sum(keep_g)))
