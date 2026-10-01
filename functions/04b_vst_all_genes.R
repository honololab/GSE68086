# 04b -- vst on ALL genes, no filterByExpr and no removal of all-zero genes
#
# Same steps as 04, minus the gene filter: counts -> DESeqDataSet -> size factors -> vst.
# Written only so 05 can draw the PCA without any gene filter and compare it with the
# filtered one. Nothing downstream (06, 07) should read this file.

library(DESeq2)

IN      <- "data/derived/subset_lung_gbm_hc.rds"       # written by 02
OUT_VSD <- "data/derived/vsd_blind_allgenes.rds"       # 05 can read this instead of vsd_blind.rds

stopifnot("subset not found -- run 02 first" = file.exists(IN))

d      <- readRDS(IN)
counts <- d$counts
ss     <- d$ss
stopifnot(identical(colnames(counts), rownames(ss)))

cat(sprintf("genes: %d in total, %d with zero reads in every sample (kept)\n",
            nrow(counts), sum(rowSums(counts) == 0)))
    # genes: 57736 in total, 31559 with zero reads in every sample (kept)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData   = ss,
                              design    = ~ batch + group)

# Median-of-ratios (NOT TMM, which 04 uses). It only uses genes with no zero in any sample; here there
# are 1387 of them, so the default works without switching to "poscounts".
dds <- estimateSizeFactors(dds)
cat(sprintf("genes with no zero in any sample (used for size factors): %d\n",
            sum(rowSums(counts == 0) == 0)))

# blind = TRUE, as in 04. The all-zero genes get the same value in every sample, so
# they have zero variance and add nothing to the PCA -- but they stay in the object.
vsd <- vst(dds, blind = TRUE)

saveRDS(vsd, OUT_VSD)
cat(sprintf("written %s (%d genes x %d samples)\n", OUT_VSD, nrow(vsd), ncol(vsd)))
