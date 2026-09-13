# 4. Build the DESeq2 object and remove low-count genes (Normalization)
library(DESeq2)
library(ggplot2)
library(edgeR)

# reload the subset from 02. gives counts and ss.
d <- readRDS("data/derived/subset_lung_gbm_hc.rds")     # output of 02
counts <- d$counts; ss <- d$ss
ss$lib_size <- colSums(counts)   # 03 only has it in-session, so recompute here


# A DESeqDataSet bundles three things that until now were floating separately: 
# the count matrix, the sample sheet, and the model formula. From here on they travel together,
# which is what prevents the labels from detaching from the columns.
# "for every gene, explain its counts by which batch the sample was in and which group it belongs to."
dds <- DESeqDataSetFromMatrix(countData = counts, colData = ss, design = ~ batch + group)

# filter genes
# keep genes with enough counts (CPM scale, group-aware) -- edgeR::filterByExpr
keep_g <- filterByExpr(counts(dds), group = dds$group)
print(table(keep_g))
dds <- dds[keep_g, ]
print(nrow(dds))      # 6,994

## 5. Normalization
dds <- estimateSizeFactors(dds) # library size factor 
ss$size_factor <- sizeFactors(dds)

print(ggplot(ss, aes(lib_size / 1e6, size_factor, colour = group)) + geom_point() +
  labs(x = "Library size (M reads)", y = "DESeq2 size factor"))

# Before vs after: distribution of log counts per sample
par(mfrow = c(1, 2))
boxplot(log2(counts(dds) + 1), las = 2, col = as.integer(ss$group), main = "Raw log2 counts", outline = FALSE, xaxt = "n")
boxplot(log2(counts(dds, normalized = TRUE) + 1), las = 2, col = as.integer(ss$group), main = "Normalized log2 counts", outline = FALSE, xaxt = "n")
par(mfrow = c(1, 1))

# export data
saveRDS(dds, "data/derived/dds_filtered.rds")

## 6. Variance-stabilizing transform (needed for PCA and heatmaps)
vsd <- vst(dds, blind = TRUE)     # blind = does not use the design; good for exploration
head(assay(vsd)[, 1:4])

# saveRDS(vsd, "data/derived/vsd_blind.rds")

dds