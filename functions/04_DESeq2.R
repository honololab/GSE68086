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
# keep genes with >= 10 reads in at least 38 samples (38 = size of the smallest group)
keep_g <- rowSums(counts(dds) >= 10) >= 38 # a gene survives if it is reasonably expressed in at least one whole group's worth of samples. 57,736 → 6,109
print(table(keep_g))
dds <- dds[keep_g, ]
print(nrow(dds))         # 6,109 genes

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