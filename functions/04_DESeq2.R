# 4. Build the DESeq2 object and remove low-count genes (Normalization)
library(DESeq2)
library(ggplot2)

d <- readRDS("data/derived/subset_lung_gbm_hc.rds")     # output of 02
counts <- d$counts; ss <- d$ss
ss$lib_size <- colSums(counts)   # 03 only has it in-session, so recompute here


# A DESeqDataSet bundles three things that until now were floating separately: 
# the count matrix, the sample sheet, and the model formula. From here on they travel together,
# which is what prevents the labels from detaching from the columns.
dds <- DESeqDataSetFromMatrix(countData = counts, colData = ss, design = ~ batch + group)

# keep genes with >= 10 reads in at least 38 samples (38 = size of the smallest group)
keep_g <- rowSums(counts(dds) >= 10) >= 38
table(keep_g)
dds <- dds[keep_g, ]
nrow(dds)          # 6,109 genes
## 5. Normalization
dds <- estimateSizeFactors(dds)
ss$size_factor <- sizeFactors(dds)

print(ggplot(ss, aes(lib_size / 1e6, size_factor, colour = group)) + geom_point() +
  labs(x = "Library size (M reads)", y = "DESeq2 size factor"))

# Before vs after: distribution of log counts per sample
par(mfrow = c(1, 2))
boxplot(log2(counts(dds) + 1), las = 2, col = as.integer(ss$group), main = "Raw log2 counts", outline = FALSE, xaxt = "n")
boxplot(log2(counts(dds, normalized = TRUE) + 1), las = 2, col = as.integer(ss$group), main = "Normalized log2 counts", outline = FALSE, xaxt = "n")
par(mfrow = c(1, 1))
