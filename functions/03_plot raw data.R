#  3. plotting raw data

library(ggplot2)   # 02 is base R only; the plots live here

d <- readRDS("data/derived/subset_lung_gbm_hc.rds")
counts <- d$counts; ss <- d$ss
print(dim(counts))

ss$lib_size  <- colSums(counts)          # total reads per sample
ss$n_detected <- colSums(counts > 0)     # genes with at least 1 read

print(ggplot(ss, aes(group, lib_size / 1e6, fill = group)) + geom_boxplot() +
  labs(y = "Library size (millions of reads)", title = "Sequencing depth per group"))

print(ggplot(ss, aes(group, n_detected, fill = group)) + geom_boxplot() +
  labs(y = "Genes detected (count > 0)", title = "Library complexity per group"))

# Which genes dominate? (composition)
top10 <- head(sort(rowSums(counts), decreasing = TRUE), 10)
print(round(100 * sum(top10) / sum(counts), 1))     # % of all reads in 10 genes
 barplot(colSums(counts[names(top10), ]) / colSums(counts), col = as.integer(ss$group),
        border = NA, main = "Fraction of reads in the top-10 genes, per sample")