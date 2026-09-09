# libraries setup (uncomment only if missing)
# renv::install(c("ggplot2", "pheatmap", "ggrepel",
#                 "bioc::DESeq2", "bioc::apeglm", "bioc::org.Hs.eg.db",
#                 "bioc::edgeR"))
# install.packages("ashr")   # CRAN, not Bioconductor

# base R only -- ggplot2 moved to 03 along with the plots

# not needed yet -- uncomment as you add those steps:
# library(DESeq2); library(pheatmap); library(ggrepel)
# library(org.Hs.eg.db); library(limma)
# =================================

# 1. load counts and sample information
counts <- as.matrix(read.csv("GSE68086_TEP_data_matrix.csv", row.names = 1, check.names = FALSE))
ss     <- read.csv("GSE68086_sample_sheet_claude.csv", check.names = FALSE)

ss <- data.frame(sample = ss$source,
                 group  = ss$`cancer type`,
                 batch  = ss$batch,
                 row.names = ss$source)

stopifnot(setequal(rownames(ss), colnames(counts)))   # must be TRUE
counts <- counts[, rownames(ss)]                        # same order
print(dim(counts)); print(table(ss$group))

# 2. choose the substet
keep_s <- ss$group %in% c("Lung", "GBM", "HC") & ss$batch %in% c("Batch02", "Batch03", "Batch04")
ss     <- droplevels(ss[keep_s, ])
counts <- counts[, rownames(ss)]

ss$group <- relevel(factor(ss$group), ref = "HC")      # HC = reference level
ss$batch <- factor(ss$batch)
print(table(ss$group, ss$batch))

dir.create("data/derived", recursive = TRUE, showWarnings = FALSE)
saveRDS(list(counts = counts, ss = ss), "data/derived/subset_lung_gbm_hc.rds")

