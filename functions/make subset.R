# libraries setup (uncomment only if missing)
# renv::install(c("ggplot2", "pheatmap", "ggrepel",
#                 "bioc::DESeq2", "bioc::apeglm", "bioc::org.Hs.eg.db",
#                 "bioc::edgeR"))
# install.packages("ashr")   # CRAN, not Bioconductor

library(DESeq2); library(ggplot2); library(pheatmap); library(ggrepel)
library(org.Hs.eg.db); library(limma)

# load counts and sample information
counts <- as.matrix(read.csv("GSE68086_TEP_data_matrix.csv", row.names = 1, check.names = FALSE))
ss     <- read.csv("GSE68086_sample_sheet_claude.csv", check.names = FALSE)

ss <- data.frame(sample = ss$source,
                 group  = ss$`cancer type`,
                 batch  = ss$batch,
                 row.names = ss$source)

stopifnot(setequal(rownames(ss), colnames(counts)))   # must be TRUE
counts <- counts[, rownames(ss)]                        # same order
dim(counts); table(ss$group)

# choose the substet
keep_s <- ss$group %in% c("Lung", "GBM", "HC") & ss$batch %in% c("Batch02", "Batch03", "Batch04")
ss     <- droplevels(ss[keep_s, ])
counts <- counts[, rownames(ss)]

ss$group <- relevel(factor(ss$group), ref = "HC")      # HC = reference level
ss$batch <- factor(ss$batch)
table(ss$group, ss$batch)