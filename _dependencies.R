# Not meant to be run. This file exists so renv::snapshot() can see the
# packages used in the .md pipeline docs -- renv does not scan plain .md.
# Add a library() line here whenever you start using a new package.

library(GEOquery)       # functions/sample sheet creation.R
library(BiocManager)

library(DESeq2)         # GSE68086 DESeq2 / PCA pipeline
library(limma)
library(apeglm)
library(ashr)

library(AnnotationDbi)  # annotation
library(org.Hs.eg.db)

library(ggplot2)        # plotting
library(ggrepel)
library(pheatmap)

library(glmnet)          # lasso / elastic net
