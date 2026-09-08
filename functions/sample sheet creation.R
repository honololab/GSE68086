library(GEOquery)

gse <- getGEO("GSE68086", GSEMatrix = TRUE)[[1]]
raw <- pData(gse)

# check this first if anything still looks wrong:
names(raw)

ss <- data.frame(
  sample   = raw$source_name_ch1,
  group    = raw$`cancer type:ch1`,
  batch    = raw$`batch:ch1`,
  subclass = raw$`mutational subclass:ch1`,
  patient  = raw$`patient id:ch1`,
  row.names = raw$source_name_ch1
)

dim(ss)            # expect 285 x 5
table(ss$group)     # expect HC 55, Lung 60, CRC 42, GBM 40, Breast 39, Pancreas 35, Hepatobiliary 14

counts <- as.matrix(read.csv("GSE68086_TEP_data_matrix.csv", row.names = 1, check.names = FALSE))
stopifnot(setequal(rownames(ss), colnames(counts)))
counts <- counts[, rownames(ss)]

write.csv(ss, "GSE68086_sample_sheet.csv", row.names = TRUE)