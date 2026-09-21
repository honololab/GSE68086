# 01 -- sample sheet creation
# One row per sample: group, batch, subclass, patient. No genes here.
#
# Metadata comes from the series matrix file already in the repo, NOT from a download:
# getGEO(filename=) parses it offline, so this runs without network and gives the same
# result every time. GEOquery splits the characteristics by KEY, which is why the 40
# second-wave samples whose keys arrive in a rotated order still land in the right
# columns -- see 00b_series_df_dictionary.ipynb for the evidence.

library(GEOquery)

SERIES <- "GSE68086_series_matrix.txt"
COUNTS <- "GSE68086_TEP_data_matrix.csv"
OUT    <- "GSE68086_sample_sheet.csv"

# Paths are relative to the project root, which is what TEP.Rproj sets the wd to.
stopifnot(
  "series matrix not found -- is the working directory the project root?" = file.exists(SERIES),
  "counts matrix not found -- is the working directory the project root?" = file.exists(COUNTS)
)

# The series matrix table is empty (!Sample_data_row_count is "0" for every sample --
# the counts live in a separate supplementary file), so this gives 0 features and a
# full pData(). GEOquery's warning about the empty table is expected.
gse <- getGEO(filename = SERIES)
raw <- pData(gse)

ss <- data.frame(
  sample   = raw$source_name_ch1,
  group    = raw$`cancer type:ch1`,
  batch    = raw$`batch:ch1`,
  subclass = raw$`mutational subclass:ch1`,   # NA where the group has no subclass
  patient  = raw$`patient id:ch1`,
  # note     = raw$description,                 # authors' QC flags -- recorded, not acted on
  stringsAsFactors = FALSE
)

# What the comments used to only claim. 285 = 39+42+40+55+14+60+35.
stopifnot(nrow(ss) == 285)
stopifnot(anyDuplicated(ss$sample) == 0)
stopifnot(identical(
  as.integer(table(ss$group)[c("Breast", "CRC", "GBM", "HC",
                               "Hepatobiliary", "Lung", "Pancreas")]),
  c(39L, 42L, 40L, 55L, 14L, 60L, 35L)
))

# Check the sheet against the counts matrix using the HEADER ONLY: the sample names
# are all we need here, so there is no reason to read 57,736 rows of counts.
counts_cols <- colnames(read.csv(COUNTS, nrows = 1, row.names = 1, check.names = FALSE))
stopifnot(setequal(ss$sample, counts_cols))

write.csv(ss, OUT, row.names = FALSE)
