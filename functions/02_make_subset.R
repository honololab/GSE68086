# 02 -- make subset
# Cuts the 285 samples down to the 133 that support a clean comparison, then saves the
# counts and the sample sheet together, already aligned, so 03 and 04 just load them.
#
# Why these batches: every one of the 9 group x batch cells has to contain samples, or
# the model in 04 cannot tell a batch effect apart from a group effect. Batch02/03/04 is
# where Lung, GBM and HC all appear together -- GBM and HC have no samples in Batch05 or
# Batch06, and Batch01 holds only 2 samples in total. This matters because sequencing
# depth here differs a lot between batches.
#
# Why these groups: the batches rule out Breast, which has no samples in Batch02. They do
# not rule out CRC, which has 4, 21 and 13 samples across Batch02/03/04 and would give a
# working design too. Pancreas and Hepatobiliary have only 1 sample each in Batch02, too
# few to be useful. So picking Lung and GBM over CRC is a choice, not something the data
# forced -- add "CRC" to GROUPS below if you want it in.
#
# Base R only -- ggplot2 and the plots live in 03.

COUNTS <- "GSE68086_TEP_data_matrix.csv"
SHEET  <- "GSE68086_sample_sheet.csv"                      # written by 01
OUTDIR <- "data/derived"
OUT    <- file.path(OUTDIR, "subset_lung_gbm_hc.rds")      # 03 and 04 read this path

GROUPS  <- c("Lung", "GBM", "HC")
BATCHES <- c("Batch02", "Batch03", "Batch04")
REF     <- "HC"                                            # baseline for every group contrast

# Paths are relative to the project root, which is what TEP.Rproj sets the wd to.
stopifnot(
  "counts matrix not found -- is the working directory the project root?" = file.exists(COUNTS),
  "sample sheet not found -- run 01 first" = file.exists(SHEET)
)

# 1. load counts and sample information
counts <- as.matrix(read.csv(COUNTS, row.names = 1, check.names = FALSE))
ss     <- read.csv(SHEET, check.names = FALSE)

# 01 writes these three columns. Check them here: if one were renamed, ss$<name> would
# be NULL, data.frame() would quietly drop it, and the error would surface much later.
stopifnot(all(c("sample", "group", "batch") %in% names(ss)))

ss <- data.frame(
  sample = ss$sample,
  group  = ss$group,
  batch  = ss$batch,
  row.names = ss$sample,
  stringsAsFactors = FALSE
)

stopifnot(setequal(rownames(ss), colnames(counts)))
counts <- counts[, rownames(ss)]                           # same order
print(dim(counts)); print(table(ss$group))

# 2. choose the subset
keep_s <- ss$group %in% GROUPS & ss$batch %in% BATCHES
ss     <- ss[keep_s, , drop = FALSE]                       # drop=FALSE: stay a data.frame
counts <- counts[, rownames(ss), drop = FALSE]

# Factors come after the subsetting, so the groups we just dropped leave no empty level
# behind. HC goes first because R treats the first level as the baseline, which is what
# makes every group coefficient in 04 read "<group> vs HC".
# sort() fixes the order of the other two (HC, GBM, Lung -- same order relevel() gave).
# That order matters: results(dds) with no arguments in 04 returns the LAST group
# coefficient, so this line is what makes Lung vs HC the default result there.
ss$group <- factor(ss$group, levels = c(REF, sort(setdiff(GROUPS, REF))))
ss$batch <- factor(ss$batch, levels = BATCHES)

# 3. the property the whole analysis rests on: no empty group x batch cell.
cells <- table(ss$group, ss$batch)
print(cells)
stopifnot(
  "a group x batch cell is empty -- 04 could not separate batch from group" = all(cells > 0),
  "counts and sample sheet disagree on sample count" = ncol(counts) == nrow(ss),
  "counts columns are not in sample sheet order" = identical(colnames(counts), rownames(ss)),
  "HC is not the reference level" = levels(ss$group)[1] == REF,
  setequal(levels(ss$group), GROUPS),
  setequal(levels(ss$batch), BATCHES)
)
# Expected for this GROUPS x BATCHES configuration; change the constants, change this.
stopifnot("unexpected subset size" = ncol(counts) == 133)

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(counts = counts, ss = ss), OUT)

# Same subset again as plain text, only so it can be opened and looked at (Excel, a
# Python notebook, less). The .rds above stays the file 03 and 04 read: reading these
# back would lose the level order, so HC would no longer be the group compared against.
# Two files because the .rds holds two objects. Both are gitignored.
# write.csv(counts, sub("\\.rds$", "_counts.csv", OUT))                     # genes x samples
# write.csv(ss, sub("\\.rds$", "_samples.csv", OUT), row.names = FALSE)     # 133 x 3
