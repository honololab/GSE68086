# 04 -- DESeq2 object, gene filter, normalisation
# Takes the 133-sample subset that 02 saved and turns it into the object that 05 will
# test. Four steps, in this order:
#   1. counts + sample sheet + the model formula (~ batch + group) go into one object
#   2. genes with too few reads to be tested are removed (edgeR::filterByExpr)
#   3. one size factor per sample is estimated (DESeq2), and its "beyond depth" part is
#      set next to edgeR's TMM factor, group by group -- they turn out to disagree
#   4. counts are put on a log-like scale (vst) for one PCA, to see what structure is
#      left once depth has been taken out
# No differential-expression test runs here; that is 05. Figures are saved next to the
# ones from 03, numbered 05-08.
#
# Numbers quoted in the comments were measured on this 133-sample subset (GBM 38). The
# subset still contains VU398-GBM, a sample the authors excluded; the review of 00b
# suggests dropping it in 02, which would change those numbers slightly (6,994 -> 7,082
# genes kept) without changing any conclusion here.
#
# ---- what the objects are ----------------------------------------------------
#   d, counts, ss   exactly as in 03: the .rds written by 02, the 57,736 x 133 count
#                   matrix, and the 133 x 3 sample sheet (sample, group, batch) whose
#                   rows are in the same order as the columns of counts.
#   dds             a DESeqDataSet. It holds the counts (rows = genes, columns =
#                   samples), the sample sheet (colData(dds), one row per column of
#                   counts) and the design formula, all in one object. Subsetting
#                   dds[rows, ] or dds[, cols] moves all three together, so the labels
#                   cannot come apart from the columns. counts(dds) gives the raw matrix
#                   back; counts(dds, normalized = TRUE) gives it divided by the size
#                   factors.
#   keep_g          one TRUE/FALSE per gene: passes the filter or not.
#   vsd             dds after vst(): assay(vsd) is a matrix on a log2-like scale where
#                   the spread of a gene no longer grows with its mean. For plots only.
#                   05 tests the raw counts in dds, never these values.
# -----------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(edgeR)          # filterByExpr, DGEList, normLibSizes
  library(ggplot2)
})

IN      <- "data/derived/subset_lung_gbm_hc.rds"     # written by 02
FIGDIR  <- "data/derived/figures"                    # same folder as 03
OUT_DDS <- "data/derived/dds_filtered.rds"           # 05 starts from this
OUT_VSD <- "data/derived/vsd_blind.rds"              # PCA / heatmaps read this

stopifnot("subset not found -- run 02 first" = file.exists(IN))
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(IN)
counts <- d$counts # raw counts
ss     <- d$ss # sample sheet

# 02 checked these before saving. Checking again costs nothing, and it means an edited or
# stale .rds fails here instead of giving a wrong answer three steps later.
stopifnot(
  "counts columns are not in sample sheet order" = identical(colnames(counts), rownames(ss)),
  "HC is not the reference level"                = levels(ss$group)[1] == "HC",
  "a group x batch cell is empty"                = all(table(ss$group, ss$batch) > 0)
)
print(dim(counts))

# Three per-sample numbers from 03, worked out again because 03 keeps them only in its
# own session. All use ALL genes, before any filtering: a sample's depth is its whole column.
ss$lib_size    <- colSums(counts)              # total reads in that sample
ss$n_detected  <- colSums(counts > 0)          # genes with at least 1 read
top10          <- names(head(sort(rowSums(counts), decreasing = TRUE), 10))
ss$top10_share <- colSums(counts[top10, ]) / ss$lib_size

# ---- plot style: identical to 03 so the figures match --------------------------
GROUP_COL <- c(HC = "#2a78d6", GBM = "#eb6834", Lung = "#1baf7a")
INK  <- "#0b0b0b"      # titles
INK2 <- "#52514e"      # axis text, labels
SURF <- "#fcfcfb"      # page colour, also the ring around each dot

theme_qc <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = SURF, colour = NA),
    panel.background = element_rect(fill = SURF, colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "#e6e5e1"),
    strip.text       = element_text(colour = INK, size = 10),
    axis.title       = element_text(colour = INK2, size = 9.5),
    axis.text        = element_text(colour = INK2),
    plot.title       = element_text(colour = INK, face = "bold", size = 12),
    plot.subtitle    = element_text(colour = INK2, size = 9.5),
    legend.title     = element_blank(),
    legend.text      = element_text(colour = INK2),
    legend.position  = "top",
    plot.margin      = margin(12, 14, 10, 12)
  )

show_and_save <- function(p, name, w = 7.2, h = 4.4) {
  print(p)
  ggsave(file.path(FIGDIR, name), p, width = w, height = h, dpi = 150, bg = SURF)
}

# ---- 1. one object: counts + sample sheet + model ------------------------------
# ~ batch + group says: for every gene, explain its counts by which batch the sample was
# processed in and which group it belongs to. Batch is in the formula because 03 showed
# that depth follows batch and that the healthy controls sit mostly in the shallow batch;
# leaving batch out would let that processing difference read as a disease difference.
# The order of the two terms does not change the fit. HC is the first level of group (set
# in 02), so every group coefficient in 05 will read "<group> vs HC".
dds <- DESeqDataSetFromMatrix(countData = counts, colData = ss, design = ~ batch + group)

# ---- 2. remove genes with too few reads to test --------------------------------
# 03 found that 54.7% of genes have no read in any sample, and that HC libraries are about
# half as deep as Lung (median 1.04 M against 1.97 M reads). Both shape the filter.
#
# filterByExpr keeps a gene if it reaches a cutoff in enough samples. The cutoff is set in
# counts per million (CPM), not raw reads: 10 reads in the median library, about 5.9 CPM
# here. That is about 6 reads in a median HC library and 12 in a median Lung library, so
# every sample is asked the same question relative to its own depth. A raw cutoff (say
# ">= 10 reads") is twice as hard for an HC sample to pass, and 00c measured the cost:
# 886 genes expressed mainly in HC would be dropped, exactly the genes that go DOWN in
# cancer. So the cutoff is on the CPM scale.
#
# "Enough samples" comes from group=: the smallest group (38), reduced by edgeR's
# tolerance to 10 + 0.7 * 28 = 29.6 samples. With design= instead, the minimum would come
# from the smallest batch x group cell (about 19 samples), and a gene expressed in one
# cell alone would pass. That is a batch-driven gene, the kind we want out. Our contrasts
# are between groups, so the group rule matches the question. A gene must also reach 15
# reads in total.
#
# On this subset the rules compared in 00c keep: filterByExpr(group=) 6,994; the raw rule
# 6,109; filterByExpr(design=) 8,079; the looser CPM > 1 and ">= 5 reads in 10% of
# samples" rules about 9,700-9,900, mostly genes with a few reads in a minority of samples
# that cannot be tested. Details: GSE68086_review_00b_filtering_and_paper.md, section 2.

# The numbers filterByExpr works with, printed so the rule is not a black box. They follow
# the function's defaults (min.count = 10, min.total.count = 15, large.n = 10, min.prop = 0.7).
cpm_cut <- 10 / median(ss$lib_size) * 1e6
n_min   <- min(table(ss$group))
n_need  <- 10 + (n_min - 10) * 0.7
cat(sprintf(paste0(
  "filterByExpr: CPM cutoff %.2f = about %.0f reads in a median HC library, %.0f in a median Lung library\n",
  "              needed in at least %.1f samples (smallest group %d, after tolerance); total >= 15 reads\n"),
  cpm_cut,
  cpm_cut * median(ss$lib_size[ss$group == "HC"])   / 1e6,
  cpm_cut * median(ss$lib_size[ss$group == "Lung"]) / 1e6,
  n_need, n_min))

keep_g <- filterByExpr(counts(dds), group = dds$group)
print(table(keep_g))                       # 6,994 kept on this 133-sample subset

# What the filter costs in reads rather than rows: it drops most rows but almost no reads,
# because the rows it drops are the near-empty ones (03, finding 4).
kept_share <- colSums(counts(dds)[keep_g, ]) / ss$lib_size
cat("share of each sample's reads kept by the filter, median by group:\n")
print(round(tapply(kept_share, ss$group, median), 4))

dds <- dds[keep_g, ]

# ---- 3. size factors: DESeq2 next to edgeR's TMM -------------------------------
# DESeq2 gives every sample one number, its size factor, and divides the sample's counts
# by it. How it is found: for each gene, the sample's count is divided by that gene's
# geometric mean across all samples; the size factor is the median of those ratios over
# the genes that have a count in every sample. Two things follow.
#  - It mostly tracks library size. A sample with twice the reads has about twice the
#    count for every gene, so the ratios double and so does their median.
#  - It need not be only library size. If a sample's reads pile into a few very abundant
#    genes (03, finding 3), the ordinary gene gets fewer of them, and the size factor
#    should come out below what depth alone would give. Whether that happens here, and
#    whether the two normalisation methods see it the same way, is measured below.
dds <- estimateSizeFactors(dds)
ss$size_factor <- sizeFactors(dds)

# Only genes counted in every sample enter the median above. There should be well over a
# hundred; if this number were tiny the size factors would be unreliable and
# estimateSizeFactors(type = "poscounts") would be the alternative.
cat("\ngenes with a count in every sample (the ones the size factors are based on):",
    sum(rowSums(counts(dds) == 0) == 0), "of", nrow(dds), "\n")

# edgeR measures the same "beyond depth" part with TMM. Its factor is already relative to
# library size, so to compare, take the DESeq2 size factor and divide out the relative
# library size. Everything is scaled to a geometric mean of 1 first. A third column uses
# DESeq2's other method, poscounts, which builds each gene's geometric mean from all
# samples with a count instead of using only genes counted in every sample.
gm <- function(v) exp(mean(log(v)))
y  <- normLibSizes(DGEList(counts(dds), lib.size = ss$lib_size))
sf_pos <- sizeFactors(estimateSizeFactors(dds, type = "poscounts"))
ss$rel_lib      <- ss$lib_size / gm(ss$lib_size)                         # depth alone
ss$bd_ratio     <- (ss$size_factor / gm(ss$size_factor)) / ss$rel_lib    # DESeq2 default
ss$bd_poscounts <- (sf_pos / gm(sf_pos)) / ss$rel_lib                    # DESeq2 poscounts
ss$bd_tmm       <- y$samples$norm.factors                                # edgeR TMM

cat("\nmedian by group (bd_ = the part of the factor beyond depth):\n")
print(round(do.call(rbind, lapply(split(ss, ss$group), function(s) c(
  lib_size_M   = median(s$lib_size) / 1e6,
  size_factor  = median(s$size_factor),
  bd_ratio     = median(s$bd_ratio),
  bd_poscounts = median(s$bd_poscounts),
  bd_tmm       = median(s$bd_tmm)))), 3))

# Measured on this subset: bd_ratio is about 1.0 in all three groups, while bd_tmm is
# about HC 1.22, GBM 0.89, Lung 0.82 (as 00c found) and bd_poscounts about HC 1.25,
# GBM 0.84, Lung 0.78. So DESeq2's default sees no group-wide composition effect; TMM and
# poscounts see a large one. They disagree because they anchor on different genes. The
# default ratio method only uses genes counted in every sample, about 1,400 of the most
# abundant ones, and among those HC and Lung sit at the same level per million reads. TMM
# and poscounts use all kept genes, and among the other 5,600 genes HC is about twice as
# high per million reads as Lung. That fits 03, finding 3: cancer libraries spend more of
# their reads on the abundant genes, so less is left for the rest. The two lines below
# print that split.
allnz   <- rowSums(counts(dds) == 0) == 0
cpm_raw <- cpm(counts(dds), lib.size = ss$lib_size)          # depth-scaled only, no TMM
hc_lung <- log2(rowMeans(cpm_raw[, ss$group == "HC"])   + 0.5) -
           log2(rowMeans(cpm_raw[, ss$group == "Lung"]) + 0.5)
cat(sprintf(paste0(
  "median log2(HC / Lung) of mean CPM: %+.2f over the %d genes counted in every sample,\n",
  "                                    %+.2f over the other %d kept genes\n"),
  median(hc_lung[allnz]), sum(allnz), median(hc_lung[!allnz]), sum(!allnz)))
# This is a choice for 05 and the tool comparison to make, not something to settle
# quietly here. Under the default ratio method most of the 5,600 less-abundant genes will
# come out DOWN in cancer; under TMM or poscounts that shift is normalised away and the
# abundant genes come out UP in cancer instead. Either way the cancer-vs-HC lists will be
# long and will partly reflect library make-up, not only tumour biology; Lung-vs-GBM is
# the cleaner contrast. dds keeps the DESeq2 default. To switch:
# dds <- estimateSizeFactors(dds, type = "poscounts").

# ---- figure 5: what the size factors do to each sample -------------------------
# One box per sample: the 5th, 25th, 50th, 75th and 95th percentile of its log2 counts
# over the kept genes, before and after dividing by its size factor. Samples are split by
# batch and ordered by depth inside each panel, so the Batch04 dip from 03 (figure 1)
# should be visible in the top row and gone in the bottom row. Percentiles are computed
# first and drawn as-is; drawing 930,000 points twice would say the same thing slower.
five <- function(m, stage) {
  q <- apply(log2(m + 1), 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  data.frame(sample = colnames(m), group = ss$group, batch = ss$batch, stage = stage,
             ymin = q[1, ], lower = q[2, ], middle = q[3, ], upper = q[4, ], ymax = q[5, ])
}
box <- rbind(five(counts(dds),                    "Raw counts"),
             five(counts(dds, normalized = TRUE), "Divided by size factor"))
box$stage  <- factor(box$stage, levels = c("Raw counts", "Divided by size factor"))
box$sample <- factor(box$sample, levels = ss$sample[order(ss$batch, ss$lib_size)])

p5 <- ggplot(box, aes(sample, fill = group)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
               stat = "identity", colour = INK2, linewidth = 0.2, width = 0.75) +
  facet_grid(stage ~ batch, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = GROUP_COL) +
  labs(title = "Size factors line the samples up across batches",
       subtitle = paste0(
         "One box per sample: 5th to 95th percentile of log2 counts over the kept genes, ordered by depth within batch.\n",
         "The batch difference goes; the spread in box length stays. That is per-sample read concentration (03, figure 4)."),
       x = NULL, y = "log2(count + 1)") +
  theme_qc +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.grid.major.x = element_blank())
show_and_save(p5, "05_normalisation_before_after.png", w = 9, h = 5.2)

# ---- figure 6: size factor against depth ---------------------------------------
# The dashed line is what the size factor would be if it were depth alone. Distance from
# the line is the "beyond depth" part tabulated above, sample by sample.
above <- tapply(ss$bd_ratio > 1, ss$group, mean)
p6 <- ggplot(ss, aes(lib_size / 1e6, size_factor)) +
  geom_abline(intercept = 0, slope = gm(ss$size_factor) / (gm(ss$lib_size) / 1e6),
              linetype = "22", linewidth = 0.4, colour = "#9a9994") +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.5) +
  scale_fill_manual(values = GROUP_COL) +
  labs(title = "DESeq2's default size factors are depth and little else",
       subtitle = sprintf(paste0(
         "Dashed line = depth alone. Samples above it: HC %.0f%%, GBM %.0f%%, Lung %.0f%%.\n",
         "No group sits off the line. TMM sees it differently (table printed by the script)."),
         100 * above["HC"], 100 * above["GBM"], 100 * above["Lung"]),
       x = "Library size (million reads)", y = "DESeq2 size factor") +
  theme_qc
show_and_save(p6, "06_size_factor_vs_depth.png")

# ---- 4. vst and one PCA --------------------------------------------------------
# Raw counts spread more the higher a gene's mean, so a PCA on them would be a PCA of the
# few most abundant genes. vst() puts every gene on a scale where that no longer happens,
# and applies the size factors on the way. blind = TRUE means it does not look at the
# design, which is right when the question is what structure the data has on its own.
vsd <- vst(dds, blind = TRUE)

# PCA on all kept genes, centred, not scaled. No "top N most variable genes" cut: it is
# one fewer choice to justify, and with 6,994 genes it is not needed for speed.
pca <- prcomp(t(assay(vsd)))
pct <- 100 * pca$sdev^2 / sum(pca$sdev^2)
cat(sprintf("\nvariance carried by PC1-5: %s\n",
            paste0(round(pct[1:5], 1), "%", collapse = ", ")))
pc <- cbind(ss, pca$x[, 1:3])

# How much of each PC is group, batch, depth, genes detected, or read concentration?
# R-squared for the two factors, Spearman correlation for the three per-sample numbers.
# Group and batch overlap (HC sits mostly in Batch04), so their R-squared values are not
# additive. The last row is the check 03 left open: do the samples with very concentrated
# reads (top10_share from 12% to 50%) still stand out once depth has been normalised?
# Measured: yes, strongly. PC1 (39% of the variance) follows top10_share with Spearman
# 0.82 and genes detected with -0.87 (two faces of the same thing: concentrated reads
# leave more genes at zero), far more than it follows group (R-squared 0.36). PC2 is batch
# (R-squared 0.46), PC3 is depth (Spearman 0.69). Concentration is mostly a per-sample
# property, only weakly tied to group (03, finding 4), so it is noise that ~ batch + group
# does not remove. Worth considering as a covariate in 05 -- to review, not decided here.
r2 <- function(p, f) summary(lm(pc[[p]] ~ pc[[f]]))$r.squared
assoc <- sapply(c("PC1", "PC2", "PC3"), function(p) c(
  group_R2          = r2(p, "group"),
  batch_R2          = r2(p, "batch"),
  depth_spearman    = cor(pc[[p]], pc$lib_size,    method = "spearman"),
  detected_spearman = cor(pc[[p]], pc$n_detected,  method = "spearman"),
  top10_spearman    = cor(pc[[p]], pc$top10_share, method = "spearman")))
cat("\nwhat each PC lines up with:\n")
print(round(assoc, 2))

p7 <- ggplot(pc, aes(PC1, PC2)) +
  geom_point(aes(fill = group, shape = batch), colour = SURF, stroke = 0.5, size = 2.6) +
  scale_fill_manual(values = GROUP_COL) +
  scale_shape_manual(values = c(Batch02 = 21, Batch03 = 22, Batch04 = 24)) +
  guides(fill  = guide_legend(override.aes = list(shape = 21)),
         shape = guide_legend(override.aes = list(fill = INK2, colour = INK2))) +
  labs(title = "PCA on vst values, all kept genes",
       subtitle = sprintf(
         "Colour = group, shape = batch. PC1 follows read concentration (Spearman %.2f), PC2 follows batch (R-squared %.2f).",
         assoc["top10_spearman", "PC1"], assoc["batch_R2", "PC2"]),
       x = sprintf("PC1 (%.0f%% of variance)", pct[1]),
       y = sprintf("PC2 (%.0f%% of variance)", pct[2])) +
  theme_qc
show_and_save(p7, "07_pca_vst.png", w = 8, h = 5)

# ---- figure 8: the open question from 03, answered -----------------------------
# PC1 against the share of a sample's reads in the ten most abundant genes (the quantity
# figure 4 in 03 plotted per group). If the points fall on a line, the main axis of
# variation after normalisation is how concentrated a sample's reads are, not which group
# it is in. Groups would then show as a vertical offset at a given concentration, not as
# separate clouds.
p8 <- ggplot(pc, aes(top10_share, PC1)) +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.5) +
  scale_fill_manual(values = GROUP_COL) +
  scale_x_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "PC1 is read concentration",
       subtitle = sprintf("Spearman correlation %.2f over %d samples.",
                          assoc["top10_spearman", "PC1"], nrow(pc)),
       x = "Share of the sample's reads in the top 10 genes",
       y = sprintf("PC1 (%.0f%% of variance)", pct[1])) +
  theme_qc
show_and_save(p8, "08_pc1_vs_top10_share.png")

# ---- export --------------------------------------------------------------------
# dds carries the filtered counts, the sample sheet (colData now also has lib_size,
# n_detected, top10_share and sizeFactor) and the design. 05 starts from it. vsd is saved
# so plots do not have to recompute it. Both are gitignored, like everything under
# data/derived.
# saveRDS(dds, OUT_DDS)
# saveRDS(vsd, OUT_VSD)
# cat(sprintf("\nwritten %s (%d genes x %d samples) and %s\nfigures written to %s\n",
#             OUT_DDS, nrow(dds), ncol(dds), OUT_VSD, FIGDIR))
print(dds)
