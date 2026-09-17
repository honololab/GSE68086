# 03 -- plot raw data
# Looks at the 133-sample subset before any normalising or gene filtering, so these
# plots are about data quality, not biology. Saves each figure to data/derived/figures
# as well as drawing it, so they can be compared between runs.
#
# ---- what the objects are ----------------------------------------------------
# Everything comes out of the one .rds file that 02 wrote.
#
#   d        the subset.rds file itself, read back as a list with exactly two items:
#            d$counts and d$ss. Nothing else is in it.
#
#   counts   (= d$counts) the read counts. A table of 57,736 rows by 133 columns:
#            one ROW per gene, one COLUMN per sample. Each cell is how many reads
#            landed on that gene in that sample. Row names are Ensembl gene ids
#            (ENSG...), column names are sample names (VU280-GBM, HD-27-2, ...).
#            Numbers only -- it does not know which sample is a patient or a
#            healthy control.
#
#   ss       (= d$ss) the sample sheet, which holds those labels. 133 rows, one per
#            sample, with three columns: sample, group (HC / GBM / Lung) and batch
#            (Batch02 / 03 / 04). Its rows are in the same order as the columns of
#            counts, which is what 02 checked before saving, so sample i in ss is
#            column i in counts.
#
# So: counts has the measurements, ss says what each measurement belongs to. Below,
# this script adds three more columns to ss (lib_size, n_detected, top10_share), each
# one a number worked out from counts and attached to the right sample.
# -----------------------------------------------------------------------------

library(ggplot2)

IN     <- "data/derived/subset_lung_gbm_hc.rds"
FIGDIR <- "data/derived/figures"

stopifnot("subset not found -- run 02 first" = file.exists(IN))
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(IN)
counts <- d$counts
ss     <- d$ss
print(dim(counts))

# One colour per group, tied to the group NAME. A group keeps its colour in every plot,
# even one that leaves a group out, so "blue is HC" stays true throughout. These three
# stay apart for colour-blind readers; they are the first three slots of a palette that
# was checked as a set rather than judged by eye.
GROUP_COL <- c(HC = "#2a78d6", GBM = "#eb6834", Lung = "#1baf7a")
INK  <- "#0b0b0b"      # titles
INK2 <- "#52514e"      # axis text, labels
SURF <- "#fcfcfb"      # page colour, also the ring around each dot

# Grid lines are thin, solid and pale so the data sits in front of them.
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

# ---- per-sample numbers added to ss ------------------------------------------
ss$lib_size    <- colSums(counts)              # total reads in that sample
ss$n_detected  <- colSums(counts > 0)          # genes with at least 1 read
top10          <- names(head(sort(rowSums(counts), decreasing = TRUE), 10))
ss$top10_share <- colSums(counts[top10, ]) / ss$lib_size

# ---- 1. depth by batch AND group ---------------------------------------------
# Depth is plotted against batch and group together, because depth here tracks batch
# and the groups are not spread evenly across batches. A plot of depth by group alone
# would show a difference that is really a batch difference.
# Cells hold 7-24 samples, too few for a box plot to say anything honest, so every
# sample is drawn and the bar is just the median.
p1 <- ggplot(ss, aes(batch, lib_size / 1e6)) +
  # an errorbar with min = max = median draws a plain flat median line
  stat_summary(aes(colour = group), fun = median, fun.min = median, fun.max = median,
               geom = "errorbar", width = 0.5, linewidth = 0.5,
               position = position_dodge(width = 0.75)) +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.1,
             position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.75,
                                             seed = 1)) +
  scale_fill_manual(values = GROUP_COL) +
  scale_colour_manual(values = GROUP_COL, guide = "none") +
  labs(title = "Sequencing depth by batch and group",
       subtitle = "One dot per sample, bar = median. Depth follows the batch more than the group.",
       x = NULL, y = "Library size (million reads)") +
  theme_qc
show_and_save(p1, "01_depth_by_batch_group.png")

# ---- 2. does more depth buy more genes? --------------------------------------
# Two separate box plots of depth and of genes detected cannot answer this; the
# relationship between them needs both on one pair of axes.
rho <- suppressWarnings(cor(ss$lib_size, ss$n_detected, method = "spearman"))
p2 <- ggplot(ss, aes(lib_size / 1e6, n_detected)) +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.5) +
  scale_fill_manual(values = GROUP_COL) +
  labs(title = "Genes detected against sequencing depth",
       subtitle = sprintf(
         "Spearman correlation %.2f over %d samples. Deeper libraries are not finding many more genes.",
         rho, nrow(ss)),
       x = "Library size (million reads)", y = "Genes with at least 1 read") +
  theme_qc
show_and_save(p2, "02_genes_vs_depth.png")

# ---- 3. how concentrated are the reads? --------------------------------------
# Ranks the genes by total reads, then plots the running share of all reads. The x axis
# is on a log scale because almost everything happens in the first few hundred genes.
ranks <- unique(round(10^seq(0, log10(nrow(counts)), length.out = 400)))

running_share <- function(m) {
  tot <- sort(rowSums(m), decreasing = TRUE)
  (cumsum(tot) / sum(tot))[ranks]
}

# One curve per group, PLUS a curve for all 133 samples pooled. The pooled one is drawn
# because the headline numbers quoted elsewhere (top 10 = 29% of reads, 50% of reads in
# 57 genes) are pooled figures. Without it on the chart those numbers match none of the
# lines, which is misleading. It is grey, not a fourth group colour, because it is a
# reference line rather than another group.
ALL <- "All 133 pooled"
cum <- do.call(rbind, c(
  lapply(levels(ss$group), function(g) {
    data.frame(group = g, rank = ranks,
               share = running_share(counts[, ss$group == g, drop = FALSE]))
  }),
  list(data.frame(group = ALL, rank = ranks, share = running_share(counts)))
))
cum$group <- factor(cum$group, levels = c(levels(ss$group), ALL))
CURVE_COL <- c(GROUP_COL, setNames("#52514e", ALL))

# The two numbers the annotations point at, both worked out from the pooled counts so
# they describe the grey line.
pooled     <- sort(rowSums(counts), decreasing = TRUE)
pooled_cum <- cumsum(pooled) / sum(pooled)
top10_all  <- pooled_cum[10]
genes_half <- sum(pooled_cum < 0.5) + 1

p3 <- ggplot(cum, aes(rank, share, colour = group)) +
  geom_hline(yintercept = 0.5, linewidth = 0.25, colour = "#c9c8c3") +
  geom_vline(xintercept = 10, linewidth = 0.25, colour = "#c9c8c3") +
  geom_line(linewidth = 0.7) +
  annotate("text", x = 11, y = 0.05, hjust = 0, size = 3, colour = INK2,
           label = sprintf("top 10 genes: %.0f%% of all reads", 100 * top10_all)) +
  annotate("text", x = 1.15, y = 0.55, hjust = 0, size = 3, colour = INK2,
           label = sprintf("half of all reads: %d genes", genes_half)) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000),
                labels = c("1", "10", "100", "1,000", "10,000")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, 1)) +
  scale_colour_manual(values = CURVE_COL) +
  labs(title = "A handful of genes hold most of the reads",
       subtitle = "Genes ranked by total reads, then the running share. Grey = all 133 samples pooled.",
       x = "Number of genes (ranked most abundant first)",
       y = "Share of all reads") +
  theme_qc
show_and_save(p3, "03_read_concentration.png")

# ---- 4. the same concentration, per sample -----------------------------------
# The curve above averages over samples and hides how much they differ. This shows the
# spread. No legend: the x axis already names the groups.
p4 <- ggplot(ss, aes(group, top10_share)) +
  stat_summary(aes(colour = group), fun = median, fun.min = median, fun.max = median,
               geom = "errorbar", width = 0.45, linewidth = 0.5) +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.1,
             position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  scale_fill_manual(values = GROUP_COL, guide = "none") +
  scale_colour_manual(values = GROUP_COL, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "Read concentration varies widely between samples",
       subtitle = "Share of each sample's reads sitting in the 10 most abundant genes overall",
       x = NULL, y = "Share of reads in the top 10 genes") +
  theme_qc
show_and_save(p4, "04_top10_share_per_sample.png")

cat(sprintf("\nfigures written to %s\n", FIGDIR))
