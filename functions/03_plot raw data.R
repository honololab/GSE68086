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
# this script adds four more columns to ss (lib_size, n_detected, top10_share,
# zero_share), each one a number worked out from counts and attached to the right
# sample.
#
#   all_ss   the one exception to "everything comes from the .rds": the full sample
#            sheet 01 wrote, all 285 samples of the series and all seven groups. Only
#            the last plot uses it, to show the three groups against the four the
#            subset leaves out. No counts are read for those other samples.
# -----------------------------------------------------------------------------

library(ggplot2)

IN     <- "data/derived/subset_lung_gbm_hc.rds"
SHEET  <- "GSE68086_sample_sheet.csv"          # all 285 samples, written by 01
FIGDIR <- "data/derived/figures"
FIGW   <- 7.6                                  # one width for every figure, in inches

stopifnot("subset not found -- run 02 first" = file.exists(IN))
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(IN)
counts <- d$counts
ss     <- d$ss
print(dim(counts))

# One colour per group, tied to the group NAME. A group keeps its colour in every plot,
# even one that leaves a group out, so "blue is HC" stays true throughout.
#
# These are soft pastels, chosen for how they look. Measured in OKLab, the distance
# between every pair comes out between 8.2 and 14.1 for normal vision, under the 15
# that is usually taken as the floor for telling two categories apart on sight, and
# none of them reaches a 3:1 contrast against the page (1.5 to 2.2). Under simulated
# colour blindness the three groups stay apart, with one exception: HC against the
# neutral in plot 9 falls to 4.3 for deuteranopia.
#
# So colour is a hint here, not the thing carrying the reading. That is fine because
# nothing in this script asks colour to work alone -- every plot names its groups on an
# axis or in a legend, the bars carry their counts, and the same counts are printed to
# the console. Anyone who needs the colours to separate on their own should raise the
# three away from each other and away from MUTE.
GROUP_COL <- c(HC = "#e79897", GBM = "#b7cbdb", Lung = "#fcc88a")
INK  <- "#0b0b0b"      # titles
INK2 <- "#52514e"      # axis text, labels
SURF <- "#ffffff"      # page colour, also the ring around each dot
RED  <- "#768e78"      # reference lines only, never a group
MUTE <- "#c6c09c"      # a bar that carries no group meaning

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

# ggplot does not wrap a title or a subtitle. A long one is drawn on one line and the
# end of it simply runs off the right edge of the image, with no warning.
#
# A title has to be there, so a long one is broken into lines. A subtitle does not: one
# that needs two lines crowds the plot, so it is dropped instead of wrapped. That is the
# rule -- a subtitle either fits on one line or it does not appear. Nothing is lost by
# dropping one, because the numbers they carry are printed to the console too, and the
# console line below names any that went, so a subtitle never disappears quietly.
#
# The two divisors are how many points wide a character is at each size. They are
# measured, not guessed: strwidth() on the same font and sizes the theme uses, over the
# actual titles and subtitles below, gives at most 6.25 points per character for the
# 12pt bold title and 4.17 for the 9.5pt subtitle. Both are rounded up a little so a
# line of unusually wide words still fits.
#
# Measured against that, none of the subtitles below is anywhere near the edge -- the
# longest, on plot 7, fills 79% of the drawable width. So nothing is dropped today; the
# rule is here to catch the next one that is written too long.
fit_lines <- function(s, chars) if (is.null(s)) NULL else strwrap(s, width = chars)

# Every figure is also kept in this list as it is drawn: the plot itself, its file name
# and its size. 03b reads the list to write the same set again without titles, for
# LaTeX. Doing it this way rather than repeating the nine calls over there means the
# two sets cannot drift apart -- a plot added below turns up in both, and a height
# changed here changes in both.
FIGURES <- list()

# Every figure is the same width, so they can be put side by side or dropped into a
# document together without one of them coming out at a different scale. Only the
# height changes, and only where a plot needs a squarer shape.
show_and_save <- function(p, name, w = FIGW, h = 4.4) {
  avail <- w * 72 - 26                  # the figure, less the left and right margins
  title <- fit_lines(p$labels$title,    floor(avail / 6.6))
  sub   <- fit_lines(p$labels$subtitle, floor(avail / 4.4))

  p$labels$title <- if (is.null(title)) NULL else paste(title, collapse = "\n")
  if (!is.null(sub) && length(sub) > 1) {
    cat(sprintf("  subtitle dropped from %s -- needs %d lines\n", name, length(sub)))
    p$labels$subtitle <- NULL
  }

  FIGURES[[length(FIGURES) + 1]] <<- list(p = p, name = name, w = w, h = h)

  print(p)
  out <- file.path(FIGDIR, name)
  ggsave(out, p, width = w, height = h, dpi = 150, bg = SURF)
  cat(sprintf("  saved %s\n", out))
}

# The median labels sit inside the cloud of dots, where plain text is unreadable the
# moment a dot lands behind a digit. geom_label puts an opaque patch behind each number
# instead. label.size = 0 drops the border the patch would otherwise get, so what is
# left is the number on a clean piece of page. It does cover a dot or two; the median
# bar underneath stays visible, and an unreadable number would be worth less than the
# dots it hides.
# x and y are not passed: they are inherited from the plot, which works because each
# median table carries the same column names as ss. Only the text and, where the bars
# are dodged, the grouping have to be named here.
median_label <- function(data, mapping, dodge = NULL) {
  pos <- if (is.null(dodge)) "identity" else position_dodge(width = dodge)
  geom_label(data = data, mapping = mapping, position = pos, vjust = -0.35,
             size = 2.7, colour = INK, fill = SURF, linewidth = 0,
             label.padding = unit(0.08, "lines"))
}



# ---- per-sample numbers added to ss ------------------------------------------
ss$lib_size    <- colSums(counts)              # total reads in that sample
ss$n_detected  <- colSums(counts > 0)          # genes with at least 1 read
top10          <- names(head(sort(rowSums(counts), decreasing = TRUE), 10))
ss$top10_share <- colSums(counts[top10, ]) / ss$lib_size
ss$zero_share  <- colMeans(counts == 0)        # genes with no reads, as a share of all

# ---- 1. depth by batch AND group ---------------------------------------------
# Depth is plotted against batch and group together, because depth here tracks batch
# and the groups are not spread evenly across batches. A plot of depth by group alone
# would show a difference that is really a batch difference.
# Cells hold 7-24 samples, too few for a box plot to say anything honest, so every
# sample is drawn and the bar is just the median.
# The number each median bar stands for. The bars on their own let you compare one cell
# against another; the number lets you quote one. Two decimals, rounded.
med_depth     <- aggregate(lib_size ~ batch + group, ss, median)
med_depth$lab <- sprintf("%.2f", med_depth$lib_size / 1e6)

p1 <- ggplot(ss, aes(batch, lib_size / 1e6)) +
  # an errorbar with min = max = median draws a plain flat median line
  stat_summary(aes(colour = group), fun = median, fun.min = median, fun.max = median,
               geom = "errorbar", width = 0.5, linewidth = 0.5,
               position = position_dodge(width = 0.75)) +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.1,
             position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.75,
                                             seed = 1)) +
  # dodged by the same width as the bar it belongs to, so it lands over its own cell
  median_label(med_depth, aes(label = lab, group = group), dodge = 0.75) +
  scale_fill_manual(values = GROUP_COL) +
  scale_colour_manual(values = GROUP_COL, guide = "none") +
  labs(title = "Profundidad de secuenciación por batch y grupo",
       subtitle = "Un punto por muestra, la barra es la mediana. La profundidad sigue al batch más que al grupo.",
       x = NULL, y = "Tamaño de la librería (millones de lecturas)") +
  theme_qc
show_and_save(p1, "01_depth_by_batch_group.png")

# ---- 2. does more depth buy more genes? --------------------------------------
# Two separate box plots of depth and of genes detected cannot answer this; the
# relationship between them needs both on one pair of axes.
rho <- suppressWarnings(cor(ss$lib_size, ss$n_detected, method = "spearman"))
p2 <- ggplot(ss, aes(lib_size / 1e6, n_detected)) +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.5) +
  scale_fill_manual(values = GROUP_COL) +
  labs(title = "Genes detectados frente a la profundidad de secuenciación",
       subtitle = sprintf(
         "Correlación de Spearman %.2f sobre %d muestras. Más profundidad no encuentra muchos más genes.",
         rho, nrow(ss)),
       x = "Tamaño de la librería (millones de lecturas)", y = "Genes con al menos 1 lectura") +
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
ALL <- "Las 133 juntas"
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
           label = sprintf("los 10 genes más altos: %.0f%% de las lecturas", 100 * top10_all)) +
  annotate("text", x = 1.15, y = 0.55, hjust = 0, size = 3, colour = INK2,
           label = sprintf("la mitad de las lecturas: %d genes", genes_half)) +
  scale_x_log10(breaks = c(1, 10, 100, 1000, 10000),
                labels = c("1", "10", "100", "1.000", "10.000")) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     limits = c(0, 1)) +
  scale_colour_manual(values = CURVE_COL) +
  labs(title = "Unos pocos genes concentran casi todas las lecturas",
       subtitle = "Genes ordenados por lecturas totales y su parte acumulada. Gris = las 133 muestras juntas.",
       x = "Número de genes (del más abundante al menos)",
       y = "Parte de todas las lecturas") +
  theme_qc
show_and_save(p3, "03_read_concentration.png")

# ---- 4. the same concentration, per sample -----------------------------------
# The curve above averages over samples and hides how much they differ. This shows the
# spread. No legend: the x axis already names the groups.
med_top10     <- aggregate(top10_share ~ group, ss, median)
med_top10$lab <- sprintf("%.2f%%", 100 * med_top10$top10_share)

p4 <- ggplot(ss, aes(group, top10_share)) +
  stat_summary(aes(colour = group), fun = median, fun.min = median, fun.max = median,
               geom = "errorbar", width = 0.45, linewidth = 0.5) +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.1,
             position = position_jitter(width = 0.15, height = 0, seed = 1)) +
  # no dodge here: one bar per group, so the label sits straight above it
  median_label(med_top10, aes(label = lab)) +
  scale_fill_manual(values = GROUP_COL, guide = "none") +
  scale_colour_manual(values = GROUP_COL, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "La concentración de lecturas varía mucho entre muestras",
       subtitle = "Parte de las lecturas de cada muestra que cae en los 10 genes más abundantes",
       x = NULL, y = "Parte de las lecturas en los 10 genes más altos") +
  theme_qc
show_and_save(p4, "04_top10_share_per_sample.png")

# ---- 5. mean against variance, per gene --------------------------------------
# One dot per GENE this time, not per sample: its average count over the 133 samples on
# x, the variance of those same 133 counts on y. The reference line is where variance
# equals the mean. That line is what counts would do if the only thing moving them was
# the random chance of a read landing on one gene rather than another. Dots above the
# line are more spread out than that, which is the thing a count model has to allow for.
#
# Both axes are log10 because the means run from 0.008 to 148,000 reads and the
# variances over twelve powers of ten; on straight axes every dot but a few hundred
# would pile up in one corner.
#
# 31,559 of the 57,736 genes have no reads at all in any of the 133 samples (54.7%).
# Their mean and variance are both exactly 0, and 0 has no place on a log axis, so they
# are dropped from this plot only -- nothing downstream is filtered here. The 26,177
# genes left are every gene seen at least once.
#
# The variance is taken across all 133 samples together, so it also holds the real
# differences between HC, GBM and Lung and the differences in sequencing depth seen in
# plot 1. It is not measurement noise on its own.
n_samp    <- ncol(counts)
gene_mean <- rowMeans(counts)
# Variance written out by hand rather than apply(counts, 1, var): the same number, from
# two passes over the matrix instead of 57,736 separate function calls.
gene_var  <- (rowSums(counts^2) - n_samp * gene_mean^2) / (n_samp - 1)

mv       <- data.frame(mean = gene_mean, var = gene_var)[gene_mean > 0, ]
over_pct <- 100 * mean(mv$var > mv$mean)

# Tick labels as powers of ten. Spelling out 10,000,000,000 would not fit, and the same
# style on both axes keeps that line readable as the place where the two axes agree.
pow10 <- function(x) parse(text = sprintf("10^%d", round(log10(x))))

p5 <- ggplot(mv, aes(mean, var)) +
  # 26,177 dots overlap heavily, so they are small and mostly transparent: where the
  # cloud looks solid, that is many genes, not one.
  geom_point(colour = INK2, alpha = 0.12, size = 0.35, shape = 16) +
  geom_abline(slope = 1, intercept = 0, colour = RED, linewidth = 0.6) +
  # label placed just under the line, in the gap between it and the cloud, so it reads
  # as a name for the line rather than as a note about the dots
  annotate("text", x = 2e3, y = 1.2e2, hjust = 0, size = 3, colour = RED,
           label = "varianza = media") +
  scale_x_log10(breaks = 10^(-2:5), labels = pow10) +
  scale_y_log10(breaks = 10^(seq(-2, 10, by = 2)), labels = pow10) +
  labs(title = "Los conteos varían mucho más de lo que predice la media",
       subtitle = sprintf(
         "Un punto por gen, %s genes con al menos una lectura. El %.0f%% queda por encima de la línea.",
         # Spanish number: dot for thousands, so the decimal mark has to move too
         format(nrow(mv), big.mark = ".", decimal.mark = ","), over_pct),
       x = "Media de los conteos en las 133 muestras",
       y = "Varianza de esos conteos") +
  theme_qc
show_and_save(p5, "05_mean_vs_variance.png", h = 5.6)   # taller: twelve powers of ten

# ---- 6. how many samples of each kind ----------------------------------------
# The plainest fact about the subset, and the one every later plot depends on: a median
# taken over 38 samples is a steadier number than one taken over 7. The count is printed
# as well as drawn so it can be read off exactly.
freq <- data.frame(group = factor(levels(ss$group), levels = levels(ss$group)),
                   n     = as.integer(table(ss$group)))
freq$share <- freq$n / sum(freq$n)

cat("\nsamples per group\n")
print(table(ss$group))

p6 <- ggplot(freq, aes(group, n, fill = group)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("%d  (%.0f%%)", n, 100 * share)),
            vjust = -0.55, size = 3.1, colour = INK2) +
  scale_fill_manual(values = GROUP_COL, guide = "none") +
  # headroom at the top so the labels are not cut off by the panel edge
  scale_y_continuous(expand = expansion(mult = c(0, 0.13))) +
  labs(title = "Los controles son el grupo más grande, el 41% del subconjunto",
       subtitle = sprintf(
         "%d muestras en total; los dos cánceres juntos son %d. Ningún grupo es demasiado pequeño.",
         sum(freq$n), sum(freq$n[freq$group != "HC"])),
       x = NULL, y = "Muestras") +
  theme_qc
show_and_save(p6, "06_samples_per_group.png", h = 4.2)

# ---- 7. which groups sit in which batch --------------------------------------
# A batch is a processing and sequencing run, not a biological thing. If one batch were
# mostly one group, any difference between groups would also be a difference between
# runs and the two could not be told apart. So the question is not "is there a batch
# effect" -- plot 1 already shows depth follows the batch -- but "can group and batch be
# separated at all".
#
# Each bar is one batch, cut up by group, so the bars answer that directly. The fourth
# bar is all 133 samples pooled, the mix each batch would show if the split were even.
# It is drawn in the same group colours rather than a fourth colour because it is the
# same three groups, just counted over everything; an empty slot on the axis keeps it
# from being read as a fourth batch.
tb <- table(group = ss$group, batch = ss$batch)

cat("\nsamples per group and batch\n")
print(tb)
cat("\ncomposition of each batch, % of its samples\n")
print(round(100 * prop.table(tb, 2), 1))

# Every expected count is above 10, so the chi-squared approximation holds here and the
# p value can be read as it stands.
p_batch <- chisq.test(tb)$p.value

ALLB <- "Todas las muestras"
GAP  <- " "        # a level with no data in it, which leaves an empty slot on the axis
comp <- as.data.frame(tb, responseName = "n")
comp <- rbind(comp, data.frame(group = freq$group, batch = ALLB, n = freq$n))
comp$batch <- factor(comp$batch, levels = c(levels(ss$batch), GAP, ALLB))
comp$share <- ave(comp$n, comp$batch, FUN = function(x) x / sum(x))

# n on the axis because shares alone would hide that the batches are different sizes.
# The empty slot has no samples, so tapply gives NA there and it gets a blank label.
bn        <- tapply(comp$n, comp$batch, sum)
batch_lab <- setNames(ifelse(is.na(bn), "", sprintf("%s\nn = %d", names(bn), bn)),
                      names(bn))

# reverse = TRUE puts HC at the bottom, so the stack reads in the same order as the
# legend. The count labels use the same setting or they land on the wrong segment.
p7 <- ggplot(comp, aes(batch, share, fill = group)) +
  geom_col(width = 0.6, position = position_stack(reverse = TRUE),
           colour = SURF, linewidth = 0.9) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5, reverse = TRUE),
            colour = INK, size = 3.1) +
  scale_fill_manual(values = GROUP_COL) +
  # drop = FALSE keeps the empty slot; the default would throw the unused level away
  scale_x_discrete(labels = batch_lab, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Todos los grupos están en todos los batches, así que se pueden separar",
       subtitle = sprintf(
         "Los números son muestras. La mezcla cambia entre batches, pero no más de lo que daría el azar (p = %.2f).",
         p_batch),
       x = NULL, y = "Parte del batch") +
  theme_qc
show_and_save(p7, "07_group_by_batch.png", h = 4.6)

# ---- 8. how much of each sample is nothing at all ----------------------------
# The share of genes with no reads in a sample. This is the same measurement as plot 2's
# y axis turned upside down -- zero_share is exactly 1 - n_detected / 57,736 -- so the
# new thing here is not the quantity but the split: it is drawn against batch AND group
# together, like plot 1, because that is the only way to see which of the two it follows.
# The answer: Lung has the highest median inside every one of the three batches (89.0%,
# 88.3%, 85.8%), so this one tracks the group. Batch02 does sit a little above the other
# two batches, but that is the batch holding the most Lung samples, 19 of its 47.
zero_med  <- median(ss$zero_share)
zero_lung <- median(ss$zero_share[ss$group == "Lung"])
zero_hc   <- median(ss$zero_share[ss$group == "HC"])

med_zero     <- aggregate(zero_share ~ batch + group, ss, median)
med_zero$lab <- sprintf("%.2f", 100 * med_zero$zero_share)

p8 <- ggplot(ss, aes(batch, zero_share)) +
  stat_summary(aes(colour = group), fun = median, fun.min = median, fun.max = median,
               geom = "errorbar", width = 0.5, linewidth = 0.5,
               position = position_dodge(width = 0.75)) +
  geom_point(aes(fill = group), shape = 21, colour = SURF, stroke = 0.6, size = 2.1,
             position = position_jitterdodge(jitter.width = 0.14, dodge.width = 0.75,
                                             seed = 1)) +
  median_label(med_zero, aes(label = lab, group = group), dodge = 0.75) +
  scale_fill_manual(values = GROUP_COL) +
  scale_colour_manual(values = GROUP_COL, guide = "none") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
  labs(title = "En todas las muestras casi toda la lista de genes está vacía",
       subtitle = sprintf(
         "Mediana del %.0f%% de genes sin lecturas. Lung es el más alto en cada batch (%.0f%% frente a HC %.0f%%).",
         100 * zero_med, 100 * zero_lung, 100 * zero_hc),
       x = NULL, y = "Parte de los 57.736 genes sin lecturas") +
  theme_qc
show_and_save(p8, "08_zero_share_by_batch_group.png", h = 4.6)

# ---- 9. how many samples of every group in the series -------------------------
# Plot 6 counts the 133 samples this analysis uses. This one counts all 285 in the
# series, so those three groups can be read against the four that were left out. It is
# the only plot here that looks outside the subset file, and it reads the sample sheet
# only -- no counts for the other 152 samples are loaded.
#
# One bar per group and one number on each: a frequency plot answers one question, how
# many, and the height plus the count answer it. The only other thing the bars carry is
# colour, and it means what it means everywhere else in this script -- HC blue, GBM
# orange, Lung green. The four groups the analysis does not use stay neutral, so the
# three it does use can be picked out without reading the names.
#
# For the record, since the plot no longer shows it: 02 keeps Batch02/03/04, the only
# batches where HC, GBM and Lung all appear, and that batch choice is what removes 20
# Lung and 2 GBM samples from the 133. Of the rest, Breast has no Batch02 samples at
# all, Pancreas and Hepatobiliary have one each, too few to hold up a group x batch
# cell. CRC has 4, 21 and 13, so CRC would have worked -- leaving it out was a choice.
stopifnot("full sample sheet not found -- run 01 first" = file.exists(SHEET))
all_ss <- read.csv(SHEET, stringsAsFactors = FALSE)

# A check, not a drawing: the subset rule read back off the subset itself rather than
# 02's constants copied out by hand. If it stops reproducing the 133 samples, the two
# files have drifted apart and the count below is about a series this analysis no
# longer comes from.
kept <- all_ss$group %in% levels(ss$group) & all_ss$batch %in% levels(ss$batch)
stopifnot("the kept rule no longer reproduces the subset" = sum(kept) == ncol(counts))

dis <- data.frame(group = names(table(all_ss$group)),
                  n     = as.integer(table(all_ss$group)))
dis <- dis[order(-dis$n), ]
dis$group <- factor(dis$group, levels = dis$group)   # biggest group first

# A group keeps its colour from GROUP_COL; the rest fall back to the neutral. Built off
# the names in GROUP_COL rather than a list typed out again here, so a group added to
# the analysis later is coloured in this plot too without anyone remembering to do it.
BAR_COL  <- c(GROUP_COL, otros = MUTE)
dis$fill <- ifelse(as.character(dis$group) %in% names(GROUP_COL),
                   as.character(dis$group), "otros")
dis$fill <- factor(dis$fill, levels = names(BAR_COL))

cat("\nsamples per group across the whole series\n")
print(setNames(dis$n, as.character(dis$group)))

# Sideways. Seven group names on an upright x axis is where labels start touching each
# other, and "Hepatobiliary" is 13 characters; turned on its side every name gets a
# whole line to itself, stays horizontal and cannot run into its neighbour however
# long it is. The count sits just outside the end of each bar, clear of the bar and of
# the name.
p9 <- ggplot(dis, aes(n, group, fill = fill)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = n), hjust = -0.35, size = 3.2, colour = INK2) +
  scale_fill_manual(values = BAR_COL, guide = "none") +  # the axis already names them
  scale_x_continuous(expand = expansion(mult = c(0, 0.09))) +
  scale_y_discrete(limits = rev(levels(dis$group))) +   # biggest at the top
  labs(title = "Muestras por grupo en toda la serie",
       x = "Muestras", y = NULL) +
  theme_qc +
  theme(panel.grid.major.y = element_blank())           # nothing to line up along a bar

show_and_save(p9, "09_samples_per_group_all.png", h = 4.2)

cat(sprintf("\nfigures written to %s\n", FIGDIR))
