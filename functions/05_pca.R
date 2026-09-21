# 05 -- PCA on the vst values, drawn twice: coloured by group, then by batch

library(DESeq2)
library(ggplot2)

IN_DDS <- "data/derived/dds_filtered.rds" # written by 04
IN_VSD <- "data/derived/vsd_blind.rds" # written by 04
FIGDIR <- "data/derived/figures"
TEXDIR <- "data/derived/figures_latex"
FIGW   <- 7.6

stopifnot("dds not found -- run 04 first" = file.exists(IN_DDS),
          "vsd not found -- run 04 first" = file.exists(IN_VSD))
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TEXDIR, recursive = TRUE, showWarnings = FALSE)
 
dds <- readRDS(IN_DDS)
vsd <- readRDS(IN_VSD)
ss  <- as.data.frame(colData(dds))           # sample sheet: group, batch, ...

# --- style

GROUP_COL <- c(HC = "#e79897", GBM = "#b7cbdb", Lung = "#fcc88a")
# Batch colours are deliberately NOT the group colours, so the two legends are never
# mistaken for each other.
BATCH_COL <- c(Batch02 = "#7b8fa6", Batch03 = "#a88bb5", Batch04 = "#8fae8b")
INK  <- "#0b0b0b"
INK2 <- "#52514e"
SURF <- "#ffffff"
 
theme_qc <- theme_minimal(base_size = 11) +
  theme(
    plot.background  = element_rect(fill = SURF, colour = NA),
    panel.background = element_rect(fill = SURF, colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "#e6e5e1"),
    axis.title       = element_text(colour = INK2, size = 9.5),
    axis.text        = element_text(colour = INK2),
    plot.title       = element_text(colour = INK, face = "bold", size = 12),
    plot.subtitle    = element_text(colour = INK2, size = 9.5),
    legend.title     = element_blank(),
    legend.text      = element_text(colour = INK2),
    legend.position  = "top",
    plot.margin      = margin(12, 14, 10, 12)
  )

# Saves the figure twice: with titles for reading, without titles for LaTeX.
save_both <- function(p, name, w = FIGW, h = 5.2) {
  print(p)
  ggsave(file.path(FIGDIR, name), p, width = w, height = h, dpi = 150, bg = SURF)
 
  p_tex <- p
  p_tex$labels$title    <- NULL
  p_tex$labels$subtitle <- NULL
  ggsave(file.path(TEXDIR, name), p_tex, width = w, height = h, dpi = 300, bg = SURF)
 
  cat(sprintf("  saved %s (and the copy without titles in %s)\n", name, TEXDIR))
}

# --- 1. PCA

# prcomp wants samples as rows, so the matrix is transposed. Centred, not scaled, on
# all kept genes: no "top N most variable genes" cut, one fewer choice to justify.

pca <- prcomp(t(assay(vsd)))
pct <- 100 * pca$sdev^2 / sum(pca$sdev^2) # % of variance per component

cat(sprintf("variance carried by PC1-5: %s\n", paste0(round(pct[1:5], 1), "%", collapse = ", ")))

pc <- cbind(ss, pca$x[, 1:3]) # one row per sample: labels + PC1...PC3

# --- 2. what each component lines up with

# R-squared of a regression of the component on the factor. Group and batch overlap
# (HC sits mostly in Batch04), so the two values are not additive.

r2 <- function(component, factor) summary(lm(pc[[component]] ~ pc[[factor]]))$r.squared
 
assoc <- rbind(
  group_R2 = sapply(c("PC1", "PC2", "PC3"), r2, factor = "group"),
  batch_R2 = sapply(c("PC1", "PC2", "PC3"), r2, factor = "batch"))
cat("\nhow much of each component is group, and how much is batch (R-squared):\n")
print(round(assoc, 2))

# --- 3. plot the same

# both panels share the axis limir, so every sample is in the same place in each

xl <- range(pc$PC1)
yl <- range(pc$PC2)
 
pca_plot <- function(colour_by, cols) {
  ggplot(pc, aes(PC1, PC2)) +
    geom_point(aes(fill = .data[[colour_by]]),
               shape = 21, colour = SURF, stroke = 0.6, size = 2.8) +
    scale_fill_manual(values = cols) +
    coord_cartesian(xlim = xl, ylim = yl) +
    labs(x = sprintf("PC1 (%.0f%% de la varianza)", pct[1]),
         y = sprintf("PC2 (%.0f%% de la varianza)", pct[2])) +
    theme_qc
}
 
p_group <- pca_plot("group", GROUP_COL) +
  labs(title    = "PCA on vst values, coloured by group",
       subtitle = sprintf("%d genes, %d samples. Group: R-squared %.2f of PC1, %.2f of PC2.",
                          nrow(vsd), ncol(vsd), assoc["group_R2", "PC1"], assoc["group_R2", "PC2"]))
save_both(p_group, "10_pca_by_group.png")
 
p_batch <- pca_plot("batch", BATCH_COL) +
  labs(title    = "The same PCA, coloured by batch",
       subtitle = sprintf("Batch: R-squared %.2f of PC1, %.2f of PC2.",
                          assoc["batch_R2", "PC1"], assoc["batch_R2", "PC2"]))
save_both(p_batch, "11_pca_by_batch.png")
 
cat(sprintf("\nfigures written to %s and %s\n", FIGDIR, TEXDIR))
