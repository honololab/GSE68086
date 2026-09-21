# 03b -- the same figures again, stripped for LaTeX
# Overleaf puts the caption under the figure itself, so a title and a subtitle baked
# into the image end up saying the same thing twice, in a different font and at a size
# nobody chose. This script writes a second copy of every figure with those two lines
# removed and nothing else touched.
#
# The versions with titles are NOT replaced. They stay exactly where they were, in
# data/derived/figures, and are still the ones to look at while working. This only adds
# a second folder next to them.
#
#   data/derived/figures         with title and subtitle -- for reading
#   data/derived/figures_latex   without either          -- for \includegraphics
#
# Everything else is identical: same data, same colours, same sizes, same axis labels
# and legends, which stay because they are part of the plot rather than a caption.
#
# How it works: it runs 03 and then re-saves what 03 drew. 03 records each figure in a
# list called FIGURES as it goes, so there is no second copy of the nine plots here to
# fall out of step with the first. Add a plot to 03 and it appears here on its own.
#
# Run from the project root, like 03:  Rscript "functions/03b_plot raw data for latex.R"

source(file.path("functions", "03_plot raw data.R"))

TEXDIR <- "data/derived/figures_latex"
dir.create(TEXDIR, recursive = TRUE, showWarnings = FALSE)

stopifnot("03 recorded no figures -- has FIGURES been removed from it?" = length(FIGURES) > 0)

# 300 dpi rather than the 150 used for the on-screen copies. A figure printed across
# half a page is about 3 inches wide on paper, and 150 dpi at that size shows its
# pixels; 300 does not. The inch sizes are left alone so the text inside the figure
# comes out the same size relative to the plot as it does on screen.
TEXDPI <- 300

cat(sprintf("\nwriting %d figures without titles to %s\n", length(FIGURES), TEXDIR))

for (f in FIGURES) {
  p <- f$p
  p$labels$title    <- NULL
  p$labels$subtitle <- NULL

  out <- file.path(TEXDIR, f$name)
  ggsave(out, p, width = f$w, height = f$h, dpi = TEXDPI, bg = SURF)
  cat(sprintf("  saved %s\n", out))
}

cat(sprintf("\n%d figures written to %s\n", length(FIGURES), TEXDIR))
cat(sprintf("the versions with titles are untouched in %s\n", FIGDIR))
