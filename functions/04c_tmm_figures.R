# 04c -- filterByExpr + TMM, and figures 1a, 1b, 1c and 2a of Robinson & Oshlack (2010)
#
# Robinson MD, Oshlack A. A scaling normalization method for differential expression
# analysis of RNA-seq data. Genome Biology 2010, 11:R25.
#
# Pipeline (same as the Overleaf, "Filtro por expresión"):
#   1. DGEList with the 133-sample subset written by 02
#   2. filterByExpr(group = group), then library sizes recomputed (keep.lib.sizes = FALSE)
#   3. calcNormFactors(method = "TMM")
#
# Figures (each one compares two libraries, as in the paper):
#   1a  M values between two halves of the HC samples, library-size scaling only
#   1b  M values between Lung and HC; red = TMM factor, green = housekeeping genes
#   1c  M vs A plot of the same pair as 1b
#   2a  M vs A plot of a simulated pair: genes unique to one library (orange) and
#       asymmetric DE (blue); red = TMM estimate, dashed = true value
#
# Differences with the paper, on purpose:
#   - The paper compares single deep libraries. Here each library is the sum of the
#     samples of a group (see section 2), and 1a uses two halves of HC instead of
#     technical replicates.
#   - Only the genes kept by filterByExpr are used, because that is the pipeline.
#   - Housekeeping genes: Eisenberg & Levanon (2013) list (3,804 genes), the update of the
#     2003 list the paper used (545 genes).
#
# Run from the project root, like the other scripts:  Rscript functions/04c_tmm_figures.R
# Reads the subset only. Writes PNGs to data/derived/figures_tmm (and caches the
# housekeeping list in data/derived/HK_genes.txt).

suppressPackageStartupMessages({
  library(edgeR)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

IN     <- "data/derived/subset_lung_gbm_hc.rds"            # written by 02
FIGDIR <- "data/derived/figures_tmm"
HK_URL <- "https://www.tau.ac.il/~elieis/HKG/HK_genes.txt"
HK_TXT <- "data/derived/HK_genes.txt"                       # local copy of HK_URL

stopifnot("subset not found -- run 02 first" = file.exists(IN))
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
set.seed(2010)

d      <- readRDS(IN)
counts <- d$counts
ss     <- d$ss
stopifnot(identical(colnames(counts), rownames(ss)))

# ---- 1. pipeline: filterByExpr + TMM ------------------------------------------

y    <- DGEList(counts, group = ss$group)
keep <- filterByExpr(y, group = ss$group)
y    <- y[keep, , keep.lib.sizes = FALSE]
y    <- calcNormFactors(y, method = "TMM")

cat(sprintf("genes: %d in the matrix -> %d kept by filterByExpr\n", length(keep), sum(keep)))
cat("TMM factors, median by group:\n")
print(round(tapply(y$samples$norm.factors, ss$group, median), 3))

# ---- 2. one library per group ------------------------------------------------

# The paper compares two libraries at a time (one liver vs one kidney library, each
# sequenced deeply). Single samples here have ~1-2 million reads, so their M values
# are very noisy. Summing the samples of a group gives one deep library per group,
# which is the closest analogue. For 1a, the HC samples are split at random into two
# halves: two libraries of the same population, the analogue of technical replicates.
pool <- function(cols) rowSums(y$counts[, cols, drop = FALSE])

hc_cols   <- which(ss$group == "HC")
half      <- sample(hc_cols, length(hc_cols) %/% 2)
lib_hc_a  <- pool(half)
lib_hc_b  <- pool(setdiff(hc_cols, half))
lib_hc    <- pool(hc_cols)
lib_lung  <- pool(which(ss$group == "Lung"))

# ---- 3. helper functions -------------------------------------------------------

# M value of library k against reference r, after scaling by library size only.
# Genes with a zero in either library have no M value and are left out (NA).
m_values <- function(k, r) {
  M <- log2(k / sum(k)) - log2(r / sum(r))
  M[k == 0 | r == 0] <- NA
  M
}

# log2 TMM factor of k against r: the offset that TMM removes from the M values.
tmm <- function(k, r) {
  f <- calcNormFactors(cbind(r, k), refColumn = 1)
  log2(f[2] / f[1])
}

# Histogram of M values (figures 1a, 1b). Green line: the same for housekeeping genes.
m_hist <- function(M, hk, xlab, red = NULL) {
  h  <- hist(M, breaks = 100, plot = FALSE)
  dk <- density(M[hk], na.rm = TRUE)
  plot(h, freq = FALSE, col = "grey85", border = "grey60", main = "", xlab = xlab,
       xlim = c(-3, 3), ylim = c(0, max(h$density, dk$y)))
  lines(dk, col = "darkgreen", lwd = 2)
  abline(v = 0, lty = 3)
  if (!is.null(red)) abline(v = red, col = "red", lwd = 2)
}

# M vs A plot (figures 1c, 2a), as in the paper:
#   A = average log2 proportion of the gene in the two libraries
#   M = log2 ratio of the proportions (library k over reference r)
# Genes with a zero in one library have no M or A. They are drawn in orange at the
# left edge (the "smear"), with M computed after adding 0.5 to both counts.
# `highlight` genes are drawn on top in `hl_col`.
ma_plot <- function(k, r, highlight, hl_col, ylab) {
  zero <- k == 0 | r == 0
  A <- (log2(k / sum(k)) + log2(r / sum(r))) / 2
  M <- m_values(k, r)
  A[zero] <- min(A[!zero]) - runif(sum(zero), 0.5, 1.5)
  M[zero] <- log2((k[zero] + 0.5) / sum(k)) - log2((r[zero] + 0.5) / sum(r))

  col <- ifelse(zero, "orange", "black")
  col[highlight] <- hl_col
  top <- order(col != "black")                     # coloured points drawn last
  plot(A[top], M[top], col = col[top], pch = 19, cex = 0.3, xlab = "A", ylab = ylab)
  abline(h = 0, lty = 3)
}

save_png <- function(name, expr) {
  png(file.path(FIGDIR, name), width = 6, height = 5, units = "in", res = 300)
  expr
  dev.off()
  cat("saved", file.path(FIGDIR, name), "\n")
}

# ---- 4. housekeeping genes -----------------------------------------------------

if (!file.exists(HK_TXT)) download.file(HK_URL, HK_TXT, quiet = TRUE)
hk_symbols <- trimws(read.table(HK_TXT, sep = "\t", quote = "")[, 1])
hk_ensembl <- suppressMessages(mapIds(org.Hs.eg.db, keys = hk_symbols, column = "ENSEMBL",
                                      keytype = "SYMBOL", multiVals = "first"))
is_hk <- rownames(y) %in% hk_ensembl
cat(sprintf("housekeeping genes: %d in the list -> %d among the kept genes\n",
            length(hk_symbols), sum(is_hk)))

# ---- 5. figure 1a: HC half vs HC half ------------------------------------------

M_rep <- m_values(lib_hc_b, lib_hc_a)

save_png("fig1a_hc_vs_hc.png",
         m_hist(M_rep, is_hk, xlab = "M = log2( HC half 2 / HC half 1 )"))

# ---- 6. figure 1b: Lung vs HC --------------------------------------------------

M_tis   <- m_values(lib_lung, lib_hc)
log_tmm <- tmm(lib_lung, lib_hc)
log_hk  <- median(M_tis[is_hk], na.rm = TRUE)

cat(sprintf("\nLung vs HC: TMM factor %.2f (log2 %.2f), median M of housekeeping genes %.2f\n",
            2^log_tmm, log_tmm, log_hk))
# Same offset seen from the pipeline: median 133-sample factor of Lung over that of HC
nf <- y$samples$norm.factors
cat(sprintf("            from the pipeline factors: log2(median Lung / median HC) = %.2f\n",
            log2(median(nf[ss$group == "Lung"]) / median(nf[ss$group == "HC"]))))

save_png("fig1b_lung_vs_hc.png",
         m_hist(M_tis, is_hk, xlab = "M = log2( Lung / HC )", red = log_tmm))

# ---- 7. figure 1c: M vs A, Lung vs HC ------------------------------------------

# Housekeeping genes are transparent green because there are many of them.
save_png("fig1c_ma_lung_vs_hc.png", {
  ma_plot(lib_lung, lib_hc, highlight = is_hk, hl_col = adjustcolor("darkgreen", 0.35),
          ylab = "M = log2( Lung / HC )")
  abline(h = log_hk,  col = "darkgreen", lwd = 2)
  abline(h = log_tmm, col = "red", lwd = 2)
})

# ---- 8. figure 2a: simulated pair ----------------------------------------------

# Same idea as the paper's simulation: draw gene expression levels from the observed
# counts (here the HC library), make some genes unique to library 1 and some DE, then
# draw Poisson counts. Proportions as in the paper's Figure 3: 10% unique to library 1,
# 5% DE at 2-fold, 80% of the DE genes higher in library 1.
G       <- 10000
N       <- 10e6                                    # reads per simulated library
mu      <- sample(lib_hc[lib_hc > 0], G, replace = TRUE)
unique1 <- 1:(0.10 * G)                            # expressed in library 1 only
de      <- max(unique1) + 1:(0.05 * G)
de_up1  <- de[1:(0.8 * length(de))]                # higher in library 1
de_up2  <- setdiff(de, de_up1)                     # higher in library 2

mu1 <- mu; mu1[de_up1] <- 2 * mu[de_up1]
mu2 <- mu; mu2[de_up2] <- 2 * mu[de_up2]; mu2[unique1] <- 0

sim1 <- rpois(G, N * mu1 / sum(mu1))
sim2 <- rpois(G, N * mu2 / sum(mu2))

log_sim  <- tmm(sim1, sim2)
log_true <- log2(sum(mu2) / sum(mu1))              # where the non-DE common genes sit
cat(sprintf("\nsimulation: log2 TMM = %.3f, true = %.3f\n", log_sim, log_true))

save_png("fig2a_ma_simulation.png", {
  ma_plot(sim1, sim2, highlight = de, hl_col = "blue",
          ylab = "M = log2( library 1 / library 2 )")
  abline(h = log_sim,  col = "red", lwd = 2)
  abline(h = log_true, lty = 2)
  legend("topright", bty = "n", cex = 0.8,
         legend = c("unique to library 1", "DE", "TMM estimate", "true value"),
         col = c("orange", "blue", "red", "black"), pch = c(19, 19, NA, NA), lty = c(NA, NA, 1, 2))
})
