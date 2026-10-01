# 06 -- Batch x condition: test of independence
# Three tests on the same contingency table (rows = batch, columns = condition):
#   1. Pearson chi-square   2. G-test (likelihood ratio)   3. Fisher's exact
# Each test uses ALL the selected batches at once.
# The question is "is the condition mix the same in every batch?", which needs >= 2 batches.
#
# Three runs, from wide to narrow:
#   A. all conditions  x all batches      how the whole study was run
#   B. CONDITIONS      x all batches      effect of choosing groups only
#   C. CONDITIONS      x BATCHES          the analysis subset (same choice as in 02)
#
# Run from the project root:  Rscript "functions/06_batch_independency_tests.R"

# ================== CHOICES -- edit here ==================
CONDITIONS <- c("Lung", "GBM", "HC")               # NULL = all conditions
BATCHES    <- c("Batch02", "Batch03", "Batch04")   # NULL = all batches
B          <- 1e5                                  # Monte Carlo draws, if Fisher cannot be exact
SEED       <- 1                                    # makes the Monte Carlo p-value repeatable

SHEET <- "GSE68086_sample_sheet.csv"               # written by 01

# ---- load the full sample sheet (all 285 samples, before the subsetting in 02) ----
stopifnot("sample sheet not found -- run 01 first" = file.exists(SHEET))
sheet <- read.csv(SHEET, check.names = FALSE)

# 01 writes these columns. Check them: a missing one would be NULL and data.frame()
# would drop it without a word.
stopifnot(all(c("group", "batch") %in% names(sheet)))
ss_all <- data.frame(group = sheet$group, batch = sheet$batch, stringsAsFactors = FALSE)
stopifnot(nrow(ss_all) == 285, !anyNA(ss_all))

# A typo in CONDITIONS or BATCHES would silently select nothing, so check the names.
stopifnot(
  "unknown condition in CONDITIONS" = all(CONDITIONS %in% ss_all$group),
  "unknown batch in BATCHES"        = all(BATCHES    %in% ss_all$batch)
)

print(table(ss_all$batch, ss_all$group))

# ---- helpers ----
# Keep only the chosen conditions and batches. NULL keeps everything.
select_samples <- function(ss, conditions = NULL, batches = NULL) {
  if (!is.null(conditions)) ss <- ss[ss$group %in% conditions, ]
  if (!is.null(batches))    ss <- ss[ss$batch %in% batches, ]
  ss
}

# Batch x condition table, without levels that the selection emptied.
make_table <- function(ss) table(batch = factor(ss$batch), condition = factor(ss$group))

# G-test (not in base R)
g_test <- function(tab) {
  E  <- outer(rowSums(tab), colSums(tab)) / sum(tab)   # expected counts
  G  <- 2 * sum(ifelse(tab > 0, tab * log(tab / E), 0))
  df <- (nrow(tab) - 1) * (ncol(tab) - 1)
  c(G = G, df = df, p = pchisq(G, df, lower.tail = FALSE))
}

# Fisher: exact when feasible, Monte Carlo when the table is too large
fisher_safe <- function(tab) {
  tryCatch(fisher.test(tab, workspace = 2e8),
           error = function(e) {
             set.seed(SEED)
             fisher.test(tab, simulate.p.value = TRUE, B = B)
           })
}

# ---- main function: run the three tests on one selection ----
independence_tests <- function(ss, label = "") {
  tab <- make_table(ss)
  stopifnot("need at least 2 batches and 2 conditions" = all(dim(tab) >= 2))
  n <- sum(tab)

  cat("\n==========", label, "==========\n")
  print(tab)
  cat("\nRow percentages (condition mix within each batch):\n")
  print(round(100 * prop.table(tab, 1)))

  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  E   <- chi$expected
  cat(sprintf("\nn = %d | min expected = %.2f | cells with expected < 5: %d of %d | zero cells: %d\n",
              n, min(E), sum(E < 5), length(E), sum(tab == 0)))

  # 1. chi-square
  cat(sprintf("1. Chi-square: X2 = %.2f, df = %d, p = %.3g\n",
              chi$statistic, chi$parameter, chi$p.value))
  # 2. G-test
  g <- g_test(tab)
  cat(sprintf("2. G-test:     G  = %.2f, df = %d, p = %.3g\n", g["G"], g["df"], g["p"]))
  # 3. Fisher
  fis <- fisher_safe(tab)
  cat(sprintf("3. Fisher:     p = %.3g  (%s)\n", fis$p.value, fis$method))

  # effect size + where the imbalance is
  V <- unname(sqrt(chi$statistic / (n * (min(dim(tab)) - 1))))
  cat(sprintf("\nCramer's V = %.2f  (0 = no association, 1 = batch fully determines condition)\n", V))
  cat("Adjusted residuals (|value| > 2 = cell far from independence):\n")
  print(round(chi$stdres, 1))
  if (any(E < 5)) cat("NOTE: small expected counts -> trust Fisher, not chi-square/G.\n")

  invisible(list(table = tab, chisq = chi, g = g, fisher = fis, cramer_v = V))
}

# ---- post-hoc: compare batches two at a time (Bonferroni) ----
pairwise_batches <- function(ss) {
  pairs <- combn(sort(unique(ss$batch)), 2)
  p <- apply(pairs, 2, function(b) fisher_safe(make_table(ss[ss$batch %in% b, ]))$p.value)
  data.frame(batch_1 = pairs[1, ], batch_2 = pairs[2, ],
             p = signif(p, 3), p_bonferroni = signif(pmin(p * ncol(pairs), 1), 3))
}

# ================== RUNS ==================
ss_A <- ss_all
ss_B <- select_samples(ss_all, conditions = CONDITIONS)
ss_C <- select_samples(ss_all, conditions = CONDITIONS, batches = BATCHES)

resA <- independence_tests(ss_A, label = "A. All conditions x all batches")
resB <- independence_tests(ss_B, label = "B. Chosen conditions x all batches")
resC <- independence_tests(ss_C, label = "C. Chosen conditions x chosen batches (analysis subset)")

cat("\nPost-hoc, subset C: batches compared two at a time\n")
print(pairwise_batches(ss_C))
