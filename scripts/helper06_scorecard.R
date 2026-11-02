# Purpose: metric helpers of the frozen benchmark: ROC-based detection metrics (AUC, partial AUC, TPR at a
#          capped FPR) with donor-bootstrap intervals and paired differences. Shared by 69 and 70.
# Run:     source("scripts/helper06_scorecard.R") from the repo root (dplyr, purrr, tibble attached;
#          auc_rank() comes from scripts/helper03_realdata.R, sourced here).
# Inputs:  none
# Outputs:
#   roc_pts(score, truth)           tpr and fpr at every distinct score (call it once, pass it on)
#   tpr_at_fpr(score, truth, cap)   largest TPR among ROC points with FPR <= cap (empirical, ties kept together)
#   pauc(score, truth, fmax)        mean TPR over FPR in [0, fmax] (area / fmax; 0.5 * fmax is chance)
#   det_stats(score, truth)         named vector: auc, pauc10, tpr_fpr05, tpr_fpr10
#   boot_det(a, truth, b, B, seed)  det_stats of a with percentile CIs; with b, the paired difference a - b.
#                                   Donors are resampled within case and within non-case, same draws for a and b.

source("scripts/helper03_realdata.R")

PAUC_FMAX <- 0.1 # partial AUC range, fixed by benchmark/definition.tsv
PAUC_GRID <- 201L # FPR grid points for the partial-area average

roc_pts <- function(score, truth) {
  ok <- !is.na(score)
  score <- score[ok]
  truth <- truth[ok]
  thr <- sort(unique(score), decreasing = TRUE)
  tibble(
    tpr = vapply(thr, \(t) mean(score[truth] >= t), 0),
    fpr = vapply(thr, \(t) mean(score[!truth] >= t), 0)
  )
}

tpr_at_fpr <- function(score, truth, cap, roc = roc_pts(score, truth)) {
  max(0, roc$tpr[roc$fpr <= cap])
}

pauc <- function(score, truth, fmax = PAUC_FMAX, roc = roc_pts(score, truth)) {
  mean(vapply(
    seq(0, fmax, length.out = PAUC_GRID),
    \(f) tpr_at_fpr(NULL, NULL, f, roc),
    0
  ))
}

det_stats <- function(score, truth) {
  roc <- roc_pts(score, truth)
  c(
    auc = auc_rank(score, truth),
    pauc10 = pauc(NULL, NULL, PAUC_FMAX, roc),
    tpr_fpr05 = tpr_at_fpr(NULL, NULL, 0.05, roc),
    tpr_fpr10 = tpr_at_fpr(NULL, NULL, 0.10, roc)
  )
}

boot_det <- function(a, truth, b = NULL, B = 2000L, seed = 1L) {
  set.seed(seed)
  i1 <- which(truth)
  i0 <- which(!truth)
  est <- det_stats(a, truth) - if (is.null(b)) 0 else det_stats(b, truth)
  draws <- replicate(B, {
    i <- c(sample(i1, replace = TRUE), sample(i0, replace = TRUE))
    det_stats(a[i], truth[i]) - if (is.null(b)) 0 else det_stats(b[i], truth[i])
  })
  tibble(
    stat = names(est),
    value = unname(est),
    lower = apply(draws, 1, quantile, 0.025, na.rm = TRUE),
    upper = apply(draws, 1, quantile, 0.975, na.rm = TRUE),
    n_case = length(i1),
    n_noncase = length(i0)
  )
}

# Simulation metrics: the replicate (seed) is the unit, as in 08_summarise-benchmark.R.
# d: one row per donor with columns seed, truth (logical or numeric), a, and optionally b.
# stat: "auc" (truth logical) or "spearman" (truth numeric). Returns the statistic of a, or of a - b, with a
# percentile interval from resampling seeds.
seed_boot <- function(d, stat, use_b = FALSE, B = 1000L, seed = 1L) {
  f <- \(x, y) {
    if (stat == "auc") {
      auc_rank(x, y)
    } else {
      suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))
    }
  }
  g <- \(dd) f(dd$a, dd$truth) - if (use_b) f(dd$b, dd$truth) else 0
  rows <- split(seq_len(nrow(d)), d$seed)
  set.seed(seed)
  draws <- replicate(
    B,
    g(d[
      unlist(rows[sample(length(rows), replace = TRUE)], use.names = FALSE),
    ])
  )
  tibble(
    value = g(d),
    lower = quantile(draws, 0.025, na.rm = TRUE),
    upper = quantile(draws, 0.975, na.rm = TRUE),
    n_seeds = length(rows)
  )
}
