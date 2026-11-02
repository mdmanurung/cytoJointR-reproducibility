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

det_stats <- function(score, truth, full = TRUE) {
  if (!full) {
    return(c(auc = auc_rank(score, truth)))
  }
  roc <- roc_pts(score, truth)
  c(
    auc = auc_rank(score, truth),
    pauc10 = pauc(NULL, NULL, PAUC_FMAX, roc),
    tpr_fpr05 = tpr_at_fpr(NULL, NULL, 0.05, roc),
    tpr_fpr10 = tpr_at_fpr(NULL, NULL, 0.10, roc)
  )
}

boot_det <- function(a, truth, b = NULL, B = 2000L, seed = 1L, full = TRUE) {
  set.seed(seed)
  i1 <- which(truth)
  i0 <- which(!truth)
  est <- det_stats(a, truth, full) -
    if (is.null(b)) 0 else det_stats(b, truth, full)
  draws <- replicate(B, {
    i <- c(sample(i1, replace = TRUE), sample(i0, replace = TRUE))
    det_stats(a[i], truth[i], full) -
      if (is.null(b)) 0 else det_stats(b[i], truth[i], full)
  })
  draws <- matrix(draws, nrow = length(est))
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
    lower = unname(quantile(draws, 0.025, na.rm = TRUE)),
    upper = unname(quantile(draws, 0.975, na.rm = TRUE)),
    n_seeds = length(rows)
  )
}

# Cached comparator scores of the simulation replicates (written by 06 and 07 in `run`), one row per donor:
# scenario_id, seed, donor_id, responder_sample, true_rms_sample, compass (COMPASS_default FS, 0 when the fit
# failed, as 08), polyics (clipped net sum), polyics_mag (unclipped net sum).
sim_comparators <- function(run, scenarios, seeds) {
  reps <- tidyr::expand_grid(scen = scenarios, seed = seeds)
  read_rep <- \(scen, seed, f, sub = "") {
    p <- file.path(run, "reps", scen, sprintf("seed_%02d", seed), sub, f)
    readr::read_tsv(p, show_col_types = FALSE) |>
      dplyr::mutate(scenario_id = scen, seed = seed)
  }
  rd <- \(f, sub, m, sc, nm) {
    d <- purrr::pmap(reps, \(scen, seed) read_rep(scen, seed, f, sub)) |>
      purrr::list_rbind() |>
      dplyr::filter(method == m)
    if (!is.na(sc)) {
      d <- dplyr::filter(d, score == sc)
    }
    d$v <- if ("value" %in% names(d)) d$value else d$estimate
    d |>
      dplyr::select(scenario_id, seed, donor_id, v) |>
      dplyr::rename(!!nm := v)
  }
  truth <- purrr::pmap(reps, \(scen, seed) {
    read_rep(scen, seed, "truth_donor.tsv")
  }) |>
    purrr::list_rbind() |>
    dplyr::select(
      scenario_id,
      seed,
      donor_id,
      responder_sample,
      true_rms_sample
    )
  k <- c("scenario_id", "seed", "donor_id")
  truth |>
    dplyr::left_join(
      rd("donor_scores.tsv", "", "COMPASS_default", "FS", "compass"),
      k,
      relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      rd(
        "donor_scores.tsv",
        "polyicsflow",
        "polyICSFlow_clipped",
        "net_sum",
        "polyics"
      ),
      k,
      relationship = "one-to-one"
    ) |>
    dplyr::left_join(
      rd(
        "donor_magnitude.tsv",
        "polyicsflow",
        "polyICSFlow_unclipped",
        NA,
        "polyics_mag"
      ),
      k,
      relationship = "one-to-one"
    ) |>
    dplyr::mutate(compass = dplyr::coalesce(compass, 0))
}
