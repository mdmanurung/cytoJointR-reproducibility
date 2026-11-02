#!/usr/bin/env Rscript
# Purpose: collect the per-replicate TSVs of one pilot stage into metric tables:
#          bias, RMSE, 90% interval coverage and false sign rate per scenario x model x score,
#          fit diagnostics, replicate completion, and the COMPASS side-by-side.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/04_summarise-pilot.R <reps_run_dir> <out_dir> <n_seeds> <in_smoke|in_pilot>
# Inputs:
#   <reps_run_dir>/scenarios.tsv              [02_define-scenarios.R]
#   <reps_run_dir>/reps/*/seed_*/*.tsv        [03_run-replicate.R]
#   scripts/helper01_true-scores.R            write_session_info() [this repo]
# Outputs (in <out_dir>, new):
#   replicate_status.tsv        expected scenario x seed (stage column, n_seeds), whether it is complete
#   metrics_scores.tsv          per scenario, model, score: bias, RMSE, coverage, false sign rate, with MC SEs
#   metrics_visit_change.tsv    the same for visit change, donor level and population level
#   metrics_diagnostics.tsv     per scenario and model: worst R-hat, lowest ESS, divergences, fit seconds
#   compass_comparison.tsv      per contrast: COMPASS FS/PFS beside cytoJointR estimate and truth (RMS, poly_absolute)
#   compass_summary.tsv         per scenario: mean FS/PFS, and rank correlation with the true scores
#   session_info_summary.txt    R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})
source("scripts/helper01_true-scores.R") # write_session_info()

# Parameters -------------------------------------------------------------------
SIGN_PROB <- 0.95 # a false sign is P(change > epsilon) above this while the truth is <= 0
# Metric definitions. Per replicate, each rate is a mean over that replicate's contrasts
# (donors share hyperparameters, so contrasts are not independent). Across replicates:
#   bias, RMSE   mean error; RMSE = sqrt(mean of per-replicate MSE), SE by the delta method
#   coverage     mean of per-replicate coverage; SE = sqrt(c (1 - c) / n_replicates), binomial
#   false sign   mean over replicates that have a contrast with truth <= 0; same binomial SE

# Arguments --------------------------------------------------------------------
# Reps are read from reps_run_dir; tables go to a fresh out_dir.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 04_summarise-pilot.R <reps_run_dir> <out_dir> <n_seeds> <in_smoke|in_pilot>" = length(
    args
  ) ==
    4,
  "out_dir already exists" = !dir.exists(args[2])
)
reps_run_dir <- args[1]
out_dir <- args[2]
n_seeds <- as.integer(args[3])
stage <- args[4]
scenarios <- read_tsv(
  file.path(reps_run_dir, "scenarios.tsv"),
  show_col_types = FALSE
)

# Replicate status -------------------------------------------------------------
# Expected tasks are the stage's scenarios x seeds; only complete ones (diagnostics.tsv) are read.
replicate_status <- scenarios |>
  filter(.data[[stage]]) |>
  select(scenario_id) |>
  expand_grid(seed = seq_len(n_seeds)) |>
  mutate(
    rep_dir = file.path(
      reps_run_dir,
      "reps",
      scenario_id,
      sprintf("seed_%02d", seed)
    ),
    done = file.exists(file.path(rep_dir, "diagnostics.tsv"))
  )
read_reps <- function(file) {
  replicate_status |>
    filter(done) |>
    pull(rep_dir) |>
    file.path(file) |>
    keep(file.exists) |>
    map(\(f) {
      # donor_id is absent from some tables; coerce only where present
      read_tsv(f, show_col_types = FALSE) |>
        mutate(across(any_of("donor_id"), as.character))
    }) |>
    list_rbind()
}

# Mixing flags -----------------------------------------------------------------
# A replicate that mixes badly is flagged, not dropped: metrics keep it and report the count.
R_HAT_MAX <- 1.01 # usual thresholds, as cj_fit's own warning
ESS_MIN <- 400
mixing <- read_reps("diagnostics.tsv") |>
  mutate(
    bad_mixing = max_rhat > R_HAT_MAX | min_ess_bulk < ESS_MIN | n_divergent > 0
  )

# Metrics ----------------------------------------------------------------------
# Per-replicate rates first, then the across-replicate summary with MC SEs.
summarise_metrics <- function(x, by) {
  x |>
    left_join(
      select(mixing, scenario_id, seed, model, bad_mixing),
      join_by(scenario_id, seed, model),
      relationship = "many-to-one"
    ) |>
    summarise(
      err = mean(estimate - truth, na.rm = TRUE),
      mse = mean((estimate - truth)^2, na.rm = TRUE),
      cover = mean(covered, na.rm = TRUE),
      false_sign = if_else(
        any(truth <= 0, na.rm = TRUE),
        mean(p_positive[which(truth <= 0)] > SIGN_PROB, na.rm = TRUE),
        NA_real_
      ),
      n_unavailable = sum(is.na(covered)),
      bad_mixing = first(bad_mixing),
      .by = all_of(c(by, "seed"))
    ) |>
    summarise(
      n_reps = n(),
      n_reps_used = sum(!is.na(err)), # replicates with at least one estimate
      n_reps_no_estimate = n_reps - n_reps_used,
      bias = mean(err, na.rm = TRUE),
      bias_mcse = sd(err, na.rm = TRUE) / sqrt(n_reps_used),
      rmse = sqrt(mean(mse, na.rm = TRUE)),
      rmse_mcse = sd(mse, na.rm = TRUE) / sqrt(n_reps_used) / (2 * rmse),
      coverage = mean(cover, na.rm = TRUE),
      coverage_mcse = sqrt(coverage * (1 - coverage) / n_reps_used),
      n_reps_truth_le0 = sum(!is.na(false_sign)),
      false_sign_rate = mean(false_sign, na.rm = TRUE),
      false_sign_mcse = sqrt(
        false_sign_rate * (1 - false_sign_rate) / n_reps_truth_le0
      ),
      n_unavailable = sum(n_unavailable),
      n_reps_bad_mixing = sum(bad_mixing),
      .by = all_of(by)
    )
}
score_recovery <- read_reps("score_recovery.tsv")
stopifnot(
  "metric replicates are a subset of the done set" = nrow(anti_join(
    distinct(score_recovery, scenario_id, seed),
    filter(replicate_status, done),
    join_by(scenario_id, seed)
  )) ==
    0
)
metrics_scores <- score_recovery |>
  summarise_metrics(c("scenario_id", "model", "score"))
metrics_visit_change <- read_reps("visit_change_recovery.tsv") |>
  summarise_metrics(c("scenario_id", "model", "level", "score"))

# Diagnostics ------------------------------------------------------------------
metrics_diagnostics <- mixing |>
  summarise(
    n_reps = n(),
    max_rhat = max(max_rhat),
    min_ess_bulk = min(min_ess_bulk),
    min_ess_tail = min(min_ess_tail),
    total_divergent = sum(n_divergent),
    n_reps_divergent = sum(n_divergent > 0),
    n_reps_treedepth = sum(n_max_treedepth > 0),
    n_reps_bad_mixing = sum(bad_mixing),
    min_ebfmi = min(min_ebfmi),
    median_fit_seconds = median(fit_seconds),
    max_fit_seconds = max(fit_seconds),
    .by = c(scenario_id, model)
  )

# COMPASS ----------------------------------------------------------------------
# Same contrasts, full model, seeds with COMPASS output. FS and PFS have no truth, so
# they sit beside the true RMS and poly_absolute and are compared by rank.
compass_comparison <- read_reps("compass_scores.tsv") |>
  select(scenario_id, seed, contrast_id, compass_FS, compass_PFS) |>
  inner_join(
    score_recovery |>
      filter(model == "full", score %in% c("RMS", "poly_absolute")) |>
      pivot_wider(
        id_cols = c(scenario_id, seed, contrast_id),
        names_from = score,
        values_from = c(estimate, truth)
      ),
    join_by(scenario_id, seed, contrast_id),
    relationship = "one-to-one"
  )
compass_summary <- compass_comparison |>
  summarise(
    n_contrasts = n(),
    mean_FS = mean(compass_FS),
    mean_PFS = mean(compass_PFS),
    mean_true_RMS = mean(truth_RMS),
    mean_true_poly_absolute = mean(truth_poly_absolute),
    spearman_FS_trueRMS = suppressWarnings(cor(
      compass_FS,
      truth_RMS,
      method = "spearman"
    )),
    spearman_PFS_truePoly = suppressWarnings(cor(
      compass_PFS,
      truth_poly_absolute,
      method = "spearman"
    )),
    spearman_cj_RMS_trueRMS = suppressWarnings(cor(
      estimate_RMS,
      truth_RMS,
      method = "spearman"
    )),
    .by = scenario_id
  )

# Write outputs ----------------------------------------------------------------
dir.create(out_dir, recursive = TRUE)
write_tsv(select(replicate_status, -rep_dir), file.path(out_dir, "replicate_status.tsv"))
write_tsv(metrics_scores, file.path(out_dir, "metrics_scores.tsv"))
write_tsv(metrics_visit_change, file.path(out_dir, "metrics_visit_change.tsv"))
write_tsv(metrics_diagnostics, file.path(out_dir, "metrics_diagnostics.tsv"))
write_tsv(compass_comparison, file.path(out_dir, "compass_comparison.tsv"))
write_tsv(compass_summary, file.path(out_dir, "compass_summary.tsv"))
write_session_info(file.path(out_dir, "session_info_summary.txt"))
