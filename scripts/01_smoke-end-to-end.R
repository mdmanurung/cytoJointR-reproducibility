#!/usr/bin/env Rscript
# Purpose: one end-to-end pass of cytoJointR on generator-A data: simulate,
#          fit M1, score, visit change, COMPASS adapter; compare scores to truth.
#          A smoke test, not calibration (one seed).
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/01_smoke-end-to-end.R runs/2026-11-02_smoke
# Inputs:
#   cytoJointR/                 package source, loaded with pkgload [this repo]
# Outputs (in <run_dir>):
#   score_recovery.tsv          per contrast and score: truth, estimate, 90% interval, covered
#   scores.tsv                  cj_scores() table
#   visit_change_donor.tsv      cj_visit_change() per donor, v1 -> v2
#   visit_change_population.tsv cj_visit_change() equal-donor average, v1 -> v2
#   timings.tsv                 wall-clock seconds per step
#   fit_diagnostics.tsv         per parameter: rhat, ess_bulk, ess_tail
#   compass_scores.tsv          native COMPASS FS/PFS per donor and run
#   session_info.txt            R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
N_DONORS <- 16 # smallest donor count rung 4 recovered at
N_CELLS <- 3000 # parent cells per sample; rung-4 scale, keeps fit ~1 min
MISSING_CONTROL <- 0.1 # fraction of donor-visits without a control sample
SEED <- 20261102
EPSILON <- c(
  RMS = 0.001,
  poly_absolute = 0.001,
  poly_marginal = 0.001,
  CPS = 0.01,
  CCS = 0.001,
  breadth = 0.01,
  marker_magnitude = 0.001
)
# smoke-test thresholds on the fraction scale, not validated

# Run directory ----------------------------------------------------------------
# A fresh directory per run; never overwrite an earlier one.
run_dir <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(
  "usage: Rscript 01_smoke-end-to-end.R <run_dir>" = !is.na(run_dir),
  "run_dir already exists" = !dir.exists(run_dir)
)
dir.create(run_dir, recursive = TRUE)
pkgload::load_all("cytoJointR", quiet = TRUE)
source("scripts/helper01_true-scores.R") # true_scores(), write_session_info()

step_seconds <- c()

# Simulate ---------------------------------------------------------------------
# Generator A: data drawn from M1 itself, so true pattern probabilities are known.
sim <- cj_simulate_counts(
  n_donors = N_DONORS,
  markers = c("IFNG", "TNF", "IL2"),
  stimulations = c("CMV", "SEB"),
  visits = c("v1", "v2"),
  N = N_CELLS,
  donor_visit = TRUE,
  missing_control = MISSING_CONTROL,
  seed = SEED
)

# Fit --------------------------------------------------------------------------
# Centered + dense metric: the model agent measured 4x faster sampling at this scale.
step_seconds["fit"] <- system.time(
  fit <- cj_fit(sim$data, seed = SEED, centered = TRUE)
)[["elapsed"]]

# True scores ------------------------------------------------------------------
# Same reference algebra the score layer uses, applied once to the true probabilities.
truth_scores <- true_scores(sim) |>
  select(contrast_id, score, truth)

# Scores and recovery ----------------------------------------------------------
step_seconds["scores"] <- system.time(
  scores <- cj_scores(fit, epsilon = EPSILON)
)[["elapsed"]]
score_recovery <- scores |>
  inner_join(
    truth_scores,
    join_by(contrast_id, score),
    relationship = "one-to-one"
  ) |>
  mutate(covered = truth >= lower & truth <= upper)

# Visit change -----------------------------------------------------------------
step_seconds["visit_change"] <- system.time(
  visit_change <- cj_visit_change(
    fit,
    epsilon = EPSILON,
    from = "v1",
    to = "v2"
  )
)[["elapsed"]]

# Diagnostics ------------------------------------------------------------------
fit_diagnostics <- fit$engine_fit$summary(
  NULL,
  "rhat",
  "ess_bulk",
  "ess_tail"
) |>
  filter(!startsWith(variable, "p[") & !startsWith(variable, "log_lik"))

# COMPASS ----------------------------------------------------------------------
# Native COMPASS on the same counts, its own default filters (reported by the adapter).
step_seconds["compass"] <- system.time(
  compass <- cj_compass(sim$data)
)[["elapsed"]]

# Write outputs ----------------------------------------------------------------
write_tsv(score_recovery, file.path(run_dir, "score_recovery.tsv"))
write_tsv(scores, file.path(run_dir, "scores.tsv"))
write_tsv(visit_change$donor, file.path(run_dir, "visit_change_donor.tsv"))
write_tsv(
  visit_change$population,
  file.path(run_dir, "visit_change_population.tsv")
)
write_tsv(
  tibble(step = names(step_seconds), seconds = step_seconds),
  file.path(run_dir, "timings.tsv")
)
write_tsv(fit_diagnostics, file.path(run_dir, "fit_diagnostics.tsv"))
write_tsv(compass$scores, file.path(run_dir, "compass_scores.tsv"))
write_session_info(file.path(run_dir, "session_info.txt"))
