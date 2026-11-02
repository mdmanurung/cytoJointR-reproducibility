#!/usr/bin/env Rscript
# Purpose: one pilot replicate: simulate one scenario at one seed, fit M1 (and the
#          donor_visit = FALSE ablation where the scenario asks), score, take visit
#          change, compare with the true scores, and run COMPASS for seeds <= COMPASS_MAX_SEED.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/03_run-replicate.R <run_dir> <scenario_id> <seed>
# Inputs:
#   <run_dir>/scenarios.tsv     scenario parameters [02_define-scenarios.R]
#   cytoJointR/                 package source, loaded with pkgload [this repo]
#   scripts/lib_generator-b.R   generator-B simulator [this repo]
#   scripts/helper01_true-scores.R  true_scores(), write_session_info() [this repo]
# Outputs (in <run_dir>/reps/<scenario_id>/seed_<seed>/):
#   score_recovery.tsv          per model, contrast and score: truth, estimate, 90% interval, P(>eps), covered
#   visit_change_recovery.tsv   per model, donor or population, score: truth, estimate, interval, covered
#   diagnostics.tsv             per model: max R-hat, min ESS, divergences, treedepth hits, min E-BFMI, fit seconds;
#                               written last, so its presence marks a complete replicate
#   compass_status.tsv          ok or failed with the error message (seeds <= COMPASS_MAX_SEED only)
#   compass_scores.tsv          COMPASS FS/PFS per contrast (seeds <= COMPASS_MAX_SEED, when COMPASS ran)
#   session_info.txt            R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
SEED_OFFSET <- 20261102L # data and fit seed = SEED_OFFSET + seed; shared across scenarios
COMPASS_MAX_SEED <- 5L # COMPASS takes minutes per replicate; seeds 1-5 only
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
MARKERS <- c("IFNG", "TNF", "IL2", "CD107a") # first J are used

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 03_run-replicate.R <run_dir> <scenario_id> <seed>" = length(
    args
  ) ==
    3
)
run_dir <- args[1]
seed <- as.integer(args[3])
sc <- read_tsv(file.path(run_dir, "scenarios.tsv"), show_col_types = FALSE) |>
  filter(scenario_id == args[2])
stopifnot("scenario_id not in scenarios.tsv" = nrow(sc) == 1)
rep_dir <- file.path(
  run_dir,
  "reps",
  sc$scenario_id,
  sprintf("seed_%02d", seed)
)
stopifnot("replicate directory already exists" = !dir.exists(rep_dir))
partial_dir <- paste0(rep_dir, ".partial") # outputs land here and move to rep_dir when complete
unlink(partial_dir, recursive = TRUE) # a stale one comes from a killed task
pkgload::load_all("cytoJointR", quiet = TRUE)
source("scripts/lib_generator-b.R")
source("scripts/helper01_true-scores.R") # true_scores(), write_session_info()

# Truth ------------------------------------------------------------------------
# M1 parameters from the scenario row; column order follows cj_design()$col_info.
j <- seq_len(sc$J)
pairs <- ncol(combn(sc$J, 2))
stim_cols <- c(
  "stim_CMV",
  "stim_SEB",
  "visit_v2",
  "visit_v2:stim_CMV",
  "visit_v2:stim_SEB"
)
gamma <- rbind(
  sc$beta_cmv + sc$beta_step * (j - 1), # stim_CMV
  sc$beta_seb + sc$beta_step * (j - 1), # stim_SEB
  0 * j, # visit_v2 main effect: none
  sc$visit_cmv + 0 * j, # visit_v2:stim_CMV
  sc$visit_seb + 0 * j
) # visit_v2:stim_SEB
dimnames(gamma) <- list(stim_cols, MARKERS[j])
omega <- matrix(sc$omega_rho, 3, 3)
diag(omega) <- 1
truth <- list(
  alpha = setNames(sc$alpha0 + 0.3 * (j - 1), MARKERS[j]),
  gamma = gamma,
  theta0 = sc$theta0_start + 0.1 * (seq_len(pairs) - 1),
  delta = rbind(CMV = rep(sc$delta_cmv, pairs), SEB = rep(sc$delta_seb, pairs)),
  tau = c(sc$tau_basal, sc$tau_induced, sc$tau_induced),
  Omega = omega,
  tau_u = sc$tau_u
)

# Simulate ---------------------------------------------------------------------
# Generator A is exact M1; generator B is the logistic-normal misspecification.
data_seed <- SEED_OFFSET + seed
sim <- if (sc$generator == "A") {
  cj_simulate_counts(
    n_donors = sc$n_donors,
    markers = MARKERS[j],
    stimulations = c("CMV", "SEB"),
    visits = c("v1", "v2"),
    N = sc$N,
    truth = truth,
    donor_visit = sc$sim_donor_visit,
    missing_control = sc$missing_control,
    seed = data_seed
  )
} else {
  simulate_generator_b(sc, truth, data_seed)
}

# True scores ------------------------------------------------------------------
# Same reference algebra the score layer uses, applied to the true probabilities.
truth_scores <- true_scores(sim)
# True visit change v1 -> v2 per donor and stimulation, and its equal-donor average.
true_change_donor <- truth_scores |>
  select(donor_id, stimulation, visit, score, truth) |>
  pivot_wider(names_from = visit, values_from = truth) |>
  drop_na(v1, v2) |>
  mutate(truth = v2 - v1, .keep = "unused")
true_change_pop <- true_change_donor |>
  summarise(truth = mean(truth), .by = c(stimulation, score))

# Fit and summarise ------------------------------------------------------------
# One pass per model. Centered + dense metric, passed explicitly.
models <- c(full = TRUE, no_donor_visit = FALSE)
if (!sc$ablation) {
  models <- models["full"]
}
fit_one <- function(donor_visit) {
  elapsed <- system.time(
    fit <- cj_fit(
      sim$data,
      seed = data_seed,
      centered = sc$centered,
      donor_visit = donor_visit
    )
  )[["elapsed"]]
  dg <- fit$diagnostics
  scores <- cj_scores(fit, epsilon = EPSILON) |>
    filter(is.na(marker)) |>
    inner_join(
      truth_scores,
      join_by(contrast_id, donor_id, visit, stimulation, score),
      relationship = "one-to-one"
    )
  vc <- cj_visit_change(fit, epsilon = EPSILON, from = "v1", to = "v2")
  change <- bind_rows(
    vc$donor |>
      filter(is.na(marker), status != "unavailable") |>
      inner_join(
        true_change_donor,
        join_by(donor_id, stimulation, score),
        relationship = "one-to-one"
      ) |>
      mutate(level = "donor"),
    vc$population |>
      filter(is.na(marker)) |>
      inner_join(
        true_change_pop,
        join_by(stimulation, score),
        relationship = "one-to-one"
      ) |>
      mutate(level = "population")
  )
  list(
    scores = scores,
    change = change,
    diagnostics = tibble(
      max_rhat = dg$max_rhat,
      min_ess_bulk = dg$min_ess_bulk,
      min_ess_tail = dg$min_ess_tail,
      n_divergent = dg$n_divergent,
      n_max_treedepth = dg$n_max_treedepth,
      min_ebfmi = min(dg$ebfmi),
      fit_seconds = elapsed
    )
  )
}
fits <- models |>
  as.list() |>
  lapply(fit_one)
collect <- \(what) {
  fits |>
    lapply(\(f) f[[what]]) |>
    list_rbind(names_to = "model") |>
    mutate(scenario_id = sc$scenario_id, seed = seed, .before = 1)
}

score_recovery <- collect("scores") |>
  select(
    scenario_id,
    seed,
    model,
    contrast_id,
    donor_id,
    visit,
    stimulation,
    score,
    estimate,
    lower,
    upper,
    epsilon,
    p_positive,
    status,
    truth
  ) |>
  mutate(covered = truth >= lower & truth <= upper)
visit_change_recovery <- collect("change") |>
  select(
    scenario_id,
    seed,
    model,
    level,
    donor_id,
    stimulation,
    score,
    estimate,
    lower,
    upper,
    epsilon,
    p_positive,
    status,
    truth
  ) |>
  mutate(covered = truth >= lower & truth <= upper)
diagnostics <- collect("diagnostics")

# Write fit outputs ------------------------------------------------------------
# The fits took minutes to hours, so they are on disk before COMPASS can fail.
dir.create(partial_dir, recursive = TRUE)
write_tsv(score_recovery, file.path(partial_dir, "score_recovery.tsv"))
write_tsv(
  visit_change_recovery,
  file.path(partial_dir, "visit_change_recovery.tsv")
)

# COMPASS ----------------------------------------------------------------------
# Native COMPASS on the same counts, its own default filters; seeds 1-5 only.
# A COMPASS error is recorded, not raised: the fit outputs stay.
compass <- if (seed <= COMPASS_MAX_SEED) {
  tryCatch(
    cj_compass(sim$data, seed = data_seed)$scores |>
      mutate(scenario_id = sc$scenario_id, seed = seed, .before = 1),
    error = \(e) e
  )
}
if (inherits(compass, "error")) {
  write_tsv(
    tibble(status = "failed", message = conditionMessage(compass)),
    file.path(partial_dir, "compass_status.tsv")
  )
}
if (is.data.frame(compass)) {
  write_tsv(tibble(status = "ok", message = NA_character_), file.path(partial_dir, "compass_status.tsv"))
  write_tsv(compass, file.path(partial_dir, "compass_scores.tsv"))
}

# Complete the replicate -------------------------------------------------------
# diagnostics.tsv is the done marker, so it goes last; then the directory takes its final name.
write_session_info(file.path(partial_dir, "session_info.txt"))
write_tsv(diagnostics, file.path(partial_dir, "diagnostics.tsv"))
file.rename(partial_dir, rep_dir)
