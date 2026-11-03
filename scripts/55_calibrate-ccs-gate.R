#!/usr/bin/env Rscript
# Purpose: calibrate the default thresholds of the CCS support gate (cytoJointR, branch feat/ccs-gate). The gate
#          withholds CCS unless (1) the control and the stimulated sample each hold at least CCS_MIN_POLY_CELLS
#          observed polyfunctional cells (degree >= 2) and (2) the posterior CV of the control's polyfunctional
#          mass is at most CCS_MAX_CV. Every contrast is scored ungated, with the gate quantities beside it, so
#          one set of fits serves every threshold in the grid. A second proxy for condition 2, the IPF gain
#          (counterfactual / control polyfunctional mass), is recorded for comparison.
# Run:     R=<R4_51>/bin/Rscript; export PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib
#          $R scripts/55_calibrate-ccs-gate.R <run_dir> setup
#          $R scripts/55_calibrate-ccs-gate.R <run_dir> fit <regime_id> <seed>     (one fit, about 1 min on 4 cores)
#          $R scripts/55_calibrate-ccs-gate.R <run_dir> summarise
# Inputs:
#   runs/2026-11-02_redesign/scenarios.tsv    sparse scenarios and the effect sizes [05_define-benchmark.R]
#   <PKG_DIR>/                                package source, loaded with pkgload; env CJ_PKG_DIR overrides it [worktree]
#   scripts/helper01_true-scores.R, helper02_benchmark-sim.R   true scores, simulators [this repo]
# Outputs (in <run_dir>):
#   regimes.tsv                 one row per regime: generator, N, delta, solved alpha, beta [setup]
#   reps/<regime>/seed_NN.tsv   one row per contrast: gate quantities, CCS summary, both truths [fit]
#   contrasts.tsv               all replicates stacked [summarise]
#   gate_off.tsv                per regime: coverage and false sign with no gate
#   gate_grid.tsv               per regime x CCS_MIN_POLY_CELLS x CCS_MAX_CV: pass rate, gated coverage, false sign, SEs
#   gate_grid_gain.tsv          same with the IPF gain as condition 2
#   gate_chosen.tsv             per regime at the chosen thresholds, with the calibration verdict
#   session_info.txt            R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/ccs-gate") # package source under test
REDESIGN_DIR <- "runs/2026-11-02_redesign" # sparse scenarios and their solved effects
DENSE_BASELINE <- 0.05 # median control any-positive frequency of the dense regimes
DENSE_N <- c(1e4, 1e5) # parent cells per sample in the dense regimes
DENSE_DELTA <- 0.5 # change in each pairwise log odds ratio on stimulation, in the one dense cell with a real CCS
DENSE_BETA_FROM <- "rare_moderate_J3" # dense regimes reuse this scenario's induced effect
SPARSE_FROM <- c("moderate_J3", "rare_moderate_J3", "generator_B_J3") # existing regimes (generators A, A, B)
TUNE_DONORS <- 240L # donors in the cohort that solves alpha
TUNE_SEED <- 20261103L # seed of that cohort
SEED_OFFSET <- 20261200L # replicate data seed = SEED_OFFSET + 1000 * regime index + seed
EPSILON <- c(
  RMS = 1e-4,
  poly_absolute = 1e-4,
  poly_marginal = 1e-4,
  CPS = 0.01,
  CCS = 1e-4,
  breadth = 0.01,
  marker_magnitude = 1e-4
) # as the benchmark (06)
LEVEL <- 0.9 # interval level of the benchmark
SIGN_PROB <- 0.95 # a sign is claimed when P(CCS > eps) or P(CCS < -eps) exceeds this, as the benchmark
COVERAGE_BAND <- c(0.85, 0.97) # the gated CCS must cover in this band
FALSE_SIGN_MAX <- 0.07 # and claim a wrong or unsupported sign at most this often
MIN_PASSING <- 24 # a regime enters the verdict only if at least this many contrasts pass (one cohort's worth)
GRID_CELLS <- c(1, 5, 10, 20, 50, 100, 200) # CCS_MIN_POLY_CELLS candidates
GRID_CV <- c(0.1, 0.2, 0.3, 0.5, 1, Inf) # CCS_MAX_CV candidates (Inf: condition 2 off)
GRID_GAIN <- c(1.5, 2, 3, 5, 10, Inf) # largest IPF gain candidates
CHOSEN_CELLS <- 20 # chosen CCS_MIN_POLY_CELLS, set after reading gate_grid.tsv
CHOSEN_CV <- 0.5 # chosen CCS_MAX_CV, set after reading gate_grid.tsv

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]
stage <- args[2]
stopifnot(
  "usage: Rscript 55_calibrate-ccs-gate.R <run_dir> setup|fit|summarise" = stage %in%
    c("setup", "fit", "summarise")
)
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper01_true-scores.R") # true_table(), write_session_info()
source("scripts/helper02_benchmark-sim.R") # bench_simulate(), bench_truth_a()

# Simulator --------------------------------------------------------------------
# Generator A with a pairwise change delta (bench_simulate() fixes delta at 0), or generator B through bench_simulate().
simulate_regime <- function(sc, seed, n_donors = sc$n_donors) {
  if (sc$generator == "B") {
    return(bench_simulate(sc, seed, n_donors))
  }
  truth <- bench_truth_a(sc)
  truth$delta <- matrix(sc$delta, 1, choose(sc$J, 2))
  sim <- cj_simulate_counts(
    n_donors = n_donors,
    markers = c("IFNG", "TNF", "IL2", "CD107a")[seq_len(sc$J)],
    stimulations = "CMV",
    visits = "v1",
    N = sc$N,
    truth = truth,
    donor_visit = FALSE,
    seed = seed
  )
  sim$p_struct <- sim$p_struct %||% sim$mu %||% sim$p
  sim
}

# Setup ------------------------------------------------------------------------
if (stage == "setup") {
  stopifnot("run_dir already exists" = !dir.exists(run_dir))
  scenarios <- read_tsv(
    file.path(REDESIGN_DIR, "scenarios.tsv"),
    show_col_types = FALSE
  ) |>
    mutate(delta = 0)
  sparse <- scenarios |>
    filter(scenario_id %in% SPARSE_FROM) |>
    mutate(regime_id = paste0(scenario_id, "_sparse"))
  beta_from <- filter(scenarios, scenario_id == DENSE_BETA_FROM)
  # Dense regimes: generators A and B at each N with the effect of DENSE_BETA_FROM, plus A at the largest N with a
  # real CCS (delta). Alpha is solved below.
  dense <- scenarios |>
    filter(scenario_id %in% c("moderate_J3", "generator_B_J3")) |>
    mutate(across(starts_with("beta_"), \(x) 0)) |>
    select(-N) |>
    cross_join(tibble(N = DENSE_N, delta_cell = FALSE)) |>
    bind_rows(
      scenarios |>
        filter(scenario_id == "moderate_J3") |>
        select(-N) |>
        cross_join(tibble(N = max(DENSE_N), delta_cell = TRUE))
    ) |>
    mutate(
      across(starts_with("beta_"), \(x) 0),
      delta = if_else(delta_cell, DENSE_DELTA, 0),
      regime_id = sprintf(
        "%s_dense_N%s%s",
        generator,
        format(N, scientific = TRUE),
        if_else(delta_cell, "_delta", "")
      )
    )
  dense[paste0("beta_", 1:4)] <- beta_from[
    rep(1, nrow(dense)),
    paste0("beta_", 1:4)
  ]
  baseline <- \(sim) {
    median(rowSums(sim$p[
      na.omit(sim$data$contrasts$control_sample_id),
      rowSums(sim$data$bits) >= 1
    ]))
  }
  solved <- dense |>
    group_split(row_number()) |>
    map(\(sc) {
      alpha <- uniroot(
        \(x) {
          sc[paste0("alpha_", 1:4)] <- as.list(x + 0.3 * 0:3)
          baseline(simulate_regime(sc, TUNE_SEED, TUNE_DONORS)) - DENSE_BASELINE
        },
        c(-8, 0),
        tol = 1e-3
      )$root
      sc[paste0("alpha_", 1:4)] <- as.list(alpha + 0.3 * 0:3)
      sc
    }) |>
    list_rbind()
  regimes <- bind_rows(sparse, solved) |>
    mutate(
      regime_index = row_number(),
      seed_offset = SEED_OFFSET + 1000L * regime_index
    ) |>
    select(
      regime_id,
      regime_index,
      generator,
      J,
      n_donors,
      N,
      delta,
      theta0_start,
      tau_basal,
      tau_induced,
      omega_rho,
      starts_with("gb_"),
      starts_with("alpha_"),
      starts_with("beta_"),
      seed_offset,
      missing_control,
      n_visits,
      visit_effect,
      tau_u,
      resp_frac,
      eps_pattern
    )
  dir.create(run_dir, recursive = TRUE)
  write_tsv(regimes, file.path(run_dir, "regimes.tsv"))
  write_session_info(file.path(run_dir, "session_info.txt"), PKG_DIR)
  quit(save = "no")
}

# Fit one replicate ------------------------------------------------------------
if (stage == "fit") {
  regime <- filter(
    read_tsv(file.path(run_dir, "regimes.tsv"), show_col_types = FALSE),
    regime_id == args[3]
  )
  seed <- as.integer(args[4])
  out_file <- file.path(
    run_dir,
    "reps",
    regime$regime_id,
    sprintf("seed_%02d.tsv", seed)
  )
  stopifnot(
    "one regime_id expected" = nrow(regime) == 1,
    "replicate already exists" = !file.exists(out_file)
  )
  sim <- simulate_regime(regime, regime$seed_offset + seed)
  fit <- cj_fit(sim$data, level = "sample", seed = regime$seed_offset + seed) # package default model
  d <- fit$draws
  ct <- d$contrasts |> filter(!is.na(control_sample_id))
  poly <- rowSums(d$bits) >= 2
  # Observed polyfunctional cells per sample, and the control's posterior polyfunctional mass.
  poly_cells <- rowSums(fit$data$counts[, poly, drop = FALSE])
  eng <- cytoJointR:::score_engine(d)
  gate_q <- ct |>
    mutate(
      control_poly = poly_cells[control_sample_id],
      stim_poly = poly_cells[stim_sample_id],
      control_mass = map_dbl(control_sample_id, \(u) {
        mean(rowSums(d$p[, u, poly, drop = FALSE]))
      }),
      cv_control_mass = map_dbl(control_sample_id, \(u) {
        m <- rowSums(d$p[, u, poly, drop = FALSE])
        sd(m) / mean(m)
      }),
      # IPF gain: counterfactual polyfunctional mass over control polyfunctional mass, posterior median over draws.
      # CCS = Ps_poly - Q_poly in every draw, so Q_poly = Ps_poly - CCS.
      ipf_gain = pmap_dbl(
        list(contrast_id, control_sample_id, stim_sample_id),
        \(id, u, s) {
          ps_poly <- rowSums(d$p[, s, poly, drop = FALSE])
          median(
            (ps_poly - eng[[id]]$G[, "CCS"]) /
              rowSums(d$p[, u, poly, drop = FALSE]),
            na.rm = TRUE
          )
        }
      )
    ) |>
    select(
      contrast_id,
      donor_id,
      control_poly,
      stim_poly,
      control_mass,
      cv_control_mass,
      ipf_gain
    )
  gate_off <- if ("ccs_gate" %in% names(formals(cj_scores))) {
    list(ccs_gate = NULL)
  } else {
    list()
  }
  ccs <- do.call(cj_scores, c(list(fit, EPSILON, level = LEVEL), gate_off)) |>
    filter(score == "CCS") |>
    select(
      contrast_id,
      estimate,
      lower,
      upper,
      p_positive,
      p_negative,
      failure_fraction,
      status
    )
  truth <- true_table(sim) |>
    filter(grain == "global", score == "CCS") |>
    select(contrast_id, truth_sample, truth_struct)
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  gate_q |>
    left_join(ccs, join_by(contrast_id), relationship = "one-to-one") |>
    left_join(truth, join_by(contrast_id), relationship = "one-to-one") |>
    mutate(regime_id = regime$regime_id, seed = seed, .before = 1) |>
    write_tsv(out_file)
  quit(save = "no")
}

# Summarise --------------------------------------------------------------------
contrasts <- list.files(
  file.path(run_dir, "reps"),
  pattern = "\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
) |>
  map(\(f) read_tsv(f, show_col_types = FALSE)) |>
  list_rbind() |>
  mutate(
    covered_sample = truth_sample >= lower & truth_sample <= upper,
    covered_struct = truth_struct >= lower & truth_struct <= upper,
    claim_positive = p_positive > SIGN_PROB,
    claim_negative = p_negative > SIGN_PROB,
    # a claimed sign is false when the truth does not exceed epsilon in that direction
    false_sign = (claim_positive & truth_sample <= EPSILON[["CCS"]]) |
      (claim_negative & truth_sample >= -EPSILON[["CCS"]])
  )
write_tsv(contrasts, file.path(run_dir, "contrasts.tsv"))
n_seed <- n_distinct(contrasts$seed)

# One row per regime for the contrasts in `keep`. SE is clustered by seed (donors share a fit), over seeds that
# have at least one passing contrast; pass_rate is over all contrasts.
summarise_gate <- function(data) {
  per_seed <- data |>
    summarise(
      n = n(),
      k = sum(keep),
      cov = if_else(k > 0, mean(covered_sample[keep], na.rm = TRUE), NA_real_),
      cov_struct = if_else(
        k > 0,
        mean(covered_struct[keep], na.rm = TRUE),
        NA_real_
      ),
      fs = if_else(k > 0, mean(false_sign[keep], na.rm = TRUE), NA_real_),
      .by = c(regime_id, seed)
    )
  per_seed |>
    summarise(
      n_contrasts = sum(n),
      n_pass = sum(k),
      pass_rate = n_pass / n_contrasts,
      pass_rate_se = sd(k / n) / sqrt(n()),
      coverage = weighted.mean(cov, k, na.rm = TRUE),
      coverage_se = sd(cov, na.rm = TRUE) / sqrt(sum(!is.na(cov))),
      coverage_struct = weighted.mean(cov_struct, k, na.rm = TRUE),
      false_sign = weighted.mean(fs, k, na.rm = TRUE),
      false_sign_se = sd(fs, na.rm = TRUE) / sqrt(sum(!is.na(fs))),
      .by = regime_id
    )
}
by_gate <- \(data, cells, cv, gain) {
  keep <- coalesce(
    data$control_poly >= cells &
      data$stim_poly >= cells &
      data$cv_control_mass <= cv &
      data$ipf_gain <= gain,
    FALSE
  )
  summarise_gate(mutate(data, keep = keep)) |>
    mutate(min_cells = cells, max_cv = cv, max_gain = gain)
}
gate_off <- by_gate(contrasts, 0, Inf, Inf) |>
  select(-min_cells, -max_cv, -max_gain)
gate_grid <- expand_grid(min_cells = GRID_CELLS, max_cv = GRID_CV) |>
  pmap(\(min_cells, max_cv) by_gate(contrasts, min_cells, max_cv, Inf)) |>
  list_rbind()
gate_grid_gain <- expand_grid(min_cells = GRID_CELLS, max_gain = GRID_GAIN) |>
  pmap(\(min_cells, max_gain) by_gate(contrasts, min_cells, Inf, max_gain)) |>
  list_rbind()
# Verdict: a regime with fewer than MIN_PASSING passing contrasts is not judged. Coverage above the band is flagged
# (over-covers: conservative, not harmful); coverage below it, or too many false signs, fails.
verdict <- \(g) {
  g |>
    mutate(
      judged = n_pass >= MIN_PASSING,
      over_covers = judged & coverage > COVERAGE_BAND[2],
      calibrated = !judged |
        (coverage >= COVERAGE_BAND[1] & false_sign <= FALSE_SIGN_MAX)
    )
}
gate_chosen <- by_gate(contrasts, CHOSEN_CELLS, CHOSEN_CV, Inf) |> verdict()
write_tsv(gate_off, file.path(run_dir, "gate_off.tsv"))
write_tsv(gate_grid, file.path(run_dir, "gate_grid.tsv"))
write_tsv(gate_grid_gain, file.path(run_dir, "gate_grid_gain.tsv"))
write_tsv(gate_chosen, file.path(run_dir, "gate_chosen.tsv"))
write_session_info(file.path(run_dir, "session_info.txt"), PKG_DIR)
