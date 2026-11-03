#!/usr/bin/env Rscript
# Purpose: two follow-up checks for the CCS diagnosis (45_diagnose-ccs.R).
#          (A) Power: every benchmark scenario has a structured CCS truth of exactly 0 (delta = 0), so "coverage 1.0"
#              only says the interval holds 0. Simulate generator A with a real change in baseline dependence
#              (delta > 0) and ask whether C1 and C0 recover it.
#          (B) C2 emulator: draw p ~ Dirichlet(phi * mu + y) with the oracle mu = p_struct and the stored phi, the
#              residual model of C2, and compare its CCS with the stored C2 CCS under generator B.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/46_diagnose-ccs-power.R runs/2026-11-03_diagnose-ccs
# Inputs:
#   cytoJointR-wt/redesign-run/                 package source, frozen at 49278c9 [this repo]
#   runs/2026-11-02_redesign/scenarios.tsv      scenario rows [05_define-benchmark.R]
#   runs/2026-11-02_redesign/reps/generator_B_J3/seed_NN/{model_scores,diagnostics}.tsv   [06_run-benchmark-replicate.R]
# Outputs (in <run_dir>, which must already hold the outputs of 45_diagnose-ccs.R):
#   power_delta.tsv             per delta, seed, contrast and model: truth, estimate, interval
#   power_summary.tsv           per delta and model: mean truth, mean estimate, coverage, share of contrasts with lower > 0
#   c2_emulator.tsv             per contrast: emulated and stored C2 CCS against both truths
#   c2_emulator_summary.tsv     bias, coverage and RMSE of the emulated and the stored C2
#   session_info_46.txt         R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- "cytoJointR-wt/redesign-run" # frozen redesign package
REDESIGN_DIR <- "runs/2026-11-02_redesign" # stored benchmark replicates
POWER_SCENARIO <- "moderate_J3" # generator A, 24 donors, N = 10000, rms 0.01
DELTA_GRID <- c(0, 0.5, 1) # change in each pairwise log odds ratio on stimulation (prior SD is 0.3)
POWER_SEEDS <- 1:3 # one C1 fit takes about 20 s with 4 cores
EMULATOR_SCENARIO <- "generator_B_J3"
EMULATOR_SEEDS <- 1:10 # stored replicates to compare with
N_DRAWS <- 1000 # posterior draws per sample
EPS_CCS <- 1e-4 # CCS epsilon used by the benchmark
LEVEL <- 0.9 # interval level used by the benchmark
EPSILON <- c(
  RMS = 1e-4,
  poly_absolute = 1e-4,
  poly_marginal = 1e-4,
  CPS = 0.01,
  CCS = EPS_CCS,
  breadth = 0.01,
  marker_magnitude = 1e-4
)

# Setup ------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(
  "usage: Rscript 46_diagnose-ccs-power.R <run_dir>" = !is.na(run_dir),
  "run_dir must exist (run 45_diagnose-ccs.R first)" = dir.exists(run_dir),
  "power_delta.tsv already exists" = !file.exists(file.path(
    run_dir,
    "power_delta.tsv"
  ))
)
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper01_true-scores.R") # true_table(), write_session_info()
source("scripts/helper02_benchmark-sim.R") # bench_truth_a(), bench_simulate()
scenarios <- read_tsv(
  file.path(REDESIGN_DIR, "scenarios.tsv"),
  show_col_types = FALSE
)
score_contrast <- cytoJointR:::score_contrast

# (A) Power with a real CCS ----------------------------------------------------
# bench_simulate()'s generator-A branch, with delta set to DELTA in every pairwise term.
power_sc <- filter(scenarios, scenario_id == POWER_SCENARIO)
simulate_delta <- function(delta, seed) {
  truth <- bench_truth_a(power_sc)
  truth$delta <- matrix(delta, 1, choose(power_sc$J, 2))
  cj_simulate_counts(
    n_donors = power_sc$n_donors,
    markers = c("IFNG", "TNF", "IL2", "CD107a")[seq_len(power_sc$J)],
    stimulations = "CMV",
    visits = "v1",
    N = power_sc$N,
    truth = truth,
    donor_visit = FALSE,
    missing_control = power_sc$missing_control,
    seed = power_sc$seed_offset + seed
  )
}
power_delta <- expand_grid(delta = DELTA_GRID, seed = POWER_SEEDS) |>
  mutate(
    table = map2(delta, seed, \(delta, seed) {
      sim <- simulate_delta(delta, seed)
      truth <- true_table(sim) |>
        filter(grain == "global", score == "CCS") |>
        select(contrast_id, truth_sample)
      fits <- list(
        C0 = cj_fit_dirichlet(
          sim$data,
          a0 = 0.5,
          n_draws = N_DRAWS,
          seed = seed
        ),
        C1 = cj_fit(sim$data, model = "m1", level = "sample", seed = seed)
      )
      fits |>
        map(\(x) {
          cj_scores(x, EPSILON, level = LEVEL) |>
            filter(score == "CCS") |>
            select(contrast_id, estimate, lower, upper, p_positive, p_negative)
        }) |>
        list_rbind(names_to = "model") |>
        left_join(truth, join_by(contrast_id), relationship = "many-to-one")
    })
  ) |>
  unnest(table)
power_summary <- power_delta |>
  summarise(
    mean_truth = mean(truth_sample),
    mean_estimate = mean(estimate),
    bias = mean(estimate - truth_sample),
    coverage = mean(truth_sample >= lower & truth_sample <= upper),
    share_lower_above_zero = mean(lower > 0),
    share_p_positive_above_95 = mean(p_positive > 0.95),
    .by = c(delta, model)
  )

# (B) C2 emulator --------------------------------------------------------------
# p_U ~ Dirichlet(phi * mu_U + y_U) and p_S ~ Dirichlet(phi * mu_S + y_S), mu the oracle structured probability.
draw_dirichlet <- function(shape, n) {
  g <- pmax(
    matrix(rgamma(n * length(shape), shape = rep(shape, each = n)), n),
    1e-300
  )
  g / rowSums(g)
}
emulator_sc <- filter(scenarios, scenario_id == EMULATOR_SCENARIO)
c2_emulator <- tibble(seed = EMULATOR_SEEDS) |>
  mutate(
    table = map(seed, \(seed) {
      rep_dir <- file.path(
        REDESIGN_DIR,
        "reps",
        EMULATOR_SCENARIO,
        sprintf("seed_%02d", seed)
      )
      sim <- bench_simulate(emulator_sc, emulator_sc$seed_offset + seed)
      bits <- sim$data$bits
      phi <- read_tsv(
        file.path(rep_dir, "diagnostics.tsv"),
        show_col_types = FALSE
      ) |>
        filter(model == "C2") |>
        pull(phi_mean)
      stored <- read_tsv(
        file.path(rep_dir, "model_scores.tsv"),
        show_col_types = FALSE
      ) |>
        filter(
          model == "C2",
          draws == "sample",
          grain == "global",
          score == "CCS"
        ) |>
        select(
          contrast_id,
          stored_estimate = estimate,
          stored_lower = lower,
          stored_upper = upper,
          truth_sample,
          truth_struct
        )
      set.seed(seed)
      ct <- filter(sim$data$contrasts, !is.na(control_sample_id))
      emulated <- pmap(
        list(ct$contrast_id, ct$control_sample_id, ct$stim_sample_id),
        \(id, u, s) {
          ccs <- score_contrast(
            draw_dirichlet(
              phi * sim$p_struct[u, ] + sim$data$counts[u, ],
              N_DRAWS
            ),
            draw_dirichlet(
              phi * sim$p_struct[s, ] + sim$data$counts[s, ],
              N_DRAWS
            ),
            bits
          )$G[, "CCS"]
          tibble(
            contrast_id = id,
            emulated_estimate = mean(ccs, na.rm = TRUE),
            emulated_lower = quantile(
              ccs,
              (1 - LEVEL) / 2,
              na.rm = TRUE,
              names = FALSE
            ),
            emulated_upper = quantile(
              ccs,
              (1 + LEVEL) / 2,
              na.rm = TRUE,
              names = FALSE
            ),
            phi = phi
          )
        }
      ) |>
        list_rbind()
      left_join(
        stored,
        emulated,
        join_by(contrast_id),
        relationship = "one-to-one"
      )
    })
  ) |>
  unnest(table)
c2_emulator_summary <- c2_emulator |>
  pivot_longer(
    c(stored_estimate, emulated_estimate),
    names_to = "source",
    names_pattern = "(.*)_estimate",
    values_to = "estimate"
  ) |>
  mutate(
    lower = if_else(source == "stored", stored_lower, emulated_lower),
    upper = if_else(source == "stored", stored_upper, emulated_upper)
  ) |>
  pivot_longer(
    c(truth_sample, truth_struct),
    names_to = "truth",
    values_to = "value"
  ) |>
  summarise(
    mean_estimate = mean(estimate),
    bias = mean(estimate - value),
    rmse = sqrt(mean((estimate - value)^2)),
    coverage = mean(value >= lower & value <= upper),
    mean_phi = mean(phi),
    .by = c(source, truth)
  )

# Write outputs ----------------------------------------------------------------
write_tsv(power_delta, file.path(run_dir, "power_delta.tsv"))
write_tsv(power_summary, file.path(run_dir, "power_summary.tsv"))
write_tsv(c2_emulator, file.path(run_dir, "c2_emulator.tsv"))
write_tsv(c2_emulator_summary, file.path(run_dir, "c2_emulator_summary.tsv"))
write_session_info(file.path(run_dir, "session_info_46.txt"), PKG_DIR)
