#!/usr/bin/env Rscript
# Purpose: find why C2 (m1dm) over-covers CPS on generator A. Re-fit C2 on benchmark replicates under
#          several inv_phi prior scales and score the same fit with three kinds of p draws:
#          sample (Dirichlet(phi mu + y), as benchmarked), structured (mu), and plugin (the sample draw
#          with each Gamma(a_k) replaced by its mean a_k wherever a_k = phi mu_k + y_k < 1).
#          Diagnosis only: the package is not edited.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/40_diagnose-cps.R runs/2026-11-03_diagnose-cps rare_moderate_J3
# Inputs:
#   cytoJointR-wt/redesign-run/                    frozen package (49278c9), loaded with pkgload
#   runs/2026-11-02_redesign/scenarios.tsv         scenario rows [05_define-benchmark.R]
#   scripts/helper01_true-scores.R, helper02_benchmark-sim.R   truth tables and simulators [this repo]
# Outputs (in <run_dir>, existing dir):
#   scores_<scenario>.tsv      per prior x seed x draws x contrast x score (RMS, CPS): estimate, 90% interval, both truths
#   phi_<scenario>.tsv         per prior x seed: phi quantiles (5, 50, 95%), N, fit seconds, divergences
#   patterns_<scenario>.tsv    per prior x seed x pattern degree: median a_k of the stimulated samples and share with a_k < 1

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- "cytoJointR-wt/redesign-run"
SEEDS <- as.integer(strsplit(Sys.getenv("DIAG_SEEDS", "1,2,3,4"), ",")[[1]]) # benchmark replicate seeds; same data as reps/<scenario>/seed_NN
INV_PHI_SD <- c(default = 0.05, tight = 1e-6, small = 0.01, wide = 0.2) # half-normal scale of inv_phi; tight forces phi >> N
PLUGIN_BELOW <- 1 # a_k below this is replaced by its mean in the plugin draws
EPSILON <- c(RMS = 1e-4, poly_absolute = 1e-4, poly_marginal = 1e-4, CPS = 0.01, CCS = 1e-4, breadth = 0.01, marker_magnitude = 1e-4)
# same epsilons as 06_run-benchmark-replicate.R; only the intervals are used here

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]
scenario_arg <- args[2]
out_scores <- file.path(run_dir, paste0("scores_", scenario_arg, ".tsv"))
stopifnot(
  "usage: Rscript 40_diagnose-cps.R <run_dir> <scenario_id>" = length(args) == 2,
  "output already exists" = !file.exists(out_scores)
)
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper01_true-scores.R")
source("scripts/helper02_benchmark-sim.R")
sc <- read_tsv("runs/2026-11-02_redesign/scenarios.tsv", show_col_types = FALSE) |>
  filter(scenario_id == scenario_arg)
stopifnot(nrow(sc) == 1)

# One replicate ----------------------------------------------------------------
# Fit C2 at one prior scale; score the three kinds of draws against both truths.
run_one <- function(seed, prior_name) {
  data_seed <- sc$seed_offset + seed
  sim <- bench_simulate(sc, data_seed)
  truth <- true_table(sim) |>
    filter(grain == "global", score %in% c("RMS", "CPS")) |>
    select(contrast_id, score, truth_sample, truth_struct)
  fit <- cj_fit(sim$data, model = "m1dm", level = "sample", seed = data_seed,
                priors = list(inv_phi_sd = INV_PHI_SD[[prior_name]]))
  eng <- fit$engine_fit
  phi <- as.vector(posterior::as_draws_matrix(eng$draws("phi")))
  p_samp <- fit$draws_sample
  mu <- fit$draws_structured$p # [draw, sample, pattern]
  y <- sim$data$counts # [sample, pattern]

  # plugin draws: Gamma(a_k) -> a_k where a_k < PLUGIN_BELOW, then renormalise. a = phi mu + y per draw.
  nd <- dim(mu)[1]
  a <- mu * phi + array(rep(y, each = nd), dim(mu)) # phi recycles along the draw axis
  g <- array(rgamma(length(a), shape = a), dim(a))
  g[a < PLUGIN_BELOW] <- a[a < PLUGIN_BELOW]
  p_plug <- p_samp
  p_plug$p <- g / array(rep(apply(g, 1:2, sum), dim(g)[3]), dim(g))

  draws <- list(sample = p_samp, structured = fit$draws_structured, plugin = p_plug)
  scores <- imap(draws, \(x, nm) {
    cj_scores(x, epsilon = EPSILON) |>
      filter(is.na(marker), score %in% c("RMS", "CPS")) |>
      transmute(contrast_id, score, estimate, lower, upper) |>
      mutate(draws = nm)
  }) |>
    list_rbind() |>
    left_join(truth, join_by(contrast_id, score), relationship = "many-to-one") |>
    mutate(prior = prior_name, seed = seed)

  # a_k of the stimulated samples by pattern degree, at the posterior median phi
  stim <- sim$data$samples$condition != "control"
  deg <- rowSums(sim$data$bits)
  a_med <- apply(mu[, stim, , drop = FALSE], c(2, 3), mean) * median(phi) + y[stim, , drop = FALSE]
  pats <- tibble(degree = rep(deg, each = sum(stim)), a = as.vector(a_med)) |>
    summarise(median_a = median(a), share_below_1 = mean(a < 1), .by = degree) |>
    mutate(prior = prior_name, seed = seed)

  dg <- fit$diagnostics
  list(
    scores = scores,
    phi = tibble(prior = prior_name, seed = seed, N = sc$N,
                 phi_q05 = quantile(phi, 0.05), phi_q50 = median(phi), phi_q95 = quantile(phi, 0.95),
                 n_divergent = dg$n_divergent, max_rhat = dg$max_rhat, fit_seconds = dg$time_s),
    patterns = pats
  )
}

# Run --------------------------------------------------------------------------
grid <- expand_grid(seed = SEEDS, prior_name = names(INV_PHI_SD))
res <- pmap(grid, run_one)
write_tsv(res |> map("scores") |> list_rbind() |> mutate(scenario_id = scenario_arg), out_scores)
write_tsv(res |> map("phi") |> list_rbind(), file.path(run_dir, paste0("phi_", scenario_arg, ".tsv")))
write_tsv(res |> map("patterns") |> list_rbind(), file.path(run_dir, paste0("patterns_", scenario_arg, ".tsv")))
writeLines(capture.output(sessionInfo()), file.path(run_dir, paste0("session_info_", scenario_arg, ".txt")))
