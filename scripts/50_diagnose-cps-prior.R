#!/usr/bin/env Rscript
# Purpose: test whether the CPS bias of the m1dm structured mean comes from the theta0 / delta priors.
#          Refit m1dm on benchmark replicates with the interaction priors scaled up (x2, x5, x10) and
#          with learned interaction SDs (local Stan copies, package not edited). Score CPS and RMS from
#          the structured mean (mu) and from the sample draws (p) against both truths.
#          Diagnosis only: the package is not edited.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/50_diagnose-cps-prior.R runs/2026-11-03_cps-prior <scenario_id>
# Inputs:
#   cytoJointR-wt/redesign-run/                    frozen package (49278c9), loaded with pkgload
#   runs/2026-11-02_redesign/scenarios.tsv         scenario rows [05_define-benchmark.R]
#   scripts/helper01_true-scores.R, helper02_benchmark-sim.R   truth tables and simulators [this repo]
#   scripts/50_m1-hier-sd.stan, 50_m1-hier-mean.stan           m1.stan with hierarchical theta0 / delta [this repo]
# Outputs (in <run_dir>, existing dir):
#   scores_<scenario>.tsv      per variant x seed x draws x contrast x score (RMS, CPS): estimate, 90% interval, both truths
#   fits_<scenario>.tsv        per variant x seed: fit seconds, divergences, max R-hat, phi mean
#   theta_<scenario>.tsv       per variant x seed x parameter: posterior mean and SD of theta0 / delta, and the truth

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- "cytoJointR-wt/redesign-run"
SEEDS <- as.integer(strsplit(Sys.getenv("DIAG_SEEDS", "1,2,3,4,5"), ",")[[1]]) # benchmark replicate seeds; same data as reps/<scenario>/seed_NN
VARIANTS <- tribble(
  ~variant,    ~stan,                          ~theta0_sd, ~delta_sd, ~donor_effects,
  "default",   NA,                             0.5,        0.3,       TRUE,
  "x2",        NA,                             1,          0.6,       TRUE,
  "x5",        NA,                             2.5,        1.5,       TRUE,
  "x10",       NA,                             5,          3,         TRUE,
  "hier_sd",   "scripts/50_m1-hier-sd.stan",   1,          1,         TRUE, # scale of the half-normal hyperprior on the SDs
  "hier_mean", "scripts/50_m1-hier-mean.stan", 1,          1,         TRUE,
  "nodonor",   NA,                             0.5,        0.3,       FALSE # hypothesis 2: no donor loadings (cj_fit donor_effects = FALSE)
)
VARIANT_ARG <- Sys.getenv("DIAG_VARIANTS", "") # comma list to run a subset
EPSILON <- c(RMS = 1e-4, poly_absolute = 1e-4, poly_marginal = 1e-4, CPS = 0.01, CCS = 1e-4, breadth = 0.01, marker_magnitude = 1e-4)
# same epsilons as 06_run-benchmark-replicate.R; only the intervals are used here

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]
scenario_arg <- args[2]
tag <- paste0(scenario_arg, Sys.getenv("DIAG_TAG", ""))
out_scores <- file.path(run_dir, paste0("scores_", tag, ".tsv"))
stopifnot(
  "usage: Rscript 50_diagnose-cps-prior.R <run_dir> <scenario_id>" = length(args) == 2,
  "output already exists" = !file.exists(out_scores)
)
if (nzchar(VARIANT_ARG)) VARIANTS <- filter(VARIANTS, variant %in% strsplit(VARIANT_ARG, ",")[[1]])
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper01_true-scores.R")
source("scripts/helper02_benchmark-sim.R")
sc <- read_tsv("runs/2026-11-02_redesign/scenarios.tsv", show_col_types = FALSE) |>
  filter(scenario_id == scenario_arg)
stopifnot(nrow(sc) == 1)

# Stan swap --------------------------------------------------------------------
# cj_fit() calls m1_compile(); point it at a local Stan file (compiled once, cached next to the file).
m1_compile_orig <- cytoJointR:::m1_compile
use_stan <- function(path) {
  fn <- if (is.na(path)) m1_compile_orig else \() cmdstanr::cmdstan_model(path)
  assignInNamespace("m1_compile", fn, ns = "cytoJointR")
}

# One replicate ----------------------------------------------------------------
# Fit m1dm under one variant; score mu (structured) and p (sample) against both truths.
run_one <- function(seed, variant) {
  v <- VARIANTS[VARIANTS$variant == variant, ]
  use_stan(v$stan)
  data_seed <- sc$seed_offset + seed
  sim <- bench_simulate(sc, data_seed)
  tr <- true_table(sim)
  truth <- tr |>
    filter(grain == "global", score %in% c("RMS", "CPS")) |>
    select(contrast_id, score, truth_sample, truth_struct)
  fit <- cj_fit(sim$data, model = "m1dm", level = "sample", seed = data_seed, donor_effects = v$donor_effects,
                priors = list(theta0_sd = v$theta0_sd, delta_sd = v$delta_sd))
  scores <- list(sample = fit$draws_sample, structured = fit$draws_structured) |>
    imap(\(x, nm) {
      cj_scores(x, epsilon = EPSILON) |>
        filter(is.na(marker), score %in% c("RMS", "CPS")) |>
        transmute(contrast_id, score, estimate, lower, upper) |>
        mutate(draws = nm)
    }) |>
    list_rbind() |>
    left_join(truth, join_by(contrast_id, score), relationship = "many-to-one") |>
    mutate(variant = variant, seed = seed)

  dg <- fit$diagnostics
  sm <- dg$summary
  truth_theta <- bench_truth_a(sc)$theta0 # generator A truth; B uses the same theta0 and delta = 0
  th <- sm |>
    filter(grepl("^theta0|^delta|^s_theta0|^s_delta|^m_theta0", variable)) |>
    transmute(
      parameter = variable,
      mean,
      sd,
      truth = case_when(
        grepl("^theta0\\[", variable) ~ truth_theta[as.integer(sub("theta0\\[(\\d+)\\]", "\\1", variable))],
        grepl("^delta", variable) ~ 0
      )
    ) |>
    mutate(variant = variant, seed = seed)
  # stimulated-sample interaction theta0 + delta: what the stimulated double-positive patterns identify
  eng <- fit$engine_fit
  th_stim <- map(1:3, \(i) {
    d <- posterior::as_draws_matrix(eng$draws(c(sprintf("theta0[%d]", i), sprintf("delta[1,%d]", i))))
    tibble(parameter = sprintf("theta_stim[%d]", i), mean = mean(rowSums(d)), sd = sd(rowSums(d)),
           truth = truth_theta[i])
  }) |>
    list_rbind() |>
    mutate(variant = variant, seed = seed)
  list(
    scores = scores,
    fits = tibble(variant = variant, seed = seed, fit_seconds = dg$time_s, n_divergent = dg$n_divergent,
                  max_rhat = dg$max_rhat, min_ess_bulk = dg$min_ess_bulk,
                  phi_mean = sm$mean[grepl("^phi", sm$variable)][1]),
    theta = bind_rows(th, th_stim)
  )
}

# Run --------------------------------------------------------------------------
grid <- expand_grid(seed = SEEDS, variant = VARIANTS$variant)
res <- pmap(grid, run_one)
write_tsv(res |> map("scores") |> list_rbind() |> mutate(scenario_id = scenario_arg), out_scores)
write_tsv(res |> map("fits") |> list_rbind(), file.path(run_dir, paste0("fits_", tag, ".tsv")))
write_tsv(res |> map("theta") |> list_rbind(), file.path(run_dir, paste0("theta_", tag, ".tsv")))
writeLines(capture.output(sessionInfo()), file.path(run_dir, paste0("session_info_", tag, ".txt")))
