#!/usr/bin/env Rscript
# Purpose: simulation-based calibration of cj_fit(model = "m1dm"): draw the parameters
#          (alpha, gamma, theta0, delta, tau, Omega, inv_phi) from the priors, simulate counts
#          with cj_simulate_counts(truth = list(phi = 1 / inv_phi)), refit, and record the rank
#          of each true parameter among the posterior draws. Ranks are uniform when the
#          sampler and the model code agree. The 100-simulation run is a cluster job.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/10_sbc-m1dm.R <run_dir> [sim_from=1] [sim_to=100] [pkg_dir=cytoJointR]
#          Chunks share one run_dir: sbatch --array tasks with disjoint sim_from..sim_to.
# Inputs:
#   <pkg_dir>/                  package source, loaded with pkgload [this repo]
#   scripts/helper01_true-scores.R  write_session_info() [this repo]
# Outputs (in <run_dir>/sims/):
#   sim_<NNN>.tsv               per parameter: truth, rank among the thinned draws, number of draws,
#                               R-hat and bulk ESS; divergences and fit seconds repeat on every row
# Outputs (in <run_dir>/):
#   session_info_<sim_from>.txt R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

# Parameters -------------------------------------------------------------------
SEED_OFFSET <- 20261102L # prior-draw seed = SEED_OFFSET + simulation number
DATA_SEED_OFFSET <- 30261102L # data and fit seed; must differ from SEED_OFFSET: cj_simulate_counts() reseeds, and a shared seed ties the donor effects to alpha
N_DONORS <- 8L # small on purpose: a hundred refits must fit a day on the cluster
MARKERS <- c("IFNG", "TNF") # J = 2: 4 patterns
N_CELLS <- 2000L # per sample; with the prior on inv_phi this spans near-multinomial to strong overdispersion
CHAINS <- 4L
ITER <- 1000L # warmup and sampling each
THIN <- 4L # ranks use every THIN-th draw, so autocorrelation does not inflate rank counts
SBC_PARAMS <- "^(alpha|gamma|inv_phi)" # Stan names to rank; phi is 1 / inv_phi, so its ranks are reversed

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 10_sbc-m1dm.R <run_dir> [sim_from] [sim_to] [pkg_dir]" = length(
    args
  ) %in%
    1:4
)
run_dir <- args[1]
sim_from <- if (length(args) >= 2) as.integer(args[2]) else 1L
sim_to <- if (length(args) >= 3) as.integer(args[3]) else 100L
pkg_dir <- if (length(args) >= 4) args[4] else "cytoJointR"
dir.create(file.path(run_dir, "sims"), recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(pkg_dir, quiet = TRUE)
source("scripts/helper01_true-scores.R") # write_session_info()

# Design -----------------------------------------------------------------------
# One throwaway simulation gives the design the priors are drawn for.
sim0 <- cj_simulate_counts(
  n_donors = N_DONORS,
  markers = MARKERS,
  N = N_CELLS,
  donor_visit = FALSE,
  truth = list(phi = 100)
)
d0 <- sim0$design
pr <- cytoJointR:::cj_resolve_priors(list(), d0)

# Simulate, fit, rank ----------------------------------------------------------
for (i in sim_from:sim_to) {
  out <- file.path(run_dir, "sims", sprintf("sim_%03d.tsv", i))
  stopifnot("simulation file already exists" = !file.exists(out))
  set.seed(SEED_OFFSET + i)
  seed <- DATA_SEED_OFFSET + i
  par <- cytoJointR:::m1_prior_draw(d0, pr) # noncentered Stan scale
  sim <- cj_simulate_counts(
    n_donors = N_DONORS,
    markers = MARKERS,
    N = N_CELLS,
    donor_visit = FALSE,
    truth = list(
      alpha = par$alpha,
      gamma = par$gamma,
      theta0 = par$theta0,
      delta = par$delta,
      tau = par$tau,
      Omega = par$L_Omega %*% t(par$L_Omega),
      phi = 1 / par$inv_phi
    ),
    seed = seed
  )
  fit <- cj_fit(
    sim$data,
    chains = CHAINS,
    iter_warmup = ITER,
    iter_sampling = ITER,
    seed = seed,
    donor_visit = FALSE
  )
  # True values under their Stan names, from the simulator's truth.
  truth_vec <- cytoJointR:::m1_stan_vector(
    sim$truth[c("alpha", "gamma", "theta0", "delta", "tau", "inv_phi")],
    sim$design
  )
  truth_vec <- truth_vec[grepl(SBC_PARAMS, names(truth_vec))]
  draws <- fit$engine_fit$draws(names(truth_vec), format = "draws_matrix")
  draws <- draws[seq(1, nrow(draws), by = THIN), , drop = FALSE]
  fit_tab <- fit$diagnostics$summary |> select(variable, rhat, ess_bulk)
  tibble(
    sim = i,
    variable = names(truth_vec),
    truth = unname(truth_vec),
    rank = map_int(names(truth_vec), \(v) sum(draws[, v] < truth_vec[[v]])),
    n_draws = nrow(draws),
    n_divergent = fit$diagnostics$n_divergent,
    fit_seconds = fit$diagnostics$time_s
  ) |>
    left_join(fit_tab, by = "variable", relationship = "one-to-one") |>
    write_tsv(out)
}

write_session_info(
  file.path(run_dir, sprintf("session_info_%03d.txt", sim_from)),
  pkg_dir
)
