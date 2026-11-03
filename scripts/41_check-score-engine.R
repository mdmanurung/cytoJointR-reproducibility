#!/usr/bin/env Rscript
# Purpose: check the score engine on known draws, and show how Dirichlet noise at fixed phi spreads CPS and RMS.
#          (1) constant draws must reproduce cj_score_reference(); (2) for one contrast whose true p is the structured
#          mu, draw p ~ Dirichlet((phi + N) mu) (the posterior when y = N mu) and report the sd of RMS and CPS draws.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/41_check-score-engine.R runs/2026-11-03_diagnose-cps
# Inputs:
#   cytoJointR-wt/redesign-run/            frozen package (49278c9)
#   runs/2026-11-02_redesign/scenarios.tsv, scripts/helper02_benchmark-sim.R
# Outputs (in <run_dir>):
#   engine_check.tsv       max abs difference between score_contrast() and cj_score_reference() on constant draws
#   noise_by_phi.tsv       per scenario x phi: sd of RMS and CPS draws, and median a_k of polyfunctional patterns

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
SCENARIOS <- c("rare_moderate_J3", "moderate_J3")
PHI <- c(1e3, 1e5, 1e6, 1e8, 1e10) # 1e3-1e5: posterior range seen on generator A; 1e10: the multinomial limit
N_DRAWS <- 4000
SEED <- 1
SCORE_COLS <- c("RMS", "poly_absolute", "poly_marginal", "CPS", "CCS", "breadth")

run_dir <- commandArgs(trailingOnly = TRUE)[1]
stopifnot("run_dir must exist" = dir.exists(run_dir))
pkgload::load_all("cytoJointR-wt/redesign-run", quiet = TRUE)
source("scripts/helper02_benchmark-sim.R")
scenarios <- read_tsv("runs/2026-11-02_redesign/scenarios.tsv", show_col_types = FALSE)
set.seed(SEED)

# Engine on constant draws -----------------------------------------------------
# score_contrast() with every draw equal to the true p must give the reference scores.
one_contrast <- function(scenario_arg) {
  sc <- filter(scenarios, scenario_id == scenario_arg)
  sim <- bench_simulate(sc, sc$seed_offset + SEED)
  ct <- sim$data$contrasts |> filter(!is.na(control_sample_id)) |> slice(1)
  list(sc = sc, bits = sim$data$bits, mu_u = sim$p_struct[ct$control_sample_id, ], mu_s = sim$p_struct[ct$stim_sample_id, ])
}
cases <- set_names(SCENARIOS) |> map(one_contrast)

engine_check <- imap(cases, \(x, nm) {
  ref <- cj_score_reference(x$mu_u, x$mu_s, x$bits)$global
  G <- cytoJointR:::score_contrast(matrix(x$mu_u, 5, length(x$mu_u), byrow = TRUE),
                                   matrix(x$mu_s, 5, length(x$mu_s), byrow = TRUE), x$bits)$G
  tibble(scenario_id = nm, score = SCORE_COLS, reference = unname(ref[SCORE_COLS]),
         engine = G[1, SCORE_COLS], abs_diff = abs(reference - engine))
}) |>
  list_rbind()
write_tsv(engine_check, file.path(run_dir, "engine_check.tsv"))

# Dirichlet noise by phi -------------------------------------------------------
rdirichlet <- function(n, a) {
  g <- matrix(rgamma(n * length(a), shape = rep(a, each = n)), n)
  g / rowSums(g)
}
noise <- cross_join(tibble(scenario_id = SCENARIOS), tibble(phi = PHI)) |>
  pmap(\(scenario_id, phi) {
    x <- cases[[scenario_id]]
    N <- x$sc$N
    a_u <- (phi + N) * x$mu_u
    a_s <- (phi + N) * x$mu_s
    G <- cytoJointR:::score_contrast(rdirichlet(N_DRAWS, a_u), rdirichlet(N_DRAWS, a_s), x$bits)$G
    poly <- rowSums(x$bits) >= 2
    tibble(scenario_id, phi,
           true_CPS = cj_score_reference(x$mu_u, x$mu_s, x$bits)$global[["CPS"]],
           true_RMS = cj_score_reference(x$mu_u, x$mu_s, x$bits)$global[["RMS"]],
           sd_CPS = sd(G[, "CPS"], na.rm = TRUE), sd_RMS = sd(G[, "RMS"], na.rm = TRUE),
           median_a_poly = median(a_s[poly]), min_a_poly = min(a_s[poly]))
  }) |>
  list_rbind()
write_tsv(noise, file.path(run_dir, "noise_by_phi.tsv"))
writeLines(capture.output(sessionInfo()), file.path(run_dir, "session_info_engine.txt"))
