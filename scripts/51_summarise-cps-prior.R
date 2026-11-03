#!/usr/bin/env Rscript
# Purpose: summarise 50_diagnose-cps-prior.R: CPS and RMS bias and 90% coverage by variant, scenario and draw kind.
#          Each seed is one unit (donor contrasts are averaged within a seed); the MC SE is the SD across seeds / sqrt(seeds).
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/51_summarise-cps-prior.R runs/2026-11-03_cps-prior
# Inputs:
#   <run_dir>/scores_<scenario>.tsv, fits_<scenario>.tsv, theta_<scenario>.tsv   [50_diagnose-cps-prior.R]
# Outputs (in <run_dir>):
#   summary.tsv      variant x scenario x draws x score: n_seeds, bias, bias_se, coverage, coverage_se, mean truth
#   summary_fits.tsv variant x scenario: mean fit seconds, divergences, max R-hat
#   summary_theta.tsv variant x scenario x parameter: mean over seeds of posterior mean and SD, truth

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
})

# Read -------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[1]
read_all <- \(prefix) {
  list.files(run_dir, pattern = paste0("^", prefix, "_.*\\.tsv$"), full.names = TRUE) |>
    set_names(\(f) sub(paste0("^", prefix, "_(.*)\\.tsv$"), "\\1", basename(f))) |>
    map(read_tsv, show_col_types = FALSE) |>
    list_rbind(names_to = "scenario")
}
scores <- read_all("scores")

# Summarise --------------------------------------------------------------------
# draws "structured" is scored against the structured truth, draws "sample" against the realised-sample truth.
per_seed <- scores |>
  mutate(truth = if_else(draws == "structured", truth_struct, truth_sample)) |>
  summarise(
    bias = mean(estimate - truth),
    coverage = mean(lower <= truth & truth <= upper),
    truth = mean(truth),
    .by = c(scenario, variant, draws, score, seed)
  )
summary_tbl <- per_seed |>
  summarise(
    n_seeds = n(),
    bias_se = sd(bias) / sqrt(n_seeds),
    bias = mean(bias),
    coverage_se = sd(coverage) / sqrt(n_seeds),
    coverage = mean(coverage),
    truth = mean(truth),
    .by = c(scenario, variant, draws, score)
  )
write_tsv(summary_tbl, file.path(run_dir, "summary.tsv"))
write_tsv(
  read_all("fits") |>
    summarise(n_fits = n(), fit_seconds = mean(fit_seconds), n_divergent = sum(n_divergent),
              max_rhat = max(max_rhat), phi_mean = mean(phi_mean), .by = c(scenario, variant)),
  file.path(run_dir, "summary_fits.tsv")
)
write_tsv(
  read_all("theta") |>
    filter(!grepl("raw", parameter)) |>
    summarise(mean = mean(mean), sd = mean(sd), truth = first(truth), .by = c(scenario, variant, parameter)),
  file.path(run_dir, "summary_theta.tsv")
)
