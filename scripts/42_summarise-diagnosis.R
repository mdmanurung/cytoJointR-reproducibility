#!/usr/bin/env Rscript
# Purpose: summarise the CPS over-coverage experiments: 90% coverage and mean width per scenario x prior x draws kind x score,
#          against both truths, plus the fitted phi and the size of the Dirichlet shape parameters.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/42_summarise-diagnosis.R runs/2026-11-03_diagnose-cps
# Inputs (in <run_dir>): scores_<scenario>.tsv, phi_<scenario>.tsv, patterns_<scenario>.tsv [40_diagnose-cps.R]
# Outputs (in <run_dir>):
#   coverage.tsv    per scenario x prior x draws x score x target: n contrasts, coverage, coverage SE across seeds, mean width
#   phi.tsv         all phi_<scenario>.tsv rows
#   patterns.tsv    all patterns_<scenario>.tsv rows

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

run_dir <- commandArgs(trailingOnly = TRUE)[1]
read_all <- function(prefix) {
  list.files(run_dir, pattern = paste0("^", prefix, "_.*[.]tsv$"), full.names = TRUE) |>
    set_names(\(f) sub(paste0("^", prefix, "_(.*)[.]tsv$"), "\\1", basename(f))) |>
    map(\(f) read_tsv(f, show_col_types = FALSE) |> select(-any_of("scenario_id"))) |> # scores files carry scenario_id already
    list_rbind(names_to = "scenario_id")
}

# Coverage ---------------------------------------------------------------------
# Coverage per seed first, so the SE counts replicates, not donors.
coverage <- read_all("scores") |>
  pivot_longer(c(truth_sample, truth_struct), names_to = "target", values_to = "truth", names_prefix = "truth_") |>
  summarise(
    covered = mean(coalesce(lower <= truth & truth <= upper, FALSE)), width = mean(upper - lower, na.rm = TRUE), n = n(),
    .by = c(scenario_id, prior, draws, score, target, seed)
  ) |>
  summarise(
    n_seeds = n(), n_contrasts = sum(n), coverage = mean(covered), coverage_se = sd(covered) / sqrt(n()), width = mean(width),
    .by = c(scenario_id, prior, draws, score, target)
  ) |>
  arrange(scenario_id, score, target, prior, draws)
write_tsv(coverage, file.path(run_dir, "coverage.tsv"))
write_tsv(read_all("phi"), file.path(run_dir, "phi.tsv"))
write_tsv(read_all("patterns"), file.path(run_dir, "patterns.tsv"))
