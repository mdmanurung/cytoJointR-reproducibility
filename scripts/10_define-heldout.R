#!/usr/bin/env Rscript
# Purpose: define the held-out generator-D run (plan: "Held-out generator D"). Takes the solved rare_weak and
#          moderate scenarios of the A/B/C run and turns each into a D1 (cell-level phenotype mixture) and a D2
#          (generator A with t3 donor effects) scenario at J = 3 and 4. D is for vetoing a winner only: run it
#          after 11_apply-decision-rule.R has written decision.tsv.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/10_define-heldout.R <main_run_dir> <d_run_dir> <n_seeds>
# Inputs:
#   <main_run_dir>/scenarios.tsv   solved alpha, beta, tau of the A/B/C run [05_define-benchmark.R]
# Outputs (in <d_run_dir>, new):
#   scenarios.tsv               one row per D scenario, seed offset 900000 (pre-registered)
#   jobs.tsv                    one row per replicate: task, scenario_id, seed

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
BASE_SCENARIOS <- c("rare_weak", "moderate") # plan: scenarios rare_weak and moderate
GENERATORS <- c("D1", "D2") # plan: cell-level phenotype mixture, t3 donor effects
SEED_OFFSET_D <- 900000L # plan: seed offset 900000; replicate data seed = this + 1000 * scenario_index + seed
PACK <- 6L # replicates per SLURM task, as in 05

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 10_define-heldout.R <main_run_dir> <d_run_dir> <n_seeds>" = length(
    args
  ) ==
    3,
  "d_run_dir already exists" = !dir.exists(args[2])
)
d_run_dir <- args[2]
n_seeds <- as.integer(args[3])

# Scenarios --------------------------------------------------------------------
scenarios <- read_tsv(
  file.path(args[1], "scenarios.tsv"),
  show_col_types = FALSE
) |>
  filter(scenario %in% BASE_SCENARIOS, role == "main") |>
  cross_join(tibble(generator_d = GENERATORS)) |>
  mutate(
    scenario = paste0(generator_d, "_", scenario),
    scenario_id = paste0(scenario, "_J", J),
    generator = generator_d,
    scenario_index = row_number(),
    seed_offset = SEED_OFFSET_D + 1000L * scenario_index,
    n_seeds = n_seeds
  ) |>
  select(-generator_d) # keep J and every other column the replicate script reads
jobs <- scenarios |>
  mutate(seed = map(n_seeds, seq_len)) |>
  select(scenario_id, seed) |>
  unnest(seed) |>
  mutate(task = (row_number() - 1L) %/% PACK + 1L, .before = 1)

# Write outputs ----------------------------------------------------------------
dir.create(d_run_dir, recursive = TRUE)
write_tsv(scenarios, file.path(d_run_dir, "scenarios.tsv"))
write_tsv(jobs, file.path(d_run_dir, "jobs.tsv"))
