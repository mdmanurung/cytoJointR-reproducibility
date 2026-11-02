#!/usr/bin/env Rscript
# Purpose: polyICSFlow (real package, environment B) on the counts of one benchmark replicate.
#          Expands the counts to one logical row per cell, then assignMarkerCombinations() and
#          calcPolyfunctionality() with the control as background_val. Records the clipped
#          perc_pos_backgroundSubtracted and the unclipped perc_pos - background.
# Run:     R_LIBS=renv/lib-r44 /exports/para-lipg-hpc/mdmanurung/batchdx-cyto/renv/bin/Rscript scripts/07_run-polyicsflow.R <run_dir> <scenario_id> <seed>
# Inputs:
#   <rep_dir>/counts.tsv, pattern_map.tsv   written by 06_run-benchmark-replicate.R (rep_dir must be complete)
#   renv/lib-r44/polyICSFlow                installed from GitHub SoegaardLab/polyICSFlow [install log in logs/]
# Outputs (in <run_dir>/reps/<scenario_id>/seed_<seed>/polyicsflow/):
#   donor_scores.tsv            method, score, donor_id, value: sum over patterns of net percent / 100 (fraction)
#   donor_magnitude.tsv         method, donor_id, estimate: the same sums, clipped and unclipped
#   pattern_scores.tsv          method, score, donor_id, pattern_key, value: per-pattern net / 100
#   runtime.tsv                 method, step, seconds, cpus; written last, so its presence marks a complete run
#   session_info.txt            R and package versions, including the polyICSFlow commit

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  library(polyICSFlow)
})

# Parameters -------------------------------------------------------------------
MARKERS <- c("IFNG", "TNF", "IL2", "CD107a") # column order of the logical matrix = marker order in MarkerComb (first J used)
BACKGROUND <- "control" # value of the condition column that marks the unstimulated sample
POP <- "CD4" # one parent population; calcPolyfunctionality() needs a pop_col

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 07_run-polyicsflow.R <run_dir> <scenario_id> <seed>" = length(
    args
  ) ==
    3
)
rep_dir <- file.path(
  args[1],
  "reps",
  args[2],
  sprintf("seed_%02d", as.integer(args[3]))
)
out_dir <- file.path(rep_dir, "polyicsflow")
if (read_tsv(file.path(args[1], "scenarios.tsv"), show_col_types = FALSE) |>
  filter(scenario_id == args[2]) |>
  pull(n_visits) == 2) {
  quit(save = "no") # polyICSFlow gives no visit change; 06 records the same for COMPASS
}
stopifnot(
  "replicate from 06 is not complete" = file.exists(file.path(
    rep_dir,
    "diagnostics.tsv"
  )),
  "polyicsflow output already exists" = !dir.exists(out_dir)
)
partial_dir <- paste0(out_dir, ".partial")
unlink(partial_dir, recursive = TRUE)
counts <- read_tsv(
  file.path(rep_dir, "counts.tsv"),
  col_types = cols(pattern_key = "c", donor_id = "c")
)
pattern_map <- read_tsv(
  file.path(rep_dir, "pattern_map.tsv"),
  col_types = cols(pattern_key = "c")
)

# Expand counts to cells -------------------------------------------------------
# One row per cell: logical marker columns (the only columns assignMarkerCombinations() may see),
# then donor and condition. Counts are exact, so the expansion loses nothing.
expand_seconds <- system.time({
  cells <- counts |> filter(n > 0) |> uncount(n)
  key_bits <- pattern_map$pattern_key |>
    strsplit("") |>
    map(\(b) b == "1") |>
    reduce(rbind)
  colnames(key_bits) <- MARKERS[seq_len(ncol(key_bits))] # J = 3 or 4
  bits <- as_tibble(key_bits[
    match(cells$pattern_key, pattern_map$pattern_key),
  ])
})[["elapsed"]]

# polyICSFlow ------------------------------------------------------------------
assign_seconds <- system.time(
  comb <- assignMarkerCombinations(bits, mode = "simple") |>
    mutate(
      donor_id = cells$donor_id,
      condition = cells$condition,
      cell_type = POP
    )
)[["elapsed"]]
calc_seconds <- system.time(
  res <- calcPolyfunctionality(
    comb,
    md_cols = "donor_id",
    pop_col = "cell_type",
    resolution = "MarkerComb",
    condition_col = "condition",
    background_val = BACKGROUND
  )
)[["elapsed"]]

# Net per donor and pattern ----------------------------------------------------
# Clipped: perc_pos_backgroundSubtracted as returned. Unclipped: perc_pos minus the control's perc_pos.
net <- res |>
  filter(condition == BACKGROUND) |>
  select(donor_id, MarkerComb, perc_background = perc_pos) |>
  inner_join(
    res |> filter(condition != BACKGROUND),
    join_by(donor_id, MarkerComb),
    relationship = "one-to-one"
  ) |>
  mutate(
    clipped = perc_pos_backgroundSubtracted / 100,
    unclipped = (perc_pos - perc_background) / 100,
    polyics_name = as.character(MarkerComb)
  ) |>
  left_join(
    pattern_map |> select(pattern_key, polyics_name),
    join_by(polyics_name),
    relationship = "many-to-one"
  ) |>
  select(donor_id, pattern_key, clipped, unclipped) |>
  pivot_longer(
    c(clipped, unclipped),
    names_to = "variant",
    values_to = "value"
  ) |>
  mutate(method = paste0("polyICSFlow_", variant))
donor_net <- net |> summarise(estimate = sum(value), .by = c(method, donor_id))
donor_scores <- donor_net |>
  transmute(method, score = "net_sum", donor_id, value = estimate)
pattern_scores <- net |>
  transmute(method, score = "net", donor_id, pattern_key, value)
runtime <- tibble(
  method = "polyICSFlow",
  step = c(
    "expand_counts_to_cells",
    "assignMarkerCombinations",
    "calcPolyfunctionality"
  ),
  seconds = c(expand_seconds, assign_seconds, calc_seconds),
  cpus = 1L
)

# Write outputs ----------------------------------------------------------------
dir.create(partial_dir, recursive = TRUE)
write_tsv(donor_scores, file.path(partial_dir, "donor_scores.tsv"))
write_tsv(
  donor_net |> select(method, donor_id, estimate),
  file.path(partial_dir, "donor_magnitude.tsv")
)
write_tsv(pattern_scores, file.path(partial_dir, "pattern_scores.tsv"))
writeLines(
  c(
    paste(
      "polyICSFlow",
      as.character(packageVersion("polyICSFlow")),
      packageDescription("polyICSFlow")$RemoteSha
    ),
    capture.output(sessionInfo())
  ),
  file.path(partial_dir, "session_info.txt")
)
write_tsv(runtime, file.path(partial_dir, "runtime.tsv"))
file.rename(partial_dir, out_dir)
