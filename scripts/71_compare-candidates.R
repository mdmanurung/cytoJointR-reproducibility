#!/usr/bin/env Rscript
# Purpose: one table of the dev metrics that decide "beats COMPASS", one row per leaderboard candidate.
# Run:     R_LIBS=renv/lib Rscript scripts/71_compare-candidates.R [out.tsv]
# Inputs:  benchmark/leaderboard.tsv [70_run-scorecard.R], benchmark/definition.tsv (comparator values)
# Outputs: out.tsv (default: printed only)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
METRICS <- c(
  auc = "rv144_dev.auc.P_RMS.fit",
  pauc10 = "rv144_dev.pauc10.P_RMS.fit",
  tpr_fpr10 = "rv144_dev.tpr_fpr10.P_RMS.fit",
  fpr_rule = "rv144_dev.fpr_rule.P_RMS.fit",
  tb_s = "tb_s.auc.P_RMS.fit",
  tb_u = "tb_u.auc.P_RMS.fit",
  sim_auc = "sim.auc.rare_weak.P_RMS.fit",
  sim_spear = "sim.spearman.rare_weak.RMS_mean.fit"
)

# Table ------------------------------------------------------------------------
lb <- read_tsv("benchmark/leaderboard.tsv", show_col_types = FALSE)
wide <- lb |>
  filter(metric_id %in% METRICS) |>
  mutate(name = names(METRICS)[match(metric_id, METRICS)]) |>
  select(candidate, package_sha, name, value) |>
  pivot_wider(names_from = name, values_from = value)
guards <- lb |>
  filter(str_detect(metric_id, "^sim.cov90")) |>
  summarise(cov_pass = sum(pass == "TRUE"), cov_n = n(), .by = candidate)
auc_diff <- lb |>
  filter(metric_id == METRICS[["auc"]]) |>
  select(candidate, auc_diff_compass_lo = diff_compass_lower)
verdict <- lb |>
  filter(metric_id == "win.beats_compass") |>
  select(candidate, beats_compass = value)
tab <- wide |>
  left_join(auc_diff, join_by(candidate)) |>
  left_join(guards, join_by(candidate)) |>
  left_join(verdict, join_by(candidate)) |>
  mutate(across(where(is.numeric), \(v) round(v, 3)))
print(as.data.frame(tab), row.names = FALSE)
out <- commandArgs(trailingOnly = TRUE)
if (length(out)) write_tsv(tab, out[1])
