#!/usr/bin/env Rscript
# Purpose: apply the pre-registered decision rule of the redesign plan to the A/B/C summary, mechanically.
#          It reads nothing from the held-out generator D, which is summarised separately afterwards.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/11_apply-decision-rule.R <summary_dir> <out_dir>
# Inputs (in <summary_dir>, from 08_summarise-benchmark.R on the A/B/C run):
#   scenario_roles.tsv          scenario, J, role of each scenario
#   calibration.tsv             coverage and RMSE per scenario x model x target x grain x score
#   false_sign.tsv              false sign rate per scenario x model x target x grain x score
#   ipf.tsv                     ipf_failed rate per scenario x model
#   diagnostics.tsv             median fit seconds per scenario x model
# Outputs (in <out_dir>, new):
#   decision_cells.tsv          per scenario x model: coverage, false sign, ipf_failed and the pass flags
#   decision.tsv                per model: cells passed, calibrated, mean RMSE, gain over C0, runtime, selected
#   decision.md                 one paragraph generated from decision.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Parameters (each from the plan, "Decision rule, stated before any result") ---------
DECISION_TARGET <- "sample" # the cj_draws estimand stays the default (plan "Estimand"); the struct target is reported, not decided on
DECISION_SCORES <- c("RMS", "CPS", "CCS") # rule 1: RMS, CPS and CCS coverage
COVERAGE_MIN <- 0.85 # rule 1: coverage lower bound for a 90% interval
COVERAGE_MAX <- 0.97 # rule 1: coverage upper bound
COVERAGE_SE_MULT <- 2 # rule 1: "within 2 SE" of the band, using the seed-clustered coverage SE
NO_UPPER_BOUND_BASES <- "null" # amendment 2026-11-02, before any result: every truth is exactly 0 in null
                              # cells, where over-coverage is conservative, not harmful; false signs guard them
FALSE_SIGN_MAX <- 0.07 # rule 1: false sign rate at most this, in cells that have zero-truth rows (null, mixture)
IPF_FAILED_MAX <- 0.01 # rule 1: ipf_failed rate below this
RMSE_BASES <- c("rare_weak", "rare_moderate") # rule 3: mean RMSE on rare_weak plus rare_moderate (J = 3 and 4)
RMSE_SCORE <- "RMS" # rule 3: donor-level RMS is the magnitude score of the benchmark
MIN_GAIN <- 0.10 # rule 4: a structured model must beat C0 by at least this share of RMSE
TIE_TOL <- 0.05 # rule 5: RMSE within this share of the best is a tie
BASELINE <- "C0" # rule 4: the default if no structured model earns its place
COMPLEXITY <- c(C0 = 0, C1 = 1, C2 = 2, C2phi = 3, C3 = 4) # rule 5: fewer parameters first (C2 adds phi to C1; C2phi one per condition; C3 adds interaction SDs)

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 11_apply-decision-rule.R <summary_dir> <out_dir>" = length(
    args
  ) ==
    2,
  "out_dir already exists" = !dir.exists(args[2])
)
summary_dir <- args[1]
out_dir <- args[2]
read_summary <- \(f) read_tsv(file.path(summary_dir, f), show_col_types = FALSE)
roles <- read_summary("scenario_roles.tsv") |> filter(role == "main")
calibration <- read_summary("calibration.tsv") |>
  filter(
    scenario_id %in% roles$scenario_id,
    target == DECISION_TARGET,
    draws == "sample",
    grain == "global"
  )
false_sign <- read_summary("false_sign.tsv") |>
  filter(
    scenario_id %in% roles$scenario_id,
    target == DECISION_TARGET,
    draws == "sample",
    grain == "global"
  )
ipf <- read_summary("ipf.tsv") |> filter(scenario_id %in% roles$scenario_id, draws == "sample")
diagnostics <- read_summary("diagnostics.tsv") |>
  filter(scenario_id %in% roles$scenario_id)

# Rule 1: pass a cell ----------------------------------------------------------
# Coverage of every decision score in the band, false sign rate within the limit where it is defined,
# ipf_failed below the limit. A missing value fails the cell.
cells <- calibration |>
  filter(score %in% DECISION_SCORES) |>
  mutate(
    lower_ok = coverage + COVERAGE_SE_MULT * coverage_se >= COVERAGE_MIN,
    upper_ok = coverage - COVERAGE_SE_MULT * coverage_se <= COVERAGE_MAX |
      str_remove(scenario_id, "_J[0-9]+$") %in% NO_UPPER_BOUND_BASES
  ) |>
  summarise(
    coverage_min = min(coverage),
    coverage_max = max(coverage),
    band_ok = all(lower_ok & upper_ok),
    n_scores = n(),
    .by = c(scenario_id, model)
  ) |>
  left_join(
    false_sign |>
      filter(score %in% DECISION_SCORES) |>
      summarise(false_sign_max = max(false_sign), .by = c(scenario_id, model)),
    join_by(scenario_id, model),
    relationship = "one-to-one"
  ) |>
  left_join(
    ipf |> select(scenario_id, model, ipf_failed_rate),
    join_by(scenario_id, model),
    relationship = "one-to-one"
  ) |>
  mutate(
    coverage_ok = n_scores == length(DECISION_SCORES) & band_ok,
    false_sign_ok = is.na(false_sign_max) | false_sign_max <= FALSE_SIGN_MAX, # NA: no zero-truth rows in this cell
    ipf_ok = !is.na(ipf_failed_rate) & ipf_failed_rate < IPF_FAILED_MAX,
    pass = coalesce(coverage_ok, FALSE) & false_sign_ok & ipf_ok
  )

# Rules 2 to 5: calibrated, best, tie, default --------------------------------
models <- cells |>
  summarise(
    n_cells = n(),
    n_pass = sum(pass),
    n_fail_coverage = sum(!coverage_ok),
    n_fail_false_sign = sum(!false_sign_ok),
    n_fail_ipf = sum(!ipf_ok),
    .by = model
  ) |>
  mutate(calibrated = n_cells == nrow(roles) & n_pass == n_cells) |> # rule 2: passes every A/B/C cell
  left_join(
    calibration |>
      filter(
        score == RMSE_SCORE,
        str_remove(scenario_id, "_J[0-9]+$") %in% RMSE_BASES
      ) |>
      summarise(mean_rmse = mean(rmse), .by = model),
    join_by(model),
    relationship = "one-to-one"
  ) |>
  left_join(
    diagnostics |>
      summarise(fit_seconds = median(fit_seconds_median), .by = model),
    join_by(model),
    relationship = "one-to-one"
  ) |>
  mutate(complexity = COMPLEXITY[model])
rmse_baseline <- models$mean_rmse[models$model == BASELINE]
structured <- models |> filter(calibrated, model != BASELINE)
# Rule 3 with rule 5: lowest mean RMSE among calibrated structured models; models within TIE_TOL of the
# lowest tie, and the tie goes to the fewer parameters, then the shorter runtime.
best <- structured |>
  filter(mean_rmse <= (1 + TIE_TOL) * min(mean_rmse, Inf)) |>
  arrange(complexity, fit_seconds) |>
  slice(1)
gain <- if (nrow(best) == 1) 1 - best$mean_rmse / rmse_baseline else NA_real_
default_model <- NA_character_ # no model passes every cell: no default is set
if (any(models$calibrated)) {
  default_model <- BASELINE # no structured model is calibrated, or none earns its place
}
if (any(models$calibrated) && nrow(best) == 1 && gain >= MIN_GAIN) {
  default_model <- best$model # rule 4: the structured model earns its place
}
models <- models |>
  mutate(
    gain_over_baseline = 1 - mean_rmse / rmse_baseline,
    selected = model %in% default_model
  ) |>
  arrange(complexity)

# Decision text ----------------------------------------------------------------
calibrated_names <- models |> filter(calibrated) |> pull(model)
fail_text <- models |>
  filter(!calibrated) |>
  mutate(
    text = sprintf(
      "%s passes %d of %d cells (coverage fails in %d, false sign in %d, ipf_failed in %d)",
      model,
      n_pass,
      n_cells,
      n_fail_coverage,
      n_fail_false_sign,
      n_fail_ipf
    )
  ) |>
  pull(text)
baseline_text <- sprintf(
  "C0 has mean %s RMSE %.2e on %s",
  RMSE_SCORE,
  rmse_baseline,
  paste(RMSE_BASES, collapse = " plus ")
)
verdict <- sprintf("No model passes every cell, so the rule sets no default.")
if (!is.na(default_model) && nrow(best) == 0) {
  verdict <- sprintf(
    "No structured model is calibrated, so the default is %s.",
    default_model
  )
}
if (!is.na(default_model) && nrow(best) == 1) {
  verdict <- sprintf(
    "%s is the best calibrated structured model (mean RMSE %.2e, %.1f%% below C0, at least %.0f%% required), so %s.",
    best$model,
    best$mean_rmse,
    100 * gain,
    100 * MIN_GAIN,
    if (default_model == best$model) {
      "it becomes the default"
    } else {
      "C0 stays the default and Stan becomes opt-in"
    }
  )
}
decision_text <- paste0(
  sprintf(
    "Decision on the %s target over %d A/B/C cells: default model %s. ",
    DECISION_TARGET,
    nrow(roles),
    coalesce(default_model, "none")
  ),
  if (length(calibrated_names)) {
    sprintf("Calibrated models: %s. ", paste(calibrated_names, collapse = ", "))
  } else {
    "No model is calibrated. "
  },
  if (length(fail_text)) {
    paste0("Not calibrated: ", paste(fail_text, collapse = "; "), ". ")
  } else {
    ""
  },
  baseline_text,
  ". ",
  verdict,
  sprintf(
    " A cell passes when RMS, CPS and CCS coverage lie in %.2f to %.2f, the false sign rate is at most %.2f and ipf_failed is below %.0f%%.",
    COVERAGE_MIN,
    COVERAGE_MAX,
    FALSE_SIGN_MAX,
    100 * IPF_FAILED_MAX
  )
)

# Write outputs ----------------------------------------------------------------
dir.create(out_dir, recursive = TRUE)
write_tsv(cells, file.path(out_dir, "decision_cells.tsv"))
write_tsv(models, file.path(out_dir, "decision.tsv"))
writeLines(decision_text, file.path(out_dir, "decision.md"))
