#!/usr/bin/env Rscript
# Purpose: collect the per-replicate TSVs of the redesign evaluation into metric tables. Candidate models
#          (C0, C1, C2, C2phi): 90% coverage with a seed-clustered SE, bias, RMSE, interval width relative to C0,
#          false sign rate, ipf_failed rate, runtime and diagnostics, each against both the realised sample
#          probabilities (target "sample") and the structured ones (target "struct"). Candidates, COMPASS and
#          polyICSFlow: donor-response detection (ROC, AUC), magnitude accuracy and pattern-level AUC.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/08_summarise-benchmark.R <run_dir> <out_dir>
# Inputs:
#   <run_dir>/scenarios.tsv, jobs.tsv            [05_define-benchmark.R]
#   <run_dir>/reps/*/seed_*/*.tsv                [06_run-benchmark-replicate.R]
#   <run_dir>/reps/*/seed_*/polyicsflow/*.tsv    [07_run-polyicsflow.R]
#   scripts/helper01_true-scores.R               write_session_info() [this repo]
# Outputs (in <out_dir>, new):
#   scenario_roles.tsv          scenario_id, scenario, generator, J, role, n_donors (read by 11_apply-decision-rule.R)
#   replicate_status.tsv        scenario x seed: whether 06 and 07 finished, package commit, COMPASS status per variant
#   model_status.tsv            scenario x model: replicates ok / failed / not_run
#   truth_summary.tsv           achieved true donor RMS quantiles and responder / responding-pattern fractions
#   calibration_replicate.tsv   per seed: coverage, bias, MSE, width per scenario x model x draws x target x grain x score (x degree)
#   calibration.tsv             across seeds: coverage with seed-clustered SE, bias, RMSE, width, width ratio to C0
#   false_sign_replicate.tsv    per seed: share of zero-truth rows with P(> eps) or P(< -eps) above SIGN_PROB
#   false_sign.tsv              across seeds: false sign rate with seed-clustered SE
#   ipf.tsv                     ipf_failed rate (mean share of failed draws) and withheld summaries per scenario x model
#   change_replicate.tsv        two-visit scenarios, per seed: coverage, bias, MSE of the v1 -> v2 change
#   change.tsv                  across seeds: the same with seed-clustered SE
#   donor_table.tsv             scenario, seed, donor, detector, value, truth: every score used for detection
#   roc_curves.tsv              pooled ROC points per scenario, target, comparison, detector
#   auc.tsv                     AUC per scenario, target, comparison, detector, with seed-bootstrap SE
#   auc_replicate.tsv           the same AUC within one seed (one point per replicate)
#   rules.tsv                   sensitivity and specificity at each method's conventional rule, with SE
#   magnitude_replicate.tsv     per seed: bias and RMSE of the donor-RMS estimate (all methods), coverage (candidates)
#   magnitude.tsv               per scenario, target and estimator: mean bias, RMSE, coverage with across-seed SE
#   pattern_auc.tsv             AUC over donor x pattern per scenario, target and detector, with seed-bootstrap SE
#   runtime_replicate.tsv       seconds per method and replicate
#   runtime.tsv                 median and range of seconds per method and scenario
#   diagnostics.tsv             candidate convergence per scenario and model
#   session_info_summary.txt    R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})
source("scripts/helper01_true-scores.R") # write_session_info()

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign") # only for the commit recorded in the session info
BOOT_N <- 1000L # seed-bootstrap resamples for AUC, sensitivity and specificity SEs
BOOT_SEED <- 20261104L # seeds the bootstrap
COMPARISONS <- c("within", "vs_null") # within: responders vs non-responders of the scenario; vs_null: vs all donors of the J-matched null scenario
TARGETS <- c(sample = "sample", struct = "struct") # truth: realised sample probability, or the structured mean
MODELS <- c("C0", "C1", "C2", "C2phi") # candidate ids written by 06
# Conventional rule per detector: score > threshold calls a donor a responder.
# COMPASS Response() > 0.5 is our choice; COMPASS gives no rule for it.
RULES <- tribble(
  ~detector                     , ~threshold ,
  "C0 P_RMS_gt_eps"             , 0.95       ,
  "C1 P_RMS_gt_eps"             , 0.95       ,
  "C2 P_RMS_gt_eps"             , 0.95       ,
  "C2phi P_RMS_gt_eps"          , 0.95       ,
  "COMPASS_default Response"    , 0.5        ,
  "COMPASS_nofilter Response"   , 0.5        ,
  "polyICSFlow_clipped net_sum" , 0
)
# Detectors that fill a missing value with 0: a COMPASS fit that failed reports no response.
COMPASS_METHODS <- c("COMPASS_default", "COMPASS_nofilter")
SIGN_PROB <- 0.95 # a false sign is P(change > eps) or P(change < -eps) above this while the truth is exactly 0, as in the pilot
ZERO_TOL <- 1e-12 # a truth this close to 0 is a null truth (null cell, mixture non-responder)
RHAT_MAX <- 1.01 # mixing flags, as in the pilot
ESS_MIN <- 400

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 08_summarise-benchmark.R <run_dir> <out_dir>" = length(
    args
  ) ==
    2,
  "out_dir already exists" = !dir.exists(args[2])
)
run_dir <- args[1]
out_dir <- args[2]
scenarios <- read_tsv(
  file.path(run_dir, "scenarios.tsv"),
  show_col_types = FALSE
)
jobs <- read_tsv(file.path(run_dir, "jobs.tsv"), show_col_types = FALSE)
one_visit_ids <- scenarios |> filter(n_visits == 1) |> pull(scenario_id)
set.seed(BOOT_SEED)

# Replicate status -------------------------------------------------------------
replicate_status <- jobs |>
  select(scenario_id, seed) |>
  mutate(
    rep_dir = file.path(
      run_dir,
      "reps",
      scenario_id,
      sprintf("seed_%02d", seed)
    ),
    done_a = file.exists(file.path(rep_dir, "diagnostics.tsv")),
    done_b = file.exists(file.path(rep_dir, "polyicsflow", "runtime.tsv")),
    # the package is loaded from source at run time, so the commit is recorded per replicate
    package_commit = map_chr(rep_dir, \(d) {
      f <- file.path(d, "session_info.txt")
      if (file.exists(f)) str_remove(readLines(f, 1), "^cytoJointR commit: ") else NA_character_
    })
  )
# Read one file of every finished replicate (done_a, or done_b for polyicsflow/) that has it.
read_reps <- function(
  file,
  subdir = "",
  need = "done_a",
  ids = scenarios$scenario_id
) {
  replicate_status |>
    filter(.data[[need]], scenario_id %in% ids) |>
    mutate(path = file.path(rep_dir, subdir, file)) |>
    select(scenario_id, seed, path) |>
    pmap(\(scenario_id, seed, path) {
      # pattern_key ("0001"), donor_id and marker must stay text; name only the columns this file has
      text_cols <- intersect(
        c("pattern_key", "donor_id", "marker", "contrast_id"),
        names(read_tsv(path, n_max = 0, show_col_types = FALSE))
      )
      read_tsv(
        path,
        show_col_types = FALSE,
        col_types = rlang::inject(cols(
          !!!rlang::set_names(
            rep(list(col_character()), length(text_cols)),
            text_cols
          ),
          .default = col_guess()
        ))
      ) |>
        mutate(scenario_id = scenario_id, seed = seed, .before = 1)
    }) |>
    list_rbind()
}
compass_status <- read_reps("compass_status.tsv")
replicate_status <- replicate_status |>
  left_join(
    compass_status |>
      select(scenario_id, seed, method, status) |>
      pivot_wider(
        names_from = method,
        values_from = status,
        names_prefix = "status_"
      ),
    join_by(scenario_id, seed),
    relationship = "one-to-one"
  )
model_status <- read_reps("model_status.tsv") |>
  count(scenario_id, model, status) |>
  pivot_wider(names_from = status, values_from = n, values_fill = 0L)

# Truth ------------------------------------------------------------------------
truth_donor <- read_reps("truth_donor.tsv")
truth_pattern <- read_reps("truth_pattern.tsv")
truth_summary <- truth_donor |>
  summarise(
    n_donors = n(),
    rms_q05 = quantile(true_rms_sample, 0.05),
    rms_q25 = quantile(true_rms_sample, 0.25),
    rms_median = median(true_rms_sample),
    rms_q75 = quantile(true_rms_sample, 0.75),
    rms_q95 = quantile(true_rms_sample, 0.95),
    rms_struct_median = median(true_rms_struct),
    frac_responders = mean(responder_sample),
    .by = scenario_id
  ) |>
  left_join(
    truth_pattern |>
      summarise(
        frac_patterns_responding = mean(responding_sample),
        frac_patterns_down = mean(delta_sample < -min(scenarios$eps_pattern)),
        .by = scenario_id
      ),
    join_by(scenario_id),
    relationship = "one-to-one"
  )

# Calibration of the candidates ------------------------------------------------
# One row per replicate x model x contrast x score (x marker, pattern) with both truths. A summary the model
# withheld (failure tolerance) has NA bounds and counts as not covered.
model_scores <- read_reps("model_scores.tsv") |>
  mutate(
    degree = if_else(
      grain == "pattern",
      str_count(pattern_key, "1"),
      NA_integer_
    )
  )
cal_keys <- c("scenario_id", "seed", "model", "draws", "grain", "score", "degree")
calibration_replicate <- TARGETS |>
  map(\(target) {
    model_scores |>
      mutate(
        truth = .data[[paste0("truth_", target)]],
        err = estimate - truth
      ) |>
      summarise(
        n = n(),
        coverage = mean(coalesce(lower <= truth & truth <= upper, FALSE)),
        bias = mean(err, na.rm = TRUE),
        mse = mean(err^2, na.rm = TRUE),
        width = mean(upper - lower, na.rm = TRUE),
        .by = all_of(cal_keys)
      )
  }) |>
  list_rbind(names_to = "target")
across_keys <- c("scenario_id", "model", "draws", "target", "grain", "score", "degree")
calibration <- calibration_replicate |>
  summarise(
    n_seeds = n(),
    coverage_se = sd(coverage) / sqrt(n()), # seed-clustered: the replicate is the unit
    bias_se = sd(bias) / sqrt(n()),
    rmse_se = sd(mse) / sqrt(n()) / (2 * sqrt(mean(mse))), # delta method on the mean MSE
    coverage = mean(coverage),
    bias = mean(bias),
    rmse = sqrt(mean(mse)),
    width = mean(width),
    .by = all_of(across_keys)
  )
calibration <- calibration |>
  left_join(
    calibration |>
      filter(model == "C0") |>
      select(scenario_id, target, draws, grain, score, degree, width_c0 = width),
    join_by(scenario_id, target, draws, grain, score, degree),
    relationship = "many-to-one"
  ) |>
  mutate(width_ratio = width / width_c0, .keep = "unused") |>
  relocate(coverage_se, bias_se, rmse_se, .after = coverage)

# False sign -------------------------------------------------------------------
# Rows whose truth is exactly 0 (null cells, mixture non-responders): a claim of a sign is a false sign.
false_sign_replicate <- TARGETS |>
  map(\(target) {
    model_scores |>
      filter(abs(.data[[paste0("truth_", target)]]) <= ZERO_TOL) |>
      summarise(
        n = n(),
        false_sign = mean(coalesce(
          p_positive > SIGN_PROB | p_negative > SIGN_PROB,
          FALSE
        )),
        .by = all_of(cal_keys)
      )
  }) |>
  list_rbind(names_to = "target")
false_sign <- false_sign_replicate |>
  summarise(
    n_seeds = n(),
    n_rows = sum(n),
    false_sign_se = sd(false_sign) / sqrt(n()),
    false_sign = mean(false_sign),
    .by = all_of(across_keys)
  ) |>
  relocate(false_sign_se, .after = false_sign)

# ipf_failed rate --------------------------------------------------------------
# Share of draws whose IPF failed, averaged over contrasts and replicates; same draws for every score, so RMS rows.
ipf <- model_scores |>
  filter(grain == "global", score == "RMS") |>
  summarise(
    n_contrasts = n(),
    ipf_failed_rate = mean(failure_fraction, na.rm = TRUE),
    n_withheld = sum(status == "failure_tolerance_exceeded", na.rm = TRUE),
    .by = c(scenario_id, model, draws)
  )

# Visit change -----------------------------------------------------------------
change_scores <- read_reps(
  "model_change.tsv",
  ids = setdiff(scenarios$scenario_id, one_visit_ids)
)
# No two-visit scenario (the held-out run) leaves both tables empty.
change_replicate <- tibble()
change <- tibble()
if (nrow(change_scores) > 0) {
  change_replicate <- TARGETS |>
    map(\(target) {
      change_scores |>
        mutate(
          truth = .data[[paste0("truth_", target)]],
          err = estimate - truth
        ) |>
        summarise(
          n = n(),
          coverage = mean(coalesce(lower <= truth & truth <= upper, FALSE)),
          bias = mean(err, na.rm = TRUE),
          mse = mean(err^2, na.rm = TRUE),
          width = mean(upper - lower, na.rm = TRUE),
          .by = c(scenario_id, seed, model, grain, score)
        )
    }) |>
    list_rbind(names_to = "target")
  change <- change_replicate |>
    summarise(
      n_seeds = n(),
      coverage_se = sd(coverage) / sqrt(n()),
      coverage = mean(coverage),
      bias = mean(bias),
      rmse = sqrt(mean(mse)),
      width = mean(width),
      .by = c(scenario_id, model, target, grain, score)
    )
}

# Donor table ------------------------------------------------------------------
# Truth rows x every detector, so a method that produced nothing for a replicate still has rows.
# A COMPASS value missing because the fit failed becomes 0 (COMPASS reports no response) and is
# flagged; polyICSFlow rows exist only for replicates where 07 finished. One-visit scenarios only.
attach_scores <- function(truth, scores, by) {
  truth |>
    filter(scenario_id %in% one_visit_ids) |>
    cross_join(
      scores |> distinct(across(any_of(c("method", "score", "detector"))))
    ) |>
    left_join(
      replicate_status |> select(scenario_id, seed, done_b),
      join_by(scenario_id, seed),
      relationship = "many-to-one"
    ) |>
    filter(!str_detect(method, "^polyICSFlow") | done_b) |>
    select(-done_b) |>
    left_join(scores, by = by, relationship = "one-to-one") |>
    mutate(
      filled = is.na(value) & method %in% COMPASS_METHODS,
      value = if_else(filled, 0, value)
    ) |>
    filter(!is.na(value))
}
read_methods <- function(file) {
  bind_rows(
    read_reps(file, ids = one_visit_ids),
    read_reps(file, "polyicsflow", "done_b", ids = one_visit_ids)
  )
}
donor_table <- attach_scores(
  truth_donor,
  read_methods("donor_scores.tsv"),
  c("scenario_id", "seed", "donor_id", "method", "score")
) |>
  mutate(detector = paste(method, score))

# Detection: ROC, AUC, rules ---------------------------------------------------
# d: rows with seed, value, y (TRUE = responder). Statistic f(value, y); SE from resampling seeds.
auc <- \(x, y) {
  n_pos <- sum(y)
  n_neg <- sum(!y)
  if (n_pos == 0 || n_neg == 0) {
    return(NA_real_)
  }
  (sum(rank(x)[y]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
boot_stat <- function(d, f) {
  rows <- split(seq_len(nrow(d)), d$seed)
  draws <- replicate(BOOT_N, {
    idx <- unlist(rows[sample(length(rows), replace = TRUE)], use.names = FALSE)
    f(d$value[idx], d$y[idx])
  })
  tibble(
    estimate = f(d$value, d$y),
    se = sd(draws, na.rm = TRUE),
    n_pos = sum(d$y),
    n_neg = sum(!d$y),
    n_seeds = length(rows)
  )
}
donor_by_detector <- split(donor_table, donor_table$detector)
# Positives are the scenario's responders under the target; negatives its non-responders (within) or all donors
# of the J-matched null scenario (vs_null).
detection_data <- function(scenario, comparison, det, target) {
  d <- donor_by_detector[[det]] |>
    mutate(y = .data[[paste0("responder_", target)]])
  mine <- d |> filter(scenario_id == scenario)
  if (comparison == "within") {
    return(mine |> transmute(seed, value, y))
  }
  null_id <- paste0(
    "null_J",
    scenarios$J[match(scenario, scenarios$scenario_id)]
  )
  bind_rows(
    mine |> filter(y) |> transmute(seed, value, y = TRUE),
    d |> filter(scenario_id == null_id) |> transmute(seed, value, y = FALSE)
  )
}
detection <- expand_grid(
  scenario_id = one_visit_ids,
  target = TARGETS,
  comparison = COMPARISONS,
  detector = names(donor_by_detector)
) |>
  filter(!(str_starts(scenario_id, "null_") & comparison == "vs_null")) |>
  mutate(
    d = pmap(list(scenario_id, comparison, detector, target), detection_data)
  )
auc_table <- detection |>
  mutate(res = map(d, \(d) if (nrow(d) == 0) NULL else boot_stat(d, auc))) |>
  select(-d) |>
  unnest(res) |>
  rename(auc = estimate)
auc_replicate <- detection |>
  mutate(
    res = map(d, \(d) {
      if (nrow(d) == 0) NULL else summarise(d, auc = auc(value, y), .by = seed)
    })
  ) |>
  select(-d) |>
  unnest(res)
roc_curves <- detection |>
  filter(map_lgl(d, \(d) sum(d$y) > 0 && sum(!d$y) > 0)) |>
  mutate(
    roc = map(d, \(d) {
      d |>
        summarise(pos = sum(y), neg = sum(!y), .by = value) |>
        arrange(desc(value)) |>
        transmute(tpr = cumsum(pos) / sum(pos), fpr = cumsum(neg) / sum(neg)) |>
        add_row(tpr = 0, fpr = 0, .before = 1)
    })
  ) |>
  select(-d) |>
  unnest(roc)
# Conventional rules: sensitivity among responders, specificity among non-responders of the scenario.
rules <- RULES |>
  expand_grid(
    scenario_id = one_visit_ids,
    target = TARGETS,
    rate = c("sensitivity", "specificity")
  ) |>
  filter(detector %in% names(donor_by_detector)) |>
  mutate(
    res = pmap(
      list(scenario_id, detector, threshold, rate, target),
      \(sc, det, thr, rate, target) {
        d <- detection_data(sc, "within", det, target)
        keep <- if (rate == "sensitivity") d$y else !d$y
        if (!any(keep)) {
          return(NULL)
        }
        boot_stat(d[keep, ], \(x, y) {
          mean(if (rate == "sensitivity") x > thr else x <= thr)
        })
      }
    )
  ) |>
  unnest(res)

# Magnitude --------------------------------------------------------------------
# Donor RMS of every method against both truths; only the candidates have an interval.
magnitude_table <- attach_scores(
  truth_donor,
  read_methods("donor_magnitude.tsv") |> rename(value = estimate),
  c("scenario_id", "seed", "donor_id", "method")
)
magnitude_replicate <- TARGETS |>
  map(\(target) {
    magnitude_table |>
      mutate(truth = .data[[paste0("true_rms_", target)]]) |>
      summarise(
        bias = mean(value - truth),
        rmse = sqrt(mean((value - truth)^2)),
        coverage = mean(lower <= truth & truth <= upper),
        .by = c(scenario_id, seed, method)
      )
  }) |>
  list_rbind(names_to = "target")
magnitude <- magnitude_replicate |>
  summarise(
    n_seeds = n(),
    bias_se = sd(bias) / sqrt(n()), # SEs first: the means below overwrite bias, rmse and coverage
    rmse_se = sd(rmse^2) / sqrt(n()) / (2 * sqrt(mean(rmse^2))), # delta method on the mean MSE
    coverage_se = sd(coverage) / sqrt(n()),
    bias = mean(bias),
    rmse = sqrt(mean(rmse^2)),
    coverage = mean(coverage),
    .by = c(scenario_id, target, method)
  ) |>
  relocate(bias_se, rmse_se, coverage_se, .after = coverage)

# Pattern-level detection ------------------------------------------------------
# AUC over donor x pattern; COMPASS patterns removed by its filter score 0 (already so in 06).
pattern_table <- attach_scores(
  truth_pattern |>
    select(
      scenario_id,
      seed,
      donor_id,
      pattern_key,
      responding_sample,
      responding_struct
    ),
  read_methods("pattern_scores.tsv"),
  c("scenario_id", "seed", "donor_id", "pattern_key", "method", "score")
) |>
  mutate(detector = paste(method, score))
pattern_auc <- pattern_table |>
  group_split(scenario_id, detector) |>
  map(\(d) {
    TARGETS |>
      map(\(target) {
        boot_stat(
          transmute(d, seed, value, y = .data[[paste0("responding_", target)]]),
          auc
        )
      }) |>
      list_rbind(names_to = "target") |>
      mutate(
        scenario_id = d$scenario_id[1],
        detector = d$detector[1],
        .before = 1
      )
  }) |>
  list_rbind() |>
  rename(auc = estimate)

# Runtime and diagnostics ------------------------------------------------------
runtime_step <- bind_rows(
  read_reps("runtime.tsv"),
  read_reps("runtime.tsv", "polyicsflow", "done_b")
)
runtime_replicate <- runtime_step |>
  summarise(
    seconds = sum(seconds),
    cpus = max(cpus),
    .by = c(scenario_id, seed, method)
  )
runtime <- runtime_replicate |>
  summarise(
    n = n(),
    cpus = max(cpus),
    median_s = median(seconds),
    q25_s = quantile(seconds, 0.25),
    q75_s = quantile(seconds, 0.75),
    min_s = min(seconds),
    max_s = max(seconds),
    .by = c(scenario_id, method)
  )
diagnostics <- read_reps("diagnostics.tsv") |>
  mutate(
    flagged = max_rhat > RHAT_MAX | min_ess_bulk < ESS_MIN | n_divergent > 0
  ) |>
  summarise(
    n_replicates = n(),
    n_flagged = sum(flagged),
    max_rhat_median = median(max_rhat),
    max_rhat_max = max(max_rhat),
    min_ess_bulk_min = min(min_ess_bulk),
    n_divergent_total = sum(n_divergent),
    n_max_treedepth_total = sum(n_max_treedepth),
    phi_median = median(phi_mean),
    fit_seconds_median = median(fit_seconds),
    .by = c(scenario_id, model)
  )

# Write outputs ----------------------------------------------------------------
dir.create(out_dir, recursive = TRUE)
write_tsv(
  scenarios |>
    select(scenario_id, scenario, generator, J, role, n_donors, n_visits),
  file.path(out_dir, "scenario_roles.tsv")
)
write_tsv(
  replicate_status |> select(-rep_dir),
  file.path(out_dir, "replicate_status.tsv")
)
write_tsv(model_status, file.path(out_dir, "model_status.tsv"))
write_tsv(truth_summary, file.path(out_dir, "truth_summary.tsv"))
write_tsv(
  calibration_replicate,
  file.path(out_dir, "calibration_replicate.tsv")
)
write_tsv(calibration, file.path(out_dir, "calibration.tsv"))
write_tsv(false_sign_replicate, file.path(out_dir, "false_sign_replicate.tsv"))
write_tsv(false_sign, file.path(out_dir, "false_sign.tsv"))
write_tsv(ipf, file.path(out_dir, "ipf.tsv"))
write_tsv(change_replicate, file.path(out_dir, "change_replicate.tsv"))
write_tsv(change, file.path(out_dir, "change.tsv"))
write_tsv(donor_table, file.path(out_dir, "donor_table.tsv"))
write_tsv(roc_curves, file.path(out_dir, "roc_curves.tsv"))
write_tsv(auc_table, file.path(out_dir, "auc.tsv"))
write_tsv(auc_replicate, file.path(out_dir, "auc_replicate.tsv"))
write_tsv(rules, file.path(out_dir, "rules.tsv"))
write_tsv(magnitude_replicate, file.path(out_dir, "magnitude_replicate.tsv"))
write_tsv(magnitude, file.path(out_dir, "magnitude.tsv"))
write_tsv(pattern_auc, file.path(out_dir, "pattern_auc.tsv"))
write_tsv(runtime_replicate, file.path(out_dir, "runtime_replicate.tsv"))
write_tsv(runtime, file.path(out_dir, "runtime.tsv"))
write_tsv(diagnostics, file.path(out_dir, "diagnostics.tsv"))
write_session_info(file.path(out_dir, "session_info_summary.txt"), PKG_DIR)
