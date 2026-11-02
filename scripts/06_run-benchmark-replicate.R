#!/usr/bin/env Rscript
# Purpose: one redesign-evaluation replicate in environment A: simulate one scenario at one seed, write the
#          counts and both truths (realised sample p and structured p) for environment B (07), then fit
#          each candidate model (C0 conjugate Dirichlet, C1 m1, C2 m1dm, C2phi) and COMPASS (default
#          filter and no filter) on the identical counts. Candidate scores are written with both truths.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/06_run-benchmark-replicate.R <run_dir> <scenario_id> <seed>
# Inputs:
#   <run_dir>/scenarios.tsv     scenario parameters [05_define-benchmark.R]
#   <PKG_DIR>/                  package source, loaded with pkgload; env CJ_PKG_DIR overrides it [this repo / worktree]
#   scripts/helper01_true-scores.R, helper02_benchmark-sim.R   true scores, simulators [this repo]
# Outputs (in <run_dir>/reps/<scenario_id>/seed_<seed>/):
#   counts.tsv                  long counts: sample_id, donor_id, condition, N, pattern_key, n  (input of 07)
#   pattern_map.tsv             pattern_key, degree, and each method's name for the pattern
#   truth_donor.tsv             contrast_id, donor_id, true RMS for p_sample and p_struct, responder flags (> eps_response)
#   truth_pattern.tsv           contrast_id, donor_id, pattern_key, degree, delta and responding flag for both truths
#   model_status.tsv            per candidate: ok, failed (with the message) or not_run (with the reason)
#   model_scores.tsv            model, draws (sample, or structured for overdispersed models), contrast, grain (global, marker, pattern), score, marker, pattern_key,
#                               estimate, lower, upper, p_positive, p_negative, failure_fraction, status,
#                               truth_sample, truth_struct: every metric of the plan is computed from this
#   model_change.tsv            two-visit scenarios only: visit change v1 -> v2, donor and population, both truths
#   donor_scores.tsv            method, score, donor_id, value: detection scores (one-visit scenarios)
#   donor_magnitude.tsv         method, donor_id, estimate, lower, upper: estimated donor RMS (one-visit scenarios)
#   pattern_scores.tsv          method, score, donor_id, pattern_key, value, removed: pattern-level scores (one-visit)
#   compass_status.tsv          per COMPASS variant: ok, failed with the message, or not_run (two-visit scenarios)
#   runtime.tsv                 method, step, seconds, cpus
#   diagnostics.tsv             per candidate: max R-hat, min ESS, divergences, treedepth hits, E-BFMI, phi; written
#                               last, so its presence marks a complete replicate
#   session_info.txt            R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign") # package source under test; the redesign branch holds C0-C2
CJ_CPUS <- 4L # cj_fit runs 4 chains in parallel; the sbatch allocation is 4 CPUs
PHI_BY_CONDITION_ARG <- "phi_by_condition" # cj_fit argument that gives one phi per condition (C2phi); not run if absent
DRAWS_C0 <- 4000L # draws of the conjugate baseline, as many as 4 chains x 1000 of the Stan fits
A0_C0 <- 0.5 # Dirichlet prior pseudo-count per pattern of C0 (plan: Dir(1/2 + y))
EPSILON <- c(
  RMS = 1e-4,
  poly_absolute = 1e-4,
  poly_marginal = 1e-4,
  CPS = 0.01,
  CCS = 1e-4,
  breadth = 0.01,
  marker_magnitude = 1e-4
) # cj_scores() wants all of them; RMS is set to eps_response so P(RMS > eps) asks the truth's question
COMPASS_VARIANTS <- c("COMPASS_default", "COMPASS_nofilter") # default category filter, and category_filter = NULL
MARKERS <- c("IFNG", "TNF", "IL2", "CD107a")

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 06_run-benchmark-replicate.R <run_dir> <scenario_id> <seed>" = length(
    args
  ) ==
    3
)
run_dir <- args[1]
seed <- as.integer(args[3])
sc <- read_tsv(file.path(run_dir, "scenarios.tsv"), show_col_types = FALSE) |>
  filter(scenario_id == args[2])
stopifnot("scenario_id not in scenarios.tsv" = nrow(sc) == 1)
rep_dir <- file.path(
  run_dir,
  "reps",
  sc$scenario_id,
  sprintf("seed_%02d", seed)
)
stopifnot("replicate directory already exists" = !dir.exists(rep_dir))
partial_dir <- paste0(rep_dir, ".partial") # outputs land here and move to rep_dir when complete
unlink(partial_dir, recursive = TRUE) # a stale one comes from a killed task
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper01_true-scores.R") # true_table(), write_session_info()
source("scripts/helper02_benchmark-sim.R") # bench_simulate(), true_pattern_delta()
data_seed <- sc$seed_offset + seed
one_visit <- sc$n_visits == 1 # donor-level detection and COMPASS / polyICSFlow are defined for one visit only
pattern_epsilon <- c(RMS = sc$eps_pattern, CPS = 0.01, CCS = sc$eps_pattern)

# Simulate and truth -----------------------------------------------------------
sim <- bench_simulate(sc, data_seed)
bits <- sim$data$bits
pattern_map <- tibble(
  pattern_key = apply(bits, 1, paste, collapse = ""),
  degree = rowSums(bits),
  cj_label = colnames(sim$data$counts),
  compass_name = compass_pattern_names(bits),
  polyics_name = apply(bits, 1, \(b) {
    paste0(MARKERS[seq_len(sc$J)], if_else(b == 1, "+", "-"), collapse = "")
  })
)
counts <- sim$data$counts |>
  as_tibble(rownames = "sample_id") |>
  pivot_longer(-sample_id, names_to = "cj_label", values_to = "n") |>
  left_join(pattern_map, join_by(cj_label), relationship = "many-to-one") |>
  left_join(
    sim$data$samples |> select(sample_id, donor_id, condition, N),
    join_by(sample_id),
    relationship = "many-to-one"
  ) |>
  select(sample_id, donor_id, condition, N, pattern_key, n)
truth <- true_table(sim)
truth_pattern <- true_pattern_delta(sim, "p") |>
  rename(delta_sample = delta) |>
  left_join(
    true_pattern_delta(sim, "p_struct") |>
      select(contrast_id, pattern_key, delta_struct = delta),
    join_by(contrast_id, pattern_key),
    relationship = "one-to-one"
  ) |>
  mutate(
    responding_sample = delta_sample > sc$eps_pattern,
    responding_struct = delta_struct > sc$eps_pattern
  )
truth_donor <- truth |>
  filter(grain == "global", score == "RMS") |>
  transmute(
    contrast_id,
    donor_id,
    true_rms_sample = truth_sample,
    true_rms_struct = truth_struct,
    responder_sample = true_rms_sample > sc$eps_response,
    responder_struct = true_rms_struct > sc$eps_response
  )
stopifnot(
  "pattern changes sum to true RMS" = isTRUE(all.equal(
    truth_pattern |>
      summarise(rms = sum(delta_sample), .by = contrast_id) |>
      arrange(contrast_id) |>
      pull(rms),
    truth_donor |> arrange(contrast_id) |> pull(true_rms_sample)
  ))
)

# Candidate models -------------------------------------------------------------
# C0 conjugate Dirichlet, C1 current M1, C2 M1 mean + Dirichlet-multinomial residual, C2phi one phi per condition.
# A candidate the loaded package cannot fit is recorded as not_run; a fit that errors is recorded as failed.
cj_formals <- names(formals(cj_fit))
fit_model <- function(model) {
  switch(
    model,
    C0 = cytoJointR:::cj_fit_dirichlet(
      sim$data,
      a0 = A0_C0,
      n_draws = DRAWS_C0,
      seed = data_seed
    ),
    C1 = if ("model" %in% cj_formals) {
      cj_fit(sim$data, model = "m1", level = "sample", seed = data_seed)
    } else {
      cj_fit(sim$data, seed = data_seed)
    },
    C2 = cj_fit(sim$data, model = "m1dm", level = "sample", seed = data_seed),
    C2phi = do.call(
      cj_fit,
      c(
        list(sim$data, model = "m1dm", level = "sample", seed = data_seed),
        setNames(list(TRUE), PHI_BY_CONDITION_ARG)
      )
    )
  )
}
not_run_reason <- c(
  C0 = if (!exists("cj_fit_dirichlet", asNamespace("cytoJointR"))) {
    "package has no cj_fit_dirichlet()"
  },
  C2 = if (!"model" %in% cj_formals) {
    "package has no model argument in cj_fit()"
  },
  C2phi = if (!PHI_BY_CONDITION_ARG %in% cj_formals) {
    paste0("cj_fit() has no ", PHI_BY_CONDITION_ARG, " argument")
  }
)
fits <- c("C0", "C1", "C2", "C2phi") |>
  set_names() |>
  map(\(model) {
    if (model %in% names(not_run_reason)) {
      return(list(status = "not_run", message = not_run_reason[[model]]))
    }
    tryCatch(
      {
        seconds <- system.time(fit <- fit_model(model))[["elapsed"]]
        list(
          status = "ok",
          message = NA_character_,
          fit = fit,
          seconds = seconds
        )
      },
      error = \(e) list(status = "failed", message = conditionMessage(e))
    )
  })
model_status <- fits |>
  map(\(f) tibble(status = f$status, message = f$message)) |>
  list_rbind(names_to = "model")
fits_ok <- fits |> keep(\(f) f$status == "ok")
diagnostics <- fits_ok |>
  imap(\(f, model) {
    dg <- if (inherits(f$fit, "cj_fit")) f$fit$diagnostics
    phi_row <- if (!is.null(dg)) {
      which(startsWith(as.character(dg$summary$variable), "phi"))[1]
    }
    tibble(
      max_rhat = dg$max_rhat %||% NA_real_,
      min_ess_bulk = dg$min_ess_bulk %||% NA_real_,
      min_ess_tail = dg$min_ess_tail %||% NA_real_,
      n_divergent = dg$n_divergent %||% NA_integer_,
      n_max_treedepth = dg$n_max_treedepth %||% NA_integer_,
      min_ebfmi = if (is.null(dg)) NA_real_ else min(dg$ebfmi),
      phi_mean = if (is.null(dg)) NA_real_ else dg$summary$mean[phi_row], # NA for models without a phi
      fit_seconds = f$seconds
    )
  }) |>
  list_rbind(names_to = "model")

# Score the candidates ---------------------------------------------------------
# Global, marker and exact-pattern scores of every contrast, joined to both truths. The truth table leads, so
# a contrast the model withheld (failure tolerance) still has a row with NA estimates. cj_patterns() runs first:
# its engine pass also serves cj_scores(). Each model is scored on its sample-level draws (draws = "sample", the
# cj_draws estimand); an overdispersed model is also scored on its structured draws (draws = "structured").
score_draws <- function(x) {
  pat <- cj_patterns(x, c("RMS", "CPS", "CCS"), epsilon = pattern_epsilon)
  glob <- cj_scores(x, epsilon = EPSILON)
  keep_cols <- c(
    "estimate",
    "lower",
    "upper",
    "p_positive",
    "p_negative",
    "failure_fraction",
    "status"
  )
  bind_rows(
    glob |>
      transmute(
        contrast_id,
        grain = if_else(is.na(marker), "global", "marker"),
        score,
        marker = as.character(marker),
        pattern_key = NA_character_,
        across(all_of(keep_cols))
      ),
    pat |>
      left_join(
        pattern_map |> select(cj_label, pattern_key),
        join_by(pattern == cj_label),
        relationship = "many-to-one"
      ) |>
      transmute(
        contrast_id,
        grain = "pattern",
        score,
        marker = NA_character_,
        pattern_key,
        across(all_of(keep_cols))
      )
  ) |>
    right_join(
      truth,
      join_by(contrast_id, grain, score, marker, pattern_key),
      relationship = "one-to-one"
    )
}
draws_levels <- function(fit) {
  if (!inherits(fit, "cj_fit")) {
    return(list(sample = fit)) # C0 is a cj_draws
  }
  out <- list(sample = fit$draws_sample %||% fit$draws) # an older package has only $draws
  if (
    !is.null(fit$draws_structured) &&
      !identical(fit$draws_structured, fit$draws_sample)
  ) {
    out$structured <- fit$draws_structured
  }
  out
}
scored <- fits_ok |>
  map(\(f) {
    draws_levels(f$fit) |>
      map(\(x) {
        seconds <- system.time(table <- score_draws(x))[["elapsed"]]
        list(table = table, seconds = seconds)
      })
  })
model_scores <- scored |>
  map(\(per_model) {
    per_model |> map("table") |> list_rbind(names_to = "draws")
  }) |>
  list_rbind(names_to = "model")
score_seconds <- scored |>
  map_dbl(\(per_model) sum(map_dbl(per_model, "seconds")))

# Visit change -----------------------------------------------------------------
# Two-visit scenarios: v1 -> v2 change per donor and for the equal-donor average, against both truths.
model_change <- tibble()
if (!one_visit) {
  change_truth <- truth |>
    filter(grain != "pattern") |>
    select(
      donor_id,
      stimulation,
      visit,
      score,
      marker,
      truth_sample,
      truth_struct
    ) |>
    pivot_wider(
      names_from = visit,
      values_from = c(truth_sample, truth_struct)
    ) |>
    transmute(
      donor_id,
      stimulation,
      score,
      marker,
      truth_sample = truth_sample_v2 - truth_sample_v1,
      truth_struct = truth_struct_v2 - truth_struct_v1
    )
  change_truth_pop <- change_truth |>
    summarise(
      across(c(truth_sample, truth_struct), mean),
      .by = c(stimulation, score, marker)
    )
  keep_cols <- c(
    "score",
    "marker",
    "estimate",
    "lower",
    "upper",
    "p_positive",
    "p_negative",
    "failure_fraction",
    "status"
  )
  model_change <- fits_ok |>
    map(\(f) {
      vc <- cj_visit_change(f$fit, epsilon = EPSILON, from = "v1", to = "v2")
      bind_rows(
        vc$donor |>
          select(donor_id, stimulation, all_of(keep_cols)) |>
          left_join(
            change_truth,
            join_by(donor_id, stimulation, score, marker),
            relationship = "one-to-one"
          ) |>
          mutate(grain = "donor"),
        vc$population |>
          select(stimulation, all_of(keep_cols)) |>
          left_join(
            change_truth_pop,
            join_by(stimulation, score, marker),
            relationship = "one-to-one"
          ) |>
          mutate(grain = "population")
      )
    }) |>
    list_rbind(names_to = "model")
}

# COMPASS ----------------------------------------------------------------------
# One run per variant, inside tryCatch so a COMPASS error is recorded and the cytoJointR outputs stay.
# SimpleCOMPASS() does not compute the posterior (COMPASS() does, at COMPASS.R:470), so
# PosteriorDiff() needs compute_posterior() first; it samples from the stored fit and leaves FS and PFS untouched.
# Pattern scores are mean_gamma; patterns the filter removed score 0 and carry removed = TRUE.
positive_names <- pattern_map |> filter(degree >= 1) |> pull(compass_name)
run_compass <- function(variant) {
  filter_fn <- list(
    COMPASS_default = COMPASS_DEFAULT_FILTER,
    COMPASS_nofilter = NULL
  )[[variant]]
  tryCatch(
    {
      seconds <- system.time(
        cc <- cj_compass(
          sim$data,
          category_filter = filter_fn,
          seed = data_seed
        )
      )[["elapsed"]]
      fit_c <- cc$fits[[1]]
      stopifnot("COMPASS did not fit" = !is.null(fit_c))
      set.seed(data_seed)
      fit_c$fit$posterior <- COMPASS:::compute_posterior(fit_c)
      # Sum of the posterior-mean differences over the kept non-null patterns, fraction scale. This is
      # rowSums(COMPASS::PosteriorDiff(fit_c)); PosteriorDiff() itself errors with one kept pattern
      # (Posterior.R:57-59, sapply() returns a vector), so the per-donor vectors are summed directly.
      diff_sum <- map_dbl(fit_c$fit$posterior, \(p) sum(p$diff))
      gamma <- fit_c$fit$mean_gamma[, -ncol(fit_c$fit$mean_gamma), drop = FALSE] # drops the null category
      list(
        seconds = seconds,
        scores = cc$scores |>
          transmute(
            donor_id,
            FS = compass_FS,
            PFS = compass_PFS,
            Response = compass_response_degree1
          ) |>
          pivot_longer(-donor_id, names_to = "score", values_to = "value"),
        magnitude = tibble(
          donor_id = names(diff_sum),
          estimate = unname(diff_sum)
        ),
        pattern = expand_grid(
          donor_id = rownames(gamma),
          compass_name = positive_names
        ) |>
          left_join(
            as_tibble(gamma, rownames = "donor_id") |>
              pivot_longer(
                -donor_id,
                names_to = "compass_name",
                values_to = "value"
              ),
            join_by(donor_id, compass_name),
            relationship = "one-to-one"
          ) |>
          mutate(
            removed = is.na(value),
            value = coalesce(value, 0),
            score = "mean_gamma"
          ),
        status = tibble(
          status = "ok",
          message = NA_character_,
          patterns_kept = ncol(gamma)
        )
      )
    },
    error = \(e) {
      list(
        status = tibble(
          status = "failed",
          message = conditionMessage(e),
          patterns_kept = NA_integer_
        )
      )
    }
  )
}
compass <- if (one_visit) {
  set_names(COMPASS_VARIANTS) |> map(run_compass)
} else {
  list() # COMPASS and polyICSFlow give no visit change
}


# Collect ----------------------------------------------------------------------
# Long tables with one row per method, donor and score (and pattern); candidates come from model_scores.
bind_compass <- \(what) {
  compass |>
    map(\(x) x[[what]]) |>
    compact() |>
    list_rbind(names_to = "method")
}
cj_rms <- model_scores |>
  filter(draws == "sample", grain == "global", score == "RMS") |>
  transmute(method = model, donor_id, estimate, lower, upper, p_positive)
cj_pat <- model_scores |>
  filter(draws == "sample", grain == "pattern", score == "RMS") |>
  transmute(method = model, donor_id, pattern_key, p_positive)
donor_scores <- bind_rows(
  cj_rms |>
    transmute(method, score = "P_RMS_gt_eps", donor_id, value = p_positive),
  bind_compass("scores")
)
donor_magnitude <- bind_rows(
  cj_rms |> select(method, donor_id, estimate, lower, upper),
  bind_compass("magnitude")
)
pattern_scores <- bind_rows(
  cj_pat |>
    transmute(
      method,
      score = "P_contribution_gt_eps",
      donor_id,
      pattern_key,
      value = p_positive,
      removed = FALSE
    ),
  if (one_visit) {
    bind_compass("pattern") |>
      left_join(
        pattern_map |> select(pattern_key, compass_name),
        join_by(compass_name),
        relationship = "many-to-one"
      ) |>
      select(method, score, donor_id, pattern_key, value, removed)
  }
)
compass_status <- if (one_visit) {
  compass |>
    map(\(x) x$status) |>
    list_rbind(names_to = "method")
} else {
  tibble(
    method = COMPASS_VARIANTS,
    status = "not_run",
    message = "two-visit scenario",
    patterns_kept = NA_integer_
  )
}
runtime <- bind_rows(
  tibble(
    method = rep(names(fits_ok), 2),
    step = rep(c("fit", "scores_patterns"), each = length(fits_ok)),
    seconds = c(map_dbl(fits_ok, "seconds"), score_seconds),
    cpus = if_else(method == "C0", 1L, CJ_CPUS)
  ),
  compass |>
    map(\(x) x$seconds) |>
    compact() |>
    map(\(s) tibble(step = "cj_compass", seconds = s, cpus = 1L)) |>
    list_rbind(names_to = "method")
)

# Write outputs ----------------------------------------------------------------
# The fits took minutes, so everything is on disk before the rename; diagnostics.tsv goes last.
dir.create(partial_dir, recursive = TRUE)
write_tsv(counts, file.path(partial_dir, "counts.tsv"))
write_tsv(pattern_map, file.path(partial_dir, "pattern_map.tsv"))
write_tsv(truth_donor, file.path(partial_dir, "truth_donor.tsv"))
write_tsv(truth_pattern, file.path(partial_dir, "truth_pattern.tsv"))
write_tsv(model_status, file.path(partial_dir, "model_status.tsv"))
write_tsv(model_scores, file.path(partial_dir, "model_scores.tsv"))
write_tsv(model_change, file.path(partial_dir, "model_change.tsv"))
write_tsv(donor_scores, file.path(partial_dir, "donor_scores.tsv"))
write_tsv(donor_magnitude, file.path(partial_dir, "donor_magnitude.tsv"))
write_tsv(pattern_scores, file.path(partial_dir, "pattern_scores.tsv"))
write_tsv(compass_status, file.path(partial_dir, "compass_status.tsv"))
write_tsv(runtime, file.path(partial_dir, "runtime.tsv"))
write_session_info(file.path(partial_dir, "session_info.txt"), PKG_DIR)
write_tsv(diagnostics, file.path(partial_dir, "diagnostics.tsv"))
file.rename(partial_dir, rep_dir)
