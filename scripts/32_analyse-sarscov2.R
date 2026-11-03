#!/usr/bin/env Rscript
# Purpose: benchmark cytoJointR (m1dm, m1, C0), COMPASS and polyICSFlow on the public SARS-CoV-2 Spike ICS counts
#          (Zenodo 8088178) for one parent population. Stage "cytojointr" (environment A) fits the three
#          cytoJointR models and COMPASS, scores Spike and PHA against the shared DMSO control, computes the visit
#          change, the m1dm phi posterior, LOO (m1dm vs m1) and posterior predictive checks. Stage "polyics"
#          (environment B) runs polyICSFlow on the same counts.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/32_analyse-sarscov2.R <run_dir> cytojointr <CD4|CD8>
#          R_LIBS=renv/lib-r44 /exports/para-lipg-hpc/mdmanurung/batchdx-cyto/renv/bin/Rscript scripts/32_analyse-sarscov2.R <run_dir> polyics <CD4|CD8>
# Inputs:
#   data/processed/sarscov2/counts.tsv   long Boolean counts [24_prepare-sarscov2.R]
#   cytoJointR-wt/redesign-run/          package, frozen at commit 49278c9 [this repo]
#   scripts/helper05_realdata-sarscov2.R, helper01_true-scores.R (session info)   [this repo]
# Outputs (in <run_dir>/<pop>/, written by stage cytojointr; the directory must not exist):
#   scores.tsv          model, eps_set, cj_scores() rows per Spike and PHA contrast (global and marker views)
#   visit_change.tsv    model m1dm, grain (donor, population), cj_visit_change() rows for each visit pair
#   compass.tsv         donor_id, visit, stimulation, FS, PFS, Response (degree 1 posterior probability)
#   compass_groups.tsv  one row per COMPASS run (stimulation x visit): fitted or failed; compass_filters.tsv lists what the filter removed
#   diagnostics.tsv     model, max R-hat, min ESS, divergences, treedepth hits, E-BFMI, phi summary
#   phi_draws.tsv       draw, phi of m1dm
#   loo.tsv             model, elpd_loo, se, p_loo, number of Pareto k > 0.7; loo_compare.tsv; loo_pointwise.tsv
#   ppc.tsv             model, sample_id, statistic, observed, replicated mean and 90% interval, mid-p
#   samples.tsv         sample_id, donor_id, visit, condition, N of the fitted cj_data
#   runtime.tsv         method, step, seconds, cpus
#   session_info.txt    package commit and sessionInfo()
# Outputs (stage polyics, in <run_dir>/<pop>/polyicsflow/):
#   donor_scores.tsv    donor_id, visit, stimulation, clipped and unclipped net sum over positive patterns (fraction)
#   runtime.tsv         step, seconds, cpus; written last
#   session_info.txt    polyICSFlow commit and sessionInfo()

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
SEED <- 20261103L # every stochastic step derives from it: Stan, C0 draws, COMPASS, PPC
EPSILON <- c(
  RMS = 1e-4,
  poly_absolute = 1e-4,
  poly_marginal = 1e-4,
  CPS = 0.01,
  CCS = 1e-4,
  breadth = 0.01,
  marker_magnitude = 1e-4
) # as scripts/06; RMS 1e-4 is eps_response of 05_define-benchmark.R (one cell in 10 000)
EPS_SENS_FACTOR <- 10 # sensitivity set: all epsilons x10, because CD8 parents can be as small as 400 cells
VISIT_PAIRS <- list(c("E01", "E02"), c("E01", "E03"), c("E02", "E03")) # first two requested; E02 -> E03 is the dose-2 step
DRAWS_C0 <- 4000L # as many as 4 chains x 1000 of the Stan fits
A0_C0 <- 0.5 # Dirichlet pseudo-count of C0, as scripts/06
CJ_CPUS <- 4L # cj_fit runs 4 chains in parallel
N_PPC <- 400L # posterior draws replicated per model; enough for 90% intervals of counts
RARE_DEGREE <- 3L # patterns with at least this many positive markers are "rare"
PKG_DIR <- "cytoJointR-wt/redesign-run" # frozen worktree, commit 49278c9
BACKGROUND <- "control" # condition value of the DMSO sample
STIM_LABELS <- c(peptide = "Spike", PHA = "PHA") # antigen in the workbook -> stimulation in the fit

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: 32_analyse-sarscov2.R <run_dir> <cytojointr|polyics> <CD4|CD8>" = length(
    args
  ) ==
    3
)
run_dir <- args[1]
stage <- args[2]
pop <- args[3]
pop_dir <- file.path(run_dir, pop)
counts <- read_tsv(
  "data/processed/sarscov2/counts.tsv",
  col_types = cols(pattern_key = "c", donor_id = "c")
) |>
  filter(population == pop)
stopifnot(
  "stage must be cytojointr or polyics" = stage %in% c("cytojointr", "polyics")
)

# Stage polyics ----------------------------------------------------------------
# Environment B. Counts expand to one row per cell, as scripts/07; net = stimulation minus the DMSO of the same
# donor and visit. Clipped is perc_pos_backgroundSubtracted; unclipped is perc_pos minus the control's perc_pos.
if (stage == "polyics") {
  suppressPackageStartupMessages(library(polyICSFlow))
  out_dir <- file.path(pop_dir, "polyicsflow")
  stopifnot(
    "stage cytojointr has not finished" = file.exists(file.path(
      pop_dir,
      "runtime.tsv"
    )),
    "polyicsflow output already exists" = !dir.exists(out_dir)
  )
  markers <- c("CD40L", "IFNg", "IL2", "TNFa")
  patterns <- tibble(pattern_key = unique(counts$pattern_key)) |>
    mutate(
      degree = nchar(gsub("0", "", pattern_key)),
      polyics_name = map_chr(pattern_key, \(k) {
        paste0(
          markers,
          if_else(strsplit(k, "")[[1]] == "1", "+", "-"),
          collapse = ""
        )
      })
    )
  expand_seconds <- system.time({
    cells <- counts |>
      mutate(
        condition = if_else(
          antigen == "DMSO",
          BACKGROUND,
          unname(STIM_LABELS[antigen])
        )
      ) |>
      filter(n > 0) |>
      uncount(n)
    key_bits <- patterns$pattern_key |>
      strsplit("") |>
      map(\(b) b == "1") |>
      reduce(rbind)
    colnames(key_bits) <- markers
    bits <- as_tibble(key_bits[
      match(cells$pattern_key, patterns$pattern_key),
    ])
  })[["elapsed"]]
  assign_seconds <- system.time(
    comb <- assignMarkerCombinations(bits, mode = "simple") |>
      mutate(
        donor_id = cells$donor_id,
        visit = cells$visit,
        condition = cells$condition,
        cell_type = pop
      )
  )[["elapsed"]]
  calc_seconds <- system.time(
    res <- calcPolyfunctionality(
      comb,
      md_cols = c("donor_id", "visit"),
      pop_col = "cell_type",
      resolution = "MarkerComb",
      condition_col = "condition",
      background_val = BACKGROUND
    )
  )[["elapsed"]]
  net <- res |>
    filter(condition == BACKGROUND) |>
    select(donor_id, visit, MarkerComb, perc_background = perc_pos) |>
    inner_join(
      res |> filter(condition != BACKGROUND),
      join_by(donor_id, visit, MarkerComb),
      relationship = "one-to-many"
    ) |>
    mutate(polyics_name = as.character(MarkerComb)) |>
    inner_join(patterns, join_by(polyics_name), relationship = "many-to-one") |>
    filter(degree >= 1) |>
    mutate(stimulation = condition)
  donor_scores <- net |>
    summarise(
      clipped = sum(perc_pos_backgroundSubtracted) / 100,
      unclipped = sum(perc_pos - perc_background) / 100,
      n_patterns = n(),
      .by = c(donor_id, visit, stimulation)
    )
  dir.create(out_dir, recursive = TRUE)
  write_tsv(donor_scores, file.path(out_dir, "donor_scores.tsv"))
  writeLines(
    c(
      paste(
        "polyICSFlow",
        as.character(packageVersion("polyICSFlow")),
        packageDescription("polyICSFlow")$RemoteSha
      ),
      capture.output(sessionInfo())
    ),
    file.path(out_dir, "session_info.txt")
  )
  write_tsv(
    tibble(
      method = "polyICSFlow",
      step = c(
        "expand_counts_to_cells",
        "assignMarkerCombinations",
        "calcPolyfunctionality"
      ),
      seconds = c(expand_seconds, assign_seconds, calc_seconds),
      cpus = 1L
    ),
    file.path(out_dir, "runtime.tsv")
  )
  quit(save = "no")
}

# Stage cytojointr: set-up -----------------------------------------------------
stopifnot(
  "run directory for this population already exists" = !dir.exists(pop_dir)
)
partial_dir <- paste0(pop_dir, ".partial")
unlink(partial_dir, recursive = TRUE) # a stale one comes from a killed job
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper01_true-scores.R") # write_session_info()
source("scripts/helper05_realdata-sarscov2.R") # ics_to_cj(), ppc_replicate()
x <- ics_to_cj(counts, pop)
sample_info <- x$samples |> select(sample_id, donor_id, visit, condition, N)
# The fitted count matrix and the sample table must list the same samples in the same order.
stopifnot(
  "sample table matches count matrix" = identical(
    sample_info$sample_id,
    rownames(x$counts)
  )
)
degree <- rowSums(x$bits)

# Fit --------------------------------------------------------------------------
# m1dm (default) and m1 by cmdstan; C0 is the conjugate Dirichlet baseline. Each fit is timed.
timed <- \(expr) {
  t0 <- proc.time()[["elapsed"]]
  value <- expr
  list(fit = value, seconds = proc.time()[["elapsed"]] - t0)
}
fits <- list(
  m1dm = timed(cj_fit(x, model = "m1dm", seed = SEED)),
  m1 = timed(cj_fit(x, model = "m1", seed = SEED)),
  C0 = timed(cytoJointR:::cj_fit_dirichlet(
    x,
    a0 = A0_C0,
    n_draws = DRAWS_C0,
    seed = SEED
  ))
)
draws <- list(
  m1dm = fits$m1dm$fit$draws_sample,
  m1 = fits$m1$fit$draws_sample,
  C0 = fits$C0$fit
)

# Scores -----------------------------------------------------------------------
# Every score of every Spike and PHA contrast, at the main epsilon and at the x10 sensitivity set.
eps_sets <- list(main = EPSILON, sens = EPSILON * EPS_SENS_FACTOR)
score_runs <- map(draws, \(d) map(eps_sets, \(e) timed(cj_scores(d, e))))
scores <- score_runs |>
  map(\(per_model) map(per_model, "fit") |> list_rbind(names_to = "eps_set")) |>
  list_rbind(names_to = "model")

# Visit change -----------------------------------------------------------------
# m1dm only: to minus from, per donor and for the equal-donor average. The paper defines the visits (24_prepare).
visit_runs <- map(VISIT_PAIRS, \(v) {
  timed(cj_visit_change(fits$m1dm$fit, EPSILON, from = v[1], to = v[2]))
})
visit_change <- visit_runs |>
  map(\(r) {
    bind_rows(
      r$fit$donor |> mutate(grain = "donor"),
      r$fit$population |> mutate(grain = "population")
    )
  }) |>
  list_rbind()

# COMPASS ----------------------------------------------------------------------
# One COMPASS run per stimulation and visit (default category filter); Response is the degree-1 posterior probability.
compass_run <- timed(cj_compass(x, seed = SEED))
compass <- compass_run$fit$scores |>
  transmute(
    donor_id,
    visit,
    stimulation,
    FS = compass_FS,
    PFS = compass_PFS,
    Response = compass_response_degree1
  )

# Diagnostics and phi ----------------------------------------------------------
diagnostics <- c("m1dm", "m1") |>
  set_names() |>
  map(\(m) {
    dg <- fits[[m]]$fit$diagnostics
    phi <- dg$summary |> filter(startsWith(as.character(variable), "phi"))
    tibble(
      max_rhat = dg$max_rhat,
      min_ess_bulk = dg$min_ess_bulk,
      min_ess_tail = dg$min_ess_tail,
      n_divergent = dg$n_divergent,
      n_max_treedepth = dg$n_max_treedepth,
      min_ebfmi = min(dg$ebfmi),
      phi_mean = phi$mean[1] %||% NA_real_,
      phi_q5 = phi$q5[1] %||% NA_real_,
      phi_q95 = phi$q95[1] %||% NA_real_,
      phi_rhat = phi$rhat[1] %||% NA_real_
    )
  }) |>
  list_rbind(names_to = "model")
phi_draws <- tibble(
  phi = as.vector(posterior::as_draws_matrix(fits$m1dm$fit$engine_fit$draws(
    "phi"
  )))
) |>
  mutate(draw = row_number(), .before = phi)

# LOO --------------------------------------------------------------------------
# PSIS-LOO on the pointwise log_lik of each sample (leave one sample out). r_eff comes from the likelihood scaled
# by its column maximum: r_eff is scale free, and the raw likelihood of a 100 000-cell sample underflows to 0.
loo_seconds <- system.time({
  loos <- c("m1dm", "m1") |>
    set_names() |>
    map(\(m) {
      ll <- fits[[m]]$fit$engine_fit$draws("log_lik")
      ll <- array(as.numeric(ll), dim(ll))
      r_eff <- loo::relative_eff(exp(sweep(ll, 3, apply(ll, 3, max))))
      loo::loo(ll, r_eff = r_eff)
    })
})[["elapsed"]]
loo_table <- loos |>
  imap(\(l, m) {
    tibble(
      elpd_loo = l$estimates["elpd_loo", "Estimate"],
      se_elpd_loo = l$estimates["elpd_loo", "SE"],
      p_loo = l$estimates["p_loo", "Estimate"],
      n_pareto_k_gt_0.7 = sum(l$diagnostics$pareto_k > 0.7)
    )
  }) |>
  list_rbind(names_to = "model")
cmp <- loo::loo_compare(loos) # first row is the better model; elpd_diff is relative to it
loo_compare <- tibble(
  model = loo_table$model[order(-loo_table$elpd_loo)], # loo_compare() sorts by elpd, best first; its rownames are positions
  elpd_diff = cmp[, "elpd_diff"],
  se_diff = cmp[, "se_diff"]
)
loo_pointwise <- loos |>
  imap(\(l, m) {
    tibble(
      sample_id = rownames(x$counts),
      elpd = l$pointwise[, "elpd_loo"],
      pareto_k = l$diagnostics$pareto_k
    )
  }) |>
  list_rbind(names_to = "model")

# Posterior predictive checks --------------------------------------------------
# Replicate every sample from the structured probabilities mu (m1dm: Dirichlet-multinomial with phi; m1: multinomial)
# for N_PPC posterior draws. Statistics per sample: cells positive for any marker, cells in rare patterns
# (>= RARE_DEGREE positive markers) and cells positive for all four markers. These are the same samples the model
# was fitted to, so this is an in-sample check. Mid-p: P(rep < obs) + P(rep == obs) / 2.
set.seed(SEED)
ppc_stat <- \(y) {
  cbind(
    any_positive = rowSums(y[, degree >= 1, drop = FALSE]),
    rare_patterns = rowSums(y[, degree >= RARE_DEGREE, drop = FALSE]),
    all_four = rowSums(y[, degree == max(degree), drop = FALSE])
  )
}
observed <- ppc_stat(x$counts)
ppc_seconds <- system.time({
  ppc <- c("m1dm", "m1") |>
    set_names() |>
    map(\(m) {
      mu <- fits[[m]]$fit$draws_structured$p
      phi <- if (m == "m1dm") {
        as.vector(posterior::as_draws_matrix(fits[[m]]$fit$engine_fit$draws(
          "phi"
        )))
      }
      idx <- sample(dim(mu)[1], N_PPC)
      reps <- map(idx, \(d) {
        ppc_stat(ppc_replicate(mu[d, , ], phi[d], sample_info$N))
      }) |>
        simplify2array() # sample x statistic x draw
      expand_grid(
        sample_id = rownames(x$counts),
        statistic = colnames(observed)
      ) |>
        mutate(
          observed = map2_dbl(sample_id, statistic, \(s, k) observed[s, k]),
          rep = map2(sample_id, statistic, \(s, k) {
            reps[match(s, rownames(x$counts)), k, ]
          }),
          rep_mean = map_dbl(rep, mean),
          rep_lower = map_dbl(rep, quantile, 0.05),
          rep_upper = map_dbl(rep, quantile, 0.95),
          mid_p = map2_dbl(rep, observed, \(r, o) {
            mean(r < o) + mean(r == o) / 2
          })
        ) |>
        select(-rep)
    }) |>
    list_rbind(names_to = "model")
})[["elapsed"]]

# Runtime ----------------------------------------------------------------------
runtime <- bind_rows(
  tibble(
    method = names(fits),
    step = "fit",
    seconds = map_dbl(fits, "seconds"),
    cpus = c(CJ_CPUS, CJ_CPUS, 1L)
  ),
  tibble(
    method = names(score_runs),
    step = "cj_scores",
    seconds = map_dbl(score_runs, \(r) sum(map_dbl(r, "seconds"))),
    cpus = 1L
  ),
  tibble(
    method = "m1dm",
    step = "cj_visit_change",
    seconds = sum(map_dbl(visit_runs, "seconds")),
    cpus = 1L
  ),
  tibble(
    method = "COMPASS",
    step = "cj_compass",
    seconds = compass_run$seconds,
    cpus = 1L
  ),
  tibble(method = "m1dm_m1", step = "loo", seconds = loo_seconds, cpus = 1L),
  tibble(method = "m1dm_m1", step = "ppc", seconds = ppc_seconds, cpus = 1L)
)

# Write outputs ----------------------------------------------------------------
dir.create(partial_dir, recursive = TRUE)
write_tsv(scores, file.path(partial_dir, "scores.tsv"))
write_tsv(visit_change, file.path(partial_dir, "visit_change.tsv"))
write_tsv(compass, file.path(partial_dir, "compass.tsv"))
write_tsv(compass_run$fit$groups, file.path(partial_dir, "compass_groups.tsv"))
write_tsv(
  compass_run$fit$filters,
  file.path(partial_dir, "compass_filters.tsv")
)
write_tsv(diagnostics, file.path(partial_dir, "diagnostics.tsv"))
write_tsv(phi_draws, file.path(partial_dir, "phi_draws.tsv"))
write_tsv(loo_table, file.path(partial_dir, "loo.tsv"))
write_tsv(loo_compare, file.path(partial_dir, "loo_compare.tsv"))
write_tsv(loo_pointwise, file.path(partial_dir, "loo_pointwise.tsv"))
write_tsv(ppc, file.path(partial_dir, "ppc.tsv"))
write_tsv(sample_info, file.path(partial_dir, "samples.tsv"))
write_session_info(file.path(partial_dir, "session_info.txt"), PKG_DIR)
write_tsv(runtime, file.path(partial_dir, "runtime.tsv"))
file.rename(partial_dir, pop_dir)
