#!/usr/bin/env Rscript
# Purpose: diagnose the 22% placebo false-positive rate of cytoJointR m1dm on real data (RV144 placebo, TB QFT-)
#          where COMPASS gives 0%. Diagnosis only; the package is not touched.
#          fit   (env A): refit m1dm, score BOTH estimands (sample level p, structured level mu) at an epsilon
#                 grid, add the model-free stim - control difference, and check the pairwise log-linear mean
#                 (observed vs fitted mu for any-positive and degree >= 2, 3 patterns). Saves the draws.
#          stats (env A): AUC with a donor bootstrap, false-positive / true-positive rate at P > 0.95, the
#                 mean-misfit summary, the raw-difference null calibration and, for RV144, the infection OR.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/60_diagnose-estimand.R <fit|stats> <run_dir> <dataset>
#          <dataset> is rv144, tb_s, tb_u or tb_simplecompass (complete 256 -> 5 marker partition, no negatives).
#          Slurm: scripts/slurm_estimand.sh. Stan is already compiled (see slurm_vaccine.sh).
# Inputs:
#   data/processed/rv144/       cj_data.rds, donor_meta.tsv [20_prepare-rv144.R]
#   data/processed/<tb>/        counts.tsv, meta.tsv, markers.tsv [22_prepare-tb-compass.R, 23_prepare-tb-simplecompass.R]
#   <PKG_DIR>/                  package source, loaded with pkgload [frozen worktree, commit 49278c9]
#   <OLD_RUN>/                  COMPASS scores of the main benchmark, read only [30_analyse-vaccine.R, 31_analyse-tb.R]
#   scripts/helper01_true-scores.R, helper03_realdata.R, helper04_realdata-tb.R
# Outputs (in <run_dir>/<dataset>/):
#   fit:   scores_levels.tsv    level (sample, structured, raw), donor_id, score, value, n_failed_draws
#          misfit_samples.tsv   sample, statistic, observed fraction, fitted mu fraction (mean, 5%, 95%)
#          fit_draws.rds        thinned p (sample level) and mu (structured) arrays, parameter draws, design
#          diagnostics.tsv, runtime.tsv, session_info.txt
#   stats: performance.tsv      level, score, AUC and donor-bootstrap CI, FPR and TPR at the rule
#          level_summary.tsv    median and quartiles of each score by group
#          misfit_summary.tsv   fitted / observed per group, condition and statistic; net stim - control
#          misfit_donor.tsv     per-donor net (stim - control) observed vs fitted
#          raw_null.tsv         calibration of the raw difference z in the negative group
#          infection.tsv        (rv144) OR per SD among vaccinees, covariates as 30_analyse-vaccine.R M1

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign-run") # frozen package copy, commit 49278c9
OLD_RUN <- "runs/2026-11-03_realdata" # COMPASS scores, read only
FIT_ARGS <- Sys.getenv("CJ_FIT_ARGS", "model=m1dm,level=sample") # cj_fit arguments of the fit stage, key=value pairs
SEED <- 20261103L # seed of the fit and the bootstrap, as 30_analyse-vaccine.R
EPS_GRID <- c(1e-4, 5e-4, 1e-3) # 1e-4 is the pre-registered rule; the others only show sensitivity, no re-tuning
RULE_P <- 0.95 # P(RMS > eps) above this calls a donor positive
Z_RAW <- qnorm(0.95) # raw difference: one-sided nominal 5% test
THIN <- 1000L # draws kept in fit_draws.rds (scoring uses all 4000)
BOOT <- 2000L # donor bootstrap replicates for AUC
INFECTION_COV <- c("gender", "age", "risk.medium", "risk.high") # covariates of M1 in 30_analyse-vaccine.R
STATS_DEF <- list(any_positive = 1L, degree2plus = 2L, degree3plus = 3L) # minimum pattern degree of each statistic

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 60_diagnose-estimand.R <fit|stats> <run_dir> <rv144|tb_s|tb_u|tb_simplecompass>" = length(
    args
  ) ==
    3,
  args[1] %in% c("fit", "stats"),
  args[3] %in% c("rv144", "tb_s", "tb_u", "tb_simplecompass")
)
stage <- args[1]
ds <- args[3]
out_dir <- file.path(args[2], ds)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- \(f) file.path(out_dir, f)
stopifnot(
  "stage output already exists; use a fresh run dir" = !file.exists(out(
    if (stage == "fit") "diagnostics.tsv" else "performance.tsv"
  ))
)
source("scripts/helper03_realdata.R")

# Stage fit --------------------------------------------------------------------
run_fit <- function() {
  pkgload::load_all(PKG_DIR, quiet = TRUE)
  source("scripts/helper01_true-scores.R") # write_session_info()
  source("scripts/helper04_realdata-tb.R") # tb_read(), tb_cj()
  x <- if (ds == "rv144") {
    readRDS("data/processed/rv144/cj_data.rds")
  } else {
    dat <- tb_read(ds)
    tb_cj(dat$counts, dat$markers)
  }
  stopifnot(inherits(x, "cj_data"))

  # Fit once; both estimands come from the same posterior.
  fit_args <- strsplit(strsplit(FIT_ARGS, ",")[[1]], "=") |>
    (\(kv) setNames(map(kv, \(p) type.convert(p[2], as.is = TRUE)), map_chr(kv, 1)))()
  sec <- system.time(
    fit <- do.call(cj_fit, c(list(x, seed = SEED), fit_args))
  )[["elapsed"]]
  draws <- list(sample = fit$draws_sample, structured = fit$draws_structured)
  ct <- x$contrasts[, c(
    "contrast_id",
    "donor_id",
    "stim_sample_id",
    "control_sample_id"
  )]

  # Score each estimand: RMS and CPS means and P(RMS > eps) over the epsilon grid.
  eps_lab <- format(EPS_GRID, scientific = TRUE)
  score_level <- \(d, level) {
    eng <- cytoJointR:::score_engine(d)
    imap(eng, \(e, id) {
      ok <- e$status == "ok"
      rms <- e$G[ok, "RMS"]
      tibble(
        level,
        contrast_id = id,
        score = c("RMS_mean", "CPS_mean", paste0("P_RMS_gt_", eps_lab)),
        value = c(
          mean(rms),
          mean(e$G[ok, "CPS"]),
          map_dbl(EPS_GRID, \(e0) mean(rms > e0))
        ),
        n_failed_draws = sum(!ok)
      )
    }) |>
      list_rbind()
  }
  sec_score <- system.time(scored <- imap(draws, score_level) |> list_rbind())[[
    "elapsed"
  ]]

  # Model-free difference: stimulated minus control positive fraction, and its binomial z.
  pos <- x$bits |> rowSums() >= 1
  q <- \(id) {
    sum(x$counts[id, pos]) / x$samples$N[match(id, x$samples$sample_id)]
  }
  n_of <- \(id) x$samples$N[match(id, x$samples$sample_id)]
  raw <- ct |>
    mutate(
      q_s = map_dbl(stim_sample_id, q),
      q_u = map_dbl(control_sample_id, q),
      se = sqrt(
        q_s *
          (1 - q_s) /
          n_of(stim_sample_id) +
          q_u * (1 - q_u) / n_of(control_sample_id)
      )
    ) |>
    transmute(
      contrast_id,
      level = "raw",
      RMS_mean = q_s - q_u,
      z = (q_s - q_u) / se,
      n_failed_draws = 0L
    ) |>
    pivot_longer(c(RMS_mean, z), names_to = "score", values_to = "value")
  write_tsv(
    bind_rows(scored, raw) |>
      left_join(
        ct |> select(contrast_id, donor_id),
        join_by(contrast_id),
        relationship = "many-to-one"
      ) |>
      select(level, donor_id, score, value, n_failed_draws),
    out("scores_levels.tsv")
  )

  # Mean-misfit check: fitted structured mu vs observed fractions, per sample.
  mu <- fit$draws_structured$p
  deg <- rowSums(x$bits)
  obs_n <- x$samples$N[match(rownames(x$counts), x$samples$sample_id)]
  misfit <- imap(STATS_DEF, \(min_deg, nm) {
    sel <- deg >= min_deg
    m <- rowSums(mu[,, sel, drop = FALSE], dims = 2) # draw x sample
    tibble(
      sample_id = rownames(x$counts),
      statistic = nm,
      observed = rowSums(x$counts[, sel, drop = FALSE]) / obs_n,
      fitted = colMeans(m),
      fitted_q05 = apply(m, 2, quantile, 0.05),
      fitted_q95 = apply(m, 2, quantile, 0.95)
    )
  }) |>
    list_rbind()
  write_tsv(
    misfit |>
      left_join(
        x$samples |> select(sample_id, donor_id, condition, N),
        join_by(sample_id),
        relationship = "many-to-one"
      ),
    out("misfit_samples.tsv")
  )

  # Draws: thinned p and mu plus the parameters, so the estimands can be re-scored without a refit.
  idx <- round(seq(1, dim(mu)[1], length.out = min(THIN, dim(mu)[1])))
  vars <- intersect(
    c(
      "alpha",
      "gamma",
      "theta0",
      "delta",
      "tau",
      "rho",
      "tau_u",
      "inv_phi",
      "phi"
    ),
    fit$engine_fit$metadata()$stan_variables
  )
  saveRDS(
    list(
      p_sample = fit$draws_sample$p[idx, , , drop = FALSE],
      mu_structured = mu[idx, , , drop = FALSE],
      draw_index = idx,
      params = posterior::as_draws_df(fit$engine_fit$draws(vars)),
      samples = x$samples,
      contrasts = x$contrasts,
      bits = x$bits,
      seed = SEED
    ),
    out("fit_draws.rds")
  )
  write_tsv(
    tibble(
      dataset = ds,
      n_samples = nrow(x$counts),
      N_median = median(x$samples$N),
      max_rhat = fit$diagnostics$max_rhat,
      min_ess_bulk = fit$diagnostics$min_ess_bulk,
      n_divergent = fit$diagnostics$n_divergent,
      phi_mean = mean(as.vector(posterior::as_draws_matrix(fit$engine_fit$draws(
        "phi"
      ))))
    ),
    out("diagnostics.tsv")
  )
  write_tsv(
    tibble(
      step = c("fit", "score_both_levels"),
      seconds = c(sec, sec_score),
      cpus = c(4L, 1L)
    ),
    out("runtime.tsv")
  )
  write_session_info(out("session_info.txt"), PKG_DIR)
}

# Stage stats ------------------------------------------------------------------
cp_ci <- \(k, n) unname(binom.test(k, n)$conf.int)

run_stats <- function() {
  # Truth: RV144 vaccine vs placebo; TB QFT+ vs QFT-. tb_simplecompass has no negatives: misfit only.
  meta <- if (ds == "rv144") {
    read_tsv("data/processed/rv144/donor_meta.tsv", show_col_types = FALSE) |>
      mutate(
        across(any_of(c("gender", "age")), as.factor),
        truth = arm == "vaccine"
      )
  } else {
    read_tsv(
      file.path("data/processed", ds, "meta.tsv"),
      show_col_types = FALSE
    ) |>
      mutate(truth = if (ds == "tb_simplecompass") NA else group == "positive")
  }
  group_lab <- c(
    `TRUE` = if (ds == "rv144") "vaccine" else "QFT+",
    `FALSE` = if (ds == "rv144") "placebo" else "QFT-"
  )
  if (ds != "tb_simplecompass") {
    meta <- meta |> mutate(group = unname(group_lab[as.character(truth)]))
  }
  sc <- read_tsv(out("scores_levels.tsv"), show_col_types = FALSE)

  # COMPASS from the main benchmark, read only.
  compass <- if (ds == "rv144") {
    read_tsv(
      file.path(OLD_RUN, "vaccine/rv144/compass_scores.tsv"),
      show_col_types = FALSE
    )
  } else if (ds != "tb_simplecompass") {
    read_tsv(
      file.path(OLD_RUN, "tb", ds, "donor_scores.tsv"),
      show_col_types = FALSE
    ) |>
      filter(method == "COMPASS")
  }
  if (!is.null(compass)) {
    sc <- bind_rows(
      sc,
      compass |>
        filter(score %in% c("FS", "PFS", "Response")) |>
        transmute(level = "COMPASS", donor_id, score, value)
    )
  }
  wide <- sc |>
    mutate(name = paste(level, score, sep = "|")) |>
    select(donor_id, name, value) |>
    pivot_wider(names_from = name, values_from = value)
  tab <- meta |> left_join(wide, join_by(donor_id), relationship = "one-to-one")
  stopifnot(
    "scores not matched one-to-one to donors" = nrow(tab) == nrow(meta),
    !anyNA(tab$donor_id)
  )

  # Misfit summaries (all datasets) ------------------------------------------
  ms <- read_tsv(out("misfit_samples.tsv"), show_col_types = FALSE) |>
    left_join(
      meta |> select(donor_id, group),
      join_by(donor_id),
      relationship = "many-to-one"
    ) |>
    mutate(
      group = coalesce(group, "all"),
      cond = if_else(condition == "control", "control", "stimulated")
    )
  write_tsv(
    ms |>
      summarise(
        n = n(),
        observed = mean(observed),
        fitted = mean(fitted),
        fitted_over_observed = mean(fitted) / mean(observed),
        frac_obs_below_q05 = mean(observed < fitted_q05),
        frac_obs_above_q95 = mean(observed > fitted_q95),
        .by = c(group, cond, statistic)
      ) |>
      arrange(statistic, cond, group) |>
      mutate(dataset = ds, .before = 1),
    out("misfit_summary.tsv")
  )
  net <- ms |>
    select(donor_id, group, cond, statistic, observed, fitted) |>
    pivot_wider(names_from = cond, values_from = c(observed, fitted)) |>
    transmute(
      donor_id,
      group,
      statistic,
      obs_net = observed_stimulated - observed_control,
      fit_net = fitted_stimulated - fitted_control
    )
  write_tsv(net |> mutate(dataset = ds, .before = 1), out("misfit_donor.tsv"))

  if (anyNA(meta$truth)) {
    return(invisible())
  } # no negatives: nothing more to test
  truth <- tab$truth

  # Performance ---------------------------------------------------------------
  # AUC for every score; FPR and TPR at the rule for P(RMS > eps), COMPASS Response and the raw z.
  scores <- grep("\\|", names(tab), value = TRUE)
  rules <- c(
    setNames(
      rep(RULE_P, sum(grepl("\\|P_RMS_gt_", scores))),
      grep("\\|P_RMS_gt_", scores, value = TRUE)
    ),
    `COMPASS|Response` = 0.5,
    `raw|z` = Z_RAW
  )
  perf <- map(scores, \(nm) {
    s <- tab[[nm]]
    called <- if (nm %in% names(rules)) s > rules[[nm]]
    r <- if (!is.null(called)) {
      k0 <- sum(called[!truth], na.rm = TRUE)
      n0 <- sum(!truth & !is.na(called))
      k1 <- sum(called[truth], na.rm = TRUE)
      n1 <- sum(truth & !is.na(called))
      tibble(
        rule_cut = rules[[nm]],
        n_neg = n0,
        n_false_pos = k0,
        fpr = k0 / n0,
        fpr_lower = cp_ci(k0, n0)[1],
        fpr_upper = cp_ci(k0, n0)[2],
        n_pos = n1,
        n_called = k1,
        tpr = k1 / n1
      )
    }
    bind_cols(auc_boot(s, truth, BOOT, SEED), tibble(score_name = nm), r)
  }) |>
    list_rbind() |>
    separate_wider_delim(score_name, "|", names = c("level", "score")) |>
    mutate(dataset = ds, .before = 1)
  write_tsv(perf, out("performance.tsv"))
  write_tsv(
    tab |>
      pivot_longer(
        all_of(scores),
        names_to = "score_name",
        values_to = "value"
      ) |>
      summarise(
        n = sum(!is.na(value)),
        median = median(value, na.rm = TRUE),
        q25 = quantile(value, 0.25, na.rm = TRUE),
        q75 = quantile(value, 0.75, na.rm = TRUE),
        .by = c(group, score_name)
      ) |>
      mutate(dataset = ds, .before = 1),
    out("level_summary.tsv")
  )
  # Raw null: among negatives the raw z should be N(0, 1) if stim and control differ only by multinomial noise.
  z0 <- tab[["raw|z"]][!truth]
  write_tsv(
    tibble(
      dataset = ds,
      n_neg = length(z0),
      z_mean = mean(z0),
      z_sd = sd(z0),
      frac_z_gt_nominal = mean(z0 > Z_RAW),
      frac_abs_z_gt_1.96 = mean(abs(z0) > 1.96)
    ),
    out("raw_null.tsv")
  )

  # Infection association among vaccinees (RV144) ------------------------------
  if (ds == "rv144") {
    v <- tab |> filter(truth, !is.na(infected))
    cols <- c(
      sample_RMS = "sample|RMS_mean",
      structured_RMS = "structured|RMS_mean",
      raw_RMS = "raw|RMS_mean",
      sample_CPS = "sample|CPS_mean",
      structured_CPS = "structured|CPS_mean",
      COMPASS_FS = "COMPASS|FS",
      COMPASS_PFS = "COMPASS|PFS"
    )
    write_tsv(
      imap(cols, \(cn, lab) {
        d <- v |>
          select(infected, all_of(c(cn, INFECTION_COV))) |>
          drop_na() |>
          mutate(z = as.numeric(scale(.data[[cn]])))
        f <- glm(reformulate(c("z", INFECTION_COV), "infected"), binomial, d)
        b <- unname(coef(f)["z"])
        se <- sqrt(vcov(f)["z", "z"])
        tibble(
          label = lab,
          n = nrow(d),
          events = sum(d$infected),
          or_per_sd = exp(b),
          lower = exp(b - 1.96 * se),
          upper = exp(b + 1.96 * se),
          p = 2 * pnorm(-abs(b / se))
        )
      }) |>
        list_rbind() |>
        mutate(dataset = ds, .before = 1),
      out("infection.tsv")
    )
  }
}

# Run --------------------------------------------------------------------------
switch(stage, fit = run_fit(), stats = run_stats())
