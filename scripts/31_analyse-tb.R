#!/usr/bin/env Rscript
# Purpose: benchmark cytoJointR (m1dm, m1, C0) against COMPASS and polyICSFlow on one prepared TB data set, in
#          three steps. fit (environment A): the three cytoJointR models and COMPASS on the same counts, donor
#          scores, model criticism (phi, LOO of m1dm vs m1, posterior predictive checks) and runtimes.
#          polyics (environment B): polyICSFlow on the same counts. summary (environment A): TB-status AUC with a
#          donor bootstrap and Spearman agreement between methods. CCS is experimental: reported, never a basis.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/31_analyse-tb.R <run_dir> <dataset> fit
#          R_LIBS=renv/lib-r44 /exports/para-lipg-hpc/mdmanurung/batchdx-cyto/renv/bin/Rscript scripts/31_analyse-tb.R <run_dir> <dataset> polyics
#          PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/31_analyse-tb.R <run_dir> <dataset> summary
#          <dataset> is tb_s, tb_u or tb_simplecompass. Each step refuses to overwrite its outputs.
# Inputs:
#   data/processed/<dataset>/   counts.tsv, meta.tsv, markers.tsv [22_prepare-tb-compass.R, 23_prepare-tb-simplecompass.R]
#   <PKG_DIR>/                  package source, loaded with pkgload [frozen worktree 49278c9]
#   scripts/helper01_true-scores.R (write_session_info), helper03_realdata.R (auc_boot), helper04_realdata-tb.R
# Outputs (in <run_dir>/<dataset>/ unless stated):
#   fit:     donor_scores.tsv   method, score, donor_id, value (cytoJointR P(RMS > eps), RMS / CPS / CCS means, COMPASS FS / PFS / Response)
#            model_status.tsv   per model: ok or failed with the message
#            diagnostics.tsv    per Stan model: R-hat, ESS, divergences, treedepth hits, E-BFMI
#            phi_summary.tsv    m1dm phi and inv_phi: mean, sd, 90% interval
#            phi_draws.tsv      m1dm posterior draws of phi and inv_phi
#            loo.tsv            per model: elpd_loo, SE, p_loo, Pareto-k counts (observation = one sample)
#            loo_compare.tsv    m1dm vs m1: elpd difference and its SE
#            ppc_sample.tsv     model, sample, statistic, observed, replicate mean and 90% interval, mid p-value
#            ppc_pattern.tsv    model, condition, pattern, statistic, observed, replicate mean and 90% interval, mid p-value
#            compass_status.tsv COMPASS run status, patterns kept
#            runtime.tsv        method, step, seconds, cpus
#            session_info.txt   R and package versions
#   polyics: polyicsflow/donor_scores.tsv, runtime.tsv, session_info.txt
#   summary: auc.tsv (method, score, AUC, 95% CI, n), agreement.tsv (Spearman between methods), donor_wide.tsv

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot(
  "usage: Rscript 31_analyse-tb.R <run_dir> <dataset> <fit|polyics|summary>" = length(
    args
  ) ==
    3,
  "step must be fit, polyics or summary" = args[3] %in%
    c("fit", "polyics", "summary")
)
run_dir <- args[1]
dataset <- args[2]
step <- args[3]

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
  if (step == "polyics") library(polyICSFlow)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign-run") # frozen worktree, commit 49278c9
SEED <- 20261103L # one seed for fits, PPC and bootstrap, as SEED_OFFSET in 05_define-benchmark.R
CJ_CPUS <- 4L # cj_fit runs 4 chains in parallel; the sbatch allocation is 4 CPUs
A0_C0 <- 0.5 # Dirichlet prior pseudo-count of C0, as 06_run-benchmark-replicate.R
DRAWS_C0 <- 4000L # draws of C0, as 4 chains x 1000 of the Stan fits
EPS_RMS <- 1e-4 # P(RMS > eps), as 05_define-benchmark.R; 1/N is 1.1e-5 to 1.1e-4 here, so eps is 1 to 9 cells
EPS_GRID <- c(1e-5, 1e-3, 1e-2) # sensitivity: P(RMS > eps) one decade below and two above; at 1e-4 every donor saturates in some data sets
EPSILON <- c(
  RMS = EPS_RMS,
  poly_absolute = 1e-4,
  poly_marginal = 1e-4,
  CPS = 0.01,
  CCS = 1e-4,
  breadth = 0.01,
  marker_magnitude = 1e-4
) # cj_scores() wants all of them; only RMS is analysed here, as in 06
PPC_DRAWS <- 1000L # posterior draws used for replicates; evenly thinned from the 4000
RARE_K <- 2L # a pattern is rare in a sample when its count is at most this (0.1 expected cells in 05)
AUC_BOOT <- 2000L # donor bootstrap replicates for the AUC interval
BACKGROUND <- "control" # condition label of the unstimulated sample
HEADLINE <- c(
  # method / score pairs compared for rank agreement
  "m1dm P_RMS_gt_eps",
  "m1dm RMS_mean",
  "m1dm CPS_mean",
  "m1 P_RMS_gt_eps",
  "C0 P_RMS_gt_eps",
  "COMPASS FS",
  "COMPASS PFS",
  "COMPASS Response",
  "polyICSFlow_clipped net_sum",
  "polyICSFlow_unclipped net_sum"
)

# Paths ------------------------------------------------------------------------
source("scripts/helper04_realdata-tb.R")
dat <- tb_read(dataset)
out_dir <- file.path(run_dir, dataset)
poly_dir <- file.path(out_dir, "polyicsflow")
if (step == "fit") {
  stopifnot("dataset output already exists" = !dir.exists(out_dir))
  partial_dir <- paste0(out_dir, ".partial") # outputs land here and move to out_dir when complete
  unlink(partial_dir, recursive = TRUE) # a stale one comes from a killed task
}
if (step == "polyics") {
  stopifnot(
    "fit step is not complete" = file.exists(file.path(
      out_dir,
      "session_info.txt"
    )),
    "polyicsflow output already exists" = !dir.exists(poly_dir)
  )
  partial_dir <- paste0(poly_dir, ".partial")
  unlink(partial_dir, recursive = TRUE)
}
if (step == "summary") {
  stopifnot(
    "polyics step is not complete" = file.exists(file.path(
      poly_dir,
      "session_info.txt"
    )),
    "summary output already exists" = !file.exists(file.path(
      out_dir,
      "auc.tsv"
    ))
  )
}

# Step fit ---------------------------------------------------------------------
# Skipped by the other steps.
if (step == "fit") {
  pkgload::load_all(PKG_DIR, quiet = TRUE)
  source("scripts/helper01_true-scores.R") # write_session_info()
  x <- tb_cj(dat$counts, dat$markers)
  keys <- apply(x$bits, 1, paste, collapse = "") # pattern_key of each column of x$counts
  allneg <- which(keys == strrep("0", length(dat$markers)))

  # Fit the three cytoJointR models ---------------------------------------------
  # A fit that errors is recorded as failed and the other methods go on.
  fit_model <- list(
    m1dm = \() cj_fit(x, model = "m1dm", level = "sample", seed = SEED),
    m1 = \() cj_fit(x, model = "m1", level = "sample", seed = SEED),
    C0 = \() {
      cytoJointR:::cj_fit_dirichlet(
        x,
        a0 = A0_C0,
        n_draws = DRAWS_C0,
        seed = SEED
      )
    }
  )
  fits <- fit_model |>
    imap(\(f, model) {
      tryCatch(
        {
          seconds <- system.time(fit <- f())[["elapsed"]]
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
    keep(\(f) inherits(f$fit, "cj_fit")) |>
    map(\(f) {
      dg <- f$fit$diagnostics
      tibble(
        max_rhat = dg$max_rhat,
        min_ess_bulk = dg$min_ess_bulk,
        min_ess_tail = dg$min_ess_tail,
        n_divergent = dg$n_divergent,
        n_max_treedepth = dg$n_max_treedepth,
        min_ebfmi = min(dg$ebfmi)
      )
    }) |>
    list_rbind(names_to = "model")

  # Donor scores of the cytoJointR models ---------------------------------------
  # Sample-level draws, as the benchmark (draws = "sample"). The contrast table carries donor_id.
  score_cj <- function(f, eps) {
    draws <- if (inherits(f, "cj_fit")) f$draws_sample else f
    cj_scores(draws, epsilon = replace(EPSILON, "RMS", eps)) |>
      filter(is.na(marker))
  }
  cj_donor <- fits_ok |>
    map(\(f) {
      seconds <- system.time(main <- score_cj(f$fit, EPS_RMS))[["elapsed"]]
      grid <- EPS_GRID |>
        set_names() |>
        map(\(eps) score_cj(f$fit, eps) |> filter(score == "RMS")) |>
        list_rbind(names_to = "eps")
      list(
        seconds = seconds,
        scores = bind_rows(
          main |>
            filter(score == "RMS") |>
            transmute(donor_id, score = "P_RMS_gt_eps", value = p_positive),
          main |>
            filter(score == "RMS") |>
            transmute(donor_id, score = "RMS_mean", value = estimate),
          main |>
            filter(score == "CPS") |>
            transmute(donor_id, score = "CPS_mean", value = estimate),
          main |>
            filter(score == "CCS") |>
            transmute(
              donor_id,
              score = "CCS_mean_experimental",
              value = estimate
            ),
          grid |>
            transmute(
              donor_id,
              score = paste0("P_RMS_gt_", eps),
              value = p_positive
            )
        )
      )
    })

  # COMPASS ---------------------------------------------------------------------
  # Default category filter on the same collapsed counts; FS, PFS and the degree-1 response probability.
  compass <- tryCatch(
    {
      seconds <- system.time(cc <- cj_compass(x, seed = SEED))[["elapsed"]]
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
        status = cc$groups |>
          select(status) |>
          mutate(patterns_filtered_out = sum(cc$filters$n_removed))
      )
    },
    error = \(e) {
      list(
        status = tibble(
          status = paste("failed:", conditionMessage(e)),
          patterns_filtered_out = NA
        )
      )
    }
  )
  donor_scores <- bind_rows(
    cj_donor |> map(\(m) m$scores) |> list_rbind(names_to = "method"),
    if (!is.null(compass$scores)) compass$scores |> mutate(method = "COMPASS")
  )

  # phi and LOO of m1dm vs m1 -----------------------------------------------------
  # phi = 1 / inv_phi is the Dirichlet-multinomial precision: large phi means the counts are multinomial.
  # log_lik is per sample, so LOO leaves one sample out; r_eff is from exp(log_lik) rescaled per sample
  # (the scale does not change r_eff, and the raw values underflow).
  phi_summary <- tibble()
  phi_draws <- tibble()
  loo_all <- list()
  if ("m1dm" %in% names(fits_ok)) {
    eng <- fits_ok$m1dm$fit$engine_fit
    phi <- as.vector(posterior::as_draws_matrix(eng$draws("phi")))
    phi_draws <- tibble(draw = seq_along(phi), phi = phi, inv_phi = 1 / phi)
    phi_summary <- fits_ok$m1dm$fit$diagnostics$summary |>
      filter(variable %in% c("phi[1]", "inv_phi[1]")) |>
      select(variable, mean, sd, q5, q95, rhat, ess_bulk)
  }
  loo_of <- \(fit) {
    arr <- fit$engine_fit$draws("log_lik")
    ll <- matrix(as.vector(posterior::as_draws_matrix(arr)), ncol = dim(arr)[3])
    chain_id <- rep(seq_len(dim(arr)[2]), each = dim(arr)[1])
    r_eff <- loo::relative_eff(
      exp(sweep(ll, 2, apply(ll, 2, max))),
      chain_id = chain_id
    )
    loo::loo(ll, r_eff = r_eff)
  }
  loo_all <- fits_ok[intersect(c("m1dm", "m1"), names(fits_ok))] |>
    map(\(f) loo_of(f$fit))
  loo_tab <- loo_all |>
    imap(\(l, model) {
      k <- l$diagnostics$pareto_k
      tibble(
        elpd_loo = l$estimates["elpd_loo", "Estimate"],
        se_elpd_loo = l$estimates["elpd_loo", "SE"],
        p_loo = l$estimates["p_loo", "Estimate"],
        n_samples = length(k),
        n_k_gt_0.7 = sum(k > 0.7),
        n_k_gt_1 = sum(k > 1),
        max_k = max(k)
      )
    }) |>
    list_rbind(names_to = "model")
  loo_compare_tab <- if (length(loo_all) == 2) {
    loo::loo_compare(loo_all) |>
      unclass() |>
      as_tibble(rownames = "model") |>
      select(model, elpd_diff, se_diff)
  } else {
    tibble()
  }

  # Posterior predictive checks ---------------------------------------------------
  # Replicates are drawn from the structured mean mu (not from p, which is conditioned on the observed y):
  # m1dm  y_rep ~ DirMult(N, phi * mu), m1  y_rep ~ Multinomial(N, mu). Statistics per sample: the
  # any-positive count and the number of rare positive patterns; per condition and pattern: the pooled
  # count and the number of samples with a zero count.
  set.seed(SEED)
  Ns <- x$samples$N
  observed <- x$counts
  replicate_counts <- function(model) {
    mu <- fits_ok[[model]]$fit$draws_structured$p
    idx <- unique(round(seq(1, dim(mu)[1], length.out = PPC_DRAWS)))
    phi <- if (model == "m1dm") phi_draws$phi[idx]
    rep_y <- array(0, c(length(idx), nrow(observed), ncol(observed)))
    for (s in seq_len(nrow(observed))) {
      m <- mu[idx, s, ]
      p <- if (is.null(phi)) {
        m
      } else {
        g <- matrix(pmax(rgamma(length(m), shape = phi * m), 1e-300), nrow(m))
        g / rowSums(g)
      }
      rep_y[, s, ] <- t(apply(p, 1, \(pr) rmultinom(1, Ns[s], pr)))
    }
    rep_y
  }
  summarise_rep <- \(obs, rep) {
    # obs: n; rep: draws x n -> observed, mean, 90% interval, mid p-value
    tibble(
      observed = obs,
      rep_mean = colMeans(rep),
      lower = apply(rep, 2, quantile, 0.05),
      upper = apply(rep, 2, quantile, 0.95),
      p_mid = colMeans(sweep(rep, 2, obs, ">")) +
        0.5 * colMeans(sweep(rep, 2, obs, "=="))
    )
  }
  ppc <- c("m1dm", "m1") |>
    intersect(names(fits_ok)) |>
    set_names() |>
    map(\(model) {
      rep_y <- replicate_counts(model)
      pos <- setdiff(seq_len(ncol(observed)), allneg)
      any_pos <- sweep(-rep_y[,, allneg], 2, Ns, "+")
      n_rare <- apply(rep_y[,, pos, drop = FALSE] <= RARE_K, c(1, 2), sum)
      by_sample <- bind_rows(
        summarise_rep(Ns - observed[, allneg], any_pos) |>
          mutate(statistic = "any_positive"),
        summarise_rep(rowSums(observed[, pos] <= RARE_K), n_rare) |>
          mutate(statistic = "n_rare_patterns")
      ) |>
        mutate(
          sample_id = rep(x$samples$sample_id, 2),
          condition = rep(x$samples$condition, 2),
          .before = 1
        )
      by_pattern <- x$samples$condition |>
        unique() |>
        set_names() |>
        map(\(cond) {
          rows <- which(x$samples$condition == cond)
          total <- apply(rep_y[, rows, , drop = FALSE], c(1, 3), sum)
          n_zero <- apply(rep_y[, rows, , drop = FALSE] == 0, c(1, 3), sum)
          bind_rows(
            summarise_rep(colSums(observed[rows, ]), total) |>
              mutate(statistic = "total"),
            summarise_rep(colSums(observed[rows, ] == 0), n_zero) |>
              mutate(statistic = "n_zero_samples")
          ) |>
            mutate(pattern_key = rep(keys, 2), .before = 1)
        }) |>
        list_rbind(names_to = "condition")
      list(sample = by_sample, pattern = by_pattern)
    })
  ppc_sample <- ppc |> map("sample") |> list_rbind(names_to = "model")
  ppc_pattern <- ppc |> map("pattern") |> list_rbind(names_to = "model")

  # Runtime ----------------------------------------------------------------------
  runtime <- bind_rows(
    fits_ok |>
      imap(\(f, model) {
        tibble(
          step = "fit",
          seconds = f$seconds,
          cpus = if_else(model == "C0", 1L, CJ_CPUS)
        )
      }) |>
      list_rbind(names_to = "method"),
    cj_donor |>
      map(\(m) tibble(step = "cj_scores", seconds = m$seconds, cpus = 1L)) |>
      list_rbind(names_to = "method"),
    tibble(
      method = "COMPASS",
      step = "cj_compass",
      seconds = compass$seconds %||% NA_real_,
      cpus = 1L
    )
  )

  # Write outputs ----------------------------------------------------------------
  dir.create(partial_dir, recursive = TRUE)
  write_tsv(donor_scores, file.path(partial_dir, "donor_scores.tsv"))
  write_tsv(model_status, file.path(partial_dir, "model_status.tsv"))
  write_tsv(diagnostics, file.path(partial_dir, "diagnostics.tsv"))
  write_tsv(phi_summary, file.path(partial_dir, "phi_summary.tsv"))
  write_tsv(phi_draws, file.path(partial_dir, "phi_draws.tsv"))
  write_tsv(loo_tab, file.path(partial_dir, "loo.tsv"))
  write_tsv(loo_compare_tab, file.path(partial_dir, "loo_compare.tsv"))
  write_tsv(ppc_sample, file.path(partial_dir, "ppc_sample.tsv"))
  write_tsv(ppc_pattern, file.path(partial_dir, "ppc_pattern.tsv"))
  write_tsv(compass$status, file.path(partial_dir, "compass_status.tsv"))
  write_tsv(runtime, file.path(partial_dir, "runtime.tsv"))
  write_session_info(file.path(partial_dir, "session_info.txt"), PKG_DIR) # last: marks a complete step
  file.rename(partial_dir, out_dir)
}

# Step polyics -----------------------------------------------------------------
# polyICSFlow on the cells expanded from the counts, with the control as background. 07 is tied to the
# benchmark layout and to four marker names, so the same three calls are made here with this panel.
if (step == "polyics") {
  expand_seconds <- system.time({
    cells <- dat$counts |> filter(n > 0) |> uncount(n)
    bits <- do.call(rbind, strsplit(cells$pattern_key, "")) == "1"
    colnames(bits) <- dat$markers
  })[["elapsed"]]
  assign_seconds <- system.time(
    comb <- assignMarkerCombinations(as_tibble(bits), mode = "simple") |>
      mutate(
        donor_id = cells$donor_id,
        condition = cells$condition,
        cell_type = "CD4"
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

  # Net per donor ------------------------------------------------------------------
  # Clipped: perc_pos_backgroundSubtracted. Unclipped: perc_pos minus the control's perc_pos. Summed over
  # the combinations polyICSFlow reports (it reports none for all-negative cells). A combination seen in only
  # one condition has 0 percent in the other, so the join is full; n_unpaired counts those rows, and a clipped
  # value missing for such a row is max(0, unclipped).
  net <- res |>
    filter(condition == BACKGROUND) |>
    select(donor_id, MarkerComb, perc_background = perc_pos) |>
    full_join(
      res |> filter(condition != BACKGROUND),
      join_by(donor_id, MarkerComb),
      relationship = "one-to-one"
    ) |>
    mutate(
      unpaired = is.na(perc_background) | is.na(perc_pos),
      unclipped = (coalesce(perc_pos, 0) - coalesce(perc_background, 0)) / 100,
      clipped = coalesce(
        perc_pos_backgroundSubtracted / 100,
        pmax(unclipped, 0)
      )
    )
  donor_scores <- net |>
    summarise(
      polyICSFlow_clipped = sum(clipped),
      polyICSFlow_unclipped = sum(unclipped),
      n_unpaired = sum(unpaired),
      .by = donor_id
    ) |>
    pivot_longer(
      starts_with("polyICSFlow_"),
      names_to = "method",
      values_to = "value"
    ) |>
    mutate(score = "net_sum")
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
    net |> select(donor_id, MarkerComb, clipped, unclipped, unpaired),
    file.path(partial_dir, "combination_net.tsv")
  )
  write_tsv(runtime, file.path(partial_dir, "runtime.tsv"))
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
  file.rename(partial_dir, poly_dir)
}

# Step summary -----------------------------------------------------------------
# AUC for the TB-status contrast with a donor bootstrap, and Spearman agreement of the headline scores.
if (step == "summary") {
  source("scripts/helper03_realdata.R") # auc_boot()
  scores <- bind_rows(
    read_tsv(
      file.path(out_dir, "donor_scores.tsv"),
      col_types = cols(donor_id = "c"),
      progress = FALSE
    ),
    read_tsv(
      file.path(poly_dir, "donor_scores.tsv"),
      col_types = cols(donor_id = "c"),
      progress = FALSE
    ) |>
      select(method, score, donor_id, value)
  ) |>
    left_join(
      dat$meta |> select(donor_id, case),
      join_by(donor_id),
      relationship = "many-to-one"
    )
  auc <- scores |>
    filter(!is.na(case)) |>
    reframe(
      auc_boot(value, case, B = AUC_BOOT, seed = SEED),
      .by = c(method, score)
    )
  wide <- scores |>
    mutate(column = paste(method, score)) |>
    select(donor_id, column, value) |>
    pivot_wider(names_from = column, values_from = value)
  agreement <- expand_grid(a = HEADLINE, b = HEADLINE) |>
    filter(a < b) |>
    mutate(
      spearman = map2_dbl(a, b, \(a, b) {
        cor(wide[[a]], wide[[b]], method = "spearman", use = "complete.obs")
      })
    )
  write_tsv(auc, file.path(out_dir, "auc.tsv"))
  write_tsv(agreement, file.path(out_dir, "agreement.tsv"))
  write_tsv(
    wide |> left_join(dat$meta, join_by(donor_id), relationship = "one-to-one"),
    file.path(out_dir, "donor_wide.tsv")
  )
}
