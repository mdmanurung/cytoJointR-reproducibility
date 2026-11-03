#!/usr/bin/env Rscript
# Purpose: vaccine trials (RV144, HVTN078) in the real-data benchmark, one stage per call.
#          fit    (env A): cytoJointR m1dm (default), m1 and C0 on the 5-marker cj_data; donor scores RMS, CPS, CCS
#                 (CCS experimental), posterior draws for the infection analysis, phi, LOO, posterior predictive check.
#          compass (env A): COMPASS on the SAME reduced tables via cj_compass() with its default filter; RV144 also
#                 reads the paper's native 6-marker COMPASS result (reference only).
#          polyics (env B): polyICSFlow clipped and unclipped net on the same counts.
#          stats  (env A): detection AUC with a donor bootstrap, false-positive rates at each method's rule,
#                 RV144 infection association (incl. posterior-draw propagation), HVTN078 regimen contrast,
#                 method agreement, runtime table.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/30_analyse-vaccine.R <fit|compass|stats> <run_dir> <rv144|hvtn078>
#          R_LIBS=renv/lib-r44 <R_B>/Rscript scripts/30_analyse-vaccine.R polyics <run_dir> <dataset>
#          Slurm: scripts/slurm_vaccine.sh. Compile Stan once first (see that file).
# Inputs:
#   data/processed/<dataset>/   cj_data.rds, counts_long.tsv, donor_meta.tsv [20_prepare-rv144.R, 21_prepare-hvtn078.R]
#   <PKG_DIR>/                  package source, loaded with pkgload [frozen worktree, commit 49278c9]
#   scripts/helper03_realdata.R, helper01_true-scores.R   auc_rank(), auc_boot(), write_session_info()
#   data/raw/compass_data_package/.../COMPASSData.rda   native RV144 COMPASS fit (compass stage, reference only)
# Outputs (in <run_dir>/<dataset>/):
#   fit:     model_status.tsv, runtime_fit.tsv, diagnostics.tsv, cj_scores.tsv (model, donor_id, score, estimate, lower,
#            upper, p_positive, status; RMS also at eps 1e-3), draws_prop.tsv (thinned per-draw RMS, CPS),
#            phi_draws.tsv, phi_summary.tsv, loo.tsv, loo_compare.tsv, ppc_samples.tsv, ppc_summary.tsv,
#            ppc_quantiles.tsv, epsilon_check.tsv, session_info_fit.txt
#   compass: compass_scores.tsv, compass_filters.tsv, compass_groups.tsv, compass_native.tsv (rv144), runtime_compass.tsv
#   polyics: polyics_scores.tsv, runtime_polyics.tsv, session_info_polyics.txt
#   stats:   donor_table.tsv, auc.tsv, auc_diff.tsv, roc_curves.tsv, fpr.tsv, infection.tsv, infection_propagated.tsv,
#            regimen.tsv, regimen_means.tsv, agreement.tsv, runtime.tsv

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign-run") # frozen package copy, commit 49278c9
PROCESSED <- "data/processed"
SEED <- 20261103L # seed of every fit, draw thinning and bootstrap
CJ_CPUS <- 4L
DRAWS_C0 <- 4000L # as 06_run-benchmark-replicate.R
A0_C0 <- 0.5
EPS_RESPONSE <- 1e-4 # as 05_define-benchmark.R; 1e-4 is above 1 / N for every sample here (checked in epsilon_check.tsv)
EPS_SENS <- 1e-3 # sensitivity: P(RMS > 1e-3), 10x the main epsilon
EPSILON <- c(RMS = EPS_RESPONSE, poly_absolute = 1e-4, poly_marginal = 1e-4, CPS = 0.01, CCS = 1e-4,
             breadth = 0.01, marker_magnitude = 1e-4) # as 06
PROP_DRAWS <- 200L # draws kept per contrast for the regression refits (200 x logistic fit is cheap, Rubin needs few)
PPC_DRAWS <- 400L # posterior predictive replicates per sample
PPC_RARE <- 3L # rare-pattern statistics: the 3 polyfunctional (degree >= 2) patterns with the largest stimulated total
LOO_PARETO_K <- 0.7
BOOT <- 2000L # donor bootstrap replicates for AUC CIs
RULE_CJ <- 0.95 # cytoJointR: P(RMS > eps) above this
RULE_COMPASS <- 0.5 # COMPASS: Response (degree 1) above this
RULE_POLY <- 0 # polyICSFlow clipped net above this
HEADLINE <- c("COMPASS_FS", "COMPASS_PFS", "m1dm_RMS", "m1dm_CPS") # Holm family for the infection analysis, M1 adjustment
COVARIATES <- list( # infection model covariates; columns of donor_meta.tsv
  M0 = character(),
  M1 = c("gender", "age", "risk.medium", "risk.high"), # baseline covariates present in meta
  M2 = c("gender", "age", "risk.medium", "risk.high", "IgAprim", "V2prim") # + the antibody correlates in meta
)
RDA <- "data/raw/compass_data_package/pkg/COMPASSDataPackage/data/COMPASSData.rda"
VACCINE_LEVEL <- list(rv144 = "vaccine", hvtn078 = "vaccine")

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
stopifnot("usage: Rscript 30_analyse-vaccine.R <fit|compass|polyics|stats> <run_dir> <rv144|hvtn078>" =
            length(args) == 3, args[1] %in% c("fit", "compass", "polyics", "stats"), args[3] %in% c("rv144", "hvtn078"))
stage <- args[1]
ds <- args[3]
in_dir <- file.path(PROCESSED, ds)
out_dir <- file.path(args[2], ds)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out <- \(f) file.path(out_dir, f)
stopifnot("stage output already exists; use a fresh run dir" = !file.exists(out(switch(stage,
  fit = "diagnostics.tsv", compass = "compass_scores.tsv", polyics = "polyics_scores.tsv", stats = "auc.tsv"))))
if (stage %in% c("fit", "compass", "stats")) {
  pkgload::load_all(PKG_DIR, quiet = TRUE)
}
source("scripts/helper03_realdata.R")

# Stage fit --------------------------------------------------------------------
run_fit <- function() {
  source("scripts/helper01_true-scores.R") # write_session_info()
  x <- readRDS(file.path(in_dir, "cj_data.rds"))
  stopifnot(inherits(x, "cj_data"))
  cj_formals <- names(formals(cj_fit))
  stopifnot("frozen package lacks model argument" = "model" %in% cj_formals)
  fit_model <- \(m) switch(m,
    m1dm = cj_fit(x, model = "m1dm", level = "sample", seed = SEED),
    m1 = cj_fit(x, model = "m1", level = "sample", seed = SEED),
    C0 = cytoJointR:::cj_fit_dirichlet(x, a0 = A0_C0, n_draws = DRAWS_C0, seed = SEED))
  fits <- c("m1dm", "m1", "C0") |>
    set_names() |>
    map(\(m) tryCatch({
      sec <- system.time(f <- fit_model(m))[["elapsed"]]
      list(status = "ok", message = NA_character_, fit = f, seconds = sec)
    }, error = \(e) list(status = "failed", message = conditionMessage(e))))
  write_tsv(imap(fits, \(f, m) tibble(model = m, status = f$status, message = f$message)) |> list_rbind(), out("model_status.tsv"))
  fits <- keep(fits, \(f) f$status == "ok")

  # Scores and thinned draws, one pass of the engine per model.
  score_one <- \(f, model) {
    d <- if (inherits(f$fit, "cj_fit")) f$fit$draws_sample else f$fit
    nd <- dim(d$p)[1]
    idx <- round(seq(1, nd, length.out = min(PROP_DRAWS, nd)))
    sec <- system.time(eng <- cytoJointR:::score_engine(d))[["elapsed"]]
    rows <- imap(eng, \(e, id) {
      if (!e$available) return(NULL)
      bad <- e$status != "ok"
      s <- cytoJointR:::summarise_columns(e$cols, e$G, bad, EPSILON[cytoJointR:::split_score_col(e$cols)$score], 0.9, 0.01) |>
        filter(is.na(marker), score %in% c("RMS", "CPS", "CCS")) |>
        select(score, estimate, lower, upper, p_positive, status)
      rms <- e$G[!bad, "RMS"]
      sens <- tibble(score = "RMS_eps1e-3", estimate = NA_real_, lower = NA_real_, upper = NA_real_,
                     p_positive = if (mean(bad) > 0.01) NA_real_ else mean(rms > EPS_SENS), status = s$status[1])
      g <- e$G[idx, c("RMS", "CPS")]
      g[bad[idx], ] <- NA
      list(scores = bind_rows(s, sens) |> mutate(contrast_id = id),
           draws = tibble(contrast_id = id, draw = seq_along(idx), RMS = g[, "RMS"], CPS = g[, "CPS"]))
    })
    ct <- d$contrasts[, c("contrast_id", "donor_id")]
    list(
      scores = map(rows, "scores") |> list_rbind() |> left_join(ct, join_by(contrast_id), relationship = "many-to-one") |>
        mutate(model = model) |> select(model, donor_id, score, estimate, lower, upper, p_positive, status),
      draws = if (model != "C0") {
        map(rows, "draws") |> list_rbind() |> left_join(ct, join_by(contrast_id), relationship = "many-to-one") |>
          pivot_longer(c(RMS, CPS), names_to = "score", values_to = "value") |> mutate(model = model) |>
          select(model, donor_id, draw, score, value)
      },
      seconds = sec
    )
  }
  scored <- imap(fits, \(f, m) score_one(f, m))
  write_tsv(map(scored, "scores") |> list_rbind(), out("cj_scores.tsv"))
  write_tsv(map(scored, "draws") |> compact() |> list_rbind(), out("draws_prop.tsv"))
  write_tsv(tibble(model = names(fits), step = "fit", seconds = map_dbl(fits, "seconds"), cpus = if_else(names(fits) == "C0", 1L, CJ_CPUS)) |>
              bind_rows(tibble(model = names(scored), step = "scores", seconds = map_dbl(scored, "seconds"), cpus = 1L)),
            out("runtime_fit.tsv"))
  write_tsv(tibble(
    N = x$samples$N,
    dataset = ds
  ) |> summarise(dataset = first(dataset), eps = EPS_RESPONSE, N_min = min(N), N_median = median(N), N_max = max(N),
                 inv_N_at_median = 1 / median(N), inv_N_at_min = 1 / min(N), eps_over_inv_N_min = EPS_RESPONSE * min(N),
                 expected_cells_at_eps_median_N = EPS_RESPONSE * median(N)), out("epsilon_check.tsv"))

  # Diagnostics, phi ---------------------------------------------------------
  stan_fits <- keep(fits, \(f) inherits(f$fit, "cj_fit"))
  write_tsv(imap(stan_fits, \(f, m) {
    dg <- f$fit$diagnostics
    tibble(model = m, max_rhat = dg$max_rhat, min_ess_bulk = dg$min_ess_bulk, min_ess_tail = dg$min_ess_tail,
           n_divergent = dg$n_divergent, n_max_treedepth = dg$n_max_treedepth, min_ebfmi = min(dg$ebfmi),
           fit_seconds = f$seconds)
  }) |> list_rbind(), out("diagnostics.tsv"))
  if ("m1dm" %in% names(stan_fits)) {
    eng <- stan_fits$m1dm$fit$engine_fit
    phi <- as.vector(posterior::as_draws_matrix(eng$draws("phi")))
    Nmed <- median(x$samples$N)
    write_tsv(tibble(draw = seq_along(phi), phi = phi)[round(seq(1, length(phi), length.out = 1000)), ], out("phi_draws.tsv"))
    write_tsv(tibble(
      dataset = ds, phi_mean = mean(phi), phi_median = median(phi), phi_q05 = quantile(phi, 0.05), phi_q95 = quantile(phi, 0.95),
      p_phi_gt_1e6 = mean(phi > 1e6), # phi this large is indistinguishable from multinomial at these N
      N_median = Nmed, # variance inflation over multinomial for a sample: (N + phi) / (1 + phi)
      var_inflation_median = median((Nmed + phi) / (1 + phi)), var_inflation_q05 = quantile((Nmed + phi) / (1 + phi), 0.05),
      var_inflation_q95 = quantile((Nmed + phi) / (1 + phi), 0.95)
    ), out("phi_summary.tsv"))
  }

  # LOO ----------------------------------------------------------------------
  # Per sample (leave one sample out) and per donor (a donor's samples summed, leave the donor out).
  loo_run <- \(f) {
    eng <- f$fit$engine_fit
    ll <- posterior::as_draws_matrix(eng$draws("log_lik"))
    chain_id <- rep(seq_len(eng$num_chains()), each = nrow(ll) / eng$num_chains())
    donor_of <- x$samples$donor_id[match(rownames(x$counts), x$samples$sample_id)]
    ll_d <- t(rowsum(t(unclass(ll)), donor_of))
    mk <- \(m) loo::loo(unclass(m), r_eff = loo::relative_eff(exp(sweep(unclass(m), 2, apply(m, 2, max))), chain_id = chain_id))
    list(sample = mk(ll), donor = mk(ll_d))
  }
  loos <- map(stan_fits, loo_run)
  if (length(loos)) {
    write_tsv(imap(loos, \(l, m) imap(l, \(o, lev) tibble(
      model = m, level = lev, n_obs = nrow(o$pointwise), elpd_loo = o$estimates["elpd_loo", "Estimate"],
      se_elpd = o$estimates["elpd_loo", "SE"], p_loo = o$estimates["p_loo", "Estimate"],
      n_pareto_k_gt_0.7 = sum(o$diagnostics$pareto_k > LOO_PARETO_K))) |> list_rbind()) |> list_rbind(), out("loo.tsv"))
  }
  if (all(c("m1dm", "m1") %in% names(loos))) {
    write_tsv(map(c("sample", "donor"), \(lev) {
      cmp <- loo::loo_compare(list(m1dm = loos$m1dm[[lev]], m1 = loos$m1[[lev]]))
      tibble(level = lev, model = rownames(cmp), elpd_diff = cmp[, "elpd_diff"], se_diff = cmp[, "se_diff"])
    }) |> list_rbind(), out("loo_compare.tsv"))
  }

  # Posterior predictive check -----------------------------------------------
  # Replicate counts from the structured mu (not conditional on the observed counts): m1dm y ~ DM(phi mu),
  # m1 y ~ Multinomial(mu). Statistics: any-positive count, polyfunctional (degree >= 3) count, and PPC_RARE rare patterns.
  bits <- x$bits
  deg <- rowSums(bits)
  stim_rows <- x$samples$sample_id[x$samples$condition != "control"]
  rare <- order(colSums(x$counts[stim_rows, , drop = FALSE]) * (deg >= 2), decreasing = TRUE)[seq_len(PPC_RARE)]
  stats_def <- c(list(any_positive = which(deg >= 1), degree3plus = which(deg >= 3)),
                 set_names(as.list(rare), paste0("pattern_", colnames(x$counts)[rare])))
  sm <- x$samples[match(rownames(x$counts), x$samples$sample_id), ]
  ppc_one <- \(f, model) {
    mu_all <- f$fit$draws_structured$p
    idx <- round(seq(1, dim(mu_all)[1], length.out = PPC_DRAWS))
    mu <- mu_all[idx, , , drop = FALSE]
    phi <- if (model == "m1dm") as.vector(posterior::as_draws_matrix(f$fit$engine_fit$draws("phi")))[idx]
    set.seed(SEED)
    rep_stat <- array(NA_real_, c(length(idx), nrow(x$counts), length(stats_def)))
    for (t in seq_along(idx)) for (s in seq_len(nrow(x$counts))) {
      m <- mu[t, s, ]
      p <- if (!is.null(phi)) { g <- stats::rgamma(length(m), shape = phi[t] * m); if (sum(g) > 0) g / sum(g) else m } else m
      y <- stats::rmultinom(1, sm$N[s], p)[, 1]
      for (k in seq_along(stats_def)) rep_stat[t, s, k] <- sum(y[stats_def[[k]]])
    }
    obs <- map(stats_def, \(i) rowSums(x$counts[, i, drop = FALSE]))
    samples <- imap(stats_def, \(i, nm) {
      k <- match(nm, names(stats_def)); r <- rep_stat[, , k]; o <- obs[[nm]]
      tibble(model = model, stat = nm, sample_id = rownames(x$counts), condition = sm$condition, observed = o,
             rep_mean = colMeans(r), rep_q05 = apply(r, 2, quantile, 0.05), rep_q95 = apply(r, 2, quantile, 0.95),
             ppp = map_dbl(seq_along(o), \(s) mean(r[, s] < o[s]) + 0.5 * mean(r[, s] == o[s])))
    }) |> list_rbind()
    probs <- c(0.1, 0.25, 0.5, 0.75, 0.9, 0.99)
    quant <- imap(stats_def, \(i, nm) {
      k <- match(nm, names(stats_def))
      map(unique(sm$condition == "control"), \(is_ctrl) {
        cols <- which((sm$condition == "control") == is_ctrl)
        rq <- apply(rep_stat[, cols, k], 1, quantile, probs)
        tibble(model = model, stat = nm, condition = if (is_ctrl) "control" else "stimulated", prob = probs,
               observed = quantile(obs[[nm]][cols], probs), rep_median = apply(rq, 1, median),
               rep_lo = apply(rq, 1, quantile, 0.05), rep_hi = apply(rq, 1, quantile, 0.95))
      }) |> list_rbind()
    }) |> list_rbind()
    list(samples = samples, quant = quant)
  }
  ppcs <- imap(stan_fits, \(f, m) ppc_one(f, m))
  ppc_samples <- map(ppcs, "samples") |> list_rbind()
  write_tsv(ppc_samples, out("ppc_samples.tsv"))
  write_tsv(map(ppcs, "quant") |> list_rbind(), out("ppc_quantiles.tsv"))
  write_tsv(ppc_samples |>
              mutate(cond = if_else(condition == "control", "control", "stimulated")) |>
              summarise(n = n(), coverage90 = mean(observed >= rep_q05 & observed <= rep_q95),
                        frac_ppp_extreme = mean(ppp < 0.025 | ppp > 0.975), mean_ppp = mean(ppp),
                        obs_total = sum(observed), rep_total = sum(rep_mean), .by = c(model, stat, cond)),
            out("ppc_summary.tsv"))
  write_session_info(out("session_info_fit.txt"), PKG_DIR)
}

# Stage compass ----------------------------------------------------------------
run_compass <- function() {
  x <- readRDS(file.path(in_dir, "cj_data.rds"))
  sec <- system.time(cc <- cj_compass(x, seed = SEED))[["elapsed"]] # default COMPASS filter
  stopifnot("COMPASS did not fit" = length(cc$fits) == 1)
  write_tsv(cc$scores |>
              transmute(donor_id, FS = compass_FS, PFS = compass_PFS, Response = compass_response_degree1,
                        Response2 = compass_response_degree2) |>
              pivot_longer(-donor_id, names_to = "score", values_to = "value") |> mutate(method = "COMPASS"),
            out("compass_scores.tsv"))
  write_tsv(cc$filters, out("compass_filters.tsv"))
  write_tsv(cc$groups, out("compass_groups.tsv"))
  write_tsv(tibble(method = "COMPASS", step = "cj_compass", seconds = sec, cpus = 1L), out("runtime_compass.tsv"))
  if (ds == "rv144") { # the paper's own 6-marker COMPASS fit (50000 iterations); reference only
    env <- new.env()
    load(RDA, env)
    r <- env$rv144
    write_tsv(tibble(donor_id = rownames(r$data$n_s), FS = unname(COMPASS::FunctionalityScore(r)[rownames(r$data$n_s)]),
                     PFS = unname(COMPASS::PolyfunctionalityScore(r)[rownames(r$data$n_s)])), out("compass_native.tsv"))
  }
}

# Stage polyics ----------------------------------------------------------------
run_polyics <- function() {
  suppressPackageStartupMessages(library(polyICSFlow))
  counts <- read_tsv(file.path(in_dir, "counts_long.tsv"), col_types = cols(pattern_key = "c", donor_id = "c"))
  keys <- unique(counts$pattern_key)
  markers <- readRDS(file.path(in_dir, "cj_data.rds"))$markers
  # assignMarkerCombinations() is a function of the pattern only, so it runs on the distinct patterns and the
  # label is repeated by the cell counts: identical to expanding every cell first (07), at a fraction of the memory.
  bits_u <- keys |> strsplit("") |> map(\(b) b == "1") |> reduce(rbind)
  colnames(bits_u) <- markers
  sec_assign <- system.time(comb_u <- assignMarkerCombinations(as_tibble(bits_u), mode = "simple"))[["elapsed"]]
  lab <- as.character(comb_u$MarkerComb)[match(counts$pattern_key, keys)]
  sec_expand <- system.time({
    i <- rep(seq_len(nrow(counts)), counts$n)
    cells <- tibble(donor_id = counts$donor_id[i], condition = counts$condition[i], cell_type = "parent",
                    MarkerComb = factor(lab[i], levels = unique(as.character(comb_u$MarkerComb))))
    rm(i)
  })[["elapsed"]]
  sec_calc <- system.time(res <- calcPolyfunctionality(cells, md_cols = "donor_id", pop_col = "cell_type",
                                                       resolution = "MarkerComb", condition_col = "condition",
                                                       background_val = "control"))[["elapsed"]]
  net <- res |>
    filter(condition == "control") |>
    select(donor_id, MarkerComb, perc_background = perc_pos) |>
    inner_join(res |> filter(condition != "control"), join_by(donor_id, MarkerComb), relationship = "one-to-one") |>
    mutate(polyICSFlow_clipped = perc_pos_backgroundSubtracted / 100, polyICSFlow_unclipped = (perc_pos - perc_background) / 100) |>
    pivot_longer(c(polyICSFlow_clipped, polyICSFlow_unclipped), names_to = "method", values_to = "value") |>
    summarise(value = sum(value), .by = c(method, donor_id)) |>
    mutate(score = "net_sum")
  write_tsv(net, out("polyics_scores.tsv"))
  write_tsv(tibble(method = "polyICSFlow", step = c("assign_patterns", "expand_to_cells", "calcPolyfunctionality"),
                   seconds = c(sec_assign, sec_expand, sec_calc), cpus = 1L), out("runtime_polyics.tsv"))
  writeLines(c(paste("polyICSFlow", as.character(packageVersion("polyICSFlow")), packageDescription("polyICSFlow")$RemoteSha),
               capture.output(sessionInfo())), out("session_info_polyics.txt"))
}

# Stage stats ------------------------------------------------------------------
cp_ci <- \(k, n) if (n == 0) c(NA_real_, NA_real_) else unname(binom.test(k, n)$conf.int)

run_stats <- function() {
  meta <- read_tsv(file.path(in_dir, "donor_meta.tsv"), show_col_types = FALSE) |>
    mutate(across(any_of(c("gender", "age")), as.factor), vaccinated = if (ds == "rv144") arm == "vaccine" else group == "vaccine")
  cj <- read_tsv(out("cj_scores.tsv"), show_col_types = FALSE) |>
    pivot_longer(c(estimate, p_positive), names_to = "kind", values_to = "value") |>
    filter(!(score == "RMS_eps1e-3" & kind == "estimate"), !is.na(value) | kind == "p_positive") |>
    transmute(method = model, donor_id, value, score = case_when(
      score == "RMS_eps1e-3" ~ "P_RMS_gt_eps1e-3",
      kind == "estimate" ~ paste0(score, "_mean"),
      .default = paste0("P_", score, "_gt_eps")))
  cps <- read_tsv(out("compass_scores.tsv"), show_col_types = FALSE) |> transmute(method = "COMPASS", score, donor_id, value)
  poly <- read_tsv(out("polyics_scores.tsv"), show_col_types = FALSE) |> select(method, score, donor_id, value)
  long <- bind_rows(cj, cps, poly) |> mutate(name = paste(method, score, sep = "_"))
  wide <- long |> select(donor_id, name, value) |> pivot_wider(names_from = name, values_from = value)
  donor_table <- meta |> left_join(wide, join_by(donor_id), relationship = "one-to-one")
  stopifnot("scores not matched one-to-one to donors" = nrow(donor_table) == nrow(meta), !anyDuplicated(long |> select(name, donor_id)))
  write_tsv(donor_table, out("donor_table.tsv"))
  # Score columns by role. CCS is experimental and never drives a conclusion.
  cj_methods <- c("m1dm", "m1", "C0")
  det_cols <- c(unlist(map(cj_methods, \(m) paste0(m, c("_P_RMS_gt_eps", "_RMS_mean", "_P_RMS_gt_eps1e-3", "_CPS_mean", "_P_CPS_gt_eps", "_CCS_mean", "_P_CCS_gt_eps")))),
                "COMPASS_FS", "COMPASS_PFS", "COMPASS_Response", "COMPASS_Response2",
                "polyICSFlow_clipped_net_sum", "polyICSFlow_unclipped_net_sum")
  det_cols <- intersect(det_cols, names(donor_table))
  truth <- donor_table$vaccinated

  # Detection AUC -------------------------------------------------------------
  auc <- map(det_cols, \(cn) auc_boot(donor_table[[cn]], truth, BOOT, SEED) |> mutate(score_name = cn)) |> list_rbind() |>
    mutate(dataset = ds, experimental = grepl("CCS", score_name), .before = 1)
  write_tsv(auc, out("auc.tsv"))
  roc <- map(det_cols, \(cn) {
    s <- donor_table[[cn]]; thr <- sort(unique(s[!is.na(s)]), decreasing = TRUE)
    tibble(score_name = cn, threshold = c(Inf, thr),
           tpr = map_dbl(c(Inf, thr), \(t) mean(s[truth & !is.na(s)] >= t)),
           fpr = map_dbl(c(Inf, thr), \(t) mean(s[!truth & !is.na(s)] >= t)))
  }) |> list_rbind() |> mutate(dataset = ds, .before = 1)
  write_tsv(roc, out("roc_curves.tsv"))
  # Paired bootstrap of AUC differences, headline contrasts.
  pairs <- tribble(~a, ~b,
    "m1dm_RMS_mean", "COMPASS_FS", "m1dm_RMS_mean", "COMPASS_PFS", "m1dm_P_RMS_gt_eps", "COMPASS_Response",
    "m1dm_RMS_mean", "m1_RMS_mean", "m1dm_RMS_mean", "C0_RMS_mean", "m1dm_RMS_mean", "polyICSFlow_clipped_net_sum")
  pairs <- pairs |> filter(a %in% names(donor_table), b %in% names(donor_table))
  set.seed(SEED)
  i1 <- which(truth); i0 <- which(!truth)
  bi <- replicate(BOOT, c(sample(i1, replace = TRUE), sample(i0, replace = TRUE)), simplify = FALSE)
  write_tsv(pmap(pairs, \(a, b) {
    d <- map_dbl(bi, \(i) auc_rank(donor_table[[a]][i], truth[i]) - auc_rank(donor_table[[b]][i], truth[i]))
    tibble(dataset = ds, a, b, diff = auc_rank(donor_table[[a]], truth) - auc_rank(donor_table[[b]], truth),
           lower = quantile(d, 0.025, na.rm = TRUE), upper = quantile(d, 0.975, na.rm = TRUE))
  }) |> list_rbind(), out("auc_diff.tsv"))

  # False positives at each method's conventional rule -------------------------
  rules <- tribble(~method, ~score_name, ~cutoff, ~rule,
    "m1dm", "m1dm_P_RMS_gt_eps", RULE_CJ, "P(RMS > eps) > 0.95",
    "m1", "m1_P_RMS_gt_eps", RULE_CJ, "P(RMS > eps) > 0.95",
    "C0", "C0_P_RMS_gt_eps", RULE_CJ, "P(RMS > eps) > 0.95",
    "COMPASS", "COMPASS_Response", RULE_COMPASS, "Response (degree 1) > 0.5",
    "polyICSFlow_clipped", "polyICSFlow_clipped_net_sum", RULE_POLY, "clipped net > 0",
    "polyICSFlow_unclipped", "polyICSFlow_unclipped_net_sum", RULE_POLY, "unclipped net > 0 (reference, not a conventional rule)") |>
    filter(score_name %in% names(donor_table))
  write_tsv(pmap(rules, \(method, score_name, cutoff, rule) {
    called <- donor_table[[score_name]] > cutoff
    n0 <- sum(!truth & !is.na(called)); k0 <- sum(called[!truth], na.rm = TRUE)
    n1 <- sum(truth & !is.na(called)); k1 <- sum(called[truth], na.rm = TRUE)
    tibble(dataset = ds, method, rule, n_placebo = n0, n_false_positive = k0, fpr = k0 / n0,
           fpr_lower = cp_ci(k0, n0)[1], fpr_upper = cp_ci(k0, n0)[2],
           n_vaccine = n1, n_called = k1, tpr = k1 / n1, tpr_lower = cp_ci(k1, n1)[1], tpr_upper = cp_ci(k1, n1)[2])
  }) |> list_rbind(), out("fpr.tsv"))

  # Agreement -----------------------------------------------------------------
  agree_cols <- intersect(c("m1dm_RMS_mean", "m1dm_P_RMS_gt_eps", "m1dm_CPS_mean", "m1_RMS_mean", "C0_RMS_mean", "COMPASS_FS",
                            "COMPASS_PFS", "COMPASS_Response", "polyICSFlow_clipped_net_sum", "polyICSFlow_unclipped_net_sum"),
                          names(donor_table))
  write_tsv(map(c(all = "all", vaccinated = "vaccinated"), \(st) {
    d <- if (st == "all") donor_table else donor_table |> filter(vaccinated)
    expand_grid(a = agree_cols, b = agree_cols) |> filter(a < b) |>
      mutate(stratum = st, n = nrow(d), spearman = map2_dbl(a, b, \(u, v) suppressWarnings(cor(d[[u]], d[[v]], method = "spearman", use = "complete.obs"))))
  }) |> list_rbind() |> mutate(dataset = ds, .before = 1), out("agreement.tsv"))

  # Regimen contrast (HVTN078): vaccinees only, T1-T4 ----------------------------
  if (ds == "hvtn078") {
    reg_cols <- intersect(c("m1dm_RMS_mean", "m1dm_CPS_mean", "C0_RMS_mean", "m1_RMS_mean", "COMPASS_FS", "COMPASS_PFS", "polyICSFlow_clipped_net_sum"), names(donor_table))
    v <- donor_table |> filter(vaccinated)
    write_tsv(map(reg_cols, \(cn) v |> summarise(n = n(), mean = mean(.data[[cn]], na.rm = TRUE), median = median(.data[[cn]], na.rm = TRUE), .by = arm) |>
                    mutate(score_name = cn)) |> list_rbind() |> arrange(score_name, arm), out("regimen_means.tsv"))
    arms <- sort(unique(v$arm))
    write_tsv(map(reg_cols, \(cn) {
      kw <- kruskal.test(v[[cn]] ~ factor(v$arm))
      pw <- combn(arms, 2, simplify = FALSE) |> map(\(ab) {
        w <- suppressWarnings(wilcox.test(v[[cn]][v$arm == ab[1]], v[[cn]][v$arm == ab[2]], exact = FALSE))
        tibble(contrast = paste(ab[2], "vs", ab[1]), p = w$p.value,
               auc_b_gt_a = auc_rank(v[[cn]][v$arm %in% ab], v$arm[v$arm %in% ab] == ab[2]))
      }) |> list_rbind() |> mutate(p_holm = p.adjust(p, "holm"))
      pw |> mutate(score_name = cn, kruskal_p = kw$p.value, .before = 1)
    }) |> list_rbind() |> mutate(dataset = ds, .before = 1), out("regimen.tsv"))
  }

  # RV144 infection association ----------------------------------------------
  if (ds == "rv144") run_infection(donor_table)

  # Runtime -------------------------------------------------------------------
  rt <- bind_rows(
    read_tsv(out("runtime_fit.tsv"), show_col_types = FALSE) |> rename(method = model),
    read_tsv(out("runtime_compass.tsv"), show_col_types = FALSE),
    read_tsv(out("runtime_polyics.tsv"), show_col_types = FALSE))
  write_tsv(rt |> mutate(dataset = ds, n_donors = nrow(meta), .before = 1), out("runtime.tsv"))
}

run_infection <- function(donor_table) {
  v <- donor_table |> filter(vaccinated, !is.na(infected))
  score_cols <- c(m1dm_RMS = "m1dm_RMS_mean", m1dm_CPS = "m1dm_CPS_mean", m1dm_CCS = "m1dm_CCS_mean", m1_RMS = "m1_RMS_mean",
                  m1_CPS = "m1_CPS_mean", C0_RMS = "C0_RMS_mean", COMPASS_FS = "COMPASS_FS", COMPASS_PFS = "COMPASS_PFS",
                  polyICSFlow_clipped = "polyICSFlow_clipped_net_sum")
  score_cols <- score_cols[score_cols %in% names(v)]
  fit_or <- \(df, xcol, adj) {
    cov <- COVARIATES[[adj]]
    df <- df |> select(infected, all_of(c(xcol, cov))) |> drop_na() |> mutate(z = as.numeric(scale(.data[[xcol]])))
    f <- stats::glm(reformulate(c("z", cov), "infected"), binomial, df)
    list(beta = unname(coef(f)["z"]), se = sqrt(vcov(f)["z", "z"]), n = nrow(df), events = sum(df$infected))
  }
  est <- expand_grid(label = names(score_cols), adj = names(COVARIATES)) |>
    mutate(r = map2(label, adj, \(l, a) fit_or(v, score_cols[[l]], a)),
           beta = map_dbl(r, "beta"), se = map_dbl(r, "se"), n = map_int(r, "n"), events = map_int(r, "events"),
           or_per_sd = exp(beta), lower = exp(beta - 1.96 * se), upper = exp(beta + 1.96 * se),
           p = 2 * pnorm(-abs(beta / se))) |>
    select(-r) |>
    mutate(family = label %in% HEADLINE & adj == "M1", experimental = label == "m1dm_CCS",
           p_holm = if_else(family, NA_real_, NA_real_))
  fam <- est$family
  est$p_holm[fam] <- p.adjust(est$p[fam], "holm")
  write_tsv(est |> mutate(dataset = ds, .before = 1), out("infection.tsv"))

  # Propagate posterior uncertainty: refit on each thinned draw, pool by Rubin's rules.
  dr <- read_tsv(out("draws_prop.tsv"), show_col_types = FALSE)
  prop <- expand_grid(model = intersect(c("m1dm", "m1"), unique(dr$model)), score = c("RMS", "CPS"), adj = c("M0", "M1")) |>
    mutate(r = pmap(list(model, score, adj), \(m, s, a) {
      dd <- dr |> filter(model == m, score == s) |>
        inner_join(v |> select(donor_id, infected, all_of(unique(unlist(COVARIATES)))), join_by(donor_id), relationship = "many-to-one")
      per <- dd |> group_split(draw) |> map(\(g) tryCatch(fit_or(g |> mutate(sc = value), "sc", a), error = \(e) NULL)) |> compact()
      b <- map_dbl(per, "beta"); u <- map_dbl(per, \(p) p$se^2); D <- length(b)
      ubar <- mean(u); bvar <- var(b); tvar <- ubar + (1 + 1 / D) * bvar
      rr <- (1 + 1 / D) * bvar / ubar; dfr <- (D - 1) * (1 + 1 / rr)^2
      tibble(draws_used = D, beta = mean(b), se_total = sqrt(tvar), fraction_missing_info = (1 + 1 / D) * bvar / tvar,
             n = per[[1]]$n, events = per[[1]]$events, p = 2 * pt(-abs(mean(b) / sqrt(tvar)), dfr),
             or_per_sd = exp(mean(b)), lower = exp(mean(b) - qt(0.975, dfr) * sqrt(tvar)), upper = exp(mean(b) + qt(0.975, dfr) * sqrt(tvar)))
    })) |> unnest(r)
  write_tsv(prop |> mutate(dataset = ds, .before = 1), out("infection_propagated.tsv"))
}

# Run --------------------------------------------------------------------------
switch(stage, fit = run_fit(), compass = run_compass(), polyics = run_polyics(), stats = run_stats())
