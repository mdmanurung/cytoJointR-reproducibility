#!/usr/bin/env Rscript
# Purpose: score one cytoJointR candidate on the frozen benchmark (benchmark/definition.tsv) and append one row per
#          metric to benchmark/leaderboard.tsv. A candidate is a package worktree, a cj_fit argument set and a label.
#          COMPASS and polyICSFlow are never rerun: their per-donor scores are read from the cached runs.
#          Stages fit and score the candidate and write raw scores into <run_dir>/<label>/:
#            rv144     fit all 262 RV144 donors once (the model never sees arm), score both levels, all donors
#            tb_s, tb_u  fit and score the TB data sets
#            sim       fit and score the simulation guards of one scenario (--scenario=), seeds of the design
#          Stage collect reads the raw scores, computes every dev metric and the win verdicts, and appends to the
#          leaderboard. Mode --test computes the RV144-test metrics, once per label, and logs the access.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/70_run-scorecard.R <stage> <run_dir> \
#              --pkg=<worktree> --fit=model=m1dm,level=sample --label=<label> [--scenario=<id>]
#          Rscript scripts/70_run-scorecard.R collect <run_dir> --pkg=<worktree> --fit=... --label=...
#          Rscript scripts/70_run-scorecard.R --test <run_dir> --pkg=<worktree> --fit=... --label=...
#          Smoke test only: --max-seeds=2 --leaderboard=<scratch file>. Slurm: scripts/slurm_scorecard.sh.
# Inputs:
#   benchmark/definition.tsv, rv144_split.tsv, sim_scenarios.tsv   frozen benchmark [69_define-scorecard.R]
#   data/processed/{rv144,tb_s,tb_u}/                              cj_data.rds, counts.tsv, meta.tsv [20, 22]
#   runs/2026-11-03_realdata/vaccine/rv144/donor_table.tsv         cached COMPASS and polyICSFlow RV144 scores [30]
#   runs/2026-11-03_realdata/tb/{tb_s,tb_u}/donor_wide.tsv         cached COMPASS and polyICSFlow TB scores [31]
#   runs/2026-11-02_redesign/reps/                                 cached COMPASS and polyICSFlow simulation scores [06, 07]
#   <pkg>/                                                         package source, loaded with pkgload [frozen worktree]
#   scripts/helper01_true-scores.R, helper02_benchmark-sim.R, helper04_realdata-tb.R, helper06_scorecard.R
# Outputs (in <run_dir>/<label>/):
#   rv144_scores.tsv, tb_s_scores.tsv, tb_u_scores.tsv   level, donor_id, score (RMS, CPS, CCS), estimate, lower, upper, p_positive, status
#   fit_<stage>.tsv                                      diagnostics and seconds of each fit
#   sim/<scenario>/seed_NN.tsv                           level, donor, score, estimate, bounds, truth_sample, truth_struct, responder_sample
#   session_info_<stage>.txt
#   benchmark/leaderboard.tsv (appended), benchmark/test_access.tsv (appended by --test)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
SEED <- 20261103L # model seed of every real-data fit, as 30_analyse-vaccine.R; simulation fits use the data seed, as 06
BOOT <- 2000L # donor bootstrap replicates, as the definition
SEED_BOOT <- 20261103L # bootstrap seed, as the definition
RULE_CJ <- 0.95 # fixed rule: call a donor when P(RMS > eps) > this
EPSILON <- c(
  RMS = 1e-4,
  poly_absolute = 1e-4,
  poly_marginal = 1e-4,
  CPS = 0.01,
  CCS = 1e-4,
  breadth = 0.01,
  marker_magnitude = 1e-4
) # as 06
COVERAGE_MIN <- 0.85 # coverage rule of the definition
COVERAGE_MAX <- 0.97
COVERAGE_SE_MULT <- 2
INFECTION_COV <- c("gender", "age", "risk.medium", "risk.high") # model M1 of 30_analyse-vaccine.R
COLUMNS <- c(
  "timestamp",
  "candidate",
  "package_sha",
  "fit_args",
  "definition_sha",
  "mode",
  "metric_id",
  "dataset",
  "tier",
  "level_used",
  "value",
  "lower",
  "upper",
  "compass_value",
  "polyics_value",
  "diff_compass",
  "diff_compass_lower",
  "diff_compass_upper",
  "diff_polyics",
  "diff_polyics_lower",
  "diff_polyics_upper",
  "pass",
  "p",
  "note"
)
SIM_RUN <- "runs/2026-11-02_redesign" # cached comparator replicates
DETECTION_SCENARIOS <- c("rare_weak_J3", "rare_weak_J4") # as 69

# Arguments --------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
flag <- \(key, default = NULL) {
  hit <- str_subset(args, fixed(paste0("--", key, "=")))
  if (length(hit)) {
    str_remove(hit[1], fixed(paste0("--", key, "=")))
  } else {
    default
  }
}
test_mode <- "--test" %in% args
pos <- args[!str_starts(args, "--")]
stage <- if (test_mode) "test" else pos[1]
stopifnot(
  "usage: Rscript 70_run-scorecard.R <rv144|tb_s|tb_u|sim|collect> <run_dir> --pkg= --fit= --label= | --test <run_dir> ..." = stage %in%
    c("rv144", "tb_s", "tb_u", "sim", "collect", "test") &&
    length(pos) == (if (test_mode) 1 else 2),
  "--pkg, --fit and --label are required" = !is.null(flag("pkg")) &&
    !is.null(flag("fit")) &&
    !is.null(flag("label"))
)
run_dir <- pos[length(pos)]
label <- flag("label")
pkg_dir <- flag("pkg")
fit_str <- flag("fit")
fit_args <- fit_str |>
  str_split_1(",") |>
  str_split_fixed("=", 2) |>
  (\(m) set_names(map(m[, 2], \(v) type.convert(v, as.is = TRUE)), m[, 1]))()
fit_level <- fit_args$level %||% "sample" # the level the candidate's own cj_fit(level = ) selects
max_seeds <- as.integer(flag("max-seeds", NA))
leaderboard <- flag("leaderboard", "benchmark/leaderboard.tsv")
cand_dir <- file.path(run_dir, label)
dir.create(cand_dir, recursive = TRUE, showWarnings = FALSE)
out <- \(f) file.path(cand_dir, f)
pkg_sha <- system2(
  "git",
  c("-C", pkg_dir, "rev-parse", "--short", "HEAD"),
  stdout = TRUE
) |>
  paste0(
    if (
      length(system2(
        "git",
        c("-C", pkg_dir, "status", "--porcelain", "--untracked-files=no"),
        stdout = TRUE
      ))
    ) {
      "-dirty"
    } else {
      ""
    }
  )
source("scripts/helper06_scorecard.R")

# Scoring helpers --------------------------------------------------------------
# Global RMS, CPS and CCS per contrast from one set of draws.
score_draws <- \(d) {
  cytoJointR::cj_scores(d, epsilon = EPSILON) |>
    filter(is.na(marker), score %in% c("RMS", "CPS", "CCS")) |>
    select(
      contrast_id,
      donor_id,
      score,
      estimate,
      lower,
      upper,
      p_positive,
      status
    )
}
# Both levels from one fit; m1 has identical draws, so score once.
score_fit <- \(fit) {
  lv <- list(sample = fit$draws_sample, structured = fit$draws_structured)
  if (identical(lv$sample, lv$structured)) {
    lv$structured <- NULL
  }
  s <- map(lv, score_draws)
  if (is.null(lv$structured)) {
    s$structured <- s$sample
  }
  list_rbind(s[c("sample", "structured")], names_to = "level")
}
fit_candidate <- \(x, seed) {
  sec <- system.time(
    fit <- do.call(
      cytoJointR::cj_fit,
      modifyList(list(x, seed = seed), fit_args)
    )
  )[["elapsed"]]
  list(fit = fit, seconds = sec)
}
fit_row <- \(f, id) {
  dg <- f$fit$diagnostics
  tibble(
    id,
    max_rhat = dg$max_rhat,
    min_ess_bulk = dg$min_ess_bulk,
    n_divergent = dg$n_divergent,
    min_ebfmi = min(dg$ebfmi),
    fit_seconds = f$seconds
  )
}

# Fit stages -------------------------------------------------------------------
run_fit_stage <- function() {
  stopifnot(
    "stage output exists; use a fresh run dir" = !file.exists(out(paste0(
      stage,
      "_scores.tsv"
    )))
  )
  pkgload::load_all(pkg_dir, quiet = TRUE)
  source("scripts/helper01_true-scores.R") # write_session_info()
  x <- if (stage == "rv144") {
    readRDS("data/processed/rv144/cj_data.rds")
  } else {
    source("scripts/helper04_realdata-tb.R")
    dat <- tb_read(stage)
    tb_cj(dat$counts, dat$markers)
  }
  f <- fit_candidate(x, SEED)
  write_tsv(fit_row(f, stage), out(paste0("fit_", stage, ".tsv")))
  write_tsv(score_fit(f$fit), out(paste0(stage, "_scores.tsv")))
  write_session_info(out(paste0("session_info_", stage, ".txt")), pkg_dir)
}

run_sim_stage <- function() {
  pkgload::load_all(pkg_dir, quiet = TRUE)
  source("scripts/helper01_true-scores.R") # true_table(), write_session_info()
  source("scripts/helper02_benchmark-sim.R") # bench_simulate()
  sc <- read_tsv("benchmark/sim_scenarios.tsv", show_col_types = FALSE) |>
    filter(scenario_id == flag("scenario"))
  stopifnot("--scenario is not in benchmark/sim_scenarios.tsv" = nrow(sc) == 1)
  dir.create(
    out(file.path("sim", sc$scenario_id)),
    recursive = TRUE,
    showWarnings = FALSE
  )
  seeds <- seq_len(min(sc$n_seeds, max_seeds, na.rm = TRUE))
  fits <- map(seeds, \(seed) {
    path <- out(file.path(
      "sim",
      sc$scenario_id,
      sprintf("seed_%02d.tsv", seed)
    )) # a finished seed is kept: a timed-out task can be resubmitted
    if (file.exists(path)) {
      return(NULL)
    }
    data_seed <- sc$seed_offset + seed
    sim <- bench_simulate(sc, data_seed)
    truth <- true_table(sim) |>
      filter(grain == "global", score %in% c("RMS", "CPS", "CCS")) |>
      select(contrast_id, donor_id, score, truth_sample, truth_struct)
    f <- fit_candidate(sim$data, data_seed)
    score_fit(f$fit) |>
      select(-donor_id) |>
      right_join(
        truth,
        join_by(contrast_id, score),
        relationship = "many-to-one"
      ) |>
      mutate(
        scenario_id = sc$scenario_id,
        seed = seed,
        responder_sample = truth_sample > sc$eps_response,
        fit_seconds = f$seconds
      ) |>
      write_tsv(path)
  })
  write_session_info(
    out(sprintf("session_info_sim_%s.txt", sc$scenario_id)),
    pkg_dir
  )
}

# Collect: metrics -------------------------------------------------------------
def_sha <- system2(
  "git",
  c("log", "-1", "--format=%H", "--", "benchmark/definition.tsv"),
  stdout = TRUE
)
append_rows <- \(rows, mode) {
  rows <- rows |>
    mutate(
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      candidate = label,
      package_sha = pkg_sha,
      fit_args = fit_str,
      definition_sha = def_sha,
      mode = mode
    )
  rows[setdiff(COLUMNS, names(rows))] <- NA
  rows <- rows |> select(all_of(COLUMNS))
  write_tsv(
    rows,
    leaderboard,
    append = file.exists(leaderboard),
    col_names = !file.exists(leaderboard)
  )
}
# Wide donor scores of one level: RMS_mean, P_RMS, CPS_mean, P_CPS, CCS_mean, P_CCS.
wide_scores <- \(scores, lv) {
  scores |>
    filter(level == lv) |>
    pivot_wider(
      id_cols = donor_id,
      names_from = score,
      values_from = c(estimate, p_positive)
    ) |>
    transmute(
      donor_id,
      RMS_mean = estimate_RMS,
      P_RMS = p_positive_RMS,
      CPS_mean = estimate_CPS,
      P_CPS = p_positive_CPS,
      CCS_mean = estimate_CCS,
      P_CCS = p_positive_CCS
    )
}
cp_ci <- \(k, n) unname(binom.test(k, n)$conf.int)
# One detection block: own, paired-vs-COMPASS and paired-vs-polyICSFlow bootstraps of one score vector.
det_block <- \(a, truth, compass, polyics, full) {
  own <- boot_det(a, truth, B = BOOT, seed = SEED_BOOT, full = full)
  dc <- boot_det(a, truth, compass, B = BOOT, seed = SEED_BOOT, full = full)
  dp <- boot_det(a, truth, polyics, B = BOOT, seed = SEED_BOOT, full = full)
  own |>
    select(stat, value, lower, upper) |>
    left_join(
      dc |>
        select(
          stat,
          diff_compass = value,
          diff_compass_lower = lower,
          diff_compass_upper = upper
        ),
      join_by(stat),
      relationship = "one-to-one"
    ) |>
    left_join(
      dp |>
        select(
          stat,
          diff_polyics = value,
          diff_polyics_lower = lower,
          diff_polyics_upper = upper
        ),
      join_by(stat),
      relationship = "one-to-one"
    )
}

# Rows of one detection data set. d: donor_id, truth, compass, polyics, and the candidate columns of wide_scores().
detection_rows <- \(def, ds_id, d, cand_by_level) {
  rows <- def |>
    filter(
      dataset == ds_id,
      role %in% c("dev", "test"),
      str_detect(metric_id, "\\.(auc|pauc10|tpr_fpr05|tpr_fpr10)\\.")
    )
  memo <- list()
  pmap(rows, \(...) {
    r <- list(...)
    lv <- if (r$level == "fit") fit_level else r$level
    stat <- str_split_1(r$metric_id, "\\.")[2]
    full <- r$tier %in% c("primary", "companion", "test_only")
    key <- paste(lv, r$score, full)
    if (is.null(memo[[key]])) {
      dd <- d |>
        inner_join(
          cand_by_level[[lv]] |> select(donor_id, a = all_of(r$score)),
          join_by(donor_id),
          relationship = "one-to-one"
        )
      memo[[key]] <<- det_block(dd$a, dd$truth, dd$compass, dd$polyics, full)
    }
    memo[[key]] |>
      filter(stat == !!stat) |>
      mutate(
        metric_id = r$metric_id,
        dataset = r$dataset,
        tier = r$tier,
        level_used = lv,
        compass_value = r$compass_value,
        polyics_value = r$polyics_value,
        note = if_else(
          is.na(value),
          "no value: every donor score is NA (CCS support gate)",
          NA_character_
        )
      )
  }) |>
    list_rbind()
}
# Placebo FPR and vaccinee TPR at the fixed rule.
rule_rows <- \(def, ds_id, d, cand_by_level) {
  def |>
    filter(dataset == ds_id, str_detect(metric_id, "_rule\\.")) |>
    pmap(\(...) {
      r <- list(...)
      dd <- d |>
        inner_join(
          cand_by_level[[fit_level]] |> select(donor_id, a = P_RMS),
          join_by(donor_id),
          relationship = "one-to-one"
        )
      fpr <- str_detect(r$metric_id, "fpr_rule")
      called <- dd$a > RULE_CJ
      k <- sum(called[if (fpr) !dd$truth else dd$truth], na.rm = TRUE)
      n <- sum(if (fpr) !dd$truth else dd$truth)
      tibble(
        metric_id = r$metric_id,
        dataset = r$dataset,
        tier = r$tier,
        level_used = fit_level,
        value = k / n,
        lower = cp_ci(k, n)[1],
        upper = cp_ci(k, n)[2],
        compass_value = r$compass_value,
        polyics_value = r$polyics_value,
        note = sprintf("%d of %d", k, n)
      )
    }) |>
    list_rbind()
}

run_collect <- function() {
  stopifnot(
    "benchmark definition or split has uncommitted changes" = !length(system2(
      "git",
      c(
        "status",
        "--porcelain",
        "benchmark/definition.tsv",
        "benchmark/rv144_split.tsv",
        "benchmark/sim_scenarios.tsv"
      ),
      stdout = TRUE
    ))
  )
  def <- read_tsv("benchmark/definition.tsv", show_col_types = FALSE)
  split <- read_tsv("benchmark/rv144_split.tsv", show_col_types = FALSE)
  read_scores <- \(ds) {
    read_tsv(out(paste0(ds, "_scores.tsv")), show_col_types = FALSE)
  }
  by_level <- \(scores) {
    set_names(c("sample", "structured")) |> map(\(lv) wide_scores(scores, lv))
  }

  # RV144 dev
  rv_cached <- read_tsv(
    "runs/2026-11-03_realdata/vaccine/rv144/donor_table.tsv",
    show_col_types = FALSE
  ) |>
    transmute(
      donor_id,
      truth = vaccinated,
      compass = COMPASS_FS,
      compass_resp = COMPASS_Response,
      polyics = polyICSFlow_clipped_net_sum
    ) |>
    inner_join(
      split |> select(donor_id, split),
      join_by(donor_id),
      relationship = "one-to-one"
    )
  rv_scores <- read_scores("rv144")
  stopifnot(
    "RV144 candidate scores do not cover the cached donors" = setequal(
      rv_scores$donor_id,
      rv_cached$donor_id
    )
  )
  rv_dev <- rv_cached |> filter(split == "dev")
  rv_by <- by_level(rv_scores)
  rows_rv <- bind_rows(
    detection_rows(def, "rv144_dev", rv_dev, rv_by),
    rule_rows(def, "rv144_dev", rv_dev, rv_by)
  )

  # TB
  rows_tb <- map(c("tb_s", "tb_u"), \(ds) {
    t <- read_tsv(
      file.path("runs/2026-11-03_realdata/tb", ds, "donor_wide.tsv"),
      show_col_types = FALSE
    ) |>
      transmute(
        donor_id,
        truth = case,
        compass = `COMPASS FS`,
        polyics = `polyICSFlow_clipped net_sum`
      )
    detection_rows(def, ds, t, by_level(read_scores(ds)))
  }) |>
    list_rbind()

  # Simulation guards
  sc <- read_tsv("benchmark/sim_scenarios.tsv", show_col_types = FALSE)
  seeds <- seq_len(min(sc$n_seeds, max_seeds, na.rm = TRUE))
  sim <- sc$scenario_id |>
    map(\(s) {
      map(seeds, \(i) {
        read_tsv(
          out(file.path("sim", s, sprintf("seed_%02d.tsv", i))),
          show_col_types = FALSE
        )
      }) |>
        list_rbind()
    }) |>
    list_rbind() |>
    filter(level == fit_level)
  cov_seed <- sim |>
    filter(score %in% c("RMS", "CPS")) |>
    mutate(truth = if (fit_level == "sample") truth_sample else truth_struct) |>
    summarise(
      coverage = mean(coalesce(lower <= truth & truth <= upper, FALSE)),
      .by = c(scenario_id, seed, score)
    )
  rows_cov <- cov_seed |>
    summarise(
      value = mean(coverage),
      se = sd(coverage) / sqrt(n()),
      .by = c(scenario_id, score)
    ) |>
    mutate(
      metric_id = sprintf("sim.cov90.%s.%s", score, scenario_id),
      dataset = scenario_id,
      tier = "guard",
      level_used = fit_level,
      lower = value - COVERAGE_SE_MULT * se, # pass: coverage + 2 SE >= min and coverage - 2 SE <= max, as 11
      upper = value + COVERAGE_SE_MULT * se,
      pass = upper >= COVERAGE_MIN &
        (lower <= COVERAGE_MAX | str_detect(scenario_id, "^null")),
      note = sprintf(
        "coverage %.3f, SE %.3f over %d seeds",
        value,
        se,
        length(seeds)
      )
    ) |>
    select(
      metric_id,
      dataset,
      tier,
      level_used,
      value,
      lower,
      upper,
      pass,
      note
    )
  cmp <- sim_comparators(SIM_RUN, DETECTION_SCENARIOS, seeds)
  cand <- sim |>
    filter(scenario_id %in% DETECTION_SCENARIOS, score == "RMS") |>
    transmute(
      scenario_id,
      seed,
      donor_id,
      a_p = p_positive,
      a_rms = estimate,
      truth_live = truth_sample
    )
  sim_tab <- cmp |>
    inner_join(
      cand,
      join_by(scenario_id, seed, donor_id),
      relationship = "one-to-one"
    )
  stopifnot(
    "simulated truth differs from the cached truth: the generator changed" = isTRUE(all.equal(
      sim_tab$true_rms_sample,
      sim_tab$truth_live
    ))
  )
  # One sim metric: own statistic, and the paired difference to a comparator column (resampling seeds).
  sim_metric <- \(id, stat, truth, a, vs_compass, vs_polyics) {
    d <- sim_tab |>
      transmute(
        seed = paste(scenario_id, seed),
        truth = {{ truth }},
        a = {{ a }},
        c = {{ vs_compass }},
        p = {{ vs_polyics }}
      )
    own <- seed_boot(d, stat, B = BOOT / 2, seed = SEED_BOOT)
    dc <- seed_boot(
      d |> mutate(b = c),
      stat,
      use_b = TRUE,
      B = BOOT / 2,
      seed = SEED_BOOT
    )
    dp <- seed_boot(
      d |> mutate(b = p),
      stat,
      use_b = TRUE,
      B = BOOT / 2,
      seed = SEED_BOOT
    )
    df <- def |> filter(metric_id == id)
    tibble(
      metric_id = id,
      dataset = df$dataset,
      tier = "guard",
      level_used = fit_level,
      value = own$value,
      lower = own$lower,
      upper = own$upper,
      compass_value = df$compass_value,
      polyics_value = df$polyics_value,
      diff_compass = dc$value,
      diff_compass_lower = dc$lower,
      diff_compass_upper = dc$upper,
      diff_polyics = dp$value,
      diff_polyics_lower = dp$lower,
      diff_polyics_upper = dp$upper
    )
  }
  rows_sim <- bind_rows(
    sim_metric(
      "sim.auc.rare_weak.P_RMS.fit",
      "auc",
      responder_sample,
      a_p,
      compass,
      polyics
    ),
    sim_metric(
      "sim.spearman.rare_weak.RMS_mean.fit",
      "spearman",
      true_rms_sample,
      a_rms,
      compass,
      polyics_mag
    )
  )

  # Win verdicts: the criteria are quoted in benchmark/definition.tsv
  rows <- bind_rows(rows_rv, rows_tb, rows_cov, rows_sim)
  g <- \(id) rows |> filter(metric_id == id)
  guards_pass <- all(rows_cov$pass)
  c_checks <- c(
    auc_lower_gt0 = g("rv144_dev.auc.P_RMS.fit")$diff_compass_lower > 0,
    pauc_gt0 = g("rv144_dev.pauc10.P_RMS.fit")$diff_compass > 0,
    tpr10_ge0 = g("rv144_dev.tpr_fpr10.P_RMS.fit")$diff_compass >= 0,
    sim_auc_ge0 = g("sim.auc.rare_weak.P_RMS.fit")$diff_compass >= 0,
    guards = guards_pass
  )
  p_checks <- c(
    auc_lower_gt0 = g("rv144_dev.auc.P_RMS.fit")$diff_polyics_lower > 0,
    sim_magnitude_lower_gt0 = g(
      "sim.spearman.rare_weak.RMS_mean.fit"
    )$diff_polyics_lower >
      0,
    sim_auc_ge0 = g("sim.auc.rare_weak.P_RMS.fit")$diff_polyics >= 0,
    guards = guards_pass
  )
  verdict <- \(id, checks) {
    tibble(
      metric_id = id,
      dataset = "all",
      tier = "win",
      level_used = fit_level,
      value = as.numeric(all(checks)),
      pass = all(checks),
      note = paste(names(checks), checks, sep = "=", collapse = "; ")
    )
  }
  rows <- bind_rows(
    rows,
    verdict("win.beats_compass", c_checks),
    verdict("win.beats_polyics", p_checks)
  )
  append_rows(rows, "dev")
  writeLines(capture.output(sessionInfo()), out("session_info_collect.txt"))
}

# Test mode --------------------------------------------------------------------
run_test <- function() {
  access <- read_tsv(
    "benchmark/test_access.tsv",
    show_col_types = FALSE,
    col_types = cols(.default = "c")
  )
  lb <- read_tsv(leaderboard, show_col_types = FALSE)
  stopifnot(
    "label already used the test set; one access per candidate" = !label %in%
      access$candidate,
    "no dev rows for this label: run collect first" = any(
      lb$candidate == label & lb$mode == "dev"
    ),
    "benchmark definition or split has uncommitted changes" = !length(system2(
      "git",
      c(
        "status",
        "--porcelain",
        "benchmark/definition.tsv",
        "benchmark/rv144_split.tsv"
      ),
      stdout = TRUE
    ))
  )
  def <- read_tsv("benchmark/definition.tsv", show_col_types = FALSE)
  split <- read_tsv("benchmark/rv144_split.tsv", show_col_types = FALSE)
  cached <- read_tsv(
    "runs/2026-11-03_realdata/vaccine/rv144/donor_table.tsv",
    show_col_types = FALSE
  ) |>
    inner_join(
      split |> select(donor_id, split),
      join_by(donor_id),
      relationship = "one-to-one"
    ) |>
    filter(split == "test") |>
    mutate(
      truth = vaccinated,
      compass = COMPASS_FS,
      polyics = polyICSFlow_clipped_net_sum
    )
  scores <- read_tsv(out("rv144_scores.tsv"), show_col_types = FALSE)
  by_level <- set_names(c("sample", "structured")) |>
    map(\(lv) wide_scores(scores, lv))
  d <- cached |> select(donor_id, truth, compass, polyics)
  rows <- bind_rows(
    detection_rows(def, "rv144_test", d, by_level),
    rule_rows(def, "rv144_test", d, by_level)
  )
  # comparator values on the test half, computed now and only now
  rows <- rows |>
    mutate(
      compass_value = case_when(
        metric_id == "rv144_test.fpr_rule.P_RMS.fit" ~ mean(
          cached$COMPASS_Response[!cached$truth] > 0.5
        ),
        metric_id == "rv144_test.tpr_rule.P_RMS.fit" ~ mean(
          cached$COMPASS_Response[cached$truth] > 0.5
        ),
        str_detect(metric_id, "\\.auc\\.") ~ auc_rank(
          cached$compass,
          cached$truth
        ),
        str_detect(metric_id, "pauc10") ~ pauc(cached$compass, cached$truth),
        str_detect(metric_id, "tpr_fpr05") ~ tpr_at_fpr(
          cached$compass,
          cached$truth,
          0.05
        ),
        str_detect(metric_id, "tpr_fpr10") ~ tpr_at_fpr(
          cached$compass,
          cached$truth,
          0.10
        )
      ),
      polyics_value = case_when(
        metric_id == "rv144_test.fpr_rule.P_RMS.fit" ~ mean(
          cached$polyics[!cached$truth] > 0
        ),
        metric_id == "rv144_test.tpr_rule.P_RMS.fit" ~ mean(
          cached$polyics[cached$truth] > 0
        ),
        str_detect(metric_id, "\\.auc\\.") ~ auc_rank(
          cached$polyics,
          cached$truth
        ),
        str_detect(metric_id, "pauc10") ~ pauc(cached$polyics, cached$truth),
        str_detect(metric_id, "tpr_fpr05") ~ tpr_at_fpr(
          cached$polyics,
          cached$truth,
          0.05
        ),
        str_detect(metric_id, "tpr_fpr10") ~ tpr_at_fpr(
          cached$polyics,
          cached$truth,
          0.10
        )
      )
    )

  # Infection association among test-half vaccinees, as 30_analyse-vaccine.R model M1
  meta <- read_tsv(
    "data/processed/rv144/donor_meta.tsv",
    show_col_types = FALSE
  ) |>
    inner_join(
      split |> select(donor_id, split),
      join_by(donor_id),
      relationship = "one-to-one"
    ) |>
    filter(split == "test", arm == "vaccine", !is.na(infected)) |>
    mutate(across(any_of(c("gender", "age")), as.factor))
  or_row <- \(score_col, v) {
    df <- meta |>
      inner_join(v, join_by(donor_id), relationship = "one-to-one") |>
      select(infected, all_of(c(score_col, INFECTION_COV))) |>
      drop_na() |>
      mutate(z = as.numeric(scale(.data[[score_col]])))
    f <- stats::glm(
      reformulate(c("z", INFECTION_COV), "infected"),
      binomial,
      df
    )
    b <- unname(coef(f)["z"])
    s <- sqrt(vcov(f)["z", "z"])
    tibble(
      or = exp(b),
      lower = exp(b - 1.96 * s),
      upper = exp(b + 1.96 * s),
      p = 2 * pnorm(-abs(b / s)),
      n = nrow(df),
      events = sum(df$infected)
    )
  }
  pfs <- or_row("pfs", cached |> transmute(donor_id, pfs = COMPASS_PFS))
  rows_or <- def |>
    filter(str_detect(metric_id, "infection_or")) |>
    pmap(\(...) {
      r <- list(...)
      o <- or_row(
        r$score,
        by_level[[fit_level]] |> select(donor_id, all_of(r$score))
      )
      tibble(
        metric_id = r$metric_id,
        dataset = r$dataset,
        tier = r$tier,
        level_used = fit_level,
        value = o$or,
        lower = o$lower,
        upper = o$upper,
        p = o$p,
        compass_value = pfs$or,
        note = sprintf(
          "n %d, events %d; COMPASS PFS on the same donors: OR %.2f (%.2f to %.2f), p %.3f",
          o$n,
          o$events,
          pfs$or,
          pfs$lower,
          pfs$upper,
          pfs$p
        )
      )
    }) |>
    list_rbind()
  a <- rows |> filter(metric_id == "rv144_test.auc.P_RMS.fit")
  verdict <- \(id, lower) {
    tibble(
      metric_id = id,
      dataset = "rv144_test",
      tier = "win",
      level_used = fit_level,
      value = as.numeric(lower > 0),
      pass = lower > 0,
      note = "paired-bootstrap 95% lower bound of the AUC difference above 0 (criterion 5 / 4 of definition.tsv)"
    )
  }
  append_rows(
    bind_rows(
      rows,
      rows_or,
      verdict("win.test_beats_compass", a$diff_compass_lower),
      verdict("win.test_beats_polyics", a$diff_polyics_lower)
    ),
    "test"
  )
  write_tsv(
    tibble(
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      candidate = label,
      package_sha = pkg_sha,
      definition_sha = def_sha,
      access_number = nrow(access) + 1L
    ),
    "benchmark/test_access.tsv",
    append = TRUE
  )
  message(
    "test-set access number ",
    nrow(access) + 1L,
    " for ",
    label,
    ": update the counter in benchmark/README.md"
  )
}

# Run --------------------------------------------------------------------------
switch(
  stage,
  rv144 = ,
  tb_s = ,
  tb_u = run_fit_stage(),
  sim = run_sim_stage(),
  collect = run_collect(),
  test = run_test()
)
