#!/usr/bin/env Rscript
# Purpose: define the redesign evaluation (candidate models C0-C2phi vs COMPASS vs polyICSFlow). Solves each
#          scenario's baseline and effect size by simulation so the stated targets hold at J = 3 and J = 4,
#          checks the achieved true-RMS distribution, and writes one row per scenario for 06 and 07, the
#          replicate list for the SLURM array (several replicates per task). Stan is compiled separately
#          before submitting (see slurm_benchmark.sh), so the tasks do not race to build it.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/05_define-benchmark.R <run_dir> <n_seeds> <runtime_seeds> [<solved_run_dir>]
#          With <solved_run_dir> the solved scenarios.tsv of that run is reused (smoke runs), and nothing is solved.
# Inputs:
#   none                        targets below [this script]
#   <PKG_DIR>/                  package source, loaded with pkgload; env CJ_PKG_DIR overrides it [this repo / worktree]
#   scripts/helper02_benchmark-sim.R, lib_generator-b.R, lib_generator-c.R, lib_generator-d.R   simulators [this repo]
#   scripts/helper01_true-scores.R   write_session_info() [this repo]
# Outputs (in <run_dir>):
#   scenarios.tsv               one row per scenario: generator, J, solved alpha / beta, thresholds, design, role, n_seeds
#   jobs.tsv                    one row per replicate: task (SLURM array index), scenario_id, seed
#   tuning_check.tsv            achieved true-RMS quantiles, baseline, and responder / responding-pattern fractions
#                               per scenario in one large simulated cohort (TUNE_DONORS donors)
#   session_info.txt            R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign") # package source under test; the redesign branch holds C0-C2
J_VALUES <- c(3L, 4L) # markers per scenario, plan "Scale"
PACK <- 6L # replicates per SLURM task: about 545 replicates / 6 = 91 tasks, under the 100-task target; 6 x ~3 min fits the 1 h short limit
MISSING_FRAC <- 0.25 # missing_control scenario: fraction of donor-visits without a control (plan)
VISIT_EFFECT <- -0.8 # visit_change scenario: visit_v2:stim_CMV log-odds shift for every marker, as the pilot
TAU_U <- 0.3 # donor x visit SD in the visit_change scenario, package default truth
RUNTIME_DONORS <- 100L # donors in the runtime-only cell (plan)
N_DONORS <- 24L # donors per replicate, as agreed
N_CELLS <- 10000 # parent cells per sample, as agreed
EPS_RESPONSE <- 1e-4 # donor responds when true RMS > this: below 1/N_CELLS a sample cannot show the change
EPS_PATTERN <- 1e-5 # pattern responds when true change > this: 0.1 expected cells per sample, the least a count can see
ALPHA_STEP <- 0.3 # marker j baseline logit = alpha_1 + 0.3 (j - 1), as in the pilot
BETA_STEP <- 0.2 # marker j effect = beta_1 + 0.2 (j - 1), as in the pilot
BETA_DOWN <- -1 # negative scenario: log-odds drop for markers 3 and 4 (about 63% of their baseline)
TAU_BASAL <- 0.5 # donor basal SD, package default truth
TAU_INDUCED <- 0.2 # donor induced-effect SD; small so a "weak" scenario is weak for most donors
OMEGA_RHO <- 0.3 # donor-effect correlation, package default truth
THETA0_START <- 0.6 # baseline pairwise interaction, package default truth
RESP_FRAC <- 0.5 # mixture: fraction of donors that respond
GB_SD <- c(donor = 0.3, sample = 0.1, stim = 0.3, scale_sd = 0.5) # generator B noise, as in the pilot
TUNE_DONORS <- 240L # donors in the cohort used to solve and check targets; median over donors
TUNE_SEED <- 20261103L # seed of that cohort
SEED_OFFSET <- 20261103L # replicate data seed = SEED_OFFSET + 1000 * scenario_index + seed (06)
# Targets: baseline = median control any-positive frequency, rms = median true donor RMS.
TARGETS <- tribble(
  ~scenario_id      , ~generator , ~baseline , ~rms   ,
  "null"            , "A"        , 0.001     , 0      ,
  "rare_weak"       , "A"        , 0.001     , 0.0005 ,
  "rare_moderate"   , "A"        , 0.001     , 0.002  ,
  "moderate"        , "A"        , 0.001     , 0.01   ,
  "mixture"         , "C"        , 0.001     , 0.002  , # responders only; beta copied from rare_moderate
  "negative"        , "A"        , 0.01      , 0.002  , # markers 1-2 up (solved), 3-4 down (BETA_DOWN)
  "generator_B"     , "B"        , 0.001     , 0.01   ,
  "missing_control" , "A"        , 0.001     , 0.002  , # alpha and beta copied from rare_moderate, 25% of controls dropped
  "visit_change"    , "A"        , 0.001     , 0.002 # alpha and beta copied from rare_moderate, two visits with a visit change
)
COPY_OF_RARE_MODERATE <- c("mixture", "missing_control", "visit_change") # scenarios that reuse its alpha and beta

# Scenario skeleton ------------------------------------------------------------
# Parameters other than alpha and beta; alpha_1..4 and beta_1..4 are filled in below.
scenarios <- TARGETS |>
  rename(scenario = scenario_id) |>
  cross_join(tibble(J = J_VALUES)) |>
  mutate(
    scenario_id = paste0(scenario, "_J", J),
    scenario_index = row_number(),
    role = "main",
    n_donors = N_DONORS,
    N = N_CELLS,
    missing_control = if_else(scenario == "missing_control", MISSING_FRAC, 0),
    n_visits = if_else(scenario == "visit_change", 2L, 1L),
    visit_effect = if_else(scenario == "visit_change", VISIT_EFFECT, 0),
    tau_u = if_else(scenario == "visit_change", TAU_U, 0),
    tau_basal = TAU_BASAL,
    tau_induced = if_else(scenario == "null", 0, TAU_INDUCED), # null: exactly zero induced effect
    omega_rho = OMEGA_RHO,
    theta0_start = THETA0_START,
    resp_frac = if_else(generator == "C", RESP_FRAC, 1),
    gb_sd_donor = if_else(generator == "B", GB_SD[["donor"]], NA_real_),
    gb_sd_sample = if_else(generator == "B", GB_SD[["sample"]], NA_real_),
    gb_sd_stim = if_else(generator == "B", GB_SD[["stim"]], NA_real_),
    gb_scale_sd = if_else(generator == "B", GB_SD[["scale_sd"]], NA_real_),
    eps_response = EPS_RESPONSE,
    eps_pattern = EPS_PATTERN,
    seed_offset = SEED_OFFSET + 1000L * scenario_index,
    alpha_1 = NA_real_,
    beta_1 = 0
  )
set_alpha <- \(sc, x) {
  sc[paste0("alpha_", 1:4)] <- as.list(x + ALPHA_STEP * 0:3)
  sc
}
set_beta <- \(sc, x) {
  beta <- x + BETA_STEP * 0:3
  if (sc$scenario == "null") {
    beta <- rep(0, 4)
  }
  if (sc$scenario == "negative") {
    beta <- c(x, x, BETA_DOWN, BETA_DOWN)
  }
  sc[paste0("beta_", 1:4)] <- as.list(beta)
  sc
}

# Simulators -------------------------------------------------------------------
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper01_true-scores.R") # write_session_info()
source("scripts/helper02_benchmark-sim.R") # bench_simulate(), true_pattern_delta()
stat_baseline <- \(sim) {
  median(rowSums(sim$p[
    na.omit(sim$data$contrasts$control_sample_id),
    rowSums(sim$data$bits) >= 1
  ]))
}
stat_rms <- \(sim) {
  true_pattern_delta(sim) |>
    summarise(rms = sum(delta), .by = contrast_id) |>
    pull(rms) |>
    median()
}
# One root: the value x of a parameter (set by `set`) at which `stat` of a TUNE_DONORS cohort hits target.
solve_target <- function(sc, set, stat, target, interval) {
  uniroot(
    \(x) stat(bench_simulate(set(sc, x), TUNE_SEED, TUNE_DONORS)) - target,
    interval,
    tol = 1e-3
  )$root
}

# Solve alpha and beta ---------------------------------------------------------
# Baseline first (beta 0), then the effect at that baseline. null has beta 0; the copies reuse rare_moderate of
# the same J. With a solved run given, its scenarios.tsv replaces the solving.
args <- commandArgs(trailingOnly = TRUE)
run_dir <- args[1]
n_seeds_arg <- as.integer(args[2])
runtime_seeds <- as.integer(args[3])
solved_dir <- args[4]
stopifnot(
  "usage: Rscript 05_define-benchmark.R <run_dir> <n_seeds> <runtime_seeds> [<solved_run_dir>]" = !is.na(
    runtime_seeds
  ),
  "run_dir already exists" = !dir.exists(run_dir)
)
if (!is.na(solved_dir)) {
  solved <- read_tsv(
    file.path(solved_dir, "scenarios.tsv"),
    show_col_types = FALSE
  )
  scenarios <- solved |> filter(role == "main")
} else {
  solved_rows <- scenarios |>
    group_split(scenario_index) |>
    map(\(r) {
      r <- set_beta(r, 0)
      if (r$scenario %in% COPY_OF_RARE_MODERATE) {
        return(r)
      }
      alpha <- solve_target(r, set_alpha, stat_baseline, r$baseline, c(-14, -3))
      r <- set_alpha(r, alpha)
      if (r$scenario == "null") {
        return(r)
      }
      set_beta(r, solve_target(r, set_beta, stat_rms, r$rms, c(0, 5)))
    })
  rare_moderate <- keep(solved_rows, \(r) r$scenario == "rare_moderate") |>
    list_rbind()
  scenarios <- solved_rows |>
    map(\(r) {
      if (!r$scenario %in% COPY_OF_RARE_MODERATE) {
        return(r)
      }
      from <- filter(rare_moderate, J == r$J)
      r[paste0("alpha_", 1:4)] <- from[paste0("alpha_", 1:4)]
      r[paste0("beta_", 1:4)] <- from[paste0("beta_", 1:4)]
      r
    }) |>
    list_rbind()
}
# The runtime-only cell: rare_moderate at J = 4 with RUNTIME_DONORS donors, not part of the decision rule.
runtime_cell <- scenarios |>
  filter(scenario == "rare_moderate", J == max(J_VALUES)) |>
  mutate(
    scenario = "runtime_100",
    scenario_id = "runtime_100_J4",
    scenario_index = max(scenarios$scenario_index) + 1L,
    role = "runtime",
    n_donors = RUNTIME_DONORS
  )
scenarios <- bind_rows(scenarios, runtime_cell) |>
  mutate(
    n_seeds = if_else(role == "runtime", runtime_seeds, n_seeds_arg),
    seed_offset = SEED_OFFSET + 1000L * scenario_index
  )

# Check achieved truth ---------------------------------------------------------
# Simulate one large cohort per scenario at the solved values and summarise the true quantities.
tuning_check <- scenarios |>
  filter(role == "main") |>
  group_split(scenario_id) |>
  map(\(sc) {
    sim <- bench_simulate(sc, TUNE_SEED, TUNE_DONORS)
    delta <- true_pattern_delta(sim)
    rms <- delta |> summarise(rms = sum(delta), .by = contrast_id)
    tibble(
      scenario_id = sc$scenario_id,
      baseline_any_positive = stat_baseline(sim),
      rms_q05 = quantile(rms$rms, 0.05),
      rms_q25 = quantile(rms$rms, 0.25),
      rms_median = median(rms$rms),
      rms_q75 = quantile(rms$rms, 0.75),
      rms_q95 = quantile(rms$rms, 0.95),
      frac_responders = mean(rms$rms > EPS_RESPONSE),
      frac_patterns_up = mean(delta$delta > EPS_PATTERN),
      frac_patterns_down = mean(delta$delta < -EPS_PATTERN)
    )
  }) |>
  list_rbind() |>
  mutate(
    scenario_id = fct(scenario_id, levels = unique(scenarios$scenario_id))
  ) |>
  arrange(scenario_id)

# Replicate list ---------------------------------------------------------------
# One row per scenario x seed. Main replicates are packed PACK to a task; each runtime replicate has its own task.
jobs <- scenarios |>
  mutate(seed = map(n_seeds, seq_len), .keep = "unused") |>
  select(scenario_id, role, seed) |>
  unnest(seed)
jobs_main <- jobs |>
  filter(role == "main") |>
  mutate(task = (row_number() - 1L) %/% PACK + 1L)
jobs <- bind_rows(
  jobs_main,
  jobs |>
    filter(role == "runtime") |>
    mutate(task = max(jobs_main$task) + row_number())
) |>
  select(task, scenario_id, seed)

# Write outputs ----------------------------------------------------------------
dir.create(run_dir, recursive = TRUE)
write_tsv(scenarios, file.path(run_dir, "scenarios.tsv"))
write_tsv(jobs, file.path(run_dir, "jobs.tsv"))
write_tsv(tuning_check, file.path(run_dir, "tuning_check.tsv"))
write_session_info(file.path(run_dir, "session_info.txt"), PKG_DIR)
