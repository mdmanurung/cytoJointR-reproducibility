#!/usr/bin/env Rscript
# Purpose: define the pilot-benchmark scenarios, one row per scenario with every
#          generator parameter, so 03_run-replicate.R reads nothing but this table.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/02_define-scenarios.R runs/2026-11-02_pilot
# Inputs:
#   none                        parameters below [this script]
# Outputs (in <run_dir>):
#   scenarios.tsv               one row per scenario; all M1 / generator-B parameters
#   session_info.txt            R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
# Truths are M1 parameters on the Stan scale. Marker j (1-based) has
# alpha = ALPHA0 + 0.3 (j - 1) and stimulation effect beta + 0.2 (j - 1).
# Stimulation names are fixed to CMV and SEB, visits to v1 and v2, in 03.
TRUTHS <- tribble(
  ~truth             , ~alpha0 , ~beta_cmv , ~beta_seb , ~beta_step , ~visit_cmv , ~visit_seb , ~delta_cmv , ~delta_seb , ~tau_induced ,
  "null"             , -4.5    ,  0.0      ,  0.0      , 0.0        ,  0.0       ,  0.0       , 0.0        , 0.0        , 0.0          ,
  "positive_margins" , -4.5    ,  2.0      ,  1.6      , 0.2        ,  0.0       ,  0.0       , 0.0        , 0.0        , 0.3          ,
  "negative_margins" , -3.5    , -1.0      , -1.0      , 0.0        ,  0.0       ,  0.0       , 0.0        , 0.0        , 0.3          ,
  "dependence_only"  , -4.5    ,  0.0      ,  0.0      , 0.0        ,  0.0       ,  0.0       , 1.0        , 1.0        , 0.0          ,
  "visit_change"     , -4.5    ,  2.0      ,  1.6      , 0.2        , -0.8       , -0.4       , 0.0        , 0.0        , 0.3
)
# null and dependence_only: tau_induced 0 so the true stimulation effect is exactly zero
# (null) or comes only from the interactions (dependence_only).
# negative_margins: higher baseline (about 3% per marker) so a drop is visible at N = 10000.
# visit_cmv, visit_seb: visit x stimulation coefficients; the other truths have none.
N_CELLS <- 10000 # parent cells per sample, as agreed
MISSING_CONTROL <- 0.1 # fraction of donor-visits without a control, as in the smoke test
TAU_BASAL <- 0.5 # donor basal SD, package default truth
TAU_U <- 0.3 # donor x visit SD, package default truth
OMEGA_RHO <- 0.3 # donor-effect correlation, package default truth
THETA0_START <- 0.6 # baseline pairwise interaction, package default truth
GB_SD <- c(donor = 0.3, sample = 0.1, stim = 0.3, scale_sd = 0.5)
# generator B logit-scale noise SDs on positive patterns; scale_sd is
# the log-SD of a per-donor multiplier (donor heterogeneity)
GB_BASE <- "positive_margins" # generator B means come from this truth plus delta +-0.4

# Scenarios --------------------------------------------------------------------
# Full generator-A grid, then the generator-B scenario. in_pilot marks the 13
# scenarios that fit the 400-task cap at 30 seeds; in_smoke marks those run at 5 seeds.
scenarios_a <- expand_grid(
  truth = TRUTHS$truth,
  J = c(3L, 4L),
  n_donors = c(24L, 100L)
) |>
  left_join(TRUTHS, join_by(truth), relationship = "many-to-one") |>
  mutate(generator = "A")
# Generator B: M1-shaped means from GB_BASE with default-size interactions, then noise (in 03).
scenarios_b <- TRUTHS |>
  filter(truth == GB_BASE) |>
  mutate(
    truth = "generator_B",
    J = 3L,
    n_donors = 24L,
    delta_cmv = 0.4,
    delta_seb = -0.4,
    generator = "B"
  )

# Mixing check: data without donor x visit variation, fitted with donor_visit = TRUE (the
# cj_fit default); centered sampling may funnel as tau_u -> 0. Run at 5 seeds via in_dv_check.
scenarios_dv <- scenarios_a |>
  filter(truth == "positive_margins", J == 3L, n_donors == 24L) |>
  mutate(truth = "positive_margins_no_dv", generator = "A")

scenarios <- bind_rows(scenarios_a, scenarios_b, scenarios_dv) |>
  mutate(
    scenario_id = sprintf("%s_J%d_D%d", truth, J, n_donors),
    N = N_CELLS,
    missing_control = MISSING_CONTROL,
    tau_basal = TAU_BASAL,
    tau_u = TAU_U,
    omega_rho = OMEGA_RHO,
    theta0_start = THETA0_START,
    gb_sd_donor = if_else(generator == "B", GB_SD[["donor"]], NA_real_),
    gb_sd_sample = if_else(generator == "B", GB_SD[["sample"]], NA_real_),
    gb_sd_stim = if_else(generator == "B", GB_SD[["stim"]], NA_real_),
    gb_scale_sd = if_else(generator == "B", GB_SD[["scale_sd"]], NA_real_),
    ablation = truth == "visit_change", # also fit donor_visit = FALSE
    in_pilot = truth != "positive_margins_no_dv" &
      (J == 3L |
        (J == 4L &
          n_donors == 24L &
          truth %in% c("positive_margins", "dependence_only"))),
    in_smoke = in_pilot,
    # Centered sampling funnels when the true induced-effect SD is 0 (smoke: R-hat 1.1-1.4,
    # ESS < 30 for null and dependence_only), so those truths are fitted noncentered.
    centered = tau_induced > 0,
    sim_donor_visit = truth != "positive_margins_no_dv", # FALSE: no donor x visit variation in the data
    in_dv_check = truth == "positive_margins_no_dv"
  ) |>
  relocate(scenario_id, generator, truth, J, n_donors)

# Stan compile -----------------------------------------------------------------
# Compile the M1 executable once here, so the SLURM array tasks do not race to build it.
pkgload::load_all("cytoJointR", quiet = TRUE)
source("scripts/helper01_true-scores.R") # write_session_info()
cytoJointR:::m1_compile()

# Write outputs ----------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(
  "usage: Rscript 02_define-scenarios.R <run_dir>" = !is.na(run_dir),
  "run_dir already exists" = !dir.exists(run_dir)
)
dir.create(run_dir, recursive = TRUE)
write_tsv(scenarios, file.path(run_dir, "scenarios.tsv"))
write_session_info(file.path(run_dir, "session_info.txt"))
