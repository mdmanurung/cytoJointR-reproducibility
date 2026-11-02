#!/usr/bin/env Rscript
# Purpose: freeze the benchmark before any candidate is scored. Writes the seeded RV144 donor split, the reduced
#          simulation design (solved scenarios copied from the cached runs), and benchmark/definition.tsv: one row
#          per metric with its role, direction, the comparator values and the rule. Comparator values come from the
#          cached COMPASS and polyICSFlow outputs only, restricted to dev donors. No cytoJointR score is read.
#          RV144-test comparators are NOT computed here: they are computed once, in `70_run-scorecard.R --test`.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/69_define-scorecard.R
#          Refuses to run if benchmark/definition.tsv exists: the definition is frozen once committed.
# Inputs:
#   data/processed/rv144/donor_meta.tsv                         arm, infection outcome [20_prepare-rv144.R]
#   runs/2026-11-03_realdata/vaccine/rv144/donor_table.tsv       cached COMPASS and polyICSFlow RV144 scores, whole cohort [30_analyse-vaccine.R]
#   runs/2026-11-03_realdata/vaccine/rv144/infection.tsv         COMPASS PFS infection OR on the whole cohort [30_analyse-vaccine.R]
#   runs/2026-11-03_realdata/tb/{tb_s,tb_u}/donor_wide.tsv       cached COMPASS and polyICSFlow TB scores [31_analyse-tb.R]
#   runs/2026-11-02_redesign/scenarios.tsv, reps/               solved A/B/C scenarios; cached COMPASS and polyICSFlow replicates [05, 06, 07]
#   runs/2026-11-02_redesign-D2/scenarios.tsv, reps/            solved D1/D2 scenarios; cached replicates [10, 06, 07]
#   scripts/helper06_scorecard.R                                metric helpers [this repo]
# Outputs (in benchmark/):
#   rv144_split.tsv     donor_id, arm, infect, stratum, split (dev or test)
#   sim_scenarios.tsv   the solved rows of the simulation guards, same columns as 05's scenarios.tsv
#   definition.tsv      the frozen definition (columns described in benchmark/README.md)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
  library(tidyr)
})
source("scripts/helper06_scorecard.R")

# Parameters -------------------------------------------------------------------
SPLIT_SEED <- 20261104L # seeds the RV144 split; new, so it shares no stream with the model seed 20261103
BOOT_SEED <- 20261103L # seeds every donor bootstrap, as 30_analyse-vaccine.R
BOOT <- 2000L # donor bootstrap replicates, as 30_analyse-vaccine.R
SIM_SCENARIOS <- c(
  # one or two scenarios per generator; reduced from 18 + 8 so a scorecard run fits in about 1 h of SLURM
  "rare_weak_J3",
  "rare_weak_J4",
  "moderate_J3", # generator A
  "generator_B_J3", # generator B
  "mixture_J3", # generator C
  "D1_rare_weak_J3",
  "D2_rare_weak_J3" # held-out generator D (cell-level phenotype mixture, t3 donor effects)
)
SIM_SEEDS <- 1:10 # replicates per scenario, of the 30 cached ones; the cached COMPASS and polyICSFlow runs exist for them
DETECTION_SCENARIOS <- c("rare_weak_J3", "rare_weak_J4") # pooled for the simulation detection and magnitude metrics
SIM_RUNS <- c("runs/2026-11-02_redesign", "runs/2026-11-02_redesign-D2") # cached A/B/C and D runs
SCORES <- c(
  RMS_mean = "RMS posterior mean",
  P_RMS = "P(RMS > 1e-4)",
  CPS_mean = "CPS posterior mean",
  P_CPS = "P(CPS > 0.01)",
  CCS_mean = "CCS posterior mean (experimental)",
  P_CCS = "P(CCS > 1e-4) (experimental)"
)
LEVELS <- c("sample", "structured")
COVERAGE_MIN <- 0.85 # 11_apply-decision-rule.R rule 1
COVERAGE_MAX <- 0.97 # rule 1; not applied to null scenarios
COVERAGE_SE_MULT <- 2 # rule 1: pass if coverage + 2 SE >= min and coverage - 2 SE <= max
stopifnot(
  "benchmark/definition.tsv exists; the definition is frozen" = !file.exists(
    "benchmark/definition.tsv"
  )
)

# RV144 split ------------------------------------------------------------------
# Stratified 50/50 by arm and infection status (placebo infection status is not recorded, so placebo is one stratum).
meta <- read_tsv("data/processed/rv144/donor_meta.tsv", show_col_types = FALSE)
set.seed(SPLIT_SEED)
split <- meta |>
  transmute(
    donor_id,
    arm,
    infect,
    stratum = case_when(
      arm == "placebo" ~ "placebo",
      infect == "INFECTED" ~ "vaccine_infected",
      .default = "vaccine_noninfected"
    )
  ) |>
  arrange(donor_id) |>
  mutate(
    split = if_else(
      row_number() %in% sample(row_number(), n() %/% 2),
      "dev",
      "test"
    ),
    .by = stratum
  )
write_tsv(split, "benchmark/rv144_split.tsv")

# Simulation design ------------------------------------------------------------
sc <- map(SIM_RUNS, \(r) {
  read_tsv(file.path(r, "scenarios.tsv"), show_col_types = FALSE)
}) |>
  list_rbind() |>
  filter(scenario_id %in% SIM_SCENARIOS) |>
  mutate(n_seeds = length(SIM_SEEDS))
stopifnot(
  "a guard scenario is missing from the cached runs" = setequal(
    sc$scenario_id,
    SIM_SCENARIOS
  )
)
write_tsv(
  sc |> arrange(match(scenario_id, SIM_SCENARIOS)),
  "benchmark/sim_scenarios.tsv"
)

# Comparator values on dev donors ----------------------------------------------
# truth: TRUE = case. CIs are donor bootstraps within case and non-case.
comp_row <- \(compass, polyics, truth, ds) {
  cs <- boot_det(compass, truth, B = BOOT, seed = BOOT_SEED)
  ps <- boot_det(polyics, truth, B = BOOT, seed = BOOT_SEED)
  inner_join(
    cs,
    ps,
    join_by(stat),
    suffix = c("_c", "_p"),
    relationship = "one-to-one"
  ) |>
    mutate(dataset = ds)
}
rv <- read_tsv(
  "runs/2026-11-03_realdata/vaccine/rv144/donor_table.tsv",
  show_col_types = FALSE
) |>
  inner_join(
    split |> select(donor_id, split),
    join_by(donor_id),
    relationship = "one-to-one"
  ) |>
  filter(split == "dev")
tbs <- map(c(tb_s = "tb_s", tb_u = "tb_u"), \(d) {
  read_tsv(
    file.path("runs/2026-11-03_realdata/tb", d, "donor_wide.tsv"),
    show_col_types = FALSE
  )
})
comp <- bind_rows(
  comp_row(
    rv$COMPASS_FS,
    rv$polyICSFlow_clipped_net_sum,
    rv$vaccinated,
    "rv144_dev"
  ),
  imap(tbs, \(t, d) {
    comp_row(t[["COMPASS FS"]], t[["polyICSFlow_clipped net_sum"]], t$case, d)
  }) |>
    list_rbind()
)
cp_ci <- \(k, n) unname(binom.test(k, n)$conf.int)
rule_row <- \(called, truth, kind) {
  k <- if (kind == "fpr") sum(called[!truth]) else sum(called[truth])
  n <- if (kind == "fpr") sum(!truth) else sum(truth)
  tibble(value = k / n, lower = cp_ci(k, n)[1], upper = cp_ci(k, n)[2])
}

# Simulation comparators: pooled over the detection scenarios and seeds --------
read_rep <- \(run, scen, seed, f, sub = "") {
  p <- file.path(run, "reps", scen, sprintf("seed_%02d", seed), sub, f)
  read_tsv(p, show_col_types = FALSE) |> mutate(scenario_id = scen, seed = seed)
}
sim_run <- SIM_RUNS[1]
reps <- expand_grid(scen = DETECTION_SCENARIOS, seed = SIM_SEEDS)
sim_truth <- pmap(reps, \(scen, seed) {
  read_rep(sim_run, scen, seed, "truth_donor.tsv")
}) |>
  list_rbind()
sim_det <- pmap(reps, \(scen, seed) {
  read_rep(sim_run, scen, seed, "donor_scores.tsv") |>
    filter(method == "COMPASS_default", score == "FS") |>
    select(scenario_id, seed, donor_id, compass = value)
}) |>
  list_rbind()
sim_poly <- pmap(reps, \(scen, seed) {
  read_rep(sim_run, scen, seed, "donor_scores.tsv", "polyicsflow") |>
    filter(method == "polyICSFlow_clipped") |>
    select(scenario_id, seed, donor_id, polyics = value)
}) |>
  list_rbind()
sim_mag <- pmap(reps, \(scen, seed) {
  read_rep(sim_run, scen, seed, "donor_magnitude.tsv", "polyicsflow") |>
    filter(method == "polyICSFlow_unclipped") |>
    select(scenario_id, seed, donor_id, polyics_mag = estimate)
}) |>
  list_rbind()
sim_tab <- sim_truth |>
  left_join(
    sim_det,
    join_by(scenario_id, seed, donor_id),
    relationship = "one-to-one"
  ) |>
  left_join(
    sim_poly,
    join_by(scenario_id, seed, donor_id),
    relationship = "one-to-one"
  ) |>
  left_join(
    sim_mag,
    join_by(scenario_id, seed, donor_id),
    relationship = "one-to-one"
  ) |>
  mutate(compass = coalesce(compass, 0)) # a COMPASS fit that failed reports no response, as 08
sim_auc_c <- seed_boot(
  sim_tab |>
    transmute(
      seed = paste(scenario_id, seed),
      truth = responder_sample,
      a = compass
    ),
  "auc",
  B = BOOT,
  seed = BOOT_SEED
)
sim_auc_p <- seed_boot(
  sim_tab |>
    transmute(
      seed = paste(scenario_id, seed),
      truth = responder_sample,
      a = polyics
    ),
  "auc",
  B = BOOT,
  seed = BOOT_SEED
)
sim_sp_c <- seed_boot(
  sim_tab |>
    transmute(
      seed = paste(scenario_id, seed),
      truth = true_rms_sample,
      a = compass
    ),
  "spearman",
  B = BOOT,
  seed = BOOT_SEED
)
sim_sp_p <- seed_boot(
  sim_tab |>
    transmute(
      seed = paste(scenario_id, seed),
      truth = true_rms_sample,
      a = polyics_mag
    ),
  "spearman",
  B = BOOT,
  seed = BOOT_SEED
)

# Infection reference ----------------------------------------------------------
inf_ref <- read_tsv(
  "runs/2026-11-03_realdata/vaccine/rv144/infection.tsv",
  show_col_types = FALSE
) |>
  filter(label == "COMPASS_PFS", adj == "M1")

# Definition rows --------------------------------------------------------------
cv <- \(ds, st, which) {
  comp |>
    filter(dataset == ds, stat == st) |>
    select(
      value = paste0("value_", which),
      lower = paste0("lower_", which),
      upper = paste0("upper_", which)
    )
}
def_row <- \(
  metric_id,
  dataset,
  role,
  tier,
  metric,
  score,
  level,
  direction,
  cmp = NULL,
  pol = NULL,
  rule = NA_character_,
  note = NA_character_
) {
  tibble(
    metric_id,
    dataset,
    role,
    tier,
    metric,
    score,
    level,
    direction,
    compass_value = cmp$value %||% NA_real_,
    compass_lower = cmp$lower %||% NA_real_,
    compass_upper = cmp$upper %||% NA_real_,
    polyics_value = pol$value %||% NA_real_,
    polyics_lower = pol$lower %||% NA_real_,
    polyics_upper = pol$upper %||% NA_real_,
    rule,
    note
  )
}
det_stats_names <- c(
  auc = "AUC",
  pauc10 = "partial AUC, FPR <= 0.1 (mean TPR over that range)",
  tpr_fpr05 = "TPR at FPR <= 0.05",
  tpr_fpr10 = "TPR at FPR <= 0.10"
)
CMP_NOTE <- "comparators: COMPASS FS, polyICSFlow clipped net sum, cached, dev donors only"

# RV144 dev: primary and its companions
rows_rv_primary <- imap(det_stats_names, \(lab, st) {
  def_row(
    paste0("rv144_dev.", st, ".P_RMS.fit"),
    "rv144_dev",
    "dev",
    if_else(st == "auc", "primary", "companion"),
    lab,
    "P_RMS",
    "fit",
    "higher",
    cv("rv144_dev", st, "c"),
    cv("rv144_dev", st, "p"),
    rule = "detection score P(RMS > 1e-4) of the candidate at its own fit level; vaccinees vs placebos",
    note = CMP_NOTE
  )
}) |>
  list_rbind()
# Secondary: every score family at both levels, RV144 dev and the TB sets
rows_auc <- expand_grid(
  dataset = c("rv144_dev", "tb_s", "tb_u"),
  score = names(SCORES),
  level = LEVELS
) |>
  mutate(
    r = pmap(list(dataset, score, level), \(d, s, l) {
      def_row(
        sprintf("%s.auc.%s.%s", d, s, l),
        d,
        "dev",
        "secondary",
        paste("AUC of", SCORES[[s]]),
        s,
        l,
        "higher",
        cv(d, "auc", "c"),
        cv(d, "auc", "p"),
        rule = if (d == "rv144_dev") {
          "vaccinees vs placebos"
        } else {
          "QFT+ vs QFT-"
        },
        note = paste(
          CMP_NOTE,
          if (str_detect(s, "CCS")) "| CCS is experimental" else ""
        )
      )
    })
  ) |>
  pull(r) |>
  list_rbind()
# TB primary-level AUC at the candidate's fit level (the headline TB numbers)
rows_tb_fit <- expand_grid(dataset = c("tb_s", "tb_u")) |>
  mutate(
    r = map(dataset, \(d) {
      def_row(
        paste0(d, ".auc.P_RMS.fit"),
        d,
        "dev",
        "secondary",
        "AUC of P(RMS > 1e-4)",
        "P_RMS",
        "fit",
        "higher",
        cv(d, "auc", "c"),
        cv(d, "auc", "p"),
        rule = "QFT+ (22) vs QFT- (18)",
        note = CMP_NOTE
      )
    })
  ) |>
  pull(r) |>
  list_rbind()
# Fixed rule on RV144 dev
fr <- \(called, kind) rule_row(called, rv$vaccinated, kind)
rows_rule <- bind_rows(
  def_row(
    "rv144_dev.fpr_rule.P_RMS.fit",
    "rv144_dev",
    "dev",
    "secondary",
    "placebo false-positive rate at the fixed rule",
    "P_RMS",
    "fit",
    "lower",
    fr(rv$COMPASS_Response > 0.5, "fpr"),
    fr(rv$polyICSFlow_clipped_net_sum > 0, "fpr"),
    rule = "candidate: P(RMS > 1e-4) > 0.95; COMPASS: Response (degree 1) > 0.5; polyICSFlow: clipped net > 0",
    note = "fixed rule epsilon = 1e-4, P > 0.95; nobody changes it during iteration"
  ),
  def_row(
    "rv144_dev.tpr_rule.P_RMS.fit",
    "rv144_dev",
    "dev",
    "secondary",
    "vaccinee true-positive rate at the fixed rule",
    "P_RMS",
    "fit",
    "higher",
    fr(rv$COMPASS_Response > 0.5, "tpr"),
    fr(rv$polyICSFlow_clipped_net_sum > 0, "tpr"),
    rule = "as rv144_dev.fpr_rule",
    note = "COMPASS and the candidate sit at different ROC corners; read with pauc10 and tpr_fpr10"
  )
)
# Test-only rows: comparators are computed at test time
rows_test <- bind_rows(
  imap(det_stats_names, \(lab, st) {
    def_row(
      paste0("rv144_test.", st, ".P_RMS.fit"),
      "rv144_test",
      "test",
      "test_only",
      lab,
      "P_RMS",
      "fit",
      "higher",
      rule = "as the dev row, on test donors",
      note = "comparators computed once in --test mode"
    )
  }) |>
    list_rbind(),
  def_row(
    "rv144_test.fpr_rule.P_RMS.fit",
    "rv144_test",
    "test",
    "test_only",
    "placebo false-positive rate at the fixed rule",
    "P_RMS",
    "fit",
    "lower",
    rule = "as rv144_dev.fpr_rule"
  ),
  def_row(
    "rv144_test.tpr_rule.P_RMS.fit",
    "rv144_test",
    "test",
    "test_only",
    "vaccinee true-positive rate at the fixed rule",
    "P_RMS",
    "fit",
    "higher",
    rule = "as rv144_dev.fpr_rule"
  ),
  map(c("RMS_mean", "CPS_mean"), \(s) {
    def_row(
      paste0("rv144_test.infection_or.", s, ".fit"),
      "rv144_test",
      "test",
      "test_only",
      paste(
        "infection OR per SD of",
        SCORES[[s]],
        "among test-half vaccinees, adjusted M1 (gender, age, risk)"
      ),
      s,
      "fit",
      "lower",
      cmp = tibble(
        value = inf_ref$or_per_sd,
        lower = inf_ref$lower,
        upper = inf_ref$upper
      ),
      rule = "logistic regression as 30_analyse-vaccine.R M1; a candidate is protective-consistent if OR < 1 and p < 0.05",
      note = sprintf(
        "comparator shown is COMPASS PFS on the WHOLE cohort (OR %.2f, p %.3f, 226 vaccinees, 38 events); the test-half COMPASS PFS value is computed in --test mode",
        inf_ref$or_per_sd,
        inf_ref$p
      )
    )
  }) |>
    list_rbind()
)
# Simulation guards
rows_cov <- expand_grid(scenario_id = SIM_SCENARIOS, score = c("RMS", "CPS")) |>
  mutate(
    r = map2(scenario_id, score, \(s, sc_) {
      def_row(
        sprintf("sim.cov90.%s.%s", sc_, s),
        s,
        "guard",
        "guard",
        sprintf("90%% interval coverage of donor %s", sc_),
        sc_,
        "fit",
        "pass",
        rule = sprintf(
          "pass if coverage + %g SE >= %g and coverage - %g SE <= %g (no upper bound in null scenarios); SE clustered by seed; truth matches the fit level (sample or structured)",
          COVERAGE_SE_MULT,
          COVERAGE_MIN,
          COVERAGE_SE_MULT,
          COVERAGE_MAX
        ),
        note = sprintf(
          "generator %s, %d seeds x %d donors",
          sc$generator[sc$scenario_id == s],
          length(SIM_SEEDS),
          sc$n_donors[sc$scenario_id == s]
        )
      )
    })
  ) |>
  pull(r) |>
  list_rbind()
rows_sim <- bind_rows(
  def_row(
    "sim.auc.rare_weak.P_RMS.fit",
    "sim_rare_weak",
    "guard",
    "guard",
    "AUC of P(RMS > 1e-4), responders vs non-responders, pooled over rare_weak J3 + J4 and seeds",
    "P_RMS",
    "fit",
    "higher",
    tibble(
      value = sim_auc_c$value,
      lower = sim_auc_c$lower,
      upper = sim_auc_c$upper
    ),
    tibble(
      value = sim_auc_p$value,
      lower = sim_auc_p$lower,
      upper = sim_auc_p$upper
    ),
    rule = "responder = true sample-level RMS > 1e-4; CI resamples seeds",
    note = "comparators: COMPASS_default FS, polyICSFlow clipped net sum; cached seeds 1-10"
  ),
  def_row(
    "sim.spearman.rare_weak.RMS_mean.fit",
    "sim_rare_weak",
    "guard",
    "guard",
    "Spearman of the donor RMS estimate with the true RMS, pooled over rare_weak J3 + J4 and seeds",
    "RMS_mean",
    "fit",
    "higher",
    tibble(
      value = sim_sp_c$value,
      lower = sim_sp_c$lower,
      upper = sim_sp_c$upper
    ),
    tibble(
      value = sim_sp_p$value,
      lower = sim_sp_p$lower,
      upper = sim_sp_p$upper
    ),
    rule = "magnitude ranking; truth is the sample-level RMS",
    note = "comparators: COMPASS_default FS (a magnitude proxy), polyICSFlow UNCLIPPED net sum (the fair magnitude comparator)"
  )
)
# Win rows, computed by 70 from the rows above
rows_win <- bind_rows(
  def_row(
    "win.beats_compass",
    "all",
    "win",
    "win",
    "candidate beats COMPASS on the dev set",
    NA,
    "fit",
    "pass",
    rule = paste(
      "ALL of: (1) rv144_dev.auc.P_RMS.fit paired-bootstrap difference vs COMPASS FS has a 95% lower bound > 0;",
      "(2) rv144_dev.pauc10.P_RMS.fit difference > 0 and rv144_dev.tpr_fpr10.P_RMS.fit difference >= 0 (point estimates);",
      "(3) sim.auc.rare_weak.P_RMS.fit difference vs COMPASS_default FS >= 0 (point estimate);",
      "(4) every sim.cov90 guard passes. The claim stands only after (5): the same (1) on RV144 test, run once with --test."
    ),
    note = "COMPASS FS AUC is about 0.98, so (1) needs a candidate AUC near 0.99 on dev"
  ),
  def_row(
    "win.beats_polyics",
    "all",
    "win",
    "win",
    "candidate beats polyICSFlow on the dev set",
    NA,
    "fit",
    "pass",
    rule = paste(
      "ALL of: (1) rv144_dev.auc.P_RMS.fit paired-bootstrap difference vs polyICSFlow clipped has a 95% lower bound > 0",
      "(clipped detection is trivially beatable: it calls 36/36 placebos);",
      "(2) sim.spearman.rare_weak.RMS_mean.fit difference vs polyICSFlow unclipped has a 95% lower bound > 0 (the meaningful comparison is magnitude ranking);",
      "(3) sim.auc.rare_weak.P_RMS.fit difference vs polyICSFlow clipped >= 0 (point estimate; on the simulation clipped detection is not trivial);",
      "(4) every sim.cov90 guard passes. (5) the same (1) on RV144 test, run once with --test."
    ),
    note = "polyICSFlow has no uncertainty or calibration claim to test"
  )
)
# Datasets that are not scored
rows_na <- bind_rows(
  def_row(
    "hvtn078.excluded",
    "hvtn078",
    "excluded",
    "info",
    "not scored",
    NA,
    NA,
    NA,
    note = "5 placebos and every method scores about 0.4 AUC: no detection metric. Magnitude agreement may be reported by hand, not scored"
  ),
  def_row(
    "sarscov2.unusable_for_fpr",
    "sarscov2",
    "excluded",
    "info",
    "not scored",
    NA,
    NA,
    NA,
    rule = "E01 Spike cannot serve as pseudo-negatives",
    note = "Ford et al., Nat Immunol 2024, PMID 38057617, doi:10.1038/s41590-023-01692-x: the title says 'persons with previous COVID-19' and the Zenodo record 8088178 says 'persons recovered from SARS-CoV-2 infection'. Participants were not SARS-CoV-2-naive at enrolment, so E01 responses are expected. Source read: titles only; the full text could not be fetched (PMC and Nature blocked automated access)"
  )
)
definition <- bind_rows(
  rows_rv_primary,
  rows_auc,
  rows_tb_fit,
  rows_rule,
  rows_test,
  rows_cov,
  rows_sim,
  rows_win,
  rows_na
) |>
  distinct(metric_id, .keep_all = TRUE)
write_tsv(definition, "benchmark/definition.tsv")
