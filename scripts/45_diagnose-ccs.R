#!/usr/bin/env Rscript
# Purpose: diagnose why CCS fails calibration under generator B. Four hypotheses: (H1) the truth is the realised
#          sample probability, whose pattern noise no count-based estimate recovers; (H2) IPF amplifies errors in
#          rare control patterns; (H3) the C0 pseudo-count (and model shrinkage) sets the CCS of cells with no
#          control counts; (H4) a bug in ipf_draws() or in the truth. Uses regenerated simulations, the stored
#          benchmark replicates and C0 (no Stan). Diagnosis only: the package is not touched.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/45_diagnose-ccs.R runs/2026-11-03_diagnose-ccs [n_seeds]
# Inputs:
#   cytoJointR-wt/redesign-run/                 package source, frozen at 49278c9 [this repo]
#   runs/2026-11-02_redesign/scenarios.tsv      scenario rows [05_define-benchmark.R]
#   runs/2026-11-02_redesign/reps/<scenario>/seed_NN/{counts,model_scores}.tsv   stored replicates [06_run-benchmark-replicate.R]
# Outputs (in <run_dir>):
#   regeneration_check.tsv      does bench_simulate() reproduce the stored counts
#   h1_truth_variants.tsv       per contrast: CCS of realised/structured control x stimulated probabilities
#   h1_truth_spread.tsv         SD of each truth variant and of realised minus structured
#   h1_model_error.tsv          model estimate vs each truth variant: bias, rmse, coverage, correlation
#   control_information.tsv     expected and observed polyfunctional counts in control and stimulated samples
#   h3_a0_sweep.tsv             C0 CCS against the Dirichlet pseudo-count a0
#   h2_sensitivity.tsv          CCS SD from control-only and stimulated-only draws, IPF gain, elasticity
#   h4_implementation.tsv       ipf_draws vs cj_score_reference, independent IPF vs stored truth, row sums
#   pooled_reference.tsv        per-donor vs pooled-control reference: coverage against each truth
#   session_info.txt            R and package versions

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- "cytoJointR-wt/redesign-run" # frozen redesign package
REDESIGN_DIR <- "runs/2026-11-02_redesign" # stored benchmark replicates
SCENARIOS <- c("generator_B_J3", "moderate_J3") # B is the failing cell; moderate (generator A) is the zero-truth control
N_SEEDS <- 30 # as in the benchmark; override with the 2nd argument for a quick test
N_SEEDS_SMALL <- 3 # seeds for the sweeps that fit C0 several times
A0_GRID <- c(0.5, 0.05, 0.005) # Dirichlet pseudo-counts; 0.5 is the benchmark C0
N_DRAWS <- 1000 # posterior draws per sample
EPS_CCS <- 1e-4 # CCS epsilon used by the benchmark
LEVEL <- 0.9 # interval level used by the benchmark
SIGN_PROB <- 0.95 # false-sign threshold used by the benchmark
SENS_STEP <- 0.01 # log-scale tilt for the finite-difference elasticity

# Setup ------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[1]
arg_seeds <- as.integer(commandArgs(trailingOnly = TRUE)[2])
if (!is.na(arg_seeds)) {
  N_SEEDS <- arg_seeds
}
stopifnot(
  "usage: Rscript 45_diagnose-ccs.R <run_dir> [n_seeds]" = !is.na(run_dir),
  "run_dir already exists" = !dir.exists(run_dir)
)
dir.create(run_dir, recursive = TRUE)
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper01_true-scores.R") # true_table(), write_session_info()
source("scripts/helper02_benchmark-sim.R") # bench_simulate()
scenarios <- read_tsv(
  file.path(REDESIGN_DIR, "scenarios.tsv"),
  show_col_types = FALSE
)
EPSILON <- c(
  RMS = 1e-4,
  poly_absolute = 1e-4,
  poly_marginal = 1e-4,
  CPS = 0.01,
  CCS = EPS_CCS,
  breadth = 0.01,
  marker_magnitude = 1e-4
)
score_contrast <- cytoJointR:::score_contrast
ipf_draws <- cytoJointR:::ipf_draws

# Helpers ----------------------------------------------------------------------
# Regenerate the simulation of one replicate exactly as 06_run-benchmark-replicate.R does.
sim_of <- function(scenario_id, seed) {
  sc <- filter(scenarios, scenario_id == !!scenario_id)
  bench_simulate(sc, sc$seed_offset + seed)
}
# Contrast table with control, and the polyfunctional pattern mask.
contrasts_of <- function(sim) {
  filter(sim$data$contrasts, !is.na(control_sample_id))
}
is_poly <- function(sim) rowSums(sim$data$bits) >= 2
# Global CCS of one control / stimulated probability pair (the package oracle).
ccs_of <- function(pu, ps, bits) {
  cj_score_reference(pu, ps, bits)$global[["CCS"]]
}
# n draws of a Dirichlet(y + a0); a tiny shape can underflow to 0, so keep full support as cj_fit_dirichlet does.
draw_dirichlet <- function(y, a0, n) {
  g <- pmax(
    matrix(rgamma(n * length(y), shape = rep(y + a0, each = n)), n),
    1e-300
  )
  g / rowSums(g)
}
# Summary of CCS draws (failed IPF draws are NA and counted).
ccs_summary <- function(draws) {
  ok <- draws[!is.na(draws)]
  tibble(
    estimate = mean(ok),
    lower = quantile(ok, (1 - LEVEL) / 2, names = FALSE),
    upper = quantile(ok, (1 + LEVEL) / 2, names = FALSE),
    p_positive = mean(ok > EPS_CCS),
    p_negative = mean(ok < -EPS_CCS),
    sd = sd(ok),
    failure_fraction = mean(is.na(draws))
  )
}
sims <- expand_grid(scenario_id = SCENARIOS, seed = seq_len(N_SEEDS)) |>
  mutate(sim = map2(scenario_id, seed, sim_of))

# Regeneration check -----------------------------------------------------------
# The stored counts must equal the regenerated ones, or every truth below is for a different data set.
regeneration_check <- sims |>
  mutate(
    identical_counts = pmap_lgl(
      list(scenario_id, seed, sim),
      \(scenario_id, seed, sim) {
        stored <- read_tsv(
          file.path(
            REDESIGN_DIR,
            "reps",
            scenario_id,
            sprintf("seed_%02d", seed),
            "counts.tsv"
          ),
          show_col_types = FALSE
        )
        regenerated <- sim$data$counts
        colnames(regenerated) <- apply(sim$data$bits, 1, paste, collapse = "")
        compared <- regenerated |>
          as_tibble(rownames = "sample_id") |>
          pivot_longer(
            -sample_id,
            names_to = "pattern_key",
            values_to = "n_new"
          ) |>
          inner_join(
            stored,
            join_by(sample_id, pattern_key),
            relationship = "one-to-one"
          )
        nrow(compared) == nrow(stored) && all(compared$n == compared$n_new)
      }
    )
  ) |>
  select(scenario_id, seed, identical_counts)
stopifnot(
  "regenerated counts differ from the stored replicates" = all(
    regeneration_check$identical_counts
  )
)

# H1: truth variants -----------------------------------------------------------
# CCS of (realised or structured control) x (realised or structured stimulated) probabilities.
# real_real is the benchmark truth_sample, struct_struct is truth_struct; the mixed ones show which side
# of the contrast carries the pattern noise.
h1_truth_variants <- sims |>
  mutate(
    truths = map(sim, \(sim) {
      ct <- contrasts_of(sim)
      bits <- sim$data$bits
      ct |>
        mutate(
          ccs_real_real = map2_dbl(control_sample_id, stim_sample_id, \(u, s) {
            ccs_of(sim$p[u, ], sim$p[s, ], bits)
          }),
          ccs_struct_struct = map2_dbl(
            control_sample_id,
            stim_sample_id,
            \(u, s) ccs_of(sim$p_struct[u, ], sim$p_struct[s, ], bits)
          ),
          ccs_ustruct_sreal = map2_dbl(
            control_sample_id,
            stim_sample_id,
            \(u, s) ccs_of(sim$p_struct[u, ], sim$p[s, ], bits)
          ),
          ccs_ureal_sstruct = map2_dbl(
            control_sample_id,
            stim_sample_id,
            \(u, s) ccs_of(sim$p[u, ], sim$p_struct[s, ], bits)
          )
        ) |>
        select(contrast_id, starts_with("ccs_"))
    })
  ) |>
  select(scenario_id, seed, truths) |>
  unnest(truths)
h1_truth_spread <- h1_truth_variants |>
  mutate(real_minus_struct = ccs_real_real - ccs_struct_struct) |>
  pivot_longer(
    c(starts_with("ccs_"), real_minus_struct),
    names_to = "variant"
  ) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    max_abs = max(abs(value)),
    .by = c(scenario_id, variant)
  )

# H1: model error against each truth variant -----------------------------------
# Stored benchmark estimates (C0, C1, C2; C2 on sample and structured draws) scored against every truth variant.
stored_ccs <- sims |>
  distinct(scenario_id, seed) |>
  mutate(
    scores = map2(scenario_id, seed, \(scenario_id, seed) {
      read_tsv(
        file.path(
          REDESIGN_DIR,
          "reps",
          scenario_id,
          sprintf("seed_%02d", seed),
          "model_scores.tsv"
        ),
        col_select = c(
          model,
          draws,
          contrast_id,
          grain,
          score,
          estimate,
          lower,
          upper,
          failure_fraction
        ),
        show_col_types = FALSE
      ) |>
        filter(grain == "global", score == "CCS")
    })
  ) |>
  unnest(scores)
h1_model_error <- stored_ccs |>
  left_join(
    h1_truth_variants,
    join_by(scenario_id, seed, contrast_id),
    relationship = "many-to-one"
  ) |>
  pivot_longer(
    starts_with("ccs_"),
    names_to = "truth_variant",
    values_to = "truth"
  ) |>
  summarise(
    bias = mean(estimate - truth),
    rmse = sqrt(mean((estimate - truth)^2)),
    coverage = mean(truth >= lower & truth <= upper),
    correlation = suppressWarnings(cor(estimate, truth)),
    mean_posterior_sd = mean((upper - lower) / (2 * qnorm((1 + LEVEL) / 2))),
    mean_failure_fraction = mean(failure_fraction),
    .by = c(scenario_id, model, draws, truth_variant)
  )

# Control information ----------------------------------------------------------
# How many polyfunctional cells does a control sample hold? Counts near 0 mean the baseline dependence is not observed.
control_information <- sims |>
  mutate(
    info = map(sim, \(sim) {
      ct <- contrasts_of(sim)
      poly <- is_poly(sim)
      n_cells <- sim$data$samples$N[match(
        ct$control_sample_id,
        sim$data$samples$sample_id
      )]
      tibble(
        contrast_id = ct$contrast_id,
        expected_control_poly = n_cells *
          rowSums(sim$p_struct[ct$control_sample_id, poly, drop = FALSE]),
        observed_control_poly = rowSums(sim$data$counts[
          ct$control_sample_id,
          poly,
          drop = FALSE
        ]),
        expected_stim_poly = n_cells *
          rowSums(sim$p_struct[ct$stim_sample_id, poly, drop = FALSE]),
        observed_stim_poly = rowSums(sim$data$counts[
          ct$stim_sample_id,
          poly,
          drop = FALSE
        ])
      )
    })
  ) |>
  select(scenario_id, seed, info) |>
  unnest(info) |>
  summarise(
    mean_expected_control_poly = mean(expected_control_poly),
    mean_observed_control_poly = mean(observed_control_poly),
    share_controls_zero_poly = mean(observed_control_poly == 0),
    mean_expected_stim_poly = mean(expected_stim_poly),
    mean_observed_stim_poly = mean(observed_stim_poly),
    mean_pooled_control_poly_per_seed = sum(observed_control_poly) /
      n_distinct(seed),
    .by = scenario_id
  )

# H3: pseudo-count sweep -------------------------------------------------------
# C0 with several a0. If CCS tracks a0 (negative at a0 = 0.5, positive as a0 -> 0), the prior sets it.
h3_a0_sweep <- expand_grid(
  sims |> filter(seed <= N_SEEDS_SMALL),
  a0 = A0_GRID
) |>
  mutate(
    table = pmap(list(sim, a0, seed), \(sim, a0, seed) {
      x <- cj_fit_dirichlet(sim$data, a0 = a0, n_draws = N_DRAWS, seed = seed)
      truth <- true_table(sim) |>
        filter(grain == "global", score %in% c("CCS", "RMS")) |>
        select(contrast_id, score, truth_sample, truth_struct)
      cj_scores(x, EPSILON, level = LEVEL) |>
        filter(score %in% c("CCS", "RMS"), is.na(marker)) |>
        select(
          contrast_id,
          score,
          estimate,
          lower,
          upper,
          p_positive,
          p_negative,
          failure_fraction
        ) |>
        left_join(
          truth,
          join_by(contrast_id, score),
          relationship = "one-to-one"
        )
    })
  ) |>
  select(scenario_id, seed, a0, table) |>
  unnest(table) |>
  pivot_wider(
    id_cols = c(scenario_id, seed, a0, contrast_id),
    names_from = score,
    values_from = c(
      estimate,
      lower,
      upper,
      p_positive,
      p_negative,
      failure_fraction,
      truth_sample,
      truth_struct
    )
  ) |>
  # a contrast whose IPF failed in more than 1% of draws is withheld by cj_scores() (NA); count it, average the rest
  summarise(
    share_withheld = mean(is.na(estimate_CCS)),
    mean_estimate_ccs = mean(estimate_CCS, na.rm = TRUE),
    mean_truth_sample = mean(truth_sample_CCS),
    mean_truth_struct = mean(truth_struct_CCS),
    bias_vs_sample = mean(estimate_CCS - truth_sample_CCS, na.rm = TRUE),
    bias_vs_struct = mean(estimate_CCS - truth_struct_CCS, na.rm = TRUE),
    coverage_sample = mean(
      truth_sample_CCS >= lower_CCS & truth_sample_CCS <= upper_CCS,
      na.rm = TRUE
    ),
    coverage_struct = mean(
      truth_struct_CCS >= lower_CCS & truth_struct_CCS <= upper_CCS,
      na.rm = TRUE
    ),
    share_sign_claimed = mean(
      p_positive_CCS > SIGN_PROB | p_negative_CCS > SIGN_PROB,
      na.rm = TRUE
    ),
    share_negative_claimed = mean(p_negative_CCS > SIGN_PROB, na.rm = TRUE),
    share_positive_claimed = mean(p_positive_CCS > SIGN_PROB, na.rm = TRUE),
    mean_rms = mean(estimate_RMS, na.rm = TRUE),
    mean_failure_fraction = mean(failure_fraction_CCS, na.rm = TRUE),
    .by = c(scenario_id, a0)
  )

# H2: IPF sensitivity ----------------------------------------------------------
# At each contrast (C0 a0 = 0.5): CCS SD from control-only draws, stimulated-only draws and both; the IPF gain
# p0_k / p_U_k of polyfunctional cells; the elasticity d CCS / d log p_U_k (finite difference) by pattern degree.
h2_sensitivity <- sims |>
  filter(seed <= N_SEEDS_SMALL) |>
  mutate(
    table = map2(sim, seed, \(sim, seed) {
      set.seed(seed)
      bits <- sim$data$bits
      deg <- rowSums(bits)
      poly <- deg >= 2
      ct <- contrasts_of(sim)
      pmap(
        list(ct$contrast_id, ct$control_sample_id, ct$stim_sample_id),
        \(id, u, s) {
          yu <- sim$data$counts[u, ]
          ys <- sim$data$counts[s, ]
          pu <- (yu + 0.5) / sum(yu + 0.5)
          ps <- (ys + 0.5) / sum(ys + 0.5)
          rep_row <- \(p) matrix(p, N_DRAWS, length(p), byrow = TRUE)
          ccs_u_only <- score_contrast(
            draw_dirichlet(yu, 0.5, N_DRAWS),
            rep_row(ps),
            bits
          )$G[, "CCS"]
          ccs_s_only <- score_contrast(
            rep_row(pu),
            draw_dirichlet(ys, 0.5, N_DRAWS),
            bits
          )$G[, "CCS"]
          ccs_both <- score_contrast(
            draw_dirichlet(yu, 0.5, N_DRAWS),
            draw_dirichlet(ys, 0.5, N_DRAWS),
            bits
          )$G[, "CCS"]
          p0 <- ipf_draws(rbind(pu), rbind(ps), bits)$Q[1, ]
          base <- sum((ps - p0)[poly])
          elasticity <- map_dbl(seq_along(pu), \(k) {
            tilted <- pu
            tilted[k] <- tilted[k] * exp(SENS_STEP)
            tilted <- tilted / sum(tilted)
            (sum((ps - ipf_draws(rbind(tilted), rbind(ps), bits)$Q[1, ])[
              poly
            ]) -
              base) /
              SENS_STEP
          })
          # posterior SD of log p_U_k is about sqrt(trigamma(y + a0)); first-order CCS SD from the control side
          sd_log_pu <- sqrt(trigamma(yu + 0.5))
          tibble(
            contrast_id = id,
            ccs_point = base,
            sd_ccs_control_only = sd(ccs_u_only, na.rm = TRUE),
            sd_ccs_stim_only = sd(ccs_s_only, na.rm = TRUE),
            sd_ccs_both = sd(ccs_both, na.rm = TRUE),
            sd_ccs_first_order = sqrt(sum((elasticity * sd_log_pu)^2)),
            gain_degree2 = mean((p0 / pu)[deg == 2]),
            gain_degree3 = if (any(deg == 3)) {
              mean((p0 / pu)[deg == 3])
            } else {
              NA_real_
            },
            p0_poly = sum(p0[poly]),
            ps_poly = sum(ps[poly]),
            pu_poly = sum(pu[poly]),
            abs_elasticity_degree2 = mean(abs(elasticity[deg == 2])),
            abs_elasticity_degree3 = if (any(deg == 3)) {
              mean(abs(elasticity[deg == 3]))
            } else {
              NA_real_
            },
            abs_elasticity_degree0 = abs(elasticity[deg == 0])
          )
        }
      ) |>
        list_rbind()
    })
  ) |>
  select(scenario_id, seed, table) |>
  unnest(table) |>
  left_join(
    h1_truth_variants,
    join_by(scenario_id, seed, contrast_id),
    relationship = "one-to-one"
  ) |>
  mutate(gap_point_vs_sample_truth = abs(ccs_point - ccs_real_real))

# H4: implementation -----------------------------------------------------------
# (a) ipf_draws-based score_contrast() equals the per-draw oracle cj_score_reference().
# (b) an independent IPF (fixed 500 sweeps, base R) reproduces the stored sample truth.
# (c) pattern CCS truths sum to the global truth; stored C0/C1 pattern estimates sum to the global estimate.
h4_scenario <- SCENARIOS[1]
h4_sim <- sims |>
  filter(scenario_id == h4_scenario, seed == 1) |>
  pull(sim) |>
  pluck(1)
h4_bits <- h4_sim$data$bits
h4_ct <- contrasts_of(h4_sim)
h4_draws <- cj_fit_dirichlet(h4_sim$data, a0 = 0.5, n_draws = 200, seed = 1)
h4_engine <- h4_ct |>
  slice_head(n = 8) |>
  mutate(
    max_abs_diff = map2_dbl(control_sample_id, stim_sample_id, \(u, s) {
      Pu <- matrix(h4_draws$p[, u, ], 200)
      Ps <- matrix(h4_draws$p[, s, ], 200)
      G <- score_contrast(Pu, Ps, h4_bits)$G[, 1:6]
      oracle <- t(vapply(
        seq_len(200),
        \(i) cj_score_reference(Pu[i, ], Ps[i, ], h4_bits)$global,
        numeric(6)
      ))
      max(abs(G - oracle), na.rm = TRUE)
    })
  )
ipf_independent <- function(pu, ps, bits) {
  target <- colSums(bits * ps)
  q <- pu
  for (sweep in 1:500) {
    for (j in seq_len(ncol(bits))) {
      a <- sum(q[bits[, j] == 1])
      q[bits[, j] == 1] <- q[bits[, j] == 1] * target[j] / a
      q[bits[, j] == 0] <- q[bits[, j] == 0] * (1 - target[j]) / (1 - a)
    }
  }
  q
}
h4_truth <- true_table(h4_sim) |>
  filter(score == "CCS")
h4_independent <- h4_ct |>
  mutate(
    independent = map2_dbl(control_sample_id, stim_sample_id, \(u, s) {
      q <- ipf_independent(h4_sim$p[u, ], h4_sim$p[s, ], h4_bits)
      sum((h4_sim$p[s, ] - q)[is_poly(h4_sim)])
    })
  ) |>
  left_join(
    h4_truth |> filter(grain == "global") |> select(contrast_id, truth_sample),
    join_by(contrast_id),
    relationship = "one-to-one"
  )
h4_pattern_sum <- h4_truth |>
  filter(grain == "pattern") |>
  summarise(pattern_sum = sum(truth_sample), .by = contrast_id) |>
  left_join(h4_independent, join_by(contrast_id), relationship = "one-to-one")
h4_stored <- read_tsv(
  file.path(REDESIGN_DIR, "reps", h4_scenario, "seed_01", "model_scores.tsv"),
  col_select = c(model, draws, contrast_id, grain, score, estimate),
  show_col_types = FALSE
) |>
  filter(score == "CCS", grain %in% c("global", "pattern")) |>
  summarise(
    pattern_total = sum(estimate[grain == "pattern"]),
    global_value = estimate[grain == "global"],
    .by = c(model, draws, contrast_id)
  )
h4_implementation <- tibble(
  check = c(
    "score_contrast() vs cj_score_reference(), max abs diff over 6 global scores, 8 contrasts x 200 C0 draws",
    "independent IPF vs stored truth_sample CCS, max abs diff over contrasts",
    "pattern CCS truths sum to global CCS truth, max abs diff",
    "stored pattern CCS estimates sum to the global CCS estimate, max abs diff over models (C0 draws)"
  ),
  value = c(
    max(h4_engine$max_abs_diff),
    max(abs(h4_independent$independent - h4_independent$truth_sample)),
    max(abs(h4_pattern_sum$pattern_sum - h4_pattern_sum$truth_sample)),
    max(abs(h4_stored$pattern_total - h4_stored$global_value), na.rm = TRUE)
  ),
  tolerance = c(1e-9, 1e-9, 1e-12, 1e-9)
) |>
  mutate(pass = value <= tolerance)

# Pooled reference -------------------------------------------------------------
# Proposed estimand: CCS against a population baseline (the control counts pooled over donors). Compared with the
# per-donor control (C0 a0 = 0.5), against three truths: sample (realised control), struct (structured), and pool
# (structured control x realised stimulated sample, the donor-specific coordination change vs the population baseline).
pooled_reference <- sims |>
  mutate(
    table = pmap(list(sim, seed), \(sim, seed) {
      set.seed(seed)
      bits <- sim$data$bits
      ct <- contrasts_of(sim)
      pooled_y <- colSums(sim$data$counts[
        unique(ct$control_sample_id),
        ,
        drop = FALSE
      ])
      pooled_u <- draw_dirichlet(pooled_y, 0.5, N_DRAWS) # one pooled baseline, shared by the donors
      pmap(
        list(ct$contrast_id, ct$control_sample_id, ct$stim_sample_id),
        \(id, u, s) {
          p_s <- draw_dirichlet(sim$data$counts[s, ], 0.5, N_DRAWS)
          own <- ccs_summary(score_contrast(
            draw_dirichlet(sim$data$counts[u, ], 0.5, N_DRAWS),
            p_s,
            bits
          )$G[, "CCS"])
          pool <- ccs_summary(score_contrast(pooled_u, p_s, bits)$G[, "CCS"])
          bind_rows(own = own, pooled = pool, .id = "reference") |>
            mutate(
              contrast_id = id,
              truth_sample = ccs_of(sim$p[u, ], sim$p[s, ], bits),
              truth_struct = ccs_of(sim$p_struct[u, ], sim$p_struct[s, ], bits),
              truth_pool = ccs_of(sim$p_struct[u, ], sim$p[s, ], bits)
            )
        }
      ) |>
        list_rbind()
    })
  ) |>
  select(scenario_id, seed, table) |>
  unnest(table) |>
  pivot_longer(
    starts_with("truth_"),
    names_to = "truth",
    values_to = "value"
  ) |>
  summarise(
    coverage = mean(value >= lower & value <= upper),
    bias = mean(estimate - value),
    rmse = sqrt(mean((estimate - value)^2)),
    mean_posterior_sd = mean(sd),
    sd_of_truth = sd(value),
    false_sign = mean((p_positive > SIGN_PROB | p_negative > SIGN_PROB)[
      abs(value) <= 1e-12
    ]), # NaN when no zero truth
    mean_failure_fraction = mean(failure_fraction),
    .by = c(scenario_id, reference, truth)
  )

# Write outputs ----------------------------------------------------------------
write_tsv(regeneration_check, file.path(run_dir, "regeneration_check.tsv"))
write_tsv(h1_truth_variants, file.path(run_dir, "h1_truth_variants.tsv"))
write_tsv(h1_truth_spread, file.path(run_dir, "h1_truth_spread.tsv"))
write_tsv(h1_model_error, file.path(run_dir, "h1_model_error.tsv"))
write_tsv(control_information, file.path(run_dir, "control_information.tsv"))
write_tsv(h3_a0_sweep, file.path(run_dir, "h3_a0_sweep.tsv"))
write_tsv(h2_sensitivity, file.path(run_dir, "h2_sensitivity.tsv"))
write_tsv(h4_implementation, file.path(run_dir, "h4_implementation.tsv"))
write_tsv(pooled_reference, file.path(run_dir, "pooled_reference.tsv"))
write_session_info(file.path(run_dir, "session_info.txt"), PKG_DIR)
