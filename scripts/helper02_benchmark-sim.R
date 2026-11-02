# Purpose: helpers shared by the benchmark scripts: build the M1 truth of a scenario row,
#          simulate one data set (generators A, B, C, and the held-out D1, D2), and get the true
#          per-pattern change for the realised (p) and the structured (p_struct) probabilities.
# Run:     source("scripts/helper02_benchmark-sim.R") from the repo root, after pkgload::load_all("cytoJointR")
# Inputs:
#   scenario row               one row of <run_dir>/scenarios.tsv [05_define-benchmark.R]
#   scripts/lib_generator-b.R, scripts/lib_generator-c.R   generators B and C [this repo]
# Outputs:
#   bench_simulate(sc, seed, n_donors)   list(data, p, p_struct): one stimulation (CMV); one visit (v1), or v1 + v2
#                                        when sc$n_visits is 2. p is the realised sample probability, p_struct the
#                                        structured one (equal to p unless the generator adds sample noise)
#   true_pattern_delta(sim, p_name)      one row per contrast x positive pattern: contrast_id, donor_id, pattern_key, degree, delta
#                                        (stimulated minus control true probability, from sim[[p_name]];
#                                        pattern_key = bits, first marker first)

source("scripts/lib_generator-b.R")
source("scripts/lib_generator-c.R")
source("scripts/lib_generator-d.R") # pre-registered held-out generator (commit fed35d9)

# M1 truth list for generators A and C: one stimulation, no visit effect, no interactions beyond theta0.
bench_truth_a <- function(sc) {
  markers <- c("IFNG", "TNF", "IL2", "CD107a")[seq_len(sc$J)]
  alpha <- unlist(sc[paste0("alpha_", seq_len(sc$J))])
  beta <- unlist(sc[paste0("beta_", seq_len(sc$J))])
  n_pairs <- choose(sc$J, 2)
  # Two visits add the columns visit_v2 (no main effect) and visit_v2:stim_CMV (the visit change).
  gamma <- rbind(stim_CMV = beta)
  if (isTRUE(sc$n_visits == 2)) {
    gamma <- rbind(
      gamma,
      visit_v2 = 0 * beta,
      `visit_v2:stim_CMV` = sc$visit_effect + 0 * beta
    )
  }
  colnames(gamma) <- markers
  list(
    alpha = setNames(alpha, markers),
    gamma = gamma,
    theta0 = sc$theta0_start + 0.1 * (seq_len(n_pairs) - 1),
    delta = matrix(0, 1, n_pairs),
    tau = c(sc$tau_basal, sc$tau_induced),
    Omega = matrix(c(1, sc$omega_rho, sc$omega_rho, 1), 2, 2),
    tau_u = if (isTRUE(sc$n_visits == 2)) sc$tau_u else 0
  )
}

# Generator B is written for CMV + SEB and visits v1 + v2 (lib_generator-b.R), so build its full
# truth, simulate, and keep the CMV v1 samples.
bench_simulate_b <- function(sc, seed) {
  markers <- c("IFNG", "TNF", "IL2", "CD107a")[seq_len(sc$J)]
  beta <- unlist(sc[paste0("beta_", seq_len(sc$J))])
  gamma <- rbind(beta, beta, 0 * beta, 0 * beta, 0 * beta) # CMV, SEB, visit_v2 and its two interactions
  dimnames(gamma) <- list(
    c(
      "stim_CMV",
      "stim_SEB",
      "visit_v2",
      "visit_v2:stim_CMV",
      "visit_v2:stim_SEB"
    ),
    markers
  )
  truth <- list(
    alpha = setNames(unlist(sc[paste0("alpha_", seq_len(sc$J))]), markers),
    gamma = gamma,
    theta0 = sc$theta0_start + 0.1 * (seq_len(choose(sc$J, 2)) - 1),
    delta = matrix(0, 2, choose(sc$J, 2)),
    tau = c(sc$tau_basal, sc$tau_induced, sc$tau_induced),
    Omega = diag(3),
    tau_u = 0
  )
  full <- simulate_generator_b(sc, truth, seed)
  d <- full$data
  keep <- d$samples$visit == "v1" & d$samples$condition %in% c("control", "CMV")
  samples <- d$samples[keep, ]
  samples$visit <- droplevels(samples$visit)
  list(
    data = new_cj_data(
      d$counts[keep, ],
      samples,
      d$contrasts[
        d$contrasts$stimulation == "CMV" & d$contrasts$visit == "v1",
      ],
      d$bits,
      provenance = list(
        generator = "B (logistic-normal saturated), CMV v1 subset",
        seed = seed
      )
    ),
    p = full$p[keep, ],
    p_struct = full$p_struct[keep, ]
  )
}

bench_simulate <- function(sc, seed, n_donors = sc$n_donors) {
  sc$n_donors <- n_donors
  sim <- switch(
    sc$generator,
    B = bench_simulate_b(sc, seed),
    C = simulate_generator_c(sc, bench_truth_a(sc), seed),
    D1 = simulate_generator_d(sc, bench_truth_a(sc), seed),
    D2 = simulate_generator_d(sc, bench_truth_a(sc), seed),
    cj_simulate_counts(
      n_donors = sc$n_donors,
      markers = c("IFNG", "TNF", "IL2", "CD107a")[seq_len(sc$J)],
      stimulations = "CMV",
      visits = c("v1", "v2")[seq_len(if (isTRUE(sc$n_visits == 2)) 2 else 1)],
      N = sc$N,
      truth = bench_truth_a(sc),
      donor_visit = isTRUE(sc$n_visits == 2),
      missing_control = sc$missing_control,
      seed = seed
    )
  )
  # A and C add no sample noise beyond p; cj_simulate_counts() with a finite phi returns mu.
  sim$p_struct <- sim$p_struct %||% sim$mu %||% sim$p
  sim
}

true_pattern_delta <- function(sim, p_name = "p") {
  ct <- sim$data$contrasts |> dplyr::filter(!is.na(control_sample_id))
  bits <- sim$data$bits
  pos <- rowSums(bits) >= 1
  dp <- sim[[p_name]][ct$stim_sample_id, ] -
    sim[[p_name]][ct$control_sample_id, ]
  dimnames(dp) <- list(ct$contrast_id, apply(bits, 1, paste, collapse = ""))
  dp[, pos, drop = FALSE] |>
    tibble::as_tibble(rownames = "contrast_id") |>
    tidyr::pivot_longer(
      -contrast_id,
      names_to = "pattern_key",
      values_to = "delta"
    ) |>
    dplyr::left_join(
      ct |> dplyr::select(contrast_id, donor_id),
      dplyr::join_by(contrast_id),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(degree = stringr::str_count(pattern_key, "1"))
}
