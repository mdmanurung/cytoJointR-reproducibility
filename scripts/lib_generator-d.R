# Purpose: held-out generator D, pre-registered before any candidate-model result was seen. It is for
#          vetoing a winner only and is never used to tune a model (plan: plan-to-fix-all-fuzzy-thacker.md).
#          D1: cell-level phenotype mixture. Cells belong to one of H = 3 phenotypes whose markers are
#              independent given the phenotype, so the pattern probabilities are not log-linear and the
#              sample noise (a Dirichlet on phenotype weights) is neither Dirichlet-multinomial nor logistic-normal.
#          D2: generator A with multivariate t (3 df) donor effects instead of normal ones.
# Run:     source("scripts/lib_generator-d.R") after pkgload::load_all("cytoJointR")
# Inputs:
#   sc                          scenario row: generator ("D1" or "D2"), n_donors, J, N, missing_control,
#                               tau_basal, tau_induced, omega_rho (D2)
#   truth                       M1 truth list from bench_truth_a(sc) [helper02_benchmark-sim.R]; alpha is the
#                               baseline logit of each marker (D1 and D2), gamma the induced shift (D2)
#   seed                        integer data seed (seed offset 900000 + 1000 * scenario index + replicate seed)
# Outputs:
#   simulate_generator_d(sc, truth, seed)   list(data = cj_data, p = realised sample probabilities,
#                                           p_struct = structured probabilities), samples by patterns.
#                                           D1 p_struct: the donor's phenotype weights (before sample noise).
#                                           D2 p_struct = p (donor effects are part of the structure).
#                                           Needs internal new_cj_data. One stimulation, one visit.

# Pre-registered constants (do not change after a candidate result exists) ---------------
D1_SHIFT <- 3 # logit shift of the markers a program turns on
D1_PROGRAMS <- list(c("IFNG", "TNF"), c("IL2", "CD107a")) # H = 3: quiescent plus these two
D1_DONOR_CONC <- 50 # donor weights ~ Dir(50 * base weights)
D1_BASE_CONTROL <- c(0.96, 0.02, 0.02)
D1_BASE_STIM <- c(0.90, 0.07, 0.03)
D1_SAMPLE_CONC <- 200 # sample weights ~ Dir(200 * donor weights)
D2_DF <- 3 # degrees of freedom of the donor effects

rdirichlet1 <- function(a) {
  g <- rgamma(length(a), a)
  g / sum(g)
}

simulate_generator_d <- function(sc, truth, seed) {
  base <- cj_simulate_counts(
    n_donors = sc$n_donors,
    markers = c("IFNG", "TNF", "IL2", "CD107a")[seq_len(sc$J)],
    stimulations = "CMV",
    visits = "v1",
    N = sc$N,
    truth = truth,
    donor_effects = FALSE, # donor effects are drawn below (D2) or absent (D1)
    donor_visit = FALSE,
    missing_control = sc$missing_control,
    seed = seed
  )
  samples <- base$data$samples
  bits <- base$data$bits
  markers <- colnames(bits)
  is_stim <- samples$condition != "control"
  if (sc$generator == "D2") {
    # Donor effects (basal, induced) with M1's covariance, as a multivariate t: normal / sqrt(W / df),
    # rescaled by sqrt(df / (df - 2)) so the SD equals tau. They shift every marker's logit equally,
    # as in M1 (eta gains degree * shift), so the mean structure is generator A's.
    donors <- unique(samples$donor_id)
    tau <- c(sc$tau_basal, sc$tau_induced)
    omega <- matrix(c(1, sc$omega_rho, sc$omega_rho, 1), 2, 2)
    set.seed(seed + 900000L)
    z <- t(chol(omega)) %*% matrix(rnorm(2 * length(donors)), 2)
    w <- rchisq(length(donors), D2_DF)
    effect <- t(z * tau) / sqrt(w / D2_DF) / sqrt(D2_DF / (D2_DF - 2))
    rownames(effect) <- donors
    shift <- effect[samples$donor_id, 1] + is_stim * effect[samples$donor_id, 2]
    p <- m1_softmax(log(base$p) + outer(shift, rowSums(bits)))
    p_struct <- p
  } else {
    # Phenotype marker probabilities: quiescent keeps the baseline logit; each program adds D1_SHIFT to its markers.
    alpha <- unname(truth$alpha)
    profile <- rbind(
      alpha,
      alpha + D1_SHIFT * (markers %in% D1_PROGRAMS[[1]]),
      alpha + D1_SHIFT * (markers %in% D1_PROGRAMS[[2]])
    ) |>
      plogis()
    # Independent markers given the phenotype: product of Bernoulli terms per pattern. Rows: phenotype.
    by_phenotype <- t(apply(profile, 1, \(pi_h) apply(bits, 1, \(b) prod(ifelse(b == 1, pi_h, 1 - pi_h)))))
    set.seed(seed + 900000L)
    donors <- unique(samples$donor_id)
    donor_w <- list(
      control = t(replicate(length(donors), rdirichlet1(D1_DONOR_CONC * D1_BASE_CONTROL))),
      stim = t(replicate(length(donors), rdirichlet1(D1_DONOR_CONC * D1_BASE_STIM)))
    )
    w_struct <- t(vapply(
      seq_len(nrow(samples)),
      \(s) donor_w[[if (is_stim[s]) "stim" else "control"]][match(samples$donor_id[s], donors), ],
      numeric(3)
    ))
    w_sample <- t(apply(w_struct, 1, \(w) rdirichlet1(D1_SAMPLE_CONC * w)))
    floor_norm <- \(p) {
      p <- pmax(p, 1e-300)
      p / rowSums(p)
    }
    p_struct <- floor_norm(w_struct %*% by_phenotype)
    p <- floor_norm(w_sample %*% by_phenotype)
  }
  counts <- t(vapply(
    seq_len(nrow(p)),
    \(s) rmultinom(1, sc$N, p[s, ])[, 1],
    integer(ncol(p))
  ))
  dimnames(counts) <- dimnames(base$data$counts)
  dimnames(p) <- dimnames(p_struct) <- dimnames(counts)
  data <- new_cj_data(
    counts,
    samples,
    base$data$contrasts,
    bits,
    provenance = list(generator = paste0(sc$generator, " (held-out)"), seed = seed)
  )
  list(data = data, p = p, p_struct = p_struct)
}
