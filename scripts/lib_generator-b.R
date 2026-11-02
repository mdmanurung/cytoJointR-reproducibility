# Generator B: counts from logistic-normal saturated pattern probabilities.
# Not in the package. Misspecification for M1: each pattern's log-probability gets
# donor, donor-by-stimulation and sample noise that M1's unit-loading effects cannot
# reproduce, with a per-donor noise scale (heterogeneous donors).
# Needs cytoJointR loaded with pkgload::load_all() (uses internals new_cj_data, m1_softmax).
#
# sc:   one row of scenarios.tsv
# gamma, truth: the M1 truth list that gives the mean structure (built in 03)
# Returns list(data = cj_data, p = realised sample probabilities, p_struct = the structured mean before the noise,
# both samples by patterns), like cj_simulate_counts().
simulate_generator_b <- function(sc, truth, seed) {
  base <- cj_simulate_counts(
    n_donors = sc$n_donors,
    markers = c("IFNG", "TNF", "IL2", "CD107a")[seq_len(sc$J)],
    stimulations = c("CMV", "SEB"),
    visits = c("v1", "v2"),
    N = sc$N,
    truth = truth,
    donor_effects = FALSE,
    donor_visit = FALSE,
    missing_control = sc$missing_control,
    seed = seed
  )
  samples <- base$data$samples
  K <- ncol(base$p)
  set.seed(seed + 500000L)
  donors <- unique(samples$donor_id)
  scale_i <- setNames(exp(rnorm(length(donors), 0, sc$gb_scale_sd)), donors)
  noise <- function(n) matrix(rnorm(n * K), n, K)
  g_donor <- noise(length(donors)) * sc$gb_sd_donor * scale_i
  rownames(g_donor) <- donors
  stim_key <- paste(rep(donors, each = 2), c("CMV", "SEB"))
  h_stim <- noise(length(stim_key)) * sc$gb_sd_stim * rep(scale_i, each = 2)
  rownames(h_stim) <- stim_key
  h_row <- paste(samples$donor_id, samples$condition)
  h <- h_stim[match(h_row, stim_key), , drop = FALSE]
  h[is.na(h)] <- 0 # control samples carry no stimulation noise
  eta <- log(base$p) +
    g_donor[samples$donor_id, , drop = FALSE] +
    h +
    noise(nrow(samples)) * sc$gb_sd_sample
  eta[, rowSums(base$data$bits) == 0] <- log(base$p)[,
    rowSums(base$data$bits) == 0
  ]
  p <- m1_softmax(eta)
  counts <- t(vapply(
    seq_len(nrow(p)),
    \(s) rmultinom(1, sc$N, p[s, ])[, 1],
    integer(K)
  ))
  dimnames(counts) <- dimnames(base$data$counts)
  data <- new_cj_data(
    counts,
    samples,
    base$data$contrasts,
    base$data$bits,
    provenance = list(generator = "B (logistic-normal saturated)", seed = seed)
  )
  list(
    data = data,
    p = `dimnames<-`(p, dimnames(counts)),
    p_struct = `dimnames<-`(base$p, dimnames(counts))
  )
}
