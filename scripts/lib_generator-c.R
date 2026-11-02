# Purpose: generator C, a responder mixture: M1 counts in which only a random fraction of
#          donors respond; every other donor has exactly zero induced effect (stimulated
#          probabilities equal the control probabilities). M1 cannot make an exact-zero mixture.
# Run:     source("scripts/lib_generator-c.R") after pkgload::load_all("cytoJointR")
# Inputs:
#   sc, truth, seed             scenario row (n_donors, J, N, resp_frac), M1 truth list, integer seed
# Outputs:
#   simulate_generator_c(sc, truth, seed)   list(data = cj_data, p = true probabilities, responder = named logical)
#                                           like cj_simulate_counts(); needs internal new_cj_data

simulate_generator_c <- function(sc, truth, seed) {
  base <- cj_simulate_counts(
    n_donors = sc$n_donors,
    markers = c("IFNG", "TNF", "IL2", "CD107a")[seq_len(sc$J)],
    stimulations = "CMV",
    visits = "v1",
    N = sc$N,
    truth = truth,
    donor_visit = FALSE,
    missing_control = sc$missing_control,
    seed = seed
  )
  ct <- base$data$contrasts
  set.seed(seed + 700000L)
  responder <- setNames(
    ct$donor_id %in% sample(ct$donor_id, round(sc$resp_frac * nrow(ct))),
    ct$donor_id
  )
  # Non-responders: stimulated probabilities = control probabilities, counts redrawn.
  p <- base$p
  off <- ct$stim_sample_id[!responder[ct$donor_id]]
  p[off, ] <- p[ct$control_sample_id[match(off, ct$stim_sample_id)], ]
  counts <- base$data$counts
  counts[off, ] <- t(vapply(
    off,
    \(s) rmultinom(1, sc$N, p[s, ])[, 1],
    integer(ncol(counts))
  ))
  data <- new_cj_data(
    counts,
    base$data$samples,
    ct,
    base$data$bits,
    provenance = list(generator = "C (responder mixture)", seed = seed)
  )
  list(data = data, p = p, responder = responder)
}
