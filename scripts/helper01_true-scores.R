# Purpose: helpers shared by the pilot scripts: true scores of a simulated data set,
#          and the session-info file that records the cytoJointR commit.
# Run:     source("scripts/helper01_true-scores.R") from the repo root, after pkgload::load_all("cytoJointR")
# Inputs:
#   sim                         output of cj_simulate_counts() or simulate_generator_b()
#   cytoJointR/                 git repo, for the commit hash
# Outputs:
#   true_scores(sim)            one row per contrast with a control and score: contrast_id, donor_id,
#                               visit, stimulation, score, truth
#   write_session_info(path)    writes the cytoJointR commit and sessionInfo() to path

true_scores <- function(sim) {
  sim$data$contrasts |>
    dplyr::filter(!is.na(control_sample_id)) |>
    dplyr::rowwise() |>
    dplyr::mutate(
      global = list(
        cytoJointR::cj_score_reference(
          sim$p[control_sample_id, ],
          sim$p[stim_sample_id, ],
          sim$data$bits
        )$global
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::select(contrast_id, donor_id, visit, stimulation, global) |>
    tidyr::unnest_longer(global, values_to = "truth", indices_to = "score")
}

write_session_info <- function(path) {
  commit <- system2("git", c("-C", "cytoJointR", "rev-parse", "HEAD"), stdout = TRUE)
  writeLines(c(paste("cytoJointR commit:", commit), capture.output(sessionInfo())), path)
}
