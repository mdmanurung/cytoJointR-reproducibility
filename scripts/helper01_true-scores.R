# Purpose: helpers shared by the pilot and benchmark scripts: true scores of a simulated data set
#          (against the realised sample probabilities and the structured ones), and the session-info
#          file that records the cytoJointR commit.
# Run:     source("scripts/helper01_true-scores.R") from the repo root, after pkgload::load_all("cytoJointR")
# Inputs:
#   sim                         output of cj_simulate_counts() or a lib_generator-*.R simulator; sim$p is the
#                               realised sample-level probability, sim$p_struct (or sim$mu) the structured one
#   cytoJointR/                 git repo, for the commit hash
# Outputs:
#   true_scores(sim)            one row per contrast with a control and score: contrast_id, donor_id,
#                               visit, stimulation, score, truth (sim$p, global scores only; pilot scripts)
#   true_table(sim)             one row per contrast with a control x score: contrast keys, grain
#                               (global, marker, pattern), score, marker, pattern_key, truth_sample, truth_struct
#   write_session_info(path, pkg_dir)   writes the package commit and sessionInfo() to path

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

# True global, marker and exact-pattern scores of one contrast.
# Pattern scores are RMS (degree >= 1), CPS and CCS (degree >= 2), as cj_patterns() reports them.
contrast_truth <- function(p_u, p_s, bits) {
  ref <- cytoJointR::cj_score_reference(p_u, p_s, bits)
  keys <- apply(bits, 1, paste, collapse = "")
  pattern <- purrr::imap(c(RMS = 1, CPS = 2, CCS = 2), \(min_degree, sc) {
    keep <- ref$patterns$degree >= min_degree
    tibble::tibble(
      grain = "pattern",
      score = sc,
      pattern_key = keys[keep],
      truth = ref$patterns[[sc]][keep]
    )
  }) |>
    purrr::list_rbind()
  dplyr::bind_rows(
    tibble::tibble(
      grain = "global",
      score = names(ref$global),
      truth = unname(ref$global)
    ),
    tibble::tibble(
      grain = "marker",
      score = "marker_magnitude",
      marker = ref$markers$marker,
      truth = ref$markers$magnitude
    ),
    pattern
  )
}

true_table_one <- function(p, sim) {
  sim$data$contrasts |>
    dplyr::filter(!is.na(control_sample_id)) |>
    dplyr::mutate(
      truth = purrr::map2(control_sample_id, stim_sample_id, \(u, s) {
        contrast_truth(p[u, ], p[s, ], sim$data$bits)
      })
    ) |>
    dplyr::select(contrast_id, donor_id, visit, stimulation, truth) |>
    tidyr::unnest(truth)
}

true_table <- function(sim) {
  p_struct <- if (!is.null(sim$p_struct)) {
    sim$p_struct
  } else if (!is.null(sim$mu)) {
    sim$mu
  } else {
    sim$p
  }
  keys <- c(
    "contrast_id",
    "donor_id",
    "visit",
    "stimulation",
    "grain",
    "score",
    "marker",
    "pattern_key"
  )
  true_table_one(sim$p, sim) |>
    dplyr::rename(truth_sample = truth) |>
    dplyr::left_join(
      true_table_one(p_struct, sim) |> dplyr::rename(truth_struct = truth),
      dplyr::join_by(
        contrast_id,
        donor_id,
        visit,
        stimulation,
        grain,
        score,
        marker,
        pattern_key
      ),
      relationship = "one-to-one"
    )
}

write_session_info <- function(path, pkg_dir = "cytoJointR") {
  commit <- system2("git", c("-C", pkg_dir, "rev-parse", "HEAD"), stdout = TRUE)
  writeLines(
    c(paste("cytoJointR commit:", commit), capture.output(sessionInfo())),
    path
  )
}
