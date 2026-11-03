# Purpose: helpers for the SARS-CoV-2 Spike ICS analysis (Zenodo 8088178): read one Boolean-count workbook,
#          build a cj_data from the processed tables, and draw posterior predictive counts.
# Run:     source("scripts/helper05_realdata-sarscov2.R") from the repo root
# Inputs:
#   data/raw/zenodo_8088178/<POP>_ICS_Boolean_Data.xlsx   Boolean counts and frequencies [Zenodo 8088178]
# Outputs:
#   ics_read_boolean(path, pop)      long counts (sample_id, donor_id, visit, antigen, N, pattern_key, n)
#   ics_to_cj(counts, pop)           cj_data of one parent population, stimulations Spike and PHA, one control
#   ppc_replicate(mu, phi, N, seed)  one replicated count matrix: Dirichlet-multinomial (phi given) or multinomial

ICS_MARKERS <- c("CD40L", "IFNg", "IL2", "TNFa") # marker order = bit order of cj patterns and of pattern_key

# Read one workbook: the 16 Boolean columns become a pattern_key (bits in ICS_MARKERS order).
ics_read_boolean <- function(path, pop) {
  raw <- readxl::read_excel(path)
  bool_cols <- names(raw)[
    stringr::str_detect(names(raw), "\\| Count$") &
      stringr::str_count(names(raw), "[+-]CD") == 3
  ]
  key_of <- \(nm) {
    paste0(
      purrr::map_chr(ICS_MARKERS, \(m) {
        if (stringr::str_detect(nm, stringr::fixed(paste0(m, "+")))) {
          "1"
        } else {
          "0"
        }
      }),
      collapse = ""
    )
  }
  raw |>
    dplyr::mutate(sample_id = paste(pop, SUBJ_ID, VISIT, ANTIGEN, sep = "_")) |>
    dplyr::select(
      sample_id,
      donor_id = SUBJ_ID,
      visit = VISIT,
      antigen = ANTIGEN,
      N = 4,
      dplyr::all_of(bool_cols)
    ) |>
    dplyr::mutate(donor_id = paste0("D", donor_id)) |>
    tidyr::pivot_longer(
      dplyr::all_of(bool_cols),
      names_to = "column",
      values_to = "n"
    ) |>
    dplyr::mutate(
      pattern_key = purrr::map_chr(column, key_of),
      population = pop
    ) |>
    dplyr::select(-column)
}

# cj_data for one parent population; DMSO is the control, stored once per donor and visit.
ics_to_cj <- function(counts, pop) {
  cp <- counts |> dplyr::filter(population == pop)
  samples <- cp |>
    dplyr::distinct(sample_id, donor_id, visit, antigen, N) |>
    dplyr::transmute(
      sample_id,
      donor_id,
      visit = forcats::fct(visit, levels = c("E01", "E02", "E03")),
      parent_population = pop,
      condition = dplyr::case_match(
        antigen,
        "DMSO" ~ "control",
        "peptide" ~ "Spike",
        "PHA" ~ "PHA"
      ),
      N
    )
  bits <- cp |>
    dplyr::mutate(
      !!!rlang::set_names(
        purrr::map(seq_along(ICS_MARKERS), \(j) {
          rlang::expr(as.integer(substr(pattern_key, !!j, !!j)))
        }),
        ICS_MARKERS
      )
    ) |>
    dplyr::select(sample_id, dplyr::all_of(ICS_MARKERS), count = n)
  cytoJointR::cj_from_counts(
    as.data.frame(bits),
    as.data.frame(samples),
    ICS_MARKERS
  )
}

# One replicated count matrix (samples x patterns). mu: draws x samples x patterns at one draw index is passed as
# a samples x patterns matrix. phi = NULL gives a multinomial (m1); a number gives a Dirichlet-multinomial (m1dm).
ppc_replicate <- function(mu, phi, N) {
  t(vapply(
    seq_len(nrow(mu)),
    \(s) {
      p <- if (is.null(phi)) {
        mu[s, ]
      } else {
        g <- stats::rgamma(ncol(mu), shape = phi * mu[s, ])
        if (sum(g) > 0) g / sum(g) else mu[s, ] # tiny shapes can underflow to all zeros
      }
      stats::rmultinom(1, N[s], p)[, 1]
    },
    numeric(ncol(mu))
  ))
}
