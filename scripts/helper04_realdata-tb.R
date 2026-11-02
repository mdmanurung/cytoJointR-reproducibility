# Purpose: helpers shared by the TB real-data scripts and notebook: read a prepared TB data set, and build a
#          cj_data from its long counts.
# Run:     source("scripts/helper04_realdata-tb.R") from the repo root (dplyr, readr attached; tb_cj() needs
#          cytoJointR loaded)
# Inputs:
#   data/processed/<dataset>/counts.tsv, meta.tsv, markers.tsv   written by 22 and 23 [this repo]
# Outputs:
#   tb_read(dataset)            list(counts, meta, markers): long counts, one row per donor, marker names in bit order
#   tb_cj(counts, markers)      cj_data built from the long counts

tb_read <- function(dataset) {
  dir <- file.path("data/processed", dataset)
  list(
    counts = read_tsv(
      file.path(dir, "counts.tsv"),
      col_types = cols(pattern_key = "c", donor_id = "c", sample_id = "c"),
      progress = FALSE
    ),
    meta = read_tsv(
      file.path(dir, "meta.tsv"),
      col_types = cols(donor_id = "c"),
      progress = FALSE
    ),
    markers = read_tsv(
      file.path(dir, "markers.tsv"),
      col_types = cols(.default = "c"),
      progress = FALSE
    )$marker
  )
}

tb_cj <- function(counts, markers) {
  bits <- do.call(rbind, strsplit(counts$pattern_key, "")) == "1"
  colnames(bits) <- markers
  samples <- counts |>
    distinct(sample_id, donor_id, visit, parent_population, condition, N) |>
    mutate(visit = factor(visit, levels = unique(visit)))
  cj_from_counts(
    bind_cols(
      counts |> select(sample_id),
      as_tibble(bits * 1L),
      count = counts$n
    ),
    samples,
    markers
  )
}
