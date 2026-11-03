#!/usr/bin/env Rscript
# Purpose: TB_s (MTB-specific antigens) and TB_u (non-specific antigens) from the COMPASS data package, collapsed to
#          the core 5 markers by summing counts over the dropped markers, with the TB status of each donor.
#          COMPASS pre-filtered these tables: patterns it removed sit in the all-negative count, so the
#          tables are not exact partitions of the cells. The log records what the collapse and the filter did.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/22_prepare-tb-compass.R
# Inputs:
#   data/raw/compass_data_package/pkg/COMPASSDataPackage/data/COMPASSData.rda   TB_s, TB_u (COMPASSResult) [md5-verified]
#   <PKG_DIR>/                  package source, loaded with pkgload for cj_from_counts() [frozen worktree 49278c9]
#   scripts/helper03_realdata.R compass_to_cj(), write_counts_long() [vaccine agent; read-only]
# Outputs (in data/processed/<tb_s|tb_u>/):
#   counts.tsv                  long counts: sample_id, donor_id, visit, parent_population, condition, N, pattern_key, n
#   meta.tsv                    donor_id, group (QFT status), case (QFT positive), qft (IU/mL), previous_tb
#   markers.tsv                 marker names in pattern_key order (first character = first marker)
#   prep_log.tsv                one row: the collapse (helper03) plus the filter and removal counts

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign-run") # frozen worktree, commit 49278c9
KEEP <- c("IFNg", "IL2", "TNFa", "CD154", "IL17a") # core 5 (user decision); cytoJointR fits at most 5 markers
STIM <- c(tb_s = "MTB_specific", tb_u = "non_specific") # condition label of the stimulated sample

# Load -------------------------------------------------------------------------
# The rda holds only the two COMPASSResult objects; load it into its own environment.
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper03_realdata.R")
env <- new.env()
load(
  "data/raw/compass_data_package/pkg/COMPASSDataPackage/data/COMPASSData.rda",
  envir = env
)
res <- list(tb_s = env$TB_s$data, tb_u = env$TB_u$data)

# Align meta -------------------------------------------------------------------
# meta is sorted differently from n_s (PUBID1, PUBID10, ...), so align by donor name.
meta <- res |>
  map(\(d) d$meta[match(rownames(d$n_s), d$meta$PATIENT.ID), ])
stopifnot(
  "meta donors match count rows" = identical(
    map(meta, \(m) as.character(m$PATIENT.ID)),
    map(res, \(d) rownames(d$n_s))
  )
)

# Collapse and write -----------------------------------------------------------
# One data set at a time; the log row also counts the cells that the collapse moved to all-negative.
prep_log <- names(res) |>
  set_names() |>
  map(\(ds) {
    d <- res[[ds]]
    out <- file.path("data/processed", ds)
    dir.create(out, recursive = TRUE, showWarnings = FALSE)
    cj <- compass_to_cj(d, KEEP, STIM[[ds]])
    only_dropped <- rowSums(d$categories[, KEEP]) == 0 &
      rowSums(d$categories[, colnames(d$categories) != "Counts"]) > 0
    write_counts_long(cj$x, file.path(out, "counts.tsv"))
    meta[[ds]] |>
      transmute(
        donor_id = as.character(PATIENT.ID),
        group = as.character(TB),
        case = TB == "positive",
        qft = QFT,
        previous_tb = PreviousDiagnosisOfTB
      ) |>
      write_tsv(file.path(out, "meta.tsv"))
    tibble(marker = KEEP) |> write_tsv(file.path(out, "markers.tsv"))
    cj$log |>
      mutate(
        dataset = ds,
        stim_cells_min_max = paste(range(d$counts_s), collapse = "-"),
        n_case = sum(meta[[ds]]$TB == "positive"),
        n_noncase = sum(meta[[ds]]$TB != "positive"),
        patterns_only_dropped_markers_positive = sum(only_dropped),
        cells_moved_to_allneg_stim = sum(d$n_s[, only_dropped]),
        cells_moved_to_allneg_ctrl = sum(d$n_u[, only_dropped]),
        control_N_equals_other_dataset = identical(
          d$counts_u,
          res[[setdiff(names(res), ds)]]$counts_u
        ),
        .before = 1
      )
  }) |>
  list_rbind()
write_tsv(
  prep_log |> filter(dataset == "tb_s"),
  "data/processed/tb_s/prep_log.tsv"
)
write_tsv(
  prep_log |> filter(dataset == "tb_u"),
  "data/processed/tb_u/prep_log.tsv"
)
print(as.data.frame(t(prep_log)))
