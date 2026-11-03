#!/usr/bin/env Rscript
# Purpose: HVTN078 (Env peptide pools, ICS) COMPASS-paper table -> validated cj_data on the 5-marker panel, plus
#          donor metadata. Panel: IFNg, IL2, TNFa, CD154 + CD107a; IL4 and MIP1B are dropped by summing pattern counts.
#          Placebo recipients are the donors with control == 1 (arms P1-P4).
#          FILTER CAVEAT: n_s / n_u were filtered by COMPASS (26 of 128 patterns kept; removed patterns' cells were
#          added to the all-negative pattern), so they are not exact partitions; positive cells in removed
#          patterns count as negative, and any-positive frequencies are a lower bound.
#          The COMPASS call used treatment = ANTIGEN %in% env_stimulations and control = "negctrl"; n_s is
#          therefore pooled over the Env peptide pools and n_u over the two negctrl wells (counts_u ~ 2x one well).
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/21_prepare-hvtn078.R
# Inputs:
#   data/raw/compass_data_package/pkg/COMPASSDataPackage/data/COMPASSData.rda   HVTN078 COMPASSResult [READ-ONLY]
#   <PKG_DIR>/                  package source, loaded with pkgload [frozen worktree, commit 49278c9]
#   scripts/helper03_realdata.R compass_to_cj(), write_counts_long()
# Outputs (in data/processed/hvtn078/):
#   cj_data.rds, counts_long.tsv, donor_meta.tsv, prep_log.tsv   as 20_prepare-rv144.R

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign-run")
RDA <- "data/raw/compass_data_package/pkg/COMPASSDataPackage/data/COMPASSData.rda"
OUT_DIR <- "data/processed/hvtn078"
PANEL <- c("IFNg", "IL2", "TNFa", "CD154", "CD107a") # CD107a is present, so it is the dataset-specific marker
STIM <- "Env"
DONOR_COLS <- c("arm", "grp", "rx_code", "rx", "control") # must be constant within a donor (meta has one row per sample)

pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper03_realdata.R")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load and build ---------------------------------------------------------------
env <- new.env()
load(RDA, env)
d <- env$HVTN078$data
reduced <- compass_to_cj(d, PANEL, STIM)

# Metadata ---------------------------------------------------------------------
# meta is per sample (987 rows, 13 antigens); collapse to one row per PUBID, checking each column is constant.
meta <- d$meta |>
  summarise(across(all_of(DONOR_COLS), \(v) n_distinct(v)), .by = PUBID)
stopifnot("donor-level meta varies within a donor" = all(meta[DONOR_COLS] == 1L))
donor_meta <- tibble(donor_id = rownames(d$n_s)) |>
  left_join(
    d$meta |> distinct(PUBID, across(all_of(DONOR_COLS))) |> rename(donor_id = PUBID),
    join_by(donor_id), relationship = "one-to-one"
  ) |>
  mutate(group = if_else(control == 1, "placebo", "vaccine"))
stopifnot(
  "metadata order differs from the count matrix" = identical(donor_meta$donor_id, rownames(d$n_s)),
  "metadata order differs from the cj_data donors" = identical(donor_meta$donor_id, reduced$x$samples$donor_id[reduced$x$samples$condition == STIM])
)
log <- reduced$log |>
  mutate(
    donors_removed_no_control = sum(!donor_meta$donor_id %in% rownames(d$n_u)),
    donors_with_na_metadata = sum(!complete.cases(donor_meta)),
    n_vaccine = sum(donor_meta$group == "vaccine"), n_placebo = sum(donor_meta$group == "placebo"),
    median_allneg_fraction_stim = median(d$n_s[, ncol(d$n_s)] / d$counts_s),
    rowsums_equal_counts = all(rowSums(d$n_s) == d$counts_s) && all(rowSums(d$n_u) == d$counts_u),
    min_N = min(c(d$counts_s, d$counts_u)), median_N = median(c(d$counts_s, d$counts_u)),
    caveat = "COMPASS-filtered tables: removed patterns' cells folded into all-negative; not exact partitions"
  )

# Write outputs ----------------------------------------------------------------
saveRDS(reduced$x, file.path(OUT_DIR, "cj_data.rds"))
write_counts_long(reduced$x, file.path(OUT_DIR, "counts_long.tsv"))
write_tsv(donor_meta, file.path(OUT_DIR, "donor_meta.tsv"))
write_tsv(log, file.path(OUT_DIR, "prep_log.tsv"))
