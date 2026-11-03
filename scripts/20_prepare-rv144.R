#!/usr/bin/env Rscript
# Purpose: RV144 (Env, CD4 ICS) COMPASS-paper table -> validated cj_data on the 5-marker panel, plus donor metadata.
#          Panel: IFNg, IL2, TNFa, CD154 + IL4 (dataset-specific); IL17a is dropped by summing pattern counts.
#          FILTER CAVEAT: n_s / n_u were filtered by COMPASS (24 of 64 patterns kept; removed patterns' cells were
#          added to the all-negative pattern), so they are not exact partitions: rowSums(n_s) == counts_s holds
#          by construction, and positive cells in removed patterns are counted as negative. Any-positive
#          frequencies are therefore a lower bound. Collapsing preserves this.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/20_prepare-rv144.R
# Inputs:
#   data/raw/compass_data_package/pkg/COMPASSDataPackage/data/COMPASSData.rda   rv144 COMPASSResult [READ-ONLY; zip md5 6a82ee0c99b3ba326c0e62b4166a98db]
#   <PKG_DIR>/                  package source, loaded with pkgload [frozen worktree, commit 49278c9]
#   scripts/helper03_realdata.R compass_to_cj(), write_counts_long()
# Outputs (in data/processed/rv144/):
#   cj_data.rds                 cj_data, 5 markers
#   cj_data_native.rds          cj_data, all 6 native markers (reference COMPASS only)
#   counts_long.tsv             long counts of cj_data, input of polyICSFlow (environment B)
#   donor_meta.tsv              one row per donor: vaccine arm, infection, covariates
#   prep_log.tsv                what was removed, the marker collapse and the filter facts

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign-run") # frozen package copy, commit 49278c9
RDA <- "data/raw/compass_data_package/pkg/COMPASSDataPackage/data/COMPASSData.rda"
OUT_DIR <- "data/processed/rv144"
PANEL <- c("IFNg", "IL2", "TNFa", "CD154", "IL4") # core 4 + IL4 for RV144 (user decision)
STIM <- "Env" # the COMPASS call used treatment = Stim == "92TH023 Env"
META_COLS <- c("PUBID", "vaccine", "infect", "sex", "gender", "age", "stratuminds", "risk.medium", "risk.high",
               "IgAprim", "V2prim", "flrstatus", "Luminex") # covariates for the infection analysis

pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper03_realdata.R")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Load and build ---------------------------------------------------------------
env <- new.env()
load(RDA, env) # data only; nothing from the file is executed
d <- env$rv144$data
reduced <- compass_to_cj(d, PANEL, STIM)
native <- compass_to_cj(d, setdiff(colnames(d$categories), "Counts"), STIM)

# Metadata ---------------------------------------------------------------------
# meta has one row per donor (PUBID); joined to the count matrix donors, one-to-one.
donors <- tibble(donor_id = rownames(d$n_s))
meta <- d$meta |> select(all_of(META_COLS)) |> rename(donor_id = PUBID)
stopifnot("meta is not one row per donor" = !anyDuplicated(meta$donor_id))
donor_meta <- donors |>
  left_join(meta, join_by(donor_id), relationship = "one-to-one") |>
  mutate(arm = tolower(vaccine), infected = if_else(infect == "INFECTED", 1L, if_else(infect == "NON-INFECTED", 0L, NA_integer_)))
stopifnot(
  "metadata order differs from the count matrix" = identical(donor_meta$donor_id, rownames(d$n_s)),
  "metadata order differs from the cj_data donors" = identical(donor_meta$donor_id, reduced$x$samples$donor_id[reduced$x$samples$condition == STIM])
)
removed_no_control <- sum(!donors$donor_id %in% rownames(d$n_u))
na_meta <- sum(!complete.cases(donor_meta |> select(-infected, -infect)))
log <- reduced$log |>
  mutate(
    donors_removed_no_control = removed_no_control,
    donors_with_na_metadata = na_meta, # NA in a covariate; kept here, dropped per model in 30
    n_vaccine = sum(donor_meta$arm == "vaccine"), n_placebo = sum(donor_meta$arm == "placebo"),
    n_vaccine_infected = sum(donor_meta$infected == 1 & donor_meta$arm == "vaccine"),
    median_allneg_fraction_stim = median(d$n_s[, ncol(d$n_s)] / d$counts_s),
    rowsums_equal_counts = all(rowSums(d$n_s) == d$counts_s) && all(rowSums(d$n_u) == d$counts_u),
    min_N = min(c(d$counts_s, d$counts_u)), median_N = median(c(d$counts_s, d$counts_u)),
    caveat = "COMPASS-filtered tables: removed patterns' cells folded into all-negative; not exact partitions"
  )

# Write outputs ----------------------------------------------------------------
saveRDS(reduced$x, file.path(OUT_DIR, "cj_data.rds"))
saveRDS(native$x, file.path(OUT_DIR, "cj_data_native.rds"))
write_counts_long(reduced$x, file.path(OUT_DIR, "counts_long.tsv"))
write_tsv(donor_meta, file.path(OUT_DIR, "donor_meta.tsv"))
write_tsv(log, file.path(OUT_DIR, "prep_log.tsv"))
