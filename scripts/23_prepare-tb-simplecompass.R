#!/usr/bin/env Rscript
# Purpose: the SimpleCOMPASS example (real TB data, Rakshit et al. 2017 Sci Rep, doi 10.1038/s41598-017-10773-5):
#          21 donors, EPTB vs IGRA+, one stimulation (ESAT6/CFP10) against the unstimulated sample, Boolean counts
#          over 8 markers. Collapsed to the 5 markers of the Th17 focus of that paper (IFNg, IL2, TNF, IL17A,
#          IL22; no CD154 in this panel) by summing counts over IL10, IL17F and MIP1b.
#          The 256 patterns are a complete partition of the parent count (checked below, logged).
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/23_prepare-tb-simplecompass.R
# Inputs:
#   system.file("extdata", "SimpleCOMPASSExamples.csv", package = "COMPASS")   42 samples x 256 pattern counts
#   <PKG_DIR>/                  package source, loaded with pkgload for cj_from_counts() [frozen worktree 49278c9]
#   scripts/helper03_realdata.R compass_to_cj(), write_counts_long() [vaccine agent; read-only]
# Outputs (in data/processed/tb_simplecompass/):
#   counts.tsv, meta.tsv, markers.tsv, prep_log.tsv   as 22_prepare-tb-compass.R; meta$case is EPTB

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
PKG_DIR <- Sys.getenv("CJ_PKG_DIR", "cytoJointR-wt/redesign-run") # frozen worktree, commit 49278c9
KEEP <- c("IFNg", "IL2", "TNF", "IL17A", "IL22") # Th17 focus of the paper; cytoJointR fits at most 5 markers
STIM <- "ESAT6_CFP10" # slash-free label of ESAT6/CFP10
NATIVE <- c("IFNg", "IL2", "IL10", "IL17A", "IL17F", "IL22", "MIP1b", "TNF") # column order of the csv names
CSV_TOKENS <- c("IFNg", "2", "10", "17A", "17F", "22", "MIP", "TNF") # the same markers as written in the csv names

# Load -------------------------------------------------------------------------
pkgload::load_all(PKG_DIR, quiet = TRUE)
source("scripts/helper03_realdata.R")
raw <- read_csv(
  system.file("extdata", "SimpleCOMPASSExamples.csv", package = "COMPASS"),
  col_types = cols(.default = "c"),
  name_repair = "minimal",
  progress = FALSE
)
raw <- raw[, -1] |> # unnamed row-number column
  mutate(across(-c(Stimulation, `PATIENT ID`, GROUP), as.numeric))
pattern_cols <- str_which(names(raw), "IFNg[+-]2[+-]")
parent <- raw[[4]]
stopifnot(
  "pattern columns follow the parent count" = identical(pattern_cols, 5:260)
)

# Samples ----------------------------------------------------------------------
# "Unstimulated" and "UnStimulated" are one level; cj_from_counts() wants the control called "control".
samples <- raw |>
  transmute(
    donor_id = `PATIENT ID`,
    group = GROUP,
    condition = if_else(
      str_detect(Stimulation, regex("^unstim", ignore_case = TRUE)),
      "control",
      STIM
    ),
    N = parent
  )
n_per_donor <- samples |> count(donor_id, condition)

# Patterns ---------------------------------------------------------------------
# bits: 0/1 per native marker, parsed from the column name; counts: samples x patterns.
tok_re <- paste0(CSV_TOKENS, "([+-])", collapse = "")
bits <- str_match(names(raw)[pattern_cols], tok_re)[, -1] == "+"
colnames(bits) <- NATIVE
counts <- as.matrix(raw[, pattern_cols])
ord_s <- samples$condition == STIM
ord_u <- samples$condition == "control"
donors <- samples$donor_id[ord_s]
stopifnot(
  "one stimulated and one control sample per donor" = identical(
    donors,
    samples$donor_id[ord_u]
  )
)
rownames(counts) <- samples$donor_id
res_data <- list(
  n_s = counts[ord_s, ],
  n_u = counts[ord_u, ],
  counts_s = setNames(samples$N[ord_s], donors),
  counts_u = setNames(samples$N[ord_u], donors),
  categories = cbind(bits * 1L, Counts = rowSums(bits))
)
colnames(res_data$n_s) <- colnames(res_data$n_u) <- apply(
  bits * 1L,
  1,
  paste,
  collapse = ""
)

# Collapse and write -----------------------------------------------------------
cj <- compass_to_cj(res_data, KEEP, STIM)
out <- "data/processed/tb_simplecompass"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_counts_long(cj$x, file.path(out, "counts.tsv"))
meta <- samples |>
  filter(condition == STIM) |>
  transmute(donor_id, group, case = group == "EPTB")
write_tsv(meta, file.path(out, "meta.tsv"))
write_tsv(tibble(marker = KEEP), file.path(out, "markers.tsv"))
prep_log <- cj$log |>
  mutate(
    dataset = "tb_simplecompass",
    n_case = sum(meta$case),
    n_noncase = sum(!meta$case),
    patterns_complete_partition = all(rowSums(counts) == parent), # 256 patterns sum to the parent count
    csv_samples = nrow(raw),
    csv_control_spellings = paste(
      unique(raw$Stimulation[ord_u]),
      collapse = " / "
    ),
    donors_with_other_than_one_stim_and_one_control = sum(n_per_donor$n != 1),
    cells_moved_to_allneg_stim = sum(res_data$n_s[,
      rowSums(bits[, KEEP]) == 0 & rowSums(bits) > 0
    ]),
    cells_moved_to_allneg_ctrl = sum(res_data$n_u[,
      rowSums(bits[, KEEP]) == 0 & rowSums(bits) > 0
    ]),
    .before = 1
  )
write_tsv(prep_log, file.path(out, "prep_log.tsv"))
print(as.data.frame(t(prep_log)))
