#!/usr/bin/env Rscript
# Purpose: check and tidy the public SARS-CoV-2 Spike ICS Boolean counts (Zenodo 8088178, Ford et al. 2023): verify
#          the downloads against the Zenodo md5, confirm every sample is a complete 16-pattern partition of its
#          parent count, confirm every peptide and PHA sample has one DMSO at the same donor and visit, and write
#          long counts for 32_analyse-sarscov2.R.
# Run:     PATH=<R4_51>/bin:$PATH R_LIBS=renv/lib Rscript scripts/24_prepare-sarscov2.R
# Inputs:
#   data/raw/zenodo_8088178/CD4_ICS_Boolean_Data.xlsx, CD8_ICS_Boolean_Data.xlsx   [Zenodo API, downloaded by hand]
#   data/raw/zenodo_8088178/record.json                                           Zenodo API record with the md5s
#   scripts/helper05_realdata-sarscov2.R                                          workbook reader [this repo]
# Outputs (in data/processed/sarscov2/):
#   counts.tsv         population, sample_id, donor_id, visit, antigen, N, pattern_key, n: 16 rows per sample
#   checksums.tsv      file, size, md5_zenodo, md5_local, match
#   data_checks.tsv    check, population, n_failed, detail: completeness and pairing results
#   unpaired.tsv       peptide or PHA samples without exactly one DMSO at the same donor and visit (header only if none)
#   citation.txt       data record and the paper it belongs to

suppressPackageStartupMessages({
  library(dplyr)
  library(forcats)
  library(here)
  library(purrr)
  library(readr)
  library(stringr)
  library(tidyr)
})

# Parameters -------------------------------------------------------------------
RAW_DIR <- here("data", "raw", "zenodo_8088178") # downloads from the Zenodo API
OUT_DIR <- here("data", "processed", "sarscov2")
POPULATIONS <- c("CD4", "CD8") # one workbook each; cj_fit takes one parent population
N_PATTERNS <- 16L # 2^4 Boolean combinations of CD40L, IFNg, IL2, TNFa

source(here("scripts", "helper05_realdata-sarscov2.R"))

# Checksums --------------------------------------------------------------------
# Zenodo lists md5 per file; the files are untrusted until they match.
record <- jsonlite::read_json(file.path(RAW_DIR, "record.json"))
checksums <- record$files |>
  map(\(f) {
    tibble(
      file = f$key,
      size = f$size,
      md5_zenodo = str_remove(f$checksum, "^md5:")
    )
  }) |>
  list_rbind() |>
  mutate(
    md5_local = unname(tools::md5sum(file.path(RAW_DIR, file))),
    match = md5_local == md5_zenodo
  )
if (!all(checksums$match)) {
  stop(
    "md5 mismatch: ",
    paste(checksums$file[!checksums$match], collapse = ", "),
    call. = FALSE
  )
}

# Read -------------------------------------------------------------------------
counts <- POPULATIONS |>
  map(\(pop) {
    ics_read_boolean(
      file.path(RAW_DIR, sprintf("%s_ICS_Boolean_Data.xlsx", pop)),
      pop
    )
  }) |>
  list_rbind()

# Completeness -----------------------------------------------------------------
# Each sample: 16 distinct patterns, no NA, nonnegative integers, summing to the parent count (column "Count").
per_sample <- counts |>
  summarise(
    n_patterns = n_distinct(pattern_key),
    n_rows = n(),
    any_na = anyNA(n),
    any_bad = any(n < 0 | n != round(n)),
    pattern_total = sum(n),
    N = first(N),
    .by = c(population, sample_id)
  )
# Sample metadata meets counts: the parent count in the sheet is the sum of the 16 patterns.
stopifnot(
  "parent count equals the sum of the 16 patterns" = identical(
    per_sample$N,
    per_sample$pattern_total
  )
)
completeness <- per_sample |>
  summarise(
    n_samples = n(),
    n_not_16_patterns = sum(n_patterns != N_PATTERNS | n_rows != N_PATTERNS),
    n_na_or_negative = sum(any_na | any_bad),
    min_N = min(N),
    .by = population
  )

# Pairing ----------------------------------------------------------------------
# Every peptide and PHA sample needs exactly one DMSO at the same donor and visit in the same population.
samples <- counts |>
  distinct(population, sample_id, donor_id, visit, antigen, N)
dmso_n <- samples |>
  filter(antigen == "DMSO") |>
  count(population, donor_id, visit, name = "n_dmso")
unpaired <- samples |>
  filter(antigen != "DMSO") |>
  left_join(
    dmso_n,
    join_by(population, donor_id, visit),
    relationship = "many-to-one"
  ) |>
  filter(is.na(n_dmso) | n_dmso != 1)
pairing <- samples |>
  summarise(
    n_stimulated = sum(antigen != "DMSO"),
    n_dmso = sum(antigen == "DMSO"),
    .by = population
  ) |>
  mutate(n_unpaired = map_int(population, \(p) sum(unpaired$population == p)))

# Write outputs ----------------------------------------------------------------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_tsv(counts, file.path(OUT_DIR, "counts.tsv"))
write_tsv(checksums, file.path(OUT_DIR, "checksums.tsv"))
write_tsv(
  bind_rows(
    completeness |>
      pivot_longer(-population, names_to = "check", values_to = "value") |>
      mutate(group = "completeness"),
    pairing |>
      pivot_longer(-population, names_to = "check", values_to = "value") |>
      mutate(group = "pairing")
  ),
  file.path(OUT_DIR, "data_checks.tsv")
)
write_tsv(unpaired, file.path(OUT_DIR, "unpaired.tsv"))
writeLines(
  c(
    "Data: Ford ES, Mayer-Blackwell K, Jing L, Laing KJ, Burrow CJ, Koelle DM (2023). CD4+ and CD8+ T cell responses",
    "  to peptides covering SARS-CoV-2 Spike in response to mRNA vaccination in persons recovered from SARS-CoV-2",
    "  infection. Zenodo, doi:10.5281/zenodo.8088178 (CC BY 4.0), published 2023-06-27.",
    "Paper (same first author, group and ICS design; E01 n=7, E02 n=14, E03 n=14 as in the record): Ford ES et al.",
    "  Repeated mRNA vaccination sequentially boosts SARS-CoV-2-specific CD8 T cells in persons with previous",
    "  COVID-19. Nat Immunol 25(1):166-177 (2024), doi:10.1038/s41590-023-01692-x, PMID 38057617. The Zenodo record",
    "  carries no related identifier; the link rests on authors, design and sample sizes.",
    "Visits (paper): E01 before vaccine dose 1; E02 about 19.5 days after dose 1; E03 about 24 days after dose 2."
  ),
  file.path(OUT_DIR, "citation.txt")
)
print(checksums)
print(completeness)
print(pairing)
