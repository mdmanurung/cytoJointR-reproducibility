# Purpose: helpers for the real-data benchmark: turn a COMPASSResult data table into a cj_data on a reduced
#          marker panel, write long counts for polyICSFlow (environment B), and rank AUC with a donor bootstrap.
# Run:     source("scripts/helper03_realdata.R") from the repo root (dplyr, purrr, tidyr, readr attached;
#          cj_from_counts() needs the package loaded for compass_to_cj()).
# Inputs:
#   res_data                    $data of a COMPASSResult: n_s, n_u, counts_s, counts_u, categories
# Outputs:
#   compass_to_cj(res_data, keep, stim, parent, visit)  list(x = cj_data, log = one-row tibble of the collapse)
#   write_counts_long(x, path)  long counts (sample, donor, visit, parent, condition, N, pattern_key, n) as TSV
#   auc_rank(score, truth)      P(score_case > score_noncase) + 0.5 P(tie), truth is logical (TRUE = case)
#   auc_boot(score, truth, B, seed)  AUC with a percentile CI from resampling donors within case and non-case

compass_to_cj <- function(res_data, keep, stim, parent = "parent", visit = "v1") {
  cats <- res_data$categories
  all_markers <- setdiff(colnames(cats), "Counts")
  stopifnot(
    "keep markers not in categories" = all(keep %in% all_markers),
    # categories rows and n_s columns are the same patterns in the same order: COMPASS names
    # ("TNFa&!IFNg") or, in HVTN078, 0/1 strings in categories column order
    "categories do not match n_s columns" = identical(
      apply(cats[, all_markers, drop = FALSE], 1, \(b) paste0(ifelse(b == 1, "", "!"), all_markers, collapse = "&")),
      colnames(res_data$n_s)
    ) || identical(apply(cats[, all_markers, drop = FALSE], 1, paste, collapse = ""), colnames(res_data$n_s)),
    "rowSums(n_s) != counts_s" = all(rowSums(res_data$n_s) == res_data$counts_s),
    "rowSums(n_u) != counts_u" = all(rowSums(res_data$n_u) == res_data$counts_u)
  )
  key <- apply(cats[, keep, drop = FALSE], 1, paste, collapse = "")
  collapse <- \(n) t(rowsum(t(n), key)) # donors x reduced patterns, exact sum over dropped markers
  one <- \(n, donor, condition, N) {
    m <- collapse(n)
    bits <- do.call(rbind, strsplit(colnames(m), "")) == "1"
    colnames(bits) <- keep
    tibble(sample_id = rep(paste(donor, condition, sep = "_"), each = ncol(m))) |>
      bind_cols(as_tibble(bits * 1L)[rep(seq_len(ncol(m)), times = length(donor)), ]) |>
      mutate(count = as.vector(t(m)))
  }
  donors <- rownames(res_data$n_s)
  stopifnot(identical(donors, rownames(res_data$n_u)))
  counts <- bind_rows(one(res_data$n_s, donors, "stim", NULL), one(res_data$n_u, donors, "ctrl", NULL))
  samples <- bind_rows(
    tibble(sample_id = paste(donors, "stim", sep = "_"), donor_id = donors, condition = stim, N = res_data$counts_s),
    tibble(sample_id = paste(donors, "ctrl", sep = "_"), donor_id = donors, condition = "control", N = res_data$counts_u)
  ) |>
    mutate(visit = factor(visit, levels = visit), parent_population = parent)
  x <- cj_from_counts(counts, samples, keep)
  list(
    x = x,
    log = tibble(
      markers_native = paste(all_markers, collapse = ","),
      markers_kept = paste(keep, collapse = ","),
      markers_dropped = paste(setdiff(all_markers, keep), collapse = ","),
      patterns_native_kept_by_compass = nrow(cats),
      patterns_native_possible = 2^length(all_markers),
      patterns_reduced_possible = 2^length(keep),
      patterns_reduced_nonzero = sum(colSums(x$counts) > 0),
      n_donors = length(donors)
    )
  )
}

write_counts_long <- function(x, path) {
  keys <- apply(x$bits, 1, paste, collapse = "")
  x$counts |>
    as_tibble(rownames = "sample_id") |>
    pivot_longer(-sample_id, names_to = "label", values_to = "n") |>
    mutate(pattern_key = keys[match(label, colnames(x$counts))]) |>
    left_join(as_tibble(x$samples) |> select(sample_id, donor_id, visit, parent_population, condition, N),
              join_by(sample_id), relationship = "many-to-one") |>
    mutate(visit = as.character(visit)) |>
    select(sample_id, donor_id, visit, parent_population, condition, N, pattern_key, n) |>
    write_tsv(path)
}

auc_rank <- function(score, truth) {
  ok <- !is.na(score)
  score <- score[ok]; truth <- truth[ok]
  n1 <- sum(truth); n0 <- sum(!truth)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
auc_boot <- function(score, truth, B = 2000L, seed = 1L) {
  set.seed(seed)
  i1 <- which(truth & !is.na(score)); i0 <- which(!truth & !is.na(score))
  b <- replicate(B, {
    i <- c(sample(i1, replace = TRUE), sample(i0, replace = TRUE))
    auc_rank(score[i], truth[i])
  })
  tibble(auc = auc_rank(score, truth), lower = quantile(b, 0.025, na.rm = TRUE), upper = quantile(b, 0.975, na.rm = TRUE),
         n_case = length(i1), n_noncase = length(i0))
}
