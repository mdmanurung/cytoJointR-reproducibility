#!/usr/bin/env Rscript
# Purpose: bonus check that cytoJointR::cj_from_gatingset() runs on a real GatingSet (flowWorkspaceData's
#          gs_manual, environment B) and that its pattern counts equal counts taken straight from the gate indices.
#          The two "marker" nodes are arbitrary CD4 sub-gates, so this tests the function, not the biology.
# Run:     R_LIBS=renv/lib-r44 /exports/para-lipg-hpc/mdmanurung/batchdx-cyto/renv/bin/Rscript scripts/09_bonus-gatingset.R runs/2026-11-02_benchmark
# Inputs:
#   renv/lib-r44/flowWorkspaceData   installed from Bioconductor 3.20 [logs/install_flowWorkspaceData.log]
#   cytoJointR/                      package source, loaded with pkgload [this repo]
# Outputs (in <run_dir>):
#   bonus_gatingset.txt              counts from cj_from_gatingset(), the independent counts, and whether they match

suppressPackageStartupMessages({
  library(flowWorkspace)
})

# Parameters -------------------------------------------------------------------
PARENT <- "CD4" # parent gate of the cells that are classified
MARKER_NODES <- c(A_38pDRp = "CD4/38+ DR+", B_CCR7p45RAp = "CD4/CCR7+ 45RA+") # arbitrary sub-gates standing in for markers

# Run --------------------------------------------------------------------------
run_dir <- commandArgs(trailingOnly = TRUE)[1]
pkgload::load_all("cytoJointR", quiet = TRUE)
gs <- load_gs(file.path(
  system.file("extdata", package = "flowWorkspaceData"),
  "gs_manual"
))
id <- sampleNames(gs)
samples <- data.frame(
  sample_id = id,
  donor_id = "D1",
  visit = factor("v1"),
  parent_population = PARENT,
  condition = "control",
  N = gh_pop_get_count(gs[[id]], PARENT)
)
from_cj <- cj_from_gatingset(gs, PARENT, MARKER_NODES, samples)

# Independent counts: cross-tabulate the two gate indices among parent events.
gh <- gs[[id]]
parent <- gh_pop_get_indices(gh, PARENT)
direct <- table(
  gh_pop_get_indices(gh, MARKER_NODES[[1]])[parent],
  gh_pop_get_indices(gh, MARKER_NODES[[2]])[parent]
)
# cj pattern order has the first marker as the low bit: (--, +-, -+, ++)
direct_vec <- c(
  direct["FALSE", "FALSE"],
  direct["TRUE", "FALSE"],
  direct["FALSE", "TRUE"],
  direct["TRUE", "TRUE"]
)

# Write outputs ----------------------------------------------------------------
writeLines(
  c(
    paste("GatingSet sample:", id, "| parent", PARENT, "| N =", samples$N),
    "cj_from_gatingset() counts:",
    capture.output(print(from_cj$counts)),
    paste("direct counts (--, +-, -+, ++):", paste(direct_vec, collapse = " ")),
    paste(
      "identical:",
      identical(unname(as.integer(from_cj$counts[1, ])), as.integer(direct_vec))
    ),
    capture.output(sessionInfo())
  ),
  file.path(run_dir, "bonus_gatingset.txt")
)
