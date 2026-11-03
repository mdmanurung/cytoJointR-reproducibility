# Donor-effect candidates (2026-11-04)

Raw scores of the leaderboard rows for the candidates of branch `feat/donor-effects` (package worktrees under `cytoJointR-wt/`). Dev mode only; test mode was not run (test-set access counter stays 0).
All candidates run `model=m1dm,level=structured` plus the flag shown. Run with `scripts/slurm_scorecard.sh` (rv144 task on `all` or `short`, the rest on `short`), collected on the login node.

| Label | Package | Extra fit arguments | Note |
|---|---|---|---|
| `c1-marker` | f02b669 (`cytoJointR-wt/donor-effects`) | `donor_loading=marker` | RV144 R-hat 1.06, ESS 87 in the PPC fit |
| `c2-responder` | aece1a4 (`de-c2`) | `responder=TRUE` | RV144 R-hat 1.73, ESS 6: two chains sat in a worse mode. Rows kept as a record; superseded by c6. |
| `c3-degree` | 139d835 (`de-c34`) | `degree_terms=TRUE` | |
| `c4-scaled` | 139d835 (`de-c34`) | `donor_loading=scaled` | Sim and TB stages ran; the RV144 fit (over 57 min) did not run in time on `short` and the `all` job never left the queue. No leaderboard rows. |
| `c5-responder-init` | dec9e56 (`de-c5`) | `responder=TRUE` | Starting values did not remove the bad mode (R-hat 1.53 in the PPC fit). Cancelled after the PPC fit; no leaderboard rows. Its PPC draws show chain 4 in the bad mode. |
| `c6-responder-screen` | 3188b90 (`de-c6`) | `responder=TRUE` | Chains far below the best chain's log posterior are dropped (`max_lp_gap`). |
| `c7-responder-degree` | 3188b90 (`de-c6`) | `responder=TRUE,degree_terms=TRUE` | |

`ppc-<label>/rv144/` holds the `60_diagnose-estimand.R` outputs (set `CJ_FIT_ARGS` and `CJ_PKG_DIR`): misfit_summary.tsv is the posterior predictive check on RV144, performance.tsv the all-donor (not dev-only) AUC and rule rates.
`ppc-c4-scaled/` is an empty directory: that PPC fit timed out at 60 min.
Table of all candidates: `Rscript scripts/71_compare-candidates.R`.
