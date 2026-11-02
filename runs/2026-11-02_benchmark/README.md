# Benchmark run 2026-11-02: cytoJointR vs COMPASS vs polyICSFlow

Detection-focused benchmark on simulated ICS counts. J = 4 markers, one stimulation (CMV), one visit (v1),
24 donors, N = 10 000 parent cells per sample, 7 scenarios x 30 seeds, 0 missing controls.

## How it was run

1. `scripts/05_define-benchmark.R` (env A) solved each scenario's baseline and effect size and wrote `scenarios.tsv`, `tuning_check.tsv`.
2. Smoke: 7 scenarios x 3 seeds in `runs/2026-11-02_benchmark-smoke` (summary in `-smoke-summary`). Not part of the results.
3. Full: `sbatch --array=1-100 scripts/slurm_benchmark.sh runs/2026-11-02_benchmark 30 0`, then `--array=1-110 ... 30 100`.
   Each task runs `06_run-benchmark-replicate.R` (env A) then `07_run-polyicsflow.R` (env B). Partition `short`, 4 CPUs, 24G.
4. `scripts/08_summarise-benchmark.R` wrote `runs/2026-11-02_benchmark-summary`; `notebooks/benchmark_report.qmd` reads it.

## Versions

| Item | Version |
|---|---|
| cytoJointR | loaded from source at run time; the commit of each replicate is in `session_info.txt` and in `replicate_status.tsv` (`package_commit`). 05 ran at `cb65cc4`. The smoke ran at `ac43752`. |
| Env A R | 4.5.1 (conda `R4_51`); cmdstanr 0.8.0, CmdStan 2.39.0, posterior 1.7.0, COMPASS 1.48.0; library `renv/lib` |
| Env B R | 4.4.3 (`batchdx-cyto/renv`), Bioconductor 3.20; flowWorkspace 4.18.0, flowCore 2.18.0 |
| polyICSFlow | 0.99.4, GitHub SoegaardLab/polyICSFlow, SHA `9dc2b374a299c8f06f8d7fa62fa2366f54d093e3` (HEAD on 2026-11-02) |
| Extra env B library | `renv/lib-r44` (first in `R_LIBS`) holds polyICSFlow, openCyto 2.18.0, flowClust 3.44.0, flowWorkspaceData 3.18.0, posterior, cmdstanr. Nothing was installed into batchdx-cyto. Install logs: `logs/install_*.log`. |

The package moved between commits during the benchmark (`cb65cc4` to `ac43752` to `77b45a2` and later). The commits
between `cb65cc4` and `77b45a2` change `R/31_fit.R` (a private CmdStan `output_dir` per `cj_fit`), CI, NEWS and one
test. They touch neither the Stan model nor the scores. The commit set actually used is `table(package_commit)` in
`replicate_status.tsv`.

## Definitions

- Responder: true donor RMS > `EPS_RESPONSE` = 1e-4. That is 1/N: a difference below one cell in 10 000 cannot show in one sample.
- Responding pattern: true change in pattern probability > `EPS_PATTERN` = 1e-5 (0.1 expected cells per sample).
- cytoJointR scores use epsilon = 1e-4 (donor RMS) and 1e-5 (pattern contribution), the same thresholds as the truth.
- Conventional rules: cytoJointR P > 0.95; COMPASS `Response()` > 0.5 (our choice, COMPASS gives no rule); polyICSFlow clipped net > 0.
- COMPASS magnitude: sum over kept non-null patterns of the posterior-mean difference p_s - p_u (fraction scale, `samplePuPs.cpp:64`).
- polyICSFlow: `calcPolyfunctionality(md_cols = "donor_id", pop_col = "cell_type", condition_col = "condition", background_val = "control")`.
  Clipped = `perc_pos_backgroundSubtracted`; unclipped = `perc_pos` minus the control's `perc_pos`; both /100.
- AUC `vs_null`: the scenario's responders against all donors of the null scenario. AUC `within`: against the scenario's own
  non-responders; defined only where both classes exist (mostly `mixture`, `negative`, `rare_weak`).
- SE of AUC and rates: resampling seeds, 1000 draws. SE of bias, RMSE and coverage: across seeds.

## Deviations and package findings

- COMPASS `SimpleCOMPASS()` does not compute `fit$posterior` (only `COMPASS()` does, `COMPASS.R:470`), and `cj_compass()`
  (`cytoJointR/R/60_compass.R:109`) calls `SimpleCOMPASS()`. So `PosteriorDiff()` fails on `cj_compass()$fits` as returned
  (IMPLEMENTATION_STATUS lists it as untested). 06 calls `COMPASS:::compute_posterior()` after `set.seed(data_seed)`. This adds
  the posterior and leaves FS, PFS and `Response()` untouched.
- COMPASS bug: `PosteriorDiff()` errors when one non-null pattern survives the filter (`Posterior.R:57-59`, `sapply()` returns a vector).
  It failed on null seed 1 in the smoke. 06 sums `p$diff` per donor directly, which equals `rowSums(PosteriorDiff())` otherwise.
- A failed COMPASS fit scores 0 for every donor, pattern and magnitude (COMPASS reports no response), flagged `filled` in `donor_table.tsv`.
  Failures are counted in `replicate_status.tsv`.
- COMPASS default filter: patterns it removes score 0 at pattern level; the `nofilter` variant (`category_filter = NULL`) is reported separately.
- Generator B (`lib_generator-b.R`) is hard-wired to CMV + SEB and visits v1 + v2. It is simulated in full and subset to CMV / v1.
- Generator C (`lib_generator-c.R`, new): 50% of donors respond at the `rare_moderate` size; the rest have stimulated probabilities
  equal to control probabilities, counts redrawn.
- cytoJointR bug (not on the benchmark path): `cj_describe()` errors on a `cj_data` with no contrasts (control-only samples),
  `cytoJointR/R/20_describe.R` in `agg()`, `order(match(out$contrast_id, ...), out$degree)` on an empty `out`.
  Found while running the GatingSet bonus (`bonus_gatingset.txt`).
- cytoJointR fits keep the default `adapt_delta`; some fits have divergences. All fits stay in every metric (see `diagnostics.tsv` of the summary).
- CPUs: cytoJointR 4 (parallel chains), COMPASS 1, polyICSFlow 1. polyICSFlow runtime includes about 0.7 s to expand counts to cells.

## Bonus: cj_from_gatingset()

Ran in env B on `flowWorkspaceData` `gs_manual` (one sample, parent `CD4`, two arbitrary CD4 sub-gates as "markers").
Counts equal an independent cross-tabulation of the gate indices. See `bonus_gatingset.txt`, `scripts/09_bonus-gatingset.R`.

## Not run

(filled in after the run)
