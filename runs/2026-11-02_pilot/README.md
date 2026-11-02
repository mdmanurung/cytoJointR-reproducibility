# Pilot benchmark, 2026-11-02

Pilot stage of the plan in `cytoJointR/docs/SCOPE.md`: generator A (data from
M1) and one generator-B misspecification, 30 seeds per scenario, run as SLURM
arrays. Results are summarised in `cytoJointR/IMPLEMENTATION_STATUS.md`
("Pilot results"); figures are in `pilot_report.html`.

## Code that produced it

- Package: `cytoJointR` loaded with `pkgload::load_all()` from `main`. `main`
  moved during the run, from `5c7b607` (score-speed fix) to `bac00be` (docs
  merge). The commits in between changed only documentation, vignettes and
  build files, not model, score or input code. The replicates did not record
  the commit; scripts now write it to `session_info.txt` from 2026-11-02 on.
- Scripts: `scripts/02_define-scenarios.R`, `03_run-replicate.R`,
  `04_summarise-pilot.R`, `lib_generator-b.R`, `slurm_pilot.sh`.
- Notebook: `notebooks/pilot_report.qmd`, rendered with
  `/home/mdmanurung/.local/bin/quarto` (the env's quarto is broken).

## Deviations from the plan

- **Parameterisation per scenario**, from the `centered` column of
  `scenarios.tsv`: null and dependence-only were fitted noncentered, all others
  centered. Centered fits of null and dependence-only funnelled in the smoke
  stage (R-hat up to 1.4, ESS < 30). Those smoke replicates are kept in
  `../2026-11-02_pilot-smoke/reps_centered_funnel/` and are not in any metric.
- **Simulation** used `donor_visit = TRUE`. A separate 5-seed cell
  (`positive_margins_no_dv_J3_D24`) simulates without donor-by-visit variation
  and fits with it; all 5 replicates mixed badly (R-hat 1.23, ESS 14).
- **`03_run-replicate.R` was edited while smoke jobs ran**, truncating 10 smoke
  tasks (positive margins and visit change, 100 donors). They were rerun.
- **`scenarios.tsv` in the smoke run was replaced** by the final version before
  those reruns.
- **Not run** (array limit and the 400-task budget): J = 4 at 100 donors; J = 4
  at 24 donors for null, negative margins and visit change.
- Flagged replicates (R-hat > 1.01, ESS < 400 or any divergence) are kept in
  every metric; counts are in `metrics_diagnostics.tsv`. Treedepth saturation
  was not a flag criterion (dependence-only J = 4 saturated in 29 of 30).
- COMPASS was run on seeds 1-5 only (about 125 s per run).

## Corrected summary (2026-11-02, after code review)

The metric tables in this directory were made by the first `04_summarise-pilot.R`,
which read every replicate directory. That leaked the 5-seed diagnostic cell
`positive_margins_no_dv_J3_D24` into `metrics_scores.tsv`. The corrected tables
are in `../2026-11-02_pilot-summary-v2/` (same replicates, no refits). Every
other scenario x score is unchanged; prefer the v2 tables.
