# Model redesign run: provenance (2026-11-02)

This run fitted three candidate models and two comparators to 545 simulated replicates.
All 545 replicates finished. The decision rule and the user's choice of `m1dm` are in [the decision record](../2026-11-02_redesign-decision/README.md).
This page records what was run and where the outputs sit. It reports no new results.

## Design

- 545 replicates = 18 scenario cells x 30 seeds (540) + `runtime_100_J4` x 5 seeds (`jobs.tsv`).
- The 18 cells are 9 scenarios at J = 3 and J = 4: `null`, `rare_weak`, `rare_moderate`, `moderate`, `mixture`, `negative`, `generator_B`, `missing_control`, `visit_change` (`scenarios.tsv`).
- Generators: A (most scenarios), B (`generator_B`), C (`mixture`). Each replicate has 24 donors (runtime scenario excepted; see `scenarios.tsv`).
- Candidates: C0 (conjugate Dirichlet), C1 (`m1`), C2 (`m1dm`). C2phi was never run: `model_status.tsv` in the summary shows 545 `not_run`.
- Comparators: COMPASS (default filter and no filter) and polyICSFlow.
- `visit_change` replicates (60) have no polyICSFlow output by design. `07_run-polyicsflow.R` exits for two-visit scenarios, so `done_b` is FALSE for them in `replicate_status.tsv`. COMPASS status is `not_run` for the same 60.

## Frozen package

| Item | Value |
|---|---|
| Worktree | `cytoJointR-wt/redesign-run`, detached at 49278c9 (`docs: NEWS for m1dm, C0 and the side fixes`) |
| Branch lineage | 49278c9 is on `feat/model-redesign` |
| Commit recorded in all 545 `reps/*/seed_*/session_info.txt` | 49278c9614959063b57da55304b1217ae67dbd7e |
| Commit recorded in the run-level `session_info.txt` | db8d20956e5a8408bf39b63d9f37e448f03bc1d5 |

db8d209 is an ancestor of 49278c9. The run-level file was probably written when the run directory was set up, before the worktree moved to 49278c9. That is inferred, not recorded.
The difference between the two commits is `m1dm`, C0, a COMPASS adapter fix and a `cj_describe` fix (five commits, `git log db8d209..49278c9`).
The replicate files and the summary (`session_info_summary.txt`) record 49278c9. Treat 49278c9 as the package that produced the results.

## Scripts

| Step | Script | Environment |
|---|---|---|
| Define scenarios, solve effect sizes, write `jobs.tsv` | `scripts/05_define-benchmark.R` | R 4.5.1 (conda `R4_51`) |
| One replicate: simulate, fit C0 to C2, run COMPASS | `scripts/06_run-benchmark-replicate.R` | R 4.5.1 |
| polyICSFlow for the same counts | `scripts/07_run-polyicsflow.R` | R 4.4.3 (`batchdx-cyto`) |
| Summaries | `scripts/08_summarise-benchmark.R` | R 4.5.1 |
| Pre-registered decision rule | `scripts/11_apply-decision-rule.R` | R 4.5.1 |
| Simulators and true scores | `helper01_true-scores.R`, `helper02_benchmark-sim.R`, `lib_generator-b.R`, `lib_generator-c.R` | |

Held-out generator D uses `lib_generator-d.R` and `scripts/10_define-heldout.R`. It was committed before any candidate ran (fed35d9).

## SLURM

- Wrapper: `scripts/slurm_benchmark.sh`. Job 25939505, array 1 to 95 (95 log files in `logs/`).
- Partition `short`, 4 CPUs, 24G, 1 h per task. Each task runs several replicates and starts no new one after 2700 s.
- A finished step is skipped, so a failed task can be resubmitted.
- Stan was compiled once before submission, so tasks did not race to build it.
- Task `i` runs the rows of `jobs.tsv` with `task == i`.

Command pattern (from the wrapper header):

    sbatch --array=1-95 scripts/slurm_benchmark.sh runs/2026-11-02_redesign

## Outputs

Directory `runs/2026-11-02_redesign/`:

| File | Holds |
|---|---|
| `scenarios.tsv` | Solved parameters of the 19 scenarios |
| `jobs.tsv` | Task, scenario and seed of the 545 replicates |
| `tuning_check.tsv` | Achieved true-RMS distribution per scenario |
| `reps/<scenario>/seed_NN/` | Per-replicate counts, truths, model scores, diagnostics, runtime, `polyicsflow/` |
| `logs/` | One log per array task |
| `session_info.txt` | Run-level session (see the commit note above) |

Related directories:

| Directory | Holds |
|---|---|
| [`../2026-11-02_redesign-summary/`](../2026-11-02_redesign-summary/) | Metric tables for this run (`calibration.tsv`, `false_sign.tsv`, `auc.tsv`, `replicate_status.tsv`, `model_status.tsv` and others), from `08_summarise-benchmark.R` |
| [`../2026-11-02_redesign-decision/`](../2026-11-02_redesign-decision/README.md) | Decision rule output (`decision.tsv`, `decision_cells.tsv`, `decision.md`) and the decision record |
| [`../2026-11-02_redesign-interim/`](../2026-11-02_redesign-interim/) | Summary of a partial state: 482 replicates complete in both steps (`replicate_status.tsv`) |
| [`../2026-11-02_redesign-smoke/`](../2026-11-02_redesign-smoke/) and `-smoke-summary` | Smoke test: 37 replicates, job 25939399 |
| [`../2026-11-02_redesign-vccheck/`](../2026-11-02_redesign-vccheck/) | Single-scenario `visit_change_J3` rerun; no notes were recorded |
| [`../2026-11-02_redesign-D/`](../2026-11-02_redesign-D/README.md) | Generator D, first attempt. It failed before any fit; no results. |
| [`../2026-11-02_redesign-D2/`](../2026-11-02_redesign-D2/) and `-D2-summary` | Generator D rerun: 240 replicates, job 25940215. Opened after the decision. |
| `notebooks/redesign_report.qmd` | Report notebook |

## What the rule and the decision said

From `../2026-11-02_redesign-decision/decision.tsv` (sample target, 18 A/B/C cells):

| Model | Cells passed of 18 | Mean RMSE |
|---|---|---|
| C0 | 5 | 6.39e-4 |
| C1 (`m1`) | 9 | 4.42e-4 |
| C2 (`m1dm`) | 6 | 4.42e-4 |

The rule returned no default. The user then chose C2 (`m1dm`) as a post hoc decision. The reasons, the pre-registration trail and the corrections are in the decision record.
