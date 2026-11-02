# cytoJointR frozen benchmark

The benchmark scores a cytoJointR candidate against COMPASS and polyICSFlow on fixed data, fixed splits and fixed rules.
`definition.tsv` and `rv144_split.tsv` were committed before any candidate was scored. Do not edit them.
A change to the definition starts a new benchmark version, with a new directory and a new baseline.

Files:

- `definition.tsv`: one row per metric. Columns: `metric_id`, `dataset`, `role` (dev, test, guard, win, excluded), `tier`, `metric`, `score`, `level`, `direction`, comparator value and 95% bounds for COMPASS and polyICSFlow, `rule`, `note`. `level = fit` means the level of the candidate's own `cj_fit(level = )`.
- `rv144_split.tsv`: donor, arm, infection status, stratum, `split` (dev or test).
- `sim_scenarios.tsv`: the solved simulation scenarios of the guards.
- `leaderboard.tsv`: one row per candidate and metric, appended by `scripts/70_run-scorecard.R`.
- `test_access.tsv`: one row per test-set access.
- Built by `scripts/69_define-scorecard.R`, `scripts/helper06_scorecard.R`. Run by `scripts/70_run-scorecard.R`, `scripts/slurm_scorecard.sh`.

## Fixed rule

epsilon = 1e-4. A donor is called positive when P(RMS > epsilon) > 0.95. Nobody changes this during iteration.
Comparator rules: COMPASS Response (degree 1) > 0.5, polyICSFlow clipped net sum > 0.

## Datasets

| Dataset | Role | Use |
|---|---|---|
| RV144, dev half | dev | Primary detection metric and its companions. 131 donors: 113 vaccinees, 18 placebos. |
| RV144, test half | test | Test-only metrics. 131 donors: 113 vaccinees, 18 placebos. Seeded 50/50 split (seed 20261104), stratified by placebo, vaccine-infected, vaccine-non-infected. |
| TB_s, TB_u | dev | Secondary AUC. QFT+ 22 against QFT- 18. Whole data sets. |
| Simulation, generators A, B, C, D | guard | Coverage, detection and magnitude. Reduced design below. |
| HVTN078 | excluded | 5 placebos and every method scores about 0.4 AUC. No detection metric. |
| SARS-CoV-2 (Ford et al.) | unusable for FPR | See below. |

RV144 is fitted once, on all 262 donors. The model never sees arm or outcome, so the fit is unsupervised with respect to the label.
Metrics then use dev donors only (test donors only in `--test` mode).
Fitting on all donors lets test-donor data shape the shared hyperparameters (phi, tau). This costs a little independence, and every candidate pays it the same way.
The scores of test donors sit in the run directory. Do not open them before `--test`.

### SARS-CoV-2 naivety

Not usable for false-positive rates. The paper is Ford et al., Nat Immunol 2024, PMID 38057617, doi:10.1038/s41590-023-01692-x.
Its title reads "Repeated mRNA vaccination sequentially boosts SARS-CoV-2-specific CD8+ T cells in persons with previous COVID-19".
The Zenodo record 8088178 (the data we use) says "persons recovered from SARS-CoV-2 infection".
Participants were therefore not naive at enrolment, and visit E01 Spike responses are expected, not false positives.
That `m1dm` called 100% of CD4 E01 donors responders is then consistent with biology.
Evidence level: titles and abstract summary only. PMC and Nature blocked the full text, so the enrolment criteria were not read.

## Metrics

Primary: RV144-dev AUC of the detection score P(RMS > epsilon), at the candidate's fit level.
Companions on the same score: partial AUC for FPR <= 0.1 (mean TPR over that range, chance is 0.05), TPR at FPR <= 0.05 and at FPR <= 0.10.
COMPASS and our model sit at different ROC corners, so one operating point misleads. With 18 dev placebos, FPR moves in steps of 0.056, so FPR <= 0.05 means no false positive.

Secondary (all on dev): AUC for every score family (RMS mean, P(RMS), CPS mean, P(CPS), CCS mean and P(CCS), the CCS pair marked experimental), at the sample and the structured level, on RV144 dev, TB_s and TB_u.
Also placebo FPR and vaccinee TPR at the fixed rule on RV144 dev.

Test-only (score once, at the end): RV144-test AUC, pAUC, TPR at matched FPR, FPR and TPR at the rule, and the infection association.
Infection association: odds ratio per SD of RMS and CPS among test-half vaccinees (about 113 donors, 19 events), logistic model with gender, age and risk (model M1 of `30_analyse-vaccine.R`).
The COMPASS PFS reference is OR 0.64, p 0.026 on all 226 vaccinees. With 19 events the test-half test is underpowered; the OR direction is the main read.

CIs are donor bootstraps (2000 draws, seed 20261103), resampling within case and within non-case. Differences use the same draws for both methods.

### Simulation guards (reduced design)

Scenarios from the 2026-11-02 redesign runs, solved parameters copied to `sim_scenarios.tsv`:
`rare_weak_J3`, `rare_weak_J4`, `moderate_J3` (A), `generator_B_J3` (B), `mixture_J3` (C), `D1_rare_weak_J3`, `D2_rare_weak_J3` (D).
Seeds 1 to 10 of the original 30, 24 donors each, 70 replicates per candidate. The cached COMPASS and polyICSFlow runs cover the same seeds, so comparators are not rerun.
The original design has 18 scenarios and 30 seeds, plus D. We cut it so one scorecard fits in about 1 h of SLURM across 7 array tasks.
The replicate uses the data seed, simulator and truth of `06_run-benchmark-replicate.R` (`helper02_benchmark-sim.R`, `helper01_true-scores.R`) with the same scoring call; it fits only the candidate and does not rerun COMPASS or polyICSFlow.

- Coverage guard: 90% interval coverage of donor RMS and CPS against the truth at the fit level (sample truth for `level = "sample"`, structured truth for `"structured"`). Pass rule is rule 1 of `11_apply-decision-rule.R`: coverage + 2 SE >= 0.85 and coverage - 2 SE <= 0.97, SE clustered by seed. No upper bound in null scenarios (none are in this reduced set). The false-sign and ipf rules are not guarded.
- Detection guard: AUC of P(RMS > epsilon) for responders against non-responders (true sample RMS > 1e-4), pooled over rare_weak J3 and J4 and over seeds, against COMPASS_default FS and polyICSFlow clipped. CIs resample seeds.
- Magnitude guard: Spearman of the donor RMS estimate with the true RMS, same pooling, against polyICSFlow unclipped.
- With 10 seeds, a coverage SE is about 0.01 to 0.03. A single scorecard cannot detect small coverage drift.

## Win criteria

A candidate "beats COMPASS" on dev when all hold:

1. `rv144_dev.auc.P_RMS.fit`: paired-bootstrap difference against COMPASS FS has a 95% lower bound above 0.
2. `rv144_dev.pauc10` difference above 0, and `rv144_dev.tpr_fpr10` difference at least 0 (point estimates).
3. `sim.auc.rare_weak` difference against COMPASS_default FS at least 0 (point estimate).
4. Every `sim.cov90` guard passes.

The claim stands only after the same test as (1) passes on RV144 test, run once. COMPASS FS has dev AUC 0.974, so (1) needs a candidate AUC near 0.99.

A candidate "beats polyICSFlow" on dev when all hold:

1. `rv144_dev.auc.P_RMS.fit`: difference against polyICSFlow clipped has a lower bound above 0. This is trivially beatable: clipped polyICSFlow calls 36 of 36 placebos.
2. `sim.spearman.rare_weak.RMS_mean.fit`: difference against polyICSFlow unclipped has a lower bound above 0. Magnitude ranking is the meaningful comparison.
3. `sim.auc.rare_weak` difference against polyICSFlow clipped at least 0 (point estimate).
4. Every `sim.cov90` guard passes.

Again confirmed on RV144 test, run once. `scripts/70_run-scorecard.R` writes both verdicts as `win.beats_compass` and `win.beats_polyics` rows (value 1 or 0, `note` lists the failed parts).

## Running

```
sbatch scripts/slurm_scorecard.sh <bench|test|collect> ...   # see the header of the script
```

`70_run-scorecard.R` stages: `rv144`, `tb`, `sim <scenario>` (each fits and writes raw scores), `collect` (reads raw scores, computes the dev metrics, appends to the leaderboard), `test` (the `--test` mode).
Candidate = package worktree path, `cj_fit` arguments such as `model=m1dm,level=structured`, a label and a run directory.
Every leaderboard row carries the package SHA, the candidate label, the git SHA of `definition.tsv` and a timestamp.

## Test-set access counter

Test-set accesses so far: **0**.

`--test` appends a row to `test_access.tsv`, prints the new count and refuses a second test run for the same candidate label.
Update this number with each access. Record why the access was made.

## Baseline (2026-11-03)

Package e54b883 (main), three candidates: `baseline-m1dm-sample` (default), `baseline-m1dm-structured`, `baseline-m1`.
Rows are in `leaderboard.tsv`, raw scores in `runs/2026-11-03_bench-baseline/`. No candidate beats COMPASS or polyICSFlow on the dev win criteria; test mode was not run.

## Known limits

- RV144 placebos are the only negatives. Placebo infection status is not recorded, so the split stratifies placebo as one group.
- Dev and test come from one trial. The test half checks selection on the dev half, not transport to another trial.
- TB_s and TB_u are used for model development too, so they have no hold-out. Treat them as secondary.
- Cached comparators were run on the whole data set. The dev-half comparator values are recomputed from the cached per-donor scores, not refitted.
