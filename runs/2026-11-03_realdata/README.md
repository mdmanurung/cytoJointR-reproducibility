# Real-data report: RV144, HVTN078, TB and SARS-CoV-2 (2026-11-03)

COMPASS ranks donors better than the default cytoJointR model (`m1dm`, sample level) on RV144, TB_s and TB_u.
The default P(RMS > eps) rule calls 8 of 36 RV144 placebos positive; COMPASS calls none.
HVTN078 and SARS-CoV-2 cannot test this. HVTN078 has 5 placebos, and the SARS-CoV-2 donors had prior COVID-19.
This page only reports tables that exist in this repo. It runs no new analysis. Each table names its source file.

Run directories (all under `runs/2026-11-03_realdata/` unless stated):

| Directory | Holds | Detail |
|---|---|---|
| `vaccine/rv144`, `vaccine/hvtn078` | RV144 and HVTN078 fits and comparisons | this page |
| `tb/` | TB_s, TB_u, SimpleCOMPASS | [tb/README.md](tb/README.md) |
| `sarscov2/` | SARS-CoV-2 CD4 and CD8 | [sarscov2/README.md](sarscov2/README.md) |
| `vaccine-dev/` | Earlier partial vaccine run (jobs 25941161 to 25941164). This page does not use it. | none |
| `../2026-11-03_estimand/` | Sample level against structured level, and mean misfit | this page |

Notebooks that draw the figures: `notebooks/realdata_vaccine.qmd`, `notebooks/realdata_tb.qmd`, `notebooks/realdata_sarscov2.qmd` (rendered `.html` files sit beside them).

## Datasets and provenance

| Data set | Source | Identifier and licence | Checksum |
|---|---|---|---|
| RV144, HVTN078, TB_s, TB_u | COMPASSDataPackage 1.0.0, Lin et al. 2015, doi:10.1038/nbt.3187, file `COMPASSData.rda` | Zenodo doi:10.5281/zenodo.17500. The record says `license: other-open`; the GitHub repo ships an Apache 2.0 `LICENSE` file. | zip md5 `6a82ee0c99b3ba326c0e62b4166a98db` (record file checksum; script `scripts/20_prepare-rv144.R`) |
| SimpleCOMPASS TB | `SimpleCOMPASSExamples.csv` in the COMPASS 1.48.0 package, Rakshit et al. 2017, doi:10.1038/s41598-017-10773-5 | No Zenodo ID recorded. No licence recorded here. | none recorded |
| SARS-CoV-2 Spike ICS | Ford et al. 2023 | Zenodo doi:10.5281/zenodo.8088178, CC BY 4.0, published 2023-06-27 | CD4 `975e31de20dd715a15e02062e4e2b312`, CD8 `da0fd5e8ebee2194e6aae54cd2ed86e1`; both match Zenodo |

Sources of the facts above: `data/raw/compass_data_package/files.json`, `data/raw/zenodo_8088178/record.json`, `data/processed/sarscov2/{citation.txt,checksums.tsv}` and the headers of `scripts/20` to `scripts/24`.
`data/` is not committed.

What each table holds (from `data/processed/*/prep_log.tsv`, `data/processed/sarscov2/data_checks.tsv`):

| Data set | Donors | Contrast | Markers kept | Note |
|---|---|---|---|---|
| RV144 | 262: 226 vaccine (38 infected), 36 placebo | vaccine against placebo | IFNg, IL2, TNFa, CD154, IL4 (IL17a dropped) | COMPASS kept 24 of 64 patterns |
| HVTN078 | 76: 71 vaccine, 5 placebo | vaccine against placebo | IFNg, IL2, TNFa, CD154, CD107a (IL4, MIP1B dropped) | COMPASS kept 26 of 128 patterns |
| TB_s, TB_u | 40: QFT+ 22, QFT- 18 | IGRA status, not active TB | IFNg, IL2, TNFa, CD154, IL17a (IL22, IL4 dropped) | COMPASS kept 22 (TB_s) and 19 (TB_u) of 128 patterns |
| SimpleCOMPASS | 21: EPTB 11, IGRA+ 10 | EPTB against IGRA+ | IFNg, IL2, TNF, IL17A, IL22 | complete 256-pattern partition |
| SARS-CoV-2 | CD4 and CD8, 105 samples each; visits E01 (n 7), E02 (n 14), E03 (n 14) | Spike and PHA against DMSO | IFNg, IL2, TNFa, CD40L | complete 16-pattern partition |

COMPASS filtered the RV144, HVTN078, TB_s and TB_u tables.
Cells in removed patterns sit in the all-negative count, so positive counts are lower bounds. Every method sees the same tables.

## Methods

- **Package:** frozen worktree `cytoJointR-wt/redesign-run`, commit 49278c9. Seed 20261103. R 4.5.1.
- **cytoJointR models:** `m1dm` (default; M1 mean structure plus a Dirichlet-multinomial sample residual), `m1` and `C0` (conjugate Dirichlet).
- **Scores:** RMS, CPS and CCS. CCS is experimental and no statement here rests on it.
- **cytoJointR rule:** a donor is positive when P(RMS > 1e-4) > 0.95, at the sample level.
- **COMPASS 1.48.0:** FS, PFS and Response. Rule: Response (degree 1) > 0.5.
- **polyICSFlow 0.99.4** (R 4.4.3): net sum over reported combinations. Rule: clipped net sum > 0.
- **AUC:** vaccine against placebo uses all 262 donors, not the benchmark dev half. CIs are donor bootstraps within group, 2000 draws.
- **Infection association:** logistic model among vaccinees with gender, age and behavioural risk (model M1), odds ratio per SD. Holm correction over COMPASS FS and PFS plus `m1dm` RMS and CPS.
- **SLURM jobs** (from log names): RV144 fit 25941346, COMPASS 25941342, polyICSFlow 25941343; HVTN078 25941262, 25941264, 25941266; TB 25941339 and 25941340; SARS-CoV-2 25942136; estimand 25942625 and 25942631 to 25942633.
- **Scripts:** `scripts/20` to `24` prepare data, `30_analyse-vaccine.R`, `31_analyse-tb.R`, `32_analyse-sarscov2.R` fit and score, `60_diagnose-estimand.R` runs the estimand test.

## RV144 and HVTN078

AUC of vaccine against placebo (source: `vaccine/rv144/auc.tsv`, `vaccine/hvtn078/auc.tsv`).
RV144 has 226 vaccinees and 36 placebos. HVTN078 has 71 vaccinees and 5 placebos.

| Score | RV144 AUC (95% CI) | HVTN078 AUC (95% CI) |
|---|---|---|
| m1dm P(RMS>eps), sample level (default) | 0.755 (0.672 to 0.828) | 0.428 (0.261 to 0.690) |
| m1dm RMS mean | 0.732 (0.638 to 0.816) | 0.386 (0.135 to 0.670) |
| m1dm CPS mean | 0.921 (0.875 to 0.958) | 0.330 (0.149 to 0.507) |
| m1 P(RMS>eps) | 0.793 (0.719 to 0.854) | 0.487 (0.338 to 0.724) |
| COMPASS FS | 0.977 (0.959 to 0.991) | 0.397 (0.124 to 0.685) |
| COMPASS PFS | 0.972 (0.952 to 0.989) | 0.392 (0.124 to 0.676) |
| polyICSFlow clipped net sum | 0.834 (0.734 to 0.920) | 0.400 (0.144 to 0.690) |

The paired AUC difference between the `m1dm` sample-level P(RMS > eps) and COMPASS Response is -0.221 (95% CI -0.303 to -0.149) on RV144 (`vaccine/rv144/auc_diff.tsv`).

False positives among placebos, each method at its own rule (source: `vaccine/rv144/fpr.tsv`, `vaccine/hvtn078/fpr.tsv`):

| Data set | Method | Rule | Placebo false positives | FPR 95% CI | Vaccinees called |
|---|---|---|---|---|---|
| RV144 | m1dm | P(RMS > eps) > 0.95 | 8/36 | 0.101 to 0.392 | 157/226 |
| RV144 | COMPASS | Response (degree 1) > 0.5 | 0/36 | 0.000 to 0.097 | 42/226 |
| RV144 | polyICSFlow_clipped | clipped net > 0 | 36/36 | 0.903 to 1.000 | 226/226 |
| HVTN078 | m1dm | P(RMS > eps) > 0.95 | 4/5 | 0.284 to 0.995 | 56/71 |
| HVTN078 | COMPASS | Response (degree 1) > 0.5 | 3/5 | 0.147 to 0.947 | 24/71 |
| HVTN078 | polyICSFlow_clipped | clipped net > 0 | 5/5 | 0.478 to 1.000 | 71/71 |

RV144 infection association among vaccinees, model M1 (source: `vaccine/rv144/infection.tsv`; the Holm family is these four rows):

| Score | OR per SD (95% CI) | p | Holm p | Events |
|---|---|---|---|---|
| m1dm_RMS | 1.281 (0.884 to 1.856) | 0.191 | 0.574 | 38 |
| m1dm_CPS | 1.099 (0.765 to 1.579) | 0.609 | 0.609 | 38 |
| COMPASS_FS | 0.779 (0.534 to 1.137) | 0.195 | 0.574 | 38 |
| COMPASS_PFS | 0.639 (0.431 to 0.948) | 0.026 | 0.104 | 38 |

No score passes Holm correction. COMPASS PFS has the only nominal p below 0.05, and its odds ratio points the opposite way from `m1dm` RMS.

## TB

AUC of QFT+ against QFT- (TB_s, TB_u) and EPTB against IGRA+ (SimpleCOMPASS).
Source: `tb/tb_s/auc.tsv`, `tb/tb_u/auc.tsv`, `tb/tb_simplecompass/auc.tsv`. Data facts and deviations: [tb/README.md](tb/README.md).

| Score | TB_s | TB_u | SimpleCOMPASS |
|---|---|---|---|
| m1dm P(RMS>eps), sample level (default) | 0.758 (0.593 to 0.896) | 0.649 (0.467 to 0.821) | 0.550 (0.500 to 0.650) |
| m1dm RMS mean | 0.705 (0.533 to 0.856) | 0.604 (0.422 to 0.780) | 0.845 (0.627 to 1.000) |
| m1dm CPS mean | 0.841 (0.697 to 0.947) | 0.725 (0.543 to 0.884) | 0.373 (0.127 to 0.636) |
| m1 P(RMS>eps) | 0.876 (0.745 to 0.967) | 0.755 (0.588 to 0.904) | 0.550 (0.500 to 0.650) |
| COMPASS FS | 0.987 (0.954 to 1.000) | 0.864 (0.730 to 0.965) | 0.673 (0.418 to 0.900) |
| COMPASS PFS | 0.955 (0.884 to 0.997) | 0.864 (0.730 to 0.965) | 0.564 (0.300 to 0.809) |
| polyICSFlow clipped net sum | 0.717 (0.556 to 0.861) | 0.631 (0.462 to 0.798) | 0.873 (0.664 to 1.000) |

COMPASS FS is higher than the `m1dm` default on TB_s and TB_u.
On SimpleCOMPASS (21 donors) the intervals overlap everywhere, and polyICSFlow clipped has the highest point estimate.
`runs/2026-11-04_donor-effects/tb_u_diagnosis.md` explains the TB_u gap. The model-free net any-positive difference gives AUC 0.58. The net fraction of cells with 3 or more markers gives 0.88. The group signal sits in polyfunctional patterns, not in the total.

## SARS-CoV-2

Source: `sarscov2/CD4/` and `sarscov2/CD8/`; the table is rendered in `notebooks/realdata_sarscov2.html` (fraction of Spike donor-visits called, pooled over E01 to E03, at each method's own rule).

| Population | m1dm | m1 | C0 | COMPASS | polyICSFlow clipped | polyICSFlow unclipped |
|---|---|---|---|---|---|---|
| CD4 (35 donor-visits) | 1.000 | 1.000 | 1.000 | 0.857 | 1.000 | 1.000 |
| CD8 (35 donor-visits; COMPASS 28) | 0.714 | 0.800 | 0.600 | 0.750 | 1.000 | 0.914 |

Per visit, the `m1dm` call at eps 1e-4 and at eps 1e-3 against COMPASS Response > 0.5.
These counts come from rows of `scores.tsv` (model m1dm, score RMS, `p_positive` > 0.95, `eps_set` main and sens) and `compass.tsv`, counted by hand from the stored tables.

| Population | Visit | Donors | m1dm, eps 1e-4 | m1dm, eps 1e-3 | COMPASS |
|---|---|---|---|---|---|
| CD4 | E01 | 7 | 7 | 6 | 4 |
| CD4 | E02 | 14 | 14 | 14 | 13 |
| CD4 | E03 | 14 | 14 | 13 | 13 |
| CD8 | E01 | 7 | 1 | 0 | no call (COMPASS failed) |
| CD8 | E02 | 14 | 11 | 9 | 10 |
| CD8 | E03 | 14 | 13 | 10 | 11 |

All PHA donor-visits (35 per population) are called by `m1dm` and by COMPASS.
COMPASS failed for CD8 Spike at E01 ("Not a matrix.", one category left after its filter), so 7 donor-visits have no COMPASS Spike call ([sarscov2/README.md](sarscov2/README.md)).
A CD8 parent can have only 398 cells (`data_checks.tsv`), so eps 1e-4 is below one cell for some samples.
CCS is not used for any conclusion.

## Sample level against structured level

`60_diagnose-estimand.R` refits `m1dm` and scores two estimands. The sample level is the realised sample probability. The structured level is the model's fitted mean.
Source: `../2026-11-03_estimand/{rv144,tb_s,tb_u}/performance.tsv`.

| Data set | Level | AUC of RMS mean | AUC of P(RMS>1e-4) | False positives at P>0.95 | TPR at P>0.95 |
|---|---|---|---|---|---|
| RV144 | sample | 0.732 | 0.755 | 8/36 | 0.695 |
| RV144 | structured | 0.948 | 0.932 | 5/36 | 0.881 |
| TB_s | sample | 0.705 | 0.758 | 3/18 | 0.409 |
| TB_s | structured | 0.919 | 0.917 | 0/18 | 0.091 |
| TB_u | sample | 0.604 | 0.649 | 1/18 | 0.136 |
| TB_u | structured | 0.662 | 0.610 | 0/18 | 0.000 |

COMPASS FS in the same files: RV144 0.977, TB_s 0.987, TB_u 0.864.

- On RV144 and TB_s the structured level ranks much better than the sample level.
- At the P > 0.95 rule the structured level calls few TB donors: 2 of 22 QFT+ on TB_s (TPR 0.091) and 0 of 22 on TB_u, with 0 false positives. Its AUC does not make it a usable rule there.
- On TB_u it does not: P(RMS > 1e-4) falls from 0.649 to 0.610, and only the RMS mean rises (0.604 to 0.662).
- The structured level still calls 5 of 36 RV144 placebos at 1e-4. At eps 5e-4 and 1e-3 it calls none (TPR 0.735 and 0.513; `rv144/performance.tsv`).
- Its placebo RMS is not near zero. The median structured RMS mean is 0.00126 in placebos and 0.00306 in vaccinees; the sample-level placebo median is 3.7e-5 (`rv144/level_summary.tsv`).
- Its infection odds ratio per SD is 0.997 (95% CI 0.700 to 1.420) for RMS (`rv144/infection.tsv`).
- A model-free stimulated-minus-control z test also calls 10 of 36 placebos at z > 1.645 (`rv144/performance.tsv`, level raw). The placebo excess is not unique to the model.

## Mean misfit

The structured mean over-predicts observed rates in the RV144 fit and in SimpleCOMPASS.
It is close for any-positive cells in TB_s (ratios 0.99 to 1.12) and TB_u (1.03 to 1.09; not in the table). It over-predicts cells with 3 or more positive markers almost everywhere.
Column `fitted_over_observed` of `../2026-11-03_estimand/<data set>/misfit_summary.tsv`; a ratio of 1 is a perfect fit:

| Data set | Group | Condition | Statistic | Fitted / observed |
|---|---|---|---|---|
| rv144 | placebo | control | any_positive | 1.26 |
| rv144 | vaccine | control | any_positive | 1.18 |
| rv144 | placebo | stimulated | any_positive | 1.61 |
| rv144 | vaccine | stimulated | any_positive | 1.48 |
| rv144 | placebo | control | degree3plus | 8.25 |
| rv144 | vaccine | control | degree3plus | 7.70 |
| rv144 | placebo | stimulated | degree3plus | 24.20 |
| rv144 | vaccine | stimulated | degree3plus | 4.46 |
| tb_s | QFT+ | control | degree3plus | 2.68 |
| tb_s | QFT- | control | degree3plus | 1.81 |
| tb_s | QFT+ | stimulated | degree3plus | 0.94 |
| tb_s | QFT- | stimulated | degree3plus | 4.24 |
| tb_u | QFT+ | control | degree3plus | 2.86 |
| tb_u | QFT- | control | degree3plus | 2.03 |
| tb_u | QFT+ | stimulated | degree3plus | 1.37 |
| tb_u | QFT- | stimulated | degree3plus | 3.58 |
| tb_simplecompass | EPTB | control | degree3plus | 29.22 |
| tb_simplecompass | IGRA+ | control | degree3plus | 23.91 |
| tb_simplecompass | EPTB | stimulated | degree3plus | 6.93 |
| tb_simplecompass | IGRA+ | stimulated | degree3plus | 3.32 |

The fitted mean is pairwise log-linear. Whether higher-order terms would remove the misfit was not tested.

## Model checks

Source: `phi_summary.tsv`, `loo.tsv`, `loo_compare.tsv`, `diagnostics.tsv` in each run directory (sample-level LOO for the vaccine data).
A positive elpd difference would favour `m1`. Every difference is negative, so LOO favours `m1dm`.

| Data set | phi (mean) | Variance inflation at median N | Max R-hat | Min bulk ESS | Divergences | elpd diff, m1 minus m1dm (SE) | Pareto k > 0.7, m1dm |
|---|---|---|---|---|---|---|---|
| RV144 | 817 | 68.1 | 1.010 | 721 | 0 | -57246 (5265) | 154/524 |
| HVTN078 | 8861 | 13.6 | 1.003 | 456 | 0 | -11580 (1258) | 54/152 |
| TB_s | 4368 | n/a | 1.003 | 868 | 0 | -2199 (343) | 40/80 |
| TB_u | 4706 | n/a | 1.023 | 329 | 0 | -2171 (262) | 14/80 |
| SimpleCOMPASS | 117 | n/a | 1.007 | 804 | 0 | -55666 (13545) | 13/42 |
| SARS-CoV-2 CD4 | 164 | n/a | 1.017 | 444 | 0 | -16653 (3390) | 23/105 |
| SARS-CoV-2 CD8 | 310 | n/a | 1.008 | 563 | 0 | -3557 (557) | 34/105 |

The Pareto k counts are high, so the LOO estimates are unreliable. The ranking of `m1dm` over `m1` is large against its SE in every row.
Variance inflation is the overdispersion `m1dm` adds over a multinomial at the median parent count.

## Limitations

- **HVTN078 AUC is uninformative.** There are 5 placebos. Every method scores about 0.4, and the m1dm CI is 0.261 to 0.690. The placebo FPR CI is 0.284 to 0.995 for `m1dm`.
- **SARS-CoV-2 cannot give a false-positive rate.** The record title says "persons recovered from SARS-CoV-2 infection" and the paper title says "persons with previous COVID-19". Evidence is titles and abstract; the enrolment criteria were not read (`benchmark/README.md`). CD4 E01 Spike is positive in 7 of 7 donors for `m1dm`, which fits prior infection.
- **COMPASS beats the cytoJointR default on real-data classification.** This holds on RV144, TB_s and TB_u. It is not shown on SimpleCOMPASS.
- **The default rule over-calls placebos.** The sample-level rule gives 8/36 placebo false positives on RV144 (COMPASS 0/36).
- **The structured level ranks better, with exceptions.** It does so on RV144 and TB_s, not on TB_u. It has no infection signal, and its placebo RMS is above zero.
- **Mean misfit.** See the table above. Fitted rates exceed observed rates, mostly for cells with 3 or more markers.
- **Resolution of eps.** At the smallest RV144 parent count (5141 cells) eps = 1e-4 is 0.51 of one cell; at the median (54905) it is 5.5 cells (`vaccine/rv144/epsilon_check.tsv`).
- **Small data.** TB has 40 donors and SimpleCOMPASS 21. The TB_u AUC SE is about 0.07 (`tb_u_diagnosis.md`). The TB contrast is IGRA status, not active TB.
- **Filtered tables.** Positive counts in the COMPASS-filtered tables are lower bounds.
- **Infection association.** There are 38 events. No score survives Holm correction.
- **No external hold-out.** Each data set is one trial or cohort, and the TB sets were also used for model development.
- **Experimental model.** A responder-mixture model on branch `feat/donor-effects` is close to COMPASS on RV144 dev. It is not merged and does not beat COMPASS on the frozen criteria. See the iteration log in [benchmark/README.md](../../benchmark/README.md).
- **CCS is experimental.** No conclusion here rests on it.

## Rerun

Never overwrite a run directory. Commands for TB are in [tb/README.md](tb/README.md) and for SARS-CoV-2 in [sarscov2/README.md](sarscov2/README.md).
Render a notebook with `quarto render notebooks/realdata_<name>.qmd` from the repo root.
