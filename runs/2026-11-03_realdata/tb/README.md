# TB real-data benchmark (cytoJointR, COMPASS, polyICSFlow)

Package: frozen worktree `cytoJointR-wt/redesign-run`, commit 49278c9. Seed 20261103. Never overwrite this directory.

## How it was run

    Rscript scripts/22_prepare-tb-compass.R            # data/processed/tb_s, tb_u
    Rscript scripts/23_prepare-tb-simplecompass.R      # data/processed/tb_simplecompass
    sbatch --array=1-3 scripts/slurm_realdata_tb.sh runs/2026-11-03_realdata/tb fit    # job 25941339
    sbatch --array=1-3 --dependency=afterok:25941339 scripts/slurm_realdata_tb.sh runs/2026-11-03_realdata/tb post   # job 25941340
    quarto render notebooks/realdata_tb.qmd

Each `<dataset>/` holds donor scores, model status, diagnostics, phi, LOO, PPC, runtime and AUC tables
(column lists are in the header of `scripts/31_analyse-tb.R`). `polyicsflow/` holds the environment-B output.

## Data facts

- TB_s and TB_u: 40 adolescents, all `StudyGroup = passive`. The contrast is `meta$TB`, which the data package
  builds from `QuantiferonInTubeResult` (QFT positive 22, negative 18). It is IGRA infection status, not active TB.
- COMPASS filter: `categories` lists only the kept 7-marker patterns (21 plus all-negative for TB_s, 18 plus
  all-negative for TB_u, of 128). Removal is by omission, and the cells sit in the all-negative count
  (`rowSums(n_s) == counts_s` holds). Positive counts are lower bounds. Every method sees the same tables.
- Collapse to the core 5 markers (IFNg, IL2, TNFa, CD154, IL17a) sums over IL22 and IL4. The two kept patterns
  positive only for IL22 or IL4 become all-negative (TB_s: 4122 stimulated and 3335 control cells).
- TB_s and TB_u share the control sample size (`counts_u` identical); the kept patterns differ.
- SimpleCOMPASS: 256 patterns sum exactly to the parent count in all 42 samples (complete partition). 21 donors:
  EPTB 11, IGRA+ 10 (the brief said 10 and 11). Markers kept: IFNg, IL2, TNF, IL17A, IL22; IL10, IL17F and MIP1b
  are summed out (151290 stimulated and 88182 control cells become all-negative).

## Deviations

- One script, `31_analyse-tb.R`, has three steps (fit, polyics, summary) because environment B must run polyICSFlow.
  `07_run-polyicsflow.R` was not reusable: it has four marker names and the benchmark folder layout.
- polyICSFlow net is summed over the combinations it reports, with a full join of stimulated and control rows
  (missing side = 0 percent); `n_unpaired` is 0 in all three data sets.
- Sensitivity in P(RMS > eps) at 1e-5, 1e-3 and 1e-2, next to the main eps = 1e-4.
