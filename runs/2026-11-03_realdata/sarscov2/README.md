# SARS-CoV-2 Spike ICS benchmark (Zenodo 8088178)

Run of `scripts/32_analyse-sarscov2.R` on `data/processed/sarscov2/counts.tsv` (from `scripts/24_prepare-sarscov2.R`), submitted with `scripts/slurm_realdata_sarscov2.sh` (job 25942136, medium partition).
Package: `cytoJointR-wt/redesign-run`, commit 49278c9. Seed 20261103. Report: `notebooks/realdata_sarscov2.qmd`.

- `CD4/`, `CD8/`: one fit each; Spike and PHA share the DMSO control of the same donor and visit. Tables are described in the header of script 32.
- `CD4/polyicsflow/`, `CD8/polyicsflow/`: polyICSFlow 0.99.4 (environment B) net sums per donor, visit and stimulation.
- `logs/`: job logs.

Data: Ford et al. 2023, doi:10.5281/zenodo.8088178 (CC BY 4.0); paper doi:10.1038/s41590-023-01692-x (visits: E01 before dose 1, E02 after dose 1, E03 after dose 2).
Known gaps: COMPASS failed for CD8 Spike at E01 ("Not a matrix.", one category left after its filter), so 7 donor-visits have no COMPASS Spike call in CD8.
CCS is experimental and is not used for any conclusion.
