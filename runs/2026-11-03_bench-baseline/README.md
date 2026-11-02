# Benchmark baseline (2026-11-03)

Raw scores of the first leaderboard rows. Package: `cytoJointR-wt/bench-baseline`, commit e54b883 (main). Definition commit 32d9447, scorecard commit 6013eb0.
Run with `scripts/slurm_scorecard.sh` (array jobs 25945249, 25945250, 25945251 for the fits; 25945481 to 25945483 for collect).

- `baseline-m1dm-sample/`: `model=m1dm,level=sample` (the package default).
- `baseline-m1dm-structured/`: `model=m1dm,level=structured`.
- `baseline-m1/`: `model=m1,level=sample`.

Each holds `rv144_scores.tsv`, `tb_s_scores.tsv`, `tb_u_scores.tsv`, `sim/<scenario>/seed_NN.tsv`, fit diagnostics and logs.
`rv144_scores.tsv` contains test-half donors too. Do not open it for the test half before `--test`.
Test mode was not run (test-set access counter 0).
The metrics are in `benchmark/leaderboard.tsv`.
