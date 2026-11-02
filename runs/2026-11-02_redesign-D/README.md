# Generator D, first attempt (failed, 2026-11-02)

All 240 replicates failed before simulating: `scripts/10_define-heldout.R`
used `mutate(.keep = "unused")`, which dropped the `J` column from
`scenarios.tsv`, so `bench_truth_a()` errored. No model was fitted and no
result exists here. The pre-registered generator (`scripts/lib_generator-d.R`,
commit fed35d9) was not changed. The harness fix keeps all columns; the rerun
is `../2026-11-02_redesign-D2/`.
