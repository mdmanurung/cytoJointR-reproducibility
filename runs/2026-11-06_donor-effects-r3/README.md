# Donor-effects round 3: polyfunctional response offsets

Package worktree `cytoJointR-wt/de-r3` (commit 208c17d), dev mode only, test mode not run.
Both labels are c7f-fast (`responder=TRUE,degree_terms=TRUE,iter_warmup=400`, `level=structured`) plus:
- `c11-fast-stimdeg3`: `stim_degree=TRUE` (stimulation-specific offsets for degrees 3..J).
- `c12-fast-stimdeg2`: `stim_degree=2` (degrees 2..J).
`ppc-<label>/rv144/` hold the `60_diagnose-estimand.R` outputs.
