# Donor-effects round 2 (labels as in the leaderboard)

Package worktree `cytoJointR-wt/de-r2` (commit 7f125c9) and `de-c6` (3188b90). Dev mode only; test mode not run.
`c8-resp-deg-basalmarker`: `responder=TRUE,degree_terms=TRUE,donor_loading=marker_basal,iter_warmup=400`. Result: worse than c7 (AUC 0.878).
`c7f-fast`: c7 with `iter_warmup=400`. Same metrics as c7, fit 16 min instead of 36 min.
Cancelled before finishing, no leaderboard rows: `c8-resp-deg-stimdeg` (stim_degree) and `c10-resp-deg-soft` (soft spike) were stopped when the scope changed (other agents own the soft spike); `c9-resp-deg-basalmarker` hit the 1 h limit on `short`, re-run as `c8-resp-deg-basalmarker`. Their partial sim directories are unused.
The PPC fit of c8 (`ppc-c8-resp-deg-basalmarker`) timed out at 1 h.
TB_u diagnosis: `runs/2026-11-04_donor-effects/tb_u_diagnosis.md`.
