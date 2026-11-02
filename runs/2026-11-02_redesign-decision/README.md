# Model redesign decision, 2026-11-02

**Default model: C2 (`m1dm`, M1 mean structure plus a Dirichlet-multinomial
sample residual).** This was chosen by the user after the pre-registered
rule found no calibrated model. It is a post hoc decision, recorded here with
its reasons.

## What the rule said

`scripts/11_apply-decision-rule.R` on `../2026-11-02_redesign-summary/`
(545 replicates, 18 A/B/C cells, sample target) returned **no default**.
`decision.tsv` and `decision.md` here are its unedited output. One correction:
the last sentence of `decision.md` is a fixed text that still names CCS. The
run applied the amended score set, RMS and CPS.

Pre-registration trail (outer repo commits, in order):
1. `fed35d9`: held-out generator D committed before any candidate ran.
2. `6e51462`: amendment made before any result. The band is widened by
   2 seed-clustered SEs, and null cells have no upper bound.
3. Post hoc commit: decide on RMS and CPS only; CCS is experimental.
   - Reason: no candidate covered CCS under generator B (best 0.34).

## Why C2, against the rule

After step 3, the models still failed cells, but in opposite directions:

| Model | Failing cells | Direction | CPS coverage in failing cells |
|---|---|---|---|
| C1 (M1) | 9 of 18 | under-coverage | 0.24 to 0.79 |
| C2 (m1dm) | 12 of 18 | over-coverage | 0.988 to 0.999 (nominal 0.90) |
| C0 (Dirichlet) | 13 of 18 | mixed; CCS false signs up to 0.99 | |

These facts support C2:
- RMS coverage is within the band for C1 and C2 in every cell.
- C2's failures are all conservative.
- RMS error on rare_weak plus rare_moderate is 31% below C0's for both C1
  and C2.
- C1's failures understate uncertainty, which is disqualifying.
- C2's cost is wider CPS intervals than needed.

The user chose C2 on 2026-11-02 (option 1 of 2). The alternative, no default
until C2's CPS over-coverage is fixed, was declined.

## Consequences for the package

- `m1dm` is the default; `m1` stays available.
- CPS intervals are documented as conservative.
- CCS is marked experimental: it is not calibrated under donor heterogeneity.
- Follow-up: find why C2 over-covers CPS. It is unverified, but the
  Dirichlet residual may be inflating ratio-type scores.

Held-out generator D was opened only after this decision.
