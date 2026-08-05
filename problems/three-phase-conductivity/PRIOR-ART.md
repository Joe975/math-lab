# Three-phase conductivity — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None.** This problem was onboarded (2026-08-04) and has not been worked.
There is no prior art to be informed by, so `blind` and `informed` mode are
currently equivalent here — which makes the first attempts worth running
blind, since blind costs nothing while the record is empty.

## Editorial view of the attack surface

Homogenization / optimal composites is a new attack surface for the lab,
arrived at from a physics direction (metamaterial engineering of effective
response; see `docs/PLAN-conductance-problems.md` for the scoping decision).
Unlike the lab's famous long shots this frontier moved recently and is
genuinely contested: the best three-phase bounds (Nesi 1995, Cherkaev 2009)
and best structures meet on part of parameter space and leave a documented
gap elsewhere. A careful census can add something real.

The exact-arithmetic fit is the best since mahler-4d: hierarchical laminates
are rational functions of rational data, so attainment is an equality in ℚ
and the Keller–Dykhne duality identity gives every record a free internal
cross-check. The harness computes effective tensors by two independent
algebraic routes (projection formula vs rotated-frame interface averaging)
that must agree exactly.

Two cautions. First, the improved bounds (Nesi, Cherkaev) are *not* in the
harness — they need careful transcription from the papers, which is attempt
work, must be marked as transcription per CONTRIBUTING rule 6, and should be
cross-checked against the papers' own worked examples before anything is
killed against them. Second, a bounded-rank laminate search that stalls below
a bound proves nothing about all microstructures — Milton's and Cherkaev's
optimal assemblages are exactly the warning that new structure classes can
appear; scope every census claim by its class.

Concrete lines, if you want them:

- Self-test attempt: two-phase ground truth. Rank-2 laminates attaining the
  2D two-phase HS bounds exactly in ℚ at several (f, σ) points, duality on
  nontrivial trees, series/parallel closed forms. Expected `VERIFIED`,
  scope = the checks run; its purpose is to validate the tooling against
  known answers before any open-regime work.
- Three-phase attainability map: fix rational (σ₁, σ₂, σ₃), sweep a rational
  grid over the fraction simplex, optimize bounded-rank laminates
  (float screen, exact certification of the best found), and chart the gap
  to the HS bound — where it closes, where it stalls. `MAP`/`EVIDENCE`,
  scoped by rank, direction set and grid.
- Transcribe the Nesi and Cherkaev bounds as attempt tooling and re-run the
  map against them; where laminate optima match the tighter bounds exactly,
  that is attainability evidence worth isolating as candidate exact
  structures.

Kill condition, inline: if bounded-rank optima plateau strictly inside the
bounds across the whole grid with no structure in the residual gap, the
census route yields only a negative map — recordable once, then cap the
budget rather than refining the grid indefinitely.
