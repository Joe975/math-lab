# 013 — 009's separation survives a much wider falsification search

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (executes 009 lead 1 and 010 lead 2)
- **Type:** falsification search
- **Tools:** `explore/tp_constrained_hunt.py` (`--selftest`, exhaustive
  axis-normal sweeps by rank), plus a dense sweep of the 006 family.
  Data: `data/hunt_r3_m1_3_25.json`, `data/hunt_r4_m1_3_25.json`,
  `data/hunt_r4.log`.
- **Sources:** this repo's 006, 009, 010.

## Approach

009 concluded that Cherkaev's B2 separates exactly on his constraint (4.26):
every structure found below B2 violates it, every structure respecting it lies
above. That reading rests on a single laminate family at three grid points, and
it makes a sharp falsifiable prediction — **no structure both respects (4.26)
and lies below B2**. One counterexample overturns it and reopens 006. This
attempt tries hard to find one.

Why this rather than re-deriving Cherkaev's Remark 4.6: the derivation is what
would say who is *right*, but two passes failed to reproduce it, whereas the
prediction is cheap to attack at scale and a hit would be decisive on its own.

## What was done

Three independent searches at σ = (1,2,5), m₂ = 1/4, all in exact arithmetic
after a float screen, classifying every survivor by the exact sign of
(S − ς_N)² − D² over phases 1 and 2:

| search | m₁ | survivors | below B2 | respecting (4.26) | **both** |
|---|---|---|---|---|---|
| exhaustive axis-normal **rank 3** | 3/25 | 621 | 0 | 432 | **0** |
| exhaustive axis-normal **rank 4** | 3/25 | 18 859 | 0 | 12 273 | **0** |
| dense sweep of the 006 family | 31/250 | 981 | 19 | 500 | **0** |
| dense sweep of the 006 family | 3/25 | 992 | 58 | 421 | **0** |
| dense sweep of the 006 family | 11/100 | 1 035 | 197 | 304 | **0** |

**Over 22 000 structures, not one both respects (4.26) and lies below B2.**

The dense family sweeps matter most: they are the only searches that contain
below-B2 structures at all (19, 58 and 197 of them), and every single one
violates the constraint. The minimum among respecting members stays above B2 by
+2.5e−4, +2.5e−3 and +7.6e−3 at the three points, margins that grow as m₁ falls.

Two sanity checks on the classifier itself (`--selftest`): it reports slack
exactly 0 at 004's attaining structure and strictly negative at 006's
below-B2 structure.

## Outcome

`EVIDENCE`, strengthening 009. Scope: σ = (1,2,5), m₂ = 1/4, at m₁ = 31/250,
3/25, 11/100; exhaustive over axis-normal topologies at ranks 3 and 4 at the
middle point, dense within the 006 family at all three; float screen followed by
exact classification.

No counterexample found, so 009's reading stands: **B2 behaves exactly as a
correct bound over the class its derivation constrains**, and every structure
beating it steps outside that class.

**Not claimed.** A proof. This is a finite search: it covers axis-normal
topologies to rank 4 and one family at rank 5, not all microstructures, not
non-axis normals, and not higher ranks. Nothing new about whether (4.26) is
necessary at the optimum — still the decisive question, and still untouched by
any amount of searching.

## Why it failed / what survived

The falsification attempt failed to falsify, which is the informative outcome
here. What survived is 009's reading, now on a much broader base: the
separation is not an artefact of one family, since it holds across two
exhaustive topology classes and a dense parameter sweep, and the only
structures that beat the bound are constraint violators.

The logical position is worth restating precisely, because it is easy to
overstate in either direction. Our structures put the true minimum below B2.
Cherkaev's derivation puts the minimum *over (4.26)-respecting structures*
above B2, and this search confirms that empirically. Both can hold at once only
if the true minimiser violates (4.26) — that is, only if (4.26) is not a valid
necessary condition at the optimum. That is the single remaining disagreement,
and it is a question about his structural-variation argument, not about any
computation.

## Leads generated

1. **Re-derive Remark 4.6** (carried from 008/009). Now the *only* thing that
   can move this, since searching cannot.
2. **Extend to non-axis normals and rank 5** for completeness, though the
   marginal value is low given 22 000 negatives.
3. **Test the separation at a second conductivity triple** — everything so far
   is σ = (1,2,5). If B2 separates on (4.26) there too, the reading is
   essentially secure; if not, something is special about this triple.
4. **Ask the author** (carried, for the human).

## References

- This repo: `attempts/006`, `attempts/009`, `attempts/010`,
  `data/hunt_r3_m1_3_25.json`, `data/hunt_r4_m1_3_25.json`
