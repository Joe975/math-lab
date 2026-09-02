# 003 — Attainability is upward closed in f₁ (proved), why the rank-3 shape dies, and the literature's verdict

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (follows 001 and 002)
- **Type:** proof + exact computation + literature check
- **Tools:** `explore/tp_coated.py`, `explore/tp_fields.py`, `explore/tp_shapes.py`
  (all from 001/002), `explore/tp_rank5.py` (new, delegated).
  Tests: `tests/test_three_phase_fields.py` (7). Reports:
  `data/structures-check.md` (new literature pass), `data/rank5/`.
- **Sources:** Cherkaev 2009 (Mech. Mater. 41, 411–433; preprint
  arXiv:1009.3060) §8.1 and Thm 8.1 [T]; Albin–Cherkaev–Nesi 2007 (J. Mech.
  Phys. Solids 55, 1513–1553) Thm 2, Thm 6, Remark 5 [T]; Gibiansky–Sigmund
  2000 (J. Mech. Phys. Solids 48, 461–498) via Cherkaev's Remark 8.2 [T],
  with the caveat below. Full quotations in `data/structures-check.md`.

## Approach

002 supplied diagnostics but left the central question open and flagged a
`SPECULATION` that low-rank searches are parameter-starved. This attempt
attacks the question from three sides at once: prove something structural
about the *shape* of the attainable region, explain mechanistically why 001's
rank-3 family died where it did, and settle what the literature actually
claims. The three turn out to interlock.

## What was done

### 1. Coating with the comparison medium preserves HS-optimality (proved)

**Lemma.** Let a structure attain the three-phase HS lower bound (comparison
medium σ₁) at fractions g. Coat it with extra phase 1 in two orthogonal
normals with equal fractions, leaving the core at volume fraction F. Then the
result is isotropic and attains the HS lower bound at the new fractions.

*Proof.* The coating is a two-phase isotropic coated laminate of an isotropic
core v₀ (fraction F) in matrix σ₁, which attains the two-phase HS lower bound
(established exactly in 001 §1), so

    1/(value + σ₁) = F/(v₀ + σ₁) + (1−F)/(2σ₁).

The core attains its own bound, so 1/(v₀+σ₁) = (1/F)·Σ_inner f_i/(σ_i+σ₁).
Substituting, the first term becomes Σ_inner f_i/(σ_i+σ₁); the second is the
*outer* phase-1 contribution, since its comparison-medium denominator is
σ₁+σ₁ = 2σ₁. Together they are Σ over all phases, i.e. 1/(HS_lo+σ₁) at the
new fractions. ∎

Verified exactly in ℚ on 12 instances (3 base structures × 4 coating
fractions): value − HS_lo = 0 to the last digit, isotropy exact, and the
required coating parameter comes out ρ = 1/2 every time
(`test_coating_with_comparison_medium_preserves_hs_optimality`).

**Corollary (upward closure).** Coating scales f₂ and f₃ by the same factor,
so the ratio f₂:f₃ is preserved while f₁ increases. Hence along any ray of
fixed f₂:f₃, **the set of f₁ at which HS_lo is attainable is upward closed**,
and attainability is governed by a single threshold f₁\*(f₂:f₃). It also costs
only two extra laminations to move up the ray, so minimal *rank* is a question
about the bottom of the ray, not about its interior.

This retro-explains 001: the Milton assemblage's feasible set {HS_lo ≤ σ₂} is
exactly such an upward-closed region, and its boundary is where phase 2 goes
uncoated (φ₂ = 1). Going below Milton's threshold is therefore not about
finding better structures everywhere, only about finding one better structure
at the bottom.

### 2. Why 001's rank-3 shape dies at f₁\* ≈ 0.0301

002's `Var_req = (c²S₂−1)/f₁` diverges as f₁ → 0. Tracking 001's rank-3
family down its ray (f₂:f₃ = 1:6, σ = (1,2,5)), the achieved phase-1 variance
runs the *other* way and collapses:

| f₁ | gap | Var achieved | Var required | ratio | smaller p1 leaf / f₁ |
|---|---|---|---|---|---|
| 0.125 | 2.1e−2 | 1.832 | 2.050 | 0.894 | 0.465 |
| 0.0625 | 3.8e−2 | 2.033 | 3.199 | 0.636 | 0.310 |
| 0.0333 | 5.0e−2 | 0.586 | 4.744 | 0.124 | 0.058 |
| 0.0303 | 5.1e−2 | 0.048 | 5.051 | 0.0095 | 0.0044 |
| 0.0294 | — | no isotropic member of this shape | | | |

The mechanism is in the last column. Phase 1 occupies two leaves, so its
variance is (v₁v₂/f₁²)·|ΔE|²; the field separation |ΔE| stays near-constant
(2.71 → 3.30), but **one of the two leaf volumes is squeezed out**, falling to
0.4% of f₁. The variance collapses with it, and the shape's isotropic member
disappears exactly there. So 001's feasibility floor is a volume degeneration,
not a field one: the shape runs out of the phase-1 splitting it needs to
supply a diverging required fluctuation.

Testing the natural prediction — that a shape with more freedom survives
lower — the rank-4 family of 001 does survive well past the rank-3 floor
(computed down to f₁ = 1/96 ≈ 0.0104), though its gap grows steadily
(5.46e−4 at f₁ = 1/8 up to 1.18e−1 at f₁ = 1/96) and its variance ratio sits
at 0.24–0.33 there. Note the rank-4 shape still has only *two* phase-1 leaves,
so what buys survival is the extra free parameter, not extra leaves —
correcting the guess that leaf count is what matters.

### 3. Rank 5, and the gap sequence

A delegated rank-5 search over axis-normal topologies
(`explore/tp_rank5.py`, `data/rank5/`) at f = (1,1,6)/8, σ = (1,2,5):

| rank | best gap above HS_lo | ratio to previous |
|---|---|---|
| 3 | 2.116e−2 | — |
| 4 | 5.456e−4 | 38.8 |
| 5 | 1.6e−5 | ≈34 |

The rank-5 sweep was still running when this record was written (about half
the topologies enumerated), so 1.6e-5 is the best value found so far and
therefore an upper bound on the rank-5 optimum gap, not the optimum. The
per-rank ratio is near-constant, which is geometric decay toward HS_lo
rather than a plateau. `SPECULATION` (labelled): a *geometric* decay is
consistent with the bound being approached but not attained at any finite
rank within this axis-normal family, which would sit oddly beside the
literature's finite-rank attainment claim in §4; the alternative is that the
optimizer is not finding the exact attaining topology. Distinguishing these is
lead 1.

### 4. What the literature actually claims — HS is attained here

The interval left open by 001 was m₁₁ ≤ m₁ < 2Θ(1−m₂), where Cherkaev's 2009
bound degenerates to HS and Milton's assemblage does not reach it. A dedicated
literature pass (`data/structures-check.md`) finds that **the published
position is that HS_lo is attained on this whole interval**:

- Cherkaev 2009 §8.1 [T] describes a structure "L13,2,13,1,1": a T²-structure
  carrying phase 1 in amount exactly m₁₁, then "sequentially laminated by the
  two orthogonal layers of k₁" in the amount m₁ − m₁₁. Theorem 8.1 [T]: "The
  bound … is exact in each point: There exist laminates of a finite rank that
  realize the bounds."
- ACN 2007 Thm 2 asserts an isotropic HS-optimal structure exists on
  m₁₁ ≤ m₁ ≤ 2Θ(1−m₂), with their own Remark 5 flagging it as not yet
  constructive at the time; their Thm 6 aims to supply the construction.
- No source found states or implies non-attainability there.

**The outer part of Cherkaev's structure is exactly this attempt's §1 lemma**
— extra phase 1 added as two orthogonal laminations — which is why it
preserves both isotropy and optimality, and why his threshold is a bottom-of-
ray condition. Our independently proved lemma and his construction agree, and
neither was derived from the other.

Two caveats recorded rather than smoothed over. Cherkaev's Remark 8.2 credits
Gibiansky–Sigmund 2000 with proving optimality on m₁ ∈ [m₁₁, 1], but that
reference is an *elasticity* (bulk modulus) paper; the conductivity reading is
Cherkaev's own reinterpretation, which ACN 2007 also make explicitly (their
footnote 3 [T]). And the worker could not reconstruct the T²-structure's
isotropic member from the transcribed formulas at our numeric point, so we do
not yet hold a ready-to-verify exact tree.

## Outcome

- `VERIFIED` (proved, and exact in ℚ on 12 instances): coating with the
  comparison medium preserves HS-optimality; hence attainability is upward
  closed in f₁ along fixed f₂:f₃.
- `EVIDENCE` (scope: 001's rank-3 family on the f₂:f₃ = 1:6 ray at
  σ = (1,2,5), exact at each computed f₁): the feasibility floor is a
  phase-1 leaf-volume degeneration, with the required variance diverging as
  the achieved variance collapses.
- `EVIDENCE` (scope: axis-normal topologies, seeded search, f = (1,1,6)/8,
  σ = (1,2,5)): gaps 2.1e−2, 5.46e−4, 1.6e−5 at ranks 3, 4, 5, ratio ≈ 35–39.
- `MAP`: the literature claims attainment on the whole interval, with the
  two source caveats above.

**Not claimed.** We have not built or verified an HS-attaining structure below
Milton's threshold; the rank-5 number is a small gap, not an attainment. The
literature verdict rests on transcriptions of two papers, one of whose
attribution chain runs through an elasticity result. Nothing here bears on
whether HS is attainable below m₁₁, which is a different and untouched
question.

## Why it failed / what survived

Nothing failed; the attempt's own §3 raises the sharpest open tension, between
a geometric gap decay and a finite-rank attainment claim. The strongest
survivor is §1: a short proof, exactly verified, that reduces "where is HS
attainable" from a two-dimensional region question to a threshold on each
ray, and reduces the construction problem to the bottom of the ray. Combined
with §4 it makes the next step concrete and checkable rather than exploratory.

Reusable: the coating lemma as a *construction generator* (any HS-optimal
structure yields a one-parameter family of them at larger f₁, at a cost of two
laminations); the variance-collapse diagnostic as an explanation for why a
laminate family has a feasibility floor; `tp_rank5.py`.

## Leads generated

1. **Build Cherkaev's L13,2,13,1,1 explicitly and check it in ℚ.** Work
   through Cherkaev 2009 §8.1 and the Appendix (§9, eqs. 9.1–9.5) into a
   lamination tree the harness can consume, at σ = (1,2,5), f = (1,1,6)/8, and
   test `value == HS_lo` exactly. Definite outcome. If it attains, §3's
   geometric decay is an artefact of the axis-normal restriction or the
   optimizer, and we have the object; if it does not, we have a concrete
   discrepancy with a published theorem, which is worth much more.
2. **Find the true bottom of the ray in our own search.** By §1 it suffices to
   attain at the smallest f₁ possible. Search directly for exact attainment at
   f₁ just below Milton's 2Θ(1−m₂) = 0.4375 and walk down; the first f₁ where
   exact attainment fails brackets our constructive threshold and can be
   compared with m₁₁ = 0.114277.
3. **Settle §3's tension.** Push the axis-normal family to rank 6 and 7: a
   continuing constant ratio argues the family only converges, whereas a sudden
   drop to zero identifies the attaining rank. Either resolves the
   `SPECULATION`.
4. **Does the coating lemma have a dual?** Keller–Dykhne (002 §6) should turn
   it into "coating with the *largest* phase preserves HS-upper-optimality",
   giving an upward-closure statement on the upper side for free. Cheap to
   state and verify.
5. **Test whether m₁₁ is the true threshold** by checking whether the required
   variance at m₁₁ is exactly what a T²-structure can supply — i.e. whether the
   §2 mechanism reproduces Cherkaev's threshold as a variance-feasibility
   boundary. If it does, the §2 diagnostic predicts thresholds rather than just
   explaining one.

## References

- A. Cherkaev, Mech. Mater. 41 (2009) 411–433; preprint arXiv:1009.3060. [T]
- P. Albin, A. Cherkaev, V. Nesi, J. Mech. Phys. Solids 55 (2007) 1513–1553. [T]
- L. Gibiansky, O. Sigmund, J. Mech. Phys. Solids 48 (2000) 461–498 — cited
  via Cherkaev's Remark 8.2; an elasticity paper, see caveat in §4. [T]
- G. W. Milton, Appl. Phys. A 26 (1981) 125–130 — not reached; see 001.
- This repo: `attempts/001-coated-laminate-census.md`,
  `attempts/002-hs-attainment-field-mechanism.md`,
  `data/structures-check.md`, `data/rank5/`.
