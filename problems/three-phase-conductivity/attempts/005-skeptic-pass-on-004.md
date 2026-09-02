# 005 — Skeptic pass on 004: the attaining construction survives, with three corrections

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed
- **Type:** skeptic review (of `attempts/004-closed-form-attainment-at-m11.md`)
- **Tools:** the reviewer's own third implementation and
  `explore/tp_cherkaev.py` (an independently derived construction), both
  standard-library Python; report in `data/skeptic-attain.md`. The precondition
  and identity assertion added in response are in `explore/tp_attain.py`.
- **Sources:** Cherkaev 2009 §8.1 [T] via `data/structures-check.md`; 004.

## Approach

004 claimed an exact object, and the repo's bar for a computational claim is
re-implementation from scratch by a different route, not a re-run. 004's two
agreeing routes both live in the harness it used (`laminate.py`'s projection
form and `verify_laminate.py`'s rotated-frame averaging), so neither counted.
The review was commissioned with the default stance *refute*: build a third
implementation, prove the algebraic side effects by hand rather than by running
004's selftest, hunt adversarially for failing instances, derive the two
pattern-matched parameter identities, and settle novelty.

## What was done, and what it found

**The construction stands. Nothing refuted it.**

1. **Third independent implementation** — the reviewer's own lamination algebra,
   sharing no code with either harness route. Agrees on **590+ instances**
   including extreme conductivity ratios, **zero mismatches**. This is the check
   004 was missing, and it is what lifts 004 from `LIVE` to `VERIFIED`.
2. **m₁₁ < Milton's threshold proved algebraically**, unconditionally for every
   m₂ ∈ (0,1), replacing 004's instance-wise check.
3. **Domain claims proved by hand**: a₀, a₁, a₄, a₅ ∈ (0,1) under the stated
   conditions, unconditionally.
4. **Adversarial hunt: zero failures** across 3 575+ instances plus an exact
   perturbation test, probing a₃ at both endpoints of (1−r, 1), r → 0 and
   r → 1, and near-degenerate phases.
5. **Novelty settled**: an exact, independent rediscovery of Cherkaev 2009
   §8.1's L13,2,13 structure at m₁ = m₁₁ [T]. Recorded as a rediscovery, which
   is what the repo wants from a construction reached without the source's
   formulas.

### Three corrections to 004

- **a₅ = 1 − Θ, an exact algebraic identity.** 004 carried the clumsier
  (σ₂σ₃−σ₁²)/((σ₂+σ₁)(σ₃−σ₁)) found by pattern-matching. Expanding
  1 − Θ gives that numerator identically. So the construction's parameters are
  just a₀ = 1−r, a₁ = rΘ, a₅ = 1−Θ — all in terms of Θ and r alone. This is a
  genuine simplification of the stated result, now asserted in the tool.
- **The object is really rank 4, not rank 5.** Nodes C and B laminate along the
  *same* normal e₂ in sequence, and two consecutive same-normal laminations
  collapse exactly: laminate(laminate(p2@a₄, D; e₂)@a₃, p2; e₂) equals
  laminate(D@w, p2; e₂) with w = a₃(1−a₄), checked as an exact tensor identity.
  Substituting 004's a₄ gives **w = 1−r independent of a₃** — which is exactly
  *why* a₃ is free. So 004's "one-parameter family" is a redundant five-node
  encoding of a four-parameter object, not a genuine geometric degree of
  freedom. Every member is still a bona fide distinct tree that passes every
  check, so nothing computational changes; but 004's framing of a₃ as a real
  freedom, and its billing as rank 5, are both wrong. It also tightens the
  correspondence with Cherkaev's L13,2,13, whose combined weight on D is
  likewise 1−r.
- **σ₁ < σ₂ < σ₃ is a precondition, not a convention.** 004's code never
  checked it, while `hs_bounds()` independently takes min and max of the
  conductivities, so an unordered triple would silently compare the structure
  against a bound built from a different comparison medium. No test hit this
  because every test built sorted triples. Now raised explicitly in
  `check_ordering()`.

### What remains open

The reviewer showed a₁ = rΘ and a₅ = 1−Θ agree exactly with an independently
derived closed form (`tp_cherkaev.py`, reached by direct algebraic solving of
the L13,2,13 family in ℚ(√m₂) and matching node by node after the
laminate(A@m,B) = laminate(B@(1−m),A) symmetry). That is corroboration by a
second derivation, not yet a first-principles derivation from 002's field
conditions. Upgrading it is a strengthening, not a doubt.

## Outcome

`VERIFIED` — scoped to the construction of 004 as stated: exact attainment of
the HS lower bound at f₁ = m₁₁ = 2Θ√m₂(1−√m₂), f₂ = m₂, for σ₁ < σ₂ < σ₃ and
a₃ ∈ (1−r, 1) with f₃ > 0. Established by three independent implementations
agreeing exactly, hand proofs of the algebraic side claims, and an adversarial
sweep of 3 575+ instances with zero failures. Recorded as an **independent
rediscovery** of Cherkaev 2009 §8.1.

**Not claimed.** Nothing about f₁ < m₁₁. No first-principles derivation of the
parameters. The rank-4 collapse is about *this* encoding; whether rank 4 is
minimal for attainment at m₁₁ is still open (004 lead 4, now sharper since the
object is already rank 4 — the question becomes whether rank 3 can do it).

## Why it failed / what survived

Nothing failed. The corrections are the value: a simpler parameterisation, a
lower true rank, and a silent precondition that could have produced a
false-positive attainment claim for an unordered triple. The middle one matters
for the record — 004's headline rank was wrong, and its free parameter was an
artefact of a redundant encoding rather than the geometric freedom it was
presented as.

Method note worth keeping: the free parameter that made 004's float optimizer
wander, which 004 read as a signature of an exact family, was really a
signature of a *redundant encoding*. Both readings predict a flat direction;
only the algebra distinguishes them.

## Leads generated

1. **Derive a₁ = rΘ and a₅ = 1−Θ from 002's field conditions** via the
   path-product rule, turning the construction into a proof rather than a
   verified formula.
2. **Is rank 3 enough?** Now that the object is known to be rank 4, exhaustively
   search rank 3 at exactly f₁ = m₁₁; cheap and definite.
3. **Is m₁₁ the true threshold for all microstructures?** Still the frontier.
   Within 004's family the gap is 0 at m₁₁ and strictly positive below, but that
   is evidence about one family.
4. **Assert preconditions across the other tools** — the ordering hazard found
   here may exist in `tp_coated.py` and `tp_shapes.py` too. Cheap audit.

## References

- `attempts/004-closed-form-attainment-at-m11.md` (the reviewed record)
- `data/skeptic-attain.md` (the full review, with its exact instance counts)
- `explore/tp_cherkaev.py` (the independent derivation)
- A. Cherkaev, Mech. Mater. 41 (2009) 411–433; arXiv:1009.3060, §8.1 [T]
