# 007 — Exact attainment at the original target point, via a coated T² structure in ℚ(√105)

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (follows 001–005; builds on 003's coating lemma)
- **Type:** construction (object), exact in a quadratic number field
- **Tools:** `explore/tp_cherkaev.py` (`--selftest`, `--build`), which carries its
  own two independent lamination algorithms generic over a field. Records:
  `data/cherkaev/cherkaev_rational_{lo,hi}.json`, both passing
  `verify_laminate.py`; `data/cherkaev/cherkaev_target.json` is a target
  descriptor, not a laminate record.
- **Sources:** A. Cherkaev (2009), §8.1 and Thm 8.1, via
  `data/structures-check.md` and `data/cherkaev-primary-excerpts.md` [T];
  attempts 003 (the coating lemma) and 004 (the m₁₁ construction).

## Approach

004 attained the HS lower bound at f₁ = m₁₁ exactly, and 003's coating lemma
extends that up the ray. But the point this whole investigation started from —
σ = (1,2,5), f = (1/8, 1/8, 3/4), where 001's bounded-rank searches stalled at
5.46e−4 — is **not** at m₁₁ for its own m₂, so neither result reached it
directly. This attempt closes that specific point.

The move that works: do not fix f₁ = 1/8 and search. Instead build the inner
T²-structure at *its own* m₁₁, where it attains exactly, and then choose the
inner m₂ so that coating it with extra phase 1 lands on the target fractions.
The inner m₂ is then not the target m₂, which is exactly why local searches
around f₁ = 1/8 never found this — the required inner structure is at a
different, irrational, volume fraction.

## What was done

### 1. A closed form for the T² structure, derived not transcribed

The literature passes could not reconstruct Cherkaev's own formulas
(`data/structures-check.md` records the failure). So the L13,2,13 family was
re-derived here. For σ = (1,2,5), with Θ = 1/4 and u = √m₂, the rank-4
axis-normal tree

    X = laminate(p3, p1; fraction b; e₂)
    Z = laminate(p1, p3; fraction d; e₁)
    Y = laminate(Z, p2;  fraction c; e₂)
    T² = laminate(X, Y;  fraction a; e₁)

with **a = c = 1−u, d = Θ = 1/4, b = 1−u/4** is isotropic and equals HS_lo
exactly at m₁ = m₁₁(m₂) = 2Θ·u(1−u). Verified by exact field arithmetic at five
values of m₂, by two independent lamination algorithms (the projection formula
and the rotate-and-average algebra, both reimplemented generically over a
field). This is an independent confirmation of Cherkaev's Thm 8.1 claim, and it
agrees exactly with 004's construction after the swap symmetry
laminate(A@m, B) = laminate(B@(1−m), A) — see 005 §5, where the two derivations
were shown to be the same object.

### 2. Hitting the target point

Requiring the coated structure to land on f₂:f₃ = 1:6 forces the inner
structure's own m₂ to satisfy

    13x² + x − 2 = 0,   x = √(inner m₂),   so x = (−1 + √105)/26,

and the whole construction lives in the single quadratic extension ℚ(√105).
The inner T² sits at g = ((−33+7√105)/338, (53−√105)/338, (159−3√105)/169) and
attains its own HS_lo exactly there; coating it with extra phase 1 in two
orthogonal equal-fraction layers — **003's lemma, which is what makes this
work** — lands exactly on fractions (1/8, 1/8, 3/4) and exactly on value
**37/11 = HS_lo**.

Every intermediate node fraction is irrational, yet the final tensor and the
final volume fractions are exactly rational. The check is literal equality in
the field: `tensor == 37/11 + 0·√105`.

**Rank: 6** (4-node T² core plus 2-node coating).

### 3. Verification

- Exact symbolic equality in ℚ(√105), zero floating point, with two independent
  lamination algorithms agreeing; they were also cross-checked against each
  other on an unrelated generic anisotropic tree.
- Harness cross-check: √105 snapped to two rational approximations of width
  2⁻¹⁹⁰ (one below, one above), the identical tree rebuilt as plain Fractions,
  and both run through the actual `laminate.py` and `verify_laminate.py`. Both
  pass, Keller–Dykhne holds, and the residual from 37/11 is 7.8e−58 and
  −7.7e−60 respectively — consistent with the precision of the √105
  approximation, not with a real gap.

### 4. Why 001's and 003's searches missed it

Confirmed by Newton-solving the same topology at fixed f₁ = 1/8: **no valid
solution exists in the physical domain** there. The plain T² attains only at
f₁ = m₁₁ ≈ 0.1143, never at 0.125. The coating step is structurally necessary,
and it requires the inner structure to sit at a different m₂ than the target —
invisible to any local perturbation around f₁ = 1/8. That is the precise reason
001's rank-3/4 screens and 003's rank-5 sweep converged toward the bound
without reaching it.

Also checked and discarded: coating a *Milton-type* assemblage merely
re-derives Milton's own threshold and gets nowhere below it. The coating lemma
only buys anything when the inner structure attains by a non-Milton mechanism,
i.e. the T².

## Outcome

`VERIFIED` for the construction as stated — exact attainment of HS_lo at
σ = (1,2,5), f = (1/8, 1/8, 3/4), by a rank-6 laminate, established by exact
arithmetic in ℚ(√105) with two independent lamination algorithms plus two
rational-snap witnesses passing the harness verifier.

This **resolves 003's open tension** at this point in favour of attainment: the
near-constant per-rank gap ratio in the axis-normal sweeps was not convergence
to an unattainable bound, it was a search space that did not contain the
attaining object.

**Not claimed.** Nothing below m₁₁ (that is 006's unresolved escalation, and is
untouched by this). No claim that rank 6 is minimal here. Novelty: this is an
independent constructive confirmation of what Cherkaev 2009 Thm 8.1 asserts,
reached without his formulas, and should be recorded as such.

## Why it failed / what survived

Nothing failed. The transferable lesson is about search design: the attaining
object at a target point was **not reachable by any local search at that
point**, because the construction routes through an inner structure at a
different, irrational volume fraction. Fixing the target fractions and
searching over node parameters — the obvious design, used in 001 and 003 — is
structurally incapable of finding it. What worked was building at the place
where attainment is easy (the inner structure's own threshold) and then using a
proved transformation (003's coating lemma) to move to the target.

Reusable: `tp_cherkaev.py`'s field-generic lamination algorithms, which work
over ℚ(√n) and so handle constructions whose parameters are irrational but whose
answers are not; the T² closed form; the rational-snap technique for pushing a
quadratic-field construction through a rational-only harness.

## Leads generated

1. **Generalise the target-hitting step.** Here the inner m₂ solved a quadratic.
   For an arbitrary target (σ, f) above m₁₁, does the analogous condition always
   have a solution in the physical domain, and of what algebraic degree? If yes,
   004 + 003 + this give a construction for the whole attainable region.
2. **Minimal rank at the target.** Rank 6 here versus rank 4 at m₁₁ (005). Is
   rank 5 enough at (1/8,1/8,3/4)? 003's rank-5 sweep says no within axis
   normals, but it also could not have found this object.
3. **Wrap `tp_cherkaev.py` into the pytest gate.** Its `--selftest` and
   `--build` currently stand outside `pytest tests/ -q`; the rational-snap
   witnesses are cheap enough to run there.
4. **Feed this back into 006.** The below-m₁₁ discrepancy used a family that
   cannot express coated T² structures. Re-running the below-m₁₁ search over
   *this* richer class would either lower our values further, sharpening the
   conflict with B2, or reveal that the family was the limitation.

## References

- A. Cherkaev, Mech. Mater. 41 (2009) 411–433, §8.1 and Thm 8.1 [T]
- P. Albin, A. Cherkaev, V. Nesi, J. Mech. Phys. Solids 55 (2007) 1513–1553 [T]
- This repo: `attempts/003` (the coating lemma), `attempts/004`, `attempts/005`,
  `data/structures-check.md`, `data/cherkaev/`
