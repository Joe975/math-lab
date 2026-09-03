# 018 — No tilt of the normals can beat the family curve, and 016's enrichment ties rather than loses

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (executes 015 lead 2; reviews 016)
- **Type:** exact computation + skeptic review of 016
- **Tools:** `explore/tp_enrich.py` (new; degree-2 jets in seven variables over
  ℚ, exact rational linear algebra, `--selftest --report --recheck-016`), built
  on 017's `explore/tp_curvature.py`. Tests:
  `tests/test_three_phase_curvature.py`. Deterministic, standard library only.
- **Sources:** this repo's 004, 005, 015, 016, 017.

## Approach

015's most informative open item was lead 2: **does a richer class drop below
the family's minimum below m₁₁?** 016 attacked it by re-laying node B along the
orthogonal normal and grid-searching two parameters, and reported the enriched
class landing 1–2 orders of magnitude *worse* (ratio 65.87 and 6.71 against
0.53). Two things about that answer were worth re-examining.

First, 017 makes the question local. The family's gap below m₁₁ is γ(m₁−m₁₁)²
with γ the Schur complement of the gap's Hessian at 004's attaining point. A
class that **contains** that point has the same base point, and enlarging the
parameter vector can only lower a Schur complement. So "does a richer class win?"
is finite exact linear algebra at one point — one Hessian, no search, an answer
covering every m₁ just below m₁₁ at once.

Second, 016's class does not obviously exclude the attaining structure the way
016 claims ("containing no member of the original family"), because node B's
normal stops mattering when its p₂ layer has vanishing weight.

## What was done

### 1. The largest fixed-topology enrichment: tilt every normal

Give each of the five internal nodes a tilted normal, all reducing to the
axis-aligned base at t = 0:

    root (1, t₀)   A (t_A, 1)   D (1, t_D)   C (t_C, 1)   B (t_B, 1).

Tilting breaks the automatic vanishing of σ₁₂, so isotropy becomes **two**
conditions instead of one. a₅ and t₀ are eliminated from them as series (the
2×2 Jacobian is invertible at the base), leaving six free parameters
a₀, a₃, t_A, t_D, t_C, t_B. a₃ is redundant in the axis-aligned family (005) but
not once normals tilt, so it is carried as a genuine parameter. This class
strictly contains the 006 family and contains 004's attaining structure at
t = 0.

The gap g = value − HS_lo is expanded to degree 2 in (u, a₀, a₃, t_A, t_D, t_C,
t_B) with u = m₁ − m₁₁, and γ is the Schur complement of the u-row.

```
python tp_enrich.py --selftest
python tp_enrich.py --report
```

### 2. Result: γ does not move at all

At σ = (1,2,5), m₂ = 1/4, the gap's Hessian at the attaining point is

|  | u | a₀ | a₃ | t_A | t_D | t_C | t_B |
|---|---|---|---|---|---|---|---|
| **u** | 352/27 | 20/27 | **0** | **0** | **0** | **0** | **0** |
| **a₀** | 20/27 | 1.31481 | 0 | 0 | 0 | 0 | 0 |
| **a₃** | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **t_A** | 0 | 0 | 0 | 2.59259 | 3.11111 | 0.34568 | 0.17284 |
| **t_D** | 0 | 0 | 0 | 3.11111 | 3.83333 | 0.48148 | 0.24074 |
| **t_C** | 0 | 0 | 0 | 0.34568 | 0.48148 | 0.10905 | 0.02675 |
| **t_B** | 0 | 0 | 0 | 0.17284 | 0.24074 | 0.02675 | 0.04115 |

It **block-diagonalises**: the tilt directions are PSD among themselves and
exactly decoupled from both u and a₀. So the optimal tilt stays at 0 for every
m₁, and

    γ(a₀ only) = γ(a₀, a₃, t_A, t_D, t_C, t_B),  exactly, on all five triples.

| σ | γ, a₀ only | γ, all six | c = γ/β |
|---|---|---|---|
| (1,2,5) | 6.30985915 | 6.30985915 | 189/355 |
| (1,3,7) | 21.07784431 | 21.07784431 | 385/668 |
| (1,4,9) | 49.28492849 | 49.28492849 | 665/1111 |
| (1,2,9) | 14.88299445 | 14.88299445 | 191675/277733 |
| (2,5,11) | 21.19923746 | 21.19923746 | 1209312/2270015 |

A flat direction cannot be exploited at higher order either: if the quadratic
form has no dependence on a direction, every contribution from it enters at
cubic order or beyond and cannot change the u² coefficient.

**Guard against the trivial failure mode.** The obvious way to get "tilts do
nothing" is a bug in which the tilts never enter. The selftest asserts that a
tilt of 1/10 at node A makes the tensor anisotropic *before* the (a₅, t₀)
re-solve (σ₁₂ ≠ 0 and σ₁₁ ≠ σ₂₂). Numerically, re-solving isotropy at
m₁ = m₁₁ with t_A = 0.01 / 0.05 / 0.1 moves a₅ from 0.75 to 0.749909 / 0.747740
/ 0.741180 and t₀ to −0.035 / −0.176 / −0.354, and raises the value above HS_lo
by 1.2958e−4 / 3.2115e−3 / 1.2506e−2 — against the Hessian's predictions
1.2963e−4 / 3.2407e−3 / 1.2963e−2, and symmetric in ±t_A. The tilts are real,
they produce genuinely different isotropic structures, and they are strictly
uphill.

Direct minimisation at m₁ = m₁₁ − 10⁻³ confirms it end to end: min over a₀ at
t_A = 0 gives gap/h² = 6.3429, at t_A = ±0.005 gives 37.44, at ±0.02 gives
495.3.

### 3. 016 re-checked: its enrichment ties, it does not lose

016's altered topology (node B along e₁) with a₃ swept, minimised over a₀ in
exact ℚ with a₅ from exact isotropy bisection (`python tp_enrich.py
--recheck-016`), at σ = (1,2,5), m₂ = 1/4:

| a₃ | m₁ = 3/25, ratio | m₁ = 11/100, ratio |
|---|---|---|
| 1/2 | 85.19 | 8.57 |
| 3/5 | **65.21** | **6.58** |
| 9/10 | 14.23 | 1.78 |
| 99/100 | 1.836 | 0.647 |
| 999/1000 | 0.661 | 0.541 |
| 1 − 10⁻⁵ | **0.5326** | **0.5289** |

016's reported minima (65.87 and 6.71) sit at a₃ ≈ 0.6. They are an **interior
local minimum, not the class minimum**: the value decreases monotonically as
a₃ → 1 and converges to the original family's minimum (ratios 0.531326 and
0.528734, from 015). At a₃ = 1 − 10⁻⁵ — a perfectly ordinary laminate, no
degeneracy — the altered topology reaches 3.0270079656 against the family's
3.0270075702, a difference of 4e−10.

The reason is structural. At a₃ → 1 the p₂ layer at node B has vanishing weight,
so **B's normal stops mattering** and the structure degenerates onto 005's
collapsed rank-4 object. The two classes share that boundary; 016's claim that
its class "contains no member of the original family" holds in the interior and
fails in the closure, and the optimum lives in the closure.

## Outcome

`VERIFIED`, scoped to the second-order coefficient at m₁ = m₁₁, m₂ = 1/4, five
conductivity triples, this topology:

- The gap's Hessian at 004's attaining point block-diagonalises. Every
  normal-tilt direction has **zero** coupling to m₁, so the six-parameter Schur
  complement equals the one-parameter γ exactly. No tilt of this topology's
  normals can lower the 015 curve.

`VERIFIED`, scoped to σ = (1,2,5), m₂ = 1/4, m₁ ∈ {3/25, 11/100}, exact ℚ
throughout — this **refutes part of 016**:

- The infimum of 016's altered class equals the original family's minimum, not
  a value 66× worse. 016's reported numbers are an interior local minimum of its
  own class.

**Not claimed.**

- That no class beats the family. This covers tilts and a₃ at *fixed topology*.
  Higher rank, different trees, and non-laminate microstructures are untouched —
  015 lead 2 stays open for those, but the tool to settle each one is now cheap.
- That 016's headline verdict was wrong. Its `REFUTED` — "this enrichment does
  not beat the curve" — stands. What is corrected is its *magnitude* (a tie, not
  a rout) and its *explanation*.
- Anything about m₁ far below m₁₁, or about the a₃ → 1 limit as an attained
  optimum: it is an infimum, approached but not reached inside the open class.

## Why it failed / what survived

016's moral was that "the collapse is the right structure, not a limitation" and
that recovering the redundant freedom "leaves the right subspace" at a cost of
one to two orders of magnitude. The correct statement is milder and more useful:
**the extra freedom is inert at the optimum, not punished.** 016's factor of 66
measured the distance from an interior local minimum to the bound, not the cost
of the enrichment.

Two things went wrong in 016 and both are cheap to avoid next time.

1. **A grid search over a bounded box stopped short of its boundary.** The
   optimum sat at a₃ → 1, where one layer's weight vanishes. Any parameter whose
   limit *degenerates the topology* deserves an explicit boundary check, because
   that limit is exactly where a richer class re-contains the poorer one.
2. **"Strictly richer" was asserted from the interior.** Two laminate classes
   that differ only in one node's normal always meet where that node's minority
   layer vanishes. A class comparison has to be made on closures.

What survived, and is the real content: the tilt block of the Hessian is
positive semidefinite and **exactly decoupled from m₁**. Rotating a normal costs
energy at second order in the tilt with no first-order gain in the m₁ direction,
so it can never trade against a volume-fraction change. That is a much stronger
statement than 016 could reach by searching — it holds for every m₁ near m₁₁
simultaneously, at every triple checked, exactly — and it is the first genuine
sharpness result for the 015 curve.

Reusable: `tp_enrich.py` — add any parameter that leaves the attaining structure
in the class, and read off whether its u-row entry is nonzero. That single number
decides whether the enrichment can move the curve, before any search.

## Leads generated

1. **Rank-increasing enrichments, done on the closure.** A higher-rank class
   generally contains the rank-4 base only on its boundary, so the extra
   directions are one-sided. The same Hessian machinery applies with the
   minimisation over a half-line; the question is whether any such direction has
   a nonzero u-row entry *of the right sign*. This is the sharpest remaining
   form of 015 lead 2.
2. **Different topologies that also attain at m₁₁.** 007's coated T² attains at
   f = (1/8,1/8,3/4) by a completely different route; if a variant of it attains
   at m₁₁ for the same (σ, m₂), its γ is a second, independent data point on the
   true curvature — and if it differs from the family's, the smaller one wins.
3. **Explain the decoupling.** The u-row vanishing on all five tilts is too
   clean to be accidental. 010 derived the attainment fields from interface
   continuity; the same argument should say why a rotation of any interface is
   second-order neutral in the volume-fraction direction. That would upgrade this
   from a computed identity at five triples to a structural fact.
4. Carried from 016, now answered in the negative for tilts and still open for
   the rest: enrichments that preserve the collapse.

## References

- This repo: `attempts/004`, `attempts/005`, `attempts/015`, `attempts/016`,
  `attempts/017`; `explore/tp_enrich.py`, `explore/tp_curvature.py`.
