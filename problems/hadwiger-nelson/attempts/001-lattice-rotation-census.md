# 001 — Rotations of the triangular lattice: a census, and a 4-colouring of the whole ring they generate

- **Problem:** hadwiger-nelson, `problems/hadwiger-nelson/PROBLEM.md`
- **Date:** 2026-07-27
- **Mode:** blind (source-commit `e5af9e6285ea9af0b5dba21adafe402d588df83d`)
- **Type:** computational search for an object, with an exact-arithmetic
  formalization and two proved barriers
- **Target shape:** object (a finite point set a checker can settle), so the
  rule was to build the checker first and then keep pushing while a
  structurally different construction idea remained. **Push rounds: 9.**
- **Tools:** written for this attempt, all pure-Python standard library, all
  deterministic.
  - `harness/hadwiger-nelson/exact_field.py` — multiquadratic exact arithmetic
  - `harness/hadwiger-nelson/exact_field_poly.py` — independent reference implementation of the same
  - `harness/hadwiger-nelson/unit_distance.py` — two-sided unit-distance certifier
  - `harness/hadwiger-nelson/colouring.py` — four independent exhaustive colouring algorithms
  - `harness/hadwiger-nelson/colouring_fast.py` — the scaling solver (k-core, components, propagation, heuristics)
  - `harness/hadwiger-nelson/selftest.py` — the cross-checks between all of them
  - `problems/hadwiger-nelson/explore/lattice.py` — lattice patches, exact
    rotations, compass intersections. Route-specific, so deliberately *not* in
    the tier-0 harness even though it is reusable: `PROBLEM.md` asks the harness
    for the field arithmetic, the certifier and the colouring checker, and this
    is none of the three.
  - `problems/hadwiger-nelson/explore/*.py` — eighteen route-specific scripts
  - Runtimes are quoted per command. The whole set re-runs in about an hour;
    the arena sweeps are the slow part.
- **Sources:** `problems/hadwiger-nelson/PROBLEM.md` only, whose own figures are
  marked machine-transcribed `[T]` there. **No paper was consulted and no
  literature was read**; this was a blind run with no network use. In
  particular the Moser spindle here was rebuilt from the lattice rather than
  transcribed from a coordinate list, and I make no claim about whether any
  other graph below is already known.

## Approach

The verification contract in `PROBLEM.md` names this problem's specific failure
mode — a float turning a non-edge into an edge — so the first deliverable had to
be exact arithmetic in a named number field, and the first result had to be a
published object rebuilt from scratch inside it.

**Why this family rather than the obvious alternatives.**

- *Transcribe and re-verify a known 5-chromatic graph.* Rejected: blind mode
  means I do not have the coordinate list, and re-verifying a published object
  with my own code tests my code, not the problem.
- *Search the space of point configurations directly.* A search over a
  continuum. Any discretisation reintroduces the floating-point failure the
  contract forbids.
- *Rotations of the triangular lattice.* The smallest family that is
  **calibratable** (it provably contains a 4-chromatic graph, so a wrong
  implementation shows up at once), **censusable** (one integer parameter), and
  **exact** (every coordinate lands in a named quadratic field).

The construction. Write the triangular lattice as the Eisenstein integers
E = Z + Z·ζ, ζ = exp(iπ/3), so the lattice point (a,b) = a + bζ has squared
length the integer norm N(a,b) = a² + ab + b². For a positive integer n let ρ_n
be the rotation about the origin with

    cos θ = 1 − 1/(2n),    sin θ = √(4n−1) / (2n),

the unique angle (up to sign) carrying a point at distance √n from the origin to
a point at distance exactly 1 from where it started. Its entries lie in
Q(√(4n−1)), so P ∪ ρ_n(P) lives in Q(√3, √(4n−1)) for any lattice patch P.

The attempt ran in two phases. The first censused that family and measured why
it stalls. The second — after the target was re-framed as an **object** a
checker can settle, with no ceiling on the claim — went after a 5-chromatic
graph directly, and ended up proving that the natural arena for it cannot
contain one.

## What was done

### 0. Environment

`where gcc`, `where cc`, `where clang`, `where tcc` all return nothing on this
machine, and there is no SAT solver and no nauty. So a DRAT/LRAT certificate
was not available. The contract's stated alternative — "an exhaustive search
reproduced by a second, independently written implementation" — is what every
non-colourability claim below rests on, with the second implementation named
each time.

### 1. Exact arithmetic, and its independent check

`exact_field.py` implements Q(√d₁,…,√d_k) for pairwise coprime squarefree
generators, an element being an integer coordinate vector over the subset basis
with a common denominator. Those 2^k numbers are linearly independent over Q,
so equality is a coordinate-wise integer comparison and a distance decision is
an integer comparison. Nothing branches on a float; `to_float` and `to_decimal`
exist only for reporting.

`exact_field_poly.py` is a second implementation written separately: a dict of
monomials with `Fraction` coefficients, generic polynomial multiplication, then
reduction modulo (x_i² − d_i). Different data structure, different algorithm,
no shared code.

```
python harness/hadwiger-nelson/selftest.py        # ~40 s, ALL PASS
```

checks 8000 random operations of the fast path against the reference across
five fields, plus a 140-digit `Decimal` evaluation as a third net (worst
disagreement 1.0e-136, i.e. rounding), plus exhaustive inversion.

**This immediately earned its cost.** The first run reported 3384 mismatches
out of 8000: `__add__` was forming the common denominator as `a*lb` instead of
`a*la`. Every distance decision taken before that fix would have been silently
wrong, and nothing in the outputs would have looked odd.

### 2. Certifier

`unit_distance.build_graph` tests **every** pair exhaustively and exactly — no
numeric prefilter, so the non-edges are certified by construction, not by
omission. `certify` recomputes all of it and compares in both directions, and
verifies the points are pairwise distinct as field elements.
`numeric_separation` reports, at 80 digits, the largest |d²−1| over the edges
and the smallest over the non-edges.

### 3. Colouring algorithms

`colouring.py` implements (1) DSATUR-ordered backtracking; (2) unpruned
enumeration of kⁿ; (3) the chromatic polynomial by deletion/contraction; (4) an
inclusion–exclusion oracle over independent-set counts,
cov_k(G) = Σ_{S⊆V} (−1)^{|V|−|S|} i(S)^k, nonzero exactly when G is
k-colourable, which never constructs a colouring. All four agree on χ on 60
random graphs, and (3) and (4) agree on the verdict for every k.

A correction worth recording: routine (4) was first written as a *count of
proper colourings*, which it is not — it counts covers by independent sets. The
selftest caught the discrepancy against the chromatic polynomial. It is correct
as a colourability oracle and is used only as one.

### 4. Calibration: rebuild a 4-chromatic graph from the lattice

```
python problems/hadwiger-nelson/explore/calibrate_spindle.py     # <1 s
```

The lattice patch of norm ≤ 3 (13 points) and its image under ρ₃ (cos 5/6,
sin √11/6, orthogonality checked exactly) give 25 distinct points in
Q(√3,√11), all 300 pairs decided exactly, 54 edges, χ = 4. Greedy
vertex-critical extraction gives

- **7 vertices, 11 edges, degree sequence [3,3,3,3,3,3,4]**
- P(G,k) = k⁷ − 11k⁶ + 51k⁵ − 129k⁴ + 188k³ − 148k² + 48k, so **P(G,3) = 0**
  and P(G,4) = 384
- χ = 4 by DSATUR, brute force over 3⁷, the chromatic polynomial, and the
  inclusion–exclusion oracle
- smallest |d²−1| over a non-edge: 0.4574… at 80 digits

That is the Moser spindle, recovered rather than transcribed.

### 5. A proved barrier: what a unit compass alone can build

Operation: given two existing points at distance strictly between 0 and 2, add
the two points at distance exactly 1 from both.

**Lemma.** The closure of a unit segment under this operation is contained in
the Eisenstein lattice E, whose unit-distance graph has χ = 3.

*Proof.* The only Eisenstein norms strictly between 0 and 4 are 1 and 3, so
those are the only squared distances at which the operation applies. E is
invariant under translation by E and multiplication by ζ, which act
transitively on the ordered pairs of each norm, so one pair of each suffices:
for (0, 1) the intersections are ζ and 1−ζ = ζ̄, both in E; for (0, 1+ζ) they
are 1 and ζ, both in E. Hence E is closed, and since {0,1} ⊆ E the closure is
contained in E. For the colouring, let λ = 1 + ζ, of norm 3; reduction modulo λ
sends ζ to −1, so it is (a,b) ↦ (a−b) mod 3, and adjacent lattice points differ
by a unit, never 0 mod λ. Three colours suffice, and a unit triangle forces
three. ∎ *(That the closure is all of E I did not prove; the computation reaches
all 55 lattice points inside |p|² ≤ 13.)*

```
python problems/hadwiger-nelson/explore/compass_closure.py 13      # <1 s
```

55 points after 9 rounds; the only squared distances below 4 ever encountered
are 1 and 3; zero field escapes; every point lands on the lattice; the graph
(138 edges) is 3-chromatic and (a−b) mod 3 is verified proper by the independent
linear checker. The closure step was separately checked on all 420 intersecting
pairs of the norm ≤ 21 patch, with zero non-lattice results.

**Consequence.** A unit compass can never beat three colours, so every
4-chromatic construction must import a non-unit radius. The spindle imports
exactly one: the far tip of the second rhombus is pinned by a unit circle *and a
circle of radius √3*, after which its two neighbours are ordinary unit-circle
intersections.

### 6. The census

G(n,m) = P(m) ∪ ρ_n(P(m)), P(m) the lattice patch of norm ≤ m.

```
python problems/hadwiger-nelson/explore/rotation_census.py \
       --max-n 79 --patches 79 --copies 2       # ~4 min, 30 rows
```

Over all 30 Eisenstein norms n ≤ 79 at m = 79 (graphs of 295–589 vertices,
every pair decided exactly):

> χ(G(n,79)) = 4 when 3 | n, and 3 otherwise — for every n tested, with no
> exceptions.

χ = 4 at n = 3, 9, 12, 21, 27, 36, 39, 48, 57, 63, 75; χ = 3 at
n = 1, 4, 7, 13, 16, 19, 25, 28, 31, 37, 43, 49, 52, 61, 64, 67, 73, 76, 79.
Same verdicts at m = 48 and m = 28 for the n those patches cover.

The invariant behind the split is exact: mod 3, a² + ab + b² ≡ (a−b)², so
**3 | N(p) ⟺ p ≡ 0 mod λ ⟺ p has colour 0** in the canonical 3-colouring — that
is, "3 | n" says every lattice point of norm n carries the same colour. Checked
exhaustively over |a|,|b| ≤ 200 (160801 points, 0 counterexamples).

*Why one colour class being the source of all the built-in edges should force a
fourth colour is `SPECULATION`.* I have the correlation across 30 values of n
and the algebraic characterisation of the split, and no proof connecting them.

### 7. Four certified 4-critical graphs

```
python problems/hadwiger-nelson/explore/verify_cores.py 21 21     # ~6 s
```

| n | field | core | edges | degree sequence | worst non-edge gap |
|---|---|---|---|---|---|
| 3 | Q(√3,√11) | 7 | 11 | 3⁶ 4 | 0.457427 |
| 12 | Q(√3,√47) | 13 | 21 | 3¹⁰ 4³ | 0.119764 |
| 9 | Q(√3,√35) | 17 | 29 | 3¹² 4³ 5² | 0.069275 |
| 21 | Q(√3,√83) | 19 | 33 | 3¹⁰ 4⁹ | 0.037130 |

Each carries, in `data/verified-cores.json`: exact coordinates; a two-sided
geometry certificate over all pairs; 80-digit separation figures;
non-3-colourability by DSATUR **and** the inclusion–exclusion oracle (plus brute
force and the chromatic polynomial where size allows); an explicit 4-colouring
checked by the linear checker; and vertex-criticality, verified by deleting each
vertex in turn. Being vertex-critical and larger than 7, **none of the 13-, 17-
and 19-vertex cores contains a Moser spindle**.

### 8. Phase-one pushes, and the interface bound

Four structurally different pushes, all landing on 4: powers of one rotation
(217 vertices); many rotation centres and both directions (**751 vertices, 3796
edges, 2140 cross-copy** — still χ = 4); two rotation families at once in
Q(√3,√11,√35) (289 vertices); and the commensurable case 4n−1 = 3m², where the
union collapses into one triangular lattice scaled by 1/√n and stayed
3-chromatic in all 15 configurations checked.

The measured obstruction was that **the interface between a patch and its
rotated image is bounded independently of the patch** — 12 cross edges at patch
norms 3, 7, 13, 21 and 28 alike while the patch grows from 13 to 109 points —
with a proof and an integer re-derivation given in the obstruction section
below.

### 9. A scaling solver, and the asymmetry that makes large graphs tractable

`colouring_fast.py` adds k-core reduction (a vertex of degree < k is coloured
last, so G is k-colourable iff its k-core is), component splitting, bitmask
domains with unit propagation, and DSATUR order with colour-symmetry breaking.

The two directions are handled asymmetrically **on purpose**: a colouring is
self-certifying, verified in linear time, so the YES direction may be answered
by heuristics; only the NO direction needs the complete search. `decide_k_colourable`
therefore runs the complete search under a small node budget first (with
propagation it is itself a strong heuristic, and it *settles* easy NO instances),
then randomised greedy plus min-conflicts, then the complete search under the
full budget. Reordering these gave a ~100× speedup on the gadget searches.

Every one of these is cross-checked in `selftest.py` against the four
exhaustive algorithms: verdicts for k = 1..5 on 400 random graphs, the
properness of every colouring returned, agreement between reduced and
unreduced runs, and the semantics of the contraction helper against explicit
enumeration of all colourings.

### 10. The right arena: the ring Z[ζ][ρ], and its ladder

ρ = ρ₃ satisfies 3ρ² = 5ρ − 3, so ρ is not an algebraic integer and
Z[ζ][ρ] is not a finitely generated Z[ζ]-module. It is the increasing union

    H_0 ⊂ H_1 ⊂ H_2 ⊂ …,    H_k = 3^(−k) (E + ρE).

Scaling by 3^k is a similarity, so **the unit-distance graph on H_k is the
Cayley graph on the rank-4 group E ⊕ E whose generators are the vectors of
squared length 9^k**. Same vertex group at every level; a strictly richer
generator set as k grows.

Two independent enumerations of those generators agree exactly:

- a sweep of a coefficient box, deciding each candidate in Q(√3,√11);
- an arithmetic classification: |Z + ρW|² rational forces Z̄W real, hence W ∥ Z,
  hence Z = bY, W = aY with gcd(a,b)=1 and N(Y)·Q(a,b) = 3·9^k for
  Q(a,b) = 3a² + 5ab + 3b². Since Q(a,b) ≥ (a²+b²)/2 the search is finite.

They agree at levels 0–3 (18, 42, 66, 90 generators — the box sweep takes 23 s
at level 3, the arithmetic 0.00 s), and the arithmetic one reaches level 8. The
count is **24k + 18**.

### 11. Searching inside the ring

Cayley balls V_L(k) are the *maximal* candidate at radius L: any connected
unit-distance graph on ring points, of diameter ≤ L and containing the origin,
is a subgraph. Edges come from the generator table and were cross-checked
against the pair-by-pair exact certifier on prefixes (agreement on every case
tested).

| level | L | vertices | edges | mean degree | 4-colourable |
|---|---|---|---|---|---|
| 0 | 3 | 715 | 3906 | 10.93 | yes |
| 0 | 4 | 2113 | 13296 | 12.58 | yes |
| 1 | 2 | 883 | 3768 | 8.53 | yes |
| 1 | 3 | 10315 | 77016 | 14.93 | yes |
| 2 | 2 | 2179 | 9192 | 8.44 | yes |

Coefficient boxes at level 0 reached 1849 vertices and mean degree 12.4, also
4-colourable. Gadget searches found nothing either: over 144 non-adjacent pairs
of the level-0 depth-2 ball, no pair is forced equal or forced unequal by
4-colourings, and every colour pattern on a distance-√3 terminal triple of the
883-vertex level-1 ball is realizable.

### 12. Why all of that had to fail: a 4-colouring of the whole ring

```
python problems/hadwiger-nelson/explore/ring_four_colouring.py    # ~5 s
```

**Theorem.** Write a point of H_k as (Z + ρW)/3^k with Z, W ∈ E, and set

    c(p) = (λZ + W)  mod 2E,     λ = 1 + ζ.

E/2E has four elements, so c is a 4-colouring, and it is proper on the
unit-distance graph of the entire ring Z[ζ][ρ].

*Well defined.* The inclusion H_k ↪ H_{k+1} is (Z,W) ↦ (3Z,3W) and 3 ≡ 1
mod 2, so the value is unchanged. One rule colours every level at once.

*Proper.* Let (Z,W) be a unit vector, so |Z + ρW|² = 9^k, odd. Note
|Z + ρW|² = N(Z) + N(W) + (5/3)·r always lies in (1/3)Z. Suppose
λZ + W ≡ 0 (mod 2E). N(λ) = 3 is odd, so λ is invertible in E/2E. By the
parallelism lemma bW = aZ in E for coprime a, b (the cases Z = 0 and W = 0 are
immediate: λ and 1 are invertible mod 2). Reduce mod 2E and split on parities,
which cannot both be even:
 – a, b both odd ⟹ W̄ = Z̄ ⟹ (λ+1)Z̄ = ζZ̄ = 0 ⟹ Z̄ = 0;
 – a even, b odd ⟹ W̄ = 0 ⟹ λZ̄ = 0 ⟹ Z̄ = 0;
 – a odd, b even ⟹ Z̄ = 0 ⟹ W̄ = 0.
So Z ≡ W ≡ 0 (mod 2E), hence Z + ρW = 2(Z′ + ρW′) and
9^k = 4·|Z′ + ρW′|² ∈ 4·(1/3)Z — impossible, since 9^k has no factor 2. ∎

**Consequence.** Every finite unit-distance graph whose points lie in
Z[ζ][ρ] is 4-colourable. The Moser spindle lies in H_0, so the ring's
unit-distance graph is **exactly 4-chromatic**, and no construction built from
the triangular lattice and the spindle rotation — to any depth, with any number
of rotated copies, at any level of the ladder — can ever witness χ(R²) ≥ 5.

Verification, none of which reuses the proof: the colour map kills no unit
vector at any level 0–8 (all 18…210 of them, from the arithmetic enumeration);
and applying the rule to real balls built by the group machinery gives zero
monochromatic edges at level 0 depth 4 (2113 v), level 1 depth 3
(**10315 vertices, 77016 edges**), and levels 2–5 depth 2 (up to 9523 v).

### 13. Which arenas the argument does *not* close

The same idea generalises: for a group G = Σ_j r_j E of rotated lattices,
ψ(Σ r_j Z_j) = Σ μ_j Z_j mod 2E is a 4-colouring, proper exactly when no unit
vector lies in ker ψ. Two preconditions are checked, not assumed: that G is a
direct sum (otherwise a coefficient-wise rule is not a function of the point at
all), and that the unit vectors are enumerated completely for the box used.

```
python problems/hadwiger-nelson/explore/mod2_cap_sweep.py --max-n 21
```

39 groups, one and two rotations. The injectivity check earns its keep: every
group containing a *commensurable* rotation (4n−1 = 3m², i.e. n = 7, 19, 37, …)
fails it, because ρ_n then lies in Q(√3) and the sum is not direct — those rows
are excluded rather than believed. Of the rest, **the only groups with no
capping functional are those containing ρ₄**, whose sine is √15/8 — denominator
a power of 2, exactly the modulus the argument uses. And ρ₄ alone is then closed
at a different modulus: the quotient (E/3E)² is 4-colourable (81 vertices, 729
edges), so E + ρ₄E is capped too.

That leaves the mixed groups E + ρ₃E + ρ₄E and relatives as the arenas still
open at the end of this attempt; the object hunt there is Lead 1.

## Outcome

- `VERIFIED` — **the ring Z[ζ][ρ] is exactly 4-chromatic.** Explicit colouring
  c = (λZ + W) mod 2E, proof in §12, plus: the colour map kills no unit vector
  at levels 0–8, and zero monochromatic edges on real balls up to 10315
  vertices / 77016 edges. Scope: every finite unit-distance graph with all
  points in Z[ζ][ρ]. This is the headline result and it is a *negative* one.
- `VERIFIED` — **the generator classification.** Every unit vector of H_k has
  parallel components and there are exactly 24k + 18 of them. Two independently
  written enumerations (field box sweep vs arithmetic over Z) agree exactly at
  levels 0–3.
- `VERIFIED` — **the compass barrier.** The closure of a unit segment under
  unit-circle intersection lies in the Eisenstein lattice, whose unit-distance
  graph is 3-chromatic. Proved; checked inside |p|² ≤ 13 and on all 420
  intersecting pairs of the norm ≤ 21 patch.
- `VERIFIED` — **the interface bound.** For non-commensurable n, every
  cross-copy unit edge (p, ρ_n q) has p ∥ q and N(p) ≤ 4n²/(4n−1) < n+1.
  Cross-checked against the exact geometry on 112 (n, patch) pairs, zero
  disagreements, zero bound violations.
- `VERIFIED` — **a 7-vertex, 11-edge 4-chromatic unit-distance graph** in
  Q(√3,√11), rebuilt from the lattice; χ = 4 by four independent algorithms;
  P(G,3) = 0, P(G,4) = 384. The Moser spindle: χ(R²) ≥ 4 re-derived, not
  improved.
- `VERIFIED` — **four 4-vertex-critical unit-distance graphs** on 7, 13, 17 and
  19 vertices in four quadratic fields, coordinates shipped,
  non-3-colourability confirmed by two independent complete algorithms each.
- `EVIDENCE` — **the census rule** χ(P(m) ∪ ρ_n P(m)) = 4 iff 3 | n, over all
  30 Eisenstein norms n ≤ 79 at patch norm 79, reconfirmed at 48 and 28.
- `EVIDENCE` — **the capping sweep.** Of 39 rotation groups, every one that is a
  direct sum is 4-capped by a mod-2 functional except those containing ρ₄; ρ₄
  alone is capped mod 3.
- `REFUTED` (as a route, by the theorem above) — **any construction inside
  Z[ζ][ρ]**. Not to be re-attempted as stated: enlarging patches, stacking
  rotated copies, going deeper in the ladder, or searching bigger Cayley balls
  are all provably incapable of reaching 5.

**Not claimed.** No 5-chromatic unit-distance graph, and no lower bound beyond
the classical 4. Nothing about χ(R²) itself: the theorem says the *ring* is
4-colourable, which is a statement about one countable point set, not about the
plane — the plane contains points outside every such ring. No claim that the
13-, 17- or 19-vertex cores are new; blind mode, no literature consulted. No
claim that the census rule holds beyond n ≤ 79, and none that the mod-2 sweep is
exhaustive over rotations.

## Why it failed / what survived

Two obstructions, the second of which supersedes and explains the first.

### The measured obstruction: a bounded interface

Write p = a₁ + b₁ζ, q = a₂ + b₂ζ, N₁ = N(a₁,b₁), N₂ = N(a₂,b₂),
2D = 2a₁a₂ + a₁b₂ + a₂b₁ + 2b₁b₂ and K = a₂b₁ − a₁b₂. Then

    |p − ρ_n q|² = N₁ + N₂ − ((2n−1)/n)·D − (√(3(4n−1)) / (2n))·K.

When 3(4n−1) is not a perfect square the last term is irrational unless K = 0,
so a unit distance forces **two integer conditions**: K = 0, and
2n(N₁+N₂) − (2n−1)(2D) = 2n. K = 0 says p ∥ q, so q = μp and
N₁(1 + μ² − 2cμ) = 1 with c = (2n−1)/(2n); since 1 + μ² − 2cμ ≥ 1 − c²,

    **N₁ ≤ 4n²/(4n−1) < n + 1.**

Every cross edge has both endpoints inside the disc of squared radius n,
*whatever the patch is*. The constraint count is Θ(1) in the patch while the
vertex count is Θ(area). `interface_theory.py` recomputes every cross-edge count
from the integer criterion alone — no field arithmetic — and agrees with the
exact geometry on all 112 (n, patch) pairs. Making √(3(4n−1)) rational is the
only escape and is self-defeating: the union then collapses into one triangular
lattice scaled by 1/√n and stayed 3-chromatic in all 15 cases checked.

### The real obstruction: the arena is 4-colourable

Adding rotation centres does break the interface saturation — 2140 cross edges
at 751 vertices — so the phase-one obstruction was not the end of the story. The
end of the story is §12: **every point that the lattice and the spindle rotation
can generate lies in Z[ζ][ρ], and that whole ring is 4-colourable.** The step
that breaks is not any particular construction; it is the choice of arena.

The quantity that goes the wrong way is *rigidity*. Lemma L1 — |Z + ρW|²
rational forces W ∥ Z — says the ring's unit vectors are far scarcer than the
ambient plane's: there are only 24k + 18 at level k, all lying on lines through
the origin in coefficient space. That scarcity is exactly what lets a single
F₄-linear functional avoid all of them. A construction that beats four colours
needs an arena whose unit vectors are *not* confined to that few directions, and
one rotation cannot supply it — no matter how many copies, centres or ladder
levels are stacked, because all of them stay inside the same ring.

This retro-explains every negative result in the attempt: the census, the
depth sweeps, the multi-centre sweeps, the Cayley balls to 10315 vertices, the
absence of forcing pairs, and the fully realizable pattern profiles. None of
them could have succeeded.

### What survived / is reusable

- The harness: exact multiquadratic arithmetic with an independent reference
  implementation, a two-sided certifier, four cross-checked exhaustive colouring
  algorithms, and a scaling solver whose YES/NO asymmetry is what makes
  10⁴-vertex instances practical.
- **The capping technique itself**, which is the most transferable thing here: to
  rule out a whole family of constructions at once, exhibit a group
  homomorphism from the coefficient group under which no unit vector vanishes.
  It replaces an unbounded search with a finite check, and it is what turned an
  inconclusive search into a theorem.
- The classification of the ring's unit vectors (parallel components, 24k + 18)
  and its two independent enumerations.
- The compass-closure lemma; the interface bound and its integer criterion; the
  mod-λ observation 3 | N(p) ⟺ colour 0; four certified 4-critical graphs.
- The arena sweep, which says exactly where searching can still pay.

### Methodological lesson

The two-implementation rule paid for itself twice, and neither failure was
visible from the outputs: a wrong common denominator in field addition (caught
by an independently written reference implementation) and an
inclusion–exclusion identity mislabelled as a colouring count (caught by
cross-checking against the chromatic polynomial). Both produced entirely
plausible numbers.

The second lesson is about search versus structure. Eight of the nine push
rounds were searches, and all of them returned "4-colourable" without saying
why. The ninth was a five-line algebraic check, and it settled every one of them
at once — and would have settled them before any of the searches were run. When
a search family keeps returning the same answer, the next move is to look for
the invariant that forces it, not for a bigger instance.

## Leads generated

1. **Hunt in the arena the sweep leaves open.** E + ρ₃E + ρ₄E is a direct sum
   (checked), has 30 unit vectors, and admits no mod-2 capping functional; the
   quotient (E/3E)³ did not settle under the node budget. Build Cayley balls
   there and decide 4-colourability, and separately decide whether *some*
   quotient (E/NE)³ is 4-colourable. Definite outcome either way: a capping
   colouring closes the arena, and a non-4-colourable ball is a 5-chromatic
   unit-distance graph.
2. **Characterise the capping obstruction.** Conjecture: a group Σ r_j E admits
   a mod-2 capping functional unless some ρ_n in it has an even denominator
   (4n−1's rotation written as (2n−1)/(2n) with 2n a power of 2 times …). Test
   by extending `mod2_cap_sweep.py` to all n ≤ 200 and to three-rotation groups,
   and check whether "contains ρ₄, ρ₁₆, ρ₆₄, …" predicts the failures exactly.
3. **Generalise off the diagonal.** The rotation family only ever asks
   |p − ρ_θ p| = 1. For distinct lattice points p, q of norms m, n a unit
   distance needs cos(θ + φ) = (m + n − 1)/(2√(mn)), real iff
   (m − n)² ≤ 2(m + n) − 1. Enumerate the admissible (m, n) with m, n ≤ 50 and
   ask whether the resulting rings are still mod-2 capped.
4. **Escape the ring by leaving the group.** Compass intersections introduce
   halves and new square roots, so the constructible closure is strictly larger
   than Z[ζ][ρ] — and mod 2E is exactly what halving destroys. Compute the
   closure of the Moser spindle under unit-circle intersections *within*
   Q(√3,√11), which needs a square-root routine for multiquadratic fields, and
   test whether the capping colouring extends to it.
5. **Settle the census rule.** Prove or refute χ(P(m) ∪ ρ_n P(m)) = 4 iff 3 | n
   for m ≥ n. A refutation needs one n with 3 | n and χ = 3, or 3 ∤ n and χ = 4.
6. **Are the new cores minimal, and distinct?** Run the criticality extraction
   from many vertex orders on G(12, 21) and take the minimum; test the four
   cores for pairwise isomorphism. Definite outcome: a number for the smallest
   4-chromatic subgraph of G(12, 21), currently 13.
7. **Cap the commensurable case properly.** For n = 3k²+3k+1 the union sits in
   one triangular lattice; compute χ of the *infinite* distance-n graph on the
   triangular lattice via its quotient by NE for small N, which caps every
   finite patch at once.

## References

- `problems/hadwiger-nelson/PROBLEM.md` — the statement and the published
  figures, themselves marked machine-transcribed `[T]` there. The only source
  consulted.
- No papers were read. `PROBLEM.md` cites de Grey (2018) for the 1581-vertex
  5-chromatic graph and Parts, arXiv:2010.12665, for the 509-vertex record;
  **neither was consulted**, and nothing here depends on them beyond knowing
  that the classical range is 4 ≤ χ(R²) ≤ 7 and that the lower bound is now 5.
  That published 5 is consistent with everything above: it says a 5-chromatic
  unit-distance graph exists, and §12 says it cannot live in Z[ζ][ρ].
- No prior attempt records in this repository were read: blind run from
  source-commit `e5af9e6285ea9af0b5dba21adafe402d588df83d`.
- Evidence files: `problems/hadwiger-nelson/data/` — `calibration-spindle.json`,
  `compass-closure.json`, `rotation-census-n48.json`, `rotation-census-n79.json`,
  `verified-cores.json`, `interface-growth.json`, `interface-theory.json`,
  `commensurable-collapse.json`, `multi-centre-sweep.json`,
  `spindling-depth.json`, `mixed-rotations.json`, `ring-four-colouring.json`,
  `ball-search.json`, `mod2-cap-sweep.json`.
