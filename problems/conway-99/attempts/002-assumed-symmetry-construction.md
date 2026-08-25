# 002 — Assumed-symmetry construction: orbit matrices, and the order-33 case

- **Problem:** conway-99, `problems/conway-99/PROBLEM.md`
- **Date:** 2026-07-27
- **Mode:** blind
- **Source-commit:** `e5af9e6285ea9af0b5dba21adafe402d588df83d`
- **Continues:** `attempts/001-forced-pair-model-and-partner-regular-split.md`,
  by the same agent in the same cycle. 001 stopped at a precisely-described
  obstruction; this record is what happened when that stopping point was
  overruled and the attempt was pushed at **construction** instead. Everything
  001 established — the exact feasibility harness, the forced pair model and
  its bijection, the partner-regular split, lemmas L1–L4, and the finding that
  the propagating search fails its own positive control — is **assumed here and
  not restated**.
- **Type:** attempted construction under assumed symmetry + orbit-matrix
  elimination
- **Target shape:** object. A 99×99 adjacency matrix is the entire answer and
  `harness/conway-99/srg.py::check_srg` settles it in microseconds, so the rule
  followed here was: build the checker first, then do not stop at a partial
  result while a structurally different idea remains. **Push rounds: 5.**
- **Tools:** standard-library Python 3.13, deterministic. New in this attempt,
  all under `problems/conway-99/explore/`:
  - `orbit_matrix.py` — orbit matrices for a semiregular group, with isomorph
    rejection; `orbit_general.py` — the same for mixed orbit sizes
  - `orbit_lift.py` — lifting an orbit matrix to a graph; `orbit_calibrate.py`,
    `lift_calibrate.py` — the calibration harnesses for both
  - `m33_exhaustive.py`, `m33_verify.py` — the order-33 elimination and its
    independent verification
  - `order3_lemma.py` — lemma L5 and its checks
  - `srg_anneal.py`, `pair_anneal.py`, `cyclic_anneal.py`, `group_anneal.py` —
    four construction searches; `ga_control_243.py` — the positive control
  - `make_symmetry_evidence.py` — regenerates `data/symmetry-evidence.json`
  - `tests/test_conway99_orbits.py` — 19 tests, ~4 s
  - Still no C compiler, `nauty` or SAT solver on this machine. No SAT encoding
    was used and no UNSAT claim appears below.
- **Sources:** as 001 — only `problems/conway-99/PROBLEM.md`, machine
  transcribed `[T]`. No papers consulted; L5 below is an independent
  re-derivation.

## Approach

001 concluded with a barrier and stopped. That was the wrong stopping rule for
an object-shaped target: when a checker settles the object instantly there is no
gap between "found it" and "proved it", so a construction attempt should keep
going while distinct structural ideas remain.

The route chosen is the standard one for finding graphs like this: **assume an
automorphism and search the quotient.** The published restrictions make the
candidate groups tiny, and 001's own lemmas L1–L3 already pin most of the orbit
profiles, so the assumption costs little and collapses the space by orders of
magnitude. Concretely, for orbits O_i of sizes n_i, projecting
A² = kI + λA + μ(J−I−A) onto orbits gives

    Σ_l r_il·r_lj = k·δ_ij + λ·r_ij + μ·(n_j − δ_ij − r_ij),   r_ij·n_i = r_ji·n_j

which is exact integer arithmetic and enumerable when the number of orbits is
small.

The methodological rule from 001 was kept and tightened: **every engine had to
rebuild a graph that exists before its output on 99 was read.** That is what
produced the one result here and what disqualified everything else.

## What was done

### Round 1 — orbit matrices, and the semiregular Z₃₃ elimination

`orbit_matrix.py` enumerates the conditions above for a semiregular group, with
isomorph rejection by a non-increasing diagonal (validated to preserve every
canonical class, not merely to shrink the output). Calibration first: the true
orbit matrices of Paley(9)/Z₃, rook(3)/Z₃, rook(4)/Z₄ and BvLS(243)/Z₃ all
satisfy the conditions, and for the small cases the enumeration provably
*contains* them.

**Semiregular Z₃₃ (three orbits of 33) has exactly one orbit matrix**,
R = [[2,6,6],[6,2,6],[6,6,2]] — re-derived independently by brute force over
every symmetric 3×3 candidate with even diagonal and row sums 14.

The lift then closes completely, and the mechanism is worth stating because it
generalises. Write A_U(y) for the autocorrelation of U ⊆ Z₃₃; it is symmetric,
A_U(−y) = A_U(y), since swapping a pair negates the difference. The three
*diagonal* strongly-regular equations are

    A_{S00} + A_{S01} + A_{S02} = f₀
    A_{S01} + A_{S11} + A_{S12} = f₁      (using S₁₀ = −S₀₁)
    A_{S02} + A_{S12} + A_{S22} = f₂

— three linear equations in the three unknown functions, with an invertible
coefficient matrix. So the three diagonal sets S_ii = {a,−a} **force** every
off-diagonal autocorrelation:

    A_{S01} = (f₀+f₁−f₂)/2,  A_{S02} = (f₀−f₁+f₂)/2,  A_{S12} = (−f₀+f₁+f₂)/2

Over all 16³ = 4096 diagonal choices, **4050 die immediately** on a
half-integer autocorrelation, and the surviving **46** demand autocorrelation
vectors that **no 6-subset of Z₃₃ possesses**.

Verified two independent ways in `m33_verify.py`: the forcing formula
reproduces the *true* autocorrelations of Paley(9) under its Z₃ (same t = 3
shape); and the unrealisability is re-derived by brute force over all
C(32,5) = 201376 6-subsets containing 0, tabulating the 16778 autocorrelations
that actually occur. Translation invariance, which is what makes "containing 0"
lossless, is checked rather than assumed.

### Round 2 — mixed orbit sizes, and lemma L5

A semiregular action is not the only way order 33 can act, and assuming
otherwise would be exactly the incomplete-case-split failure the verification
contract warns about. `orbit_general.py` handles mixed orbit sizes; it is
calibrated on Petersen/Z₃ (sizes [3,3,3,1]) and Paley(9)/Z₂ (sizes [1,2,2,2,2]),
whose real orbit matrices it contains — without that calibration a count of zero
would prove nothing.

By 001's L2, an order-33 σ can have no orbit of size 1 or 3 (the stabiliser
would contain the order-11 subgroup, giving an order-11 element with a fixed
point). So 99 = 33a + 11b with a ≥ 1: three profiles. Enumerating them:

| profile | orbit matrices |
|---|---|
| 3×33 (semiregular) | 1, eliminated above |
| 1×33 + 6×11 | **0** (963673 nodes) |
| 2×33 + 3×11 | 6, **never lifted** |

Closing the third needed a new lemma.

**L5 — an order-3 automorphism of an SRG(99,14,1,2) has exactly 0 or 3 fixed
points.** For v ∈ Fix(τ), τ permutes the 7 matching edges of N(v) in orbits of
size 1 or 3, so the fixed-edge count is ≡ 7 ≡ 1 (mod 3), i.e. 1, 4 or 7; a
τ-fixed edge has both endpoints fixed (odd order), and 7 would fix N(v)
pointwise and force τ = id by 001's L1. Hence |Fix ∩ N(v)| ∈ {2,8}. Fix is
closed under common neighbours (odd order fixes a set of size ≤ 2 pointwise),
and for any induced F with that closure,

    Σ_{x ∈ N_F(u)} deg_F(x) = deg_F(u)·(1+λ−μ) + μ·(f−1)  =  2f − 2   for (λ,μ)=(1,2)

for every u ∈ F. Degree 2 forces f ≤ 9; degree 8 forces f ≥ 9; mixed degrees
force f = 9 with degree-8 vertices adjacent only to degree-2 ones and vice
versa, which λ = 1 kills (the third vertex of a triangle on a (8,2) edge would
need both degrees). So F is regular: degree 2 gives f = 3, a triangle; degree 8
gives f = 33, i.e. an SRG(33,8,1,2), which the harness rejects on multiplicity
integrality. Hence **f ∈ {0,3}**.

Applied to order 33: the 11-orbit points are exactly Fix(σ¹¹), and σ¹¹ has order
3, so their number 11b must be 0 or 3 — forcing b = 0.

Two further consequences, pure arithmetic on top of L5: **no automorphism of
order 27** (neither 99 nor 96 is divisible by 27), and **an order-9 automorphism
is semiregular** with 11 orbits of 9 (96 is not divisible by 9).

L5's checks: the counting identity and the closure hold for every order-3
automorphism of Paley(9), rook(3), Petersen and rook(4). The identity is checked
in its general (λ,μ) form — the 2f−2 shape is special to λ=1, μ=2, and an
earlier version of the check that hardcoded it produced a false failure on
rook(4). And L5's k=4 analogue predicts f = 0 only, which is exactly what
|Aut| = 72 delivers for both 9-vertex graphs.

### Rounds 3–5 — four construction searches, and their ceilings

Three further structurally distinct searches were built (the fourth,
`pair_search`, is 001's):

- **`srg_anneal`** — anneal over all k-regular graphs, degree-preserving
  2-swaps, incremental energy (verified against full recompute).
- **`pair_anneal`** — the same but *inside* 001's pair model, so v₀, N(v₀) and
  every P–D edge are correct by construction and only D×D is searched.
- **`cyclic_anneal`** — anneal over Z_m-invariant graphs via connection sets
  S_ij ⊆ Z_m. Regularity is free here: c_ii(0) = Σ_l |S_il| *is* the degree.
- **`group_anneal`** — anneal over unions of orbits of an assumed group acting
  on vertex pairs. Fully general: handles fixed points, mixed orbit sizes and
  non-cyclic groups, and collapses C(99,2) = 4851 variables to 441 (Z₁₁), 693
  (Z₇ with its forced fixed point) or 1617 (Z₃).

Measured ceilings, each by rebuilding graphs that exist:

| engine | rebuilds | fails |
|---|---|---|
| `srg_anneal` | n=27 (27,10,1,5) λ=1; also 9, 15, 16 | n=25, 36, 45 |
| `cyclic_anneal` | n=27 at t=9 **in 16 s**; also 9, 10, 15, 16 | n=45 at t=9 |
| `group_anneal` | n=9, 16 | n=243 (BvLS) under an assumed order-3 map |
| `pair_search` (from 001) | n=9 | n=243: 14/220 settled after 120000 nodes |

### Where every run at 99 was stopped

All runs were terminated at wrap-up on a load decision. **None had written a
witness**, so there was no checkpoint to salvage and the stopping point is the
result. Cost is the energy defined in each module; a witness is cost 0.

| run | assumed structure | stopped at | controlled? |
|---|---|---|---|
| `cyclic_anneal` Z₁₁ | 9 orbits of 11 | 11 restarts × 400k iters, best **336** | **no** |
| `cyclic_anneal` Z₉ | 11 orbits of 9 | 8 × 400k, best **440** | **no** |
| `cyclic_anneal` Z₃ | 33 orbits of 3 | 1 × 300k, best **1448** | **no** |
| `group_anneal` Z₁₁ | 441 pair-orbits | 3 × 60k, best **2508** | **no** |
| `group_anneal` Z₇ | 693 pair-orbits, 1 fixed point | 5 × 60k, best **2653** | **no** |
| `group_anneal` Z₃ | 1617 pair-orbits | 5 × 60k, best **2733** | **no** |
| `srg_anneal` | none | 2 × 800k, best **2546** | **no** |
| `pair_anneal` | pair model, D×D only | 1 × 1.5M, best **2460** | **no** |
| `ga_control_243` | *the positive control itself* | 1 × 25k, best **24729** (target 0) | **the control failed** |
| `orbit_matrix` m=9,11,3 | complete enumeration | no first solution in minutes | n/a |

"Controlled" means the engine had first rebuilt a known graph *at comparable
scale*. None had. So none of these stopping points is evidence about 99; they
are recorded so the next attempt knows what was tried and how far it got.

### Reproducing every number

```bash
python -m pytest tests/test_conway99_orbits.py -q
python problems/conway-99/explore/orbit_calibrate.py
python problems/conway-99/explore/lift_calibrate.py
python problems/conway-99/explore/orbit_matrix.py -m 33 --sorted-diagonal
python problems/conway-99/explore/m33_exhaustive.py
python problems/conway-99/explore/m33_verify.py
python problems/conway-99/explore/order3_lemma.py
python problems/conway-99/explore/make_symmetry_evidence.py
```

## Outcome

**`MAP`**, containing one new elimination that stands on its own and one lemma
that does a lot of work.

`VERIFIED` — resting only on exhaustive computation, double-checked, and on no
lemma that has yet to face an adversary:

- **No SRG(99,14,1,2) admits a *semiregular* automorphism of order 33** (three
  orbits of 33). Scope: exhaustive over the unique orbit matrix and all 4096
  diagonal choices, verified two independent ways.
- **The orbit profile (1×33, 6×11) admits no orbit matrix at all** — standalone
  exhaustive enumeration, 963673 nodes.
- Semiregular Z₃₃ has exactly one orbit matrix (two independent derivations).

`LIVE` (proved here, no skeptic pass): **L5** — an order-3 automorphism has 0 or
3 fixed points — and its consequences *no automorphism of order 27* and *an
order-9 automorphism is semiregular*.

`LIVE`, and explicitly **not** `VERIFIED`: **no automorphism of order 33 at
all**. The full statement needs 001's L2 (to restrict the profiles) and L5 (to
kill (2×33, 3×11), whose 6 orbit matrices were never lifted). Both are
unreviewed proofs, so the combined claim inherits their status. If a skeptic
pass clears L2 and L5 this upgrades to a complete elimination of order 33;
until then it must be quoted with the dependency attached.

**What is not claimed.** No graph was found. Orders 2, 3, 6, 7, 9, 11 are
untouched, as is the trivial-automorphism case, which is where the answer most
likely lives. The eight null construction runs are **not** evidence of
nonexistence. No SAT encoding or UNSAT answer is involved anywhere.

## Why it failed / what survived

### The obstruction

**The orbit-matrix route is complete but does not scale.** It closed Z₃₃
outright and it closed it *cheaply* — 94 search nodes for the enumeration, a
4096-case analysis for the lift. But it terminated only where the number of
orbits was 3. For t ≥ 9 (Z₁₁ with 9 orbits, Z₉ with 11, Z₃ with 33, and the
15-orbit Z₇ profile) the backtracking search does not reach even a *first*
solution in minutes. The exact quantity that fails is the branching: each row
of the orbit matrix has thousands of arrangements consistent with the row-sum
and sum-of-squares constraints, and the coupling equations between rows only
bite once several rows are complete.

The reason Z₃₃ was tractable is worth naming precisely, because it says what a
successor needs: at t = 3 the diagonal equations form a **square** system in the
off-diagonal autocorrelations, so the diagonal *determines* everything else. At
t = 9 there are 9 diagonal equations and 36 unknown autocorrelations — the
system is underdetermined and the trick evaporates.

**Every construction search fails its own positive control.** `cyclic_anneal`
rebuilds the λ=1 graph at n=27, t=9 in 16 seconds and cannot do n=45 at the same
t; `srg_anneal` tops out at n=27; `group_anneal` cannot rebuild BvLS(243) even
with a genuine automorphism of it supplied. The target is n=99. So the low
energies reached under assumed symmetry (336 at Z₁₁, an order of magnitude
below the unrestricted 2546) show that **the symmetry assumption is doing real
work**, and nothing more than that.

### What the method caught

An earlier `order3_lemma.py` check "failed" on rook(4). The lemma was right; the
*test* had hardcoded the λ=1, μ=2 form 2f−2 of the counting identity instead of
the general deg_F(u)(1+λ−μ) + μ(f−1). This is the second time in this cycle that
a bug surfaced only because the same fact was derived twice, once by hand and
once by machine, and the two disagreed — the first being the propagator gap
recorded in 001. Neither would have shown up from testing outputs alone.

### What survived, and is reusable

- The **forced-autocorrelation technique** for t = 3 orbit matrices. It is what
  closed Z₃₃ and it should close any other t = 3 profile in this family.
- **L5**, the most productive lemma of the cycle: it alone kills order 27 and
  the mixed order-33 profiles, and pins order-9 actions to be semiregular.
- **Three calibrated search engines** and the general orbit machinery, including
  `orbit_general.py`, which handles fixed points and mixed orbit sizes and is
  what any Z₇ or Z₂ attack will need.
- The **calibration-ceiling table**, which converts "my search found nothing"
  into a quantified statement of how much stronger a search must be. Gate any
  future search on `lift_calibrate.py` and `ga_control_243.py` first.

## Leads generated

1. **Find the t ≥ 4 successor to the forced-autocorrelation argument.** At t = 3
   the diagonal system is square; at t = 9 it is 9 equations in 36 unknowns, but
   the 81 *off-diagonal* equations are entirely unused. Set up the full system
   (diagonal + off-diagonal, in the group ring of Z₁₁) and determine its rank.
   Definite outcome: either Z₁₁ falls the way Z₃₃ did, or a genuinely free
   parameter is exhibited and the case is known to need search.

2. **Lift the six surviving (2×33, 3×11) orbit matrices.** This is the only gap
   between the verified semiregular result and a clean "no automorphism of order
   33" that does not lean on L5. The lift model for mixed orbit sizes under Z₃₃
   is the missing piece: orbits of size 11 are Z₃₃/⟨11⟩. Small and definite.

3. **Order 7 is the smallest search this problem admits and remains untouched.**
   001's L3 pins it completely — one fixed vertex, 14 orbits of 7, the fixed
   vertex adjacent to exactly two whole orbits, and r₁₁ = r₂₂ = 0, r₁₂ = 1
   forced by the perfect-matching condition. `orbit_general.py` does not
   terminate on the 15-orbit profile; the connection-set formulation (14×14 sets
   over Z₇, energy as in `cyclic_anneal` plus a fixed point) should.

4. **Extend L5 to involutions.** Published work says an even |G| divides 6, so
   Z₂, Z₆, S₃ are the entire even story, and L5 already constrains the order-3
   part of Z₆ and S₃. The gap is precise: even order breaks the "odd order fixes
   a set of size ≤ 2 pointwise" step, so Fix is no longer closed under common
   neighbours and a different argument is needed.

5. **Raise the construction ceiling before running anything larger.** The
   binding constraint is search quality, not compute: every engine tops out
   below n=45. Add targeted move selection (bias moves toward the vertices or
   orbit-pairs carrying the largest local deviation) and re-measure on
   (45,12,3,3) and (50,7,0,1). Outcome either way: the ceiling passes n=50 and a
   run at 99 becomes meaningful, or annealing is the wrong family for λ=1 and
   effort should go to leads 1–3.

## References

- `attempts/001-forced-pair-model-and-partner-regular-split.md` — the pair
  model, the bijection, the partner-regular split, lemmas L1–L4, and the
  positive-control finding, all assumed here.
- `problems/conway-99/PROBLEM.md` — tier-0 statement, machine transcribed `[T]`.
  The only source read in either attempt.
- No papers consulted. Makhnev–Minakova, Behbahani–Lam, Cesarz–Woldar and
  Keramatipour are named in `PROBLEM.md` but were **not** read; L5 is an
  independent re-derivation and I make no claim about its relation to them.
- Data: `problems/conway-99/data/symmetry-evidence.json` (regenerated by
  `make_symmetry_evidence.py`), `orbit-matrices-m33.json`, and the raw run logs
  in `problems/conway-99/data/runs/`.
