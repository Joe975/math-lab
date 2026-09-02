# 002 — What attaining the HS bound demands of the fields: an exact mechanism, two dead conjectures, and a parameter count

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (follows 001; reads its record, its certified data, and
  `data/literature-check.md`)
- **Type:** mechanism analysis + exact computation
- **Tools:** `explore/tp_fields.py` (exact field propagation through a laminate
  tree; HS attainment fields; required phase-1 variance; exact gap
  decomposition; `--selftest`, `--report`), all standard-library Python.
  Tests: `tests/test_three_phase_fields.py` (6 tests). Data re-used:
  `data/certified/` (27 records from 001).
- **Sources:** attempt 001 in this repo; `data/literature-check.md` [T] for
  Milton's threshold as quoted by Albin–Cherkaev–Nesi 2007 and Cherkaev 2009.
  No new literature was consulted for the mathematics below.

## Approach

001 left the central question open: in the region where the classical coated
assemblage is infeasible, is the HS lower bound attained at all? 001 attacked
it by pushing lamination rank and measuring the energy gap. That is a blunt
instrument — the energy gap is one number, and it cannot distinguish
"converging to the bound" from "stalling just above it".

This attempt asks instead what attainment *demands*, and measures those
demands directly. The lever: attaining a variational bound has an equality
case, and the equality case is a statement about the fields. If the demands
are measurable on any candidate structure, then every structure becomes a
multi-channel observation rather than a single number, and a stall becomes
diagnosable.

Why this rather than more search: 001's own data already contained the
answer to several natural conjectures, unread. Extracting it costs no new
search, and it re-scopes what a rank-5 search can even establish (see §5).

## What was done

### 1. The attainment fields, and Milton's threshold restated

With comparison medium σ₁ (the smallest phase), attaining the multiphase HS
lower bound forces the field to be **uniform inside every phase of
nonvanishing polarization**, i.e. every i ≠ 1, at the value

    E_i = (σ_HS + σ₁)/(σ_i + σ₁) · E₀.                                   (*)

Phase 1 *is* the comparison medium, so its polarization (σ−σ₁)E vanishes
identically; its field is free pointwise and only its volume average is
pinned. `SPECULATION` (load-bearing, and standard in the homogenization
literature but not re-derived here): that (*) is *necessary* for attainment
is the classical equality analysis of the HS variational principle; this
attempt verifies it holds for the attaining structure but does not prove
necessity. Everything downstream that uses (*) as a *diagnostic* is
unaffected; only the claim "a structure violating (*) cannot attain" leans on
it.

Verified exactly (`tp_fields.py --selftest`, and
`tests/test_three_phase_fields.py`): the Milton assemblage of 001 reproduces
(*) to the last digit in phases 2 and 3 at every feasible point tried, its
phase-1 field is genuinely non-uniform, and its phase-1 volume average is
exactly the target. Summing (*) against the volume fractions recovers the HS
identity, so (*) is self-consistent.

Consequence, and the reason (*) is worth stating: since σ_HS > σ₁ always,
(*) gives E₁ > E₀ > E₃, and

    E_2 ≤ E₀   ⟺   σ_HS ≤ σ₂   ⟺   m₁ ≥ 2Θ(1−m₂),

the last being Milton's published threshold [T]. So **Milton's attainability
condition is a pointwise field constraint**: the middle phase must not be
required to carry more than the applied field. Checked on the whole k/8 grid
at σ = (1,2,5) (`test_middle_phase_field_condition_is_miltons_threshold`);
the three formulations agree at every point.

### 2. Two conjectures killed by 001's own data

The reformulation suggests an obvious mechanism — laminates fail in the gap
region because they cannot push the middle phase above the applied field.
**False.** Running the field tool over the certified records:

| record | phase-2 field ratio \|E₂\|/\|E₀\| | target |
|---|---|---|
| r3_1-8_2-8_5-8_lower | 1.5126 | 1.3333 |
| r3_1-8_1-8_6-8_lower | 1.3900 | 1.4545 |
| r4_1-1-6_lower | 1.4746 | 1.4545 |

Laminates exceed E₀ in phase 2 routinely, up to 1.51× measured, and even
overshoot the target. So exceeding E₀ is not the obstruction.

Second conjecture, also false. Attainment additionally requires a determined
amount of phase-1 fluctuation (§3); one might expect the gap to be driven by
a *deficit* of it. But the achieved variance **exceeds** the required value in
records that still have a strictly positive gap (ratios 1.0045 and 1.0139 at
gaps 1.9e−3 and 4.7e−4), and the gap is not monotone in the ratio (ratio
0.6804 at gap 8.1e−3, ratio 0.9878 at gap 5.5e−3). So the gap is not a
one-sided variance deficit either.

### 3. Attainment requires a *precise nonzero* fluctuation in the weakest phase

This is the new positive result. Suppose (*) holds. Then energy bookkeeping
gives σ* = Σ_i f_i σ_i |E_i^tgt|²/|E₀|² + f₁σ₁·Var₁, and setting σ* = HS_lo
and simplifying (write σ_i = (σ_i+σ₁) − σ₁ in the numerator) yields the
closed form

    Var₁ = (c² S₂ − 1)/f₁,   c = 1/Σ_i f_i/(σ_i+σ₁),   S₂ = Σ_i f_i/(σ_i+σ₁)²,

per |E₀|². By Cauchy–Schwarz, (Σ f_i/(σ_i+σ₁))² ≤ (Σ f_i)(Σ f_i/(σ_i+σ₁)²), so
c²S₂ ≥ 1 with equality only when all phases coincide. **So attaining the
bound does not merely permit field fluctuation in the comparison-medium
phase — it requires a strictly positive, exactly determined amount of it.**
The closed form is checked against independent energy bookkeeping on three
instances including an unrelated conductivity triple
(`test_required_phase1_variance_closed_form`), and the Milton assemblage
realises it exactly (e.g. Var₁ = 1/2 at f = (3,2,3)/8, σ = (1,2,5), matching
the closed form's 1/2).

### 4. The gap decomposes exactly

Combining, for any tree, with E₀ arbitrary and the effective value taken as
the quadratic form E₀·σ*E₀/|E₀|² so no exact isotropy is needed:

    gap = Σ_i f_i σ_i (|mean_i|² − |E_i^tgt|²)/|E₀|²
        + Σ_i f_i σ_i Var_i/|E₀|²
        − f₁σ₁·Var_req.

`gap_decomposition()` computes both sides independently and asserts equality;
it reproduces in ℚ on all 12 certified lower-side records at σ = (1,2,5) and
on the assemblage, where the mean term is 0 and the variance term equals
Var_req exactly. The mean term can be negative (−0.00293 at
r3_3-8_1-8_4-8_lower), so only the total is signed. This converts one number
into two independent defect channels, which is what §5 needs.

### 5. A parameter count that re-scopes the search

A rank-R tree carries R node fractions. Three are always spent: two volume
fractions and isotropy. So R−3 free parameters remain — **zero at rank 3, one
at rank 4, two at rank 5**. In every certified optimum so far phase 1 occupies
exactly two leaves and phases 2 and 3 one leaf each, so their fields are
automatically uniform and (*) reduces to hitting two target *values*: 2 scalar
conditions each, plus the variance condition, up to 5 further conditions.

`SPECULATION` (labelled, and the count is heuristic — the conditions are not
shown independent): exact attainment needs roughly rank 7–8, so **rank 5 is
parameter-starved and a stall there is not evidence of a floor.** This
matters: it is exactly the misreading 001's lead 1 invited. The replacement
experiment is to *solve* the attainment constraint system rather than minimise
energy, and report the smallest rank at which it is solvable.

Supporting tool, verified: for axis-normal trees every tensor stays diagonal
and the leaf fields are path products,

    E_leaf,1 = E₀,₁ · Π over e₁-laminations on the path of (T_parent,11/T_child,11)
    E_leaf,2 = E₀,₂ · Π over e₂-laminations on the path of (T_parent,22/T_child,22)

(tangential continuity passes one component through untouched; normal-flux
continuity fixes the other). Checked against the general propagator on three
trees at three applied-field directions each, exact agreement. It turns the
target conditions into two multiplicative equations per phase.

### 6. The upper side is not independent

Keller–Dykhne makes the upper-bound problem the lower-bound problem at
inverted conductivities: HS_hi(f, σ) = 1/HS_lo(f, 1/σ) with the (f_i, σ_i)
pairing preserved, and an isotropic tree of value s dualises to one of value
1/s. Verified exactly on all 27 certified records. Concretely, the upper-side
problem at σ = (1,2,5), f = (1,1,6)/8 **is** the lower-side problem at
σ = (1,5/2,5), f = (3/4,1/8,1/8), after rescaling by 5 and relabelling. So
001's upper-side screens were not wasted but were also not independent
evidence, and future censuses should run the lower side only.

## Outcome

`MAP`, with `VERIFIED` sub-claims scoped as stated.

- `VERIFIED` (exact in ℚ; the instances and grids named above): the
  assemblage realises (*); the equivalence of E₂ ≤ E₀, HS_lo ≤ σ₂ and
  Milton's threshold on the k/8 grid at σ = (1,2,5); the closed form for
  Var_req on three instances plus its positivity by Cauchy–Schwarz; the gap
  decomposition on 12 certified records; the duality reduction on 27 records;
  the path-product rule on the trees tried.
- `REFUTED`: "laminates cannot drive the middle phase above the applied
  field" (measured to 1.51×); "the gap is a one-sided phase-1 variance
  deficit" (achieved exceeds required at positive gap, and the gap is not
  monotone in the ratio).
- `SPECULATION`, labelled inline: the necessity of (*) (classical, not
  re-derived here); the rank 7–8 parameter count.

**Not claimed.** Nothing about whether HS is attainable in the gap region —
this attempt sharpens the question and supplies instruments, it does not
answer it. The necessity of (*) is assumed, not proved, so no structure is
ruled out here. The parameter count is heuristic. Nothing about all
microstructures; laminates only.

## Why it failed / what survived

The hoped-for obstruction did not materialise: both natural one-line
mechanisms are false against 001's own data, and that is the honest headline.
What survived is better instrumentation — attainment now has three measurable
channels (phase-2/3 uniformity, distance to target, phase-1 variance vs
Var_req) instead of one energy number, and an exact identity tying them to
the gap. The parameter count explains *why* low ranks fail without invoking
any floor, which is the most likely reason 001's rank sequence looked the way
it did.

Reusable: `tp_fields.py` in full; the closed form for Var_req; the gap
decomposition as a checker (it asserts its own identity, so it catches tree
and field bugs); the path-product rule as a cheap exact evaluator for
axis-normal trees; the duality reduction as a 2× economy on every future
census.

## Leads generated

1. **Solve the attainment system instead of minimising energy.** On a fixed
   topology, solve {volume fractions, isotropy, E₂ = target, E₃ = target}
   exactly and report the smallest rank at which a solution exists; then check
   whether σ* = HS_lo there. Definite outcome; directly tests §5's count.
2. **Prove or refute the necessity of (*)** for this problem, closing the one
   `SPECULATION` the diagnostics rest on. A textbook derivation of the
   multiphase HS equality case in 2D would settle it; the skeptic pass should
   demand a citation or a proof.
3. **Is Var₁ = Var_req automatic** given (*), correct fractions and isotropy?
   If yes the count in §5 drops by one condition and the predicted minimal
   rank falls. Falsifiable: exhibit a structure satisfying (*) with
   Var₁ ≠ Var_req, or prove the implication.
4. **Where does Var_req blow up?** It is (c²S₂−1)/f₁, so it diverges as
   f₁ → 0. Compare its growth with the feasibility floor f₁* ≈ 0.0301 found
   for 001's rank-3 shape: if the shape dies exactly where the required
   fluctuation outruns what two leaves can supply, that is the mechanism 001
   was missing. Compute both curves and check.
5. **Extend the diagnostics to the upper side for free** via §6, and check
   that the dual of a certified lower-side record reproduces the upper-side
   defects — a cheap consistency test of the whole framework.

## References

- `problems/three-phase-conductivity/attempts/001-coated-laminate-census.md`
- `problems/three-phase-conductivity/data/literature-check.md` (Milton's
  threshold as quoted by ACN 2007 and Cherkaev 2009, both [T])
- Harness: `harness/three-phase-conductivity/laminate.py`, `verify_laminate.py`
