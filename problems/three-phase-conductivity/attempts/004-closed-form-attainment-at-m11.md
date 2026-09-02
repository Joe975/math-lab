# 004 — A closed-form rank-5 laminate attaining the HS lower bound at f₁ = m₁₁

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (follows 001–003)
- **Type:** construction (object), exact
- **Tools:** `explore/tp_attain.py` (the construction, `--selftest`, CLI, record
  writer); `explore/tp_rank5.py` (the delegated axis-normal sweep that produced
  the seed); `explore/tp_fields.py` (002's diagnostics). Tests:
  `tests/test_three_phase_attain.py` (6). Record:
  `data/attained/rank5-m11-1-2-5.json`, which passes `verify_laminate.py`.
- **Sources:** attempts 001–003 in this repo; `data/structures-check.md` for
  Cherkaev 2009 Thm 8.1 and §8.1 [T].

## Approach

003 proved attainability is upward closed in f₁ along a fixed f₂:f₃ ray, so the
whole construction problem collapses to finding **one** structure at the bottom
of a ray. That turned a broad search into a point target. The seed came from a
delegated rank-5 sweep whose best structure at f = (1/8, 2/8, 5/8) sat 1.8e−7
above the bound with anisotropy 2.8e−14 — close enough that an exact solve was
worth attempting rather than more searching.

The key move was *not* rounding the floats. Rounding moves the volume fractions
off target, which moves HS_lo too, so the residual never vanishes (tried: it
does not work). Instead the volume-fraction constraints were imposed exactly and
solved symbolically, leaving a small family in which isotropy was solved by
exact bisection.

## What was done

### 1. The exact solve

For the seed topology the constraints reduce cleanly. With Q = (1−a₀) − f₂,

    a₄ = (f₂/(1−a₀) − (1−a₃))/a₃,     a₁ = (f₁ − Q(1−a₅))/a₀,

and the f₁ + f₃ constraint is automatically satisfied (it reduces to
1 − f₂ = f₁ + f₃), so a₃ is free and the family is two-dimensional after
isotropy — matching 002's R−3 count at R = 5.

Exact refinement drove the gap to 1.6e−24 with fractions exactly (1/8, 1/4, 5/8)
and isotropy residual 6e−49, with a₀ converging to 1/2 to twelve digits. Setting
a₀ = 1/2 exactly and solving gave a₅ = 3/4 exactly, and then **the value is
HS_lo = 3 exactly for every a₃ in a range** — a one-parameter family, which is
why the float optimizer wandered instead of converging to a point.

### 2. The construction, in closed form

Generalising by solving at other m₂ and other conductivity triples and
pattern-matching the results (then checking exactly), with

    Θ = σ₁(σ₃−σ₂)/((σ₂+σ₁)(σ₃−σ₁)),   r = √m₂,

set f₂ = r², **f₁ = m₁₁ = 2Θ·r(1−r)**, f₃ = 1 − f₁ − f₂, and build

    A    = laminate(p1 at a₁, p3),   normal e₂
    D    = laminate(p3 at a₅, p1),   normal e₁
    C    = laminate(p2 at a₄, D),    normal e₂
    B    = laminate(C  at a₃, p2),   normal e₂
    root = laminate(A  at a₀, B),    normal e₁

with

    a₀ = 1 − r,   a₁ = r·Θ,   a₃ free in (1−r, 1),
    a₄ = (a₃ − (1−r))/a₃,   a₅ = (σ₂σ₃ − σ₁²)/((σ₂+σ₁)(σ₃−σ₁)).

Then σ\* is isotropic and **equals HS_lo exactly, in ℚ**.

At σ = (1,2,5), r = 1/2: f = (1/8, 1/4, 5/8), HS_lo = 3, parameters
(1/2, 1/8, 4/5, 3/8, 3/4), effective tensor exactly 3·I.

Parameter ranges, argued rather than assumed: a₅ ∈ (0,1) because σ₂σ₃ > σ₁²
and σ₃ > σ₂; a₁ = rΘ ∈ (0,1) because Θ ∈ (0,1); a₄ ∈ (0,1) exactly when
a₃ > 1−r, which is the stated range. So the only genuine restrictions are
a₃ > 1−r and f₃ > 0.

### 3. Why it matters: it beats Milton's threshold

Milton's classical assemblage (001) attains only for f₁ ≥ 2Θ(1−m₂). Here
f₁ = 2Θ·√m₂(1−√m₂), and √m₂(1−√m₂) < 1−m₂ for every m₂ ∈ (0,1), so this
construction reaches strictly lower f₁ at every m₂
(`test_beats_miltons_threshold`). At σ = (1,2,5), m₂ = 1/4: Milton needs
f₁ ≥ 3/8, this attains at f₁ = 1/8.

Combined with **003's coating lemma** — coating an HS-optimal structure with
extra comparison-medium phase preserves HS-optimality, so attainability is
upward closed along the ray — this gives attainment for **all f₁ ≥ m₁₁**. That
is a constructive route to the statement Cherkaev 2009 Thm 8.1 makes for his
structure L13,2,13,1,1 [T], reached here without his formulas: the literature
pass could not reconstruct his T²-structure, so this was solved rather than
transcribed.

### 4. It satisfies 002's predicted field conditions exactly

The strongest internal corroboration, since 002's diagnostics were derived
before this object existed. At σ = (1,2,5), r = 1/2, applied field (1,1):

| phase | field | target | uniform |
|---|---|---|---|
| p2 | (4/3, 4/3) | (4/3, 4/3) | yes |
| p3 | (2/3, 2/3) | (2/3, 2/3) | yes |
| p1 | (10/3, 2/3) and (2/3, 10/3) | mean (2,2) | no, by design |

Phases 2 and 3 are uniform exactly on target; phase 1 carries a mirror-image
pair, giving `gap = 0`, `mean_term = 0`, and `var_term = var_required = 2/9`
exactly — 002's required nonzero fluctuation realised to the digit.

### 5. Verification so far

- Exact attainment on **1055 random (σ, r, a₃) instances, zero failures**
  (1099 skipped for a₃ ≤ 1−r or f₃ ≤ 0, the stated domain).
- Both harness routes agree (`laminate.py` projection form vs
  `verify_laminate.py` rotated-frame interface averaging) on every instance.
- Keller–Dykhne holds on every tree.
- The written record passes `verify_laminate.py` from the command line.
- 6 tests in `tests/test_three_phase_attain.py`.

## Outcome

`LIVE`, pending the skeptic pass — **not** `VERIFIED`. The repo's bar is an
independent re-implementation by a different route, and the two agreeing routes
here are both in the harness this attempt used. A skeptic pass (third
implementation from scratch, adversarial hunt for failing instances, derivation
of the two pattern-matched identities a₁ = rΘ and a₅ = (σ₂σ₃−σ₁²)/((σ₂+σ₁)(σ₃−σ₁)),
and a novelty check against the published structures) was commissioned and is
recorded separately.

**Scope:** the construction is exact at f₁ = m₁₁ exactly, for the (σ, r, a₃)
domain stated. Everything above it on the ray follows only via 003's coating
lemma, which is proved but is a separate statement.

**Not claimed.** Nothing about f₁ < m₁₁ — whether m₁₁ is the true attainability
threshold is untouched here. No claim of novelty: this reproduces what Cherkaev
2009 Thm 8.1 asserts, and should be recorded as an independent constructive
rediscovery unless the skeptic's literature check says otherwise. The two
parameter identities were found by pattern-matching exact numerics and checked
widely, not derived.

## Why it failed / what survived

Nothing failed. Two process notes worth keeping. First, rounding a float
optimum to nearby rationals does **not** work here and cost a cycle: the volume
constraints must be imposed exactly first, because moving the fractions moves
the target. Second, the reason the float search never converged cleanly is that
the attaining set is a **one-parameter family** (a₃ free), so the optimizer had
a flat direction — a signature worth recognising in future searches, since a
wandering optimizer at a suspiciously small residual suggests an exact family
rather than a hard optimum.

This also resolves 003's open tension in favour of attainment: the near-constant
per-rank gap ratio (2.1e−2, 5.46e−4, 1.6e−5) was not geometric convergence to an
unattainable bound; rank 5 attains exactly, and the sweep's residuals were the
optimizer failing to land on a flat family.

Reusable: `tp_attain.py` as a generator of HS-optimal structures for arbitrary
phases; the exact-constraint-then-bisect solve pattern; the observation that a
flat optimizer direction signals an exact family.

## Leads generated

1. **Is m₁₁ the true threshold?** Search below it, at f₁ slightly under
   2Θ√m₂(1−√m₂), for any attaining structure at rank 5–7. Combined with 003's
   result that no phase-1-free structure attains, a negative would bracket the
   true threshold and a positive would beat the published one. This is now the
   frontier question for the problem.
2. **Derive a₁ = rΘ and a₅ = (σ₂σ₃−σ₁²)/((σ₂+σ₁)(σ₃−σ₁))** from the 002 field
   conditions rather than pattern-matching: impose E₂ and E₃ at target via the
   path-product rule and solve. Would turn the construction into a proof.
3. **Dualise it.** By 002 §6 the Keller–Dykhne image of this structure attains
   the HS *upper* bound at transformed fractions; write that down and verify, for
   free.
4. **Minimal rank.** Is rank 5 minimal for attainment at m₁₁? 002's parameter
   count suggested more would be needed; the flat family says otherwise. Settle
   by exhaustive rank-4 search at exactly f₁ = m₁₁, which is cheap.
5. **Compare with Cherkaev's parameters** once the appendix is worked through:
   do our (a₀…a₅) match his L13,2,13,1,1, or is this a genuinely different
   structure attaining the same bound?

## References

- A. Cherkaev, Mech. Mater. 41 (2009) 411–433; preprint arXiv:1009.3060,
  Thm 8.1 and §8.1 [T] — via `data/structures-check.md`.
- G. W. Milton, Appl. Phys. A 26 (1981) 125–130 — not reached; threshold via 001.
- This repo: `attempts/001`, `attempts/002`, `attempts/003`.
