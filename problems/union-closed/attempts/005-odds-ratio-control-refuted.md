# 005 — Gap 1 attacked: Plackett odds-ratio control is false as stated, sharply

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-30
- **Mode:** informed
- **Type:** proof-gap attack (refutation with proofs + partial positive results),
  on the labeled gap `plackett-odds-ratio-control` (Gap 1) of the LIVE
  dependent-couplings route. Follows 003 §3 item 1 and 003 lead 1, using ONLY
  004's corrected statements where 003 and 004 disagree (recipe ceiling
  0.431496 in closed form; restated mini-theorem; half-mixing coupling).
- **Tools:** `explore/uc_odds_ratio.py` (written and run here; standard library
  only; deterministic, fixed seeds; runtime ~6 s; checkpoints
  `data/005_part[A-G].json`, log `data/005_full_run.log`). Command reproducing
  every number below:
  `python problems/union-closed/explore/uc_odds_ratio.py`
- **Sources:** 003/004 (this repo). Karlin–Rinott, *Classes of orderings of
  measures and related correlation inequalities I: Multivariate totally
  positive distributions*, J. Multivariate Anal. 10 (1980) — used for the
  standard facts that MTP₂ is preserved under restriction and marginalization;
  cited from general knowledge of the standard statement, not re-checked
  against the paper this cycle [T]. No new external fetches.

Notation as in 003: tilt coupling `π_λ(A,B) ∝ u(A)u(B)·2^{λ|A∩B|}` with both
marginals μ (Sinkhorn potentials; existence and uniqueness of the coupling on
`supp(μ)²` are classical for the strictly positive kernel, and symmetry of
kernel + marginals + uniqueness lets one take equal potentials). Coordinates
revealed in the fixed order 1..n; histories `a = A_{<i}`, `b = B_{<i}`;
`OR_i(a,b)` is the odds ratio `P₁₁P₀₀/(P₁₀P₀₁)` of the conditional table of
`(A_i, B_i)` given `(a,b)` (defined at histories where all four cells have
positive mass). `m = n − i` is the number of future coordinates.
`ψ = (3−√5)/2`, `q = 2p−p²`.

## Approach

**The target.** Gap 1 as recorded in 003 (SPECULATION there, label confirmed
accurate by 004): *for λ ≥ 0, the Sinkhorn tilt on general μ has all
conditional odds ratios in `[1, 2^λ]`* — the bit-level lemma the
tax-controlled assembly (003 lead 1) needs, "FKG-flavored, plausibly
provable". The instruction for this cycle: prove it, prove a clean partial
statement, or exhibit a precise obstruction.

**Why structure theory rather than the obvious alternative.** The obvious
moves are the FKG/four-function induction 003 sketched, or a blind numerical
sweep hunting a violation. Instead I first reduced `OR_i(a,b)` to a closed
algebraic form in the potential's future slices, because the reduction
decides the question *in both directions at once*: it hands you (i) the exact
inequality that would have to hold, (ii) the mechanism that breaks it, and
(iii) the sharp constants. The reduction is three lines; everything else
falls out of it. The key fact making counterexamples constructible: **any**
positive function u on a support S is the Sinkhorn potential of its own
marginal (define π from u, read off μ; uniqueness of the scaling makes π THE
tilt coupling for that μ) — so the four future slices entering a given
`OR_i(a,b)` with `a ≠ b` are four *unconstrained* nonnegative vectors, and
the lemma, if true, would have to be a fact about the kernel alone. It is
not.

**The reduction** (rigorous; used throughout). Split `A = (a, α, x)`,
`B = (b, β, y)` with `x, y ∈ {0,1}^m`. Since
`|A∩B| = |a∧b| + αβ + ⟨x,y⟩`, the conditional table is

    P(α,β | a,b) ∝ 2^{λαβ} · ⟨F_α, W G_β⟩,
    F_α(x) = u(a,α,x),  G_β(y) = u(b,β,y),  W(x,y) = 2^{λ⟨x,y⟩},

so, with the prefix factor cancelling,

    OR_i(a,b) = 2^λ · R,   R = ⟨F₁,WG₁⟩⟨F₀,WG₀⟩ / (⟨F₁,WG₀⟩⟨F₀,WG₁⟩).

The claimed control `OR ∈ [1, 2^λ]` is exactly `R ∈ [2^{−λ}, 1]`.

## What was done

### Proposition 1 (last coordinate; VERIFIED, trivial)
At `i = n` (no future, `m = 0`): F, G are scalars and cancel, so
`OR_n(a,b) = 2^λ` exactly at every nondegenerate history, for every μ.
(Engine check: part A reproduces the known product-μ identity `OR ≡ 2^λ` at
all 85 histories of a `Bern(0.3)^4` fit, to 1.4e-15.)

### Proposition 2 (diagonal histories; VERIFIED — kills the upper half)
**For λ ≥ 0, every μ, every i, every nondegenerate diagonal history
`a = b`: `OR_i(a,a) ≥ 2^λ`, with equality iff the two future slices
`u(a,1,·)` and `u(a,0,·)` are proportional.**

*Proof.* `W = ⊗_{j>i} [[1,1],[1,2^λ]]` is positive semidefinite (each factor
has det `2^λ − 1 ≥ 0`), strictly definite for λ > 0. On the diagonal
`G_β = F_β`, so `R = ⟨F₁,WF₁⟩⟨F₀,WF₀⟩ / ⟨F₁,WF₀⟩² ≥ 1` is exactly
Cauchy–Schwarz in the W-inner product; equality iff `F₁ ∥ F₀`. ∎

Consequences. Product μ has proportional slices — equality everywhere,
recovering the identity 003/004 verified. Any other μ has `OR > 2^λ`
strictly wherever its potential slices fail to be proportional; diagonal
histories always have positive mass (positive kernel). So the claimed upper
bound `OR ≤ 2^λ` fails *in the opposite direction*: on the diagonal the tilt
over-diagonalizes, it never under-tilts. Computation (part B): 60 random
full-support potentials at n = 5, λ ∈ {0.5, 1, 2} — minimum of `OR/2^λ` over
all diagonal histories = 1.0000000011 (≥ 1 as proved; the ~1e-9 case is one
chance near-proportional slice pair), every trial had ALL its `i < n`
diagonal histories strictly above (median per-trial minimum excess 5.4e-4,
max observed ratio 3.84), and the `i = n` rows sit at exactly `2^λ` per
Prop 1.

### Proposition 3 (the crash family; VERIFIED — kills the lower half)
**For every n ≥ 4 and λ > 0 there is an explicit μ with all element
marginals ≤ 0.37 < 0.38271 and H(μ) ≈ 1.76 bits whose tilt coupling has a
positive-mass history with `OR = 2^{λ(3−n)}` — arbitrarily far below 1.**

Take μ on four atoms, `S₁ = {2}`, `S₂ = {3..n}`, `S₃ = [n]`, `S₄ = {1}`,
with masses (0.33, 0.33, 0.04, 0.30). At `i = 2`, history
`(a,b) = (A₁=0, B₁=1)`: the four conditional cells each contain exactly one
atom pair (coordinates 1, 2 classify the four atoms), so all four potential
values appear once in the numerator and once in the denominator of OR and
**cancel identically**. What remains is kernel-only:

    OR = 2^{λ(|S₁∩S₃| + |S₂∩S₄| − |S₁∩S₄| − |S₂∩S₃|)} = 2^{λ(1+0−0−(n−2))}
       = 2^{λ(3−n)}.

This is exact for whatever potentials Sinkhorn produces — no numerics in the
proof at all. Marginals: element 1: 0.34; element 2: 0.37; elements ≥ 3:
0.37 — all below the Liu record 0.38271, so μ is squarely inside the
hypothesis class of (S-coup) at every threshold of interest, and inside the
stated scope of the Gap-1 lemma ("general μ"). Computation (part C):
Sinkhorn fit from μ (residuals < 1e-14) reproduces `log₂ OR = 3 − n` to
9e-16 at n = 5, 8, 11, 13; the crash history carries mass 0.13–0.16.
**Full-support robustness** (part C2): mixing with `ε·Unif(2^[n])` at n = 8
gives `log₂ OR = −4.62, −4.96, −4.996` at ε = 1e-2, 1e-3, 1e-4 (4-atom value
−5), so the crash is not a support-degeneracy artifact — EVIDENCE at these
ε; the exact statement is the 4-atom one.

Where the crash bites: at a crashed coordinate with conditional zero-margins
`x = y = 1−p` and `p ≥ ψ`, OR < 1 pushes the union-bit zero-probability
below the iid value `xy`, i.e. the union marginal above
`q = 2p−p² ≥ 1−p ≥ 1/2`, so `h(m) < h(p)`: the coordinate contributes
**strictly negative** gain, precisely in the record-relevant regime
`p ∈ [ψ, 0.383]`.

### Proposition 4 (sharp universal bounds; VERIFIED)
**For λ ≥ 0, all μ, all nondegenerate histories:
`OR_i(a,b) ∈ [2^{λ(1−m)}, 2^{λ(1+m)}]`, and both ends are attained** (the
lower end by Prop 3's family, where `m = n−2`; the upper end by its mirror,
part D).

*Proof.* R is, in each of the four slice vectors separately, a ratio of two
linear forms with positive coefficients; over the nonnegative orthant
`max_t (Σaᵢtᵢ)/(Σbᵢtᵢ) = maxᵢ aᵢ/bᵢ`, so the extrema of R are attained with
all four slices point masses `δ_{x₁}, δ_{x₀}, δ_{y₁}, δ_{y₀}`, where

    log₂ R = λ·(⟨x₁,y₁⟩+⟨x₀,y₀⟩−⟨x₁,y₀⟩−⟨x₀,y₁⟩) = λ·⟨x₁−x₀, y₁−y₀⟩
           ∈ λ·[−m, m]. ∎

So the deviation of `log₂ OR` from λ is λ times a **future-overlap
discrepancy** `⟨x₁−x₀, y₁−y₀⟩` (exactly, in the extremal point-mass case; a
mixed second-order form in general), and the only universally valid
two-sided control is exponentially wide in the number of remaining
coordinates — useless for the assembly. Computation (part D): 40 random
potentials, no bound violation (minimum slack ~0 attained at the `i = n`
boundary, as it must be); top end `2^{λ(1+m)}` attained exactly by the
mirror 4-atom family.

### Proposition 5 (MTP₂ partial positive result; VERIFIED, one-sided only)
**If the potential u is strictly positive and log-supermodular (MTP₂) on
{0,1}^n, then for all λ ≥ 0, ALL histories — diagonal and off-diagonal —
satisfy `OR_i(a,b) ≥ 2^λ`.**

*Proof.* `Ψ(α,x,β,y) = F_α(x)G_β(y)W(x,y)` is MTP₂ on `{0,1}^{2m+2}`:
log-supermodularity is additive, slices of an MTP₂ function are MTP₂ (F and
G), and W is a product of pairwise log-supermodular factors. Marginalizing
out (x, y) preserves MTP₂ [Karlin–Rinott 1980, standard], so
`Φ(α,β) = ⟨F_α, WG_β⟩` is TP₂, i.e. `R ≥ 1`, i.e. `OR ≥ 2^λ`. ∎

Computation (part F): 15 random ferromagnetic (log-supermodular) potentials
at n = 4, λ = 1: minimum of `OR/2^λ` over every history = 1.000000000
(attained at `i = n`), max diagonal ratio 1.055. So even under the strongest
natural structural hypothesis the control is **one-sided and on the wrong
side** of the claimed interval: `[2^λ, ∞)` meets the claimed `[1, 2^λ]` in
the single point `{2^λ}`.

Probe (part F, EVIDENCE): μ MTP₂ does **not** imply the fitted potential is
MTP₂ — 1 of 24 Sinkhorn fits of random ferromagnetic μ (n = 4) produced a u
violating log-supermodularity, worst pair ratio 0.9922 at λ = 1.5 with
Sinkhorn residual 2.1e-15 (violating instance checkpointed in
`005_partF.json`, key `worst_case_instance`). Yet all 24 fits still had
every `OR ≥ 2^λ` (min ratio 1.000000000). So "Sinkhorn preserves MTP₂" is
false at the potential level, while "μ MTP₂ ⟹ OR ≥ 2^λ" survives all tests
here and is open — it is SPECULATION if asserted; a candidate replacement
hypothesis (lead 2).

### Proposition 6 (perturbative regime; proved modulo one standard step)
**Along any smooth family μ_δ with μ₀ product (fixed full support), every
conditional log-odds-ratio satisfies `log₂ OR_i(a,b) = λ + O(δ²)` — both
sides.** *Proof sketch.* Write `u_δ = u₀(1 + δv + O(δ²))` with u₀ product —
the differentiability of the Sinkhorn potential in μ (fixed finite support,
strictly positive kernel) is standard implicit-function-theorem material and
is the one step taken as known rather than re-proved here [SPECULATION only
in that narrow sense]. Expand each of the four `log⟨F, WG⟩` terms; in the
alternating sum defining `log R`, every term additive in (A-slice) +
(B-slice) contributions cancels, and the first surviving term is the mixed
second-order bilinear form evaluated at the differences `(F₁−F₀, G₁−G₀)`,
which is `O(δ²)`. ∎

Computation (part G): mixture `(1−δ)Bern(0.30)^5 + δBern(0.05)^5`, λ = 1:
max over all histories of `|log₂ OR − λ|` = 1.5e-5, 6.0e-5, 2.3e-4, 8.6e-4,
3.0e-3, 9.4e-3 at δ = 0.005, 0.01, 0.02, 0.04, 0.08, 0.16; successive ratios
3.93, 3.86, 3.73, 3.49, 3.10 → 4, i.e. clean δ² scaling. EVIDENCE for this
family; Prop 6 gives the mechanism.

### Behavior on a realistic instance (part E; EVIDENCE)
Full OR census for the Sawin-genre mixture `0.7·Bern(0.30)^6 +
0.3·Bern(0.05)^6` at λ ∈ {0.5, 1.5}: every history has `OR ≥ 1` (no crash on
benign μ — the crash needs structured adversaries; 60 random dense
potentials in part B likewise never produced OR < 1), about 83% of history
mass has `OR > 2^λ` strictly, a small mass dips below `2^λ` by ≤ 6e-4 in
log₂, and the mass-weighted mean of `log₂ OR` per coordinate is
`λ + (0.000 … 0.043)`, decreasing to exactly λ at `i = n`. On instances the
route actually cares about, essentially the *entire* deviation budget sits
above `2^λ` — the direction the assembly's tax accounting does not cover.

## Outcome

**REFUTED — the Gap-1 lemma (`plackett-odds-ratio-control`) is false as
stated, in both directions, with sharp quantitative versions.** Scope of the
refutation: the pointwise claim "λ ≥ 0 forces all conditional odds ratios of
the Sinkhorn tilt into `[1, 2^λ]` for general μ", exactly as recorded in 003
§3 item 1 / Gap 1. (i) The upper half fails for essentially every
non-product μ: on diagonal histories `OR ≥ 2^λ` always (Prop 2,
Cauchy–Schwarz), with generic strict excess. (ii) The lower half fails as
badly as it possibly can: an explicit in-class μ (marginals ≤ 0.37, positive
entropy) achieves `OR = 2^{λ(3−n)} → 0` at a history of mass ≈ 0.14
(Prop 3), matching the sharp universal range `[2^{λ(1−m)}, 2^{λ(1+m)}]`
(Prop 4). A precise obstruction was a fully successful outcome per the
protocol; this is it.

What survives as theorems: Props 1, 2, 4, 5 (and Prop 6 modulo one standard
smoothness step). The clean partial statement with explicit hypotheses is
Prop 5: potential-level MTP₂ buys the one-sided bound `OR ≥ 2^λ` at every
history — never the claimed upper bound.

**Status of the route: the dependent-couplings interface stays LIVE.**
Nothing here touches the licensing lemma, the adversarial separations, the
no-go evasion (003, as corrected by 004), or the corrected 0.431496 recipe
ceiling. What changes is the shape of Gap 1: pointwise two-sided odds-ratio
control is dead and cannot be the bridge from the functional to a theorem;
any assembly must use history-averaged control (part E suggests the
mass-weighted mean of `log₂ OR` as the right object), or a subclass of μ
with MTP₂-type structure, or the perturbative regime (Prop 6 — which is
*good* news for 003's lead 1: near product μ the deviation is O(δ²), so the
second-order assembly is not blocked by this refutation).

**Not claimed:** no bound, no constant, no progress on Gap 2 (tax) or Gap 3
(totality); no claim that the assembly is impossible — only that this
particular bridge is; no claim about μ-level MTP₂ hypotheses (open; part F
is evidence only); Prop 6's smoothness step was not re-proved here; part C2
full-support numbers are EVIDENCE at the ε tested, not a limit theorem (the
4-atom statement is the theorem); none of this has had its independent
skeptic pass yet — per repo rules the propositions are results only after
one.

## Why it failed / what survived

The step that breaks is 003's bit-level replacement: "`z = z_ρ(x̃,ỹ)` with the
conditional odds ratio ρ drifting inside `[1, 2^λ]`". The quantity that goes
the wrong way is the **future-overlap discrepancy** `⟨x₁−x₀, y₁−y₀⟩`
(Prop 4): the conditional odds ratio at step i is `2^λ` times a cross-ratio
of future-slice interactions that the hoped-for FKG-style induction was
supposed to pin into `[2^{−λ}, 1]`, but which is actually:

- ≥ 1 on diagonal histories — Cauchy–Schwarz *reverses* the needed
  four-function inequality there (the induction is trying to prove something
  already false one step above the base case); and
- unbounded below by any constant off the diagonal — adversarial μ make the
  tilt anti-correlate a coordinate at rate `2^{−λ(m−1)}`, producing strictly
  negative per-coordinate gain for `p ≥ ψ`.

Directional summary for the assembly: realistic μ over-tilt (deviation mass
above `2^λ`, part E) — so the dependence the tax must pay for is
under-accounted by λ; worst-case μ under-tilt catastrophically (Prop 3) — so
the per-coordinate gain itself fails. One parameter λ cannot certify both
ends: the bit-level echo of the two-scale obstruction 004 found at the
recipe level (`δ_∅ ⊕ Bern(1/2+ε)`, ceiling 0.431496).

Survived / reusable:

- The slice reduction `OR = 2^λ·R` with R a cross-ratio of W-inner products
  — three lines, decides every question about conditional ORs of tilt
  couplings; plus the "any positive u is its own marginal's potential" trick
  for building counterexamples to Sinkhorn-tilt claims without running
  Sinkhorn.
- Props 1, 2, 4, 5 as permanent structure theory for the tilt class: exact
  last-coordinate identity; diagonal Cauchy–Schwarz lower bound with
  equality characterization; sharp range via the future-overlap discrepancy;
  MTP₂ ⟹ one-sided control `≥ 2^λ` (Karlin–Rinott route).
- The 4-atom crash family (Prop 3) as the hard-instance genre for ANY future
  claimed pointwise control on tilt couplings — potential-free, exact,
  tunable marginals; kin to 004's `δ_∅ ⊕ Bern(1/2+ε)` killers (both are
  "two blocks with opposite needs" constructions).
- The finding that Sinkhorn scaling does not preserve MTP₂ of the potential
  (explicit violating instance, `data/005_partF.json`) — a warning against
  hypothesizing structure on μ and assuming it transfers to u.
- `explore/uc_odds_ratio.py`: OR-census engine for arbitrary supports +
  Sinkhorn fitter, validated against the product identity (1.4e-15), the
  exact crash values (9e-16), and the closed-form range endpoints.

## Leads generated

1. **(replacement for Gap 1) History-averaged one-sided control.** Part E
   shows the mass-weighted mean `M_i = E_{(a,b)~π}[log₂ OR_i(a,b)]` sitting
   in `[λ, λ+0.05]` on realistic instances. Conjecture (SPECULATION):
   `M_i ≥ λ` for all μ and all i (true restricted to diagonal mass by
   Prop 2; true at `i = n` by Prop 1). Falsifiable first step: compute `M_i`
   on the 4-atom crash family and its mirror — if some `M_i < λ` there, the
   averaged lower bound dies too and the assembly needs the MTP₂-subclass
   route; if it holds, attempt a proof by summing Prop 2 over diagonal mass
   and bounding the off-diagonal contribution via Prop 4.
2. **μ-level MTP₂ hypothesis.** Part F: 24/24 fits of MTP₂ μ had all
   `OR ≥ 2^λ` even when u lost MTP₂. Decide: does μ MTP₂ imply `OR ≥ 2^λ`
   for the fitted tilt? A counterexample search over random MTP₂ μ at
   n = 4–5 and many λ is cheap (extend part F); a proof would need a
   Sinkhorn-compatible correlation inequality bypassing the potential (the
   failure of potential-level transfer says the direct route is closed).
3. **Perturbative assembly, now unblocked at second order.** Prop 6 says ORs
   are `λ + O(δ²)` near product μ. Redo 003 lead 1's O(ε²) expansion at the
   AHS equality point with an OR-deviation budget `c·δ²` from Prop 6's
   bilinear form, and check whether the sign of the net second-order
   coefficient survives a worst-case c. Definite outcome either way.
4. **Two-scale tilts vs the crash family.** The crash needs a
   size-inhomogeneous overlap (`|S₂∩S₃| ≫ |S₁∩S₃|`). Test whether the
   two-parameter family `2^{λ₁|A∩B| + λ₂ f(|A|,|B|)}` (003 lead 3)
   suppresses it: compute its conditional ORs on the 4-atom family (same
   engine, one-line kernel change). If yes, the restated Gap 1 should be
   attacked for the two-parameter class directly.
5. **Equality structure of Prop 2.** Characterize μ whose diagonal ORs all
   equal `2^λ` (proportional slices at every prefix). Conjecture
   (SPECULATION): exactly the product measures on product supports. A
   finite-n proof would pin the `OR ≡ 2^λ` class completely — the tilt
   family's own fixed point.

## References

- This repo: `problems/union-closed/attempts/003-dependent-couplings.md`
  (Gap 1 statement, §3 item 1; lead 1); `004-skeptic-review-of-003.md`
  (corrections adopted throughout: 0.431496 ceiling, restated mini-theorem,
  half-mixing coupling); `002-weighted-kl-ladder.md` (protocol);
  `explore/uc_odds_ratio.py`; `data/005_part[A-G].json`,
  `data/005_full_run.log`.
- S. Karlin, Y. Rinott, "Classes of orderings of measures and related
  correlation inequalities I: Multivariate totally positive distributions",
  J. Multivariate Anal. 10 (1980) 467–498 — MTP₂ closure under restriction
  and marginalization. Standard statements cited from memory, not re-fetched
  this cycle [T].
- Plackett copulas: standard 2×2 fixed-margin odds-ratio family, as in 003.
- Sinkhorn scaling existence/uniqueness for strictly positive kernels:
  classical (Sinkhorn 1967 / Csiszár I-projection), as already relied on in
  003/004.
