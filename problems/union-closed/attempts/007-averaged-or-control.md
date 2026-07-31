# 007 — Averaged odds-ratio control: well-posed, partially proved, then refuted in-regime

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed
- **Type:** formalization + proof-gap attack on the restated Gap 1
  (`averaged-odds-ratio-control`, 005 lead 1 as amended by 006's zero-mass
  bookkeeping note). Outcome: the well-posed statement is **REFUTED** by an
  explicit 10-atom measure inside the record-relevant marginal regime; four
  structured subclasses and the i-aggregated form survive, three of them with
  proofs.
- **Tools:** `explore/uc_or_avg.py` (written and run here; standard library
  only; deterministic, fixed seeds; runtime ~2 min; checkpoints
  `data/or_avg_part[WCKLSTU].json`, log `data/or_avg_run.log`). Command
  reproducing every number below:
  `python problems/union-closed/explore/uc_or_avg.py 2>&1 | tee problems/union-closed/data/or_avg_run.log`
  Cross-engine checks call `explore/uc_or_skeptic.py` (006's independent
  implementation) on the exhibited witness.
- **Sources:** 003/004/005/006 (this repo). Only 004-corrected statements
  used where 003/004 disagree. No external fetches this cycle.

Notation as in 005/006: tilt coupling `π_λ(A,B) ∝ u(A)u(B)2^{λ|A∩B|}`, both
marginals μ; histories `a = A_{<i}`, `b = B_{<i}`; `OR_i(a,b)` the
conditional odds ratio; `m = n − i`; slices `F^a_α(x) = u(a,α,x)`;
`W = ⊗_{j>i}[[1,1],[1,2^λ]]`; record threshold 0.38271 (Liu).

## Approach

**The target.** 005 lead 1, the replacement for the refuted pointwise Gap 1:
`M_i = E_{(a,b)~π}[log₂ OR_i(a,b)] ≥ λ` for all μ and all i (SPECULATION as
recorded there; survived its first falsifiable test in 006 S8). 006 flagged
that the statement is not yet well-posed: on 4-atom supports every history at
i ≥ 3 is degenerate, so the average needs stated conventions. Tasks in order:
(T1) make M_i well-defined for ALL μ, i, with each convention *forced* rather
than chosen; (T2) wide kill attempt before any proof effort; (T3) proof for a
precisely-stated class or a precise obstruction.

**Why structure theory again rather than direct FKG/induction.** 005's slice
reduction makes every question about conditional ORs a question about free
nonnegative vectors (any positive u is its own marginal's potential). Writing
`M_i − λ` as a *pairing of two matrices over prefixes* — a mass matrix m and
a deficit matrix L — decides the problem in both directions at once: the
matrix identities that hold become partial theorems (Theorems A–C below), and
the matrix inequality that fails points exactly at where a counterexample
must live. It did: the counterexample was constructed from the failure mode
of the pairing, not found by blind search (the blind searches — part K, S —
all pass).

## What was done

### 1. Well-posedness (T1): the forced definition

**Definition.** For μ on 2^[n] with H(μ) > 0 and λ ≥ 0, let π = π_λ. For
1 ≤ i ≤ n let `P_i` be the reachable (i−1)-bit prefixes of supp(μ) and

    N_i = { c ∈ P_i : both continuations c0, c1 are reachable in supp(μ) }.

For (a,b) ∈ N_i × N_i the conditional 2×2 table has all four cells positive
and `OR_i(a,b) ∈ (0,∞)`. Define, when N_i ≠ ∅ (else M_i is undefined and
any bound is vacuous at that i):

    M_i(μ,λ) = Σ_{(a,b)∈N_i²} m_ab · log₂ OR_i(a,b)  /  Σ_{(a,b)∈N_i²} m_ab ,
    m_ab = π(A_{<i} = a, B_{<i} = b).

**Degeneracy dichotomy (proved).** For the tilt coupling, cell (α,β) of the
table at (a,b) is positive iff slices `F^a_α` and `F^b_β` are both nonzero
(every atom pair contributes strictly positive kernel mass). So positivity of
a cell factorizes through the two prefixes separately: a positive-mass
history is either fully nondegenerate — exactly when a, b ∈ N_i — or has a
*deterministic conditional margin*. There is no mixed pattern with OR = 0 or
∞ at positive mass. Verified: 1016 coordinate-censuses over 200 random sparse
supports, zero exceptions (part W1).

**Why each convention is forced.**

- *Zero-mass histories:* excluded automatically by the π-weighting — no
  choice exists.
- *Margin-degenerate histories:* the odds ratio there is **unidentifiable**,
  not extreme — every ρ ∈ (0,∞] produces the same table when a margin is 0
  or 1, so the tilt makes no Plackett choice at such a history and no value
  is being hidden by exclusion. Assigning any surrogate value is refuted by
  continuity: at a degenerate history of the 4-atom crash family (n=5, i=3,
  a=b=(0,1)), the ε→0 limits of log₂ OR along three full-support mixing
  directions are +1.4716, +1.0399, +1.0000 (part W2) — direction-dependent,
  so **no assignment whatsoever is continuous**; conditioning-out is the only
  convention compatible with the ε-limits it does admit.
- *Normalized (conditional) mean rather than unnormalized sum:* forced by the
  two proven anchors. At i = n, `OR ≡ 2^λ` (005 Prop 1) must give M_n = λ
  exactly on every support; the unnormalized sum gives λ·(defined mass) < λ
  whenever any prefix is deterministic. Same for the product-μ identity.
- N_i² is a *principal block* of the history-mass matrix, which is what makes
  the definition compatible with the pairing structure below.

**Reduction to prior usage:** part C recomputes 006 S8 (crash + mirror M_i)
to 8.9e-16 and 005 part E's mass-weighted means to 2.8e-14 with this
definition — it is exactly what 005/006 computed where they computed it.

### 2. Structure theory: what is now proved (T3 content, obtained first)

All statements below are elementary given 005's reduction `OR = 2^λ·R`; each
was also checked numerically (parts L, T).

- **Lemma G (mass matrix is Gram).** `m_ab = ⟨w_a, K w_b⟩` with
  `w_a = u·1_{cylinder a}` and `K = 2^{λ|A∩B|} = ⊗_j [[1,1],[1,2^λ]] ⪰ 0`.
  So m restricted to N_i² is PSD with positive entries.
- **Lemma M (off-diagonal mass suppression).**
  `m_ab = 2^{λ|a∧b|}·S_ab` with S a Gram matrix, and
  `|a∧b| = (|a|+|b|−|aΔb|)/2`, so
  `m_ab ≤ 2^{−λ|aΔb|/2}·√(m_aa m_bb)` — exponential decay in prefix
  Hamming distance.
- **Theorem A (Frobenius form).** With `L_ab = log₂ OR_i(a,b) − λ` (the
  deficit matrix, symmetric, diagonal ≥ 0 by 005 Prop 2):
  `M_i − λ = ⟨m, L⟩_F / ⟨m, J⟩_F` over N_i². Hence **if L ⪰ 0 then
  M_i ≥ λ** (Frobenius pairing of PSD matrices is nonnegative). Identity
  verified to 9.9e-15 over 572 (instance, i) pairs (part L).
- **Theorem B (point-mass class, includes all ≤4-atom supports).** If every
  cell (a,α), a ∈ N_i, contains a single atom, then
  `L_ab = λ⟨d_a, d_b⟩` with `d_a = x^a_1 − x^a_0` the future-difference
  vectors — a Gram matrix, so L ⪰ 0 and **M_i ≥ λ**. Any support of ≤ 4
  atoms has, at every i, either |N_i| ≤ 1 (diagonal only, 005 Prop 2
  applies) or singleton cells — so **averaged control is a theorem for all
  supports of at most 4 atoms**. This is exactly why 006 S8's crash/mirror
  test passed, and both families sit ON the PSD boundary
  (`ℓ² = e_a e_b`), explaining the sharpness observed there.
- **Rigidity.** If a diagonal excess vanishes (`F^a_1 ∝ F^a_0`) the entire
  row of L vanishes: zero diagonal excess forces zero off-diagonal deficit.
- **Theorem C (first order in λ; modulo one standard step).** For fixed μ,
  `M_i(λ) − λ = λ·‖Σ_{a∈N_i} ν(a) D_a‖² + O(λ²)`, where ν is the normalized
  prefix mass and `D_a = E_μ[future | a,1] − E_μ[future | a,0]`: the
  first-order coefficient is a **perfect square**, so first-order-in-λ
  counterexamples do not exist. (Uses differentiability of the Sinkhorn
  potential in λ at λ=0 — the same implicit-function-theorem step 005 Prop 6
  takes as standard; SPECULATION only in that narrow sense.) Verified on the
  witness below: excess/λ → +4.65e-4 = the independently computed
  ‖Σν(a)D_a‖² (part T(e)).
- **MTP₂ subclass.** Potential-level MTP₂ gives pointwise `OR ≥ 2^λ`
  (005 Prop 5), hence M_i ≥ λ trivially on that class.

**Where the proof route breaks (precisely).** L is **not** PSD in general:
132 of 572 sampled (instance, i) pairs have a negative eigenvalue (min
−3.11, part L), and even the 2×2-minor control `|L_ab| ≤ √(e_a e_b)`
(e_a = diagonal excess) is false for the tilt kernel — free-slice
configurations reach `|L_ab| ≈ 6.7·√(e_a e_b)`. So every assembly that
bounds off-diagonal deficit against diagonal excess by a Cauchy–Schwarz-type
inequality (Frobenius, geometric-mean + Lemma M + Schur test) fails on
configurations with *both slice pairs nearly proportional but strongly
anti-correlated across blocks*. Yet on every sampled non-PSD instance the
mass weighting still rescued the average (min excess +5.2e-4) — the
counterexample needed the mass–deficit correlation to be attacked directly.

### 3. Kill census (T2): every known genre passes

444 instances, 0 violations (part K; tolerance 1e-9; the i = n identity
floors the all-i minimum, so margins are reported over i < n): crash + mirror
at n ∈ {4,5,6,8} and six λ including off-grid (0.05, 0.31, 0.73, 1.37, 2.6);
**ε-mixed crash** at ε down to 1e-5 — the i ≥ 3 histories 006 flagged as
zero-mass now carry positive mass and their M_i stay above λ;
δ_∅ ⊕ Bern(½+ε) (004's killer genre); Chase–Lovett slice and smoothed slice;
Sawin-genre geometric mixtures with components up to 0.92; 300 random sparse
supports. Closest approaches at i < n: Sawin mixtures +1.5e-4, the rest
≥ +6.6e-3. Hill-climbs on generic sparse supports (part S, 36 restarts)
converge onto the boundary M_i = λ from above (endpoints −4e-16, i.e.
equality configurations) and never cross. The averaged control is *sharp and
generically true* — which is what made the following counterexample
invisible to every prior test.

### 4. The counterexample (headline): REFUTED in-regime

The |N_i| = 2 case of Theorem A's pairing is, after the slice reduction, an
inequality in four free nonnegative vectors and two integer prefix weights.
Hill-climbing **that** object (not μ-space) found robust violations, first at
max marginal ≈ 1; two exact invariances then force them into the record
regime:

- *Stripping:* a coordinate equal to 1 in every atom multiplies the kernel
  by a constant and cancels — M_i unchanged, marginal-1 coordinates removed.
- *Dilution invariance (proved and verified to 1.8e-15):* adding atoms whose
  prefix differs from both active prefixes and whose response bit is
  deterministic leaves N_i, every (a,b)-table, and all relative history
  masses unchanged — M_i is **exactly** invariant — while diluting every
  marginal. Spreading the dilution over three singleton atoms keeps the
  dilution coordinates' own marginals at ≈ 0.128.

**The witness** (part T; n = 7, λ = 3.5, response coordinate i = 5;
potential u = its own marginal's Sinkhorn potential, 005's trick):

    atom      u        atom      u        atom      u
    {1,6}   3.8027    {}      0.25267    {2}     20
    {1,7}   8.0947    {6}    27.286      {3}     20
    {1,5,7} 0.0047930 {7}    12.530      {4}     20
                      {5,6}   0.36092

    M_5 − λ = −0.122033        (VIOLATION of M_i ≥ λ)
    all elementwise marginals ≤ 0.3178 < 0.38271   (in-regime)
    H(μ) = 2.72 bits; defined history mass at i=5: 0.551
    all M_i − λ:  +1.48  +1.91  +2.24  +2.77  −0.122  +3.47  0 (= i=n)

Verification battery (all in part T): (a) Sinkhorn round-trip — refitting
the tilt from μ alone reproduces every M_i to 3.4e-14 (and 006's independent
engine `uc_or_skeptic.py` gives −0.122033 identically, residual 7e-15);
(b) dilution invariance −0.122033 with/without the three dilution atoms
(diff 1.8e-15); (c) 20/20 random 3% perturbations of all ten weights remain
in-regime violations; (d) rounding every weight to **two significant
digits** still gives −0.1219 at max marginal 0.313; (e) the *same μ*
violates for every λ from ≈0.03 up through 5.0 (−0.004 at λ=0.5, −0.128 at
λ=4), and is positive below λ ≈ 0.03 exactly as Theorem C forces;
(f) sweeping the light slice's weight over [1e-6, 5e-2] moves the violation
only in the fourth decimal — this is a stable limit family, not a tuning
accident. The clamped search (part U, single-i mode) independently
rediscovers an in-regime violation (−0.0033) from random starts.

**Mechanism.** Two blocks: block a = sets containing element 1, block b =
sets avoiding it. Block b's response-slices are *nearly proportional*
(e_b small ⇒ almost no diagonal excess to spend) yet its slices
anti-correlate with block a's through the tilt (ℓ = log₂R(a,b) ≪
−√(e_a e_b) — the 2×2-minor failure), and block a's response-1 slice is
*light* (weight 0.0048), which props the off-diagonal histories open at Θ(1)
mass while the anti-correlated cross-ratio sets their OR. It is 005's crash
mechanism (future-overlap discrepancy) resurfacing one level up: what the
crash did to a pointwise OR, the light-slice two-block geometry does to the
mass-weighted average.

**A float-artifact detour worth recording** (methodological lesson). The
first "witnesses" the unclamped hill-climbs produced were IEEE artifacts: an
atom of weight ~1e-299 made u(A)u(B) underflow to exact 0 in one *diagonal*
cell, silently reclassifying a hugely positive diagonal history as
degenerate and dragging the average down. **Both engines agreed on the fake
number** (they share IEEE arithmetic), and the Sinkhorn round-trip
reproduced it too — cross-implementation agreement is no defence against a
shared floating-point pathology. It was caught by physically sweeping the
light atom's mass: the honest family was positive. Guards now in the
engine: a cell-floor diagnostic (min cell / history mass; the witness's is
3.0e-6, far from underflow — atom weights span 4.8e-3..27) and a
dynamic-range clamp in all searches.

### 5. What survives the refutation (checked)

- **The i-aggregated form** `Σ_{i<n} w_i(M_i − λ) / Σ w_i ≥ 0` (w_i =
  defined mass): the witness itself gives +1.84 (its other coordinates
  over-perform), and the clamped aggregate search (part U) found no
  in-regime violation (lowest endpoint +0.0125). Unrefuted — now the
  natural restatement, since the chain-rule assembly sums over i anyway.
- **i = 1 and i = n:** M_1 ≥ λ (single diagonal history, 005 Prop 2) and
  M_n = λ (005 Prop 1) are theorems; violations can only live at
  2 ≤ i ≤ n−1.
- **≤ 4-atom supports** (Theorem B), **potential-MTP₂ μ** (via Prop 5), and
  **first order in λ at fixed μ** (Theorem C, modulo the standard smoothness
  step): proved.

## Outcome

**REFUTED — the restated Gap 1 (`averaged-odds-ratio-control`,
`M_i ≥ λ for all μ and all i`) is false, including in the record-relevant
regime.** Scope: the well-posed statement of §1 (which provably matches
005/006's usage wherever they computed) fails at an explicit μ on 10 atoms
of 2^[7] with all elementwise marginals ≤ 0.3178 < 0.38271 and H(μ) = 2.72
bits: `M_5 = λ − 0.122` at λ = 3.5, with the same μ violating for every
λ ≳ 0.03. The witness is verified by two independent engines to 1e-14,
stable under 3% perturbations and 2-digit rounding, and reproducible via
`python problems/union-closed/explore/uc_or_avg.py` (part T).

**Not claimed:** no claim that the dependent-couplings route dies — as with
005, this kills a candidate *bridge*, not the interface (licensing lemma,
separations, 0.431496 ceiling untouched); no claim that the i-aggregated
form is true (only that it survived this cycle's searches and the witness);
no claim about μ-level MTP₂ (still open, 005 lead 2); Theorem C rests on the
same unproved Sinkhorn-smoothness step as 005 Prop 6; none of this has had
its independent skeptic pass yet — per repo rules these are results only
after one. Calibration note: the refutation is the part I trust most (two
engines, exact invariances, perturbation-stable); a skeptic should
nevertheless attack the witness first (below), then Theorem B's cell-count
argument, then the forcing arguments of §1, which are arguments of
inevitability rather than formal theorems.

## Why it failed / what survived

**The quantity that goes the wrong way.** In the Frobenius form
`M_i − λ = ⟨m, L⟩/⟨m, J⟩`, the mass matrix m is always PSD (Lemma G) with
exponentially suppressed off-diagonal (Lemma M) — but the deficit matrix L
is not PSD, and the pairwise bound `|L_ab| ≤ √(e_a e_b)` that would let the
suppression win fails by unbounded factors for mixture slices. The
counterexample genre is exactly the certificate of that failure: *small
diagonal excesses with large anti-aligned cross-ratios, propped open at Θ(1)
history mass by a light slice*. Any future averaged-control claim for this
tilt class must either exclude that geometry or exploit the m–L correlation
beyond Cauchy–Schwarz; pure spectral/PSD reasoning cannot close it.

**Why every earlier test missed it** (and 006 S8 passed): the genre needs
≥ 5 atoms with a specific two-block slice structure; all ≤ 4-atom supports
are safe by Theorem B, product mixtures and slice families have too much
symmetry (their D-vectors align), and random sparse supports hit the safe
generic mass-L correlation.

Survived / reusable:

- The **well-posed M_i definition** with the degeneracy dichotomy — any
  future averaged statement on tilt couplings should reuse it verbatim.
- **Lemmas G, M, rigidity; Theorems A, B, C** — permanent structure theory
  for the tilt class, and the partial positive results the route can still
  build on (≤4 atoms; MTP₂; first order in λ; i ∈ {1, n}).
- The **10-atom witness family** (with its 2-digit tidy version) as the
  hard instance for ANY future averaged or smoothed OR-control claim; the
  light-slice + two-block construction generalizes.
- The **strip + fresh-prefix-dilution invariance** trick: converts any
  two-block violation into an in-regime one — marginal constraints cost
  nothing against M_i-type functionals, which should temper enthusiasm for
  marginal-restricted restatements generally.
- `explore/uc_or_avg.py`: census engine with degeneracy bookkeeping, L-matrix
  spectral probe, two-prefix reduction searches, cell-floor artifact
  detector; cross-validated against both prior engines at 1e-14.
- The **underflow lesson**: shared-IEEE agreement between independent
  implementations is not independence; parameter sweeps toward the suspect
  boundary are the cheap discriminator.

## Leads generated

1. **(new Gap 1 candidate) The i-aggregated control.** Prove or kill
   `Σ_{i<n} w_i (M_i − λ) ≥ 0`. It survived part U's search and the witness.
   Falsifiable first steps: extend the two-prefix search space to two
   *response* coordinates sharing one light slice (the natural way to make
   two coordinates crash at once); on the proof side, sum Theorem A over i
   and look for cross-i cancellation in the total pairing — the diagonal
   funding Σ_i Σ_a m_aa e_a is exactly 005 part B's strict-excess budget.
2. **Margin-modulated control (what the assembly actually needs).** The
   chain-rule cost of an OR-deviation at a history is weighted by the
   sensitivity of `h(z_ρ(x̃,ỹ))` to ρ, which vanishes as the conditional
   margins degenerate — precisely where the witness hides its deficit.
   Compute the h-sensitivity-weighted census on the witness: if the weighted
   mean is ≥ λ there and on the census, conjecture
   `E[σ(x̃,ỹ)·(log₂OR − λ)] ≥ 0` with σ the Plackett sensitivity, and
   re-run this cycle's whole battery against it. Definite outcome either
   way, and it is the assembly-relevant statement.
3. **Characterize when 2×2-minor control holds.** It holds for point-mass
   slices and failed for mixtures; the boundary is a statement about
   log-Gram matrices of nonnegative vectors in the infinitely-divisible
   kernel `2^{λ⟨x,y⟩}`. Even a sufficient condition (e.g. slice supports
   laminar) would give M_i ≥ λ for a nontrivial mixture class via Lemma M +
   Schur, quantifying how far Theorem B extends.
4. **Witness vs the actual functional.** Run 003/004's gain accounting on
   the witness μ: does the crashed coordinate produce negative
   per-coordinate entropy gain for the tilt recipe, or does the recipe's
   sup over λ route around it? If the latter, averaged OR-control was
   never necessary for the assembly and Gap 1 should be re-restated at the
   functional level directly.
5. **(bookkeeping) 005 lead 2 (μ-level MTP₂ ⟹ OR ≥ 2^λ) is untouched** and
   remains the cleanest open positive candidate; the witness μ is not MTP₂
   (e.g. μ({2})·μ({3}) > 0 but μ({2,3}) = 0 while ∅ ∈ supp), consistent
   with lead 2.

## References

- This repo: `attempts/003-dependent-couplings.md` (interface, Gap 1);
  `attempts/004-skeptic-review-of-003.md` (corrected ceiling 0.431496,
  half-mixing, restated mini-theorem — used exclusively where 003/004
  differ); `attempts/005-odds-ratio-control-refuted.md` (slice reduction,
  Props 1–6, crash family, lead 1 = the target);
  `attempts/006-skeptic-review-of-005.md` (S8, the bookkeeping note this
  record resolves, λ>0 correction to Prop 2's equality case).
- Tools/data: `explore/uc_or_avg.py`; checkpoints
  `data/or_avg_part[WCKLSTU].json`; log `data/or_avg_run.log`; reviewed
  engines `explore/uc_odds_ratio.py`, `explore/uc_or_skeptic.py`.
- Ahlswede–Daykin and Karlin–Rinott facts are used only through 005/006's
  already-verified statements; no new external sources.
