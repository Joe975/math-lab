# 009 — Gap 2 probed: the mutual-information tax, made precise and tested

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed
- **Type:** proof-gap formalization + first falsifiable computational tests
  (deliberately a small probe, not a campaign), on the labeled gap
  `mutual-information-tax` (Gap 2) of the LIVE dependent-couplings route.
  Follows 003 §3 item 2 / "Why it failed" Gap 2, using 004's corrected
  statements throughout (half-mixing coupling, 0.431496 ceiling) and 005/006
  only for the crash family and the state of Gap 1.
- **Tools:** `explore/uc_mitax.py` (written and run here; standard library
  only; deterministic, no RNG; runtime ~25 s; checkpoints
  `data/mitax_part[A-E].json`, log `data/mitax_full_run.log`). Command
  reproducing every number below:
  `python3 problems/union-closed/explore/uc_mitax.py`
- **Sources:** 002/003/004/005/006 (this repo). No external fetches; all
  outside facts (Gilmer/AHS structure, Sawin gadgets, Liu record) are used
  exactly as already transcribed and verified in 001–004.

Notation as in 003/005: coordinates revealed in fixed order 1..n, histories
`a = A_{<i}`, `b = B_{<i}`; `π` a coupling with BOTH marginals μ;
`ψ = (3−√5)/2`; tilt coupling `π_λ ∝ u(A)u(B)·2^{λ|A∩B|}` (Sinkhorn), pure
tilt `∝ 2^{λ|A∩B|}` on slice×slice (marginals automatic).

## Approach

Gap 2 is recorded in 003 as one sentence: the tax
`Σᵢ I(A_i; B_{<i} | A_{<i})` "must be beaten by the per-coordinate gains; no
analogue of the (KEY) inequality with tax is known". 004 left it untouched;
005/006 explicitly claimed no progress on it. Before anyone attempts the
assembly (008's line), two cheap questions decide whether the gap as stated
is even alive, and both are answerable by computation:

1. Is the tax-paid functional smoothing-sensitive? If the chain-rule
   assembly value (below) obeyed the 002 no-go's `O(log 1/δ)` smoothing
   bound, the gap would be dead on arrival — the proof scheme, not just the
   functional, would have ceiling ψ.
2. Does the tax already eat the route's separations on the known adversary
   genres? A single genre where the coupling gain is positive but the
   tax-paid assembly value is not would restate Gap 2 as an obstruction.

Why this rather than the obvious alternative (starting the perturbative
O(ε²) expansion of 003 lead 1): the expansion is worthless if either check
fails, and both checks needed the same tool — an exact full-history
chain-rule evaluator — which did not exist yet. Cost of the general
alternative avoided: nothing; this is the prerequisite either way.

## What was done

### 1. The statement under test (constructed here — 003/004 leave it vague)

003 names the tax but never writes an inequality. The following is the
sharpest precise candidate consistent with the route's proof scheme; I
constructed it in this attempt and record that explicitly. Define, for a
coupling π of (A, B) with both marginals μ:

    z_i(a,b)  = Pr_π[A_i = 0, B_i = 0 | a, b]        (union bit U_i = 1−z)
    x̃_i(a,b) = Pr_π[A_i = 0 | a, b]
    CR(μ,π)   = Σᵢ E_π[h(z_i)] − H(μ)               (chain-rule assembly value)
    T_A(π)    = Σᵢ I(A_i; B_{<i} | A_{<i})          (the mutual-information tax)
    ST(π)     = Gain(μ,π) − CR(μ,π)                 (the second tax, new here)

Four one-line facts, each re-derived and verified numerically in part A:

- (i) `H_π(U) ≥ Σᵢ E_π[h(z_i)]`: chain rule for U plus "conditioning on the
  finer σ-field (a,b) ⊇ σ(U_{<i}) reduces entropy". Hence **CR ≤ Gain** and
  `ST = Σᵢ I(U_i; (A_{<i},B_{<i}) | U_{<i}) ≥ 0`.
- (ii) `Σᵢ E_π[h(x̃_i)] = H(μ) − T_A(π)` (identity), so
  **CR = Σᵢ E_π[h(z_i) − h(x̃_i)] − T_A(π)**: "per-coordinate gains minus
  tax". Gap 2's sentence is exactly the assertion CR > 0.
- (iii) `T_A(π) ≤ I_π(A;B)` (chain rule: `I(A_i;B|A_{<i}) ≥
  I(A_i;B_{<i}|A_{<i})`, summed).
- (iv) CR is what the route's proof scheme actually produces: every
  (KEY)-style per-coordinate argument works with the (a,b)-conditional bit
  laws, so ST is the exact price of certifying the gain by chain rule
  rather than evaluating H_π(U) directly. ST was not named in 003; it is a
  second, distinct loss channel of Gap 2.

  **(TAX at p):** for every μ with H(μ) > 0 and all element marginals < p,
  the route's recipe coupling π(μ) has CR(μ, π(μ)) > 0.

Here "recipe" is per-genre, exactly as the route currently stands (Gap 3
open): λ-swept tilt on generic/slice μ, block-adaptive shared-component
coupling on product mixtures, half-mixing on the {0}∪[1/2,1) genre.
(TAX at p) ⟹ (S-coup at p) by (i), ⟹ Frankl at p by the licensing lemma
(003, verified in 004).

The functional `μ ↦ CR(μ, π(μ))` is NOT a functional of `(law(U), μ)` — it
depends on the joint conditional structure of π — so the 002 no-go does not
formally cover it; whether it is *quantitatively* smoothing-sensitive is
part D's question. Both tax functionals respond to the coupling's history
structure, not to log-likelihoods of μ, so no smoothing-insensitivity is
built in — but that had to be demonstrated on the no-go's own certificates,
not asserted.

### 2. The evaluator, and its validation (part A)

`cr_eval` computes every quantity above exactly for an explicit coupling
(list of atom-pairs with weights), by grouping pairs on the history pair
(a,b) at each coordinate. Validation:

- Bernoulli products under per-coordinate Plackett(ρ), n = 8, p = 0.4,
  ρ ∈ {1, 1.5, 3}: CR matches the closed form `n(h(z_ρ(x,x)) − h(p))` to
  1e-9, and T_A = ST = 0 to 1e-12 (as they must be: the coupling
  factorizes over coordinates). This ties the evaluator to 003 part A.
- Pure slice tilt n=8, w0=3, λ=0.7: gain +1.585790908 matches the
  brute-force value in the checked-in `data/003_partB2.json` to 3e-13.
- 4-atom crash family: the fitted Sinkhorn coupling's conditional table at
  coordinate 2, history (A₁=0, B₁=1), reproduces 005's exact crash identity
  `log₂ OR = λ(3−n)` to 1.8e-15 at every λ and n tested — an independent
  check of the evaluator's conditioning logic against a known exact value.
- Slice DP (part E) against the atom evaluator at n=8: |ΔCR| = 4.7e-14;
  internal consistency (state masses sum to 1; cell sums reproduce the
  future count identity) ≤ 1.1e-13.

### 3. Smoothing-sensitivity: the 002-no-go check (part D)

**Small lemma (rigorous).** For the 003 part-D block coupling on a mixture
`μ = Σ_k P_k Bern(p⃗_k)^⊗` (share the latent k; per coordinate, union
marginal m_k): `E[h(z_i)] = H(U_i|a,b) ≥ H(U_i|a,b,k) = Σ_k P_k h(m_k)`,
because given k the pair (A_i,B_i) is independent of both histories. Hence

    CR_block(μ) ≥ n·Σ_k P_k h(m_k) − H(μ),

with H(μ) computed exactly (profile arithmetic). On 002's certificate
instances — the exact gadgets that killed every smoothing-insensitive
functional — with H(μ) re-derived here and matching 002's table (1927.70 /
19239.59 / 57578.99):

    n      ū       θ      max marg   H(P)    CR lower bound
    2000   0.390   0.05   0.402269   0.301      +66.58
    20000  0.386   0.02   0.390799   0.144     +741.33
    60000  0.3823  0.001  0.382536   0.011    +2418.55

(Each within H(P) of 003's part-D gain values, as the lemma predicts.) The
tax-paid chain-rule functional is Θ(n) POSITIVE on exactly the instances
where every functional covered by the 002 no-go is Θ(n) negative: it
responds at probability scale. **Not dead on arrival; the probe continues.**
Note the tax never had to be estimated here — conditioning on k routes
around it; this is the chain-rule shadow of Liu's "pay the tax with
auxiliary-mixture bookkeeping" observation in 003 §2.

### 4. Adversarial tests of (TAX at p)

**(a) 4-atom crash family (005 Prop 3), Sinkhorn tilt, λ sweep** (part B;
n ∈ {5, 8, 11}, λ ∈ [−1, 3]; marginals ≤ 0.37, H(μ) = 1.76). Headline:
**sup_λ CR > 0 at every n — the candidate survives its most dangerous
instance.** Detail worth recording:

    n=5 : best CR +0.1307 at λ=0.5;  CR > 0 for λ ∈ [−0.5, 2.0]
    n=8 : best CR +0.1281 at λ=0.25; CR > 0 for λ ∈ [−0.25, 0.5]
    n=11: best CR +0.1154 at λ=0.25; CR > 0 for λ ∈ [0, 0.5]

  The iid point CR(λ=0) = +0.1114 and Gain(λ=0) = +0.9390 are n-independent
  (coordinates ≥ 3 are deterministic given the first two, so only two
  coordinates ever contribute). The positive-CR window in λ **narrows
  toward 0 as n grows** — the OR crash `2^{λ(3−n)}` makes any fixed λ > 0
  fail the chain-rule certificate at large n, while λ = 0 always certifies
  this family. Gain itself stays positive at every λ tested (+0.28 to
  +1.00): on an O(1)-entropy family the second tax ST ≈ 0.3–0.9 bits is
  the binding loss, and it is the crash mechanism (not T_A ≤ 0.05, which
  stays tiny) that moves CR. A fixed-positive-λ recipe therefore fails
  (TAX) on crash-type μ; the λ-sweep recipe does not.

**(b) Half-mixing genre (004 R1): μ = P₁δ_∅ ⊕ (1−P₁)Bern(p_h)^n** (part C;
(p_h, P₁) = (0.6, 0.5) — 004's unshielded μ_n, marginals 0.30 — and
(0.51, 0.163) — the recipe-killer shape at cap 0.427; n = 6..14).
**ST = 0 to machine precision (≤ 1.6e-12) in all 10 runs**, i.e. CR = Gain
exactly, and T_A ≈ 0.108 / 0.126 constant in n. This is an identity:

  **Mini-lemma HM (proved here; pending skeptic pass).** For the
  half-mixing coupling, `σ(U_{<i})` and `σ(A_{<i}, B_{<i})` induce the same
  conditional law of U_i, so ST = 0 identically. *Proof.* With A = ε_A s,
  B = ε_B s: if U_{<i} = u ≠ 0 then ε_A ∨ ε_B = 1 and s_{<i} = u are
  forced, and in every cell consistent with any (a,b) refining u, U_i =
  s_i ~ Bern(p_h) independent of the past; if U_{<i} = 0 then (a,b) = (0,0)
  — the same σ-field atom. ∎

  Consequence: 004's closed-form gain `(P₁/2)·n·h(p_h) − [h(P₁) − h(P₁/2)]`
  transfers verbatim to CR (numerically confirmed: CR matches it to 4
  decimals by n = 12): **on the {0}∪[1/2,1) 2-block genre the tax-paid
  statement holds with Θ(n) margin, by proof, not just by sweep.**

**(c) Chase–Lovett pure slice, pure overlap tilt** (part E; exact DP on the
sufficient prefix state (w_a, w_b); n ∈ {24, 40, 60, 80},
p ∈ {0.42, 0.3823}, λ ∈ [0, 3]). Two findings:

  **Mini-lemma SL (proved here; pending skeptic pass): the pure slice tilt
  is tax-free — T_A = T_B = 0 identically.** *Proof.* For A extending
  prefix a inside slice(w0), the future part A′ = A ∩ {i..n} has fixed size
  w0 − w_a, and `Σ_{B ext b} ρ^{|A∩B|} = ρ^{|a∧b|} · G(n−i+1, w0−w_a,
  w0−w_b)` by exchangeability of the future block — the same constant for
  every A extending a. So the conditional law of A given (a,b) is uniform
  over slice-completions of a, independent of b: `A_i ⊥ B_{<i} | A_{<i}`
  exactly. ∎ (Numerically: |T_A| ≤ 4e-12 across all 42 slice runs, every λ.)
  The entire assembly loss on slices is therefore the second tax.

  **ST grows like 0.40·log₂n** (EVIDENCE, this range): at λ = 3, p ≈ 0.42,
  successive slopes ΔST/Δlog₂n = 0.3997, 0.4065, 0.4007 across
  n = 24→40→60→80; ST ∈ [1.74, 2.53] bits over the whole sweep. sup_λ CR >
  0 at every (n, p) tested — at p ≈ 0.3823 comfortably (+2.63 at n=60,
  +3.21 at n=80), at p ≈ 0.42 by an increasing margin in n (+0.56 at n=24
  to +1.02..1.19 at n=80/60) though the needed λ sits at the top of the
  grid (λ ≈ 2–3; iid CR is −4.75 at n=80, so the tilt genuinely does the
  rescuing here). Since CR = Gain − ST with ST = O(log n) (evidence) and
  the tilt gain on slices is Θ(n) at these p (verified to n = 2·10⁶ in
  003/004), CR inherits the Θ(n) separations wherever the ST scaling holds.

**Union-closed control:** none needed beyond 004's observation that
controls are tautological here — CR ≤ Gain ≤ 0 on exactly union-closed
uniforms by the licensing lemma, for every coupling; a positive CR on one
would be a bug in (i), not new information. The part-A closed-form checks
play the code-check role instead.

## Outcome

**EVIDENCE.** Scope of what was checked, precisely:

- The candidate (TAX at p) — constructed in §1 because 003/004 state Gap 2
  without an inequality — **survives its first falsifiable tests on all
  four adversary genres of record**: 4-atom crash family (n ≤ 11, λ ∈
  [−1,3]), half-mixing genre (n ≤ 14, two parameter points, exact identity),
  Sawin certificate gadgets (n ≤ 60000, block recipe, rigorous bound),
  Chase–Lovett pure slice (n ≤ 80, λ ∈ [0,3], p ∈ {0.42, 0.3823}).
- The tax-paid functional is smoothing-sensitive: Θ(n)-positive on the 002
  no-go's own certificate instances (part D, rigorous given float
  arithmetic; H(μ) cross-matched against 002). The gap is not dead on
  arrival.
- Two structural identities were found and proved (pending skeptic pass):
  the pure slice tilt pays zero mutual-information tax (SL), and the
  half-mixing coupling pays zero second tax (HM). Fact (iii) T_A ≤ I(A;B)
  is proved trivially. The block-coupling CR bound of part D is proved.
- The second tax ST — a loss channel Gap 2 did not name — is identified,
  measured (0 / O(1) / ≈ 0.40·log₂n by genre), and located as the binding
  loss on both the crash family and the slices.

**Not claimed:** no bound, no constant, no proof of (TAX at p) for any p —
the tests are finite-n and finite-λ-grid; "survives" is EVIDENCE about the
instances listed, never about the statement. The ST ≈ 0.40·log₂n scaling
is a 4-point fit at one λ and one p-band. The Sawin genre was tested only
with the block recipe (the tilt recipe's CR on large-n gadgets is
unevaluated — the O(n⁵) potential-DP was out of probe budget). The
mini-lemmas have had no independent refutation attempt yet and per repo
rules are not results until they do. Nothing here touches Gaps 1 or 3.

## Why it failed / what survived

Nothing failed: the probe's expected outcome was a precisely-described
obstruction, and instead every genre passed — but the probe sharpened where
the real difficulty sits. The honest map of what a proof of (TAX at p)
still needs, link by link:

- **L1 (known).** Licensing lemma: (S-coup at p) ⟹ Frankl at p. Verified
  in 004.
- **L2 (known, this attempt).** CR ≤ Gain, so (TAX at p) ⟹ (S-coup at p):
  chain-rule certification is sound, and parts B–E say it is not vacuous
  on the adversary genres.
- **L3 (open — Gap 3).** A total recipe μ ↦ π(μ). Sharpened here: the
  crash family forces the recipe's λ-sweep to include λ ≈ 0 (fixed λ > 0
  fails CR at large n even where it separates in gain), so recipe totality
  and tax control cannot be decoupled.
- **L4 (open — Gap 1 restated).** Per-coordinate control: some averaged
  odds-ratio statement (`M_i ≥ λ` flavor, 005 lead 1 / 007's line) strong
  enough to lower-bound `E[h(z_i)] − E[h(x̃_i)]` coordinate-wise. Both
  directions of deviation matter: part B shows the crash (OR below λ)
  destroying CR through ST, while 005 part E's over-tilting (OR above λ)
  under-credits the tax.
- **L5 (partially known, this attempt).** Tax control by genre: T ≡ 0 on
  pure slices (SL, proved); T bypassed entirely on shared-latent block
  couplings (part-D lemma, proved); T_A ≤ I(A;B) always (proved);
  T_A = O(1) observed on half-mixing. **Open:** T_A of the Sinkhorn tilt
  on general μ — SPECULATION (from 005 Prop 6 mechanics): near product μ
  it is O(λ²) + O(δ²) per coordinate; this is the input 008's perturbative
  assembly needs, and nothing here proves it.
- **L6 (open — the core of Gap 2).** The (KEY)-with-tax inequality: a
  per-coordinate lower bound on `E[h(z_i)] − E[h(x̃_i)]` minus tax terms
  whose sum is positive for all μ with marginals < p, p > 0.38271. No
  analogue known; the perturbative expansion at the AHS equality point is
  the concrete next step (008).
- **L7 (new — named here).** The second tax ST must be afforded: proofs
  via CR give away `Σᵢ I(U_i; (a,b) | U_{<i})`. Affordable when gains are
  Θ(n) and ST = O(log n) (slices, evidence), free on half-mixing (HM),
  but the binding loss on O(1)-entropy families (crash: ST ≈ 0.83 vs gains
  ≤ 1.0). Any assembly aiming at the extremal (near-zero-gain) regime must
  either bound ST or interleave the chain rule more cleverly.

Survived / reusable:

- `explore/uc_mitax.py`: exact full-history CR/tax evaluator for arbitrary
  explicit couplings (validated four independent ways), the slice
  (w_a, w_b) DP with its future-count identity check, and the Sawin
  profile-arithmetic bound — all stdlib, all deterministic.
- The CR/T_A/ST accounting frame itself: Gap 2 now reads "prove
  sup-recipe CR > 0", with the loss decomposition CR = (per-coordinate
  gains) − T_A and Gain − CR = ST, each term computable exactly on any
  explicit instance.
- Mini-lemmas SL and HM as permanent structure: two of the route's three
  per-genre couplings provably lose nothing to one of the two tax channels
  each; the "any positive u is its own marginal's potential" trick (005)
  now has a tax-side sibling: exchangeable-future arguments kill T_A on
  single-slice supports.
- The crash-family λ-window narrowing as the hard-instance phenomenon for
  fixed-tilt chain-rule certification.

## Leads generated

1. **Skeptic pass on this attempt** (next in queue by repo rules):
   re-derive SL and HM independently, re-implement CR with a different
   grouping (e.g. natural-log, B-side histories) and re-run parts B/C/E;
   extend part B in n (n = 14, 17 at λ ∈ {0.1, 0.25}) and check the
   positive-CR window scales like λ_max ≈ c/n against the crash exponent
   λ(3−n) — a definite scaling law either holds or does not.
2. **Prove ST = O(log n) for the pure slice tilt** (currently a 4-point
   fit, slope ≈ 0.40·log₂n). Candidate route: j-concentration plus a local
   CLT for (w_a, w_b) around their conditional means — the state
   fluctuation is O(√n), the bit-probability shift O(1/√n), the per-
   coordinate entropy loss O(1/n), summing to O(log n) only through the
   variance-tracking terms. If instead ST = Θ(n^α), α > 0, chain-rule
   certification dies on slices in the extremal regime — either outcome is
   decisive for L7.
3. **Tax of the Sinkhorn tilt near product μ** (feeds 008 / L5): expand
   T_A to second order in (λ, δ) at the AHS equality point and check
   whether the O(λ²) tax coefficient is beaten by the O(ε) diagonal-
   calibration slack of 003 lead 1 under the worst-case OR-deviation
   budget that 006 showed is family-dependent. Sign of one net
   coefficient; definite either way.
4. **Smoothed slice tax:** SL's proof breaks the moment the support leaves
   a single slice. Fit the Sinkhorn tilt on the 002-part-D smoothed slice
   (t = 2^{−8}) at n ≈ 30 and measure T_A(n, t): if T_A = Θ(n·t) rather
   than O(1), interpolating between the slice and product genres is where
   the tax first bites, and L6 must handle it there.
5. **Close the one untested cell:** tilt-recipe CR on Sawin gadgets at
   adversarial n (the exchangeable-potential (w_a, w_b) DP with g-weights;
   O(n⁵) naive, likely reducible). If tilt-CR is negative at large n where
   block-CR is positive, the recipe for the mixture genre must stay
   block-adaptive and Gap 3's totality problem inherits a genre boundary —
   worth knowing before 008 commits to the tilt.

## References

- This repo: `problems/union-closed/attempts/003-dependent-couplings.md`
  (Gap 2 statement; part-D block coupling; §2 tax remark on Liu);
  `004-skeptic-review-of-003.md` (half-mixing coupling R1; corrected
  ceiling; licensing-lemma verification); `002-weighted-kl-ladder.md`
  (no-go hypotheses and certificate instances; adversarial-test-first
  protocol); `005-odds-ratio-control-refuted.md` (crash family Prop 3,
  used as adversary; OR identity used as evaluator cross-check);
  `006-skeptic-review-of-005.md` (state of the averaged-control gap);
  `explore/uc_mitax.py`; `data/mitax_part[A-E].json`,
  `data/mitax_full_run.log`; prior data cross-checked against:
  `data/003_partB2.json`, `data/003_partD.json`, 002's certificate table.
- Gilmer arXiv:2211.09055; Alweiss–Huang–Sellke arXiv:2211.11731; Sawin
  arXiv:2211.11504; Chase–Lovett arXiv:2211.11689; Liu arXiv:2306.08824 —
  all as transcribed and verified in 001–004; no new fetches this cycle.
