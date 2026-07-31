# 011 — Skeptic review of 009 (mutual-information tax): adversarial verification

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed
- **Type:** adversarial verification of `009-mutual-information-tax.md`
  (default stance: refute). Both mini-lemmas re-derived by hand from the
  003/004 definitions; every headline number re-computed with an independent
  implementation (`explore/uc_mitax_skeptic.py`, no imports from
  `uc_mitax.py`): natural-log internals vs 009's log2; exchangeable couplings
  by a 4-cell one-coordinate recursion plus forward state filtering vs 009's
  closed-form binomial sums and prefix×future mass products; asymmetric
  alternating Sinkhorn with symmetry *checked* vs 009's damped symmetric
  iteration; taxes cross-checked by direct conditional-MI computation and the
  U-side chain rule, not only via 009's entropy-difference identities.
  Additionally, 009's admitted untested cell (tilt-recipe CR on Sawin
  gadgets) and both admitted grid edges were run here, plus 009's lead 1
  (crash-window scaling law).
- **Outcome in one line:** 009 holds. Both mini-lemmas survive independent
  re-derivation (one hidden hypothesis in the part-D block bound located and
  verified to hold on every instance used); all headline numbers reproduce
  (worst case 1.4e-10, at the Sinkhorn tolerance); the λ-window narrowing is
  sharpened to a clean law λ_max·(n−3) → ≈ 4.85; the ST ≈ 0.40·log₂n fit
  survives a 3× range extension and two new p values; the untested
  tilt-on-Sawin cell closes POSITIVE at n ≤ 300. Corrections found are
  reporting-level: sup_λ CR values at p ≈ 0.42 slices are grid-edge values
  (true suprema are larger, at λ ≈ 4–5); "increasing margin in n" mixes two
  p values; the Outcome's "p ∈ {0.42, 0.3823}" overstates which slice ratios
  were actually run; part D's "exactly" has a measured ~4e-6 absolute
  precision floor from lgamma.
- **Tools:** `explore/uc_mitax_skeptic.py` (parts SK-A..SK-F; stdlib only;
  deterministic, no RNG anywhere; runtime ≈ 4 min; checkpoints
  `data/mitaxsk_part[A-F].json`, log `data/mitaxsk_run.log`). Command
  reproducing every 011 number below:
  `python3 problems/union-closed/explore/uc_mitax_skeptic.py`
- **Sources:** 002/003/004/005/006/009 (this repo). No external fetches; no
  new outside facts used.

Notation as in 009: histories `a = A_{<i}`, `b = B_{<i}`, coupling π with
both marginals μ; `z_i(a,b) = Pr[A_i=0, B_i=0|a,b]`,
`x̃_i(a,b) = Pr[A_i=0|a,b]`; `CR = Σᵢ E[h(z_i)] − H(μ)`;
`T_A = Σᵢ I(A_i; B_{<i}|A_{<i})`; `ST = Gain − CR`. Bits throughout.

## Claims attacked

1. The four one-line facts (i)–(iv) of 009 §1, the well-posedness of
   (TAX at p) as a formalization of Gap 2, and the implication chain
   (TAX at p) ⟹ (S-coup at p) ⟹ Frankl at p.
2. **Mini-lemma SL** (pure slice tilt is tax-free, T_A = T_B = 0
   identically) — the most load-bearing structural claim.
3. **Mini-lemma HM** (half-mixing coupling has ST = 0 identically) and the
   part-D block-coupling bound `CR_block ≥ nΣP_k h(m_k) − H(μ)`.
4. Part A evaluator validation (closed forms; 003_partB2 tie-in).
5. Part B crash-family numbers, the claimed reproduction of 005's crash-OR
   identity, the λ-window narrowing claim, and both λ grid edges.
6. Part C half-mixing numbers and the transfer of 004's closed-form gain
   to CR.
7. Part D: H(μ) of the 002 certificate gadgets re-derived from 002's
   definitions, the CR lower bound, and the "each within H(P) of 003's
   part-D gain values" sandwich.
8. Part E slice numbers, the "sup_λ CR > 0 at every (n, p) tested" and
   "increasing margin in n" claims, the p grid edge, and the
   ST ≈ 0.40·log₂n four-point fit.
9. The untested cell 009 admits (tilt-recipe CR on Sawin gadgets at large
   n) — run here, since a negative sign there is the nearest available kill.

## Refutations found

Nothing load-bearing. Four genuine but small defects, all reporting-level:

### R1. The slice "sup_λ CR" values at p ≈ 0.42 are grid-edge values, and
### "increasing margin in n" mixes two different p
009 §4(c) reports sup_λ CR = +0.56 (n=24) to +1.02/+1.19 (n=80/60) with
"the needed λ sits at the top of the grid (λ ≈ 2–3)". Extending the grid
(SK-E2): CR keeps rising past λ = 3 and has its interior maximum near
λ ≈ 4–5 — at (n, w₀) = (80, 34) the true supremum is ≈ +1.4655 at λ = 5.0
(vs the reported +1.0199 at λ = 3), falling again by λ = 6 (+1.29) and
λ = 8 (+0.79). Understatement only — the sign conclusion is unchanged and
strengthened — but the reported values are not suprema, exactly the class
of slip 004 flagged in 003 (its R3 note). Second, "by an increasing margin
in n (+0.56 at n=24 to +1.02..1.19 at n=80/60)" is confounded: the n
sequence 24, 40, 60, 80 alternates w₀/n between 0.4167 and 0.425, and the
apparent 60 → 80 *decrease* (1.19 → 1.02) is the p-alternation, not an n
effect. At genuinely fixed p the margin is cleanly increasing (SK-E4,
λ = 3: p = 5/12: +0.56, +0.99, +1.38, +1.76, +2.14, +2.87, +3.95 at
n = 24..240; p = 0.425: +0.61 → +2.10 at n = 40..200).

### R2. The Outcome's "p ∈ {0.42, 0.3823}" overstates what was run
w₀ must be an integer: the slice ratios actually tested are
{0.4167, 0.425} and {0.3833, 0.3875}. In particular p = 0.3823 (the Liu
record threshold, and the label 009's Outcome uses in set notation) was
never itself tested. §4(c)'s own "p ≈" phrasing is accurate; the Outcome
bullet and the range field are not. Cosmetic, but scope lines are what the
index filters on.

### R3. Part D's "H(μ) computed exactly" has a measured ~4·10⁻⁶ floor
"Exactly (profile arithmetic)" and "rigorous given float arithmetic"
overstate the precision: profile arithmetic at n = 60000 accumulates
lgamma error. Measured directly (SK-D investigation): computing
H(Bern(1/2)^60000) — exactly 60000 bits — by the same profile route gives
an error of +4.35e-6, and the part-D sandwich (below) shows an excess of
the same size and sign. Both 009's code and mine share this floor.
Irrelevant at the +2418-bit scale of the conclusion (12 significant
figures of headroom), but "exactly" should read "to ≈ 1e-5 absolute".

### R4. The part-D block lemma is applied with an unstated legality
### hypothesis (which holds on every instance used)
009's small lemma takes "per coordinate, union marginal m_k" as given by
003's part-D coupling, with m_k = 1/2 whenever p_k ≤ 1/2. A
conditionally-iid per-coordinate coupling with margins p can only reach
union marginals in [p, 2p−p²], so m_k = 1/2 requires
`p_k ≥ 1 − 1/√2 ≈ 0.29289`. Neither 003 nor 009 states this. Checked on
all three certificate gadgets and all SK-F instances (SK-D `mk_legal`):
the only components with p_k ≤ 1/2 are the k = 0 ones, p₀ = ū ≥ 0.3823 >
0.29289 — legal everywhere used. A future reuse of the bound on mixtures
with a component marginal below 0.293 would silently overstate CR_lb.

## Claims that survive (and what was done to break them)

### Facts (i)–(iv), (TAX at p), and the implication chain (claim 1) — re-derived by hand
- (i) `H_π(U) = Σᵢ H(U_i|U_{<i}) ≥ Σᵢ H(U_i|A_{<i},B_{<i}) = Σᵢ E[h(z_i)]`:
  chain rule plus conditioning on the finer σ-field (U_{<i} is a function
  of (A_{<i}, B_{<i})). Hence CR ≤ Gain and
  `ST = Σᵢ I(U_i; (A_{<i},B_{<i}) | U_{<i}) ≥ 0`. ✓
- (ii) `Σᵢ E[h(x̃_i)] = Σᵢ H(A_i|A_{<i},B_{<i}) = Σᵢ [H(A_i|A_{<i}) −
  I(A_i;B_{<i}|A_{<i})] = H(μ) − T_A`. ✓ Numerically cross-checked by
  computing T_A *directly* as summed conditional MI from the joint tables:
  identity error ≤ 1.7e-12 on every SK-A instance.
- (iii) `I(A;B) = Σᵢ I(A_i;B|A_{<i}) ≥ Σᵢ I(A_i;B_{<i}|A_{<i}) = T_A`
  (B_{<i} is a function of B). ✓ Numeric check on the SK-A4 mixture
  coupling: T_A = 0.000617 ≤ I(A;B) = 0.074481.
- Implication chain: CR > 0 ⟹ Gain > 0 (by (i)) ⟹ (S-coup at p) for any
  class containing the recipe coupling; with the licensing lemma
  (re-derived in 004) this gives Frankl at p, |F| ≥ 2 edge cases included. ✓
- **Definitional audit, two nuances recorded, neither a defect.** (a) As a
  universally quantified statement, (TAX at p) is well-formed only relative
  to a *total* recipe μ ↦ π(μ); 009's recipe is per-genre with Gap 3
  (totality) open, which 009 says explicitly. Since the λ-swept Sinkhorn
  tilt exists for every μ (strictly positive kernel), a total default
  recipe exists and the implication to Frankl is sound; but (TAX at p) is
  really a family of statements indexed by the recipe, and its truth value
  can differ by genre assignment — exactly the L3/Gap-3 coupling 009
  records. (b) CR conditions the U-chain on the *full* (a,b) history — the
  finest choice, hence by (i)-type monotonicity the SMALLEST (most-taxed)
  member of the family of interleaved assembly values between CR and Gain
  (e.g. conditioning on (A_{<i}, U_{<i}) gives a value ≥ CR). So (TAX at p)
  is the most conservative chain-rule certificate, not the sharpest
  available; Gap 2 could in principle be provable for an intermediate
  assembly while (TAX) fails. 009's L7 ("interleave the chain rule more
  cleverly") already names this; recorded here as a definitional fact, not
  a correction.

### Mini-lemma SL (claim 2) — re-derived by hand; the proof is correct and proves more
Setup: π ∝ 2^{λ|A∩B|} on slice(w₀)², both marginals uniform on the slice
(the marginal weight of A is `Σ_B ρ^{|A∩B|}`, slice-constant by
exchangeability — 003's observation, re-checked). Fix history (a,b) at
coordinate i and any completion A′ ⊇ a in the slice with future part x′,
|x′| = w₀ − w_a fixed. Then

    Σ_{B ext b} ρ^{|A′∩B|} = ρ^{|a∧b|} · Σ_y ρ^{⟨x′,y⟩}
      = ρ^{|a∧b|} · Σ_j C(w₀−w_a, j) C(m−(w₀−w_a), (w₀−w_b)−j) ρ^j,

where y runs over the future completions of b (fixed size w₀ − w_b in the
remaining m coordinates): the sum depends on x′ only through |x′|, which
is the same for every completion of a. So P(A = A′ | a, b) is the same for
all A′ extending a — the conditional law of A given (a,b) is UNIFORM over
slice-completions of a, independent of b. This yields A_i ⊥ B_{<i} | A_{<i}
(hence T_A = 0, and T_B = 0 by symmetry) and is strictly stronger: the
whole A-future is conditionally independent of the B-history. Hypothesis
audit: needs only finite exchangeability of the future block and a single
slice as support; λ may be any real; no positivity subtleties (all pair
weights > 0). Numerically: |T_A| ≤ 1e-10 across every slice run here,
including n = 240 (float accumulation, not structure). **VERIFIED.**

### Mini-lemma HM (claim 3a) — re-derived by hand; direct σ-field check added
With cells (∅,∅), (∅,s), (s,∅), (s,s) of masses P₁/2, P₁/2, P₁/2,
1−1.5P₁ per drawn s ~ Bern(p_h)^n (this needs P₁ ≤ 1/2 — satisfied by
both 009 points and both new ones; 004's "with probability 2P₁" already
implies it): U = ∅ on the first cell, U = s otherwise. If U_{<i} = u ≠ 0
then U = s and s_{<i} = u are forced; every (a,b) refining u pins down the
cell (e.g. (0,u) ⟹ (∅,s)) but the cell choice is independent of s, so
conditionally on any such (a,b), U_i = s_i ~ Bern(p_h) — the same law as
given U_{<i} = u alone. If u = 0 the only refining history is (a,b) =
(0,0): the two σ-fields have the *same atom*, so the conditional laws
trivially agree. Hence ST = Σᵢ I(U_i; (a,b)|U_{<i}) = 0. ∎ Numerically
(SK-C): ST = 0 to ≤ 3.2e-11 over 16 runs — the two 009 parameter points
plus (p_h, P₁) = (0.75, 0.40) and (0.51, 0.50) which 009 never ran, n up
to 16 — and a DIRECT check of the proof's mechanism: max over histories
(a,b) with a∨b ≠ 0 of |P(U_i=1|a,b) − p_h| ≤ 8.2e-14. Also confirmed:
T_A ≈ 0.108 / 0.126, constant in n, as 009 reports. One nuance: 004's
closed-form gain is exact only up to the exponentially small overlap of
δ_∅ with Bern(p_h)^n (measured: |gain − closed| = 4.8e-3 at n=6 falling
to 1.9e-6 at n=16, ≈ (1−p_h)^n) — consistent with 009's "matches to 4
decimals by n = 12", which is accurate as stated. **VERIFIED.**

### The part-D block bound (claim 3b) — re-derived; hypothesis R4 added
Given k, the block coupling draws (A_i, B_i) per coordinate via common
randomness independent of both histories, with union marginal m_k; so
`E[h(z_i)] = H(U_i|a,b) ≥ H(U_i|a,b,k) = Σ_k P_k h(m_k)` (conditioning
reduces entropy; the conditional is history-free given k), summing to
`CR_block ≥ nΣ_k P_k h(m_k) − H(μ)`. Correct — subject to the m_k
achievability hypothesis of R4, verified on every instance used.
**CONFIRMED** (with R4's hypothesis now stated).

### Part A validation (claim 4) — reproduced with independent code
Product+Plackett closed form matched to < 1e-9 at ρ ∈ {1, 1.5, 3} (my
Plackett solver re-derived from the cross-ratio quadratic and verified by
recomputing the cross-ratio, ≤ 5.3e-15); pure-slice n=8 gain matches the
checked-in `003_partB2.json` to 2.4e-13 and 009's CR/ST to 2.1e-14;
recursion engine vs atom engine ≤ 2.2e-13 (slice) and ≤ 1.1e-12
(Sinkhorn-fitted two-product mixture, where my profile-Sinkhorn and
atom-Sinkhorn must and do converge to the same coupling; scaling-vector
symmetry deviation ≤ 1.1e-13, *checked* not assumed). **CONFIRMED.**

### Part B crash family (claim 5) — reproduced, extended past both edges, and the window law found
- All shared grid points reproduce to ≤ 1.4e-10 (Sinkhorn-tolerance scale,
  not structure). The n-independence of the λ = 0 point re-proved by hand
  via the 2-coordinate reduction and closed form: CR(0) = +0.111394,
  Gain(0) = +0.939005 at every n ∈ {5, 8, 11, 14, 17} — 009's +0.1114 /
  +0.9390 ✓.
- 005's crash-OR identity `log₂OR = λ(3−n)` reproduced from MY fitted
  couplings at every (n, λ ≠ 0) pair on the grid, worst error 1.4e-14 —
  009's claimed evaluator cross-check is real.
- Grid edges: CR stays negative out to λ = −2 (≈ −0.19) and λ = 5
  (≈ −0.013) at every n — no sign resurrection past 009's [−1, 3] window.
- New n: at n = 14 and 17 the supremum over λ is attained AT λ = 0 and
  equals the iid value +0.1114 — the candidate still survives, and the
  tilt component of the recipe has become useless on this family (λ-sweep
  must include 0, as 009's L3 says).
- **The window law (009 lead 1, resolved positively):** bisecting the
  upper zero λ_max(n) of CR: 2.52593, 0.97336, 0.60650, 0.44080, 0.34627,
  0.28515 at n = 5, 8, 11, 14, 17, 20, giving λ_max·(n−3) = 5.0519,
  4.8668, 4.8520, 4.8488, 4.8478, 4.8475 — a clean convergent scaling law
  λ_max ≈ 4.847/(n−3), matching the crash exponent's 1/(n−3) shape. The
  lower edge also narrows: λ_min·(n−3) drifts −1.44 → −2.12 (n = 5 → 20),
  slower to converge. EVIDENCE at these n; the law is now a concrete
  target for proof. **CONFIRMED and sharpened.**

### Part C half-mixing numbers (claim 6) — see HM above. **CONFIRMED.**

### Part D Sawin certificates (claim 7) — re-derived from 002's definitions
Own truncation (tail < 1e-30 vs 009's 1e-28), own natural-log lse: H(μ) =
1927.6997 / 19239.5926 / 57578.9854 — matches 002's table and 009's values
to ≤ 1.7e-8; max marginal Σ P_k p_k = 0.402269 / 0.390799 / 0.382536 ✓;
CR_lb = +66.5794 / +741.3309 / +2418.5522 ✓ (all diffs ≤ 1.7e-8, pure
truncation choice). The "within H(P)" sandwich: computing H(U_block)
independently (law(U) = Σ_k P_k Bern(m_k)^n, profile arithmetic),
gain_block matches `003_partD.json` to ≤ 1.5e-7 and
0 ≤ gain_block − CR_lb ≤ H(P) holds at n = 2000 and 20000; at n = 60000
the difference *exceeds* H(P) by 4.1e-6 — investigated and traced to the
lgamma floor of R3 (the identity test H(Bern(1/2)^n) = n errs by the same
+4.35e-6), so the apparent violation is numerical, not real. The
smoothing-sensitivity conclusion (Θ(n)-positive CR on the no-go's own
certificates) stands with ~12 significant figures of margin.
**CONFIRMED** (modulo R3's wording fix).

### Part E slices (claim 8) — reproduced exactly; scaling law extended 3× in n
- All 42 (n, w₀, λ) points reproduce to ≤ 5.2e-12 from the independent
  recursion engine. sup_λ CR > 0 at every (n, p) tested ✓ (and understated
  per R1).
- p edge: at n = 80 the certificate does not die past 009's grid — sup_λ
  CR = +1.05 (w₀=35, p=0.4375), +0.69 (w₀=36, p=0.45), +0.088 (w₀=38,
  p=0.475), with the needed λ growing to 5–6. No sign change anywhere near
  the tested boundary; slices are CR-certified even above the 0.4315
  2-block recipe ceiling (different genre — no contradiction).
- **ST ≈ 0.40·log₂n survives every extension tried** (SK-E4): at λ = 3,
  least-squares ST = a + b·log₂n gives b = 0.4048 (p = 5/12, n = 24..240,
  7 points, rms 0.0010), b = 0.4048 (p = 0.425, n = 40..200, rms 0.0002),
  b = 0.4090 (p = 0.3833, n = 60..240, rms 0.0001) — and the intercepts
  are ≈ 0 (−0.03..−0.06), so 009's "0.40·log₂n" is the right two-parameter
  summary, not a small-range accident. The alternative power-law fit
  ST = c·n^α is 13–40× worse in rms everywhere (α ≈ 0.22–0.25 with
  visible curvature). Successive slopes at p = 5/12 are 0.402, 0.405,
  0.406, 0.406, 0.406, 0.406 — no drift up to n = 240, in particular no
  drift toward 1/2. Caveat kept: the slope is mildly λ-dependent (b =
  0.3802 at λ = 1.5), so "0.40" is a λ = 3 constant, not universal; and a
  fixed λ = 1.5 has CR turning negative at n ≥ 72 (p = 5/12) — on slices
  too, the recipe needs the λ-sweep (here toward larger λ, opposite to
  the crash family's λ → 0: the two genres pull the sweep in opposite
  directions, which is Gap 3's totality problem in miniature).
  **CONFIRMED as EVIDENCE, range extended to n ≤ 240.**

### The untested cell (claim 9) — run here; it closes POSITIVE
Tilt-recipe CR on Sawin gadgets via profile Sinkhorn + the recursion
engine (O(n³); n ≤ 300 reachable, the 002 certificate sizes 2000–60000
are not, so this is the large-n frontier available to any exact method
here). Tie-in to verified numbers: my tilt gains reproduce 003's part (a)
table — +2.6120 at (n=200, ū=0.40, θ=0.05, λ=2.5) vs 003's +2.61, and
+10.3574 at (n=300, ū=0.3823, θ=0.02, λ=2.5) vs 003's +10.36 (004
verified +10.357). Results (sup over the λ grids run):

    shape (ū, θ)      n     max marg   sup_λ CR_tilt (λ)     block CR_lb
    (0.40, 0.05)      60    0.4124      +0.750  (2.5)          +1.402
    (0.40, 0.05)     100    0.4124      +1.250  (2.5)          +2.473
    (0.40, 0.05)     200    0.4124      +4.055  (3.5)          +5.220
    (0.3823, 0.02)   100    0.3871      +3.432  (2.5)          +3.818
    (0.3823, 0.02)   200    0.3871      +6.871  (2.5)          +7.765
    (0.3823, 0.02)   300    0.3871     +11.083  (3.5)         +11.719

sup_λ CR_tilt > 0 at every instance, growing ≈ linearly in n at fixed
shape, trailing the block bound by O(1). The assembly losses on this genre
are strikingly small — T_A ≤ 0.092 and ST ≤ 0.22 bits across the whole
table (compare gains of 2–11 bits): chain-rule certification is nearly
lossless on Sawin mixtures under the tilt, not just under the block
coupling. The λ ≲ 0.25 end is CR-negative at both shapes, so the sweep
must reach moderate λ (0.5–1) here while the crash family forces λ → 0 —
the recipe survives only because it sweeps. 009's leads-5 worry (tilt-CR
negative at large n forcing a genre boundary in Gap 3) is NOT realized at
n ≤ 300; whether it appears at the 002-certificate scale n ≥ 2000 remains
open (out of reach of exact evaluation here). **EVIDENCE, n ≤ 300: the
candidate survives its last untested genre-recipe cell.**

## Verdict

| # | Claim | Verdict |
|---|-------|---------|
| 1 | Facts (i)–(iv); (TAX at p) well-posed; ⟹ (S-coup) ⟹ Frankl | **CONFIRMED** (hand re-derivation; direct-MI and I(A;B) numeric checks; two definitional nuances recorded — recipe-totality dependence, CR = most conservative assembly member) |
| 2 | Mini-lemma SL (slice tilt tax-free) | **VERIFIED** (independent hand proof — conditional law is fully uniform, stronger than stated; \|T_A\| ≤ 1e-10 through n = 240) |
| 3 | Mini-lemma HM (half-mixing ST = 0) | **VERIFIED** (independent hand proof; direct σ-field check ≤ 8.2e-14; two new parameter points; P₁ ≤ 1/2 hypothesis made explicit) |
| 4 | Block bound CR_block ≥ nΣP_k h(m_k) − H(μ) | **CONFIRMED**, with unstated legality hypothesis p_k ≥ 1−1/√2 found and verified on all instances used (R4) |
| 5 | Part A evaluator validation | **CONFIRMED** (independent engines agree ≤ 1.1e-12) |
| 6 | Part B crash: numbers, OR-identity reproduction, λ-window narrowing | **CONFIRMED** (≤ 1.4e-10; OR ≤ 1.4e-14) **and sharpened**: λ_max(n)·(n−3) → 4.8475 (n = 20); no sign change past either grid edge; sup at λ = 0 exactly for n ≥ 14 |
| 7 | Part C half-mixing numbers; closed-form transfer | **CONFIRMED** (ST ≤ 3.2e-11 over 16 runs) |
| 8a | Part D: H(μ) vs 002, CR_lb, smoothing-sensitivity | **CONFIRMED** (own truncation, diffs ≤ 1.7e-8); "exactly" **CORRECTED** to ~1e-5 absolute floor (R3, measured) |
| 8b | Part E: 42 slice points; sup_λ CR > 0; ST ≈ 0.40·log₂n | numbers **CONFIRMED** (≤ 5.2e-12); scaling **CONFIRMED as EVIDENCE** with 3× range and two new p (b = 0.405 ± 0.004 at λ = 3, log-fit beats power-law 13–40×); sup_λ prose **CORRECTED** (R1 grid-edge values, mixed-p "increasing margin"); Outcome scope **CORRECTED** (R2: actual p ∈ {0.4167, 0.425, 0.3833, 0.3875}) |
| 9 | Untested tilt-on-Sawin cell (run here) | closes **POSITIVE** at n ≤ 300; assembly losses ≤ 0.31 bits on this genre; λ-sweep necessity confirmed from both directions |

**Net assessment:** 009's EVIDENCE verdict stands and is strengthened. The
probe's two structural identities are now independently proved (SL, HM);
the numerics reproduce from genuinely different code at or below Sinkhorn
tolerance everywhere; the one cell 009 could not afford closes in the
candidate's favor; and its ST-scaling fit turns out to be more robust than
its own record claims. The corrections are all in prose precision, not in
substance — grid-edge suprema reported as suprema (R1) being the only one
with any teeth, and it errs in the conservative direction. The gap list
for the route is unchanged: (TAX at p) remains unproved for every p, and
the hard open items are exactly 009's L3/L4/L6/L7.

## Residual risk

- **Shared-reduction risk.** My recursion engine and 009's DP both rest on
  the same mathematical reduction: for exchangeable-potential tilt
  couplings the prefix state (w_a, w_b) is sufficient for the conditional
  future. I re-derived this (the prefix tilt factor ρ^{|a∧b|} cancels in
  every conditional; the future weight depends on the prefix only through
  its two weights) and both engines were validated against brute-force
  atom enumeration — but only at n = 8. An error in the *reduction itself*
  would evade both codes at large n; the ≤ 5.2e-12 agreement is evidence
  of implementation independence, not of the reduction.
- **ST slope asymptotics.** The empirical slope 0.405 is constant to
  n = 240 with tiny residuals, but nothing here proves ST = Θ(log n); a
  slow drift (e.g. toward 0.5·log₂n, or an eventual n^α takeover) beyond
  n = 240 is not excluded. 009's lead 2 (prove it via local CLT) is still
  the decisive step for L7.
- **The tilt-on-Sawin cell is closed only to n = 300.** The O(n³) exact
  engine cannot reach the 002-certificate sizes (n ≥ 2000). The trend
  (linear-in-n CR, shrinking losses) points the right way, and the block
  recipe covers the genre rigorously at all n, but "tilt-CR > 0 on Sawin
  gadgets" is EVIDENCE for n ≤ 300 only.
- **λ grids remain finite.** Part F's suprema at (200, 0.40) and (300,
  0.3823) sit at my grid top λ = 3.5; the sign conclusion is safe (CR is
  already large and positive) but the reported suprema are again
  grid-edge values, same caveat as R1 in the other direction.
- **Sinkhorn tolerance.** All coupling-level agreement is bounded by the
  fitting tolerances (~1e-13 residuals propagating to ~1e-10 in CR on the
  crash family); nothing structural, but sub-1e-10 claims about fitted
  couplings should not be made from either code.

## References

- 009-mutual-information-tax.md (the record under review);
  003-dependent-couplings.md; 004-skeptic-review-of-003.md;
  005-odds-ratio-control-refuted.md; 006-skeptic-review-of-005.md;
  002-weighted-kl-ladder.md (this repo).
- `explore/uc_mitax_skeptic.py`; checkpoints `data/mitaxsk_part[A-F].json`,
  log `data/mitaxsk_run.log`; reviewed artifacts `explore/uc_mitax.py`,
  `data/mitax_part[A-E].json`, `data/mitax_full_run.log`; prior data used
  as cross-checks: `data/003_partB2.json`, `data/003_partD.json`, 002's
  certificate table (H values 1927.70 / 19239.59 / 57578.99).
- No external sources consulted this cycle; all outside facts enter only
  through 001–004 as already transcribed there.
