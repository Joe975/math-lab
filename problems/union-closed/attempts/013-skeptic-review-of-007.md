# 013 — Skeptic review of 007 (averaged OR-control refuted): exact-arithmetic verification

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed
- **Type:** adversarial verification of `007-averaged-or-control.md` (default
  stance: refute). 007's central vulnerability, by its own confession, is that
  its two verifying engines share IEEE float arithmetic and its first
  "witnesses" were shared-underflow artifacts. The decisive check here is
  therefore **exact rational arithmetic end-to-end**: `fractions.Fraction`
  couplings, conditional tables, odds ratios and marginals, with every log₂
  replaced by a certified rational enclosure (atanh series, directed
  rounding, explicit tail bounds ≤ 1e-40). Every sign reported by the exact
  parts is a theorem about the stated rational inputs, not a float. All
  structural claims (§1 forcing, Lemmas G/M, Theorems A/B/C, both
  invariances) were additionally re-derived by hand.
- **Outcome in one line:** 007 holds in full. The 10-atom witness is REAL —
  the violation is certified in exact rational arithmetic at rational tilt
  parameters spanning λ ≈ 0.0144 to 5 (positive exactly where the
  first-order theorem forces it, violating above the onset, including at
  exactly rational λ = 1, 2, 3 and 4), the two invariances hold as exact
  rational identities, the aggregate survives with a certified positive
  sign, all four partial theorems re-derive, and — a strengthening 007 did
  not claim — the refutation survives even the most conjecture-friendly
  alternative bookkeeping convention. No corrections found; one cosmetic
  runtime discrepancy noted.
- **Tools:** `explore/uc_or_avg_skeptic.py` (parts A–G, W2; stdlib only,
  `fractions`-based on the critical path; no imports from `uc_or_avg.py`;
  own float census + own Sinkhorn — symmetric geometric-damped scaling, not
  007's row/column alternation; deterministic, fixed seeds; runtime ~3 min;
  checkpoints `data/oravgsk_part[ABCDEFWG].json`, log
  `data/oravgsk_run.log`). Command reproducing every number below:
  `python problems/union-closed/explore/uc_or_avg_skeptic.py 2>&1 | tee problems/union-closed/data/oravgsk_run.log`
- **Sources:** 005/006/007 (this repo); 008 consulted only for the
  consistency check on the Sinkhorn-smoothness step. No external fetches.

Notation as in 005/006/007: tilt coupling `π ∝ u(A)u(B)·t^{|A∩B|}` with
`t = 2^λ`; `M_i` the mass-weighted mean of `log₂ OR_i` over `N_i²` (007 §1's
well-posed definition); record threshold 0.38271. At λ = 3.5 the tilt
`t = 8√2` is irrational, so exact verification works at exact rational t
(the float-rationalization of `2^3.5`, `181/16`, and powers of 2 where
λ = log₂ t is itself a rational integer); the conjecture at tilt t is
`M_i ≥ log₂ t`, i.e. `Σ mass·log₂(OR/t) ≥ 0` with `OR/t` exactly rational —
a certifiable sign.

## Claims attacked

1. **The headline witness** (007 §4): the 10-atom μ on 2^[7], λ = 3.5, with
   `M_5 − λ = −0.122033`, all elementwise marginals ≤ 0.3178 < 0.38271,
   H(μ) = 2.72 bits — with the specific suspicion that it is a third
   shared-float artifact; plus the claims that the same μ violates for every
   λ ≳ 0.03 and is positive below (Theorem C consistency), the 2-digit tidy
   version, and the cell-floor/no-underflow claim.
2. **Well-posedness forcing** (§1): degeneracy dichotomy; direction-dependent
   ε-limits (spread 0.47) ruling out continuous surrogate assignment; the
   normalized mean forced by the i = n and product anchors.
3. **Theorem B and the cell-count argument** (≤4-atom supports safe), and
   the **two exact invariances** (marginal-1 stripping; fresh-prefix
   dilution) — 007's own named prime skeptic targets.
4. **The other partial positives:** potential-MTP₂ subclass; i ∈ {1, n}
   anchors; Theorem C (first-order-in-λ perfect square, modulo the Sinkhorn
   smoothness step); Lemmas G/M; the Frobenius form (Theorem A); the
   L-not-PSD / 2×2-minor-failure mechanism.
5. **The kill census's negative space** (444 instances, 0 violations) —
   would the battery have caught a violation had one existed in its
   families?
6. **The i-aggregated survival claim** (+1.84 on the witness; clamped
   aggregate search empty) — the proposed replacement gap.
7. **Consistency with 005/006**: that nothing there needs correction, and
   that Theorem B explains 006 S8.

## Refutations found

**None load-bearing.** Every kill attempt failed; the exact-arithmetic
attack in particular *strengthened* the record (see A below). Cosmetic
findings only:

- **R1 (cosmetic).** 007's Tools bullet says runtime ~2 min;
  `uc_or_avg.py`'s docstring says ~4–6 min. One of the two is wrong;
  nothing depends on it. (Recorded because 006's R3 showed docstring drift
  is what the next agent trips on.)
- **R2 (precision nuance, not an error).** The headline "the same μ
  violating for every λ ≳ 0.03" is float-supported in 007 (Sinkhorn refits
  of one μ). The exact-arithmetic form established here is: certified
  violations at exactly rational λ = 1, 2, 3 for measures within 9.2e-16
  (elementwise) of that μ — rationalized refit potentials, each exactly its
  own marginal's potential — plus 007's float profile reproduced
  independently. The distinction (same μ vs. an exact μ̃ within 1e-15 per λ)
  is stated here so nobody mistakes the certificate's scope; it does not
  weaken the refutation, which needs only one admissible μ per λ.

## Claims that survive (and what was done to break them)

### 1. The witness (claim 1) — CERTIFIED in exact rational arithmetic; the shared-float loophole is closed

The one attack 007 could not run on itself, run here (part A; all
enclosure widths < 1e-40):

- **At t = the float-rationalization of 2^3.5** (an exact rational,
  λ = 3.5 to 16 digits): `M_5 − log₂ t ∈ [−0.122033410223667 ± 1e-40]` —
  **certified violation**, agreeing with both float engines to 15 digits.
  Degeneracy dichotomy checked *exactly* at every coordinate (zero cell ⟺
  zero conditional margin ⟺ prefix outside N_i — no mixed pattern).
  Max elementwise marginal is an exact rational < 38271/100000 (printed in
  the log as an integer-fraction comparison). Defined history-mass fraction
  at i = 5: 0.5510878757 exactly. H(μ) = 2.7218 ✓.
- **At t = 181/16** (clean rational near 2^3.5): all seven per-i enclosures
  match 007's list (+1.4834, +1.9055, +2.2436, +2.7702, **−0.1220**,
  +3.4737, exactly 0 at i = n), in-regime, violation certified.
- **Exact-λ instances.** On the fixed-u family, certified signs at ten
  rational t: positive at t = 101/100, 103/100, 21/20, 11/10, 3/2, 2 and
  violating at t = 4, 8, 16 (in-regime) and 32 (marginal 0.397, out of
  regime). At t = 4, 16 the statement "M_5 < λ, λ = 2 (resp. 4) exactly, at
  a μ with all marginals < 0.38271" is a fully rational, assumption-free
  refutation of the restated gap — no irrational number appears anywhere in
  it. Sinkhorn-refit exact witnesses at λ = 1, 2, 3 (R2 above) cover the
  refit path; small-t positivity is exactly where the first-order theorem
  forces it (see claim 4).
- **Path subtlety worth recording:** on the *fixed-u* path the sign flips
  between λ = 1 (+0.0125) and λ = 2 (−0.0401), while on 007's *fixed-μ
  refit* path it flips at λ ≈ 0.03 (my independent refit engine reproduces
  007's entire part T(e) profile, small-λ positives included). Both are
  correct — the two paths move through different μ. Do not read the sweep
  table as contradicting T(e).
- **Tidy witness:** the exact-decimal 2-significant-digit measure violates
  by a certified −0.121887 at t = 181/16, max marginal 0.312983.
- **Convention robustness (new, strengthens 007).** Even under the most
  conjecture-friendly alternative bookkeeping — score every degenerate
  history at exactly λ instead of conditioning it out — the witness still
  violates: `M_5 − λ = −0.067173` (certified). The refutation therefore
  does not hinge on §1's exclusion convention at all; only a convention
  awarding degenerate histories *more* than λ (indefensible: the tilt makes
  no Plackett choice there) could rescue the conjecture on this μ.
- Underflow is structurally impossible in the exact parts; independently,
  the float cell-floor claim (3.0e-6) and atom-weight span reproduce.

**Verdict: CONFIRMED — the refutation is real and is now float-free.**

### 2. Well-posedness forcing (claim 2) — re-derived; CONFIRMED as the record scopes it

- *Dichotomy:* re-derived in two lines (cell (α,β) mass is
  `2^{λ|a∧b|}·2^{λαβ}·⟨F^a_α, W F^b_β⟩` with W strictly positive, so cell
  positivity factorizes through the two slices; a zero cell zeroes a full
  row or column of the table). Exact check at every coordinate of the
  witness (part A) and float check across the battery: zero exceptions.
- *ε-limits:* own engine reproduces 007's three directional limits
  (+1.4716 / +1.0399 / +1.0000, spread 0.4716; part W2). The inference —
  the per-history OR at a degenerate history has no continuous extension,
  so any surrogate value is arbitrary — is sound. One nuance the record's
  phrasing skates over: conditioning-out does not make `M_i` continuous in
  μ at support-degenerate points either (no convention does — as ε → 0 the
  repopulated history retains Θ(1) mass with a direction-dependent OR). The
  honest forcing is unidentifiability + the anchors, which is how 007's
  Outcome section itself calibrates §1 ("arguments of inevitability rather
  than formal theorems"). And by the convention-robustness computation
  above, the witness kills the conjecture under either resolution, so
  nothing downstream turns on this nuance.
- *Normalization forced by anchors:* re-derived — at i = n, `OR ≡ 2^λ`
  (005 Prop 1) gives `M_n = λ` on every support only for the normalized
  mean; the unnormalized sum gives `λ·(defined mass) < λ` whenever any
  prefix is deterministic. Anchors verified numerically: max |M_n − λ| =
  1.3e-15, min M_1 − λ = +9.4e-4 ≥ 0 over the battery (part F4).

**CONFIRMED** (with the continuity nuance recorded above).

### 3. Theorem B, cell-count, and the two invariances (claim 3) — re-derived by hand; invariances hold as exact rational identities

- *Cell-count:* if |N_i| ≥ 2 then each active prefix needs both response
  cells nonempty (≥ 2 atoms each), so ≤ 4 atoms force exactly two active
  prefixes with exactly one atom per cell; |N_i| ≥ 3 needs ≥ 6 atoms. With
  singleton cells the potential values cancel in the cross-ratio and
  `L_ab = log₂ OR − λ = λ⟨d_a, d_b⟩`, a Gram matrix, so L ⪰ 0 and Theorem A
  gives `M_i ≥ λ`. If |N_i| ≤ 1 the single diagonal history is ≥ λ by 005
  Prop 2. Re-derivation clean. Kill attempts: 500 random ≤4-atom instances
  (own seeds, n up to 7, λ up to 3.5): min excess −7e-15 (part F1). The
  crash-family cancellation was checked *exactly*: `OR = t^{3−n}` as a
  rational identity at random rational potentials and random rational t
  (part F2) — the singleton-cell mechanics are potential-free, as claimed.
  Boundary sharpness spot-checked: for the crash at i = 2,
  `|L_ab|² = (λm)² = e_a e_b` exactly (d-vectors are ∓1/±1), which is 007's
  "sits ON the PSD boundary" and why 006 S8 saw equality on the mirror.
- *Stripping:* appending an always-1 coordinate multiplies the kernel by
  exactly t; checked as an exact identity on the witness — every OR equal
  as a rational, every unnormalized history mass scaled by exactly t
  (part B). `M_i` exactly invariant.
- *Dilution:* with the three dilution atoms removed, N_5, the violating
  histories, and every unnormalized table cell and OR are **identical
  rationals** (not merely equal to 1.8e-15 — exactly equal; part B), while
  the max marginal moves 0.5076 → 0.3178. The invariance is exact for the
  stated reason: fresh-prefix atoms enter no cylinder of N_5, so they touch
  only the normalization, which cancels in the weighted mean.

**CONFIRMED.**

### 4. The remaining partial positives (claim 4)

- *Potential-MTP₂:* trivial given 005 Prop 5 (pointwise `OR ≥ 2^λ` ⟹ every
  average ≥ λ); Prop 5 itself carries 006's verification. Own probe: 30
  random ferromagnetic potentials, min pointwise `log₂OR − λ` = −3.6e-15
  (part F3). CONFIRMED.
- *i ∈ {1, n}:* i = 1 has the single diagonal history (∅,∅), Prop 2
  applies; i = n is Prop 1. Numerics above. CONFIRMED.
- *Theorem C:* structure re-derived — at λ = 0 the coupling is μ⊗μ; the
  O(λ) term of each `log⟨F_α, W G_β⟩` splits into (a,α)-only and
  (b,β)-only pieces (which the alternating cross-ratio sum annihilates,
  including any potential-drift contribution) plus the mixed term
  `λ⟨x̄_{aα}, ȳ_{bβ}⟩`, leaving `L_ab = λ⟨D_a, D_b⟩ + O(λ²)` and
  `M_i − λ = λ‖Σ ν(a)D_a‖² + O(λ²)`. The formula was then tested as a
  formula: Richardson-extrapolated finite-difference slopes vs the directly
  computed `‖Σν(a)D_a‖²` agree to 2.2e-5 relative on random μ (part F5),
  and on the witness μ the coefficient computes to **4.6530e-4**, matching
  007's 4.65e-4 (and the profile's small-λ curvature fit). The Sinkhorn
  λ-differentiability step remains the labeled SPECULATION; note it is
  *not* literally 008's Theorem P6′ (P6′ is smoothness in μ at fixed λ;
  Theorem C needs d/dλ at fixed μ — the same IFT machinery, but anyone
  citing 008 to discharge Theorem C's step must restate it). Along the
  fixed-u path no smoothness step exists at all (π is explicit in t), and
  there the exact small-t signs of part A certify first-order positivity on
  the witness family outright. CONFIRMED modulo the correctly-labeled step.
- *Lemmas G/M, Theorem A, mechanism:* Lemma G and Lemma M re-derived
  (`m_ab = 2^{λ|a∧b|}⟨F^a, W′F^b⟩`, Cauchy–Schwarz in the W′-inner product
  gives the `2^{−λ|aΔb|/2}` decay); Theorem A is definitional given L
  symmetric (from π = πᵀ) with the PSD-pairing corollary. Frobenius
  identity max err 8.9e-16, normalized mass-matrix min eigenvalue +2.8e-4
  (PSD) over the battery (part F4). The claimed failure mode reproduces on
  the witness itself: L at i = 5 has eigenvalues {−1.4785, +1.9990} (not
  PSD) and 2×2-minor ratio `|L_ab|/√(e_a e_b)` = **6.68**, matching 007's
  "≈ 6.7" (part F6). CONFIRMED.

### 5. The census's negative space (claim 5) — the battery would have seen it

Own 96-instance spot-battery (own Sinkhorn, own seeds: crash, mirror,
ε-mixed crash to 1e-5, δ₀⊕Bern, slices, three Sawin mixtures, 60 random
sparse): **0 violations**, and the closest-approach ordering matches 007
(Sawin/part-E mixtures nearest at +9.4e-4 on my sparser λ-grid — 007's
+1.5e-4 sits on its denser grid — then δ₀⊕Bern, slices, crash; mirror and
random-sparse at the −1e-16 equality floor). Then the witness genre was
planted into the *same pipeline*: witness, tidy witness, five own-seed 3%
perturbations, two light-weight variants — **9/9 flagged** at −0.12
(part E). So the pipeline detects the genre when present; the battery's
families genuinely do not contain it, exactly as 007's "why every earlier
test missed it" argues (≤4-atom families are safe by Theorem B; symmetric
mixtures align their D-vectors). CONFIRMED.

### 6. The i-aggregated form (claim 6) — survives this review's attacks too

- Exact: on the witness at t = 181/16, `Σ_{i<n} w_i(M_i − λ)/Σw_i ∈
  [+1.844669 ± 1e-40]` — certified positive (part C), matching +1.84.
- Own hunt (part G, fresh parametrization and seeds): hill-climbs seeded
  **at the witness potential** and at two-light-slice variants built to
  crash two coordinates at once ({1,5,6} light, {6,7} light), plus 12
  random sparse restarts at λ ∈ {2, 3.5, 5}: 21 endpoints, **0 in-regime
  violations**, lowest endpoint +0.472. The witness-seeded climbs could not
  trade its i = 5 deficit against the other coordinates' surplus. The
  aggregate remains unrefuted — and remains only EVIDENCE, as 007 labels
  it. CONFIRMED (as a survival claim, not as a truth claim).

### 7. Consistency with 005/006 (claim 7)

- 007's definition provably reduces to what 005/006 computed (007 part C:
  8.9e-16 vs 006 S8, 2.8e-14 vs 005 part E); my independent engine
  reproduces both S8 directional numbers and the part-E-mixture excesses
  inside the battery. Nothing in 005/006 is contradicted: 005 lead 1 was
  labeled SPECULATION and is now settled negatively, which is the label
  working as intended; 006 S8's "survives its first falsifiable test" was
  a true statement about two 4-atom families — and Theorem B shows those
  families *could not have* failed, which upgrades S8 from evidence to a
  corollary and explains the sharpness (mirror at exact equality) 006
  observed. 006's bookkeeping note is answered by §1. CONFIRMED.

## Verdict

| # | Claim | Verdict |
|---|-------|---------|
| 1 | 10-atom witness: in-regime violation `M_5 = λ − 0.122` | **CONFIRMED — certified in exact rational arithmetic** (enclosure −0.122033410223667 ± 1e-40 at rationalized 2^3.5; exact marginals < 38271/100000; violations certified at exactly rational λ = 2, 3, 4 on the family and λ = 1, 2, 3 on rationalized refits; shared-float loophole closed) |
| 1b | Same μ violates ∀λ ≳ 0.03, positive below | **CONFIRMED** (own refit engine reproduces the full profile; exact certificates at λ = 1, 2, 3 for μ̃ within 9.2e-16; small-λ positivity certified exactly on the fixed-u path — see R2 for the precise certificate scope) |
| 2 | §1 well-posedness forcing | **CONFIRMED** as scoped by 007 (dichotomy re-derived + exact-checked; ε-limits reproduced; anchors force normalization; nuance: no convention is continuous in μ — and the new surrogate-λ computation shows the refutation is convention-robust anyway, M_5 − λ = −0.067 certified) |
| 3 | Theorem B + cell-count; stripping & dilution invariances | **CONFIRMED** (hand re-derivation; 500 random ≤4-atom instances; exact potential-free crash identity; both invariances hold as **exact rational identities**) |
| 4 | MTP₂ subclass; i ∈ {1,n}; Theorem C; Lemmas G/M/A; mechanism | **CONFIRMED** (Theorem C modulo its correctly-labeled Sinkhorn-λ-smoothness step — which is analogous to but not identical with 008's P6′; coefficient formula verified to 2.2e-5, witness coefficient 4.6530e-4 ✓; witness L non-PSD with minor ratio 6.68 ✓) |
| 5 | 444-instance census clean; genre invisible to it | **CONFIRMED** (own 96-instance battery clean with matching ordering; planted genre detected 9/9) |
| 6 | i-aggregated form survives (+1.84; searches empty) | **CONFIRMED** (+1.844669 certified positive exactly; own witness-seeded and two-coordinate attacks found nothing, lowest +0.472) |
| 7 | Nothing in 005/006 needs correction; Theorem B explains S8 | **CONFIRMED** (S8 upgraded from evidence to corollary of Theorem B, including the observed equality boundary) |

**Net assessment:** 007's headline — the restated Gap 1
(`averaged-odds-ratio-control`, `M_i ≥ λ` for all μ, i) is **REFUTED
in-regime** — should be treated as VERIFIED, and is now certified by exact
rational arithmetic, closing the one attack class (shared IEEE pathology)
that cross-engine agreement could not exclude. The four partial positives
stand as stated (Theorem C with its labeled gap). The record's
self-calibration was accurate: the part it trusted most survived the
strongest attack available, and its named skeptic targets (Theorem B,
§1 forcing) survived hand re-derivation. The route's gap list correctly
advances to `aggregated-or-control` / `margin-modulated-or-control`.

## Residual risk

- **Shared-definition risk.** All engines (005's, 006's, 007's, and this
  one) formalize "conditional odds ratio of the tilt coupling given
  histories" the same way from 003's definitions. Exact arithmetic removes
  numerical error, not a conceptual error in the shared census definition
  itself. Mitigated by 006 S1 (census vs the independently-derived
  slice/inner-product formula) and by the hand re-derivations here, but a
  from-scratch third formalization has still never been written.
- **Theorem C's smoothness step** (∂/∂λ of the Sinkhorn potential at fixed
  μ) remains open here; 008 claims the μ-direction analogue with
  quantitative constants but that statement does not literally discharge
  this one. All numerics are consistent with it.
- **The aggregated form** has now survived two independent search
  campaigns and one exact certificate — it is still only EVIDENCE on
  n ≤ 8-scale instances. My witness-seeded attacks are hill-climbs, not
  exhaustive; the two-response-coordinate construction space (007 lead 1)
  is larger than what 20+21 restarts explore.
- **Certified enclosures** rest on the correctness of the atanh-series
  tail bound implemented in `log2_enclosure` (~15 lines, elementary, but
  itself not independently reviewed).

## References

- Reviewed record: `problems/union-closed/attempts/007-averaged-or-control.md`;
  its tool `explore/uc_or_avg.py` and checkpoints
  `data/or_avg_part[WCKLSTU].json`, `data/or_avg_run.log`.
- Context: `attempts/005-odds-ratio-control-refuted.md` (Props 1–6, slice
  reduction, crash family); `attempts/006-skeptic-review-of-005.md` (S8,
  bookkeeping note, λ>0 equality correction);
  `attempts/008-perturbative-assembly.md` (Theorem P6′, consulted for the
  smoothness-step consistency check only); `attempts/004-*` (review shape).
- This review's tool/data: `explore/uc_or_avg_skeptic.py`;
  `data/oravgsk_part[ABCDEFWG].json`; `data/oravgsk_run.log`.
- No external sources; the atanh series for ln and its geometric tail bound
  are standard calculus, re-derived inline in the tool's docstrings.
