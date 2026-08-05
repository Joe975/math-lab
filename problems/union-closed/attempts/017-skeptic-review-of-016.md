# 017 — Skeptic review of 016 (aggregated control killed at large n): independent exact recomputation

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed
- **Type:** adversarial verification of `016-aggregated-control-probe-large-n.md`
  (default stance: refute the refutation — a wrong kill would close a live
  proof branch incorrectly). 016's certificates reuse 013's exact engine and
  013's log₂ enclosures; the decisive check here is therefore a **from-scratch
  third implementation of the aggregate**, derived from 007 §1's prose (not
  transcribed from any engine), in scaled-integer exact arithmetic with a
  **different certified log₂ enclosure** (interval squaring with directed
  dyadic rounding + the exact series ln 2 = Σ 1/(k·2^k), not 013's atanh
  series). Per the verification standard (2026-07-31), exact rational
  arithmetic IS the skeptic pass for claims decidable over ℚ; this review
  additionally removes the shared-enclosure-code and shared-census-code risks
  013's Residual-risk section named.
- **Outcome in one line:** **016 holds — the kill stands.** All three kill
  certificates (n = 96, 128, 160 at t = 4) re-derive from scratch with
  matching digits and certified negative sign; the positive controls
  (+1.844669005 witness anchor, P3-best n = 32, the raw-weight ladder at
  n = 96) re-certify positive; the construction, marginals, dichotomy, and
  convention-robustness claims all check exactly. Two corrections, neither
  load-bearing for the kill: a wrong atom count for the n = 128 certificate
  (252, not 253), and P4's claim that the violating λ-window sits "exactly
  where 009/011's λ-window law puts the workable range" — it does not: at
  n = 96 that law puts the assembly's workable λ at ≈ 0.052, and the
  violation lives at λ ∈ [2, 2.5], far outside it. The formal ∀λ statement
  is dead; the λ-restricted restatement (016's own lead 4) is more alive
  than the record's P4 framing suggests.
- **Tools:** `explore/uc_agg_ctrl_skeptic.py` (written here; stdlib only;
  imports NOTHING from `uc_or_avg*.py` / `uc_agg_ctrl_probe*.py`; instance
  parameters are read as data from 016's checkpoints, and the MU(n,r)
  builder is re-implemented from 016 §1's prose then cross-checked
  atom-by-atom against the probe builder in a separate snippet;
  deterministic, no randomness; full run with `--n128 --n160` ≈ 15 s).
  Command reproducing every number below:
  `python problems/union-closed/explore/uc_agg_ctrl_skeptic.py --n128 --n160`
  Checkpoints: `data/aggsk15_selftest.json`, `data/aggsk15_results.json`.
- **Sources:** 007/013 (definition of M_i, w_i, the aggregate, the anchors),
  016 (the record under attack and its checkpoints), 005 (Prop 1 anchor),
  009/011 (λ-window law, for the design audit), STATUS.md verification
  standard. No external fetches.

Notation as in 007/013: tilt coupling π ∝ u(A)u(B)·t^{|A∩B|}, t = 2^λ;
M_i the mass-weighted mean of log₂ OR_i over N_i²; w_i the defined history
mass; A(μ,λ) = Σ_{i<n, N_i≠∅} w_i(M_i − λ) / Σ w_i; record threshold
0.38271. At rational t the sign of A is the sign of
Σ_i Σ_rows m·log₂(OR/t) with every OR/t exactly rational.

## Claims attacked

1. **Quantity fidelity** (highest value): is 016's A(μ,λ) exactly the
   aggregate 013 certified positive — same M_i (007 §1's well-posed
   definition), same w_i weights, same N_i² block, same normalization
   pinned by the i = n and product anchors, same exclusion of i = n and of
   empty N_i? Killing a subtly different aggregate is not a kill.
2. **Family validity:** is the certified n = 96 instance really what the
   record says — 188 atoms, the MU construction of §1, all elementwise
   marginals < 0.38271 (exactly), H(μ) ≈ 3 bits, degeneracy dichotomy clean
   — and does anything hinge on MU being union-closed?
3. **The kill certificates themselves** (n = 96 the load-bearing one;
   n = 128, 160 the depth claims): recompute A from scratch in exact
   arithmetic with an independent log₂ enclosure and confirm or refute the
   certified signs.
4. **The positive controls:** the +1.844669 (013 part C) anchor
   reproduction, the padding invariance MU(n,1) ≡ witness, the P3-best
   n = 32 re-certification (and its phantom-zero-atom story), and the
   "surviving positive space" claim that the RAW-witness-weight ladder is
   still positive at n = 96.
5. **Robustness claims:** the one-line convention-robustness proof (sign
   invariant under 013's score-degenerate-at-λ bookkeeping); the meaning
   and scope of the 20/20 perturbation-stability claim.
6. **Design audit:** does the λ = 2 (t = 4) kill close the proof branch as
   claimed given 009/011's λ-window law; is the scope/"not claimed" section
   accurate; does the record overclaim anywhere?
7. **Enclosure smell test** (verification standard): could the tiny n = 96
   value (−7.6e-4) be an artifact of the enclosure method itself?

## Refutations found

**None load-bearing.** The kill survives every attack. Two corrections:

- **C1 (reporting).** The n = 128 tidy-rational certificate instance has
  **252 atoms** (2r + 8 with r = 122), not the "253 atoms" printed in 016's
  §3-P5/P7 block and Outcome table. Verified by reconstructing the instance
  from `aggprobe2_partP2.json` (my builder and the probe's agree
  atom-by-atom); `limit_denominator` dropped nothing (minimum θ* weight
  0.0017). Nothing depends on the count; the certified value re-derives on
  the 252-atom instance to all printed digits.
- **C2 (design framing, the one that matters downstream).** P4's sentence
  "the family's violating window sits exactly at moderate λ where 009/011's
  λ-window law puts the workable range" is **wrong at the n where the
  family violates**. The law is λ_max ≈ 4.847/(n−3) with sup CR at λ = 0
  for n ≥ 14 (011): at n = 96 the assembly's workable window is
  λ ≲ 0.052, and at every n ≥ 8 it excludes λ = 2. The violating window
  ([≈2, ≈2.5] at n = 96, positive at λ ≤ 1.5 and λ ≥ 3 per P4's own data)
  and the workable window are already separated by ~40× at the smallest
  violating n. Consequence, stated precisely: the certified kill refutes
  the recorded ∀μ∀λ statement (`aggregated-or-control`, exactly as queue
  item 1 posed it) — that much stands — but it does **not** touch the
  λ-restricted regime the assembly actually uses at these n, so 016's lead
  4 (λ ≲ c/n restatement) is not merely "possible", it is untouched by
  this family: my own certificate shows the same μ at n = 96 is positive
  at small λ (P4 float data; and the λ-profile positivity at λ ≤ 1.5 means
  no kill exists there for this family). The Ledger-consequence paragraph's
  list of live candidates (margin-modulated AND fixed-n/decaying-λ
  restatements) is accurate; P4's "sits exactly at the workable range"
  framing should not be quoted.

## Claims that survive (and what was done to break them)

### 1. Quantity fidelity (claim 1) — CONFIRMED; the aggregate is 013's

Re-derived the definition from 007 §1's prose alone (degeneracy dichotomy,
N_i from support reachability, conditional 2×2 tables, mass weighting,
normalized mean forced by the i = n and product-μ anchors; aggregate per
007 §5/013 part C: i < n, N_i ≠ ∅, weights w_i) and implemented it
from scratch in scaled-integer arithmetic (weights × lcm of denominators,
kernel × den(t)^n — no rounding anywhere outside the log₂ enclosure).
Fidelity evidence, all from my engine:

- **Anchor:** witness at t = 181/16 → A ∈ [+1.844669005341 ± 6e-21],
  matching 013 part C's certified +1.844669005 digit-for-digit. Since 013
  certified this via a third census implementation that 007's float engine
  agreed with, all four implementations now agree on the anchor.
- **i = n anchor (005 Prop 1) checked exactly:** at i = n every
  nondegenerate history has OR/t = 1 as a rational identity — verified on
  every instance run (`prop1_ok`), so the excluded i = n term is exactly 0
  and its exclusion is normalization, not choice. Including it anyway only
  adds w_n to the denominator: the n = 96 aggregate becomes
  −0.000759233..., still certified negative (computed).
- **Normalization robustness:** w_i is the globally-normalized defined
  mass; any common normalization cancels in A. Σ w_i at n = 96 is
  47.9006 (matches 016's "≈ 47.9"); per-i census structure matches the
  record's construction exactly (|N_i| = 1 at i ∈ {1,2,3,4}, |N_i| = 2 at
  the 90 response coordinates and at i = n−1; 368 nondegenerate rows).
- **Padding invariance as an exact identity:** MU(20,1) and MU(32,1) at
  t = 181/16 give enclosures with **zero** difference from the witness
  (identical rational bounds, identical W/Z) — 016's part-A anchor holds
  not to 0.00e+00 float but as exact rational equality.

### 2. Family validity (claim 2) — CONFIRMED

- **Construction:** my MU(n,r) builder, written from 016 §1's prose
  (coordinate j = bit j−1; block a {1,n−1}/{1,n}/light {1,c_j,n}; block b
  ∅/{n−1}/{n}/responder {c_j,n−1}; dilution {2},{3},{4}; c_j = coordinate
  4+j+1), reproduces the probe's `mu_ladder_theta` **atom-for-atom with
  exactly equal dyadic weights** at (7,1), (96,90), (128,122); the n = 160
  instance was taken directly from `aggprobe2_partP7.json`'s stored atom
  list (316 atoms = 2·154 + 8 ✓). The prose and the certified data agree.
- **Marginals, exactly:** every elementwise marginal satisfies
  marg·100000 < 38271·Z as an integer inequality at n = 96 (max
  0.30867, coordinate 95), n = 128 (0.30793), n = 160 (0.30719) — the
  record's "marginals ≤ 0.309 < 0.38271" is certified, not float.
- **Not degenerate:** H(μ) = 2.97–3.03 bits; defined mass ≈ 0.50 of
  history mass per coordinate; dichotomy holds exactly at every coordinate
  of every instance (nondegenerate table ⟺ both prefixes active — checked
  as an integer statement, no tolerance).
- **Union-closedness is not a hypothesis and its absence is not a flaw:**
  the statement under test quantifies over ALL μ on 2^[n] with H(μ) > 0
  ("for all μ and all i", 005 lead 1 as restated by 007; the control is a
  lemma about arbitrary tilt marginals, applied downstream to
  union-closed-supported measures). MU's support is not union-closed —
  and neither was 007's 10-atom witness that 013 certified. Checked and
  reported by the engine (`union_closed_support: false` for every
  instance including the witness) so no future reader mistakes this for
  an unexamined gap.

### 3. The kill certificates (claim 3) — CONFIRMED from scratch; the kill is real

My engine's certified enclosures (interval-squaring log₂, width ≤ 1.3e-20
per instance; every OR/t handled as an exact rational; dedup by value —
277 distinct log arguments at n = 96):

    n = 96,  188 atoms, t = 4:  A ∈ [−0.000759865185844166092,
                                     −0.000759865185844166079]  KILL (1.1s)
    n = 128, 252 atoms, t = 4:  A ∈ [−0.015405726918 ± 7e-21]   KILL (3.2s)
    n = 160, 316 atoms, t = 4:  A ∈ [−0.024225909315 ± 7e-21]   KILL (7.2s)

Each interval **contains 016's certified enclosure** (n = 96:
[−...087, −...086]; n = 128: −0.0154057269180687; n = 160:
−0.02422590931474056) and is strictly negative, with in-regime marginals
and exact dichotomy. Two structurally different enclosure methods (013's
atanh-series with 1e-50 directed rationalization vs my dyadic
interval-squaring with the ln 2 series bound) and two independently
written censuses agree to 18+ digits. The sign at n = 96 is a theorem
about the stated 188 rational weights. λ = 2 exactly (t = 4): no
irrational number appears in the statement.

### 4. The positive controls (claim 4) — CONFIRMED; the probe is not broken-negative

- Witness anchor: above (+1.844669005341).
- **P3-best n = 32** (tidy, zero-dropped, 36 atoms): my engine certifies
  +0.297679573398 at t = 4, dichotomy exact — matching 016's re-cert and
  confirming the phantom-zero-atom postmortem (the `limit_denominator`
  zero-drop is the correct fix; with zeros dropped the dichotomy flag is
  clean in my fully independent bookkeeping too).
- **"Surviving positive space" is real:** the RAW-witness-weight ladder at
  n = 96 with P6's adversarial dilution weight (wdil = 41.6935…) certifies
  **+0.157131854664** at t = 4 — my certificate upgrades 016's float
  +0.157132 (P6 was float-only) to a certified positive. The kill is
  genuinely a re-weighting result, not pure replication, exactly as the
  record's "calibrates how much of the kill is re-weighting" lead assumes.
- Round-1 exact positives cross-checked against `aggprobe_partE.json`
  (witness +1.8446690053408297, MU(22)=MU(32) +1.7307254294897931,
  MU(32)_t4 +1.0184774056855224, n = 8 endpoint +0.1661008241): the
  record quotes them faithfully.

### 5. Robustness claims (claim 5) — CONFIRMED as scoped

- **Convention robustness:** re-derived the one-line proof. Under 013's
  alternative bookkeeping (every degenerate positive-mass history scored
  at exactly λ), each such history contributes m·(λ−λ) = 0 to the
  numerator, and the per-i denominator grows from w_i to the full history
  mass; so A_alt = N/((n−1)·Z) with the SAME numerator N. Sign invariance
  is rigorous, and computed: A_alt at n = 96 ∈ [−0.000383136592916716685,
  −0.000383136592916716678], certified negative (128: −0.007751;
  160: −0.012168). No bookkeeping convention rescues the statement on
  these μ (only a convention awarding degenerate histories MORE than λ
  could, which 013 already ruled indefensible).
- **Perturbation claim, scope pinned:** what was perturbed is the 316
  potential weights (3% multiplicative noise, 20 trials, n = 160),
  evaluated in FLOAT (`aggprobe2_partP7.json` rows, range −0.0253..−0.0208,
  all in-regime). This is EVIDENCE that the violation is an open set, not
  part of the certificate — correctly so, since the certified statements
  are about the exact stated weights and need no stability. No issue; a
  reader should just not cite the 20/20 as certified.

### 6. Design audit (claim 6) — kill valid for the stated gap; one framing corrected

- The recorded gap (`aggregated-or-control`, 007 lead 1 = 013 claim 6 =
  queue item 1) is the ∀μ ∀λ>0 statement; one certified (μ, λ) pair
  refutes it, exactly as 007/013's per-i kill worked. The kill closes
  queue item 1's "if it survives n ≳ 20, attack the proof" branch as
  claimed: do not invest in the Gram/Frobenius proof of the unrestricted
  aggregate.
- But see **C2**: the violating λ-window and the assembly's workable
  window (λ_max ≈ 4.847/(n−3)) are disjoint at every violating n, so the
  practical bite on the route is narrower than P4's framing suggests, and
  lead 4 (λ-restricted control) is the natural next Gap-1 candidate
  alongside the margin-modulated one. The Outcome/Not-claimed section
  itself is accurate (it claims no ∀λ violation and flags the window);
  the record's calibration survives except for the one P4 sentence.
- The float-crossing claim "(88, 96]" matches P7's n-scan data (+0.00456
  at n = 88, −0.00076 at 96) — float EVIDENCE, correctly labeled.
- Search-completeness claims (n ≤ 32 clean, 58+6 climbs, permutation
  orbits) were spot-audited against checkpoints for consistency but NOT
  re-run — they are scoped as EVIDENCE in 016 and carry no weight in the
  kill; the certified n ≤ 32 positives (witness, P3-best, ladder points)
  are what backs "a fixed-small-n control is untouched", and those I
  re-certified.

### 7. Enclosure smell test (claim 7) — the sign is 16 orders above the noise floor

My log₂ enclosure self-test (`aggsk15_selftest.json`): exact on powers of
2 (width 0), correct on log₂ 3, log₂ 10, log₂(181/16), (7/5)^10,
1 + 1e-40 (width ≤ 1.5e-20 each), additivity-consistent, monotone. The
n = 96 aggregate enclosure width is 1.3e-20 against a certified value of
−7.6e-4: the sign has ~16 orders of margin. Underflow/cell-floor
pathologies are structurally impossible in the scaled-integer census
(every kernel cell is a positive integer; the 007-genre IEEE artifact
cannot occur). The dynamic-range worry from 007's detour does not apply
to the certificates at all — it applies only to the float searches that
found the instances, and their endpoints were then certified exactly.

## Verdict

| # | 016 claim | Verdict |
|---|-----------|---------|
| 1 | A(μ,λ) is 013's certified aggregate (007 §1 M_i, w_i weights, i < n) | **CONFIRMED** — from-scratch third implementation matches the +1.844669005 anchor and every certificate; i = n exclusion verified as an exact identity (Prop 1), inclusion variant still negative |
| 2 | n = 96/128/160 instances valid: MU prose = data, marginals < 0.38271, dichotomy, H | **CONFIRMED** (atom-for-atom builder match; marginals certified as integer inequalities, max 0.309; dichotomy exact; union-closure correctly not a hypothesis — the statement quantifies over all μ) |
| 3 | Certified kills: −0.000759865…, −0.0154057…, −0.0242259… at t = 4 | **CONFIRMED — the kill stands.** Independent census + independent log₂ enclosure re-derive all three signs with matching digits (widths ≤ 1.3e-20) |
| 4 | Positive controls (+1.844669, P3-best +0.29768, raw ladder positive) | **CONFIRMED**; raw-witness ladder at n = 96 upgraded from float to certified +0.157131854664 |
| 5 | Convention robustness; perturbation stability | **CONFIRMED** (one-line proof re-derived; A_alt certified negative at all three n; perturbation claim is float EVIDENCE, correctly scoped) |
| 6 | Ledger consequence: `aggregated-or-control` (∀λ) dies; branch closed | **CONFIRMED with C2**: the kill is valid for the stated gap, but the violating λ-window is disjoint from the assembly's workable window at every violating n — the λ-restricted restatement (lead 4) is untouched by this family and joins the margin-modulated control as a live Gap-1 candidate |
| 7 | "253 atoms" at n = 128 | **REFUTED (cosmetic)** — 252 atoms (C1) |

**Net assessment:** 016's headline — the i-aggregated odds-ratio control
(the restated Gap 1) is REFUTED, including in the record-relevant marginal
regime, by the θ*-optimized replicated-witness ladder at n ≥ 96, certified
float-free at λ = 2 exactly — should be treated as **VERIFIED**. The
certificates are now backed by two independently written exact censuses
and two structurally different certified log₂ enclosures, closing the
shared-code residual risks 013 recorded. The record's self-calibration was
good: the parts it flagged for attack first (the exact certificates, the
definition match, the anchors) all survived; the one genuine correction is
in the downstream framing (C2), not the mathematics.

## Residual risk

- **Shared-definition risk, further reduced but not zero:** all five
  engines formalize "conditional OR of the tilt coupling given histories"
  from 003/005's definitions. My implementation was written from 007 §1's
  prose, not from code, and hits the forced anchors (M_n = λ exactly,
  padding invariance exactly) — but the prose itself is the shared
  ancestor. A conceptual error in the *definition* of the route's Gap 1
  (as opposed to its computation) would pass every check here; 013's
  mitigation (006 S1's slice/inner-product derivation) still carries that
  weight.
- **Search-completeness (n ≤ 32 positive space)** was audited for
  checkpoint consistency, not re-run. It is EVIDENCE in 016 and remains
  EVIDENCE.
- **Float-level side claims** (crossing in (88, 96], perturbation 20/20,
  P4 window edges, P6 sweep minima other than the one certified here)
  were reconciled against checkpoints but not certified; only the n = 96
  raw-ladder point was upgraded to certified.
- My log₂ enclosure's correctness argument (interval squaring + crude
  endpoint bounds divided by 2^64, ln 2 series with geometric tail) is
  ~40 lines and self-tested, but has itself not been reviewed by anyone
  else. It agrees with 013's independent enclosure on every instance both
  have touched, which is the practical mitigation.

## References

- Reviewed record: `problems/union-closed/attempts/016-aggregated-control-probe-large-n.md`;
  its tools `explore/uc_agg_ctrl_probe.py`, `explore/uc_agg_ctrl_probe2.py`;
  checkpoints `data/aggprobe2_partP[1-7].json`, `data/aggprobe2_n96cert.json`,
  `data/aggprobe2_p3best_cert.json`, `data/aggprobe_part[ABCDE]*.json`, logs
  `data/aggprobe2_*.log`.
- Context: `attempts/007-averaged-or-control.md` (§1 definition, §5
  aggregate, witness); `attempts/013-skeptic-review-of-007.md` (part C
  anchor, exact standard, alternative bookkeeping, residual risks);
  `attempts/005-odds-ratio-control-refuted.md` (Prop 1, "for all μ"
  quantifier); `attempts/009/011` (λ-window law λ_max ≈ 4.847/(n−3), 011's
  measured λ_max table); `attempts/012` (large-n lesson, cited by 016).
- This review's tool/data: `explore/uc_agg_ctrl_skeptic.py`;
  `data/aggsk15_selftest.json`, `data/aggsk15_results.json`.
- No external sources; the ln 2 series Σ 1/(k·2^k) and the bounds
  1 − 1/v ≤ ln v ≤ v − 1 are standard calculus, re-derived in the tool's
  docstrings.
