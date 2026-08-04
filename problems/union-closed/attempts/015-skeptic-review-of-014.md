# 015 — Skeptic review of 014 (aggregated OR probe to n = 32): independent re-certification

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed
- **Type:** adversarial verification of `014-aggregated-or-probe-large-n.md`
  (default stance: refute). 014's highest-value confessed vulnerability is
  that its 15 exact certificates rest on its own ~60-line directed-rounding
  enclosure kit, reviewed by nobody. The decisive check here is therefore an
  **independent exact certifier built on a different algorithm**: binary
  digit-extraction (square-and-compare with directed integer shifts, sound-
  ness proved in the docstring) instead of 014's atanh series + dyadic
  compaction, integer-scaled row builders with an *exact* aggregate
  denominator, own orbit-formula re-derivation, own Sinkhorn (0.4-damped,
  not 014's sqrt damping), own direct full-support evaluator. Nothing on
  the critical path is shared with `uc_or_agg_probe.py`; it is imported
  only to reproduce the exact rational refit measures its orbit
  certificates are *about* (the measure spec is data, not verification
  logic) and as a comparison target.
- **Outcome in one line:** 014 holds: all five re-certified certificates
  (including the softest, +5.99e-7 at n = 32, and the fully rational t = 4
  point) reproduce inside 014's enclosures at my fixed width ~1.1e-14
  (~200x tighter on the n = 32 orbit certificates),
  the enclosure kit survives hand re-derivation and 437-case adversarial
  cross-bracketing, the orbit formulas re-derive and cross-check at n = 8, 9
  with exact within-class OR agreement, the dip row reproduces digit-for-
  digit and is not a grid artifact, and ~116 new instances in regions 014
  could not see (witness-direction off the product boundary, marginals
  pushed to 0.3827, lambda below the grid floor) found 0 violations.
  Corrections are reporting-level: the part-P "d^1.85–d^2.06" order range
  contradicts its own checkpoint (true range d^1.39–d^2.06, grid-edge
  artifact), and the "products are the exact equality set" phrase claims an
  unproven converse (false outright at lambda = 0).
- **Tools:** `explore/uc_or_agg_skeptic2.py` (parts K/X/C/T/V; stdlib only;
  deterministic; runtime ~3 min plus a ~1 min T-addendum inline script;
  checkpoints `data/oragg_sk2_part[KXCTV].json`, log
  `data/oragg_sk2_run.log`). Commands reproducing every number:

      python problems/union-closed/explore/uc_or_agg_skeptic2.py 2>&1 | \
          tee -a problems/union-closed/data/oragg_sk2_run.log
      # per part: K (kit attack), X (cross-checks), C (re-certification),
      # T (trend attack), V (new adversaries); the small-d order study is
      # the inline script logged at the end of oragg_sk2_run.log
      python problems/union-closed/explore/uc_or_agg_probe.py X   # re-run

- **Sources:** 005/006/007/012/013/014 (this repo). No external fetches.

Notation as in 007/013/014: tilt coupling `pi ∝ u(A)u(B) t^{|A∩B|}`,
`t = 2^lambda`; `M_i` = 007 §1's well-posed per-coordinate average;
`A(mu, lambda) = sum_{i<n} w_i (M_i − lambda) / sum w_i` with `w_i` the
defined mass; record threshold 0.38271.

## Claims attacked (014's own attack-order list)

1. **The certification kit** (`_dround` / `_ln_1_to_2_fp` / compacted
   `agg_enclosure`) and the 15 exact certificates, especially the softest
   (+5.988e-7, n = 32) and the fully rational t = 4 point.
2. **The orbit-engine formulas**: `Ncl`, `JW`, the cell formula, the
   j-independence of class ORs, the structural N_i reach test — cross-
   checked by 014 only at n = 6, 7.
3. **"Product measures are the exact equality set of A."**
4. **The surrogate-lambda sign-robustness argument.**
5. **The float trend claims**: growing margin in n; the n ≈ 48 dip;
   whether census minima are grid artifacts.

Plus the standard sweep: construct an in-regime violation in regions the
probe could not see.

## Refutations found

**None load-bearing.** Two corrections and two cosmetic findings:

- **R1 (correction — misreported order range).** Part P's prose and the
  record's one-line claim "all perturbative departures are quadratically
  upward … fitted order d^1.85–d^2.06 at small d" contradict 014's own
  checkpoint: over the 72 fitted cells in `or_agg_probe_partP.json` the
  order `log2(agg(0.01)/agg(0.005))` ranges **d^1.39–d^2.06**, with 19
  cells below 1.85 (worst: n = 32, full direction, p = 0.30, lambda = 0.5,
  in-regime, order 1.388). Re-computing that cell with my own engine at
  d down to 0.000625 gives successive-halving orders 1.388 → 1.510 →
  1.637 → 1.749 (and 1.624 → 1.734 → 1.820 for the crash direction):
  consistent with a quadratic leading term plus a large third-order
  correction, so the *qualitative* "quadratic upward departure" reading
  survives in the d → 0 limit — but the quoted range is a grid-edge
  artifact of exactly the genre 011 corrected in 009 (grid-edge values
  quoted as asymptotics). Consequence for lead 1 recorded below.
- **R2 (correction — overstated equality-set claim).** "Product measures
  are the **exact equality set**" (part O heading, one_line) claims a
  biconditional. Only the forward direction is proved — and I re-derived
  it by hand and verified it independently (product potential `x^|A|`,
  my direct evaluator: A = +4.1e-15 at n = 6; Sinkhorn-fitted
  Bern(0.3823)^8: A = +1.3e-14, max |M_i − lambda| = 1.1e-13). The
  converse (A = 0 ⟹ product) is nowhere proven or tested in 014, and at
  **lambda = 0 it is false outright** (OR ≡ 1, so A ≡ 0 for *every* mu —
  the same lambda > 0 proviso 006 had to add to 005 Prop 2's equality
  case). Downstream uses are unaffected: part E's exclusion of bern rows
  from certification and the "proof must be second-order at the boundary"
  framing need only the forward inclusion. Converse probe: the natural
  non-product candidate (the mirror family, which 006 S8 saw at per-i
  equality) fails by hand computation — its support is block-product
  {∅,{1}} × {∅,{2..n}}, the block target is product, so Sinkhorn gives a
  block-product potential, M_1 = lambda exactly but M_2 = (n−1) lambda
  > lambda, hence A > 0. The honest statement: products ⊆ equality set
  (all lambda); equality set = products is open for lambda > 0.
- **R3 (cosmetic — docstring drift, the 006-R3/013-R1 genre).**
  `uc_or_agg_probe.py`'s docstring says "O ~ 40-80 min … E ~ 10-30 min,
  X+S+H ~ 15 min"; the actual run (log + checkpoint `seconds` fields) is
  O ≈ 146 s, E ≈ 113 s, H ≈ 208 s. The record's *own* per-part times
  (~3 min, ~2 min) are correct; the docstring is stale.
- **R4 (cosmetic — instance-count slip).** The Outcome's "~1780 float
  instances" understates its own listed components: 840 + 432 + 40 + 426
  + 82 = 1820.

## Claims that survive (and what was done to break them)

### 1. The certification kit and the 15 certificates — CONFIRMED; five re-certified independently, exact agreement

*Hand re-derivation of the kit (part K docstrings carry the proofs):*
`_dround` floors/ceils onto ≤120-bit dyadics via integer floor division —
directionally sound including negatives (2000-case fuzz: 0 failures).
`_ln_1_to_2_fp`: the atanh partial sum is a lower bound (positive terms);
the tail bound `2 z^{2J+1}/((2J+1)(1−z²))` is the correct geometric
majorant; for the upper pass every intermediate is rounded up and the tail
uses the up-rounded z², which only enlarges it — outward everywhere. In
`log2_enclosure` the range reduction is exact (Fraction halving), and the
division by ln 2 picks the correct bound (LN2D_LO under a positive
numerator for the upper bound, etc.; the negative-numerator branches are
dead code but also correct). `agg_enclosure`'s compaction rounds every
running sum outward; the num/den interval ratio takes worst cases
independently, which is conservative (valid, slightly loose).
**One-sided verdict: sound.**

*Independent adversarial test (part K1):* my own enclosure —
digit-extraction `log2_bracket`, working mantissa 80 bits, 48 extracted
bits, ±(1,2)-ulp safety margins, soundness proof in the docstring, no
shared code or algorithm — against their `log2_enclosure` on 437
rationals: 400 random with numerators/denominators up to 200 bits, plus
2^k(1 ± 1e-25), 2^k(1 ± 1/3) for k ∈ {−100..100}, 181/16, 271/256.
Since my bracket provably contains the truth at width 1.1e-14, any
their-interval disjoint from mine convicts them: **0 failures** (their
max width 1.5e-13).

*Independent re-certification (part C):* five of the fifteen, on the
**identical exact rational measures** — sparse witnesses rebuilt from the
stated atom weights and layout; orbit refit measures mu-tilde rebuilt by
re-running 014's deterministic float Sinkhorn + 36-bit quantization as an
input spec, then everything downstream (row builder, integer scaling,
log2, marginals) mine:

    certificate                        mine [lo, hi]                          014's [lo, hi]                        overlap
    witness_n10_wdil16 (t=181/16)      [+1.70954246371, +1.70954246371]       [+1.70954246371, +1.70954246371]      yes
    witness_n20_t4 (lambda=2 exact)    [+0.483059568557, +0.483059568557]     same to all printed digits            yes
    witness_n32_wdil20 (t=181/16)      [+0.803594797389, +0.803594797389]     same to all printed digits            yes
    crashmix_383_05_n32 (t=271/256)    [+1.20186543655e-3, +1.20186543656e-3] [+1.2018654345e-3, +1.2018654366e-3]  yes
    softest: mixdir_slice p=.3 n=32    [+5.98823635652e-7, +5.9882364631e-7]  [+5.9881776522e-7, +5.9882370424e-7]  yes

  My width is 1.1e-14 everywhere (their orbit widths 2–6e-12; their
  sparse widths ~1e-15 are finer still); on the orbit certificates my
  intervals sit *inside* theirs and on the sparse ones the intervals
  agree to all printed digits; all signs positive, exact in-regime
  marginals confirmed
  (crashmix max marginal exactly 0.38235 < 0.38271 as a rational
  comparison). **No exact disagreement anywhere.** In addition a fresh
  certificate at the softest point on **my own** mu-tilde-prime (my
  Sinkhorn damping, my 40-bit quantization — a different exact rational
  measure near the same family point): [+5.9882e-7, +5.9882e-7],
  positive. The aggregate's insensitivity (10 digits) to which ~1e-11
  refit is used says the certified value is a property of the family
  point, not of the quantization.

### 2. The orbit-engine formulas — re-derived by hand; CONFIRMED at n = 8, 9

- *Re-derivation:* for tail-exchangeable mu and i ≥ 3, conditioning on
  prefix classes gives cell(α,β) = `t^{oa·ob + j} · t^{αβ} ·
  S[ja+α][jb+β]` with S the future pair-kernel contraction; the `t^j` and
  `t^{oa·ob}` factors cancel in the cross-ratio, so **OR = t · S00·S11 /
  (S01·S10) is j-independent** and the j-sum collapses into
  `JW[ja][jb] = Σ_j Ncl t^j` — 014's formulas are forced. `Ncl(q,ja,jb,j)
  = C(q,j)C(q−j,ja−j)C(q−ja,jb−j)` and `PK_f[s][l] = C(f,s) Σ_w
  C(s,w)C(f−s,l−w) t^w` were validated by brute-force subset enumeration
  (part K4: all q ≤ 6, f ≤ 5, rational t = 7/3 exact — 0 failures).
- *Cross-checks at n = 8 AND 9* (014 stopped at 6, 7), part X2: my orbit
  engine vs my own direct full-support evaluator vs 014's orbit engine on
  the same fitted potential, four families (crashmix, sawin_geo_a,
  slice-direction mixture, d0bern) × lambda ∈ {0.5, 1}: max disagreement
  across all 16 instances **1.1e-13**.
- *j-independence tested directly* (part X3): at n = 9, i = 6, per-history
  ORs from the direct evaluator grouped into (oa,ob,ja,jb) classes: 256
  classes, within-class relative spread **exactly 0.0**.
- *N_i reach test:* implied by the direct-engine agreement (my evaluator
  derives N_i from the actual support) and by the exact dichotomy asserts
  in both exact paths never firing.

### 3. The equality boundary — forward direction CONFIRMED (re-proved and re-measured); converse is R2

Two-line forward proof re-derived: a product potential factorizes the
kernel, so at every history both response slices are proportional, OR =
2^lambda pointwise (005 Prop 2 equality case), M_i ≡ lambda, A ≡ 0 — at
every lambda, every n. Verified with my independent evaluator (numbers in
R2). The *census-noise reading* of part O (bern rows at ±1e-14 are
boundary noise, not near-violations) is therefore right, and part E's
refusal to "certify" a boundary point is the correct move. The
biconditional phrasing is R2.

### 4. Surrogate-lambda sign-robustness — CONFIRMED (re-derived)

Under 013's alternative bookkeeping every margin-degenerate history
contributes `w·(lambda − lambda) = 0` to the numerator and `w > 0` to the
denominator; this also covers coordinates whose N_i is empty under
conditioning-out (they enter the surrogate denominator with zero
numerator). So the numerator of A is literally unchanged, the denominator
only grows: sign preserved, magnitude shrinks toward 0. One boundary case
014 skates over, harmless: if *no* coordinate has defined mass the
conditioning-out A is undefined while the surrogate A is 0 — no sign to
disagree about. Sound as stated.

### 5. The float trends — CONFIRMED; the dip and the minima are real features, not grid artifacts

- *Dip row reproduced digit-for-digit* by my engine (part T1): +6.150e-7 /
  +2.217e-7 / +1.408e-7 / +2.330e-7 / +4.305e-7 at n = 32..64 — matches
  014's checkpoints to 4 digits at every n, with independent Sinkhorn.
- *Not a lambda-grid artifact* (part T2): at the dip point (n = 48) lambda
  swept 0.005–0.5 including below 014's grid floor 0.079: A > 0
  throughout; A/lambda decreases monotonically to **+6.4e-7 as
  lambda → 0** — the first-order-in-lambda coefficient (007 Theorem C's
  perfect square) stays bounded away from 0. Not a d artifact (d ∈
  {0.01..0.2} positive), and p = 0.30 is the local minimizer over p ∈
  {0.26..0.34} at n = 48 (so 014's softest point is a genuine local
  extremum of the sampled family, not where its grid stopped).
- *Census-minimum growth in n confirmed off-grid* (part T3): crashmix_383_05
  at ten lambdas spanning 0.005–4.4 (including values between 014's big
  grid points): A monotone increasing in lambda at fixed n, and the
  margin **grows with n at every lambda tested** — the growth claim is
  not a property of the lambda ≈ 0.083 grid floor. The census minimum
  sitting at the smallest grid lambda is the O(lambda) scaling of A, as
  014 itself says, not a hidden interior extremum.
- Spot-verified from checkpoints against the record: per-n census minima
  (+4.88e-4 … +1.21e-3 at crashmix, smallest grid lambda), family
  ordering at n = 32, 660/840 in-regime, Sinkhorn residuals ≤ 1e-11 on
  all 840, max |M_n − lambda| = 1.39e-13, enclosure widths ≤ 5.94e-12,
  H best endpoint +0.005602, witness n·A → 1.7555 (my evaluator,
  independent). 014's part X re-run reproduces its checkpoint.

### New-adversary sweep (part V) — 0 in-regime violations in ~116 new instances

Chosen to hit exactly what 014's search could not see:

- **Witness-direction departure from the product boundary** (n = 10):
  `u = x^{|A|} + d · u_witness` in potential space — a non-tail-symmetric
  direction; 014's part P sampled only 4 symmetric rays and its record
  itself names the Hessian's direction space as the attack surface.
  24 cells (bulk marginal targets 0.30 / 0.375 × d ∈ {0.03, 0.3, 3} ×
  lambda ∈ {0.084, 0.5, 2, 3.5}): all positive, min +3.2e-8 at the
  smallest departure — consistent with quadratic upward curvature in this
  direction too.
- **Boundary-of-regime marginals** (max marginal pushed to 0.38264–0.38270
  by tuning crashmix p, d0bern components, and a scaled Sawin mixture, at
  n = 20, 32): all positive, and with *larger* margins than the interior
  census — the regime boundary is not where the aggregate is tight.
- **Tiny lambda** (0.002–0.04, below 014's grid floor) at the two tightest
  families, n = 20, 32: positive, cleanly linear in lambda.

## Verdict

| # | 014 claim | Verdict |
|---|-----------|---------|
| 1 | 15 exact certificates + enclosure kit | **CONFIRMED** — kit sound by hand re-derivation + 437-case adversarial cross-bracketing (0 failures); 5/15 re-certified on identical exact inputs by an independent algorithm, exact agreement, my tighter intervals nested inside theirs; fresh independent certificate at the softest point also positive |
| 2 | Orbit-engine formulas, j-independence, N_i test | **CONFIRMED** — hand re-derivation; combinatorics brute-forced; 3-way engine agreement ≤ 1.1e-13 at n = 8, 9; within-class OR spread exactly 0 |
| 3 | Products are the exact equality set | **CONFIRMED as an inclusion, corrected as stated** (R2): forward re-proved and re-measured; converse unproven and false at lambda = 0; mirror-candidate fails to break it |
| 4 | Surrogate-lambda sign-robustness | **CONFIRMED** (re-derived, one harmless undefined-A boundary case noted) |
| 5 | Trend claims (margin growth, dip-and-recover) | **CONFIRMED** — dip reproduced digit-for-digit; shown robust off-grid in lambda, d, p; growth-in-n holds at every off-grid lambda tested; part-P order range corrected (R1) |
| — | 0 violations in ~1780 instances | **CONFIRMED and extended** — count is 1820 (R4); +~116 new instances in unseen regions, 0 violations |

**Net assessment:** 014's headline — the i-aggregated OR control survives
to n = 32 certified and to n = 64 on trends, with growing margin and no
012-style reversal — should be treated as **VERIFIED at its stated scope**
(finite n, finite grids, tested genres; certificates are sign theorems
about the stated exact rational measures). The queue-item verdict (proof
effort now justified, lead 1 = product-boundary Hessian) stands, with one
sharpening from R1: at lambda ≳ 0.5 the third-order term in d is large
enough to halve the fitted order at d = 0.005, so **a Hessian-PSD result
alone will not control the boundary layer at finite d — lead 1 needs a
quantitative cubic remainder bound**, or a restriction to the small-d
regime it certifies.

## Why it failed / what survived

Nothing failed in the reviewed record beyond reporting-level slips; the
review's own kill attempts failed for identifiable reasons worth keeping:

- The kit attack failed because every rounding in 014's kit really is
  outward; the compaction idea (directed dyadic fixed-point) is sound and
  is *the* reason exact certification is cheap at n = 32. My digit-
  extraction bracket holds a fixed 1.1e-14 width at similar cost (~200x
  tighter on the n = 32 orbit certificates) — worth adopting for future
  certificates where margins are ~1e-7.
- The witness-direction boundary attack failed with clean quadratic
  positivity: even the one genre known to produce per-i deficits
  (M_i − lambda = −0.122) enters the aggregate quadratically upward when
  diluted into a product bulk. The bulk floods every prefix, so the
  gadget's deficit rows carry vanishing relative mass — same mechanism as
  014's Theta(1/n) witness-genre decay, seen from the boundary side.
- The regime-boundary attack failed because pushing marginals toward
  0.38271 *raises* the aggregate on every family tried: the tight
  direction is small lambda near products, not large marginals.

Residual risk, labelled: all engines (005–014 and this one) formalize the
same census definition from 003; exact arithmetic and re-derivation remove
numerical and algebraic error, not a shared conceptual error in that
definition (already flagged in 013; still true). The orbit certificates
remain statements about quantized refit measures (013 R2 scope), which is
the right scope for the universally-quantified claim.

Reusable: `explore/uc_or_agg_skeptic2.py` — the digit-extraction log2
bracket (`log2_bracket`, soundness proof inline, width 1.1e-14), the
integer-scaled exact orbit/sparse certifiers with exact denominators, the
from-scratch direct evaluator, and the boundary/witness-direction
instance builders. Checkpoints carry every number
(`data/oragg_sk2_part[KXCTV].json`).

## Leads generated

1. **(sharpens 014 lead 1)** The Hessian program at the product boundary
   needs a cubic remainder bound: measured fitted orders at d = 0.005 drop
   to 1.39 (full direction, lambda = 0.5, n = 32) before recovering
   toward 2 as d → 0. Compute the third derivative's growth in (n,
   lambda) along the full/crash rays with the orbit engine before
   trusting any finite-d positivity radius from the Hessian alone.
2. **Settle the equality-set converse for lambda > 0** (from R2): is
   A(mu, lambda) = 0 only at (block-degenerate) products? For ≤ 4-atom
   supports Theorem B reduces it to per-i equality, where 005 Prop 2's
   proportional-slices characterization applies; the full-support case
   looks like a rank-1-slices induction. A non-product equality point
   would be a second tangency of the conjectured inequality that the
   Hessian program would silently miss; a proof closes the boundary
   description in 014's proof plan.
3. The first-order coefficient at the dip point (A/lambda → +6.4e-7 at
   n = 48) is the smallest Theorem-C perfect-square value seen anywhere;
   014's lead 2 (map the dip in (p, layer, n)) should track this
   lambda → 0 coefficient rather than A at a fixed small lambda — it is
   the quantity that would actually cross zero first if the control fails
   along this family.

## References

- Reviewed record: `problems/union-closed/attempts/014-aggregated-or-probe-large-n.md`;
  its tool `explore/uc_or_agg_probe.py` and checkpoints
  `data/or_agg_probe_part[XOPTSHE]*.json`, `data/or_agg_probe_run.log`.
- Context: `attempts/007-averaged-or-control.md` (M_i definition, witness,
  dilution invariance); `attempts/013-skeptic-review-of-007.md` (exact
  certification standard, refit scope, surrogate-lambda);
  `attempts/012-skeptic-review-of-008.md` (orbit symmetry, reversal
  warning); `attempts/005-*`/`006-*` (Props 1–2 and the lambda > 0
  equality proviso reused in R2).
- This review's tool/data: `explore/uc_or_agg_skeptic2.py`;
  `data/oragg_sk2_part[KXCTV].json`; `data/oragg_sk2_run.log`.
- No external sources; the digit-extraction log2 soundness argument is
  elementary and written out in the tool's docstrings.
