# 012 — Skeptic review of 008 (perturbative assembly): adversarial verification

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed
- **Type:** adversarial verification of `008-perturbative-assembly.md`
  (default stance: refute). Every lemma and the two theorems re-derived by
  hand from the 003/005/006 definitions; every numeric family re-computed
  with independent implementations (`explore/uc_pert_skeptic.py`, nothing
  imported from `uc_pert.py`): a from-scratch direct engine (alternating
  scaling, both-margin residuals, coupling symmetry *checked* not assumed,
  stable conjugate-root `z_ρ`), plus an **orbit engine** — profile-space
  Sinkhorn exploiting the `S_{n−2}` symmetry of crash-direction mixtures
  (states `(a₁,a₂,k)`, exact combinatorial census) — structurally unrelated
  to 008's dense atom engine and able to run 008's budget census to
  **n = 32**, far past 008's n ≤ 9. The two engines agree to 1.2e-14 at
  n = 6, 7 (SK0).
- **Outcome in one line:** the proof content of 008 (Lemmas A–C, Theorem
  P6′) survives line-by-line re-derivation and a 120-instance numerical
  stress test, and the 78-run assembly table is genuine (spot-agreement
  ≤ 5e-6; no capped negative found even on a new adversarial direction) —
  but the conditional theorem's constants are wrong in three distinct ways
  (τ_half is not box-uniform at 008's own s₀, analytically; the stated δ₀
  formula was not the one computed; δ₀ ≈ 0.022 should be ≤ 0.0079 by 008's
  own formula and ≈ 0.004 corrected), and the census extension **reverses
  008's flattening reading: past 008's window the averaged downward
  budget's increments stop shrinking and start growing (both regimes;
  claimed limits passed by n = 16 / n = 10), so budget (B1) with an
  n-uniform constant is dead at evidence level by 008's own lead-1
  criterion, and the tax budget (B2) has the wrong δ-shape entirely
  (first-order, not δ²)**.
- **Tools:** `explore/uc_pert_skeptic.py` (parts SK0–SK7; stdlib only;
  deterministic, no RNG; runtime ~90 s total; checkpoints
  `data/pertsk_SK0.json` … `pertsk_SK7.json`). Commands reproducing every
  number below: `python problems/union-closed/explore/uc_pert_skeptic.py
  --part SKn` for SKn in {SK0, SK1, SK2, SK3, SK4, SK5, SK5B, SK5C, SK6,
  SK7} (checkpoints refuse overwrite by design; delete `data/pertsk_*` from
  a scratch copy to re-run).
- **Sources:** 008 (the record under review), 003/004/005/006 (definitions
  and corrected statements), `explore/uc_pert.py` and `data/pert_*.json` /
  `pert_full_run.log` (reviewed artifacts), STATUS.md (queue framing for
  claim 4). No external fetches this cycle.

Notation as in 008: tilt `π_λ(A,B) ∝ u(A)u(B)2^{λ|A∩B|}`, `ρ = 2^λ`,
`x₀ = 1−p`, `m₀ = min(p,1−p)^n`, capped δ-ball `C_n(p,δ)`,
`G(x,y,τ) = h(z_{2^{λ+τ}}(x,y)) − ½h(x) − ½h(y)`, budgets (B1)–(B4) as in
008's conditional theorem.

## Claims attacked

1. **Theorem P6′ and Lemmas A–C** — every step re-derived by hand: strict
   convexity of the scaling functional (K ⪰ 0 via tensor factors),
   coercivity along rays, the analytic IFT and the resolvent bound
   `‖(D_μ+Π)^{−1}‖ ≤ 1/min μ`, the second-derivative bound of Lemma C, the
   ω₀-measure-identical-across-cells step, the alternating-sum first-order
   cancellation, the Taylor remainder, and the 924 constant chase.
2. **The constant `924·min(p,1−p)^{−3n}`** — re-derived, then stress-tested
   with the independent direct engine on 120 instances (six (n,p,λ)
   triples × five directions incl. δ_full and a new half-empty direction ×
   δ up to the theorem's boundary 0.5), plus Lemma B's differential at an
   *interior* δ (008 checked δ = 0 only) and Lemma C's sup-bounds.
3. **Lead-3 / lead-4 distinctness** (against 005's text) and the lead-4
   closed form `log₂OR_crash = λ(3−n)/n` for the centered kernel — by hand
   and by two independent computations (potential-free identity + Sinkhorn).
4. **ρ*(0.383) = 1.0422 vs the queue's ρ ≈ 1.03; g(1.031, 0.383) < 0.**
5. **The budget census claims:** pointwise +0.3 bits/coordinate upward-only
   non-uniformity; max-coordinate downward ≈ 0.083(n−2)δ²; the averaged
   downward flattening (→ ≈ 0.07δ² assembly regime, ≈ 0.3δ² at (0.30, 1.0))
   — extended from n ≤ 9 to **n = 32** with the orbit engine, including
   δ → 0 scaling checks at n = 32 (δ down to 0.003125); tax and
   margin-variance shapes for (B2)/(B3).
6. **The conditional theorem's assembly:** the (B1)–(B4) ⟹ (S-coup)
   implication step by step; the box-uniform τ_half step (008's own flag);
   the P2d constants M_τ, M₂, L₂, dip; the δ₀ ≈ 0.022 arithmetic.
7. **The 78-run assembly table** and the "only negatives are cap-violating"
   claim — spot re-computation plus a 104-run hunt for a capped negative in
   cells 008 never ran: the half-empty direction `½δ_∅ + ½Bern(2p)^n`
   (marginals exactly p — the 004 killer genre inside the class), ρ = 1.06
   cap-stress cells, ρ = 1.045 just above threshold, off-grid δ = 0.15.

## Refutations found

Three genuine defects in the conditional theorem (T4), one reversal of a
census reading, and two reporting slips. Nothing touches Theorem P6′.

### R1. The box-uniform τ_half step is false at 008's own (p, ρ, s₀) — analytically
The proof outline uses "G nondecreasing in τ up to τ_half", justified by
"z(x,y,τ) < ½ on the whole box for τ ≤ τ_half, which holds for s₀ = 0.10 at
this (p,ρ) by the monotonicity of z in (x,y,τ) — checked numerically on the
P2d grid". Both halves are wrong:

- **Analytically:** at the box corner, `xy = (x₀+s₀)² = 0.717² = 0.5141 >
  ½`, and `z > xy` for every ρ > 1, so `z > ½` already **at τ = 0**:
  τ_half(corner) = 0, not 3.017. Monotonicity of z in (x,y) is precisely
  why the corner, not the center, is the binding point — the cited
  monotonicity implies the *opposite* of what it is invoked for.
- **The claimed numerical check does not exist:** `uc_pert.py` part P2d
  computes sups of |∂G/∂τ|, curvature, and dip; it contains no sign or
  z < ½ check. On my grid (SK4), 266 of 625 box points have z > ½ at
  τ = 0.999·τ_half(center), with ∂G/∂τ as negative as −0.0173 there.

Consequence: on part of Γ, *upward* τ-deviations lose margin and the loss
is charged to no budget — the theorem as stated has a hole. Repair (SK7):
the step is saved by shrinking the box; z_1.2(x,x) crosses ½ at
x_c = 0.7013, so s₀ ≤ x_c − x₀ = **0.0843** is the ceiling at (0.383, 1.2)
(and τ_half, M_τ, M₂, ‖G‖ must then be re-measured on the smaller box —
the constants are coupled). The wrong-way slope is small, so the repaired
theorem plausibly holds with adjusted constants; as printed it does not.

### R2. δ₀ ≈ 0.022 does not follow from 008's own theorem — three constant errors
- **The formula computed is not the formula stated.** The theorem's
  `δ₀ = sqrt(g / 2(M_τc₁ + c₂ + (M₂ + ‖G‖_box/s₀²)c₃ + …))` includes a
  `‖G‖_box/s₀²·c₃` term; part P2d's code computes
  `sqrt(0.5g/(M_τc₁ + c₂ + M₂c₃))`, silently dropping it, and never
  computes ‖G‖_box at all. Measured (SK4): ‖G‖_box = 0.1456, so 008's own
  stated formula gives **δ₀ = 0.0079**, not 0.0216.
- **The off-box Markov charge needs a global bound, not ‖G‖_box.** The
  margin-tail mass lies *outside* the box, where sup(−G) over
  (0,1)² × [0, τ_half] is 0.487 (SK4), not ≤ ‖G‖_box = 0.146.
- **M_τ is measured on τ ∈ [0, τ_half] but charged against τ < 0 losses**
  ((B1) is a downward budget). sup|∂G/∂τ| on the box × τ ∈ [−6, 0] is
  **0.0838**, 42% above the 0.0592 used (attained at the lower corner,
  τ ≈ −3).
- Additionally c₃ = 2.5 is under-measured: crash-direction
  `E[(x̃−x₀)²]/δ²` reaches **4.70** by n = 16 at δ = 0.05 (SK5;
  n-saturating, so (B3)'s *shape* is fine — the constant isn't; 008's 2.16
  was read at δ = 0.1, where the δ² law has not converged).

With all corrections, δ₀ ≈ 0.0045 (0.0035 with c₃ = 4.7). The theorem's
*structure* (split over Γ, Markov the tails, budgets subtract) is sound
once R1's box repair is made; the headline number 0.022 is not, off by
3–6×. (008 did hedge: "constants rest on dense finite differences, not
interval arithmetic" — but the errors above are formula-level, not
rounding.)

### R3. The averaged-downward-budget "flattening" reverses beyond 008's window — (B1) as stated is dead at evidence level
008's part P3c (n ≤ 9, δ = 0.05, crash direction): increments of
`avg_i E[(λ−log₂OR_i)_+]/δ²` decreasing ~25%/step, read as "consistent
with a bounded limit (≈ 0.07–0.08)" at (0.383, ρ=1.2) and "flattening
toward ≈ 0.3δ²" at (0.30, λ=1). The orbit engine reproduces every one of
008's n ≤ 9 values to all printed digits (SK5: 0.0414…0.0652 and
0.1288…0.2739), then extends them (SK5/SK5B/SK5C/SK7):

- **(0.383, ρ=1.2):** avg/δ² = 0.0652 (n=9) → 0.0759 (16) → 0.0824 (22) →
  0.0889 (28) at δ = 0.05, with increment *ratios rising monotonically
  through 1* (0.677 … 0.948 … 1.011): past n ≈ 22 the increments
  **increase**. This is not a fixed-δ artifact: the δ → 0 coefficient is
  stable at n = 32 (0.1015 / 0.1032 / 0.1043 at δ = 0.0125 / 0.00625 /
  0.003125 — clean δ², ratios → 4) and the δ = 0.0125 sequence 0.0832,
  0.0886, 0.0947, 0.1015 at n = 20, 24, 28, 32 has *increasing* increments
  — growth worse than logarithmic everywhere in the measured window. The
  claimed limit 0.07–0.08 is already passed by n = 16–20.
- **(0.30, λ=1):** the average blows through the claimed ≈ 0.3 limit at
  n = 10 and reaches ≥ 1.24 by n = 32. The apparent turnover at fixed
  δ = 0.05 (peak ≈ 0.48 near n ≈ 24, SK5B) is a **δ-saturation artifact**:
  at δ = 0.0125 the sequence still climbs (0.956, 1.140, 1.226, 1.242 at
  n = 20, 24, 28, 32), and the δ-scan ratio (δ=0.0125 vs 0.025) worsens
  with n, so even these are lower bounds on c₁(n).

008's lead 1 named the falsification criterion itself: "extend part P3c to
n = 10..12 … if the increments stop decreasing, (B1) as stated dies and
the theorem needs a log n factor". The increments do not merely stop
decreasing — they turn increasing. **(B1) with an n-uniform constant is
REFUTED at EVIDENCE level (n ≤ 32, crash direction, δ² coefficient
verified stable in δ at n = 32); even a log n repair is not supported in
the measured window** (log growth would have decreasing increments).
Whether c₁(n) is eventually bounded remains open; nothing in n ≤ 32 points
that way. EVIDENCE, not proof, on both sides.

### R4. The tax budget (B2) has the wrong δ-shape: first-order, not δ²
At n = 32, (0.383, ρ=1.2), halving δ from 0.05 down to 0.003125 multiplies
the average tax by 2.67, 2.49, 2.34, 2.23 — tending to **2**, i.e. the tax
is asymptotically **linear in δ** in every observable window (tax/δ² grows
like 1/δ: 0.041 → 0.301 across that δ-range; same at (0.30, λ=1)). At
fixed δ = 0.0125 the coefficient also grows ~linearly in n (0.0657 →
0.0981 over n = 20 → 32), 10× the instantiated c₂ = 0.01 already.
Mechanism (SPECULATION, consistent with the scaling): histories whose
prefix is crash-dominated carry ~δ total mass with Θ(1)-non-product
conditional tables once the product mass of those prefixes (~x₀^i) drops
below δ·0.33 — the true δ² window at n = 32 is δ ≲ 1e-6, unobservable, and
shrinks exponentially in n (echoing 006 S7 and P6′'s exp(n) constant). The
*absolute* tax stays tiny (~1e-4 at δ = 0.05, n = 32, vs g = 0.0052), so
the end-to-end assembly is not threatened — but (B2) as stated (quadratic,
n-uniform c₂) is doubly unsupported; the budget must be restated
first-order, `avg tax ≤ c₂′δ`, which the theorem's accounting absorbs
(taxes subtract directly) at the price of δ₀ ~ g/c₂′ bookkeeping instead
of sqrt.

### R5. Reporting slip: the "worst in-class net" run is out of class
Outcome 2's "(worst in-class net = +0.0003, slice n=5 δ=0.2, still
positive)": that run is flagged `capped=False` in 008's own log (slice at
n = 5 has marginals 2/5 = 0.4 > 0.383). The true worst in-class
per-coordinate net across the 78 runs is ≈ +0.00105 (ρ = 1.06, crash,
n = 6, δ = 0.02) — i.e. the sentence simultaneously mislabels the class
and *understates* 008's own result (no capped run's min net ever dropped
below its base margin, verified in SK3).

### R6. Two constant-chase nits in P6′'s write-up (bound still valid)
- The proof's final step "dividing by ln 2 and using m₀ ≤ ½" gives
  (1024 + 128)/(2 ln 2) = **831**, not 924; the printed 924 corresponds to
  the weaker m₀ ≤ 1. The stated bound 924·m₀^{−3} is therefore still a
  correct upper bound — the cited step just doesn't produce it.
- "measured c ≈ 0.6 vs bound 6.4e10 at n = 5, p = 0.3" quotes the gentlest
  direction (005's part-G family); 008's own P3 census measured c ≈ 20.6
  (crash, same n, p, λ, δ = 0.025). The bound is wildly loose either way;
  the illustration understates the measured worst case by 34×.

## Claims that survive (and what was done to break them)

### Theorem P6′ and Lemmas A–C (claim 1) — the heaviest content of 008 stands
Hand re-derivation, from the 003/005 definitions, of every step:

- **Lemma A.** ∇Ψ = 0 ⟺ scaling equations ✓; Hessian identity
  `HessΨ = D(v) + D_{e^v}KD_{e^v}` re-computed entry-wise ✓; each kernel
  factor `[[1,1],[1,2^λ]]` has det 2^λ−1 ≥ 0, tensor products and principal
  submatrices of PSD are PSD, so Q ⪰ 0 — **valid for λ ≥ 0 only**, and P6′
  correctly claims λ ≥ 0 only (for λ < 0 the factors are indefinite and the
  convexity proof genuinely fails; noted as scope, not error). Coercivity
  case split ✓ (w with a positive entry: quadratic term beats linear; all
  entries ≤ 0: the linear term −Σμv grows, using μ > 0); ray-coercive
  convex ⟹ coercive is standard ✓.
- **Lemma B.** `D(v) = D_μ` and `Q = Π` *at the solution* ✓; analytic IFT
  applies since ∂F/∂v = D_μ + Π ≻ 0 ✓; the resolvent bound is operator
  monotonicity of the inverse applied to D_μ + Π ⪰ D_μ ✓ (needs Π ⪰ 0,
  i.e. again λ ≥ 0). At the simplex boundary the bound degrades as
  1/min μ, exactly as used — the mixture segment keeps min μ_δ ≥ m₀/2, and
  directions ν putting mass on μ₀'s smallest atoms cannot break this
  (μ_δ ≥ (1−δ)μ₀ pointwise). Exponent −3n in the final constant traced to
  the 1/m₀³ term ✓ (tightness not claimed by 008; not contested).
- **Lemma C.** F linear in μ kills F_μμ and F_vμ ✓; `μ̈ = 0` ✓; the
  three-term formula for (∂²F/∂v²)[h,h] re-derived and each term bounded
  by ‖h‖∞²μ_δ(A) (middle term 2×) ✓; ‖ν−μ₀‖₂ ≤ 2 ✓; 128/m₀³ chase ✓.
- **P6′.** The probability-measure trick: L″ = E_ω[v̈⊕v̈] + Var_ω(v̇⊕v̇)
  re-derived from s″/s − (s′/s)² ✓; Var ≤ osc²/4 ≤ 4‖v̇‖∞² ✓. The
  **ω₀-identity step** (task attack (a)): at δ = 0 uniqueness forces the
  product potential, so all four slices are proportional to the same future
  product H and `ω_{αβ,0} ∝ H(x)W(x,y)H(y)` for every (α,β) — the scalars
  cancel in the normalization; hence L′(0) = f(α) + g(β) and the
  alternating sum annihilates it ✓. R(0) = 1 by the product identity ✓;
  integral-remainder Taylor ✓. Full support of μ_δ makes every history
  nondegenerate ✓.

Numerical kill attempts (SK2): 120 instances — six (n,p,λ) including
p = 0.2 and λ = 2, five directions including δ_full (mass onto the
smallest-μ₀ corner) and half-empty, δ up to the theorem's boundary 0.5 —
**zero violations**; the bound is never approached (max dev/bound
4.7e-8). Lemma B's differential formula verified at an *interior* point
δ = 0.25 against finite differences (3.3e-4, FD-limited), not just at
δ = 0 as in 008's P1; Lemma C's sup-bounds hold with slack ~500× (|v̇|∞ =
3.2 vs 1646) — loose, as expected, never violated. CONFIRMED (with R6's
cosmetic constant-chase note). This closes 005's Prop-6 smoothness gap:
P6′ should now be regarded as VERIFIED at fixed n.

### The calibration block (claim 4) — 008 right, queue framing wrong
Independent recomputation (SK1): ρ*(0.383) = 1.042205, ρ*(0.38271) =
1.030222, g(1.031, 0.383) = −0.000417 < 0, g(1.06)/g(1.20) = +0.000649 /
+0.005252 — all match 008. STATUS.md's queue line (c) says "perturbative
assembly at ρ≈1.03" with no p; 1.03 is the *threshold at p = 0.38271*
(003's own number), infeasible as a running point at p = 0.383. 008's
recalibration to ρ ∈ {1.06, 1.20} was necessary and correct. The AHS
coefficient: implicit differentiation of the Plackett quadratic gives
ż·1 = xy(1−x)(1−y) exactly (re-derived), so dg/dρ|₁ at p = ψ is
h′(ψ)φ⁶ = 0.0386888, matched by central differences to 3.3e-10. CONFIRMED.

### Leads 3/4 distinctness and the lead-4 closed form (claim 3)
005's text: lead 3 = "Perturbative assembly, now unblocked at second
order" (sign of the net second-order coefficient); lead 4 = "Two-scale
tilts vs the crash family" (kernel change). Distinct as 008 says ✓ (the
"task line conflated them" refers to the cycle prompt, which is not in the
repo — unverifiable here; the queue line names only lead 3). The centered
kernel exponent arithmetic re-done by hand: e(S₁,S₃) = 0, e(S₂,S₄) =
−(n−2)/n, e(S₁,S₄) = −1/n, e(S₂,S₃) = 0, so log₂OR = λ(3−n)/n ✓, and
008's intermediate `(3−n) − (4n−3−n²)/n` simplifies to the same ✓; the
potential-free identity and an independent Sinkhorn fit both give it
exactly at n = 5, 8 (SK6), the mirror gives λ(n−1)/n = 0.8λ at n = 5 ✓,
and the product bands [0.8265, 0.8333] (n=6), [0.8695, 0.8750] (n=8)
reproduce to 4 decimals. CONFIRMED.

### The 78-run assembly table and the cap story (claims 2, 7)
Four spot rows (P2b/P2c, both ρ, incl. the negative slice row) reproduce
with the independent engine to ≤ 5e-6 (SK3a). The hunt (SK3b, 104 runs):
the **half-empty direction ½δ_∅ + ½Bern(0.766)^n** — marginals exactly p,
the "two blocks with opposite needs" genre that produced 004's 0.431496
ceiling, never tested by 008 — plus pm03/slice/embedded-crash at ρ = 1.06
(008 ran ρ = 1.06 only at n = 6 with three gentle directions), ρ = 1.045
(margin 1.3e-4, just above threshold), off-grid δ = 0.15, n to 9: **zero
capped instances with min_net ≤ 0 or gain_lb ≤ 0, zero OR < 1 mass
anywhere** (worst capped min_net at ρ = 1.045: +0.000115, still above the
base margin g(1.045, 0.383) = 0.000103; every capped run's min net ≥ its
base margin, as 008's marginal-slack argument predicts). Meanwhile out-of-class slice
(n = 7, marginals 0.4286) is negative even in *actual* gain at ρ = 1.045
and 1.06 (−0.034, −0.030 at δ = 0.1) — the cap is load-bearing exactly as
008 says. "Every in-class run is positive" survives a strictly harsher
test than 008 ran (R5's mislabeled sentence aside). CONFIRMED.

### The pointwise budget refutations (claim 5, parts that survive)
- Naive n-uniform two-sided pointwise budget REFUTED — confirmed and
  extended: max|log₂OR−λ| at δ = 0.05 reaches 8.3 (crash, (0.30,1), n=28)
  and 1.80 ((0.383,1.2), n=28); upward-only ✓ (OR < 1 mass identically 0
  in every run of every part, both engines).
- Max-coordinate downward ≈ 0.0827(n−2)δ²: confirmed at n ≤ 9 to all
  printed digits and extended to n = 28 — with the correction that the
  slope *drifts up* (increments 0.0825 → 0.0882 by n = 16, reaching 2.36
  vs the linear fit's 2.16 at n = 28): "linear" is the right first-order
  reading, slightly optimistic in level. Always at i = 2 ✓.
- (B4) vacuity in the measured range ✓ (max upward τ ≈ 1.8 < 3.0 at
  n ≤ 28 — though R1 makes the τ_half *definition* box-dependent, so (B4)
  needs restating on the repaired box).
- (B3) shape (n-saturating margin-variance) ✓, constant CORRECTED to
  ≈ 4.7, not 2.5 (R2).

## Verdict

| # | 008 claim | Verdict |
|---|-----------|---------|
| 1 | Lemmas A–C + Theorem P6′ (the Prop-6 closure) | **CONFIRMED** (hand re-derivation of every step; 120-instance stress grid incl. δ = 0.5, 0 violations; interior-δ IFT check; λ ≥ 0 scope correctly load-bearing). Constant chase **CORRECTED**: stated derivation yields 831; printed 924 valid but corresponds to m₀ ≤ 1 (R6) |
| 2 | `924·min(p,1−p)^{−3n}δ²` bound, wildly loose | **CONFIRMED** (loosest measured gap 4.7e-8 of bound); looseness illustration **CORRECTED** (own census's worst c is 20.6, not 0.6) (R6) |
| 3 | Leads 3/4 distinct; centered-kernel crash OR = 2^{λ(3−n)/n} | **CONFIRMED** (hand + potential-free identity + independent Sinkhorn, exact; mirror and product band too) |
| 4 | ρ*(0.383) = 1.0422, queue's 1.03 miscalibrated, g(1.031, .383) < 0 | **CONFIRMED** (queue framing was wrong; 008's recalibration necessary) |
| 5a | Pointwise budget not n-uniform, +0.3 bits/coord, upward-only | **CONFIRMED**, extended to n = 28 |
| 5b | Max-coordinate downward ≈ 0.083(n−2)δ² | **CONFIRMED** to n = 28 with upward slope drift noted |
| 5c | Averaged downward budget flattens (→ ≈ 0.07δ² / ≈ 0.3δ²) | **REFUTED** (orbit census to n = 32, δ²-coefficient verified stable in δ at n = 32: growth resumes with increasing increments in both regimes; claimed limits passed at n ≈ 16 / n = 10; (B1) as stated dead at EVIDENCE level by 008's own lead-1 criterion) (R3) |
| 5d | Tax budget tiny, `tax/δ² ≲ 0.01`, shape undetermined | **CORRECTED**: shape now determined and it is not δ² — tax is first-order in δ in every observable window at large n; fixed-δ coefficient exceeds 0.01 by 10× at n = 32; absolute tax still tiny, so no assembly threat (R4) |
| 6 | Conditional theorem (B1)–(B4) ⟹ (S-coup) on δ₀-ball; δ₀ ≈ 0.022 | Structure **CONFIRMED** after repairs; box-uniform τ_half step **REFUTED as stated** (analytic corner counterexample at 008's own s₀ = 0.10; s₀ ceiling 0.0843; claimed numerical check absent from the code) (R1); δ₀ **CORRECTED** to ≤ 0.0079 (008's own formula) / ≈ 0.004 (corrected constants) (R2) |
| 7 | 78 runs all positive in class; only negatives cap-violating | **CONFIRMED** (spot ≤ 5e-6; 104-run hunt incl. new half-empty direction, 0 capped negatives); "worst in-class net" sentence **CORRECTED** (that run is out of class) (R5) |

**Net assessment.** 008's proof deliverable — the closure of 005 Prop 6's
smoothness gap — is real and survives a line-by-line hostile pass; it
should stand as VERIFIED (fixed n, λ ≥ 0). The lead-3 execution and the
cap story are genuine and survive harder adversaries than 008 ran. What
does not survive is the record's *forward-looking* reading: the
conditional theorem's instantiated constants are wrong (R1, R2), and the
two budget hypotheses it depends on go the wrong way once n leaves 008's
window — (B1)'s coefficient resumes growth by n ≈ 22 and (B2) is not even
quadratic in δ (R3, R4). The n-uniform perturbative assembly is thus in
strictly worse shape than 008 reports: gap a′'s perturbative form now has
direct counter-evidence to n = 32, and any future budget statement must be
first-order in δ for the tax and super-logarithmic-tolerant in n for the
downward OR term. The positive residue is unchanged and real: end-to-end
positivity on every capped instance tested by anyone (δ ≤ 0.2, n ≤ 9,
ρ ∈ {1.045, 1.06, 1.20}) — whatever certifies it will not be budgets
(B1)–(B4) as stated.

## Residual risk

- The census growth findings (R3, R4) are crash-direction, fixed
  (p, λ)-pair EVIDENCE to n = 32; a bounded c₁(n) with onset beyond 32, or
  a different behavior on other in-class directions, is not excluded. The
  orbit engine only handles S_{n−2}-symmetric directions; asymmetric
  adversaries at n > 9 remain unexplored.
- My δ → 0 stability check (SK7) covers n = 32 down to δ = 0.003125;
  the (B1) coefficient there moves by < 3% per δ-halving (down), while the
  tax has demonstrably not converged — conclusions about (B2) rest on the
  ratio trend (→ 2), not on a reached limit.
- Both engines share the atoms-as-bitmasks / cells-by-prefix census
  *definition* with 005/006/008; 006 S1 checked that definition against
  the slice/inner-product formula, and SK0's orbit-vs-direct agreement is
  a genuinely different decomposition of the same object, but no third
  formalization of "conditional OR given histories" exists yet.
- The repaired conditional theorem (s₀ ≤ 0.0843, re-measured constants,
  first-order (B2)) is sketched here, not written out; nobody has
  re-proved it end to end with the new box.

## References

- `problems/union-closed/attempts/008-perturbative-assembly.md` (under
  review); `005-odds-ratio-control-refuted.md` (Prop 6, leads 3–4, crash
  family); `006-skeptic-review-of-005.md` (S1 reduction, S7
  family-dependent onset, S8); `003-dependent-couplings.md` /
  `004-skeptic-review-of-003.md` (definitions, ρ*, g, capped class);
  STATUS.md (queue line for claim 4).
- Reviewed artifacts: `explore/uc_pert.py`; `data/pert_P0_validation.json`
  … `pert_P4_twoscale.json`; `data/pert_full_run.log`.
- This review's tool and data: `explore/uc_pert_skeptic.py`;
  `data/pertsk_SK0.json` … `pertsk_SK7.json`.
- Gilmer arXiv:2211.09055; Alweiss–Huang–Sellke arXiv:2211.11731; Liu
  arXiv:2306.08824 — context only, all used via 003's closed forms as in
  008. No external fetches; the IFT material is checked inline against
  008's from-scratch proof, not against a citation.
