# 014 — Large-n probe kills the i-aggregated odds-ratio control

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed
- **Type:** computational search (falsifiable probe of the restated Gap 1 at
  n ≫ 7, run BEFORE proof effort per queue item 1); adversarial families +
  exact rational certification. Outcome: the aggregated control is **REFUTED
  in the record-relevant regime** by an explicit unit-replicated family at
  n ≈ 96–128, certified float-free.
- **Tools:** `explore/uc_or_agg_probe.py` (round 1: anchors, ladder trend to
  n = 32, seeded free-support climbs, permutation attack, exact certs;
  ~3 min) and `explore/uc_or_agg_probe2.py` (round 2: ladder extension to
  n = 128, unit-symmetric θ-optimization, deep annealed climbs, adversarial
  dilution sweep `--sweep`, kill drill `--kill`; ~30 min + drill).
  Both standard-library only, deterministic (fixed seeds, fixed step counts).
  Float census is 007's engine (`uc_or_avg.py`) imported verbatim; exact
  rational census + certified log₂ enclosures are 013's engine
  (`uc_or_avg_skeptic.py`) imported verbatim; a second float engine
  (013's own `f_census`) is used as an independent cross-check on the
  violating instance. Commands reproducing every number:

      python problems/union-closed/explore/uc_or_agg_probe.py 2>&1 | tee problems/union-closed/data/aggprobe_run.log
      python problems/union-closed/explore/uc_or_agg_probe2.py 2>&1 | tee problems/union-closed/data/aggprobe2_run.log
      python problems/union-closed/explore/uc_or_agg_probe2.py --sweep 2>&1 | tee problems/union-closed/data/aggprobe2_sweep.log
      python problems/union-closed/explore/uc_or_agg_probe2.py --kill 2>&1 | tee problems/union-closed/data/aggprobe2_kill.log
      python problems/union-closed/explore/uc_or_agg_probe2.py --cert96 2>&1 | tee problems/union-closed/data/aggprobe2_n96cert.log
      python problems/union-closed/explore/uc_or_agg_probe2.py --certp3   # re-cert of P3-best, writes aggprobe2_p3best_cert.json

  Checkpoints: `data/aggprobe_part[ABDE]*.json`, `data/aggprobe_partC_n*.json`,
  `data/aggprobe2_partP[1-7]*.json`.
- **Sources:** 007/013 (definition of M_i, the aggregate, the 10-atom
  witness, exact-arithmetic standard), 012 (large-n reversal warning, orbit
  symmetrization as the way to buy n), 005/006 (families), STATUS.md queue
  item 1. No external fetches this cycle.

## Approach

**The statement under test** (007 lead 1 = 013 claim 6, the restated Gap 1
`aggregated-or-control`; stated here exactly as certified positive in 013
part C). For μ on 2^[n] with H(μ) > 0 and λ > 0, with M_i(μ, λ) the
well-posed history-averaged control of 007 §1 (mass-weighted mean of
log₂ OR_i(a,b) over the nondegenerate history block N_i × N_i of the tilt
coupling π_λ ∝ u(A)u(B)2^{λ|A∩B|}) and w_i = the defined history mass at
coordinate i:

    A(μ, λ)  =  Σ_{i<n, N_i≠∅} w_i · (M_i − λ)  /  Σ_{i<n, N_i≠∅} w_i   ≥  0 ?

(i = n is excluded; M_n = λ identically by 005 Prop 1, so it contributes 0.)
**Kill condition:** an instance with A(μ, λ) < 0 certified in exact rational
arithmetic at rational tilt t (the aggregate numerator is then
Σ Σ m·log₂(OR/t), every OR/t an exact rational, so the sign is a theorem
about the stated inputs), with all elementwise marginals < 38271/100000.
Prior state: certified positive on the 10-atom witness (+1.844669 at
t = 181/16, 013) and unrefuted by 007 part U (lowest endpoint +0.0125) and
013 part G (21 witness-seeded endpoints, lowest +0.472) — all at n ≤ 8.

**Why a large-n probe rather than more n ≤ 8 search or proof effort.** Queue
item 1's own instruction, grounded in 012's lesson: the budget census
"flattening" reversed past n ≈ 22, so small-n survival is not evidence of
n-uniformity. The aggregate has a specific large-n vulnerability that no
n ≤ 8 search can see: **unit replication**. The 007 witness's deficit lives
at ONE response coordinate propped open by a light slice; the aggregate
survives there because the O(1) remaining coordinates carry large surpluses.
But the two-block geometry can be replicated across r = Θ(n) response
coordinates sharing one block marker, two future coordinates and three
dilution coordinates, so the deficit budget grows with n while the surplus
budget stays O(1) coordinates wide. The probe therefore (a) measures the
aggregate along that ladder as n grows, (b) re-optimizes the ladder's unit
weights (an 8-parameter, unit-symmetric search — the orbit-symmetrization
idea in search form, which is how a Θ(n)-atom family stays searchable at
n = 128), and (c) keeps the free-support and permutation attacks running as
controls.

## What was done

### 1. The multi-unit ladder MU(n, r) (the killer genre, scaled)

Coordinate 1 = block marker; coordinates 2–4 = dilution singletons (weight
shared, tuned); response coordinates 5..4+r (one per unit); shared future
coordinates n−1, n. Block a: {1, n−1}, {1, n}, plus per unit j the light
slice {1, c_j, n}; block b: ∅, {n−1}, {n}, plus per unit j the responder
{c_j, n−1}. Weights per unit copied from the 007 witness (so **MU(7,1) is
exactly the 10-atom witness**). 2r + 8 atoms, all marginals in-regime after
dilution tuning.

### 2. Round 1 (`uc_or_agg_probe.py`): anchors, trend to n = 32, controls

- **Anchors (part A).** Witness aggregate reproduces: float +1.844754 at
  λ = 3.5, and the exact engine re-certifies **+1.844669** at t = 181/16
  (013's number to all printed digits). Padding invariance: MU(20,1) and
  MU(32,1) equal the witness aggregate to 0.00e+00 — unused coordinates
  carry w_i = 0, an exact engine check. Crash/mirror embedded at
  n ∈ {8,16,24,32} (≤4-atom supports are Theorem-B-safe): 16 instances, min
  aggregate +0.71. Degeneracy dichotomy flag clean everywhere.
- **Ladder trend (part B).** With witness weights, aggregate falls
  monotonically in n at every λ: e.g. λ = 3.5: +2.60 (n=7) → +1.28 (n=32);
  λ = 2: +1.35 → +0.81. Per-unit deficit constant (per-i min −0.1220 at
  λ = 3.5 at every n). Best 2-parameter fit over n = 7..32 is a − b·log n
  (max residual 0.033–0.051 vs 0.06–0.20 for a + b/n and a + b/√n),
  extrapolating to a sign change at n ≈ 320 (λ=2) / 140 (λ=3.5) / 78 (λ=5).
  Extrapolation is not evidence; it set round 2's target range.
- **Seeded free-support climbs (part C).** 7–8 climbs × 300–500 steps per
  n ∈ {8,12,16,20,24,28,32}; seeds = MU ladders, witness-spread, sparsified
  Sawin-genre mixtures, δ_∅⊕Bern genre, random sparse; mutations = weight
  jiggle, light-slice push, bit toggle, atom birth/death; dynamic-range
  clamp 1e-9 per 007's underflow guard. **0 in-regime violations**; best
  in-regime endpoints per n: +0.166 (n=8, witness-spread genre) and
  +0.26..+0.35 for n = 12..32 (all MU-seeded at λ = 2). Degenerate
  out-of-regime collapses (max marginal → 1, H → 0.5) reach +0.004 but are
  not admissible.
- **Permutation attack (part D).** Coordinate order is part of μ, and the
  aggregate is not permutation-invariant. 41 orders per target: witness
  +1.84 → min +1.69; MU(13/22/32) minima stay ≥ +1.83; but the optimized
  n = 8 endpoint drops +0.166 → **+0.0183** (9× margin cut, still
  positive). Order optimization is a real attack channel the earlier
  records never used.
- **Exact certs (part E).** Witness +1.844669005 (enclosure width < 1e-40);
  MU(22)=MU(32) at r=14, +1.730725429; MU(32) at t = 4, +1.018477406; the
  n = 8 search endpoint +0.166100824. All CERTIFIED POSITIVE, in-regime,
  dichotomy exact.

### 3. Round 2 (`uc_or_agg_probe2.py`): n up to 128 and the kill

- **P1/P6, honest ladder trend.** The dilution weight is the adversary's
  parameter too, so the per-n ladder margin is min over the dilution weight
  subject to in-regime marginals (the P1/part-B bisection tuners land on
  different feasible branches at different n — visible as jumps at
  n = 64/80 — so P6's 16-point sweep per (n, λ) is the honest curve).
  With the original witness weights the honest curve keeps falling but stays
  positive through n = 128: λ = 2: +0.784 (n=10) → +0.281 (n=32) → +0.157
  (n=96), +0.220 at n=128 (feasibility-window artifacts of the 16-point
  dilution grid at n = 48 and 128 make the curve locally non-monotone);
  λ = 3.5 similar (+0.185 at n=96). So the RAW witness ladder does not kill
  by n = 128 — re-optimizing the unit weights is what does.
- **P2, unit-symmetric optimization.** The ladder is determined by 8 shared
  weights θ = (A1, A2, light, B0, B1, B2, resp, dil). Climbing θ at n = 48
  (260 steps, marginal penalty) gives +0.0342 at λ = 2 (θ*: A1 14.0,
  A2 22.6, B0 13.6, B1 43.0, B2 24.6, dil 21.0, light 0.0079,
  **resp 0.0017** — the optimizer made block b's responder even lighter
  than the light slice, shrinking the block-coordinate surplus while
  keeping every unit's deficit at −0.059). Transferring θ* to larger n
  (dilution re-tuned only):

      λ = 2.0:  n = 64: +0.0285   n = 96: −0.0008   n = 128: −0.0154
      λ = 3.5:  n = 64: +0.0443   n = 96: −0.0069   (n = 128 tuner branch-jumped to +1.09)

  **First in-regime float-negative aggregates**, max marginal 0.309–0.365.
- **P3, deep annealed free-support climbs** at n = 24, 32 (2500 steps,
  round-1 best endpoints + fresh ladders as seeds, occasional uphill
  acceptance): no improvement — 6 endpoints, 0 violations, best in-regime
  +0.2803 (n = 32). Free-atom search does not find the replication channel
  at n ≤ 32 even when seeded next to it; consistent with the crossing
  genuinely needing n ≈ 90.
- **P4, λ-profiles.** The violating θ* ladder at n = 96 is negative on a
  λ-window, not a knife edge: +0.015..+0.025 for λ ≤ 1.5, **−0.001 at
  λ = 2.0, −0.006 at λ = 2.5**, positive again from λ = 3 (where the
  dilution re-tune jumps branch). The P3-best n = 32 instance is positive
  at every λ in [0.25, 6]. The window matters downstream: the assembly
  chooses λ per instance, and the family's violating window sits exactly
  at moderate λ where 009/011's λ-window law puts the workable range.
- **P5/P7, certification of the kill** (exact rational, t = 4 so λ = 2
  exactly; enclosure widths ≤ 1e-18; degeneracy dichotomy re-checked
  exactly at every coordinate; marginals certified < 38271/100000):

      n = 96, 188 atoms, float-exact dyadic weights:
        A ∈ [−0.000759865185844166087, −0.000759865185844166086]   KILL (142s)
      n = 128, 253 atoms, tidy rational weights (limit_denominator 1e4):
        A ∈ [−0.015405727, −0.015405727]                            KILL (237s)
      n = 160, 316 atoms, float-exact dyadic weights:
        A ∈ [−0.024225909, −0.024225909]                            KILL (427s)
      n = 160, tidy rational weights:
        A ∈ [−0.024225909, −0.024225909]                            KILL (404s)

  Positive controls certified the same way: witness +1.844669005
  (= 013 part C), MU(22/32, r=14) +1.730725429 at t = 181/16, MU(32) at
  t = 4 +1.018477406, MU(64) tidy +0.180142390, round-1 n = 8 endpoint
  +0.166100824, P3-best n = 32 +0.297679573. One bookkeeping artifact
  caught and fixed: `limit_denominator` rounded one 1e-9-weight atom of
  the P3-best instance to exactly 0, leaving a phantom zero-weight atom
  that tripped the exact dichotomy flag on the first P5 pass; weights are
  now zero-dropped after rationalization and the re-certified value
  (`data/aggprobe2_p3best_cert.json`) is unchanged at +0.297679573 with
  the dichotomy exact. No kill certificate involved weights small enough
  to be affected (θ* weights ≥ 0.0017).

### 4. Validity checks on the violating family

- **Independent engine.** The n = 96 instance evaluates to −0.0007598652
  identically (diff 1.5e-13) on 007's census and on 013's own from-scratch
  float engine (`f_census`); both share IEEE arithmetic, which is why the
  exact certification below is the actual judge (verification standard,
  2026-07-31).
- **Not a degenerate corner.** H(μ) ≈ 3.0 bits; defined mass healthy at
  Θ(n) coordinates (Σ w_i ≈ 47.9 at n = 96, i.e. ≈ 0.50 of history mass
  per coordinate); max elementwise marginal 0.307–0.309, well inside the
  record regime; cell floors far from underflow (weights span 1.7e-3..43,
  cell-floor diagnostic recorded in `aggprobe2_partP7.json`); exact
  degeneracy-dichotomy check clean at every coordinate.
- **Perturbation stability.** 20/20 independent 3% multiplicative
  perturbations of ALL 316 weights at n = 160 remain in-regime violations
  (range −0.0253 .. −0.0208) — a stable open set, not a tuning accident.
  The tidy-rational (4-significant-digit) version certifies at the same
  −0.024225909, the analogue of 007's 2-digit tidy witness.
- **Convention robustness (proved, one line).** Under 013's
  conjecture-friendliest alternative bookkeeping — score every degenerate
  history at exactly λ instead of conditioning it out — the aggregate's
  numerator is unchanged (degenerate histories contribute M − λ = 0) and
  its denominator can only grow, so the SIGN of every violation here is
  invariant under the bookkeeping choice. No convention rescues the
  statement on these μ.
- **λ is not tuned to a knife edge, but it is windowed.** At n = 96 the
  θ* family is negative on a λ-window around [2, 2.5] (P4 grid; −0.006 at
  λ = 2.5 is deeper than the λ = 2.0 optimum was tuned for) and positive
  outside it. The certified kills are at λ = 2 exactly (t = 4). A
  violation for every λ is NOT claimed and is false for this family —
  the refuted statement is the ∀μ ∀λ control, which one (μ, λ) pair
  refutes, exactly as in 007/013's per-i kill.

## Outcome

**REFUTED — the i-aggregated odds-ratio control (`aggregated-or-control`,
the restated Gap 1) is false, including in the record-relevant regime.**
The probe's kill condition is met: the unit-replicated two-block family
with re-optimized shared weights θ* violates A(μ, λ) ≥ 0 for every tested
n ≥ 96, and the violation is **certified float-free in exact rational
arithmetic at the exactly-rational tilt t = 4 (λ = 2)**:

    n = 96  (188 atoms):  A ∈ [−0.00075986518584416609, −0.00075986518584416608]
                          max marginal 0.30867 < 0.38271     (data/aggprobe2_n96cert.json)
    n = 128 (253 atoms):  A ∈ [−0.015405727, −0.015405727]
                          max marginal 0.30793               (data/aggprobe2_partP5.json)
    n = 160 (316 atoms):  A ∈ [−0.024225909, −0.024225909]   (dyadic AND tidy weights)
                          max marginal 0.30719               (data/aggprobe2_partP7.json)

both with the degeneracy dichotomy checked exactly at every coordinate, so
each sign is a theorem about the stated rational inputs (013's standard).
At t = 4 no irrational number appears anywhere in the statement. Scope of
the refutation and of the surviving positive space:

- **Violations found:** n ≥ 88–96 (float crossing between n = 88 and 96
  for θ* at λ = 2; float-negative also at λ = 3.5, n = 96), monotone
  deepening to −0.0242 at n = 160 with no sign of reversal.
- **Searched without violation (EVIDENCE, scoped):** n ≤ 32 free-support
  climbs (58 seeded climbs across n ∈ {8,...,32}, λ ∈ {2, 3.5, 5}, seeds =
  ladder/witness-spread/Sawin-sparse/δ_∅⊕Bern/random, 300–500 steps each,
  round 1) + deep annealed climbs at n ∈ {24, 32} (2500 steps, round 2);
  permutation orbits (41 orders × 5 targets); the RAW-witness-weight
  ladder for all n ≤ 128 at any dilution weight (P6 sweep, min +0.157 at
  λ = 2, n = 96); ≤4-atom families at n ≤ 32 (Theorem-B-safe anchors).
  λ grid {0.5, 2, 3.5, 5} in searches; exact certificates at t = 181/16
  and t = 4.

**Not claimed:** no claim that the dependent-couplings route dies — as
with 005/007 this kills a candidate bridge, not the interface (licensing
lemma, separations, 0.4315 ceiling untouched); no claim about the
margin-modulated control (007 lead 2 — untested against this family, now
the leading Gap-1 candidate; lead 1 below); no claim that a fixed-small-n
aggregated control fails (everything searched at n ≤ 32 is positive; the
minimal violating n for THIS family sits in (88, 96] but lead 2's direct
optimization could push lower); no claim about μ-level MTP₂ (still open).
Per repo rules nothing here is a result until an independent skeptic pass;
the exact certificates are the part to attack first (they are two-line
claims about explicit rational data), then the aggregate definition's
match to 007/013 (anchored by the +1.844669 reproduction), then the
search-completeness claims, which are weaker and scoped as EVIDENCE.

**Ledger consequence (for STATUS.md, pending skeptic):** queue item 1's
"if it survives n ≳ 20, attack the proof" branch is CLOSED NEGATIVE — do
not invest in the Gram/Frobenius proof effort for the aggregate as stated;
the live Gap-1 candidates are now the margin-modulated control and
fixed-n/decaying-λ restatements. The `aggregated-or-control` gap tag dies
with this record.

## Why it failed / what survived

**The quantity that goes the wrong way: the surplus is O(1) coordinates
wide, the deficit is Θ(n) coordinates wide.** In the ladder, every unit
contributes a deficit coordinate whose defined mass does not decay with the
number of units (the per-i min is exactly the same at every unit — the
units are exchangeable given the block), while every surplus term the
aggregate can spend lives at the O(1) structural coordinates (block marker
i = 1 with w_1 = 1, three dilution coordinates, the final future
coordinate). The witness weights leave each unit's deficit small relative
to those surpluses, so the raw ladder decays only slowly (still +0.16 at
n = 96); but the unit weights are the adversary's too, and re-optimizing
them (8 shared parameters) shrinks the block-coordinate surplus faster
than the per-unit deficit, moving the crossing into reach: the aggregate
is then a monotone decreasing function of n that crosses zero between
n = 88 and n = 96 and keeps falling (−0.024 at n = 160, still falling).
This is invisible at n ≤ 8 by construction — with one unit the surplus
dominates ~15:1 — which is exactly the small-n-optimism failure mode 012
warned about, now realized on the aggregate itself. 007's own mechanism
note ("the chain-rule assembly sums over i anyway") identified the sum
over i as the reason to trust the aggregate; the sum over i is precisely
what the replication attack weaponizes.

**Why the earlier searches missed it:** 007 part U and 013 part G searched
free supports at n ≤ 8 (where no replication exists) and hill-climbed the
witness's own atoms (which can re-weight but not replicate units). Round
1's free-support climbs at n ≤ 32 also found nothing — the basin around
the replicated geometry is apparently narrow in free-atom space and the
crossing needs n ≈ 90 even for optimized weights; the kill required
searching the unit-symmetric quotient (8 parameters at any n), which is
the searchable-at-scale parametrization 012's orbit engine prefigured.

Survived / reusable:

- **The n ≤ 32 positive space stands.** Nothing searched at n ≤ 32 —
  free supports, all recorded genres, permutation orbits, ladders with any
  dilution weight — violates; certified positives at n ≤ 32 include the
  witness (+1.844669) and ladder points (+1.7307 at t = 181/16, +1.0185
  at t = 4). A FIXED-small-n aggregated control is untouched by this
  attempt; what died is the n-uniform statement the route needed.
- **The MU(n, r) ladder family + unit-symmetric θ parametrization** — the
  reusable adversary for ANY future control statement on this route: any
  candidate must now survive unit replication with adversarial unit
  weights, not just the 10-atom witness.
- **The permutation attack** (part D): coordinate order is part of μ; it
  cut an optimized margin 9× at n = 8. Cheap, never used before this
  record; should join the standard battery.
- **Engines:** `uc_or_agg_probe.py` / `uc_or_agg_probe2.py` re-use 007's
  census and 013's exact machinery unchanged and add: aggregate evaluation,
  ladder builders, θ-space (orbit-quotient) search, dilution sweeps,
  n-scans, and exact aggregate enclosures at rational tilt on ~300-atom
  instances in minutes.
- **The padding-invariance anchor** (MU(n,1) ≡ witness at every n,
  diff 0.0) — a free exactness check for any future engine touching
  this definition.

## Leads generated

1. **Test the margin-modulated control against the ladder (007 lead 2, now
   urgent).** Compute the h-sensitivity-weighted aggregate
   `E[σ(x̃,ỹ)·(log₂OR − λ)]` on θ*-ladder instances across n = 48..160.
   Definite outcome: if it stays positive while the plain aggregate
   crosses, the margin-modulated statement graduates to the sole Gap-1
   candidate and inherits the probe machinery; if it crosses too, the
   whole OR-control corridor of the route is dead and Gap 1 must be
   re-stated at the functional level (007 lead 4).
2. **Find the minimal certified-violation n.** Optimize θ directly at the
   target n (not transferred from n = 48) with the crossing-n as the
   objective; certify the smallest n with a rational-tilt violation. The
   transferred θ* crosses at n ≈ 96; direct optimization plausibly pushes
   below n = 64, which would make skeptic re-certification cheap. Definite
   outcome either way (a certified n*, or evidence the crossing resists
   optimization below some n).
3. **Does the raw-witness-weight ladder cross?** P6's honest
   (dilution-minimized) trend still decreases at n = 96 (+0.157 at λ = 2)
   but was only run to n = 128 with a 16-point dilution grid. Extend to
   n ≈ 320 with a finer grid: crossing or asymptote is a definite outcome,
   and calibrates how much of the kill is re-weighting vs pure replication.
4. **λ-window law of the violation.** At n = 96 the θ* family violates on
   a window around λ ∈ [2, 2.5] and not outside it (P4). Measure the
   window edges λ_lo(n), λ_hi(n) along the n-scan and check them against
   009/011's assembly window λ_max ≈ 4.847/(n−3): if the violating window
   and the assembly's workable window separate as n grows, a
   λ-RESTRICTED aggregated control (λ ≤ c/n, the regime the assembly
   actually uses at large n per the λ-window law) may survive the
   replication attack — a precise, falsifiable restatement candidate the
   next cycle can test with this record's engines unchanged.
5. **Consequence for 008/012's assembly budgets.** The δ-linear tax and
   the corrected budgets (queue item 2) were to be paired with an
   aggregated OR control; that pairing is now dead as stated. Re-derive
   what budget object the conditional theorem needs if the OR control is
   margin-modulated (lead 1) — a precise restatement is a definite
   deliverable and belongs with queue item 2's rebuild.

## References

- This repo: `attempts/007-averaged-or-control.md` (M_i definition §1, the
  witness, lead 1 = the target statement); `attempts/013-skeptic-review-of-007.md`
  (exact-arithmetic standard, part C aggregate certification, part G seeded
  attacks); `attempts/012-skeptic-review-of-008.md` (large-n reversal
  lesson, orbit-symmetrized engine precedent); `attempts/005/006` (families
  and census conventions); `attempts/008` (assembly context for the
  aggregate's downstream role).
- Tools/data: `explore/uc_or_agg_probe.py`, `explore/uc_or_agg_probe2.py`;
  checkpoints `data/aggprobe_*.json`, `data/aggprobe2_*.json`; logs
  `data/aggprobe_run.log`, `data/aggprobe2_run.log`,
  `data/aggprobe2_sweep.log`, `data/aggprobe2_kill.log`.
- No external papers consulted this cycle (record threshold 0.38271 per
  prior records, Liu).
