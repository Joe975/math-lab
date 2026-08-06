# 018 — The two surviving Gap-1 candidates vs the ladder adversary

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-06
- **Mode:** informed
- **Type:** computational search + candidate-statement triage (queue items 1
  and 20; 016 leads 1 and 4). Outcome: the **margin-modulated control is
  REFUTED in both signed readings at n = 7** — by 007's own witness, the
  instance it was invented to survive — certified float-free; the |σ|-weighted
  unsigned reading and the **λ-window-restricted control survive** the full
  ladder adversary with re-optimized weights (`EVIDENCE`, scoped below).
- **Tools:** `explore/uc_gap1_candidates.py` (written here; standard-library
  only, deterministic — fixed seeds, fixed step counts). Parts A–D (default,
  ~35 min), `--ext` (part E, ~25 min), `--cert` (part F, exact rational
  certification), `--robust` (part R, kill-robustness battery), `--fast`
  (smoke). Float census logic is a fresh implementation cross-anchored
  against 007's engine (`uc_or_avg.py`, imported for the plain aggregate);
  exact machinery is 013's (`uc_or_avg_skeptic.py`) plus a new certified
  interval path for the margin-modulated quantity (dyadic bisection for the
  Plackett root, directed log₂ enclosures for the binary entropy). Commands
  reproducing every number:

      python problems/union-closed/explore/uc_gap1_candidates.py           2>&1 | tee problems/union-closed/data/gap1c_run.log
      python problems/union-closed/explore/uc_gap1_candidates.py --ext     2>&1 | tee problems/union-closed/data/gap1c_ext.log
      python problems/union-closed/explore/uc_gap1_candidates.py --robust  2>&1 | tee problems/union-closed/data/gap1c_robust.log
      python problems/union-closed/explore/uc_gap1_candidates.py --cert    2>&1 | tee problems/union-closed/data/gap1c_cert.log

  Checkpoints: `data/gap1c_part[A-F,R].json`.
- **Sources:** 007 (lead 2 — the margin-modulated candidate; §1 M_i
  definition; the 10-atom witness), 013 (exact-arithmetic standard, aggregate
  certification), 016/017 (the MU(n,r) ladder, θ\*, the λ-window
  disjointness correction C2), 009/011 (the λ-window law
  λ_max ≈ 4.847/(n−3)), 008 (h(z_ρ) gain machinery). No external fetches
  this cycle.

## Approach

**Context.** 016/017 killed the i-aggregated OR control at n ≥ 96 and left
exactly two Gap-1 candidates standing: the margin-modulated control (007
lead 2) and the λ-window-restricted control (017 C2). The queue's standing
policy is probe-before-proof, and 016 supplies the one adversary genre known
to kill at scale: the unit-replicated ladder MU(n, r) with adversarially
re-optimized shared weights θ. This attempt runs that adversary — with θ
re-optimized against each NEW statement, since a survival against weights
tuned for a different functional would be worthless — against both
candidates, plus free-support and (θ, λ)-joint attacks.

**Making candidate I precise (the fork in 007 lead 2).** "E[σ(x̃,ỹ)·(log₂OR
− λ)] ≥ 0 with σ the Plackett sensitivity" admits three readings, which this
attempt formalizes and tests separately. At a nondegenerate history (a, b)
of coordinate i, write the conditional table's zero-margins x = P(A_i = 0),
y = P(B_i = 0), realized both-zero probability z̃, and let z_ρ(x, y) be the
Plackett both-zero probability at odds ratio ρ (the unique root of
z(1−x−y+z) = ρ(x−z)(y−z) in the coupling range). Then, mass-weighted over
all nondegenerate histories of all i < n:

- **MM_sec** (secant / assembly-exact): E[h₂(z̃) − h₂(z_{2^λ}(x, y))]. By
  the mean-value theorem this IS E[σ·(log₂OR − λ)] with σ the sensitivity
  ∂h₂(z_ρ)/∂λ evaluated on the secant — and h₂(z̃) = h₂(z_OR(x, y)) exactly,
  since a 2×2 table is determined by its margins and odds ratio. This is
  the reading with chain-rule meaning: it compares each history's realized
  Plackett gain against the gain the λ-target would deliver.
- **MM_der** (signed derivative at target): E[σ_λ·(log₂OR − λ)] with
  σ_λ = ∂h₂(z_ρ(x, y))/∂λ at ρ = 2^λ.
- **MM_abs** (unsigned weight): E[|σ_λ|·(log₂OR − λ)], normalized by
  E[|σ_λ|] — a weighted mean of (log₂OR − λ) with a nonnegative weight that
  vanishes as the margins degenerate, which is the literal "downweight the
  witness's hiding place" reading.

The distinction matters because σ_λ changes sign at z = 1/2, and in-regime
margins straddle that point (z_1(x, x) = x² crosses 1/2 exactly at the
(3−√5)/2 barrier — the sensitivity sign boundary is the record threshold's
golden-ratio structure again).

**Candidate II** is the plain aggregate A(μ, λ) of 014/016 restricted to
λ ≤ λ_win(n) = 4.847/(n−3) (the 009/011 window law — the regime the
assembly can actually use at large n). 016's kill lives at λ ∈ [2, 2.5],
~40× above the window at n = 96 (017 C2), so this candidate was genuinely
untested. The probe measures the ladder INSIDE the window with θ
re-optimized at window λ, lets the adversary choose λ within the window,
and fits the small-λ expansion A ≈ a₁(n)λ + a₂(n)λ² — the structural
question being whether a₁(n) (Theorem C's perfect-square first-order
coefficient) decays fast enough for the quadratic term to cross inside a
window that itself shrinks like 1/n.

**Why this rather than proof effort on either candidate:** 016 is the
second demonstration (after 012) that small-n survival licenses nothing;
both candidates had zero adversarial exposure in their stated regimes.

## What was done

### A. Anchors (all pass; `gap1c_partA.json`)

Witness plain aggregate reproduces 013 (+1.844754 float at λ = 3.5, and
+1.844669005 exactly in the part-F certification); the fresh census's
internal plain aggregate matches 007's engine to 2.2e-16; 016's θ\* kill
instance reproduces (−0.000759865 at n = 96, λ = 2); the Plackett
round-trip z_ρ(x, y, OR) = z̃ holds to 1.6e-15 on every witness row (the
2×2-table-determination fact the secant identity rests on); secant→derivative
consistency at Δλ = 1e-3 (rel 2e-4); MM padding invariance MU(20,1) ≡
witness, diff 0.0 exactly. Product-measure anchor (part R): MM_sec vanishes
identically on Bern(p)^n (|MM_sec| < 3e-16; OR = 2^λ at every history so
z̃ = z_target row by row — an exact structural zero the float path
reproduces).

### B. Candidate I on the witness and the ladders (`gap1c_partB.json`)

The signed readings fail immediately, and not at large n:

    007 witness (n = 7, λ = 3.5, marginals ≤ 0.318):
      plain aggregate +1.8448   (the value 013 certified positive)
      MM_sec  −0.022255         VIOLATION
      MM_der  −0.895542         VIOLATION
      MM_abs  +2.354328         survives

λ-profile on the witness (part R): MM_sec < 0 for every λ ∈ [1, 5]
(minimum −0.0224 near λ = 4), positive only for λ ≲ 0.7 — exactly the
first-order-dominated region where Theorem C forces positivity. MM_der
flips negative from λ ≈ 1.5. So the margin modulation does not rescue the
control on the very instance that killed the unmodulated version — it makes
it fail *harder* (the plain aggregate is +1.84 there; the modulated ones
are negative).

On the ladders at λ = 2 (raw witness weights and 016's θ\*), MM_sec is
negative at every n ∈ {16, …, 160} (θ\*: −2.4e-3 at n = 16 shrinking to
−2.2e-4 at n = 160; the replication that kills the plain aggregate
*shrinks* the MM_sec deficit), MM_der negative throughout, MM_abs positive
throughout (≈ +0.97 on θ\*, including on the three certified kill
instances of 016). Sign context (part R): MM_sec is NOT generically
negative — it is +0.0004…+0.006 on Bernoulli mixtures, smoothed slices,
and random dense potentials, ≈ 0 at products, and negative on the
witness/ladder genre and (out-of-regime) on δ₀⊕Bern(0.54) at λ = 1. The
violation is structural to the two-block geometry, not a Jensen artifact
of the definition.

### C. The adversary vs MM_abs (`gap1c_partC.json`)

θ re-optimized against MM_abs (260 steps, marginal-penalized, n = 48,
λ ∈ {2, 3.5}), then transferred to n = 64…160 with per-n dilution re-tune;
free-support climbs at n = 8 against MM_abs (8 × 600 steps, witness-seeded
and random-seeded, penalty 300); results:

- PENDING-C

### D. Candidate II inside the window (`gap1c_partD.json`)

Ladder scans at λ = λ_win(n), λ_win/2 (and λ_win/4 for n < 224) for θ\*
and raw weights, n ∈ {48, 96, 160, 224, 320}; θ re-optimized AT window λ
(n = 48 and 96, 300 steps), transferred with per-n window λ up to n = 320;
joint (θ, λ ≤ λ_win) climb at n = 96; small-λ coefficient fits:

- PENDING-D

### E. Raw-weight ladder extension, 016 lead 3 (`gap1c_partE.json`)

Honest (dilution-swept minimum) raw-witness-weight ladder at λ = 2:

- PENDING-E

### F. Exact rational certification (`gap1c_partF.json`)

The decisive point is certified float-free by a new interval path: exact
census tables (013's engine), z_target enclosed by 80-step dyadic bisection
of the exact Plackett quadratic (per distinct margin pair — the ladder's
unit exchangeability collapses thousands of rows to a handful of distinct
(x, y) pairs), h₂ enclosed via directed log₂ enclosures at dyadic interval
endpoints (straddle-guarded at z = 1/2), all sums exact rationals:

    witness, t = 181/16 (λ = log₂(181/16), exactly rational tilt):
      MM_sec numerator/denominator enclosure:
      MM_sec ∈ [−0.0222546582056…, −0.0222546582056…]  CERTIFIED NEGATIVE
      plain aggregate ∈ [+1.844669005…, +1.844669005…] (= 013's certificate
      digit-for-digit, an internal anchor computed from the same tables)
      max marginal < 38271/100000 exactly; degeneracy dichotomy exact.

- PENDING-F

Robustness of the kill (part R, `gap1c_partR.json`): 20/20 independent 3%
weight perturbations keep MM_sec < 0 in-regime (range −0.0228…−0.0218);
the 2-significant-digit tidy witness gives −0.0223 at max marginal 0.313.
Convention robustness is a one-line proof this time: a margin-degenerate
history has z̃ equal to its boundary value and z_target pinned to the same
boundary (z_ρ(0, y) = 0, z_ρ(x, 1) = x, etc. for every ρ > 0), so its
secant contribution is exactly 0 — including degenerate histories in the
average changes neither the numerator nor the sign, only the denominator,
exactly as in 016's robustness note.

## Outcome

**REFUTED (certified): the margin-modulated odds-ratio control in both
signed readings — the assembly-exact secant form and the
derivative-at-target form — fails at n = 7, in the record-relevant regime,
on 007's own 10-atom witness.** MM_sec ∈ [−0.02225465820566…,
−0.02225465820566…] at the exactly-rational tilt t = 181/16, all marginals
< 0.38271, dichotomy exact — a theorem about the stated rational inputs by
013's standard. The h-sensitivity weighting was conjectured (007 lead 2) to
neutralize the witness's light-slice hiding mechanism; instead the witness
kills the weighted statement while its plain aggregate is certified
POSITIVE (+1.8447) on the same tables. No large n, no replication, no
re-optimization was needed.

**EVIDENCE (scoped): the two surviving statements resist the full ladder
adversary.**

- **MM_abs** (the unsigned |σ|-weighted control): positive on everything
  tried — PENDING-SCOPE-ABS
- **λ-window-restricted control**: positive at every window point tried —
  PENDING-SCOPE-WIN

**Not claimed:** no claim that MM_abs or the window-restricted control
holds (finite battery, one adversary genre family plus free-support climbs
at small n; 016 showed exactly how such survivals can be small-n or
wrong-genre artifacts); no claim about the dependent-couplings route
dying — the interface, licensing lemma, separations and the 0.4315 ceiling
are untouched; the MM_der violation is float EVIDENCE (−0.90 with the
certified secant identity tying the two readings; the derivative form
itself was not separately certified); no claim that the ladder is the
worst family for either survivor. Per repo rules nothing here is a result
until an independent skeptic pass; the exact MM certificate is the thing
to attack first (it is a finite list of rational-table claims plus one
dyadic bisection per distinct margin pair), then the census definition's
match to 007 §1 (anchored by the +1.844669005 reproduction), then the
search-completeness claims, which are the weakest.

## Why it failed / what survived

**Why the margin modulation backfires: the weight is signed, and the
witness's surplus lives on the wrong side of z = 1/2.** The sensitivity
σ = ∂h₂(z_ρ)/∂λ is positive only while z_ρ < 1/2. In the record regime the
diagonal surplus histories (the ones whose large positive log₂OR − λ the
plain aggregate spends) sit at zero-margins x, y ≈ 0.6–0.7 with z̃ close
to min(x, y) — i.e. z̃ > 1/2, where h₂ is *decreasing*: their realized
gain h₂(z̃) is LESS than the target gain, so the modulated functional
converts the plain aggregate's biggest positive terms into negatives. The
deficit history (the light slice) contributes almost nothing in either
direction — its margins are nearly degenerate, which is what 007 lead 2
correctly anticipated — but the modulation kills the *surplus*, not the
deficit. Numerically, on the witness at λ = 3.5 the plain census has one
−0.122 deficit against several multi-bit surpluses; the secant census has
those same surpluses contributing ≈ −0.02 net. The golden-ratio
coincidence (z crosses 1/2 exactly at marginal (3−√5)/2 for equal margins)
means this sign flip is intrinsic to the record regime: any
sensitivity-weighted control inherits it, which is why the failure needs
no adversarial tuning at all.

**What the |σ|-weighted form dodges, and what it costs.** Taking |σ|
restores positivity on everything tried — but it decouples the functional
from the chain rule: |σ|·(log₂OR − λ) is no longer the derivative of any
gain, so a proof of MM_abs ≥ 0 would not by itself feed the 008/012
assembly. PENDING-WHY-ABS

**Why the window survives the ladder (so far).** PENDING-WHY-WIN

Reusable:

- **The three-way disambiguation of 007 lead 2** — the margin-modulated
  candidate is not one statement but three, and only the unsigned one is
  alive. Any future proof effort must target MM_abs specifically, knowing
  it lacks direct assembly meaning.
- **The exact MM certification path** (`exact_mm` in the engine): certified
  interval evaluation of Plackett-root functionals over exact census
  tables, with margin-pair caching that makes ladder-scale instances
  affordable. First tool in the library that certifies a *nonlinear*
  functional of the census (everything before was linear in log₂OR).
- **The product-measure structural zero** (MM_sec ≡ 0 at products, exact):
  a free calibration anchor for any future engine touching these
  quantities.
- The part-R friendly-genre battery: MM_sec sign map showing the violation
  is specific to the two-block geometry.

## Leads generated

PENDING-LEADS

## References

- This repo: `attempts/007-averaged-or-control.md` (lead 2 — the candidate
  under test; §1 definitions; the witness), `attempts/013-*` (exact
  standard, part C), `attempts/016-*` / `017-*` (ladder, θ\*, C2 window
  correction, leads 1/3/4), `attempts/009-*` / `011-*` (λ-window law),
  `attempts/008-*` / `012-*` (h(z_ρ) gain, assembly budgets).
- Tools/data: `explore/uc_gap1_candidates.py`;
  `data/gap1c_part[A-F,R].json`, `data/gap1c_*.log`.
- No external papers consulted this cycle (record threshold 0.38271 per
  prior records, Liu).
