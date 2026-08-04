# 014 — Aggregated OR control probed to n = 32 (trend runs to n = 64): survives, certified

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed
- **Type:** computational probe (queue item 1, the mandated first step:
  extend the certified probe of the restated gap to n ≳ 20 BEFORE any
  proof effort). Deliberately NOT a proof attempt.
- **Tools:** `explore/uc_or_agg_probe.py` (written and run here; standard
  library only; deterministic, fixed seeds; parts X/O/P/T/S/H/E, runnable
  independently). Checkpoints `data/or_agg_probe_part[XOPTSHE]*.json`, log
  `data/or_agg_probe_run.log`. Commands reproducing every number below:

      python problems/union-closed/explore/uc_or_agg_probe.py 2>&1 | \
          tee -a problems/union-closed/data/or_agg_probe_run.log
      # or per part, in the order run here:
      python .../uc_or_agg_probe.py X   # engine cross-checks (~15 s)
      python .../uc_or_agg_probe.py O   # orbit family census (~3 min)
      python .../uc_or_agg_probe.py P   # perturbative d-sweep (~80 s)
      python .../uc_or_agg_probe.py T   # trend runs to n = 64 (~2 min)
      python .../uc_or_agg_probe.py S   # sparse structured census (~5 s)
      python .../uc_or_agg_probe.py H   # aggregate hill-climbs (~3.5 min)
      python .../uc_or_agg_probe.py E   # exact certification (~2 min)

- **Sources:** 005/006/007/008/009/011/012/013 (this repo); no external
  fetches. Engines consulted for conventions and the orbit idea:
  `explore/uc_or_avg.py` (007), `explore/uc_or_avg_skeptic.py` (013),
  `explore/uc_pert_skeptic.py` (012's orbit engine). All code here is
  fresh; nothing is imported from them.

## The exact object probed

Notation as in 005/006/007: tilt coupling `pi_lambda(A,B) prop
u(A) u(B) 2^{lambda|A&B|}` with both marginals mu; `t = 2^lambda`;
histories `a = A_{<i}`, `b = B_{<i}`; record threshold 0.38271.

With 007 §1's well-posed per-coordinate average (reused verbatim):
`P_i` = reachable (i−1)-bit prefixes of supp(mu); `N_i` = prefixes in
`P_i` with BOTH continuations reachable; for `(a,b) in N_i x N_i` the
conditional 2×2 table is entirely positive (degeneracy dichotomy, proved
in 007) and

    M_i(mu, lambda) = sum_{(a,b) in N_i^2} m_ab log2 OR_i(a,b)
                      / sum_{(a,b) in N_i^2} m_ab ,
    m_ab = pi(A_{<i} = a, B_{<i} = b).

The **i-aggregated odds-ratio control** (007 §5 item 1, the restated Gap 1
after 007/013 killed the per-i form `M_i >= lambda` with the 10-atom
witness; certified positive only at n <= 7 before this attempt, 013 part C:
+1.844669) is the claim whose sign this attempt probes:

    A(mu, lambda) = sum_{i < n, N_i != empty} w_i (M_i - lambda)
                    / sum_{i < n, N_i != empty} w_i   >=  0
    for all mu with every elementwise marginal < 0.38271, all lambda >= 0,

with `w_i` = the pi-mass of `N_i x N_i` ("defined mass"; a common
normalization cancels in the ratio). `i = n` is excluded because
`M_n = lambda` identically (005 Prop 1); coordinates with `N_i` empty
contribute nothing.

**Convention robustness of the sign (new observation, two lines).** Under
013's conjecture-friendliest alternative bookkeeping — score every
margin-degenerate history at exactly lambda instead of conditioning it
out — every degenerate history adds `w * (lambda - lambda) = 0` to the
numerator and positive mass to the denominator. The numerator of A is
unchanged and the denominator only grows, so **the SIGN of A is identical
under either convention**. A sign probe of the conditioning-out aggregate
is therefore convention-robust for free (unlike the per-i magnitude
statements 013 had to re-check separately).

## Approach

Two engines, then exact certification of the tightest points:

- **Sparse engine** (any n, m atoms, O(m² n)): direct atom-pair census
  from a potential, using 005's trick (any positive u is its own
  marginal's Sinkhorn potential). This is what makes the witness genre,
  crash/mirror, gadget stacks and adversarial searches reachable at
  n = 32 with no approximation — the cost of the obvious alternative
  (full-support Sinkhorn + census) is 4^n and dies at n ≈ 13.
- **Orbit engine** (n to 64+): for measures invariant under permutations
  of coords {3..n} — every structured adversary genre on file (Bernoulli
  cells, Sawin geometric mixtures, delta_0+Bern, Chase–Lovett slice and
  smoothed slice, eps-mixed crash, 012's crash-mixture cells). Histories
  collapse to orbit classes `(oa, ob, ja, jb)` (special-coordinate
  patterns, tail-prefix weights) with the prefix-overlap j summed into a
  mass factor `JW[ja][jb] = sum_j Ncl(ja,jb,j) t^j`; cells per class are
  exact orbit sums through tail pair-kernels
  `PK[s][l] = C(f,s) sum_w C(s,w) C(f-s,l-w) t^w`. The OR of a class does
  not depend on j (it cancels), which is what keeps the census O(n^4)
  per instance. Same S_{n-2} symmetry as 012's budget engine; fresh code,
  cross-checked against the direct engine at n = 6, 7 (part X2).
- **Exact rational certifier** (the 2026-07-31 standard: float agreement
  is NOT certification): both engines re-run end-to-end in
  `fractions.Fraction` at exact rational tilts (181/16, 4, 271/256), all
  potentials dyadic rationals, every log2 replaced by a certified
  enclosure — atanh series with explicit tail bound and **directed dyadic
  fixed-point compaction** (every intermediate rounded toward the safe
  side onto 120-bit dyadics, so long certified sums stay gcd-free and
  cheap; this is what makes exact certification affordable at n = 32).

Why these adversaries: the queue names them. The 10-atom witness genre
(near-proportional light slices, anti-aligned cross-ratios) is the only
known killer of the per-i form; Sawin mixtures were 007's closest census
approach; crash/eps-crash are 005/006's genres; the crash-mixture (p,d)
cells are exactly the genre whose budget growth reversed past n ≈ 22 in
012 (the standing warning this probe exists to answer); the lambda grid
includes the 011 window law `lam_max ~ 4.847/(n-3)` (half/at/twice) plus
the tilt-recipe cells log2(1.06), log2(1.20) and O(1) values up to 5.

## What was done

### X — engine cross-checks (all PASS)

- Sparse engine reproduces 007's witness at lambda 3.5: `M_5 - lambda =
  -0.122033`, aggregate +1.844754, max marginal 0.3178.
- Own exact path reproduces 013 part C identically: aggregate at
  t = 181/16 in [+1.844669005341 ± 1.4e-15] — same number 013 certified.
- Orbit vs direct engine on full supports (slice, Sawin, eps-crash,
  d0+Bern at n = 6, 7): max |M_i or aggregate diff| = 8.0e-13.
- Exact orbit rows vs float orbit rows: diff 7.9e-15.

### O — orbit family census: 840 instances, 0 violations

14 families × 10 lambdas × n in {10, 14, 20, 24, 28, 32} (660 of 840
in-regime; Sinkhorn residual <= 1e-11 on all; `|M_n - lambda| <= 1.4e-13`
sanity on all). Key structure in the result:

- **Product measures are the exact equality set.** Every Bernoulli cell
  gives aggregate = 0 to float noise (|A| <= 9.3e-15 across all n,
  lambda): product mu has a product potential, hence proportional
  response slices at every history, hence `OR = 2^lambda` pointwise
  (005 Prop 2's equality case), hence `M_i = lambda` for every i and
  A = 0 **identically**. The aggregated control, if true, is tight on
  the entire product family — so any proof must be second-order at that
  boundary, and any violation only needs a downward curvature direction.
- Tightest non-product in-regime point per n, always at the smallest
  grid lambda (~0.083–0.084) and always the 012 genre `crashmix_383_05`
  = 0.95·Bern(0.383)^n + 0.05·crash:

      n                10        14        20        24        28        32
      min aggregate  +4.88e-4  +6.49e-4  +8.50e-4  +9.76e-4  +1.10e-3  +1.21e-3

  The minimum margin **grows monotonically with n** — no 012-style
  reversal in the aggregate on this genre up to n = 32.
- Family ordering at the tight lambda (n = 32): crashmix_383_05
  +1.2e-3 < crashmix_30_05 +1.4e-3 < partE_mix +3.5e-3 < sawin_geo_a
  +4.3e-3 < smoothed slice +6.2e-3; slices, d0bern, eps-crash all >= 1e-2.

### P — perturbative sweep around the product equality boundary: 432 instances, 0 negative

Since products are the equality set, the kill direction is perturbative:
`mu_d = (1-d) Bern(p)^n + d nu`, nu in {crash, empty, full, slice at the
record layer}, p in {0.3823, 0.30}, d in [0.005, 0.2], lambda in
{log2 1.06, 0.5, 4.847/(n-3)}, n in {10, 20, 32}. Result: the aggregate
leaves the boundary **quadratically upward in every direction tested**
(fitted order d^1.85–d^2.06 at small d; positive at every point).
Contrast with 012 recorded explicitly: on the same crash-mixture genre
whose per-coordinate averaged DOWNWARD budget grows with n (killing
008's B1), the aggregate's margin grows with n too — the aggregate
absorbs the growing downward budget because the compensating upward
deviations grow faster in the mass-weighted sum.

### T — trend extension to n = 64 on the softest directions: 0 negative

The one decaying trend in part P — the slice-direction curvature at
p = 0.30 (aggregate at d = 0.05, lambda 0.084: 1.35e-5 at n = 10 →
6.2e-7 at n = 32) — chased to n = 64 with the orbit engine:

      n              32         40         48         56         64
      slice p=.30  +6.15e-7   +2.22e-7   +1.41e-7   +2.33e-7   +4.31e-7

The margin **dips to +1.4e-7 near n = 48 and then recovers** — a genuine
non-monotonicity (the 012 lesson about small-n reading is live), but no
zero crossing to n = 64. All other tracked directions (matched-layer
slice, both p; crash direction at fixed lambda and at the moving window
lambda; crashmix_383_05) grow or plateau positive.

### S — sparse structured census at n in {10..32}: 426 instances, 0 violations

- Pure crash / mirror (4 atoms, Theorem-B-safe: aggregate >= 0 is a
  theorem there; run for the margin trend).
- **Witness genre re-laid-out at every n** (marker coord 1, fresh-prefix
  dilution singletons 2..n-3, response n-2, futures n-1, n; by 007's
  dilution invariance the response coordinate keeps `M_i - lambda =
  -0.122` exactly at lambda 3.5, independent of n), with the dilution
  weight swept 0.25..40 and to 1e6 for the limit: the in-regime aggregate
  decays like Theta(1/n) — at lambda 0.5: +0.163 (n=10) → +0.0549 (n=32),
  with n·A → ~1.76 — because a bounded numerator (the gadget's fixed
  deficit plus fixed surpluses) is spread over Theta(n) coordinates'
  defined mass. Positive at every n, lambda, wdil; the wdil → infinity
  limit converges to a positive constant (dilution coordinates are
  single-diagonal-history coordinates, `M_i >= lambda` by 005 Prop 2).
- **Multi-gadget stacks** (007 lead 1 generalized): K = floor((n-2)/5)
  copies of the witness's light-slice gadget on disjoint coordinate
  blocks, each atom carrying its gadget's marker coordinate so the
  blocks do not pollute each other's active prefixes. All positive
  (tightest +1.83, at n = 20).

### H — adversarial hill-climbs on the aggregate: 82 endpoints, 0 violations

Marginal-penalized (target 0.375), dynamic-range-clamped (007's underflow
guard) multiplicative climbs with atom add/drop moves, minimizing the
aggregate itself: free sparse supports (8–18 atoms) at n = 10, 14;
witness-seeded at n = 10, 14, 20; multi-gadget-seeded at n = 14, 20, 24;
small-lambda climbs (lambda in {0.084, 0.15, 0.3, 4.847/(n-3)}) at
n = 10, 14. Lowest in-regime endpoint +0.005602 (n = 10, lambda 0.084 —
consistent with the O(lambda) scaling of everything at small lambda; on
the lambda >= 0.5 climbs the lowest endpoint is +0.20).

### E — exact rational certification: 15 certificates, all positive

Every certificate is a sign theorem about a stated exact rational measure
(dyadic potential, exact rational tilt), computed in Fraction end-to-end
with certified log2 enclosures; enclosure widths <= 5.9e-12. Scope
nuance, as in 013 R2: the orbit certificates are for quantized refit
measures mu-tilde (36-bit dyadic potentials, Sinkhorn residual <= 1e-11
from the named family), which suffices for the universally-quantified
claim; the sparse witness-genre certificates are for the exactly-stated
atom weights themselves.

    point                                   n    tilt t        aggregate (certified)
    witness genre, tightest in-regime wdil  10   181/16        +1.709542464
      (same, per n)                         14   181/16        +1.474040147
                                            20   181/16        +1.162947590
                                            24   181/16        +1.013509091
                                            28   181/16        +0.896763710
                                            32   181/16        +0.803594797
    witness genre, lambda = 2 EXACTLY       20   4             +0.483059569
    crashmix_383_05 (tightest census fam)   10   271/256       +0.000476924
      (same, per n; exact max marginal      14   271/256       +0.000634269
       0.382350 < 0.38271 checked as an     20   271/256       +0.000829968
       exact rational at every n)           24   271/256       +0.000952638
                                            28   271/256       +0.001075588
                                            32   271/256       +0.001201865
    softest direction (slice-dir of         20   271/256       +0.000004637
      Bern(.30), d = 0.05)                  32   271/256       +0.000000599

At t = 4 the n = 20 statement is fully rational with lambda = 2 exactly —
no irrational number appears anywhere in it. The Bernoulli equality
boundary is deliberately NOT in this table: an enclosure there straddles
0 by construction and certifies nothing (recorded as boundary, not as a
near-violation).

## Outcome

**EVIDENCE — the i-aggregated odds-ratio control survives everything this
probe threw at it, now far beyond its previous n <= 7 support.** Scope:
0 in-regime violations across ~1780 float instances (840 orbit-census at
n in {10, 14, 20, 24, 28, 32} × 14 families × 10-point lambda grids in
[0.079, 5]; 432 perturbative instances around the product boundary; 40
trend instances to n = 64; 426 sparse structured instances including the
witness genre and gadget stacks at every n; 82 adversarial climb
endpoints at n <= 24), and 15 exact-rational certificates (all positive,
widths <= 5.9e-12) at n up to 32 including the tightest point found
anywhere (+5.99e-7, certified). The minimum margin over the census GROWS
with n on every structured family; the single decaying direction turns
around at n ≈ 48 without crossing zero.

**Not claimed:** no claim that the aggregated control is true — this is
finite-n, finite-grid, finite-family EVIDENCE, and part T's dip-and-
recover is a reminder that trends here are not asymptotics; no claim
about lambda outside [0.079, 5] or off-grid; no claim for adversaries
outside the tested genres (the climbs are hill-climbs, not exhaustive —
in particular the space of K-gadget interactions is sampled thinly); the
orbit certificates certify quantized refit measures, not the ideal
family members (013 R2 scope); the equality-boundary observation is a
corollary of known facts (005 Props 1–2), not new mathematics; per repo
rules nothing here is a result until an independent skeptic pass.

**Queue-item verdict:** the mandated probe condition is met — the
aggregate is now certified positive at n = 32 (vs n <= 7 before), the
margin trend is growing-in-n on every known adversary genre, and the
012-style reversal was looked for and not found (including on 012's own
genre). Proof effort on the aggregated control is now justified by the
queue's own criterion.

## Why it failed / what survived

Nothing failed — the probed statement survived; what the probe exposed is
WHERE it is hard and where it is soft:

- **The obstruction to any easy proof is the product equality boundary.**
  A = 0 identically on all product measures (proportional slices force
  `OR = 2^lambda` pointwise), so the aggregated control has an
  (n+1)-parameter equality set inside the regime, and any proof must
  show the aggregate is a NONNEGATIVE-CURVATURE functional at every
  product point in every admissible direction, plus something global.
  The probe measured that curvature to be positive (d^1.85–d^2.06
  departure, part P) in crash/empty/full/slice directions — but the
  slice-direction curvature at p = 0.30 decays to +1.4e-7 by n = 48
  before recovering. **The quantity a proof must control is the Hessian
  of A at product measures; the quantity an adversary should attack is
  its smallest eigenvalue over directions, which this probe only sampled
  along four rays.**
- The witness genre's per-i deficit is REAL at every n (M_{n-2} - lambda
  = -0.122 exactly, dilution-invariant) but the aggregate absorbs it with
  Theta(1/n) room to spare, and no stacking of gadgets found here beats
  the bookkeeping: each gadget brings its own positive coordinates.
- 012's warning did not materialize for the aggregate: the growing
  downward budget on crash mixtures is outrun by the growing upward
  deviations in the mass-weighted sum. (Measured, not explained — see
  lead 3.)

Survived / reusable:

- `explore/uc_or_agg_probe.py`: the tail-exchangeable ORBIT CENSUS for
  M_i/aggregate (n = 64 reachable in seconds, arithmetic-generic so the
  same code path runs float and Fraction), the sparse any-n census, and
  the **directed dyadic fixed-point enclosure kit** (`_dround`,
  `_ln_1_to_2_fp`, compacted `agg_enclosure`) that makes exact
  certification ~60x cheaper than 013's Fraction-naive path — reusable
  for any future certified census in this route.
- The equality-boundary framing: aggregated control = "second-order
  nonnegativity at products + global". This is the statement to prove or
  kill, and it connects directly to 007's Theorem C (first order in
  lambda is a perfect square) and to 008/012's perturbative machinery.
- The witness-at-n layout and multi-gadget builder as hard-instance
  generators for any future aggregate claim.
- Trend tables (parts O/P/T checkpoints) for margin-vs-n on every genre.

## Leads generated

1. **Prove or kill the product-boundary curvature.** Compute the Hessian
   quadratic form of `A(mu_d)` at d = 0 for product mu analytically (sum
   007's Theorem-A pairing over i; the first order vanishes on the
   boundary, so the d² coefficient is a computable quadratic form in the
   direction nu). Falsifiable both ways: numerically minimize the form
   over nu (an eigenproblem per (n, p, lambda); the orbit engine gives
   the matrix for tail-exchangeable nu) — a negative eigenvalue with an
   in-regime eigendirection kills the aggregated control; a proof that
   the form is PSD for all p, lambda settles the boundary layer and
   names the remaining global gap.
2. **Map the slice-direction dip.** The +1.4e-7 minimum at (p = 0.30,
   record-layer slice, d = 0.05, lambda ~ 0.084, n ~ 48) is the softest
   point known. Sweep p in [0.25, 0.38] × layer ratio in [p, 0.383] ×
   n to 96 (orbit engine, hours at most) and fit where the minimum moves:
   if some cell drives it below 0 the gap dies at EVIDENCE level; if the
   dip is bounded away from 0 uniformly, record the floor.
3. **Explain the 012 contrast.** On crash mixtures the per-coordinate
   downward budget grows with n (012) yet the aggregate margin grows
   too. Decompose the aggregate numerator per coordinate into down/up
   parts (the engine already has the rows): is the compensation local
   (same i) or cross-i? A local identity would be a provable lemma;
   cross-i compensation would say the aggregate is the WRONG restatement
   and 007 lead 2 (margin-modulated, h-sensitivity-weighted) is the one
   to probe next.
4. **Margin-modulated control (007 lead 2) is still untested** and is the
   assembly-relevant statement; the census infrastructure here evaluates
   it with a one-line weight change. Run the same battery.
5. **(proof side, now unblocked by the queue's own criterion)** Attack
   `sum_i` of 007's Theorem-A pairing directly: the diagonal funding
   `sum_i sum_a m_aa e_a` is 005 part B's strict-excess budget; the probe
   says the inequality has growing room at every tested family, so a
   lossy bound may suffice away from products, with lead 1 handling the
   boundary layer.

## References

- This repo: `attempts/005-odds-ratio-control-refuted.md` (Props 1–2, 5–6,
  crash family, slice reduction); `attempts/006-skeptic-review-of-005.md`
  (S8, lambda > 0 equality correction); `attempts/007-averaged-or-control.md`
  (well-posed M_i, the witness, dilution invariance, aggregated form =
  §5 item 1, leads 1–2); `attempts/008-perturbative-assembly.md` +
  `attempts/012-skeptic-review-of-008.md` (orbit-engine idea, budget-growth
  reversal warning, crash-mixture cells); `attempts/009-mutual-information-tax.md`
  + `attempts/011-skeptic-review-of-009.md` (lambda-window law 4.847/(n-3),
  tilt-recipe cells); `attempts/013-skeptic-review-of-007.md` (exact
  certification standard, t = 181/16, R2 refit-scope convention,
  surrogate-lambda bookkeeping).
- Tools/data: `explore/uc_or_agg_probe.py`;
  `data/or_agg_probe_part[XOPTSHE]*.json`; `data/or_agg_probe_run.log`;
  conventions cross-checked against `explore/uc_or_avg.py`,
  `explore/uc_or_avg_skeptic.py`, `explore/uc_pert_skeptic.py`.
- No external sources; the atanh-series ln bound is standard calculus,
  re-derived in the tool's docstrings (as in 013).

## For the skeptic (load-bearing claims, in attack order)

1. The 15 exact certificates (part E), especially crashmix_383_05 at
   n = 32 (+0.001201865435) and the softest point (+5.988e-7 at n = 32):
   they rest on my ~60-line directed-rounding enclosure kit
   (`_dround` / `_ln_1_to_2_fp` / compacted `agg_enclosure`), which is
   NOT independently reviewed. Re-derive the directed-rounding soundness
   and re-run at least one certificate with an independent exact path
   (013's slower enclosure code is on disk and takes ~1 min at n = 10).
2. The orbit-engine formulas: class count `Ncl`, mass factor `JW`, cell
   `= t^{oa.ob + j + alpha.beta} S(ja+alpha, jb+beta)`, the j-independence
   of the class OR, and the structural N_i reach test. Hand-derivable;
   the direct-engine cross-check (part X2) is only at n = 6, 7.
3. The claim that Bernoulli products give A = 0 identically (equality
   boundary) — two-line corollary of 005 Props 1–2 but it carries the
   whole "census minima at bern rows are noise" reading of part O.
4. The aggregate's sign convention-robustness argument (numerator
   unchanged, denominator grows under surrogate-lambda scoring).
5. The trend claims (growing margin in n; dip-and-recover at n ~ 48) are
   float-only and grid-limited — checkpoints carry every number.
