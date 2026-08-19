# 035 — Canonical-order HU: the rule that rescues cap 0.45 dies at 0.49, certified

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** construction + adversarial attack + exact certification
  (033 lead 4 / queue 1(b1) follow-up). Same author as 031/033; the
  fresh-context reviewer batch covering 030/031/033 was commissioned
  separately and does not cover this record.
- **Tools:** `explore/uc_hu_canon.py` (new; the order rule, the witness
  test, weight-descent attack, fixed-seed sweeps; deterministic — RNG
  only in start shapes, seeds 917000+/921000+; checkpoint
  `data/hu_canon.json`); `explore/uc_hu_canon_skeptic.py` (new;
  independent nats-based conditional entropies, independent greedy
  returning the revelation sequence, CR via the 025/031 skeptic stack;
  exit 0); `explore/uc_hu_canon_certify.py` (new; exact-rational
  certification, each kit alone per 029; checkpoint
  `data/hu_canon_certify.json`). Reproduce: run the three in that order.
- **Sources:** none.

## Approach

033 killed the for-all-orders form of (HU-TAX) at cap 9/20 with a
witness whose bad order reveals its deterministic coordinate LAST — a
surplus-sterile slot. That dissection suggested its own defence, so this
attempt turns the heuristic into a definition and then attacks it.

**Canonical order π_can(μ):** greedily pick, among unrevealed
coordinates, one minimising the conditional entropy H(A_i | A_S), S the
already-chosen set; ties (detected with an explicit tolerance) broken by
lowest index. It reads only μ and the chosen *set* — never realised
values — so both copies use the same order and canonical-HU is a genuine
coupling, still total and closed-form.

Why this rather than optimising the order per instance: a best-order
rule needs an oracle over n! orders and cannot be a recipe; the greedy
rule is O(n²) evaluations and is the one the dissection actually points
at.

## What was done

**A. The 033 killer witness.** Canonical order (2,1,3,0) gives
CR/H = **+0.10142**, against the bad order's −0.00319 and the best
order's +0.10158 — rank 23 of 24. The rule lands essentially on the best
order of the very instance that broke for-all-orders. Certified: the
rescue is CR_HU ∈ [+2.562812450e-1, ·] > 0 at cap 9/20, both kits alone.

**B. Weight-descent against canonical-HU** (17 starts per cap; the order
is re-derived by the rule at *every* candidate, so the adversary cannot
exploit a stale order):

    cap 0.38271: global floor +0.041983  [product extremal +0.04174]  0 violations
    cap 0.45   : global floor +0.001703  [product extremal +0.00728]  0 violations
    cap 0.49   : global floor −0.014704  [product extremal +0.00029]  2 VIOLATIONS

**C. Fixed-seed sweeps** (120 in-regime instances per cap, n ≤ 4, full
order enumeration alongside): at caps 0.45 and 0.49 the canonical order
is negative on **zero** instances, and so is every other order; mean gap
between canonical and the worst order +0.076 / +0.079. The sweep's
random shapes are far easier than the descent's endpoints — the kills in
B are descent products, not typical measures.

**D. The two cap-0.49 kills, certified.** Both reproduce under the
independent builder/evaluator before certification, then certify < 0
under each kit alone:

    n=4 (6 atoms, max marginal 0.48990, H = 2.5811):
        CR_HU ∈ [·, −3.795264395e-2]   canonical order (1,3,0,2)
        16 of 24 orders negative; best order +0.1507
    n=5 (5 atoms, max marginal 0.48816, H = 2.2514):
        CR_HU ∈ [·, −2.888851676e-2]   canonical order (3,0,2,1,4)
        60 of 120 orders negative; canonical IS the worst order here;
        best order +0.0350

**E. Well-posedness — a defect the skeptic caught, and the fix.** The
first implementation broke ties by float comparison; the skeptic's
nats-based implementation then chose a *different* order on the 033
witness (+0.101583 vs +0.101419), because coordinates 0 and 2 have
*exactly* equal conditional entropy at step 2 and rounding decided it.
Ties are real (symmetric measures), so the rule now detects them with a
tolerance and breaks by lowest index; both implementations then agree
exactly. Residual, stated honestly: **lowest-index tie-breaking is not
relabel-equivariant** — it cannot be, since at a tie the rule reads the
labelling. Verified: equivariance holds on all 24 tie-free relabels of a
tie-free instance, and **the n = 4 cap-0.49 kill is entirely tie-free**
(tie trace [1,1,1,1]), so that kill is labelling-independent and
intrinsic. The n = 5 kill has one tie at step 0 and is therefore
labelling-dependent in principle.

## Outcome

- **REFUTED (certified): the canonical-order form of (HU-TAX) at
  marginal cap 49/100.** Two witnesses with CR_HU certified negative
  under each kit alone; the n = 4 one is tie-free, hence not an artefact
  of the tie-break. Since 49/100 < 1/2, **the half-union coupling — even
  with the canonical order — cannot by itself deliver Frankl**: any
  proof through HU must stop below p̄ = 0.49, or add something HU does
  not have.
- **VERIFIED (certified): the canonical rescue at cap 9/20** — the 033
  order-kill witness is positive (+0.2563) under the canonical rule, so
  the rule does exactly what its dissection predicted at that cap.
- **EVIDENCE: canonical-HU survives at the working cap** (0.38271,
  floor +0.041983, tracking the product extremal +0.04174 from above)
  and at 0.45 (floor +0.001703 — where the for-all-orders form died).
- **Not claimed:** that 0.49 is the exact threshold (the descent floors
  are upper bounds on the true inf; the crossing is bracketed only by
  these families, and 033 put the for-all-orders crossing in
  (0.44, 0.45]); anything about best-order HU, which stays positive on
  both kills and is unrefuted everywhere tested — but is an oracle, not
  a recipe; any n > 6 statement.

## Why it failed / what survived

The rule works exactly as designed and still is not enough, which is the
informative outcome. Canonical ordering buys roughly 0.04–0.05 of
marginal cap — it moves the failure threshold from ≈ 0.44 up past 0.45 —
and then hits a wall well short of 1/2. The mechanism is visible in the
n = 5 kill: there the canonical rule picks the *worst* of 120 orders.
Minimising conditional entropy step-by-step is greedy in the wrong
currency near p̄ = 1/2 — it banks predictability early, but what CR
needs is surplus cells later, and those two objectives come apart once
marginals approach 1/2. That is a sharper statement of the 031 averaging
problem than 033 could make: any order rule that fixes this must be
non-greedy, or must optimise the surplus/deficit ledger directly rather
than a per-step entropy proxy.

Consequence for the route: HU's reach is now bounded from above by
measurement (< 0.49) as well as from below by proof (the [1/4,1/2]
sub-genre, 030). Recipe v3 should carry canonical-HU only up to a cap
where it is defended — 0.45 on current evidence — and the interval
[0.45, 1/2) needs a different object.

Reusable: the order rule and its tolerance-tie-break discipline (any
future "canonical" rule needs it — the failure mode is silent and
implementation-dependent); the tie-trace diagnostic; the two 0.49 kill
witnesses as the hardest instances of record; `uc_hu_canon_certify.py`
as the template for certifying a *negative* CR under a derived order.

## Leads generated

1. **Non-greedy order rules**: optimise the surplus/deficit ledger
   (031) directly — e.g. a two-step lookahead, or ordering by expected
   cell surplus rather than conditional entropy. The n = 5 kill is the
   discriminating test: any rule that scores above the worst order there
   is an improvement.
2. **Bracket the canonical threshold** properly: descents at caps
   0.46/0.47/0.48 with n ≥ 5 full-order enumeration, to see whether the
   canonical crossing is a fixed constant or drifts with n.
3. **Best-order HU as a benchmark, not a recipe**: it survives both
   kills; measuring how far above 1/2−ε it survives bounds what ANY
   order rule could achieve, and that ceiling is worth knowing before
   more rule design.
4. Unchanged from 033: prove all-orders CR_HU ≥ 0 at p̄ ≤ 0.38271; the
   witness dissection remains the minimal concrete instance of the
   averaging problem.

## References

- This repo: 033 (the order-kill and its dissection, the heuristic this
  formalises), 031 (HU, HU-notax, the cell ledger, (HU-TAX)), 030
  (HU-mix sub-genre proof), 029 (single-kit certification standard),
  026/022/016 (kits), `data/hu_canon.json`,
  `data/hu_canon_certify.json`.
- No external sources.
