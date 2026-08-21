# 060 — Correction to 059: the caps constrained only the start, so the boundary was never tested — and it holds when it is

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** self-correction + constrained rerun (same-day review of
  059; the repo forbids editing it).
- **Tools:** `explore/uc_hu_nearboundary.py` (new; the 059 campaign
  with a marginal FLOOR added, plus the drift table from 059's own
  checkpoint; deterministic, seeds 8900+; checkpoint
  `data/hu_nearboundary.json`). Reproduce: run it.
- **Sources:** none.

## Approach

059 reported "at the hardest cap tested — 0.499 — the anneal lands
exactly on the equality family". Certifying that endpoint (059 lead 2)
exposed the problem immediately: **its own max marginal is 0.036**,
not 0.499. The cap constrains only the instances the campaign *starts*
from; the admissibility test is max marginal < 1/2, so the adversary
is free to wander to small marginals — and it does, because the
equality family is reachable there and gives margin exactly 0.

## What was done

**C. The drift, from 059's own checkpoint.** The sharpest endpoint of
four of five blocks had left its cap entirely:

    block            sharpest floor    own max marginal    verdict
    n=4 cap 0.495      +1.091e-03          0.4992          stayed
    n=4 cap 0.497      +1.416e-05          0.4846          drifted
    n=4 cap 0.499      −3.220e-15          0.0360          drifted
    n=5 cap 0.495      +9.907e-05          0.0930          drifted
    n=5 cap 0.499      +8.054e-05          0.0390          drifted

  Endpoints that *stayed* near their cap have floors three orders of
  magnitude larger (+1.2e-02 at marginal 0.4998).

**A/B. The real near-boundary floors.** Re-running with a marginal
floor (max marginal ≥ cap − 0.005, so the adversary must stay near the
boundary), 3 anneals per block, full order enumeration:

    n=4 cap 0.495:  floor +0.000e+00  marginals in [0.4918, 0.5000]
    n=4 cap 0.497:  floor +5.346e-04  marginals in [0.4978, 0.5000]
    n=4 cap 0.499:  floor +8.925e-04  marginals in [0.4954, 0.4988]
    n=5 cap 0.499:  floor +1.597e-02  marginals in [0.4943, 0.5000]

  **Zero violations**, with every endpoint verified to have stayed in
  the band. The cap-0.495 block reaches **exactly 0** — its endpoint
  is the diagonal {∅: 0.5082, full: 0.4918}, an equality-family member
  *inside the band*. That is the correct reading of the other blocks
  too: diagonals exist at every p, so the true near-boundary floor is
  0 at every cap, and the positive numbers at 0.497/0.499 are the
  anneal not finding the family within budget — the 055/057 lesson
  again, now inside the band.

## Outcome

- **CORRECTION to 059:** its headline conflates the campaign's cap
  with the endpoint's own marginal. The equality endpoint it reports
  at "cap 0.499" has marginal 0.036, so **059 did not test the
  boundary**; what it established is the (still valuable) claim that
  the adversary finds the equality family somewhere in regime, with
  zero violations anywhere. Per the repo rule 059 is left as written.
- **EVIDENCE (new, and what 059 intended to show): the promoted
  conjecture holds near the boundary.** With marginals pinned in
  [0.492, 0.500] — the route's first genuinely near-boundary
  best-order numbers — there are zero violations, and at cap 0.495 the
  adversary reaches the equality family *within the band* (floor
  exactly 0, at the diagonal {∅: 0.5082, full: 0.4918}). The positive
  floors at 0.497/0.499 (+5.3e-04, +8.9e-04) and at n = 5 (+1.6e-02)
  are budget artifacts, not bounds: diagonals exist at every p, so the
  true near-boundary floor is 0 everywhere.
- **The same lesson as 044, one level up.** 044 established that
  margins must be judged against the instance's *own* constant rather
  than the campaign's cap. This is the identical error in the search
  region: a campaign labelled by a cap says nothing about the boundary
  unless the adversary is *held* there.
- **Not claimed:** certificates (float-level); anything at n ≥ 6;
  that 12 anneals exhaust the near-boundary region.

## Why it failed / what survived

The sixth cheap escape of the window, and the one that had already
been named: 044's own-constant lesson is exactly this mistake applied
to margins, and it still slipped through in the search region a dozen
records later. The generalisation worth recording: **every constraint
in these campaigns must be checked at the endpoint, not just at the
seed** — essentiality (051/052), dimension (052), and now the cap.

What survives is stronger than 059's version: the conjecture is not
merely unviolated in the in-regime space at large, it is unviolated
*at the boundary*, which is where any counterexample would have to
live.

## Leads generated

1. **Add endpoint-constraint verification to the standing harness.**
   Three constraints have now leaked (dimension, essentiality, cap);
   a single check that every recorded endpoint satisfies the campaign's
   own stated region would have caught all three.
2. **Push the band toward 1/2** (caps 0.4995, 0.4999 with the band):
   c\*(p̄) → 0 there, so the scale-free form of 046 §G is mandatory and
   the measurement discipline gets its real test.

## References

- This repo: 059 (the record corrected here), 044 (the own-constant
  lesson this repeats in the search region), 058/055/051 (the
  adversary), 056 (the certified equality identity), 052 (the earlier
  endpoint-constraint leak). `data/hu_nearboundary.json`.
- No external sources.
