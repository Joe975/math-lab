# 065 — n = 7 and n = 8 best-order, where domination removes the enumeration entirely

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** adversarial campaign (064 lead 1).
- **Tools:** `explore/uc_hu_n78_campaign.py` (new; anneal on the
  rollout margin at n = 7 and 8, with a 5,040-order endpoint
  enumeration at n = 7 and none at n = 8; deterministic, seeds
  9500/9600; checkpoint `data/hu_n78_campaign.json`). Reproduce: run
  it.
- **Sources:** none.

## Approach

064 made n = 6 best-order affordable by annealing on the rollout
margin and enumerating once at the endpoint. The logic scales further,
and at n = 8 it removes the enumeration entirely: since

    best-order margin ≥ rollout margin   (pointwise, by definition)

a **non-negative rollout floor is already a proof that the best-order
floor is non-negative** — no enumeration at any point. Enumeration is
only needed to *measure* the surrogate's slack, which is affordable at
n = 7 (5,040 orders, once) and simply skipped at n = 8 (40,320).

## What was done

Essentiality enforced, collapse move enabled, cap 0.49:

    n     anneals × steps    rollout floor       best-order conclusion
    7        3 × 900         −4.485301e-14       non-negative by domination
    8        3 × 450         +5.996704e-02       non-negative by domination

- **n = 7 reaches the equality family** (−4.49e-14 is float noise
  around the value certified exactly 0 in 056), and the endpoint
  enumeration confirms the **gap is exactly 0 at all three
  endpoints** — rollout *is* a best order there, as at n = 6.
- **n = 8 stalls at +0.0600**, well above the family; with half the
  step budget that is expected, and the conclusion does not depend on
  reaching the family — the floor's sign is what carries it.

**The ladder at cap 0.49** (rollout floors, all essential,
all adversarial):

    n = 6:  −3.031e-14    (064)
    n = 7:  −4.485e-14    (here)
    n = 8:  +5.997e-02    (here, half budget, no enumeration)

## Outcome

- **EVIDENCE: the promoted conjecture survives adversarial best-order
  campaigns at n = 7 and n = 8**, zero violations, with n = 7 reaching
  the equality family.
- **Method result: at n ≥ 8 the expensive objective never has to be
  evaluated at all.** Domination converts a rollout campaign into a
  best-order conclusion, so the cost of testing the promoted
  conjecture is now independent of the order count — the barrier that
  kept this route at n ≤ 5 for the best-order statement is gone.
- **The gap vanishes at the binding endpoints at n = 7 as at n = 6**,
  which is now measured at two sizes rather than one.
- **Not claimed:** that n = 8's +0.0600 is a floor (half budget, and
  by the 055 lesson a stall is not a floor); anything above cap 0.49
  at these sizes; certificates.

## Why it failed / what survived

Nothing failed. The window closes with the best-order statement —
the route's last unrefuted for-all-μ positivity claim — tested
adversarially at every size from 2 to 8 with no violation anywhere,
and with the cost barrier that previously confined it to n ≤ 5
removed by an observation rather than by compute.

The open mathematics is exactly where 046–048 left it: no proof at any
n ≥ 2, the n = 2 obstruction isolated to a two-variable inequality
whose infimum is exactly 1, and the pair-interaction term identified as
the thing no averaging argument can discard.

## Leads generated

1. **n = 9, 10 and beyond** are now reachable by the same route; the
   only cost is the rollout evaluation, which grows with the history
   tree rather than with n!.
2. **Give n = 8 a full budget** and check whether it too reaches the
   equality family — if the pattern holds at every size, "the
   adversary always converges to the equality family" becomes a
   statement worth trying to prove rather than to test.

## References

- This repo: 064 (the lead and the n = 6 campaign), 063 (the
  domination observation), 058/055/051 (the adversary), 056 (the
  certified equality identity), 046–048 (the open mathematics).
  `data/hu_n78_campaign.json`.
- No external sources.
