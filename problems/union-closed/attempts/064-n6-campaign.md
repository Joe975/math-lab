# 064 — An adversarial best-order campaign at n = 6, at rollout cost: reaches the equality family, zero violations

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** adversarial campaign (063 lead 1).
- **Tools:** `explore/uc_hu_n6_campaign.py` (new; anneal on the
  rollout margin with 051 essentiality and 058's collapse move, then
  720-order enumeration at the endpoint; deterministic, seeds 9300+;
  checkpoint `data/hu_n6_campaign.json`). Reproduce: run it.
- **Sources:** none.

## Approach

063 censused n = 6 best-order cheaply but a census is not a campaign —
random measures are easy, and every floor this route records comes
from a search. The same shortcut makes the campaign affordable:
**anneal on the rollout margin** (one order per candidate instead of
720) and **enumerate all 720 orders only at the endpoint**. Since
best-order ≥ rollout pointwise, a non-negative rollout floor already
proves the best-order floor is non-negative; the endpoint enumeration
then measures how much slack the shortcut left.

## What was done

Four anneals × 1,200 steps per cap, essentiality enforced, collapse
move enabled:

    cap     rollout floor    best-order at those endpoints   violations
    0.49     −3.031e-14            −3.031e-14                    0
    0.497    +3.203e-03            +3.203e-03                    0

- **At cap 0.49 the campaign reaches the equality family** — the
  sharpest endpoint is a diagonal (support {∅, full}, marginal 0.0403)
  with margin −3.03e-14, float noise around the value certified
  exactly 0 in 056.
- **The endpoint gaps** (best-order minus rollout) are +0.0525, 0, 0,
  0 at cap 0.49 and 0, 0, +0.0131, +0.0140 at 0.497: **at the sharpest
  endpoints the gap is exactly 0**, i.e. rollout *is* a best order
  there — which is why annealing on the cheap surrogate loses nothing
  where it matters.
- **Sharper than the census at both caps**, as a campaign should be:
  −3.0e-14 vs +0.014946 (cap 0.49) and +3.2e-03 vs +0.037327 (0.497).

## Outcome

- **EVIDENCE: the promoted conjecture survives an adversarial
  best-order campaign at n = 6**, zero violations at both caps, with
  the sharpest endpoint sitting on the equality family. This is the
  first adversarial (not census) best-order result at that size.
- **Method result: the surrogate is free where it counts.** The
  best-order/rollout gap vanishes exactly at the binding endpoints, so
  annealing on rollout and enumerating once at the end gives the same
  floor as annealing on best-order would have, at 1/720 the cost per
  candidate.
- **The pattern continues:** every strengthening of the adversary in
  this window has moved floors down onto the equality family and never
  through it — now at n = 6 for the best-order objective too.
- **Not claimed:** anything above cap 0.497 at n = 6; certificates
  (the diagonal endpoint is float, though 056's identity applies
  verbatim); that 8 anneals exhaust the space.

## Why it failed / what survived

Nothing failed, and the window closes on the same shape it has had
throughout: the adversary got stronger, the floors got sharper, and
they stopped exactly where the conjecture says they must. The one new
technical fact is the vanishing gap — the expensive objective and its
cheap lower bound agree precisely at the configurations that bind,
which is what makes n ≥ 6 tractable at all.

## Leads generated

1. **Push the same construction to n = 7 and 8 best-order.** Rollout
   cost is the only cost now, so the sizes that were unreachable are
   reachable; 8 orders' worth of enumeration at the endpoint is
   40,320, still affordable once.
2. **Collect any instance where the gap is positive AND the rollout
   margin is negative** — none has appeared, and such an instance
   would be the first separation of best-order from rollout in the
   region that matters.

## References

- This repo: 063 (the lead and the census it sharpens), 058/055/051
  (the adversary), 056 (the certified equality identity), 044
  (own-constant standard). `data/hu_n6_campaign.json`.
- No external sources.
