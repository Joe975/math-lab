# 063 — Best-order at n = 6, made affordable: 6,000 instances, zero violations, no enumeration needed

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** census with a cost-reducing observation (059 lead 1).
- **Tools:** `explore/uc_hu_n6_bestorder.py` (new; the rollout
  lower-bound census with an n = 5 control; deterministic, seeds
  9100–9102; checkpoint `data/hu_n6_bestorder.json`). Reproduce: run
  it.
- **Sources:** none.

## Approach

n = 6 was the size where the promoted conjecture had almost no
adversarial evidence: the best-order objective costs 720 CR
evaluations per candidate, so 054 could afford only a single-start
spot check and 059 stopped at n = 5.

The observation that makes it cheap: **rollout is a lower bound on
best-order** — the best order is by definition at least as good as the
one rollout picks — so whenever the *rollout* margin is already ≥ 0,
the best-order bound holds with **no enumeration at all**. The 720
orders need only be enumerated on instances where rollout leaves the
margin negative.

## What was done

Essentiality enforced (051), margins against each instance's own
constant (044):

    block             instances   lower bound sufficed   enumerations   violations   min margin
    n=6 cap 0.49         3000        3000 (100.0%)             0             0        +0.014946
    n=6 cap 0.497        3000        3000 (100.0%)             0             0        +0.037327
    n=5 control          1000        1000 (100.0%)             0             0        +0.043781

**The lower bound sufficed on every one of 7,000 instances** — the
720-order enumeration was never needed. The n = 5 control enumerates
all 120 orders regardless and confirms the shortcut never disagrees
with the enumerated answer (1,000 of 1,000 agreements).

## Outcome

- **EVIDENCE: the promoted conjecture holds across 6,000 essential
  in-regime instances at n = 6**, caps 0.49 and 0.497, with zero
  violations and min own-constant margins +0.0149 and +0.0373. This is
  the size that previously had one spot check.
- **Method result: best-order positivity is cheap to census.** Because
  rollout lower-bounds it and rollout's margin was non-negative on
  every instance tried, the expensive objective never had to be
  evaluated. Any future best-order census at n ≥ 6 should use this
  shortcut and report how often it sufficed — if that fraction ever
  drops below 100%, the instances where it fails are exactly the
  interesting ones.
- **Not claimed:** anything about *descent or anneal floors* at n = 6
  best-order (this is a random census, not an adversarial campaign —
  059's caveat stands); caps above 0.497; any certificate.

## Why it failed / what survived

Nothing failed. The useful content is that the expensive objective was
never the obstacle — a lower bound that happens to be tight enough
retires the cost entirely. That the fraction is exactly 100% across
7,000 instances is itself informative: on random essential measures,
rollout is not merely close to best-order, it is already good enough
to clear the conjectured bound, which is consistent with 038's census
finding that rollout is often but not always the best order.

The route's evidence now covers n = 2 through 7 with no violation
anywhere, and the open mathematics is unchanged: no proof at any
n ≥ 2, with the n = 2 obstruction mapped precisely in 046–048.

## Leads generated

1. **Point the adversarial campaign (not just a census) at n = 6
   best-order** using the same shortcut: anneal on the *rollout*
   margin, and enumerate orders only at the endpoint. That gets an
   adversarial n = 6 best-order floor at roughly rollout cost.
2. **Watch the sufficiency fraction.** If a campaign ever finds
   instances where rollout's margin is negative but best-order's is
   not, those are the first examples separating the two statements at
   n ≥ 6 and are worth collecting.

## References

- This repo: 059 (the lead), 054 (the spot check this replaces), 038
  (rollout-vs-best census at n ≤ 5), 051/044 (essentiality and
  own-constant standards). `data/hu_n6_bestorder.json`.
- No external sources.
