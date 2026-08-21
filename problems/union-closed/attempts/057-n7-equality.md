# 057 — At n = 7 the equality family is stable but not found: the gap is search convergence, not structure

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** adversarial probe with a seeded control (056 lead 2).
- **Tools:** `explore/uc_hu_n7_equality.py` (new; the 042 identity
  check, anneals seeded at the family, and long unseeded runs;
  deterministic, seeds 8300+; checkpoint `data/hu_n7_equality.json`).
  Reproduce: run it.
- **Sources:** none.

## Approach

056 certified that the constrained anneal reaches the equality family
exactly at n = 6, but 055's n = 7 run stopped at +2.484e-04. Two
readings with very different consequences: **(i)** the search did not
converge at n = 7, so "a constrained adversary can reach equality"
survives; or **(ii)** the family is genuinely harder to reach at
n = 7, which would mean a size-dependent positive floor the route has
never seen.

The separation is a **seeded control**: put the anneal *on* the family
and see whether it stays. If the family is a local minimum of the
constrained objective at n = 7 the anneal sits there and (i) is right;
if it drifts off, the framing rather than the budget is wrong.

## What was done

**C. The target.** The 042 diagonal identity holds at both sizes: the
margin of a diagonal at p ∈ {0.3826, 0.45, 0.487} is +0.000e+00 at
n = 6 **and** n = 7 (floats; certified exactly 0 at n = 6 in 056). So
the family exists at n = 7 and is a legitimate target.

**A. Seeded control — the family is stable at n = 7.** Seeded at a
diagonal and annealed for 900 steps:

    n=6 p=0.3826:  −3.331e-15   stayed, effective n=6
    n=6 p=0.487:   −3.220e-15   stayed, effective n=6
    n=7 p=0.3826:  −2.109e-15   stayed, effective n=7
    n=7 p=0.487:   −2.665e-15   stayed, effective n=7

  (The tiny negatives are float noise on a quantity certified exactly
  0 in 056, not violations.) **The family is a stable attractor of the
  constrained objective at n = 7 exactly as at n = 6.**

**B. Unseeded, triple budget — still not found.** Three runs at 2,700
steps (3× the 055 budget) at n = 7, cap 0.49, all essential: floors
+3.473e-02, +4.904e-02, **+1.840e-04**. The best is a modest
improvement on 055's +2.484e-04 at a third of the cost per step, and
none reaches the family.

## Outcome

- **Reading (i) is correct: the n = 7 gap is search convergence, not
  structure.** The equality family is present, legitimate and stable
  at n = 7 — an anneal placed on it stays — so the essential floor at
  n = 7 is 0, and 055's and this record's positive n = 7 numbers are
  properties of the search, exactly as 056 concluded for n = 6.
- **EVIDENCE: the unseeded anneal does not find the family at n = 7
  within 2,700 steps**, improving only from +2.5e-04 to +1.8e-04. The
  search cost of reaching a fixed target grows sharply with n — 900
  steps sufficed at n = 6.
- **Not claimed:** that no budget reaches it (only that 3× does not);
  any certificate at n = 7 (the seeded endpoints are float-level
  here, though 056's identity argument applies verbatim); anything
  above cap 0.49 at n = 7.

## Why it failed / what survived

The seeded control is the technique worth keeping: when an adversary
fails to reach a known target, seeding it *at* the target separates
"cannot get there" from "did not get there", and those two have
opposite consequences for the conjecture. Without it, 055's n = 7
number would have stayed ambiguous — and by the pattern of this window
it would probably have been quoted as a floor.

This is the sixth time in the window that a number near zero, or a
number that failed to reach zero, needed a second measurement to
interpret. The through-line is now explicit: **on this route, no
extremal number means anything until something independent says which
of its possible meanings applies.**

## Leads generated

1. **Seed every future high-n campaign at the equality family as a
   control run**, alongside the unseeded search — it costs one anneal
   and tells you whether a positive floor is real.
2. **Why does reaching the family get so much harder from n = 6 to
   n = 7?** The family is one point in a space whose dimension doubles
   each step, so this may be pure volume; if it is, the anneal's
   move set (which never proposes a diagonal directly) is the fixable
   part. Adding a "collapse toward a diagonal" move would test it.

## References

- This repo: 056 (the certified n = 6 equality endpoint and this
  record's lead), 055 (the n = 7 endpoint at issue), 042 (the diagonal
  identity), 051/052 (the essentiality standard used throughout).
  `data/hu_n7_equality.json`.
- No external sources.
