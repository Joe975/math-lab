# 061 — One check for all three leaked constraints: 34% of recorded endpoints leave the region their campaign claimed

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** cross-record audit + standing harness check (060 lead 1).
- **Tools:** `explore/uc_hu_region_audit.py` (new; audits every stored
  endpoint in `data/hu_*.json` against the region its block name
  declares; deterministic, no search; checkpoint
  `data/hu_region_audit.json`). Reproduce: run it.
- **Sources:** none.

## Approach

Three constraints have leaked on this route, each caught by hand after
the fact: **dimension** (052 — 11 of 27 headline floors quoted at a
higher n than the instance uses), **essentiality** (051 — descents
saturating a bound by escaping to a lower-dimensional face), and **the
cap** (060 — an "equality at cap 0.499" endpoint whose own max
marginal is 0.0360). 060's lead asked for the single check that would
have caught all three. This is it.

Per stored endpoint, four rules:

    R1  no coordinate absent (marginal exactly 0)
    R2  every coordinate marginal >= 0.03    (051's essentiality bar)
    R3  max marginal < 1/2                   (in regime at all)
    R4  if the block name declares a cap, the endpoint's own max
        marginal is within 0.02 of it        (060's bar, deliberately
                                              loose so it flags only
                                              real drift)

Blocks predating a rule are still reported: the output is a map of
which recorded numbers describe the region they claim, not a verdict
on records written before the rule existed.

## What was done

**393 endpoints audited across every committed checkpoint; 135
flagged (34.4%):**

    R1  absent coordinate        96
    R2  below essentiality       31
    R3  out of regime             0
    R4  drifted off its cap       8

**R3 never fires** — the one constraint that never leaked is the one
the engines enforce inside their objective functions rather than at
the seed. That is the whole lesson in one line: **a constraint checked
only when a run starts will leak; a constraint checked inside the
objective does not.** 055's essentiality-in-acceptance and 060's
marginal floor were both built that way after the fact, and neither
leaks in this audit.

The flagged blocks are concentrated in the pre-051 campaigns
(`hu_attack`, `hu_canon`, `hu_order2`, `hu_rollcensus`,
`hu_roll_anneal`), exactly as expected; the post-051 checkpoints
(`hu_essential_highn`, `hu_anneal_essential`, `hu_collapse_move`,
`hu_nearboundary`) are clean.

## Outcome

- **VERIFIED (audit): 34.4% of the route's recorded endpoints leave
  the region their campaign declared**, broken down by rule above. No
  positivity or kill claim changes — every endpoint is still a valid
  in-regime instance (R3 is clean) — but a third of the route's
  recorded *floors* describe a different region than their label.
- **The diagnostic principle, now evidenced:** constraints belong
  **inside the objective**, not at the seed. R3 is enforced that way
  and has zero violations across 393 endpoints; R1/R2/R4 were seed-
  level and leak at 24%, 8% and 2%.
- **The check is committed and cheap** (a few seconds over all
  checkpoints), so future campaigns can run it as a gate.
- **Not claimed:** that the flagged endpoints are wrong — they are
  valid instances, just not evidence about the region their block name
  implies; nor that R4's 0.02 band is canonical (it is loose by
  design).

## Why it failed / what survived

This is the window's seventh measurement-discipline finding and its
last, and it generalises the other six: every one of them was a
quantity that meant something other than its label — a margin against
the wrong constant (044), a vanishing constant (046 §G), a vanishing
coverage (046 §G′), a vanishing dimension (051, 052), a stall read as
a floor (055), a cap that bounded only the start (060). The unifying
repair is not more care; it is **putting the constraint where the
search cannot avoid it, and auditing endpoints rather than intentions.**

## Leads generated

1. **Gate future campaigns on this check** — run it after every
   campaign and refuse to record a floor whose endpoint is flagged.
2. **Re-derive the flagged floors that still matter.** Most flagged
   endpoints are from superseded campaigns, but `hu_order2`'s
   `D_bestorder_*` blocks feed 037/045's headline claims and have
   R1×5 and R1×4 flags; those two are worth a constrained rerun.

## References

- This repo: 060 (the lead and the cap leak), 052 (the dimension
  leak), 051 (essentiality), 055/058 (constraints enforced inside the
  objective, which do not leak), 044 (the original own-constant
  lesson). `data/hu_region_audit.json`.
- No external sources.
