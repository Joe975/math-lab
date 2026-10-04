# 062 — The two live headline floors, recomputed from clean endpoints: one sound, one from a degenerate endpoint

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** targeted correction (061 lead 2).
- **Tools:** `explore/uc_hu_region_audit.py` (061, for the flags) plus
  the recomputation recorded here; no new engine — this reads
  `data/hu_order2.json` and re-minimises over the clean subset.
  Reproduce: the ten-line script is quoted in full below.
- **Sources:** none.

## Approach

061 flagged 34.4% of the route's endpoints as leaving their declared
region, and noted that most sit in superseded campaigns — with two
exceptions that still feed live claims: `hu_order2`'s
`D_bestorder_cap_0.49` and `D_bestorder_cap_0.497`, which supply
037's and 045's best-order floors. This recomputes both headline
numbers using **only** the endpoints that pass 061's rules (every
coordinate marginal ≥ 0.03, max marginal < 1/2).

    for key in ("D_bestorder_cap_0.49", "D_bestorder_cap_0.497"):
        rows = json.load(open("data/hu_order2.json"))[key]["rows"]
        clean = [r for r in rows
                 if 0.03 <= min(marginals(r)) and max(marginals(r)) < 0.5]
        print(key, min(r["floor"] for r in clean))

## What was done

    block                    endpoints  flagged  recorded floor   clean-only floor
    D_bestorder_cap_0.49        14         6     +6.942017e-04    +6.942017e-04
    D_bestorder_cap_0.497       14         5     +2.190730e-04    +2.632841e-04

- **cap 0.49: sound.** The recorded floor already comes from a clean
  endpoint (`product4`, n = 4, all marginals 0.4888). Nothing changes.
- **cap 0.497: the recorded floor comes from a flagged endpoint**
  (`seed0`, min marginal exactly 0 — a coordinate absent, so the
  instance is effectively smaller than n = 4). The honest floor over
  clean endpoints is **+2.632841e-04**, from `product4` at marginals
  0.4969.

## Outcome

- **CORRECTION (small, to 037/045 via `hu_order2`):** the
  `D_bestorder_cap_0.497` floor of **+2.190730e-04** is produced by a
  degenerate endpoint; the clean floor is **+2.632841e-04**. Both are
  positive and both exceed c\*(0.497) = 2.6e-05, so **no qualitative
  claim moves** — the recorded number was simply 20% sharper than the
  region it claimed could justify. Per the repo rule those records are
  left as written.
- **VERIFIED: the cap-0.49 headline floor is clean** and needs no
  restatement.
- **061 lead 2 is discharged:** of the two flagged blocks feeding live
  claims, one was sound and one is corrected here.
- **Not claimed:** anything about the other flagged blocks (they
  belong to superseded campaigns); any certificate.

## Why it failed / what survived

The interesting part is how small this is. After seven
measurement-discipline findings in one window, the live headline
numbers turn out to be off by 20% in one place and correct in the
other — the leaks were real and worth fixing, but they were not
hiding anything structural. That is itself the useful calibration: the
route's *conclusions* were robust to a third of its endpoints being
mislabelled, because the conclusions never rested on a single sharpest
number.

## Leads generated

1. **None new.** 061's gate plus this recomputation closes the
   endpoint-hygiene thread opened by 051. The open mathematics is
   unchanged: no proof at any n ≥ 2, the n = 2 obstruction mapped in
   046–048, and the best-order objective untested at n ≥ 6 beyond a
   spot check.

## References

- This repo: 061 (the audit and its lead), 052/051 (the leaks), 045
  and 037 (the records whose floors these are), 044 (own-constant
  standard). Data: `data/hu_order2.json`, `data/hu_region_audit.json`.
- No external sources.
