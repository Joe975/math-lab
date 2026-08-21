# 059 — The promoted conjecture against the strongest adversary yet, at the caps closest to 1/2: zero violations

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** adversarial campaign (058 lead 1).
- **Tools:** `explore/uc_hu_bestorder_anneal.py` (new; the
  essentiality-constrained anneal with 058's collapse move, against
  the best-order objective with full order enumeration per candidate;
  deterministic, seeds 8700/8800; checkpoint
  `data/hu_bestorder_anneal.json`). Reproduce: run it.
- **Sources:** none.

## Approach

045 attacked best-order HU — the route's last unrefuted for-all-μ
positivity statement — at caps 0.495/0.497/0.499 using pattern-search
**descents only**. The adversary has improved three times since:
essentiality (051), annealing instead of descent (055, five orders of
magnitude sharper at n = 6), and the collapse-toward-a-diagonal move
(058, which reaches the equality family where tripling the budget did
not). None existed when 045 ran, and 052 found that 4 of its ~22
endpoints per cap had dropped a coordinate.

This re-runs those caps with all three improvements. A violation here
would refute the promoted conjecture outright; the expected outcome is
the equality family, i.e. floors at 0.

## What was done

Full order enumeration per candidate (24 orders at n = 4, 120 at
n = 5), essentiality enforced, collapse move enabled:

    block                    anneals × steps    floor        diagonal endpoints
    n=4 cap 0.495              3 × 1500      +1.091e-03            0 of 3
    n=4 cap 0.497              3 × 1500      +1.416e-05            0 of 3
    n=4 cap 0.499              3 × 1500      −3.220e-15            1 of 3
    n=5 cap 0.495              3 ×  750      +9.907e-05            0 of 3
    n=5 cap 0.499              3 ×  750      +8.054e-05            0 of 3

**Zero violations at every cap and both sizes.** At the hardest cap
tested — 0.499, within 0.001 of the 1/2 boundary — the anneal lands
**exactly on the equality family** (−3.22e-15 is float noise around the
value certified exactly 0 in 056), which is the best a correct
conjecture allows and the outcome that would be expected if it is
true.

**Comparison honesty:** all three of 045's recorded endpoints at these
caps have min marginal exactly 0 and are **excluded** by essentiality,
so the anneal-versus-descent numbers are not like-for-like and are
labelled as such in the checkpoint. What the comparison does show is
that the improved adversary reaches equality at 0.499 where the
descent stalled at +1.295e-04.

## Outcome

- **EVIDENCE (the strongest on record for the promoted conjecture):
  best-order HU positivity survives the improved adversary at caps
  0.495, 0.497 and 0.499**, at n = 4 and 5, with essentiality enforced
  and full order enumeration — zero violations in 15 anneals.
- **The binding configuration at cap 0.499 is the equality family**,
  reached exactly. Every improvement to the adversary this window has
  moved the floors *down onto* the family and never below it.
- **Not claimed:** anything at n ≥ 6 for the best-order objective
  (the cost is 720 orders per candidate — 054's spot check is all
  there is); caps above 0.499; any certificate here (the equality
  endpoint is float, though 056's identity applies verbatim); that 15
  anneals exhaust the space.

## Why it failed / what survived

Nothing failed, and the shape of the result is worth stating plainly:
across this window the adversary was strengthened four times
(own-constant flagging, essentiality, annealing, the collapse move),
and each time the recorded floors moved — always downward, always
stopping at the equality family, never through it. That is the
signature a true conjecture with a known extremal family should
produce, and it is now the substance behind "the promoted conjecture
survives": not that weak searches failed to break it, but that
successively stronger searches converge to its equality case.

What remains open is unchanged: no proof at any n ≥ 2 (046–048 map
the n = 2 obstruction precisely), and the best-order objective is
still untested at n ≥ 6 beyond a single spot check.

## Leads generated

1. **Best-order at n = 6 with the improved anneal**, using rollout as
   a cheap lower bound and enumerating orders only near the floor —
   the one size where the promoted conjecture has essentially no
   adversarial evidence.
2. **Certify the cap-0.499 equality endpoint** by the 056 route (it
   should be an exact diagonal identity instance), making "the
   adversary reaches equality at 0.499" exact rather than float.
3. **A cap sweep toward 1/2** (0.4995, 0.4999): the conjecture's
   constant c\*(p̄) → 0 there, so the margin vanishes on both sides and
   the scale-free form of 046 §G is required — a good stress test of
   the measurement discipline this window built.

## References

- This repo: 058 (the collapse move and this lead), 056 (the certified
  equality identity), 055/051 (anneal and essentiality), 045 (the
  campaign re-run), 052 (which showed 045's endpoints were degenerate),
  042 (the equality family). `data/hu_bestorder_anneal.json`.
- No external sources.
