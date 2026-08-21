# 058 — The n = 7 difficulty was a missing move, not volume: one added move reaches the equality family

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** adversary improvement with a matched control (057 lead 2).
- **Tools:** `explore/uc_hu_collapse_move.py` (new; 040/055's anneal
  plus one move, with a same-code same-seed control; deterministic,
  seeds 8500/8600; checkpoint `data/hu_collapse_move.json`).
  Reproduce: run it.
- **Sources:** none.

## Approach

057 established that the equality family is a stable attractor at
n = 7 but that the unseeded anneal never finds it, and offered two
readings: pure volume (the space doubles each step) or a missing move.
Its move set — weight kick, atom add, atom drop, pairwise transfer —
**never proposes a diagonal directly**: getting from a ten-atom
support to {∅, full} needs eight coordinated drops, each of which the
temperature has to accept in sequence.

The falsifiable test: add exactly one move — **collapse toward a
diagonal**, transferring a slice of one atom's mass onto ∅ or the full
set — and re-run the same campaign. Crucially with a **control**: the
same code, same seeds, same budget, move disabled, so the comparison
isolates the move rather than the reimplementation.

## What was done

Three runs × 2,700 steps per block, essentiality enforced throughout
(055's standard):

    block                                   floor        diagonal endpoints
    n=7 cap 0.49, WITH the collapse move   −1.227e-14         3 of 3
    n=7 cap 0.49, control (move disabled)  +2.300e-05         0 of 3
    n=6 cap 0.49, WITH the collapse move   −1.008e-14         2 of 3

**With the move, every n = 7 run reaches the equality family**
(margins at float noise around the exactly-0 value certified in 056);
**without it, none does**, at identical seeds and budget. The move
also works at n = 6 without breaking what already worked there.

## Outcome

- **REFUTES 057's volume reading, in favour of its alternative: the
  n = 7 difficulty was the move set.** One move closes a gap that
  tripling the budget did not (057: +2.484e-04 → +1.840e-04 at 3×
  steps; here: → −1.2e-14 at the same 3× steps with the move).
- **EVIDENCE: the essential floor at n = 7 is 0, now attained by an
  unseeded search** — which upgrades 057's seeded-control argument
  from "the family is stable there" to "a reasonable adversary finds
  it there".
- **Method result:** an adversary that cannot construct the known
  extremal family in one move will systematically overstate floors.
  The collapse move belongs in the standing anneal, and any future
  extremal family discovered on this route should be checked against
  the move set the same way.
- **Not claimed:** that the augmented anneal is a *strong* adversary
  in general (it is tested at one cap and two sizes); that reaching
  the family means no lower point exists (the family is the conjectured
  extremal, and nothing here tests below it); any certificate — the
  endpoints are float, though 056's identity applies verbatim to the
  diagonals found.

## Why it failed / what survived

Nothing failed. This is the cleanest instance of the window's
recurring lesson, and the first where the fix was constructive rather
than a caveat: five times a number needed a second measurement to
interpret, and here the second measurement (a matched control) not
only interpreted the number but repaired the instrument that produced
it.

The standing state of the route is unchanged and slightly better
supported: no violation at any size or cap tested, the equality
families are attained in essential form at n = 6 and 7, and every
positive "floor" on record above those families is a property of a
search whose weaknesses are now catalogued.

## Leads generated

1. **Re-run 045's best-order caps (0.495–0.499) with the augmented
   anneal.** That campaign used descents only; the anneal with the
   collapse move is demonstrably stronger and those caps are where the
   conjecture is closest to its boundary.
2. **Audit the move set against the other equality families.** 042's
   block tensors are not reachable in one collapse move either
   (they need coordinated mass on several block-diagonal atoms), so
   the same failure mode may still hide there.

## References

- This repo: 057 (the two readings this decides), 056 (the certified
  n = 6 equality endpoint), 055 (the constrained anneal), 042 (the
  equality family), 051/052 (the essentiality standard).
  `data/hu_collapse_move.json`.
- No external sources.
