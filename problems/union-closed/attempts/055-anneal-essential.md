# 055 — The constrained anneal reaches the equality family at n = 6 and 7, and 054's "floors" were stalls

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** adversarial rerun + correction to 054 (054 lead 2).
- **Tools:** `explore/uc_hu_anneal_essential.py` (new; 040's anneal
  move set with essentiality enforced inside move acceptance;
  deterministic, seeds 8100/8149/8200; checkpoint
  `data/hu_anneal_essential.json`);
  `explore/uc_hu_anneal_essential_skeptic.py` (new; CR via the 050
  skeptic's nats history recursion with the rollout order re-derived
  from the records' rules; exit 0). Reproduce: run the two in that
  order.
- **Sources:** none.

## Approach

054 re-ran 038's descents under 051's essentiality constraint and left
040's anneal — the source of the window's worst dimension
overstatement (an "n = 7" floor that is a two-dimensional instance) —
undone. An anneal is a different search: it accepts uphill moves, so
it can wander into a degenerate face and stay there where a descent
would not. Here the constraint is enforced **inside move acceptance**,
so no accepted state ever loses a coordinate.

## What was done

    block              anneals × steps   floor        040's recorded floor
    n=6 cap 0.49          3 × 900      +6.26e-13          +0.000488
    n=6 cap 0.497         3 × 900      +2.079e-04         +0.000180
    n=7 cap 0.49          3 × 900      +2.484e-04         +0.000488

  Zero violations, zero endpoints lost a coordinate, and 2,933 moves
  were rejected for essentiality across the three blocks — the
  constraint was actively binding, not decorative.

**The n = 6 cap-0.49 endpoint is a diagonal.** It converges to
{∅: 0.6174, full: 0.3826} — every coordinate carrying marginal
0.3826, so genuinely 6-dimensional, and an **equality-family member**
(042's diagonals). Margin +6.26e-13 is therefore *expected*: the
adversary found the equality case, which is exactly what a correct
conjecture predicts and is not a violation.

**Comparison honesty.** Two of 040's three recorded endpoints have min
marginal exactly 0 and are **excluded** by this constraint, so those
two rows are not like-for-like comparisons and are labelled as such in
the checkpoint. The one admissible recorded endpoint (n = 6 cap 0.49,
min marginal 0.487) had floor +0.000488, while the constrained anneal
reaches +6.26e-13 — so the constrained search is **stronger**, not
weakened by the constraint.

## Outcome

- **CORRECTION to 054:** its essentiality-constrained figures
  (+0.091924 at n = 6 cap 0.49, +0.079948 at 0.497, +0.135386 at
  n = 7) are **descent stall points, not floors**. The anneal, under
  the same constraint, reaches +6.3e-13 at n = 6 cap 0.49 — five
  orders of magnitude lower. 054's claim that "the honest floor is
  higher than the degenerate one" holds only against the *degenerate*
  comparison; against a strong constrained adversary the essential
  floor at n = 6 is 0, attained on the diagonal. Per the repo rule 054
  is left as written.
- **EVIDENCE: no violation at n = 6 or 7 under the strongest
  constrained adversary tried** — floors +6.3e-13, +2.1e-04, +2.5e-04,
  all ≥ 0, with the sharpest sitting exactly on the equality family.
- **The equality family is reachable in essential form at n ≥ 6**,
  which 042 predicted but no adversary had demonstrated at that size:
  diagonals are essential (all marginals equal and nonzero) and
  extremal, so a good adversary should find them — and this one does.
- **Not claimed:** any certificate this round; a campaign larger than
  3 anneals × 900 steps per block; anything above cap 0.497.

## Why it failed / what survived

The correction is the content, and it is the fifth instance this
window of the same lesson in a new dress: **a weak adversary's stall
is not a floor.** 052 caught floors that were degenerate; this catches
floors that were merely un-converged. The two failure modes push in
opposite directions — degeneracy made floors look *too small*,
stalling makes them look *too large* — which is exactly why both need
to be checked before any floor is quoted as evidence.

What survives, and is now better supported than before: the conjecture
is not violated at n = 6 or 7 by anything this route can build, and
the binding configurations there are the equality families rather than
anything new.

## Leads generated

1. **Quote floors only with an adversary-strength note.** Every floor
   in this route's records is an upper bound on the true infimum
   produced by a *particular* search; the records should say which
   (descent / anneal / enumeration) and for how long, because 054 and
   055 differ by five orders of magnitude on the same quantity.
2. **Re-run the caps above 0.497 with the constrained anneal.** 045's
   best-order attack went to 0.499 with descents only; the anneal is
   demonstrably stronger and has never been pointed there.
3. **Certify the n = 6 diagonal endpoint** as an exact equality
   instance (it should be exactly c\*(p)·H by 042's identity at
   p = 0.3826) — cheap, and it would make "the equality family is
   attained in essential form at n = 6" exact rather than float.

## References

- This repo: 054 (the record corrected here, and its lead 2), 052/051
  (the degeneracy audit and the constraint), 042 (the diagonal
  equality family this endpoint lands on), 040 (the anneal re-run),
  050 (skeptic evaluator), 044 (own-constant standard).
  `data/hu_anneal_essential.json`.
- No external sources.
