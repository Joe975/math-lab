# 054 — The high-n floors, re-run with no degenerate escape: every one rises, none crosses

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** adversarial rerun (052 lead 1).
- **Tools:** `explore/uc_hu_essential_highn.py` (new; essentiality-
  constrained rollout descents at n = 6 and 7 and a best-order spot
  check at n = 6; deterministic, seeds 7100–7300; checkpoint
  `data/hu_essential_highn.json`);
  `explore/uc_hu_essential_highn_skeptic.py` (new; CR via the 050
  skeptic's nats history recursion, the rollout order re-derived here
  from the records' rules rather than imported, plus its own
  constrained descents at n = 6; exit 0). Reproduce: run the two in
  that order.
- **Sources:** none.

## Approach

052 showed the route's n = 6 and n = 7 claims rest on endpoints that
had dropped coordinates — most severely an "n = 7" anneal floor that
is a two-dimensional instance. This re-runs those campaigns under
051's essentiality constraint (every coordinate keeps marginal ≥ 0.03)
so the adversary cannot leave the dimension it was launched in, and
asks the falsifiable question 052 posed: **do the floors stay
positive?**

The n = 3 precedent said they should *rise* — there the constraint
moved the floor from +2.6e-09 (degenerate) to +1.8e-05 — because the
degenerate escape was making the problem look harder than it is. The
alternative outcome, a genuinely n-dimensional violation, would have
been the most consequential result this route has produced.

## What was done

Rollout descents (038's P2/P3 reruns) and a best-order spot check,
all own-constant flagged, with every endpoint checked for lost
coordinates:

    block                          descents   floor      previously recorded
    n=6 cap 0.49   (038 P2)            6     +0.091924       +0.030430
    n=6 cap 0.497  (038 P2)            6     +0.079948       +0.053740
    n=7 cap 0.49   (038 P3)            4     +0.135386       +0.099692
    n=6 best-order cap 0.49 (spot)     1     +0.074824            —

  Zero violations, and **zero endpoints lost a coordinate** — the
  constraint held throughout, which is the point of the rerun. Every
  essential floor is **above** the corresponding recorded one, by
  factors of 1.4× to 3×, and the skeptic confirms each comparison
  independently.

  The part-C best-order figure is a **spot check, not a campaign**:
  at n = 6 the objective costs 720 CR evaluations per candidate, so it
  runs one start on a short budget. Stated as such rather than dressed
  up as a floor.

## Outcome

- **EVIDENCE: the high-n floors survive the constraint and rise.**
  Rollout at n = 6 and n = 7 and best-order at n = 6 are all strictly
  positive with no coordinate loss, and the degenerate escapes were
  making the recorded floors look *tighter* than the truth — the same
  direction as the n = 3 precedent.
- **052's falsifiable question is answered in the safe direction:**
  the previously-degenerate high-n claims are not hiding a violation.
  The route's evidence at n ≥ 6, thin as 052 showed, is now genuinely
  n-dimensional at four points.
- **Confirmed independently here (053's correction (ii) does not
  disturb L1):** with t = min(1/2−q, q), the boundary ratio is
  t/(q·c\*(q)), which equals 1/c\*(q) only for q ≤ 1/4 and then
  **grows** above it (4.95 at q = 0.3, 15.3 at 0.45, 70.7 at 0.49), so
  the branch infimum is still approached only as q → 0 and 048's
  conclusion stands as its reviewer said.
- **Not claimed:** a campaign at n = 6 best-order (one short start);
  anything above cap 0.497 at n ≥ 6; that 0.03 is a canonical
  essentiality threshold (it is a choice, stated so the numbers are
  reproducible); any certificate this round.

## Why it failed / what survived

Nothing failed. The useful content is the direction of the shift: in
every case the honest floor is *higher* than the degenerate one, so
the four cheap escapes this window catalogued (a vanishing constant, a
vanishing coverage, and twice a vanishing dimension) were all making
the conjecture look harder to defend than it is, never easier. That is
worth stating because the opposite would have been far more serious —
a degeneracy that flattered the conjecture would have hidden
violations rather than manufacturing near-misses.

What remains genuinely open is unchanged: no proof at any n ≥ 2, and
the n ≥ 6 evidence is now honest but still thin (four constrained
endpoints, one of them a single-start spot check).

## Leads generated

1. **A real best-order campaign at n = 6** — the spot check needs to
   become 8–17 starts. The cost driver is the 720-order enumeration;
   a branch-and-bound over orders, or rollout as a lower bound with
   enumeration only near the floor, would make it affordable.
2. **Re-run 040's anneal under the constraint.** This record covers
   038's descents; 040's anneal (the source of the worst
   overstatement, n = 7 → effective n = 2) still has no
   essentiality-constrained counterpart.
3. **Fold essentiality into the standing harness** so future campaigns
   cannot regress: the check is three lines and belongs next to the
   own-constant flag.

## References

- This repo: 052 (the audit and its lead), 051 (the essentiality
  constraint), 050 (the n = 3 precedent and the skeptic evaluator),
  044 (own-constant standard), 038/040 (the campaigns re-run),
  053 (the reviewer correction re-checked here).
  `data/hu_essential_highn.json`.
- No external sources.
