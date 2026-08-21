# 052 — Effective-dimension audit: 11 of 27 recorded floors are quoted at a higher n than the instance uses

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** audit across records (051 lead 1), pure re-analysis of
  committed checkpoints.
- **Tools:** `explore/uc_hu_effdim.py` (new; audits every stored
  descent endpoint in `data/hu_*.json`; deterministic, no search;
  checkpoint `data/hu_effdim.json`). Reproduce: run it.
- **Sources:** none.

## Approach

051 found that 050's n = 3 descents saturated their bound by escaping
to a degenerate face. No earlier campaign on this route had ever been
checked for support collapse, so this audits all of them.

Criterion, and why it is the only one used: **a coordinate whose
marginal is exactly 0 is absent**, so the instance lives at a smaller
n and any "n = k floor" quoting it overstates the dimension. Atom
count is *not* a valid criterion — a diagonal {∅, S} has two atoms and
is genuinely n-dimensional (it is one of 042's equality families) — a
mistake this audit made on its first pass and corrected before
recording. Instances with min marginal below 0.05 but nonzero are
reported separately as near-degenerate, not folded in.

## What was done

**A. All 328 stored endpoints.** **96 (29.3%) have at least one absent
coordinate**; a further 33 are near-degenerate. The affected blocks
span the whole window and the campaigns before it: `hu_attack`
(033), `hu_canon` (035), `hu_order2` (037), `hu_rollcensus` (038),
`hu_roll_anneal` (040), `hu_blocks` (042), `hu_bestorder` (045),
`hu_n3_census` (050).

**B. The headline floors** — the minimum-floor endpoint of each block,
which is what the records actually quote. **11 of 27 are quoted at a
higher n than the instance uses:**

    block                                   claimed n   effective n   floor
    hu_attack:cap_0.38271                       4           2       +0.041983
    hu_canon:B_cap_0.38271                      4           3       +0.041983
    hu_order2:C_surp_cap_0.49                   5           3       −0.012831
    hu_order2:C_surp_cap_0.497                  5           3       −0.020503
    hu_order2:C_roll_cap_0.497                  4           3       +0.000037
    hu_order2:C2_roll_cap_0.497                 4           3       +0.000037
    hu_roll_anneal:A_n6_cap_0.497               6           4       +0.000180
    hu_roll_anneal:A_n7_cap_0.49                7           2       +0.000488
    hu_rollcensus:P2_cap_0.49                   6           5       +0.030430
    hu_rollcensus:P2_cap_0.497                  6           4       +0.053740
    hu_n3_census:D_cap_0.45                     3           2       +0.000000

  The sharpest overstatement is **040's n = 7 anneal floor: effective
  n = 2**. That also explains an oddity noticed at the time — the
  n = 7 anneal reported *exactly* the same floor (+0.000488) as the
  n = 6 run. It was the same two-dimensional instance both times.

**C. What this does and does not change.**

- **Unaffected:** every positivity result and every kill. Each
  endpoint is still a valid in-regime instance, and a certified
  negative CR at effective n = 3 is still a certified negative CR —
  044's and 035's certified kills, 050's certified order-failure and
  051's certified positive floor all stand exactly as recorded.
- **Affected:** dimension labelling. A floor quoted as "the n = k
  floor" may be a smaller instance embedded in k, so it is not
  evidence *about* dimension k. In particular the route's claims to
  have probed n = 6 and n = 7 rest on fewer genuinely high-dimensional
  endpoints than the records imply.

## Outcome

- **CORRECTION (spanning 033/035/037/038/040/050):** 11 of 27 recorded
  headline floors are quoted at a higher n than their instance uses;
  most severely, 040's "n = 7, floor +0.000488" is a two-dimensional
  instance. Per the repo rule those records are left as written and
  the correction lives here.
- **VERIFIED (no change): all positivity and kill claims**, which do
  not depend on the dimension label.
- **Method result:** the effective-dimension diagnostic is three lines
  and should run on every future descent campaign, alongside 044's
  own-constant flag, 046's scale-free margin and 051's essentiality
  constraint.
- **Not claimed:** that the high-n floors are *wrong* — only that they
  are not evidence at the n they are labelled with; genuinely
  n-dimensional floors at n ≥ 6 remain largely unmeasured.

## Why it failed / what survived

Nothing failed computationally; this is a labelling audit, and the
pattern behind it is the one this window kept meeting: a descent finds
the cheapest way to make its objective small, and dropping a
coordinate is cheap. The four cheap escapes now on record are a
vanishing constant (046 §G), a vanishing coverage (046 §G′), a
vanishing dimension in a single campaign (051), and — here — a
vanishing dimension across most of the route's campaigns at once.

What survives, and is the honest state: the conjecture's positivity
has never been violated by any endpoint at any effective dimension,
but the *evidence at n ≥ 6 is thinner than the record implies*, and
recovering it needs essentiality-constrained reruns.

## Leads generated

1. **Re-run the n ≥ 6 campaigns with 051's essentiality constraint**
   (040's anneal and 038's P2/P3 first, since those carry the route's
   only high-n claims). Expect the floors to rise, as they did at
   n = 3 (+2.6e-09 → +1.8e-05); the question is whether they stay
   positive, and that is the falsifiable part.
2. **Ask why the adversary prefers dropping coordinates.** If the
   true infimum over genuinely n-dimensional instances is attained in
   a limit where a coordinate vanishes, then the degenerate faces are
   the extremal structure rather than an artifact — worth deciding,
   because it would say the equality families of 042 are incomplete.

## References

- This repo: 051 (the lead and the essentiality constraint), 050,
  046 (§G, §G′), 044 (own-constant standard), 042 (equality families
  and the two-atom diagonals that make atom count an invalid
  criterion), 040/038/037/035/033 (the campaigns audited).
  `data/hu_effdim.json`.
- No external sources.
