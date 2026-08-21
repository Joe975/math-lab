# 056 — The n = 6 essential equality endpoint, certified exactly

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-21
- **Mode:** informed
- **Type:** certification (055 lead 3).
- **Tools:** `explore/uc_hu_diag_certify.py` (new; exact-rational
  certification under each kit alone; checkpoint
  `data/hu_diag_certify.json`). Reproduce: run it.
- **Sources:** none.

## Approach

055's constrained anneal converged at n = 6, cap 0.49 to a diagonal
with margin +6.26e-13 — the first time any adversary on this route
landed on the equality family in **essential** form at n ≥ 6. Float
evidence of an *equality* is weak evidence: +6e-13 is equally
consistent with exact equality and with a tiny genuine gap. This
certifies which.

042's identity says a diagonal {∅: 1−p, full: p} has every HU cell at
x = y = 1−p, hence z = z\*(p) = max(1/2, 1−2p), so

    CR_HU = h(z*(p)) − h(p) = c*(p)·h(p),   H(μ) = h(p),

independently of n — as far as CR is concerned the diagonal behaves
like a single coordinate, while being genuinely n-dimensional (every
marginal equals p). The margin is therefore exactly 0.

## What was done

The 055 endpoint rationalizes to p = 2047143/5349940 ≈ 0.382648, and
all three statements certify under **each kit alone**:

- **M3 (essentiality, exact):** all six coordinate marginals equal p
  **exactly** in rationals — so this is an n = 6 statement, not a
  smaller instance embedded in six coordinates (the failure mode 052
  found across the route).
- **M1 (equality):** CR_HU − c\*(p)·H(μ) is enclosed in
  [−4.930e-32, +4.930e-32] (kit A) and [−5.047e-19, +5.047e-19]
  (kit B) — both contain 0, certifying equality to the kits'
  precision.
- **M2 (positivity):** CR_HU ∈ [+0.040109336, ·] and H > 0.9599,
  certified positive.

## Outcome

- **VERIFIED (certified, each kit alone): the equality family is
  attained in essential form at n = 6.** 042 proved diagonals are
  extremal at every n; this is the first certified instance of an
  adversary reaching one at n ≥ 6, and it confirms 055's +6.26e-13 is
  exact equality rather than a small gap.
- **Consequence for the route's floors:** a constrained adversary at
  n = 6 can reach the equality family, so the essential floor there is
  exactly 0 and cannot be improved. Any future "floor" claim at n = 6
  above 0 is a statement about the search, not about the problem.
- **Not claimed:** anything at n = 7 (055's endpoint there is
  +2.484e-04, not an equality instance); that the diagonal is the only
  essential equality instance at n = 6 (042's block tensors should
  also qualify and are untested here).

## Why it failed / what survived

Nothing failed. The value is in closing a gap the window kept
re-opening: five times a float number near zero was read as meaning
something, and three of those readings were wrong (a vanishing
constant, a vanishing coverage, two vanishing dimensions, a stall).
Here the near-zero number means exactly what it appeared to, and now
that is certified rather than assumed — which is the only way to tell
the two cases apart.

## Leads generated

1. **Certify a block-tensor equality instance at n = 6** (042's
   d2⊗d2⊗d2 at the same cap) to confirm the essential equality set at
   n = 6 is the full family and not just the diagonals.
2. **The n = 7 endpoint (+2.484e-04) is not on the family** — either
   the anneal did not converge there, or the family is harder to reach
   at n = 7. Worth one longer run to decide, since "the adversary can
   always reach equality" is a much stronger statement than "it did at
   n = 6".

## References

- This repo: 055 (the endpoint), 042 (the diagonal identity and the
  equality family), 052/051 (the essentiality standard this endpoint
  satisfies exactly), 029 (single-kit certification), 026/022/016
  (kits). `data/hu_diag_certify.json`.
- No external sources.
