# Lonely Runner Conjecture

**Statement.** k runners on a unit circular track start together and run at
pairwise distinct constant speeds. The conjecture: each runner is at some
time "lonely" — at distance ≥ 1/k from every other runner.

**Status.** Proved up to k = 7 (k = 7 by Barajas–Serra); a Sept 2025
preprint by Rosenfeld (arXiv:2509.14111) claims a proof of k = 8 —
treat k = 8 as likely settled pending peer review, and k ≥ 9 as the open
frontier. Standard reduction: speeds may be assumed integers, and one
runner stationary. Tao proved it suffices to check speeds up to roughly
k^{O(k^2)}, making each k a (huge) finite problem. NOTE (corrected
2026-07-25): the conjectured "tight ⇔ speeds {1,...,k-1}" rigidity is
FALSE as stated — Goddyn–Wong constructed additional tight instances
(accelerations and sporadics); our attempt 001 recovered them at k = 8.

**Attack surface.**
- Computational: verify k = 8 for structured/small speed sets; look for
  near-violating speed tuples (gap barely above 1/(k+ε)) and study their
  arithmetic structure (they tend to be arithmetic-progression-like).
- The "tight" cases for known k are conjecturally classified (speeds
  {0,1,...,k-1} up to scaling); test tightness rigidity computationally for
  k = 8, 9.
- View-obstruction / zonotope reformulation: covering radius of certain
  polytopes — try LP/SDP relaxations for lower bounds on the gap.

**Realism.** Full k = 8 is likely out of computational reach, but structure
mining of near-tight examples is tractable and could support a partial
result (e.g. k = 8 restricted to bounded speeds).
