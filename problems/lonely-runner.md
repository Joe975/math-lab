# Lonely Runner Conjecture

**Statement.** k runners on a unit circular track start together and run at
pairwise distinct constant speeds. The conjecture: each runner is at some
time "lonely" — at distance ≥ 1/k from every other runner.

**Status.** Open for k ≥ 8 (proved up to k = 7, the k = 7 case by Barajas
and Serra). Standard reduction: speeds may be assumed integers, and one
runner stationary. Tao proved it suffices to check speeds up to roughly
k^{O(k^2)}, making each k a (huge) finite problem.

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
