# Graceful Tree Conjecture (Ringel–Kotzig)

**Statement.** Every tree on n vertices admits a graceful labeling: an
injection V → {0, ..., n−1} such that edge labels |f(u) − f(v)| are exactly
{1, ..., n−1}.

**Status.** Open in general. Verified computationally for all trees up to at
least 35 vertices. Known for many classes (caterpillars, paths, stars,
trees of diameter ≤ 5, etc.). Ringel's conjecture itself (which graceful
labeling would imply) was proved asymptotically by Montgomery–Pokrovskii–
Sudakov (2020) — but gracefulness proper remains open.

**Attack surface.**
- Probabilistic/statistical: for random trees, how many graceful labelings
  exist? If the count concentrates well above 0, which tree features
  (many leaves? long paths? high-degree hubs) minimize it? The minimizers
  point at where a counterexample would live — or evidence none exists.
- Computational: push class-specific verification (e.g. lobsters — trees
  within distance 2 of a path — are a famous open subclass; verify lobsters
  to larger n than trees generally).
- SAT encoding of graceful labeling is clean; use hardest-instance mining to
  find structurally "nearly ungraceful" trees.

**Realism.** Counterexample extremely unlikely; value is in labeling-count
statistics (potentially novel data) and hardest-instance structure.
