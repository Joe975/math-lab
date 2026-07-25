# Erdős–Gyárfás Conjecture

**Statement.** Every graph with minimum degree 3 contains a cycle whose
length is a power of 2.

**Status.** Open. Known for planar claw-free graphs and some other classes;
cubic graphs are the natural hard case (min degree exactly 3, no slack).
Computational searches have ruled out small counterexamples (cubic graphs up
to ~modest vertex counts have been checked historically).

**Attack surface.**
- Computational counterexample search: generate cubic (and min-degree-3)
  graphs avoiding cycles of length 4, 8, 16, ... — use SAT/CP encodings or
  targeted constructions (large girth ≥ 9 kills C4 and C8, so candidates
  need many vertices; quantify the trade-off).
- Structural: a counterexample must have cycles only in the "gaps" between
  powers of 2; study cycle spectra of cubic graphs (what length sets are
  realizable?) — this subproblem is interesting and publishable in itself.
- Verify/replicate the known checked bound for cubic graphs before extending it.

**Realism.** Best pure counterexample-hunting target in the set; tooling
built here (cycle-spectrum computation, constrained graph generation)
reuses across problems.
