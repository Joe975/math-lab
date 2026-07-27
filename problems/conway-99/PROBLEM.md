# Conway's 99-Graph Problem

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** Does there exist an undirected graph on 99 vertices in which
every two adjacent vertices have exactly one common neighbour, and every two
non-adjacent vertices have exactly two common neighbours?

Equivalently: a strongly regular graph with parameters (99, 14, 1, 2) — every
edge lies in a unique triangle, and every non-edge is a diagonal of a unique
4-cycle. Conway offered $1000 for a resolution either way (posed at the 2014
DIMACS conference on identifying integer sequences).

## Published status

Open, and open for a specific reason: **the standard obstructions do not
apply.** The parameter set is feasible. The eigenvalues are 3 and −4 with
multiplicities 56 and 42, all integral; the Krein and absolute bounds are
satisfied. There is no cheap counting reason for the graph not to exist, and no
construction is known.

- The local structure is forced. With k = 14 and λ = 1, the neighbourhood of
  every vertex is exactly 7 disjoint edges (7K₂).
- (99,14,1,2) is one of five parameter sets with λ = 1, μ = 2. Two are
  realised: the 9-vertex Paley graph (9,4,1,2) and the Berlekamp–van
  Lint–Seidel graph (243,22,1,2). Three are open: (99,14,1,2),
  (6273,112,1,2) and (494019,994,1,2).
- Symmetry is heavily restricted. Such a graph cannot be vertex-transitive.
  Makhnev–Minakova (2004) and Behbahani–Lam (2011) showed the automorphism
  group G has order dividing 2·3³·7·11, and that |G| divides 42 if it is even.
  Cesarz–Woldar (arXiv:2308.02978, 2023) refined this: if 7 divides |G| then
  G ≅ Z₇, and consequently if |G| is even then |G| divides 6, so G is one of
  Z₂, Z₆ or S₃.
- Keramatipour (arXiv:2604.23037, 2026) encoded the problem for SAT solvers
  and reports that they do not settle it in reasonable time, with an analysis
  of why the encoding is hard. That is the state of the art on the direct
  computational attack: read it before spending a cycle re-deriving it.

**Sources.** Assembled from encyclopaedia summaries and abstracts rather than
the primary papers; treat as machine transcribed `[T]`. The automorphism
restrictions and the five-parameter-set list are the load-bearing claims.

## Verification contract

The search space is finite and astronomically large, so essentially every
honest claim here is about a **sub-case**. The contract is about making the
case split auditable.

- **Say exactly which case was eliminated, and show the split is exhaustive.**
  "No such graph contains configuration X" is a result only if the reader can
  check that X-free was the complement of what you enumerated. An incomplete
  case split is the failure mode of this whole genre.
- **Exhaustive searches need isomorph rejection with a stated canonical form**,
  and the canonical form must be validated on parameter sets whose verdict is
  already published. A pipeline that cannot recover known strongly-regular
  existence and nonexistence results is not evidence about 99.
- **SAT-based nonexistence needs an independently checkable proof** (DRAT/LRAT
  validated by a checker that did not produce it). A solver printing UNSAT is
  not a result, especially given the published finding that this encoding
  strains solvers.
- **Feasibility arithmetic is exact.** Eigenvalues, multiplicities and Krein
  conditions are integer/rational computations; a floating-point derivation of
  an integrality condition is not acceptable.
- **Status discipline.** Absent an actual graph or a complete nonexistence
  proof, work here is `MAP` or `EVIDENCE` with the case split stated. Nothing
  about the full parameter set is `VERIFIED` by a partial search.

## Harness (tier 0)

None yet. A contributor adding one should put the strongly-regular property
checker, the exact feasibility-condition calculator and the canonical-form
routine here, since those verify the objects themselves rather than any route.
