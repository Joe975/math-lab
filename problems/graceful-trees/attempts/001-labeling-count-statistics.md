# 001 — Exact graceful-labeling counts for all trees, n ≤ 14

**Queue item.** SAT encoding for graceful labeling plus labeling-count
statistics.

**Status of counting run:** complete for all trees on n = 4..14 vertices
(5,444 trees total); budget-capped before n = 15.

## Approach

Count, exactly, the number of graceful labelings of every non-isomorphic
tree on n vertices, then study the distribution: which structures minimize
the count (where a counterexample would have to live), which maximize it,
and how the minimum grows with n.

Design decisions:

- **Backtracking counter, not #SAT.** The queue item offered a choice; for
  *counting* (as opposed to existence) every solution must be enumerated,
  and a specialized bitmask backtracker beats a generic #SAT solver by
  orders of magnitude at this size. A SAT encoding remains the right tool
  for *existence* on much larger single instances (see Leads).
- **Counting convention.** For each isomorphism class we fix one
  representative and count graceful injections f: V → {0..n−1}, divided
  by 2 for the complementation symmetry f ↦ (n−1)−f (always a distinct
  graceful labeling, so the raw count is even). This raw count is **not**
  quotiented by tree automorphisms; we additionally store |Aut(T)| (AHU
  canonical forms) and analyze the **essential count** ess = raw/|Aut| —
  labelings up to automorphism — which is the right notion of "how many
  genuinely different graceful labelings does T have". (A star has raw
  count (n−1)! but ess = 1: its labeling is essentially unique.)
- **Search.** Vertices labeled in BFS order from a max-degree root, so each
  new vertex's edge difference is fixed immediately; used-label and
  used-difference bitmasks; complementation halved at the root; pruning by
  realizability of the largest unused difference via an "open vertex"
  label mask (a future edge must join an unplaced vertex to an open placed
  vertex or another unplaced vertex — three mask operations). C core
  (`tools/graceful_core.c`, compiled on demand), 4-way parallel over
  interleaved tree chunks; the pruning + bit-iteration gave ~6× over the
  naive backtracker with **identical** counts on all of n = 10..12.

## What was done

- `tools/graceful.py` — tree generation (networkx `nonisomorphic_trees`),
  the counters (C core, pure-Python fallback, brute-force |V|!-permutation
  checker), |Aut| via AHU, structural features (max degree, diameter,
  leaves, caterpillar / lobster / spider flags), a budgeted production
  runner with per-n checkpointing, and the analysis. Subcommands:
  `validate`, `run`, `analyze`. Fully re-runnable.
- **Validation (all passed):**
  - tree counts n = 1..14 match OEIS A000055 (…, 1301, 3159);
  - brute force == Python backtracking == C core for **all** trees
    n = 4..8;
  - stars K_{1,n−1} give exactly (n−1)! raw (n = 5, 9, 11, 13) — i.e.
    ess = 1; (this check initially *failed* because I had wrongly expected
    raw count 1; the validation suite caught the conceptual error and
    forced the raw/essential distinction that shaped the whole analysis);
  - AHU |Aut| == VF2 self-isomorphism count for all trees n = 4..9;
  - path raw counts ×2 = 4, 8, 24, 32, 40, 120, 296, 648, 1328 (n = 4..12)
    match OEIS A006967 (graceful permutations) — independent external
    corroboration;
  - optimized core reproduces the unoptimized core's counts exactly on all
    892 trees n = 10..12.
- **Production run:** all trees n = 4..14, per-n JSON checkpoints in
  `attempts/graceful-trees/data/graceful_counts_nNN.json` (fields: graph6,
  count, maxdeg, diam, leaves, caterpillar, lobster, spider, |Aut|).
  n = 14 (3,159 trees) took 593 s on 4 cores; the runner's projection
  honestly stopped before n = 15 (~2.2 h projected). The hard bound is
  intrinsic: enumeration cost ≥ sum of counts, and the star alone
  contributes (n−1)!/2 ≈ 4.4·10^10 at n = 15.
- **Every tree tested has count ≥ 1** — no counterexample (as expected;
  known verified far beyond n = 14).

## Outcome

Distribution of **essential** counts (raw/|Aut|):

| n | #trees | min | 2nd-min | min non-caterpillar | min non-lobster | median | max |
|---|--------|-----|---------|--------------------|-----------------|--------|-----|
| 8 | 23 | 1 | 3 | 8 | — | 15 | 52 |
| 10 | 106 | 1 | 5 | 12 | 42 | 82.5 | 367 |
| 12 | 551 | 1 | 4 | 19 | 291 | 783 | 5,249 |
| 13 | 1,301 | 1 | 3 | 20 | 172 | 2,551 | 21,107 |
| 14 | 3,159 | 1 | 4 | 25 | 553 | 8,638 | 84,746 |

- **Minimizers.** The global essential minimum is pinned at **1 for every
  n** — always the star (essentially unique labeling), and the runner-ups
  (ess 3–7, flat in n) are high-symmetry diameter-3/4 "double brooms"
  (two hubs carrying most leaves). So the minimum does **not** grow — but
  every minimizer lies deep inside classes *proven* graceful (stars,
  caterpillars; diameter ≤ 5 trees are also settled). The flat floor is a
  statement about rigid known-graceful trees, not about danger to the
  conjecture.
- **Restricted minima grow.** Excluding caterpillars, min ess =
  5, 8, 11, 12, 19, 19, 20, 25 (n = 7..14), fit ≈ 1.18^n (slow, and the
  minimizers are still diameter-4/5 broom-like *lobsters*). Excluding all
  lobsters — i.e. trees outside every well-understood near-path family —
  min ess = 42, 209, 291, 172, 553 (n = 10..14), fit ≈ **1.64^n**:
  robustly exponential. Median ess ≈ 3.0^n, max ess ≈ 3.75^n.
- **Raw-count minimizer is the path** for every n = 7..14 (raw = 2, 4, 12,
  16, 20, 60, 148, 324, 664, 1600, 4956 for n = 4..14; sole exception
  n = 6, a 3-leg spider with raw 6). Under the raw convention nothing is
  remotely close to 0 and the minimum grows ≈ 2.7–3× per vertex.
- **Folklore test — half confirmed, half refuted.** Low end confirmed:
  small-diameter, high-max-degree, high-symmetry trees (stars, near-star
  spiders, brooms) have the fewest essential labelings. High end
  *refuted* for caterpillars: from n = 10 the maximizers are
  **non-caterpillar lobsters** with max degree 4–5, diameter ≈ n/2, and
  few leaves (near-paths with several length-2 branches). At n = 14 the
  top-20 contains 0 caterpillars and 16 lobsters; the path (ess 2,478) is
  a factor 34 below the maximum (84,746). Branching slightly, at depth 2,
  beats both the path and all caterpillars.
- **Novelty is limited:** after the run, a literature check found D. Anick,
  *Counting graceful labelings of trees: a theoretical and empirical
  study*, Discrete Appl. Math. 198 (2016) 65–81, which already built the
  full labeling database for trees on ≤ 16 edges and likewise found the
  minimal-count trees fall in two known caterpillar families (consistent
  with our brooms/double brooms). Our additions relative to Anick: the
  |Aut|-normalized essential counts, the class-restricted minima
  (non-caterpillar ≈ 1.18^n, non-lobster ≈ 1.64^n growth), and re-runnable
  in-repo tooling.

## Why it failed / what survived

- **Failed / capped:** n = 15..18 unreachable for exact counting here —
  cost is bounded below by the counts themselves, which contain (n−1)!
  (star). Any push further needs symmetry-adapted counting (count
  essential classes directly, dividing work by |Aut|) or per-class closed
  forms; brute enumeration is star-bound regardless of pruning quality.
- **Failed (conceptual, caught by validation):** the initial expectation
  that a star has "1 labeling" conflated labelings-up-to-automorphism with
  labelings of a fixed representative. The suite's star check failed,
  which is exactly what it was for.
- **Survived:** the counter and its validation chain (brute force, OEIS
  A006967 paths, A000055 tree counts, old-vs-optimized core agreement);
  the full exact dataset n ≤ 14 with structure features and |Aut|; the
  headline observations (flat essential floor inside known-graceful
  classes; exponential growth outside lobsters; lobster maximizers).
- **Interpretation, honestly stated:** the non-growing global minimum is
  *not* evidence for a nearby counterexample — the minimizers are exactly
  the rigid families where gracefulness is a theorem. For trees outside
  all known-graceful/near-path classes the count grows ≈ 1.6^n, which is
  quantitative *evidence* (not proof) that no counterexample lives at
  small n.

## Leads generated

1. **Hardest-instance mining target:** the non-lobster minimizers are the
   right seeds for SAT-based hardest-instance search at n ≫ 14, where
   existence-SAT scales far beyond counting. The n = 13 dip to ess = 172
   is the graph6 tree `LpH?GCAO??_@?A` (|Aut| = 48, maxdeg 3, diam 6,
   degree sequence 3,3,3,3,2,2,2,1^6): a symmetric ternary-branching
   "double spider" — all vertices of degree ≤ 3 but branching at *both*
   depth 1 and depth 3, i.e. exactly the shape excluded from every proven
   class. This motif (bounded degree, nested branching, high symmetry) is
   the concrete profile to mine.
2. **Prove the flat floor:** diameter-4 double brooms appear to have O(1)
   essential graceful labelings (3–7, non-monotone in n). Diameter-4 trees
   are proven graceful; an exact formula for double-broom labeling counts
   looks tractable and would turn an empirical floor into a theorem.
3. **Symmetry-adapted counter:** counting essential classes directly
   (orbit-level backtracking with symmetry breaking) removes the |Aut|
   factor from the cost and would likely reach n = 16–17, matching and
   extending Anick's database within this repo's tooling.
4. **Maximizer structure:** "maximally graceful trees" are sparse-branch
   lobsters (maxdeg 4–5, diam ≈ n/2). Characterizing them (cf. Anick's
   criterion for exceptionally many labelings; MathWorld's "maximally
   graceful tree") is a clean small open question the lab could attack.
5. **Lobster verification push** (from the problem file's attack surface):
   since lobsters are simultaneously the *maximizers* and (their broom-like
   end) the *nontrivial minimizers*, the famous open lobster subclass is
   where count statistics are most informative; a lobster-only generator
   would let both tails be tracked to larger n.

## Reproduction

```
python3 tools/graceful.py validate           # full validation suite
python3 tools/graceful.py run --nmin 4 --nmax 15 --budget-min 38
python3 tools/graceful.py analyze
```

Data: `attempts/graceful-trees/data/graceful_counts_n{04..14}.json`.

References: Anick, DAM 198 (2016) 65–81
(sciencedirect.com/science/article/pii/S0166218X15002814); OEIS A000055,
A006967; MathWorld "Maximally Graceful Tree".
