# Ledger

## TL;DR (updated 2026-07-25, cycle 0 wrap-up)

Cycle 0 mostly complete. Erdős–Gyárfás: baseline verified — all 556,471
connected cubic graphs on ≤ 20 vertices contain a power-of-2 cycle; bridged
counterexamples need ≥ 38 vertices. Union-closed: entropy barrier mapped
precisely; four candidate routes past 0.382 recorded (weighted-KL ladder is
the most testable). Erdős–Straus agent lost twice to API overload (529);
retry in progress — if absent from attempts/, it's the top queue item.

## Problem status

| Problem | Status | Budget | Active line |
|---|---|---|---|
| Erdős–Gyárfás | active | high | C4∧C8-free targeted search beyond n=20; cycle-spectrum realizability census |
| Union-closed (Frankl) | active | high | test the weighted-KL ladder (idea B) on Sawin's counterexample; extremal family structure |
| Erdős–Straus | retrying | medium | residue-class identity coverage map (agent lost to 529s twice) |
| Singmaster | queued | medium | binomial collision search design |
| Lonely runner | queued | medium | near-tight speed-tuple mining, k=8 |
| Graceful trees | queued | low | graceful-labeling count statistics; lobster verification |
| Collatz | queued | low (long shot) | failed-approach taxonomy; cycle-bound frontier |

## Attempt queue (next cycles pull from the top)

1. [erdos-straus] If attempts/erdos-straus/ is still empty: residue-class identity coverage map (see problems/erdos-straus.md; the seed agent was killed by server overload — re-run it).
2. [union-closed] Idea B from 001: compute the KL-divergence profile D(U‖Unif F) of Sawin's counterexample family to determine which weight c survives; even c=0.1 would beat the 0.38271 record if admissible.
3. [erdos-gyarfas] Targeted C4-free ∧ C8-free search beyond n=20 (bridgeless only, per the ≥38-vertex bridge bound); girth ≥ 5 pruning makes this far smaller than full enumeration.
4. [singmaster] Design the multiplicity-≥8 collision search (small-k strategy, hash of C(n,k) values); estimate reachable range.
5. [lonely-runner] Implement gap computation for integer speed tuples; scan k=8 small speeds for near-tight cases.
6. [graceful-trees] SAT encoding for graceful labeling; count labelings for all trees n ≤ 18.
7. [erdos-gyarfas] Cycle-spectrum realizability census from the n≤20 data (which length-sets occur?) — standalone interest.

## Verified results

- **[erdos-gyarfas] Cubic baseline n ≤ 20** (2026-07-25, attempt 001):
  every connected cubic graph on 4–20 vertices contains a cycle of length
  4, 8, or 16. 556,471 graphs; counts match OEIS A002851 at every n;
  spectra cross-checked against an independent implementation (networkx
  simple_cycles) on n=10,12 and generation cross-checked vs brute force
  (n ≤ 8). Computational evidence, cubic case only.
- **[erdos-gyarfas] Bridge bound** (2026-07-25, attempt 001): no near-cubic
  bridge-side block on ≤ 17 vertices avoids both C4 and C8 ⇒ a bridged
  cubic counterexample needs ≥ 38 vertices.
- **[union-closed] n ≤ 4 extremal check** (2026-07-25, attempt 001):
  exhaustive — minimum max element frequency over all union-closed families
  on ground sets of size ≤ 4 is exactly 1/2 (2/12/120/4958 families;
  optima triple-checked for closure).

## Insights / cross-problem notes

- Union-closed: the entropy method's ONLY use of closure is H(A∪B) ≤ log|F|
  for iid uniform A,B — an average-case fact, tight at (3−√5)/2 by
  Chase–Lovett's approximate family. Any advance must use worst-case
  closure of atypical/overlapping pairs, dependent couplings, or counting
  structure. See attempts/union-closed/001 §candidate-ideas (A–D).
- Erdős–Gyárfás: "{8}-only" near-misses share a bridged two-block anatomy;
  bridgeless assumption is now safe below 38 vertices. Heawood and
  Möbius–Kantor graphs are the canonical high-girth near-misses.
- Constrained graph generation (nauty geng, installed via apt) + SAT tooling
  is shared infrastructure for Erdős–Gyárfás and graceful trees.
- Ops: parallel subagents can die to 529 Overloaded during API load spikes;
  resume via SendMessage, and don't record a queue item as done until its
  files exist on disk.

## Dead ends

- **[union-closed] k-wise unions**: strictly worsen the entropy constant
  (0.382 → 0.318 → 0.276 for k = 2,3,4); recorded in attempt 001.
- **[union-closed] Gilmer's Conjecture 1** (strengthened entropy inequality):
  refuted by Sawin's construction — do not re-attempt as stated.
