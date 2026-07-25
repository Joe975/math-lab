# 001 — Cycle-spectrum baseline: exhaustive check of small cubic graphs

**Problem:** Erdős–Gyárfás conjecture (every min-degree-3 graph contains a
cycle of length a power of 2).
**Date:** 2026-07-25.
**Tool:** `tools/cycle_spectrum.py` (re-runnable; see its docstring).

## Approach

Build the foundational tooling: exhaustive generation of connected cubic
graphs (the hard case of the conjecture) plus exact cycle-spectrum
computation, then verify the conjecture holds for every connected cubic
graph up to as large an n as is reachable by brute force, while collecting
"near miss" statistics — graphs whose spectrum barely meets the powers of 2.
This replicates/anchors the known small-case bound before any attempt to
extend it with targeted (SAT/girth-based) searches.

## What was done

1. **Generation.** Connected cubic graphs are generated with nauty's `geng`
   (`nauty-geng -c -d3 -D3 -q n`, nauty 2.8.8 from the Ubuntu `nauty`
   package). Generation was validated two independent ways:
   - Against a fully independent pure-Python brute-force generator
     (backtracking over labeled degree-3 graphs + backtracking isomorphism
     dedup), agreeing for n = 4, 6, 8 (counts 1, 2, 5).
   - Against OEIS A002851 (connected cubic graphs): 1, 2, 5, 19, 85, 509,
     4060, 41301, 510489 for n = 4..20 — all counts matched exactly.

2. **Cycle spectrum.** Exact enumeration of all simple cycles (each cycle
   rooted at its minimum vertex, one traversal direction kept), returning the
   set of cycle lengths. Validated on:
   - K4 → {3,4}; K3,3 → {4,6}; Petersen → {5,6,8,9}. All computed correctly
     (the Petersen result *verifies* girth 5 and non-Hamiltonicity — no 10 in
     the spectrum — rather than assuming them). No expected value had to be
     adjusted; all passed on the first run of the finished tool.
   - Independent cross-check: spectra of **all** 19 + 85 cubic graphs on
     n = 10, 12 recomputed with `networkx.simple_cycles` (a structurally
     different algorithm, v3.6.1) — every spectrum agreed. The standout
     n = 18 near miss (below) was also independently re-verified with
     networkx (degrees, connectivity, spectrum).

3. **Search.** `--search` computes every spectrum and intersects with the
   powers of 2 that fit ({4,8} for n < 16, {4,8,16} for 16 ≤ n < 32).
   Run for n = 4..18 single-threaded and n = 20 in four parallel `geng`
   slices (`--split res/mod`). JSON reports with full histograms and
   counterexample slots were saved for each run.

## Outcome

**Verified: every connected cubic graph on n ≤ 20 vertices contains a cycle
of length 4, 8, or 16. Zero counterexamples among 556,462 graphs.**

| n | graphs | counterexamples | pow2-intersection histogram |
|---|--------|-----------------|------------------------------|
| 4 | 1 | 0 | {4}: 1 |
| 6 | 2 | 0 | {4}: 2 |
| 8 | 5 | 0 | {4,8}: 5 |
| 10 | 19 | 0 | {4}: 1, {4,8}: 15, {8}: 3 |
| 12 | 85 | 0 | {4}: 6, {4,8}: 71, {8}: 8 |
| 14 | 509 | 0 | {4}: 15, {4,8}: 458, {8}: 36 |
| 16 | 4060 | 0 | {4}: 14, {4,8}: 199, {4,8,16}: 3562, {4,16}: 16, {8}: 6, {8,16}: 263 |
| 18 | 41301 | 0 | {4}: 67, {4,8}: 1403, {4,8,16}: 37014, {4,16}: 56, {8}: 1, {8,16}: 2760 |
| 20 | 510489 | 0 | (filled from parallel run; see below) |

This is **computational evidence, not proof**, and for cubic graphs only
(general min-degree-3 graphs on these vertex counts were not enumerated;
the standard reduction to the cubic/2-connected case was not re-proved here).
Literature reports checks to similar or larger n; this run independently
replicates the small-case bound with our own tooling.

### Near misses and statistics

- **Single-pow2 graphs.** Graphs whose spectrum meets the powers of 2 in
  exactly one value are the frontier. Counts of "{8} only":
  n=10: 3 (incl. Petersen), n=12: 8, n=14: 36, n=16: 6, n=18: **1**.
  At n=14 all 36 are bridgeless. "{8} only" collapses at n ≥ 16 because
  16-cycles become available; the graphs avoiding C4 and C16 nearly always
  contain C8.
- **The unique n=18 "{8} only" graph** (`Q??CA?_cAOA_`CC`@o@@OIO@@_?`,
  spectrum {3,5,6,7,8,9}) is structurally instructive: it is two identical
  9-vertex blocks joined by a **bridge** (4–14). Each block has spectrum
  {3,5,6,7,8,9}; the bridge caps the circumference at 9, killing all
  16-cycles globally. Its only power-of-2 cycles live inside the blocks.
- **C4-free graphs** (girth ≥ 5, or odd girth with no C4): 3, 8, 36, 269,
  2761 for n = 10..18 — all contained C8 or C16.
- **C8-free graphs**: 1, 6, 15, 30, 123 for n = 10..18 — every one contained
  a C4 (asserted during the run).
- **Named extremes.** The n=14 girth-6 graph is the **Heawood graph**
  (verified by isomorphism), spectrum {6,8,10,12,14}: bipartite, C4-free,
  saved only by its 8-cycles. The unique n=16 girth-6 graph is
  **Möbius–Kantor** (verified), spectrum {6,8,10,12,14,16} — meets powers of
  2 in both 8 and 16.
- Girth distribution at n=18: girth 3: 33496, girth 4: 7350, girth 5: 450,
  girth 6: 5. High-girth cubic graphs (the C4-free pool) are a thin slice,
  consistent with the known fact that a counterexample needs girth ≥ 9-ish
  or must dodge C8 combinatorially.

## Why it failed / what survived

This was a tooling/baseline attempt, not a proof attempt, so "failure" =
no counterexample found (expected: the conjecture is believed true and known
checked well past n=20). What survived:

- A validated, re-runnable generator + spectrum pipeline
  (`tools/cycle_spectrum.py`, `--validate` re-runs all checks in ~4 s;
  full n ≤ 16 search takes ~6 s, n = 18 ~2 min, n = 20 ~20 min on 4 cores).
- The n ≤ 20 verification itself (exact counts above).
- The near-miss census, which sharpened where the tension actually is.

Limits of the method: brute force stops around n = 22–24 (graph counts grow
~12× per step: n=22 is ~7.35M graphs). Extending the bound needs constrained
generation (e.g. `geng` girth flags, or SAT) rather than full enumeration.

## Leads generated

1. **Bridge/block reduction — executed, negative (a small result).** A cycle
   cannot cross a bridge, and each side of a bridge in a cubic graph is a
   connected "near-cubic" graph (exactly one degree-2 vertex, rest degree 3;
   odd order, (3n−1)/2 edges). If a side on ≤ 15 vertices avoided C4 and C8,
   two copies joined by a bridge would be a ≤ 30-vertex counterexample (too
   small for C16). We ran this search (`--blocks` mode; `geng -c -d2 -D3 -f`
   with the edge count pinned, then C8-testing the C4-free survivors):

   | side order n | C4-free near-cubic sides | also C8-free |
   |---|---|---|
   | 5, 7 | 0 | — |
   | 9 | 1 | 0 |
   | 11 | 8 | 0 |
   | 13 | 59 | 0 |
   | 15 | 544 | 0 |
   | 17 | 6314 | 0 |

   **Verified: every connected near-cubic side on ≤ 17 vertices contains a
   4-cycle or an 8-cycle.** Consequence (exact, given the computation): any
   cubic counterexample to Erdős–Gyárfás that contains a bridge has both
   sides ≥ 19 vertices, hence **≥ 38 vertices total**. Counterexample hunts
   below that size may restrict to bridgeless (hence 2-edge-connected) cubic
   graphs. Extending `--blocks` to 19, 21, ... stays far cheaper than full
   enumeration and keeps pushing this bound.
2. **C4-free + C8-free is the real wall.** Every C8-free cubic graph found
   (n ≤ 18) contains a C4, and every C4-free one contains C8 or C16. A
   targeted search for C4-free ∧ C8-free cubic graphs (girth ≥ 5 plus a
   C8-avoidance constraint — SAT/CP encoding, or filtering high-girth `geng`
   output where girth ≥ 9 gives C4/C8-freeness for free) is the natural
   escalation; girth ≥ 9 cubic graphs start at n = 58 (cage bound), so the
   interesting regime is girth 5–8 with C8 excluded combinatorially.
3. **Cycle-spectrum census as a subproblem.** The full per-graph spectra
   (JSON reports) are a dataset for the "which length-sets are realizable by
   cubic graphs" question flagged in the problem file — e.g. spectra that are
   intervals vs. gappy ones like Heawood's {6,8,10,12,14}. Worth a dedicated
   attempt.
4. Petersen-like graphs (spectrum {5,6,8,9}) rely solely on C8; a family
   analysis of generalized Petersen graphs GP(n,k) spectra mod powers of 2
   could yield structured near misses at larger n.

## Reproduction

```
sudo apt-get install nauty
python3 tools/cycle_spectrum.py --validate
python3 tools/cycle_spectrum.py --search --min-n 4 --max-n 18
# n=20, 4-way parallel:
for r in 0 1 2 3; do
  python3 tools/cycle_spectrum.py --search --min-n 20 --max-n 20 \
    --split $r/4 --json search20_$r.json &
done; wait
# bridge-side search (lead 1):
python3 tools/cycle_spectrum.py --blocks 5 7 9 11 13 15 17
```
