# Erdős–Gyárfás — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`. Full records: `attempts/`.

## Editorial view of the attack surface

- Computational counterexample search: generate cubic (and min-degree-3)
  graphs avoiding cycles of length 4, 8, 16, … — SAT/CP encodings or targeted
  constructions. Large girth ≥ 9 kills C4 and C8, so candidates need many
  vertices; quantify the trade-off.
- Structural: a counterexample must have cycles only in the gaps between
  powers of 2. Study cycle spectra of cubic graphs (which length sets are
  realizable?) — interesting and publishable in itself.

Assessed as the best pure counterexample-hunting target in the problem set;
its tooling (cycle-spectrum computation, constrained graph generation) reuses
across problems.

## Attempts

### 001 — Cycle-spectrum baseline, n ≤ 20 · `VERIFIED` (range only)

Every connected cubic graph on 4–20 vertices contains a cycle of length 4, 8
or 16. 556,471 graphs, zero counterexamples. Counts match OEIS A002851 at
every n; generation validated against an independent pure-Python brute force
(n ≤ 8) and spectra against networkx `simple_cycles` (n = 10, 12).

**Bridge bound.** No near-cubic bridge-side block on ≤ 17 vertices avoids both
C4 and C8, so a bridged cubic counterexample needs ≥ 38 vertices.

### 002 — Girth-≥5 restriction, n = 22 and 24 · `VERIFIED` (range only)

Full enumeration stops being affordable near n = 22, but the C4-free pool is
thin, so generation was restricted to girth ≥ 5. All 90,938 girth-≥5
connected cubic graphs at n = 22, and all 1,620,479 at n = 24, contain a C8.
Zero candidates, zero C8-free graphs.

Cross-checks passed: totals match A014372 exactly; kernel totals match geng's
independent stderr counts slice-by-slice; the girth-6 count 385 at n = 22
matches the known sequence; the unique girth-7 graph at n = 24 (McGee)
appears exactly once as predicted.

## Insights

- The "{8}-only" near-misses share a bridged two-block anatomy; the
  bridgeless assumption is safe below 38 vertices. Heawood and Möbius–Kantor
  are the canonical high-girth near-misses.
- Eight `{8}`-only graphs at n = 24 versus one at n = 22.

## Open lines

- 2-connected C16-free test (lead 2 of attempt 002).
- n = 26 is roughly a 16-hour run.
- Cycle-spectrum realizability census from the n ≤ 24 data — which length
  sets actually occur? Standalone interest.
