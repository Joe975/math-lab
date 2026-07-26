# Erdős–Gyárfás Conjecture

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** Every graph with minimum degree 3 contains a cycle whose
length is a power of 2.

## Published status

Open. Known for planar claw-free graphs and some other classes. Cubic graphs
are the natural hard case (minimum degree exactly 3, no slack). Computational
searches have historically ruled out small counterexamples.

Relevant published enumerations: OEIS A002851 (connected cubic graphs) and
A014372 (connected cubic graphs of girth ≥ 5) give the ground truth that any
exhaustive search must reproduce.

## Verification contract

- An **exhaustive search claim** at a given n must reproduce the published
  graph counts (A002851 / A014372) exactly, and must state the generator and
  its version. A count that disagrees with the OEIS sequence invalidates the
  run regardless of its conclusion.
- Cycle-spectrum computation must be cross-checked against an independent
  implementation on at least one overlapping range.
- A finite search is `EVIDENCE` for the conjecture and `VERIFIED` only for
  the range actually checked. Say which range.

## Harness (tier 0)

- `harness/erdos-gyarfas/cycle_spectrum.py` — exact cycle spectrum of a
  graph; the reference implementation.
- `harness/erdos-gyarfas/cycle_filter.c` — fast C kernel for exact girth and
  presence of C4/C8/C16. Build: `cc -O2 -o cycle_filter cycle_filter.c`.
- `harness/common/run_slices.sh`, `harness/common/aggregate_slices.py` —
  restartable parallel slice runner and aggregator.

Graph generation uses nauty's `geng` (`nauty-geng`, Ubuntu `nauty` package).
