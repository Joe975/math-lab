# Ledger

## TL;DR (updated 2026-07-25, cycle 0)

Lab initialized. Six problems curated. First seed cycle running: tooling +
survey agents on Erdős–Gyárfás (cycle-spectrum search), union-closed
(entropy-barrier analysis), and Erdős–Straus (residue-class identity map).
No results claimed yet.

## Problem status

| Problem | Status | Budget | Active line |
|---|---|---|---|
| Erdős–Gyárfás | seeding | high | cubic-graph cycle-spectrum tooling + counterexample constraints |
| Union-closed (Frankl) | seeding | high | map the (3−√5)/2 entropy barrier; small-family extremal search |
| Erdős–Straus | seeding | medium | residue-class identity coverage map |
| Singmaster | queued | medium | binomial collision search design |
| Lonely runner | queued | medium | near-tight speed-tuple mining, k=8 |
| Graceful trees | queued | low | graceful-labeling count statistics; lobster verification |
| Collatz | queued | low (long shot) | failed-approach taxonomy; cycle-bound frontier |

## Attempt queue (next cycles pull from the top)

1. [erdos-gyarfas] Build + verify cycle-spectrum tool on known small cubic graphs; establish the verified "no counterexample below n vertices" baseline independently.
2. [union-closed] Write up the Gilmer→0.38 chain and the Chase–Lovett obstruction precisely; identify what exact union-closure gives that approximate closure doesn't.
3. [erdos-straus] Program the classical identity families; compute exact uncovered residue classes mod 840 and beyond.
4. [singmaster] Design the multiplicity-≥8 collision search (small-k strategy, hash of C(n,k) values); estimate reachable range.
5. [lonely-runner] Implement gap computation for integer speed tuples; scan k=8 small speeds for near-tight cases.
6. [graceful-trees] SAT encoding for graceful labeling; count labelings for all trees n ≤ 18.

## Verified results

(none yet — nothing enters here without surviving adversarial verification)

## Insights / cross-problem notes

- Constrained graph generation + SAT tooling (Erdős–Gyárfás, graceful trees)
  should be shared infrastructure in `tools/`.

## Dead ends

(none recorded yet)
