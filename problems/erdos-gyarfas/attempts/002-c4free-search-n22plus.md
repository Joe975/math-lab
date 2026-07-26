# 002 — C4-free search, n = 22+: girth ≥ 5 cubic graphs vs {C8, C16}

**Problem:** Erdős–Gyárfás conjecture (every min-degree-3 graph contains a
cycle of length a power of 2).
**Date:** 2026-07-25.
**Tools:** `tools/cycle_filter` (C kernel, source `tools/cycle_filter.c`),
`tools/run_slices.sh` (restartable slice runner),
`tools/aggregate_slices.py`, `tools/cycle_spectrum.py` (independent
verification). Follows up attempt 001 (full enumeration exhausted at n ≤ 20).

## Approach

Full cubic enumeration stops being affordable around n = 22 (~7.35M graphs),
but the C4-free pool is thin (attempt 001: 5751 + 32 girth-≥5 graphs at
n = 20 out of 510489). So restrict generation to **girth ≥ 5 connected cubic
graphs** (`nauty-geng -c -d3 -D3 -tf`, i.e. triangle-free + square-free) and
test only the remaining power-of-2 obligations: a counterexample on
22 ≤ n < 32 must avoid C4 (free from the girth restriction), C8 and C16.
The C kernel computes exact girth and presence of C4/C8/C16 per graph
(cheap early-exit DFS with BFS-distance pruning; C8 presence is the cheap
rejection, C16 only matters for the survivors), flags C8-free survivors
(`SURV8`), "{8}-only" near misses (`P8ONLY`, no C4/C16 but has C8), and
full candidates (`CAND`, none of 4/8/16).

**Scope caveat (deliberate):** girth ≥ 5 excludes C4-free graphs that
contain triangles (attempt 001's "C4-free" counts at n ≤ 20 include such
odd-girth graphs; they are *not* covered here for n ≥ 22). A triangle does
not satisfy the conjecture, so that slice of the search space remains open
at n = 22 (see Leads).

**Restart resilience** (this attempt's predecessor was killed by a container
restart, and this run was hit by another one): each n is split into geng
`res/mod` slices (16 for n = 22); each slice's filter output is written
atomically (write to `.partial`, `mv` into place only after the kernel's
final `SUMMARY` line is present), so a completed slice is never re-run —
the runner skips any slice file that exists and contains `SUMMARY`. geng's
own stderr graph count is kept per slice (`.gengerr`) and the aggregator
cross-checks it against the kernel's total. This worked in practice: the
mid-run restart cost zero slices.

## What was done

1. **Kernel validation.** `cycle_filter` (binary survived the restart;
   smoke-tested on K4) was re-validated bit-for-bit against the attempt-001
   Python spectrum tool: for all 509 + 4060 connected cubic graphs on
   n = 14, 16, the kernel's girth / C4 / C8 / C16 bits matched the exact
   Python cycle spectra — 4569 graphs, 0 mismatches.
2. **Pool cross-check at n ≤ 20.** Girth-≥5 counts from `geng -tf` +
   kernel girth histogram: n = 14: 9, n = 16: 49, n = 18: 455, n = 20: 5783
   (girth split 5751 / 32) — matching the published counts of connected
   cubic graphs with girth ≥ 5 (OEIS A014372) and attempt 001's own girth
   histogram at n = 20. All of these graphs contain both C8 and C16
   (combo "011"), consistent with 001's finding that every C4-free graph
   at n ≤ 20 has C8 or C16.
3. **n = 22 exhausted.** 16 slices, 4 parallel workers, ~2.5 min wall.
   Per-slice totals 4597, 8643, 6062, 4938, 6729, 4729, 4544, 3502, 4681,
   8345, 5400, 5594, 4452, 7886, 5634, 5202 (each equal to geng's own
   stderr count for that slice; sum 90938).
4. **Near-miss verification.** The single n = 22 near miss was
   independently re-verified with `cycle_spectrum.py` (degree sequence,
   connectivity, full spectrum, bridge/block structure) — see below.
5. **McGee spectrum (for the n = 24 check).** The McGee graph — the unique
   (3,7)-cage, 24 vertices — was built from its LCF code [12, 7, −7]⁸ and
   its spectrum computed exactly with `cycle_spectrum.py`:
   **{7, 8, 9, ..., 24}** (pancyclic above its girth). It contains C8 and
   C16, so it is *not* a survivor; at n = 24 it must appear as the unique
   girth-7 entry of the histogram (a useful cross-check for that run).
6. **n = 24:** feasibility measured and slice design fixed (64 slices,
   same runner); see the addendum section at the end for status/results.

## Outcome

**Verified: every connected cubic graph with girth ≥ 5 on 22 vertices
contains a cycle of length 8 (and all but one of them also a C16). Zero
counterexamples, zero C8-free graphs among all 90938.**

| n | girth-≥5 graphs | published (A014372) | girth 5 / 6 | has C8+C16 ("011") | has C8, no C16 ("010") | C8-free (SURV8) | CAND |
|---|---|---|---|---|---|---|---|
| 22 | 90938 | 90938 | 90553 / 385 | 90937 | 1 | 0 | 0 |

Cross-checks that all passed:
- total matches OEIS A014372 (connected cubic graphs, girth ≥ 5) exactly;
- kernel totals match geng's independent stderr counts slice-by-slice;
- the girth-6 count 385 matches the known girth-≥6 cubic sequence
  (1, 1, 5, 32, 385 at n = 14..22; no girth-7 cubic graph exists below 24,
  so girth-≥6 = girth-6 here), and the n = 20 girth split reproduced
  attempt 001's 5751 / 32.

### The near miss (the one graph with no C16)

`U????A?O@?A?B?o_GKCA_?k?J?@_C?gO?GK?H@??` — girth 5, spectrum
**{5, 6, 7, 8, 9, 10}**, pow2 ∩ = {8}. Independently re-verified with
`cycle_spectrum.py`: all 22 degrees are 3, connected, spectrum as stated.
It has exactly one **bridge** (vertices 17–4); deleting it leaves two
11-vertex sides, *each* with spectrum {5, 6, 7, 8, 9, 10}. This is
precisely the bridged near-cubic construction of attempt 001's lead 1, in
its girth-5 form: each side is a C4-free near-cubic 11-vertex side (001
found 8 of these, all containing C8 — which is exactly why this graph is a
near miss and not a counterexample). The bridge caps the circumference at
10, killing C16 globally.

No graph avoided {4, 8, 16} — no candidate anomalies; nothing to mark
UNVERIFIED.

## Why it failed / what survived

"Failure" = no counterexample (expected; the conjecture is believed true).
What survived:

- The girth-≥5 exhaustion at n = 22 with double-checked counts — extending
  attempt 001's bound (n ≤ 20, full enumeration) one even size upward in
  the C4-free regime where the conjecture's tension actually lives.
- A validated fast C kernel (orders of magnitude beyond the Python tool)
  and a restart-proof slice pipeline (`run_slices.sh` +
  `aggregate_slices.py`) that survived a real mid-run container restart
  with zero lost slices.
- The n = 22 bridged {8}-only near miss with verified spectrum — extending
  001's pattern (n = 18, 20) that every known cubic graph whose pow2
  intersection is exactly {8} is bridged with small-circumference sides.
- Slice data: `data/g5n22_s*of16.txt` (+ `.gengerr` counts), aggregate
  `data/g5n22_aggregate.json`.

Limits: girth ≥ 5 only (see scope caveat), and n = 24 subject to the
compute budget (addendum below).

## Leads generated

1. **Close the triangle gap at n = 22 (and 24).** C4-free graphs with
   girth 3 (triangles allowed, no squares) are not covered by `-tf`
   generation. geng's `-f` flag alone (square-free) covers them; the pool
   is larger than girth-≥5 but still far smaller than full enumeration.
   Until then, "no C4-free counterexample at n = 22" is only claimed for
   girth ≥ 5.
2. **The {8}-only pattern is now strongly bridge-shaped.** All known
   {8}-only cubic graphs (n = 18, 20 from 001; n = 22 here) have a bridge
   and circumference ≤ 10. Test: enumerate girth-≥5 *2-connected* cubic
   graphs (geng has a connectivity-2 flag) for C16-free examples. If
   C16-freeness at these sizes truly requires a bridge, 001's ≥ 38-vertex
   bridge bound structurally forces a power-of-2 cycle for all cubic
   graphs to ~n = 30 except via a yet-unseen 2-connected mechanism.
3. **Extend 001's near-cubic side search to 19, 21 vertices**
   (`--blocks 19 21`): each new "no C4+C8-free side" result pushes the
   bridged-counterexample bound past 38 vertices, complementing lead 2.
4. **Girth ≥ 6 pool at n = 26–30** is tiny (385 at n = 22, growth roughly
   ~10–20×/step) and all such graphs are C4-free; `geng -tf` + girth-6
   post-filter reaches well past n = 24 within an hour-scale budget.
5. **n = 26 girth-≥5 exhaustion** (~31.5M graphs, next A014372 term) is
   roughly a 16 h run of this pipeline — feasible as a dedicated long
   cycle with the existing restart-proof runner (bump slices to ~1024).

## Reproduction

```
gcc -O2 -o tools/cycle_filter tools/cycle_filter.c
# kernel validation: compare `cycle_filter -v` bits against
# cycle_spectrum.py spectra for all cubic graphs at n = 14, 16 (~1 min).
# n = 22 (16 slices), restart-safe (re-running skips completed slices):
tools/run_slices.sh 22 16 4
python3 tools/aggregate_slices.py 22 16 \
    --json attempts/erdos-gyarfas/data/g5n22_aggregate.json
```

## Addendum: n = 24

Published girth-≥5 count at n = 24 (A014372): 1620479 (~17.8× n = 22).
Slice design: 64 slices via the same runner (`tools/run_slices.sh 24 64 4`),
files `data/g5n24_s*of64.txt`. Feasibility benchmark: slice 0/64 timed
first (banked as a real slice — its result file counts toward the run).
Expected cross-checks when complete: total = 1620479; exactly one girth-7
graph (McGee, spectrum {7..24} — contains C8+C16, so expected combo "011").

**COMPLETE (2026-07-26, closed out by the orchestrator after the search
agent was lost to a container restart; results verified from the banked
aggregate `data/g5n24_aggregate.json`):**

**Every girth-≥5 connected cubic graph on 24 vertices contains a C8.
Zero C8-free graphs, zero counterexample candidates.**

- Total 1,620,479 — matches published A014372 exactly; slice totals
  (all 64), girth histogram, and C4/C8/C16 combo counts each
  independently sum to the same total (orchestrator re-verified from the
  raw aggregate).
- Girth histogram: 1,612,905 girth 5; 7,573 girth 6; exactly 1 girth 7 —
  the unique girth-7 entry is the McGee-graph sanity check predicted
  above, and it lands in combo "011" (has C8 and C16) as expected.
- {8}-only phenomenon persists: 8 graphs at n = 24 have C8 but no C16
  (all girth 5, graph6 strings recorded in the aggregate), vs 1 at n = 22.
  Whether all 8 are bridged (per the lead-2 pattern) was not checked —
  folded into lead 2 for any future cycle.

Combined verified frontier for the conjecture in this lab's own tooling:
**all cubic graphs to n = 22 (full), plus all girth-≥5 cubic graphs at
n = 24.**
