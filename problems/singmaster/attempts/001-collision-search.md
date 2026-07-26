# Attempt 001 — Binomial-coefficient collision search (multiplicity ≥ 8)

**Problem:** Singmaster's conjecture. **Line:** computational search for entries
of Pascal's triangle with high multiplicity, up to a value bound B.
**Tool:** `tools/singmaster_search.py`. **Status:** COMPLETE (run finished
2026-07-25, 2066 s wall on 4 cores; results in `data/`).

## Approach

**Multiplicity convention (explicit).** mult(V) = #{(n,k) : 0 ≤ k ≤ n,
C(n,k) = V}, counting symmetric cells C(n,k) and C(n,n−k) separately and
counting the trivial cells C(V,1), C(V,V−1). Equivalently, for V ≥ 3:

    mult(V) = 2 + Σ over canonical cells (n,k), 2 ≤ k ≤ n/2, of (1 if n = 2k else 2)

Under this convention 3003 has multiplicity 8:
(3003,1),(3003,3002),(78,2),(78,76),(15,5),(15,10),(14,6),(14,8).

**Completeness argument.** For fixed k, C(n,k) is strictly increasing in
n ≥ k, so V has at most one canonical cell per k. One canonical cell gives
mult ≤ 4 (or 3 if central). Hence **mult(V) ≥ 5 requires ≥ 2 canonical cells
with distinct k**, at most one of which has k = 2. So every V ≤ B with
mult ≥ 5 has a canonical cell with k ≥ 3, and *all* of its canonical cells
have value V ≤ B. It therefore suffices to enumerate all canonical cells with
k ≥ 3 and value ≤ B, and for each value to test k = 2 membership in O(1)
(V triangular ⇔ 8V+1 a perfect square, via `math.isqrt`). The output census
of mult ≥ 5 values is thus **complete up to B**. Values with exactly one
canonical cell (mult 4 generically; mult 3 for central binomials) are generic
infinite families and are not listed individually.

**Data-structure / feasibility analysis.** Cell counts scale as
n_max(k) ≈ (k!·B)^(1/k): the k = 3 stream dominates with ~(6B)^(1/3) values
(≈ 1.1×10^10 at B = 2.5×10^29 — far too many to table), k = 4 has
~(24B)^(1/4) ≈ 5×10^7, k ≥ 5 only ~2×10^6 in total. So:

- **No k=3/k=4 tables at all.** Both are generated as sorted streams by the
  integer recurrences C(n+1,3) = C(n,3)+C(n,2), C(n+1,4) = C(n,4)+C(n,3)
  (pure big-int additions) and **two-pointer merged**, which catches
  C(a,3) = C(b,4) collisions exactly. This is the "per-k sorted arrays +
  merge" idea taken to its limit: the arrays are never materialized.
- **k ≥ 5 values go in a hash set** (~2×10^6 entries ≈ 200 MB), built once
  before forking workers (copy-on-write shared). Collisions inside k ≥ 5 and
  triangularity of k ≥ 5 values are detected at build time.
- Every streamed value is tested for (a) triangularity — with a two-level
  quadratic-residue prefilter (squares mod 693 and mod 4160; ~4.5% pass rate)
  before the isqrt, which took the hot loop from 2.0M to 2.9M values/s/core —
  and (b) membership in the k ≥ 5 set.
- **Parallelism:** the k=3 n-range is split into 4 equal chunks; chunk
  boundaries are value cutoffs V_i = C(n_i,3) so the k=4 sub-stream and all
  equality detection stay chunk-local.
- **Memory:** ~200 MB shared set + small per-worker state. Measured < 1.5 GB
  total for all 5 processes — well under the 4 GB budget.
- **Checkpointing:** each worker atomically writes {n3, n4, hits} every 30 s;
  a restarted run resumes from checkpoints (`--no-resume` to override). The
  k ≥ 5 table is regenerated deterministically in seconds rather than stored.

**Bound choice.** Benchmark: 2.9M values/s/core (Python 3.11 big-int hot
loop, values ~95 bits). With 4 cores and a 20–25 min budget this supports
N3 ≈ 1.15×10^10, i.e. **B = 2.5×10^29**. This deliberately exceeds
C(104,39) = 61218182743304701891431482520 ≈ 6.12×10^28, the first member of
the Fibonacci/Pell-family collision C(104,39) = C(103,40) beyond 3003's
C(15,5) = C(14,6) — so that family's next entry is *in* range and serves as a
deep validation target. (The task brief anticipated it being out of range;
the streaming design made it reachable.)

**Verification discipline.** The scan uses incremental addition chains; every
hit is re-verified in post-processing by an independent path: for each
k = 2..kmax, solve C(n,k) = V by binary search using `math.comb` (fresh
big-int products, no shared state with the scan), collect all canonical
cells, and recompute the multiplicity from scratch. A `--selftest` mode
checks all known facts, and `--brute N` is a fully independent row-by-row
brute force used to cross-check the pipeline below 10^7.

## What was done

1. Designed the search as above (convention, completeness proof, feasibility
   analysis) before writing code.
2. Implemented `tools/singmaster_search.py` (re-runnable; `--bound`,
   `--procs`, `--outdir`, `--no-resume`, `--selftest`, `--brute`).
3. Validated before the big run:
   - `--selftest` (B = 10^8): 3003 found with multiplicity exactly 8 and
     cells {(78,2),(15,5),(14,6)}; 120, 210, 1540, 7140, 11628, 24310 all
     found with multiplicity exactly 6 with the correct cells.
   - Independent brute force to 10^7 (direct row-by-row enumeration with
     `math.comb`): the set of values with mult ≥ 5 and their multiplicities
     agree exactly with the pipeline.
   - Intermediate runs at B = 3×10^15 and B = 10^20 (parallel path):
     no values with mult ≥ 5 beyond the classical seven.
4. Production run at B = 2.5×10^29 on 4 workers (~20 min), results and
   checkpoints in `attempts/singmaster/data/`.

## Outcome

**Bound reached: B = 2.5×10^29** (~4× beyond C(104,39) ≈ 6.12×10^28).
Scanned every canonical cell (n,k), 3 ≤ k ≤ n/2, with value ≤ B:
11,447,142,421 cells at k=3 (n ≤ 11,447,142,426), 49,492,314 at k=4
(n ≤ 49,492,321), 2,304,208 at k ∈ [5,50]. Runtime 2066 s (34 min) on
4 workers, peak memory ~1.0 GB. One container restart occurred during the
session; checkpoints and the harness-tracked waiter carried the run through
without re-scanning completed work.

**Complete census of multiplicities ≥ 5 for all values ≤ 2.5×10^29:**

| multiplicity | count | values |
|---|---|---|
| 8 | 1 | 3003 = C(78,2) = C(15,5) = C(14,6) |
| 7 | 0 | — |
| 6 | 7 | 120, 210, 1540, 7140, 11628, 24310, and 61218182743304701891431482520 = C(104,39) = C(103,40) |
| 5 | 0 | — |

- **No value beats multiplicity 8.** No mult ≥ 9, and no second mult-8 value,
  up to 2.5×10^29. Nothing UNVERIFIED to escalate.
- The six small mult-6 values are the known C(m,2) = C(n,k) coincidences:
  120=(16,2)=(10,3), 210=(21,2)=(10,4), 1540=(56,2)=(22,3),
  7140=(120,2)=(36,3), 11628=(153,2)=(19,5), 24310=(221,2)=(17,8).
- The search independently rediscovered the second member of the
  Fibonacci/Pell infinite family, C(104,39) = C(103,40) (found via the k ≥ 5
  set collision at k=40) — the deepest known-value validation available in
  range. The task brief expected this value to be out of range; the
  streaming design reached past it.
- **Multiplicity 5 and 7 are empty in range.** Both would require a central
  binomial C(2m,m) to coincide with another nontrivial cell; no such
  coincidence exists with value ≤ 2.5×10^29.
- Every reported hit was re-verified twice: once by the tool's
  post-processing (per-k binary search with `math.comb`, independent of the
  scan's incremental arithmetic), and once by a separate throwaway script
  re-deriving all cells and multiplicities from scratch. All match.
  Validation status: selftest PASSED, brute-force cross-check to 10^7
  identical, all known values reproduced with correct cells and
  multiplicities.

Run artifacts: `data/results_b327cb2734119d400000000000.json` (full census),
`data/run_b2.5e29.log`, `data/ck_b327cb2734119d400000000000_w{0..3}.json`
(worker checkpoints, all `done`).

## Why it failed / what survived

It did not find a multiplicity-≥9 entry — the expected outcome; this
extends the empirical support for Singmaster's N = 8 (indeed N = 8 is
attained only by 3003 up to 2.5×10^29). What survived:

- A re-runnable, checkpointed, validated search tool
  (`tools/singmaster_search.py`) whose census of mult ≥ 5 is *provably
  complete* up to the bound (argument in Approach), not merely heuristic.
- The streaming-merge design removes the memory wall entirely: the limit is
  purely CPU-time on the k=3 stream, which scales as B^(1/3). Pure Python
  caps this machine at roughly B ~ 10^30 per hour-class run.
- Empirical structure note: in range, *every* extra coincidence beyond the
  trivial pair involves either k=2 (six times) or the single Pell-family
  curve C(n,k)=C(n−1,k+1) (3003 and the C(104,39) value). High multiplicity
  in range is entirely explained by known infinite families plus the k=2
  sporadics.

## Leads generated

1. **C/u128 port** of the same design (k=3 stream + bit-filter + sorted k≥5
   array) would run ~20-50× faster, reaching B ~ 10^31–10^33 in the same
   wall time; values < 2^127 fit u128 up to B ~ 1.7×10^38. Natural next
   attempt if more depth is wanted; expected yield is still "no new value"
   (a genuinely new mult-8 needs a new integer point on some C(a,j)=C(b,k)
   curve, j<k — see lead 2).
2. **The equation table** (problem file's second attack line) is the better
   use of cycles: mult ≥ 8 up to any bound is governed by simultaneous
   solutions of pairs among {C(n,2)=C(m,k)}, {C(n,3)=C(m,4)}, etc. Several
   of these are individually resolved in the literature (e.g. Avanesov for
   C(n,2)=C(m,3); de Weger's list of all C(a,j)=C(b,k) with value ≤ ~10^30
   matches our census exactly, which is an external cross-check of this
   run). Building the resolved/open table would convert "search deeper"
   into "which finitely many curves could still produce a 9th cell".
3. The mult-7 emptiness observation: a mult-7 value needs a central binomial
   C(2m,m) equal to a non-central cell. C(2m,m) grows like 4^m, so only
   ~50 central values exist below any practical bound — a *targeted* search
   testing central binomials against all k up to huge bounds (say 10^1000,
   via per-k inversion, no enumeration) is nearly free and would be a cheap
   follow-up with a clean statement ("no central binomial below 10^X is any
   other binomial").
4. The residue prefilter + streamed-merge pattern is reusable for any
   "collisions among sparse polynomial sequences" search in this lab.

