# 001 — Census of centrally symmetric {0,±1}^4-vertex polytopes vs the Mahler bound 32/3

- **Problem:** mahler-4d
- **Date:** 2026-07-30
- **Mode:** blind
- **source-commit:** ac7a82245ad734c0e190f2aa49d6316e80935559
- **Type:** exhaustive computational census (exact-integer combinatorics + float filter in C, exact-rational verification layer in Python over the tier-0 harness)
- **Tools:** `cc -O2` (int64 exact combinatorics, `problems/mahler-4d/explore/census.c`), Python 3 stdlib (`fractions`), tier-0 harness `harness/mahler-4d/polytope.py` (divergence-theorem volumes) and `harness/mahler-4d/verify_product.py` (independent Fubini slicing)
- **Sources:** tier-0 `problems/mahler-4d/PROBLEM.md` and harness only (blind working copy); no attempt records read

## Approach

Question: in an exact census of centrally symmetric polytopes with vertices in
{0,±1}^4, is anything within ε of the conjectured minimum volume product
32/3 that is not a Hanner polytope?

**Generator universe, stated precisely.** The nonzero points of {0,±1}^4 form
40 antipodal pairs. A centrally symmetric polytope with vertices in {0,±1}^4
is exactly `conv(±S)` for a set `S` of `k` pairs whose 2k points are all in
convex position ("proper") and span R^4 ("rank 4"); 0 is then automatically
interior. The census enumerates **every subset of pairs with 1 ≤ k ≤ 11, up to
the hyperoctahedral group B4** (signed coordinate permutations, |B4| = 384,
acting on pairs through its faithful quotient B4/{±I} of order 192), and
computes the volume product of every proper rank-4 orbit representative with
4 ≤ k ≤ 11. Equivalently: **every centrally symmetric polytope in R^4 with at
most 22 vertices, all in {0,±1}^4, up to signed coordinate permutation** —
which is a linear map, so the volume product is constant on each orbit and the
reduction loses nothing. The full universe of pair-subsets is 2^40; the
symmetry- and level-capped census is 18,637,214 orbits (k ≤ 11), of which
18,637,066 have k ≥ 4. Every Hanner polytope in R^4 has 8–16 vertices
(k = 4–8), so the cap k ≤ 11 contains all of them with three levels of
headroom. What is NOT covered: symmetric {0,±1}^4-polytopes with 12–20 pairs
(24–40 vertices, e.g. the 24-cell at k = 12), and any polytope whose vertices
need coordinates outside {0,±1}.

**Why this design rather than the obvious alternatives.**
- *Naive exhaustion* of 2^40 subsets is out of reach; *random sampling* would
  produce no census claim at all. Orbit enumeration by levels gives an exact,
  auditable scope line ("everything with ≤ 22 vertices").
- Enumeration is by **canonical forms** (minimum of the 40-bit mask over the
  384 group images) with level-by-level augmentation. Augmentation is
  complete: any (k+1)-subset minus an element is B4-equivalent to a canonical
  k-representative, so extending every canonical k-representative by every
  pair and re-canonicalizing reaches every (k+1)-orbit.
- A mask whose hull drops a point ("improper") is skipped, because its
  polytope already appears at its own vertex count; a mask of rank < 4 is
  degenerate (not a body). This is what makes "orbits evaluated" equal
  "polytopes counted once each".
- Arithmetic is split into two layers so that speed never touches the claim:
  the C kernel makes **every combinatorial decision (facet identification,
  sidedness, vertex tests, incidence) in exact int64 integer arithmetic**
  with overflow guards that downgrade a body to `punt` (= must be recomputed
  exactly) rather than ever approximating; only the final volume assembly is
  double (relative error ~1e-13, vs. a spectral gap of ~0.15). Every body
  with float P < 10.9 — and independently, stratified samples — was then
  recomputed **exactly over Q** with the tier-0 harness, and the key bodies
  were additionally re-verified by the harness's independent Fubini-slicing
  implementation.

**Volume computation in the kernel** (same decomposition as the harness,
reimplemented independently in C): divergence theorem recursing on dimension,
`vol_d = (1/d) Σ_F b_F · vol_{d-1}(π_j F)/|a_{F,j}|` with integer-primitive
outward normals from brute-force d-subset enumeration. The polar body needs
no second vertex enumeration: the vertices of K° are `a_F/b_F` over facets F
of K, and the facet of K° supported by `⟨v,y⟩ = 1` (v a vertex of K) has
vertex set `{a_F/b_F : F ∋ v}`, scaled to integers by `L = lcm(b_F)` (guarded;
`punt` on overflow — zero punts occurred in the entire run).

## What was done

All commands run from `problems/mahler-4d/` in this working copy;
`census` = `explore/census` built with `cc -O2 -o explore/census
explore/census.c -lm`. Everything is deterministic and restartable
(slice files are written atomically and skipped when complete).

1. **Harness self-tests** (established baseline, not re-litigated):
   `python3 harness/mahler-4d/polytope.py --selftest` and
   `python3 harness/mahler-4d/verify_product.py --selftest` both pass;
   the harness generates 22 Hanner vertex sets in R^4 (8–16 vertices),
   all with P = 32/3 exactly.

2. **Orbit enumeration.** `./explore/census seed data/level01.masks`, then
   `./explore/census expand data/level{k}.masks data/level{k+1}.masks` for
   k = 1..10. Orbit counts per level:

   | k | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 |
   |---|---|---|---|---|---|---|---|---|---|----|----|
   | orbits | 4 | 21 | 123 | 756 | 4310 | 22567 | 103649 | 415920 | 1456466 | 4478554 | 12154844 |

   **Cross-check:** `python3 explore/burnside_check.py 11` counts orbits of
   k-subsets by Burnside's lemma over the 384 group elements — an independent
   Python implementation sharing no code with the C canonicalization. All
   eleven counts match exactly. (It also verifies the action on pairs has
   kernel exactly {±I}.)

3. **Census evaluation.** `./explore/run_census.sh 4 9 200000` and
   `./explore/run_census.sh 10 11 200000` evaluate every orbit representative
   (checkpointed slices `data/eval{kk}_{start}.txt`, each slice < 2 min).
   `python3 explore/aggregate.py 4 11 --near 10.9 --json data/summary.json`
   verifies completeness (every slice present, expected line counts, masks in
   level-file order) and aggregates. Per-level results:

   | k | orbits | proper ("ok") | improper | degenerate | punt | min P (float) |
   |---|--------|----|----------|-------|------|---------------|
   | 4 | 756 | 517 | 0 | 239 | 0 | 10.666666667 |
   | 5 | 4310 | 3423 | 538 | 349 | 0 | 10.666666667 |
   | 6 | 22567 | 15850 | 6298 | 419 | 0 | 10.666666667 |
   | 7 | 103649 | 56520 | 46731 | 398 | 0 | 11.000000000 |
   | 8 | 415920 | 156078 | 259545 | 297 | 0 | 10.666666667 |
   | 9 | 1456466 | 330941 | 1125348 | 177 | 0 | 11.097222222 |
   | 10 | 4478554 | 538861 | 3939612 | 81 | 0 | 11.083333333 |
   | 11 | 12154844 | 671525 | 11483291 | 28 | 0 | 11.053240741 |

   (Levels 10–11 were evaluated by the same `census eval` slices, driven
   partly by `explore/run_census.sh 10 11 200000` and partly by three
   concurrent `explore/eval_worker.sh 11 <first> <last>` workers on disjoint
   slice ranges; slice content is deterministic and every slice was
   completeness-checked by `aggregate.py`, so the parallelism cannot change
   any result.)

   Zero `punt`s: no overflow guard ever tripped, so the exact-integer
   combinatorics covered every body.

4. **Exact verification of everything near the bound.** All 1176 orbits with
   float P < 10.9 (all at k ≤ 8; there are zero with float P < 10.9 at
   k = 7, 9, 10, 11) were recomputed exactly over Q with the harness
   (`python3 explore/exact_check.py file ...`):
   - **1113 orbits have P = 32/3 exactly** (517 at k=4, 414 at k=5, 163 at
     k=6, 19 at k=8; none at k=7, 9, 10, 11);
   - **63 orbits lie strictly above**, the smallest being
     P = 336311/31104 = 32/3 + 4535/31104 ≈ 10.8124678 (one orbit, mask
     `81101001`, 10 vertices, 22 facets);
   - **nothing came out below 32/3.**

5. **Every attainer is a Hanner polytope.** `python3
   explore/hanner_match.py matchfile ...` classifies each of the 1113
   attainers: 10 are B4-equivalent to a harness-generated Hanner vertex set
   outright, and the remaining 1103 have an **explicit invertible rational
   matrix T with T(vertices) = (harness Hanner vertices), found by exhaustive
   search over images of a fixed linear basis of antipodal pairs** (T maps
   pairs to pairs; ±T identified; invertibility and set-bijection both
   checked — an earlier version omitted the invertibility check and is
   superseded). Zero attainers are NOT-HANNER. Certificate maps are in
   `data/hanner_match_attainers.txt`.

6. **Adversarial second implementation on key bodies.** Claim records
   (`data/claims/claim_*.json`) for the nearest non-attainer `81101001`
   (P = 336311/31104), the k=7 minimum `16b40` (P = 11 exactly), the k=6
   runner-up `1680180` (P = 98/9), and attainers `4c1000`, `6016`, `6036`,
   `361b0` (P = 32/3) were re-verified by
   `python3 harness/mahler-4d/verify_product.py --check data/claims/claim_X.json`
   — Fubini slicing with exact polynomial interpolation, plus the duality
   identities. **All CONFIRMED**, none escalated.

7. **Stratified float-vs-exact samples** (guarding against a systematic
   kernel error away from the near-bound region):
   `python3 explore/exact_check.py sample data/eval07_00000000.txt 500`
   (114 bodies, k=7), `... data/eval08_00200000.txt 3000` (26, k=8),
   `... data/eval09_00600000.txt 20000` (4, k=9),
   `... data/eval10_00000000.txt 8000` and `... data/eval10_00400000.txt
   2500` (11, k=10), `... data/eval11_00000000.txt 4000` (1, k=11):
   **156 bodies, 0 mismatches > 1e-9**.

8. **Structural cross-checks.**
   - The canonical masks of all 22 harness Hanner vertex sets (10 B4-orbits)
     appear in the census attainer lists at their vertex counts; every k=4
     proper orbit (necessarily a linear image of the cross-polytope) has
     P = 32/3 exactly, as it must.
   - **Brute-force k=4 re-enumeration**: all C(40,4) = 91,390 subsets
     canonicalized by the independent Python group action give 756 orbits,
     517 of rank 4 — matching the C expand/eval pipeline exactly.
   - **Improper bookkeeping spot-check**: for sampled improper masks, the
     exact P of the redundant point set equals the exact P of its extreme
     submask, whose canonical form is present in its own level file (e.g.
     `207f` → `205a` at k=5, P = 32/3; `20bf` → `20ba` at k=6, P = 221/18).

## Outcome

**EVIDENCE.** Scope: every centrally symmetric polytope in R^4 with at most
22 vertices, all vertices in {0,±1}^4 — equivalently all `conv(±S)`, S a set
of ≤ 11 of the 40 antipodal pairs of {0,±1}^4∖{0} — enumerated exhaustively
up to signed coordinate permutations (a volume-product-preserving reduction).
Of the 18,637,066 orbit representatives evaluated (4 ≤ k ≤ 11), **1,773,715
are proper rank-4 bodies — distinct centrally symmetric 4-polytopes, each
counted exactly once**; 1,988 are degenerate (rank < 4) and 16,861,363 are
improper (their hulls appear at their true vertex count). Zero overflow
punts: exact integer combinatorics covered every body.

Within that universe:

1. **min P = 32/3 exactly; nothing lies below.** (Float filter over all
   orbits; everything with float P < 10.9 re-proved exactly in Q; float error
   bound ~1e-11 ≪ the 0.146 gap.) **VERIFIED (exact, in the harness's sense)
   for the 1176 near-bound bodies individually; EVIDENCE at census scope.**
2. **The attainment set is exactly the Hanner polytopes**: 1113 orbits attain
   32/3 and each carries an explicit linear equivalence to a harness Hanner
   polytope. No non-Hanner body attains the bound in this universe.
3. **Spectral gap:** the smallest volume product of a non-Hanner body is
   336311/31104 = 32/3 + 4535/31104 ≈ 32/3 + 0.1458 (a single 10-vertex
   orbit). So for every ε < 4535/31104 the answer to the census question is
   **no — nothing within ε of 32/3 except Hanner polytopes.** Runner-up
   values: 29165/2688 (k=5, two orbits), 6251/576 (k=5, two), 56399/5184
   (k=5, four), 98/9 (k=6, 54 orbits), 85211/7776 (k=6, one),
   299/27 (k=8 runner-up); the k=7 level minimum is exactly 11
   and the k=9 level minimum exactly 799/72, the k=10 minimum exactly
   133/12 (mask `32c028b000`), the k=11 minimum exactly 4775/432 (mask
   `72c028b000`) — all exact values recomputed in Q.
4. Levels without Hanner-compatible vertex counts (k = 7, 9, 10, 11) have
   minima exactly 11, 799/72, 133/12, 4775/432 respectively — all ≥ 11 =
   32/3 + 1/3: the near-bound region is populated ONLY at Hanner vertex
   counts, and the level minima beyond k = 8 hover just above 11 without
   a monotone trend (799/72 ≈ 11.097 > 133/12 ≈ 11.083 > 4775/432 ≈ 11.053).

**NOT claimed:** anything about polytopes with more than 22 vertices (the
24-cell at 24 vertices is outside; its P = 16 comes from the harness
self-test, not this census); anything about vertex coordinates outside
{0,±1}; anything about general convex bodies; any support for the conjecture
beyond this finite universe. A refutation of Mahler's conjecture was neither
found nor expected here — a counterexample, if it exists, is known (published
local-minimality results) not to be near a Hanner polytope, and this universe
is heavily Hanner-adjacent.

## Why it failed / what survived

This was a census, not a proof attempt; it "fails" in the sense that no
new mathematical object appeared: the universe is exhausted and the bound is
attained only by Hanner polytopes, exactly as conjectured. Specific
obstructions and what remains usable:

- **The obstruction to going further out**: orbit counts grow ~×3 per level
  (k=12 is ~29.7M orbits, ~2.4× the entire k≤11 census) while the interesting
  region empties out — the k = 9, 10, 11 minima sit at 799/72, 133/12,
  4775/432, i.e. slightly DECREASING toward 11 but stuck a full 1/3 above
  the bound, and the near-bound region is empty away from Hanner vertex
  counts. Widening the coordinate set to
  {0,±1,±2}^4 multiplies the pair universe from 40 to 312 and kills this
  enumeration strategy outright; a different reduction (e.g. by facet count
  or by fixing combinatorial type) would be needed.
- **Survived, reusable**: `explore/census.c` (exact-integer 4D polytope
  kernel: facets, vertex tests, volume, polar volume, ~0.2–0.5 ms/body;
  the guard/`punt` mechanism never fired but is load-bearing for trust);
  the canonical-augmentation enumerator with the Burnside cross-check
  pattern; `explore/hanner_match.py`'s pair-basis GL-equivalence search
  (note the singular-map pitfall it now guards against: membership
  `T(V) ⊆ W` is NOT equivalence — an early buggy run produced singular
  "equivalences" and was caught because printed matrices had repeated rows
  up to sign; the fix checks rank and set-bijection).
- **A negative worth keeping**: the float layer never disagreed with exact
  arithmetic (0 mismatches in 156 stratified samples + 1176 near-bound
  bodies + 8 exact level minima/spectrum points), and no int64 guard
  tripped at k ≤ 11 — the {0,±1} lattice is
  numerically tame at this scale, so future censuses here can trust the same
  two-layer design without a wider rational kernel.

## Leads generated

1. **Close the {0,±1}^4 universe completely: k = 12..20.** ~29.7M orbits at
   k=12 and shrinking headroom after that (the total is bounded by orbits of
   all subsets, ~2^40/192 ≈ 5.7e9, but proper masks die out fast — improper
   fraction is already 77% at k=9). Falsifiable claim to test: *no proper
   mask with k ≥ 12 has P < 11*. The kernel as-is can do k=12 in a few
   hours; k=13+ needs either the improper-detection shortcut (test convex
   position before full facet enumeration) or a C hash-free streaming
   canonicalizer.
2. **The gap body `81101001`** (vertices ±{(0,0,0,1),(0,1,1,1),(1,−1,1,0),
   (1,0,−1,1),(1,1,−1,−1)}, P = 336311/31104): check whether it is a linear
   image of a known near-minimizer family; its combinatorial type is
   simplicial (10 vertices, 22 facets, every facet a simplex; vol = 73/12,
   vol° = 4607/2592). Falsifiable: is it combinatorially a "cross-polytope
   with one split vertex pair"? Its uniqueness at its value (single orbit
   vs. 2–4 orbits for the next spectrum values) is curious; note also that
   the k=6 third-lowest body `81101101` (P = 85211/7776) is this mask plus
   one more pair.
3. **Spectrum rigidity conjecture (SPECULATION, falsifiable):** in this
   universe every volume product is rational with small denominator and the
   value set near the bound is very sparse ({32/3} ∪ [336311/31104, ...]).
   Does the minimum non-Hanner value over {0,±1}^n lattices grow or shrink
   with n? Compare with the known n=3 spectrum by running this exact
   pipeline on {0,±1}^3 (trivially cheap, 13 pairs) — if the n=3 gap is
   larger, the lattice census is getting *denser* near the bound as n grows,
   which would say something about where a counterexample search should live.
4. **{0,±1,±2}^4 spot-widening.** Not the full census (312 pairs), but the
   targeted sub-universe of masks obtained from the 63 near-bound
   non-attainers and 1113 attainers by replacing one pair with a weight-2
   point — a local search in the lattice around the near-minimizers.
   Falsifiable: *no such one-step widening goes below 98/9*.
5. **Equality-case data for the equality characterization.** The 1113
   attainer certificates (explicit T per orbit) are a dataset: how many
   GL(4,Q)-images of each Hanner combinatorial type embed in {0,±1}^4?
   (Counts by target: see `data/hanner_match_attainers.txt`.) A closed-form
   count would be a small standalone lattice-geometry result and a
   consistency check on any future claimed equality-case proof.

## References

- Tier-0 `problems/mahler-4d/PROBLEM.md` (statement, published status:
  Mahler 1939; Iriyeh–Shibata for n=3; Nazarov–Petrov–Ryabogin–Zvavitch
  2010 and J. Kim arXiv:1212.2544 on local minimality near Hanner
  polytopes; Bourgain–Milman 1987; G. Kuperberg GAFA 2008).
- Tier-0 harness: `harness/mahler-4d/polytope.py`,
  `harness/mahler-4d/verify_product.py` (this attempt's exact reference
  implementations and the independent verifier).
- This attempt's tooling: `problems/mahler-4d/explore/census.c`,
  `explore/run_census.sh`, `explore/aggregate.py`, `explore/exact_check.py`,
  `explore/hanner_match.py`, `explore/burnside_check.py`; data under
  `problems/mahler-4d/data/`.

---

*Filing note (orchestrator): this record was written in the blind working copy
as `ATTEMPT.md` and filed verbatim apart from this note and the numbered
title. The checkpointed evaluation slices (`data/eval{04..11}_*.txt`, ~1.05 GB)
and the orbit mask lists (`data/level*.masks`, ~155 MB) are NOT committed to
the repository for size reasons; they are deterministic and regenerable with
the commands recorded above (`explore/run_census.sh`, `explore/eval_worker.sh`
— total ~2 h wall on 4 cores). Every derived file the claims rest on
(candidates, exact values, certificates, claims, summaries) is committed.*
