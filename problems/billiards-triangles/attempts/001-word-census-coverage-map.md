# 001 — Word census by length for stable periodic orbits in obtuse triangles

- **Problem:** billiards-triangles
- **Date:** 2026-07-30
- **Mode:** blind
- **source-commit:** ac7a82245ad734c0e190f2aa49d6316e80935559
- **Type:** computation — enumeration + certified coverage map + growth measurement
- **Tools:** `harness/billiards-triangles/unfold.py` (exact certificates), `harness/billiards-triangles/verify_cover.py` (independent simulation re-verification), new `problems/billiards-triangles/explore/census.py` (standard-library Python 3.11.15)
- **Sources:** tier-0 working copy only (`PROBLEM.md`, harness, `CONTRIBUTING.md`, `AGENTS.md`)

**SCOPE LIMIT, stated up front:** everything below is about *open-region*
(unfolding-corridor) certificates, which can only see periodic orbits that
are stable under perturbation of the triangle. This attempt says **nothing**
about unstable/isolated periodic orbits: a triangle in the "uncovered" set
may well have a periodic orbit that no open-region certificate can see.

## Approach

The conjecture, restricted to the constructive side, is a covering problem:
each bounce word whose unfolding closes to a translation certifies an open
tile of the parameter triangle, and the question is whether some word family
tiles the obtuse region. Rather than hunting single deep orbits (the obvious
alternative), I ran a **census by word length**: enumerate *every* candidate
word up to length L, map which parts of the obtuse region each length
certifies, and measure where coverage stalls as the largest angle grows.
The point of census-over-hunt is that negative space is a deliverable: the
geometry of the uncovered set and the growth rate of the minimal certificate
length are exactly what decides whether the finite-census route is viable at
all (the stated kill condition).

Parameter space: every obtuse triangle is similar to exactly one with the
obtuse vertex at C, A=(0,0), B=(1,0), apex C=(x,y) in the open half-disk
(x−1/2)²+y²<1/4, y>0, and (up to the A↔B reflection, i.e. congruence)
x ≤ 1/2. The census region R is that half of the half-disk, area π/16. The
largest angle γ (at C) is 90° on the arc and → 180° as y → 0.

**Word universe and pruning (exact rules).** U(L) = words over sides {0,1,2}
with: (W1) no two cyclically adjacent letters equal; (W2) even length (odd
unfoldings are never translations); (W3) unfolding linear part is a rotation
with p=q=0 and r even in p·α+q·β+r·π — the harness's integer translation
condition, re-implemented incrementally for DFS and cross-checked against
`unfold.is_translation_word` on all words of length ≤ 10 (census selftest);
(W4) canonical under cyclic rotation and reversal — sound because a rotated
word's gate projections onto the corridor normal are the same set (gate n+1
is gate 1 translated by τ, and τ projects to 0), and a reversed word's chain
is congruent; (W5) w = u^k (k≥2) dropped when u is itself a translation word
(identical corridor), but kept when u is not (Fagnano = (012)² survives as
it must). The 0↔1 letter swap is *not* quotiented: it maps a tile to its
x→1−x mirror, and R is already the x ≤ 1/2 half. One measured subtlety: on a
positive-width box the harness's affine *enclosure* is order-dependent, so
rotations of the same word can differ TRUE vs UNKNOWN; the certifier
therefore retries up to 4 rotation/reversal variants per candidate
(`variants` in census.py; selftest checks no TRUE/FALSE contradiction ever
occurs and that all variants agree at zero-width boxes).

**Two-layer census (and why).** The harness enclosure of the translation
vector blows up with box width — measured: for a length-14 word the verdict
at a fixed apex goes TRUE only at box size ≈ 2⁻¹², with enclosure width ~470
at box 1/64 — and the certifiable box size shrinks ≈ 2 bits per letter
(TRUE at half-width 2⁻¹³ for length 14, 2⁻¹⁹ for 20, 2⁻²⁵ for 26). So
whole-cell certification at map resolution is impossible for the lengths that
matter, and the census has two layers:
- **layer 1 (map, float, never claimed):** corridor sign of every canonical
  word at the midpoint of each 1/64 grid cell, lengths ascending → minimal
  *candidate* length per cell. A word whose midpoint corridor is empty
  cannot certify any box containing the midpoint, so this is a sound screen
  for layer 2 up to float rounding.
- **layer 2 (anchor, exact):** for each candidate cell, a harness
  certificate in exact arithmetic on a small centred dyadic box, starting at
  half-width 2^-(len+2) (retry 2^-(len+6)) and doubling while TRUE, cap
  1/128. These anchor boxes are the only certified regions claimed, and
  every one was then attacked by `verify_cover.py`'s independent rational
  simulation.

## What was done

All commands run from the working-copy root with
`MATHLAB_OUT=$PWD/out`; every stage checkpoints to JSONL and is re-runnable
and restartable (completed cells/arcs are skipped on re-run). Validation
first: `python3 harness/billiards-triangles/unfold.py --selftest`,
`python3 harness/billiards-triangles/verify_cover.py --selftest`, and
`python3 problems/billiards-triangles/explore/census.py selftest` all pass.

1. **Universe.** `census.py words --max-len 26 --out
   problems/billiards-triangles/data/words` (~3.3 min). Canonical counts by
   length: 6:1, 10:3, 12:4, 14:15, 16:30, 18:109, 20:318, 22:1101, 24:3544,
   26:12402 — total 17,527; growth ratio ≈ ×3.2–3.5 per +2 letters.

2. **Layer-1 map** (grid 32, cell 1/64, 739 non-boundary cells + 94
   boundary): `census.py map --words .../data/words --max-len 26 --grid 32
   --row-lo 0 --row-hi 16 --out .../data/cover/map_g32.jsonl` and the same
   with `--row-lo 16 --row-hi 32` (~2 min total). Result: 264 cells with a
   candidate ≤ 26, 475 with none.

3. **Layer-2 anchors.** `census.py anchor --words .../data/words --max-len
   26 --src .../data/cover/map_g32.jsonl --row-lo 0 --row-hi 20 --out
   .../data/cover/anchor_g32.jsonl` and `--row-lo 20 --row-hi 32` (~1 min
   total): **260 of 264 candidate cells anchored** with an exact TRUE
   certificate (4 failures at cand_len 22/26, thin corridors under the
   exact_cap=4 / two-start-exponent budget).

4. **Adversarial verification.** `census.py simcheck --src
   .../data/cover/anchor_g32.jsonl --full 5 --out
   .../data/verify/simcheck_anchors.jsonl`: all 260 anchor orbits re-derived
   by exact rational billiard simulation at the anchor centre (independent
   code path, no unfolding), 52 of them additionally at 9 sample apexes
   across the whole certified box — **0 contradictions**. Same for all 107
   arc-sweep anchors below (`simcheck --src <4 arcsweep files> --full 7`,
   0 contradictions).

5. **Growth measurement (kill condition).** `census.py arcsweep` samples
   16–32 apexes per constant-γ arc and records the minimal candidate length,
   anchoring the best point exactly:
   - 90.5°–130° step 0.5° (16 pts/arc), 112.05°–112.45° step 0.05°
     (32 pts), 90.05°–90.45° step 0.05° (32 pts), 130.5°–150° step 0.5°
     (16 pts). Files under `data/growth/arcsweep_*.jsonl`.
   - `census.py bisect_death --len {14,18,22,26}` bisects the γ at which the
     last length-n candidate dies (float, 400 pts/arc, full-universe confirm
     scan at death+0.25°).
   - `census.py pointsweep` re-checks frontiers in **exact arithmetic**: all
     15 canonical length-14 words at 80 exact rational apexes on the
     γ≈112.6° arc → 0 TRUE (and 1 TRUE at 112.4°); all 109 length-18 words
     at 60 apexes at γ≈120.3° → 0 TRUE (7 TRUE at 119.7°).

6. **Family extraction and out-of-universe tests.** The last-surviving word
   at every length is one family: W(a,b) = (0(12)^a(02)^b)², length 4a+4b+2
   ... specifically W(2,1) len 14, W(2,2) len 18, W(3,2) len 22, W(3,3)
   len 26 — and W(1,0) is Fagnano. `census.py family --a A --b B` scans,
   bisects birth/death, and exactly anchors each member; run for (2,1),
   (3,2), (3,3), (4,3), (4,4), (5,5), (6,6), (8,8).

7. **Aggregation.** `census.py summarize --map ... --anchor ... --out
   .../data/summary/coverage_map.json` (machine-readable coverage map) and
   `census.py ascii ... --out .../data/summary/ascii_map.txt`.

### Headline numbers (all reproducible by the commands above)

Coverage of the interior cells (area-weighted; interior area 0.1804 of the
π/16 ≈ 0.1963 region; 0.0159 unresolved boundary fringe at this grid):

| max word length | candidate fraction | anchored fraction |
|---|---|---|
| 14 | 7.6% | 7.2% |
| 18 | 15.7% | 14.9% |
| 20 | 22.7% | 21.9% |
| 22 | 29.9% | 29.0% |
| 24 | 30.7% | 30.0% |
| 26 | 35.7% | 35.2% |

**64.3% of the obtuse region has no stable-orbit word of length ≤ 26 at its
cell midpoint.** By γ band (candidate fraction): 92–95°: 54%, 95–100°: 76%,
100–105°: 61%, 105–112.3°: 59%, 112.3–120°: 60%, 120–135°: 39%,
**135–180°: 0%**.

Minimal candidate length along each γ arc (staircase, at the stated
sampling): ℓ=20 on (90°, 92°) [but see family note below], ℓ=14 on
[92°, 112.5°), ℓ=18 on [112.5°, 120°), ℓ=22 on [120°, 130°), ℓ=26 on
[130°, 135°), **nothing ≤ 26 at any sampled point for γ ≥ 135°**.

Measured family deaths (float bisection, 400 arc points, ±~0.005°
sampling noise, deaths only ever *under*-estimated by finite sampling):

| word | len | birth | death | nice value |
|---|---|---|---|---|
| W(2,1) | 14 | 90.017 (K=1600) | 112.4989 | 112.5 = 5π/8 |
| W(2,2) | 18 | — | 119.9991 | 120 = 2π/3 |
| W(3,2) | 22 | 120.0595 | 129.9881 | 130 = 13π/18 |
| W(3,3) | 26 | 127.5064 | 134.9941 | 135 = 3π/4 |
| W(4,3) | 30 | 135.0486 | 140.6161 | 140.625 |
| W(4,4) | 34 | 139.5134 | 143.9964 | 144 |
| W(5,5) | 42 | 147.0320 | 149.9963 | 150 |
| W(6,6) | 50 | 152.1523 | 154.2842 | 154.2857 |
| W(8,8) | 66 | 158.7699 | 159.9981 | 160 |

SPECULATION (fits the nine data points above to ~10⁻²–10⁻³ degrees, no
proof): death(W(a,a)) = 180·a/(a+1)°, death(W(a,a−1)) =
180·(1 − (2a−1)/(2a²))°; W(1,0) = Fagnano fits with death exactly 90°.
The a=4..8 rows were *predictions* made from the a≤3 rows and then measured.

Every family member was **exactly anchored and independently simulated**,
e.g. W(8,8) (length 66) certifies the box of half-width 2⁻⁵⁷ centred at
apex (15381/32768, 2989/32768), obtuse angle γ ≈ 159.25° — re-derived by
exact rational simulation, closing after 33 bounces (orientation-reversing
return, period = half the word, like Fagnano's 3-of-6).

**Growth rate (kill-condition verdict).** Along the *best* point of each
arc, minimal stable word length grows only ~linearly in 1/(180°−γ):
ℓ(W(a,a)) = 8a+2 with death 180a/(a+1)° gives ℓ(γ) ≈ 1440/(180−γ) − 6.
That is mild — NOT super-exponential. What actually kills census-by-length
is measured elsewhere: (i) the canonical universe grows ×~3.4 per +2
letters, so *finding* the needles at γ=150° (ℓ≈42) means ~10¹¹ canonical
words; (ii) the tiles collapse: the certifiable anchor half-width shrinks
≈2 bits per letter (2⁻¹³ at ℓ=14 → 2⁻⁵⁷ at ℓ=66; part harness enclosure,
but the float tile x-width visibly collapses too — W(2,1)'s alive window at
γ=91° is x∈[0.491, 0.5) and pinches to the right-isoceles corner as γ→90⁺),
so certified area per word shrinks exponentially in ℓ even though ℓ(γ)
doesn't blow up.

### Uncovered-set geometry (the MAP deliverable)

From `data/summary/coverage_map.json` + `ascii_map.txt`:
1. **Everything with γ > 135°** — no candidate ≤ 26 anywhere sampled; the
   family curve continues (certified anchors to 159.25°) but as isolated
   thin tiles hugging a curve, not area.
2. **A band along the right-angle arc**: the 91–92° band has candidate
   fraction 0 at cell-midpoint resolution — near-90° tiles exist (W(2,1)
   down to γ=90.017°) but only in slivers near x=1/2 thinner than a cell.
3. **Vertical funnels below arc points** at x ≈ 0.27–0.40, including
   directly below the 30-60-90 apex (1/4, √3/4): a persistent no-candidate
   pocket at all lengths ≤ 26 (the published account of tiles failing to
   cover a neighbourhood of (π/6, π/3) is consistent with this; found here
   blind).
4. **Pinch gaps at the nice angles**: at γ = 120° and 135° (=180a/(a+1))
   the covering family dies exactly there while its successor is born
   *above* it — measured gap [135.000°, 135.049°] between death(W(3,3)) and
   birth(W(4,3)) in which nothing ≤ 30 is known alive. The frontier angles
   of the census are accumulation points, not soft edges.
5. Scattered mid-region holes at every length (see `ascii_map.txt` '#'
   cells interleaved with certified cells at 92–135°).

## Outcome

- **VERIFIED** (each independently re-derived by exact rational simulation,
  which is this problem's verification contract): 260 anchor certificates
  from the map + 107 from arc sweeps + 9 family anchors — each an exact
  statement "word w certifies the axis-aligned box of stated dyadic centre
  and half-width", including a length-66 stable orbit at γ ≈ 159.25°. Exact
  boxes, words and half-widths in `data/cover/anchor_g32.jsonl`,
  `data/growth/arcsweep_*.jsonl`, `data/growth/family_W*.json`.
- **EVIDENCE** (bounded by the stated design: canonical universe ≤ 26,
  1/64 cell midpoints, 16–400 points per arc, float layer-1 screen with
  exact spot-confirmation): the coverage table, the staircase of minimal
  lengths, the 112.5°/120°/130°/135° frontier angles, and the death-angle
  measurements. The exact pointsweeps at 112.4/112.6° and 119.7/120.3° are
  exact-arithmetic confirmations at the stated finite apex samples.
- **MAP**: the uncovered-set geometry (items 1–5 above) and the
  machine-readable coverage map (`data/summary/coverage_map.json`).
- **SPECULATION** (labelled above): the closed-form death-angle laws; the
  identification of the census frontier with the published ~112.5° frontier.
- **NOT claimed:** anything about unstable orbits; any positive or negative
  statement about triangles outside the certified boxes; coverage of the
  boundary fringe (0.016 area) or of γ ∈ (90°, 91°) (no interior cells at
  this grid); that the uncovered 64.3% has no stable orbit of length ≤ 26
  anywhere in each cell (only at the 739 midpoints, float-screened, spot-
  confirmed exactly); minimality of anchor word lengths beyond the stated
  caps (top-4 candidates per length, 4 variants, 2 starting exponents).

The kill condition asked whether minimal certificate length grows
super-exponentially in γ. Measured answer: **no** — along the family curve
it grows hyperbolically, ℓ ≈ 1440/(180−γ). The route dies differently: the
*area* certified per word collapses exponentially in word length (tile
pinching + 2-bits-per-letter enclosure loss), and the 30-60-90-type funnels
and nice-angle pinch gaps mean no finite length covers even the 92–135°
band. Census-by-length is the wrong axis: length buys angle reach but not
area.

## Why it failed / what survived

**Obstruction, specifically.** Three separate mechanisms, all measured:
(1) tile-width collapse — the certified/alive windows of the covering words
shrink to curves as γ grows or approaches the arc, so area coverage stalls
at ~36% even with a 17,527-word universe; (2) accumulation funnels — below
the 30-60-90 arc point and at the family's own death angles 180a/(a+1)°,
minimal length is locally unbounded (nothing ≤ 26, and at 135° nothing ≤ 30
known), reproducing blind the published picture of tile patterns stopping
short of certain parameter points; (3) practically, the harness's affine
enclosure loses ≈ 2 bits of certifiable box size per letter, so even where
tiles exist, exact certification at map scale is impossible for ℓ ≥ 14 —
certified area from 376 anchors totals ~2.2·10⁻⁶.

**What survived.** The W(a,b) = (0(12)^a(02)^b)² family: a single
two-parameter generalization of Fagnano that (measured, partly predicted-
then-confirmed) organizes the entire constructive frontier — the shortest
stable words at every length 14–26, death angles on clean rational-multiple-
of-π arcs, certified members to γ = 159.25°, and (SPECULATION) plausibly the
same object as the published 112.3° construction at its (2,1) member. The
census tooling itself (word DFS with incremental translation bookkeeping,
two-layer map/anchor design, death-bisection, exact pointsweep refutation)
is reusable for any future word-family study.

## Leads generated

1. **Prove the death-angle laws.** Conjecture: the corridor of W(a,a)
   degenerates exactly on the arc γ = 180·a/(a+1)°, and of W(a,a−1) at
   180·(1−(2a−1)/(2a²))°. The unfolding is explicit; the corridor endpoints
   are rational functions of the apex, so this is a finite symbolic
   computation per a — falsifiable by `census.py family` at higher
   precision or by exact certification along parametrized apexes. Test
   prediction first: W(10,9) dies at 162.90°, W(12,12) at 166.154°.
2. **Close the pinch gaps.** Measure what covers [death(W(a,a)),
   birth(W(a+1,a))] (e.g. [135.000°, 135.049°]): scan words of length 28–34
   (universe DFS is ~1–4 h in C, or restrict to words with the family's
   letter statistics). If nothing bounded covers a neighbourhood of
   γ=180a/(a+1) on the arc-minimum curve, these angles are genuine
   accumulation points of the constructive problem — a precise, falsifiable
   statement.
3. **The 30-60-90 funnel.** Quantify minimal candidate length vs distance
   to apex (1/4, √3/4) using layer 1 at grid 128 around that point; compare
   against a power law. Decides whether the funnel is log-shallow (finite
   census can close it) or polynomial-deep (it cannot).
4. **Beat the enclosure, not the search.** The harness loses ~2 bits/letter
   because first-order forms square away correlation in the foot-parameter
   products. A second-order (Taylor-model) variant or interval-Newton on
   the two corridor endpoint functions would plausibly certify at cell
   scale for ℓ ≤ 26, turning this census's 36% candidate map into ~36%
   *certified area*. Concrete and testable on the existing anchor data —
   but it is a harness change, so it must be built as a new tier-0 tool and
   validated against the current one.
5. **Near-90° structure.** W(2,1) survives to γ = 90.017° in a window
   pinching to the right-isoceles corner (x→1/2). Identify the limiting
   orbit in the 45-45-90 triangle (it should be an unstable right-triangle
   orbit) and check whether a perturbation family exists on the *other*
   side (x slightly > 1/2, i.e. the mirrored half) — if yes, the 91–92°
   uncovered band is sliver-covered from both sides and the census grid just
   cannot see it; falsifiable at grid 256 restricted to that band.
6. **Letter-statistics pruning.** All last-survivors have the form
   0(12)^a(02)^b doubled — one 0-1 contact, everything else alternating
   with 2. If provable that high-γ tiles require words with ≤ c contacts
   between sides 0 and 1 (the two sides adjacent to the obtuse vertex meet
   at an ever-flatter angle, SPECULATION), the universe at length ℓ drops
   from ~1.85^ℓ to polynomial, and the census becomes feasible to ℓ ~ 60.

## References

- Working copy `problems/billiards-triangles/PROBLEM.md` (tier-0 published
  background: Fagnano 1775; Masur 1986; Schwartz, Experimental Math. 18
  (2009) 137–171; Garber–Marinov–Moore–Tokarsky arXiv:1808.06667; Forni
  arXiv:2606.10102).
- `harness/billiards-triangles/unfold.py`, `verify_cover.py` (tier-0
  reference implementations; all certificates and re-verifications above go
  through them unmodified).
- New code: `problems/billiards-triangles/explore/census.py`. Data:
  `problems/billiards-triangles/data/{words,cover,growth,verify,summary}/`.
- Note: this blind working copy contains no `tests/` directory and is not a
  git repository; the three selftests (`unfold.py --selftest`,
  `verify_cover.py --selftest`, `census.py selftest`) all pass as of this
  record, and CI/pytest must be run at merge time by the integrating agent.
