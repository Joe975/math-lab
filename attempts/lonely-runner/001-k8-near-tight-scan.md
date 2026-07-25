# 001 — k=8 exact gap tool + near-tight speed-tuple scan to V=72

**Problem:** Lonely Runner Conjecture, k = 8 case (7 moving runners).
**Date:** 2026-07-25.
**Tools:** `tools/lonely_runner.py` (exact ML + scan; re-runnable),
`tools/lonely_runner_analyze.py` (independent re-verification of scan output).
**Data:** `attempts/lonely-runner/data/k8_scan_V72.json`,
`attempts/lonely-runner/data/k8_spectrum_V40_thr1-7.json`,
`attempts/lonely-runner/data/k8_scan_V72.log`.

## Approach

**Reduction (standard).** Loneliness of a fixed runner depends only on speed
differences, so subtract that runner's speed: one runner stationary at 0,
the others at pairwise distinct nonzero speeds. A limit/approximation
argument (Bohl; Bienia et al.) reduces to integer speeds, and ||−vt|| = ||vt||
lets us take them positive. So for k = 8: distinct positive integers
v₁ < ... < v₇, and the conjecture asserts

    ML(v) := sup_t min_i ||v_i t|| ≥ 1/8,     ||x|| = distance to nearest integer.

ML is invariant under scaling v → cv (substitute t → t/c), so primitive
(gcd = 1) tuples suffice.

**Exact algorithm.** f(t) = min_i ||v_i t|| is continuous, 1-periodic,
piecewise linear; every linear piece of every ||v_i t|| has the form
v_i t − a or b − v_i t (slopes ±v_i, all nonzero). At a global maximizer t*
the left derivative of f is ≥ 0 and the right is ≤ 0, so some active rising
piece v_i t − a (achieving f just left of t*) meets some active falling piece
b − v_j t (just right of t*) at t*, the case i = j (peak kink of a single
||v_i t||) included. Solving v_i t* − a = b − v_j t* gives t* = (a+b)/(v_i+v_j).
Hence

    ML(v) = max { f(m/(v_i+v_j)) : 1 ≤ i ≤ j ≤ 7, 0 ≤ m < v_i+v_j },

a finite computation; at t = m/d, ||v t|| = min(vm mod d, d − vm mod d)/d,
so everything is exact integer arithmetic (output as `Fraction`).

**Plan.** Validate the implementation on the known tight cases, then scan all
7-subsets of {1..V} for the largest feasible V, recording every tuple with
ML < 1/8 + 1/200 = 13/100 ("near-tight") exactly, and mine the arithmetic
structure of the hits. Speed came from a sound early-bail: any single value
f(t) is a certified lower bound on ML, so a tuple is discarded the moment
any probe or candidate time witnesses f(t) ≥ 13/100 (a handful of small-
denominator probe times kill almost every tuple in O(1); full exact scan
only runs for survivors). Non-primitive tuples were scanned too (cheaper
than per-tuple gcd) and deduplicated in post — which doubles as a
consistency check, since every scaled copy must reproduce the primitive's ML.

## What was done

1. **Validation suite** (`--validate`, seeded, all PASSED):
   - (1,...,k−1) gives ML exactly 1/k for k = 3..8; scaled copy 3·(1..7)
     gives exactly 1/8.
   - 150 random tuples: the sums-only candidate set agrees exactly with a
     strictly larger "paranoid" candidate set (sums, |differences|, 2v_i).
   - 25 random 7-tuples: exact ML vs 10⁵-point float grid agree within the
     Lipschitz slack (max |diff| 1.95e−4).
   - 400 random 7-tuples: the early-bail scanner agrees with full exact ML
     on both the near-tight decision and the exact value.

2. **Main scan** — all C(72,7) = 1,473,109,704 7-subsets of {1..72},
   threshold 13/100, exact arithmetic, 2027 s on 4 cores (~727k tuples/s).

3. **Auxiliary spectrum scan** — all 18,643,560 7-subsets of {1..40} with
   threshold 1/7, to map the ML spectrum in the band (1/8, 1/7) just above
   the conjectured bound (15 s); plus a sharper **band test**: all
   73,629,072 7-subsets of {1..48} with threshold 3/23 (102 s), see
   `data/k8_band_V48_thr3-23.json`.

4. **Independent re-verification** (`lonely_runner_analyze.py`): every
   primitive hit recomputed three ways (sums-only exact, paranoid exact,
   2·10⁵/4·10⁵-point float grid); scaled-copy completeness and ML-equality
   across scales checked. All OK; zero mismatches.

5. **Literature cross-check** (web): the two nontrivial tight tuples found
   are exactly the known Goddyn–Wong tight instances; the off-pattern
   spectrum values match the Fan–Sun amended spectrum conjecture (see below).

## Outcome

**Near-tight census (ML < 13/100), V = 72: exactly 3 primitive tuples, all
with ML = 1/8 exactly — no counterexamples, and an empty open band.**

| primitive tuple | ML | maximizing times | structure |
|---|---|---|---|
| (1,2,3,4,5,6,7) | 1/8 | 1/8, 3/8, 5/8, 7/8 | canonical tight instance |
| (1,2,3,4,5,7,12) | 1/8 | 1/8, 3/8, 5/8, 7/8 | Goddyn–Wong acceleration: 6 → 2·6 in (1..7) |
| (1,4,5,6,7,11,13) | 1/8 | 1/8, 3/8, 5/8, 7/8 | Goddyn–Wong sporadic tight instance |

- The 21 raw hits are precisely the scaled copies (g·v with g·max(v) ≤ 72:
  10 + 6 + 5), each independently computed with identical ML — internal
  consistency check passed.
- **No tuple with ML < 1/8** among all 1.47·10⁹ tuples: consistent with the
  conjecture (and with Rosenfeld's 2025 preprint proof of k = 8, see Leads).
- **No tuple with 1/8 < ML < 13/100**: the spectrum has an empirical gap
  above 1/8, for all speeds ≤ 72. The sharper band test extends this: for
  all speeds ≤ 48, **no ML lies in the open interval (1/8, 3/23)** at all.
- The GW acceleration criterion (replace r by mr in (1..7); tight iff
  gcd(r,x) > 1 for all x in [8−r, m(8−r)−1]) admits **only** r=6, m=2 with
  max speed ≤ 72, matching the scan exactly. So up to speed 72, tight = the
  three known instances; the scan proves this census is complete in that range.

**Structure of the band (1/8, 1/7), V ≤ 40 (24 further primitives):**
ML values observed, with multiplicities:

    3/23 (×2), 2/15 (×11), 5/37 (×1), 3/22 (×2), 4/29 (×5), 5/36 (×2), 6/43 (×1)

- Every value has the form **s/(7s+k) with k ∈ {1,2}**: s/(7s+1) for
  s = 1..6, plus s/(7s+2) for odd s = 3, 5. (For even s, s/(7s+2) reduces
  into the k=1 family, so odd s is the only way k=2 produces new values.)
  Kravitz's original spectrum conjecture allows only k = 1; the value
  3/23 < 2/15 at (1,2,3,4,5,7,18) is a Kravitz-violating value of exactly
  the type Fan–Sun found for n = 4, and everything observed fits their
  amended conjecture ML ∈ {s/(ns+k), k ≤ n} — with only k ≤ 2 appearing.
- **Acceleration families dominate.** Accelerating one runner of a tight
  instance explains almost all near-tight tuples:
  (1,2,3,4,5,6,7m) has ML = m/(7m+1) for m = 2,3,4,5 (exact ladder);
  (1,2,3,4,5,7,6m): m=2 tight, m=3 → 3/23, m=4 → 4/29;
  (1,3,4,5,6m,7,11): m=1..4 → 2/15, m=5 → 5/37, m=6 → 6/43;
  (1,4,5,6,7,11,x) for x ∈ {10,16,26,39} → 2/15, x=23 → 4/29, x=29 → 5/36.
- All 27 primitives below 1/7 contain a pair v_i + v_j ≡ 0 (mod 8) and no
  speed ≡ 0 (mod 8); all three tight ones attain their max exactly on the
  odd multiples of 1/8 and contain both residues ±1 mod 8 (speeds 1 and 7).
- Two sporadic near-tight primitives contain no speed 1:
  (2,6,7,8,10,13,14) → 2/15 and (2,7,9,11,12,13,20) → 4/29.

## Why it failed / what survived

Nothing failed computationally; the honest calibration is about what this
can and cannot say:

- This is **evidence, not proof**: bounded-speed exhaustion says nothing
  about speeds > 72 (Tao-type bounds needed to check k = 8 are astronomically
  larger), and mid-scan a literature check surfaced Rosenfeld's Sept 2025
  preprint (arXiv:2509.14111, revised Oct 2025, not yet peer-reviewed)
  claiming a computer-assisted proof of the full k = 8 case — so the
  counterexample-hunting aspect of this line is likely moot.
- What survives regardless of that proof: **tightness rigidity and spectrum
  structure remain open**, and this scan gives an exact, adversarially
  re-checkable census: up to speed 72, tight instances are exactly the three
  known ones, and the spectrum in (1/8, 13/100) is empty.
- The problem file's claim that tight cases are conjecturally "speeds
  {0,1,...,k−1} up to scaling" is **wrong as stated** — Goddyn–Wong tight
  instances (both recovered here from scratch) are the known correction.
  `problems/lonely-runner.md` needs a status update (Rosenfeld + GW).

## Leads generated

1. **Spectrum-gap conjecture at n = 7, concrete and falsifiable:** the data
   suggests spectrum ∩ (1/8, 1/7) = {s/(7s+1) : s ≥ 2} ∪ {s/(7s+2) : s ≥ 3 odd},
   with an empty band (1/8, 3/23). The band test **confirmed the sharper gap
   directly: for all speeds ≤ 48, no ML lies in (1/8, 3/23)** — the only
   values below 3/23 are the three tight instances at exactly 1/8, so 3/23
   (attained by (1,2,3,4,5,7,18)) is the true second spectrum point in this
   range. Extending the band test to V = 72 costs ~35 min; proving the gap
   (1/8, 3/23) for k = 8 even restricted to bounded speeds would be a
   publishable-style partial result.
2. **Acceleration calculus:** the exact ladders ML((1,2,3,4,5,6,7m)) =
   m/(7m+1) and the GW-criterion match suggest a clean general formula for
   ML of "one accelerated runner" tuples; worth proving directly (finite
   check per m via the candidate-set method, or a covering argument mod 7m+1).
   This would turn most of the near-tight census into theorems.
3. **Rigidity beyond V = 72:** GW's criterion bounds where new acceleration-
   type tight instances can appear; combining it with this scan's sporadic-
   free census up to 72 could support a "no new tight instances below V"
   statement with a much larger V, since the scan bail gets *faster* per
   tuple as V grows for non-near-tight tuples.
4. **Housekeeping:** update `problems/lonely-runner.md` — status (Rosenfeld
   preprint for k = 8; conjecture open for k ≥ 9), and the tight-instance
   classification (Goddyn–Wong). The tooling here works for any k: the k = 9
   near-tight scan (8 moving runners, conjectured bound 1/9) is the natural
   next queue item and is genuinely open territory.

**References** (checked via web search this cycle):
[Goddyn–Wong tight instances / Fan–Sun amended spectrum](https://arxiv.org/html/2306.10417),
[Kravitz, Barely lonely runners](https://arxiv.org/pdf/1912.06034),
[Rosenfeld, The lonely runner conjecture holds for eight runners](https://arxiv.org/pdf/2509.14111),
[Rusza–et-al. survey: The Lonely Runner Conjecture turns 60](https://www.sciencedirect.com/science/article/pii/S1574013725000747),
[Barajas–Serra, The lonely runner with seven runners](https://arxiv.org/pdf/0710.4495).
