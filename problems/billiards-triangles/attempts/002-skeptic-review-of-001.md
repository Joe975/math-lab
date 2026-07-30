# 002 — Skeptic review of 001 (word census by length)

- **Problem:** billiards-triangles
- **Date:** 2026-07-30
- **Mode:** blind (this review was done entirely inside the tier-0 blind
  working copy; no other checkout or prior-art file was read)
- **Type:** skeptic review of `ATTEMPT.md` (census attempt, same working copy)
- **Tools:** new from-scratch code under `problems/billiards-triangles/explore/`
  (`skeptic_orbit.py` — exact-rational orbit checker; `skeptic_enum.c` —
  independent word enumeration; `skeptic_family.py` — independent float
  corridor, death bisection, sliver probes; `skeptic_aggregate.py` —
  independent aggregation). The tier-0 harness `unfold.py` was used only as a
  cross-reference on 5 sampled boxes and for the short-length translation-test
  cross-check; `verify_cover.py` and the attempt's `census.py` were **not**
  re-run and none of their code was reused.
- **Sources:** `ATTEMPT.md`, `PROBLEM.md`, `harness/billiards-triangles/*`,
  `problems/billiards-triangles/explore/census.py` (read, not executed),
  `problems/billiards-triangles/data/**` (read-only inputs)
- **Stance:** adversarial; default REFUTE. Every load-bearing number below was
  re-derived by code written for this review, not by re-running the attempt's.

My data lives under `problems/billiards-triangles/data/skeptic/`. No file of
the attempt was modified.

---

## Claims attacked

All six load-bearing claims of attempt `001-word-census-coverage-map.md`
(referred to below as ATTEMPT.md, its filename in the blind working copy),
each attacked with from-scratch re-implementations. Per-claim findings:

### Claim 1 — 376 exact certificates. **CONFIRMED**

**Independent method.** `skeptic_orbit.py` is a from-scratch exact-rational
checker with deliberately different implementation choices from both tier-0
tools and `census.py`: it unfolds by reflecting **all three vertices** across
the gate line (no foot-parameter shortcut), accepts the translation property
only **geometrically** (all three vertex displacements equal after the word),
intersects the corridor from scratch, and then re-derives the orbit with its
own billiard simulator (own Cramer ray/segment solve, own direction
reflection), requiring: struck sequence = word rotation, exact closure of
position *and* direction, every bounce strictly interior to its side.
Calibration (`selfcheck`): Fagnano accepted on an acute apex, rejected on
obtuse, junk word rejected, single traversal rejected. (One sign error in my
own Cramer solve was caught by this selfcheck and fixed — the checker can
fail.)

**Sample (seeded, seed 20260730):** 36 certificates spanning all three
sources — 15 of 260 map anchors, 12 of 107 arc-sweep anchors, all 9 family
anchors including the longest, W(8,8), length 66, half-width 2⁻⁵⁷ at apex
(15381/32768, 2989/32768), γ ≈ 159.25°. Each verified at **5 exact rational
apexes** (centre + 4 corners pulled to 3/4 half-width of the certified box):
180 exact verifications.

```
python3 .../skeptic_orbit.py check-anchors --src data/cover/anchor_g32.jsonl --sample 15 --seed 20260730 --out data/skeptic/check_map_anchors.jsonl
python3 .../skeptic_orbit.py check-anchors --src data/growth/arcsweep_*.jsonl --sample 12 --seed 20260730 --out data/skeptic/check_arc_anchors.jsonl
python3 .../skeptic_orbit.py check-anchors --src data/growth/family_W*.json --out data/skeptic/check_family_anchors.jsonl
```

**Result: 36/36 confirmed, 0 failures.** Cross-reference: the box claim
itself (harness `unfold.certify` = TRUE on the full dyadic box) re-checked on
5 sampled boxes including the length-66 one — all TRUE. Count arithmetic
confirmed from the raw files: 260 (map) + 107 (arcs) + 9 (family) = 376.
The 4 NO_ANCHOR cells are, as stated, at cand_len 22/26.

*Bounded by sampling:* 340 of 376 certificates were not re-checked by me;
they were all attacked by the attempt's own simcheck (a different code path
from the certifier), and my 36-sample drew from every source file with a
published seed.

### Claim 2 — universe = 17,527 canonical words. **CONFIRMED (exact match)**

**Independent method.** `skeptic_enum.c` (cc -O2): full DFS over all bounce
words to length 26 with its own hand-derived direction bookkeeping
(reflection d → 2t−d on integer triples p·α+q·β+r·π), own
rotation+reversal canonicalization, own u^k pruning (drop only when u is
itself a translation word). Runtime 2.2 s.

- Counts by length: **6:1, 10:3, 12:4, 14:15, 16:30, 18:109, 20:318,
  22:1101, 24:3544, 26:12402, total 17,527** — identical to `counts.json`
  and to the record. No words at lengths 4 or 8, matching the record's gap.
- Stronger than counts: the **word lists themselves are identical** —
  `sort` + `diff` of my ten files against `data/words/words_L*.jsonl`:
  no differences at any length.
- My translation test vs the tier-0 harness (`unfold.is_translation_word`):
  all **8,190** bounce words of length ≤ 12, **0 disagreements**.

### Claim 3 — canonicalization soundness. **CONFIRMED**

**Derivation re-done by hand.** (a) *Rotation:* the chain of
rot(w) = s₂…sₙs₁ is congruent (via the reflection taking T₀ to T₁) to the
sub-chain T₁…Tₙ₊₁ of w's doubled chain, whose gates are g₂,…,gₙ, gₙ₊₁ with
gₙ₊₁ = g₁ + τ. Projection onto the normal of τ kills the τ shift, so the
projected gate intervals — hence the corridor — are the same set. The
corridor test is isometry-invariant, so rot(w) certifies at an apex iff w
does. (b) *Reversal:* running the chain backwards from Tₙ = τ(T₀) unfolds
reversed(w) with the same gate set, translation −τ, same normal line; the
chains are congruent by the translation −τ. (c) *W5:* gates(u^k) =
⋃ᵢ τᵤⁱ(gates(u)) and τ_{u^k} = k·τᵤ, so the projected interval set equals
u's — dropping u^k when u is a translation word loses nothing. Fagnano
(u = 012 not a translation word) is correctly kept.

**Empirical exact check** (`skeptic_orbit.py invariance`): 8 length-14 words
× 4 random rational apexes and 5 length-20 words × 3 apexes, comparing the
**exact Fraction corridor width across every rotation and every rotation of
the reversal** (28–40 variants per word): **0 violations** — widths
identical as exact rationals, not merely same-signed. The order-dependence
the record flags is real but only afflicts the affine *enclosure* on
positive-width boxes (TRUE vs UNKNOWN), never the point corridor, and the
attempt's variant-retry handles exactly that.

### Claim 4 — frontier at 112.5°/120°. **CONFIRMED (as scoped)**

**Independent method:** my exact corridor (Fractions, all-vertex unfolding)
at apexes taken from the attempt's pointsweep files; verdicts mine.

- γ ≈ 112.6°: 16 of the 80 exact apexes (seeded sample) × all 15 canonical
  length-14 words → **0 positive corridors** (their file: 0 TRUE over 80).
- γ ≈ 112.4°: all 80 apexes × 15 words → **exactly 1 positive**, the word
  `01212020121202` = W(2,1), at (4899/16384, 9699/32768) — matching their
  1 TRUE.
- γ ≈ 120.3°: 6 apexes × 109 length-18 words → 0 positives; γ ≈ 119.7°:
  10 apexes → 1 positive (W(2,2) word) — consistent with their 7 TRUE over
  60 apexes.
- **Sliver probe** (the record's own caveat): my own float corridor at
  **4000 x-points** per arc (50× the attempt's 80), all 15 words:
  γ = 112.6° → 0 positive, max width −2.2·10⁻⁵; γ = 112.51° → 0 positive
  (max −7.4·10⁻⁶); γ = 112.4989° → 0 positive at this grid; γ = 112.49° →
  3 positive samples (W(2,1)). So the length-14 family's alive window closes
  between 112.49° and 112.51° on my measurements too, and no sliver at
  1/8000 x-resolution survives at 112.6°. A thinner sliver remains logically
  possible; the record says so itself, and the negative claim is correctly
  stated as sample-bounded EVIDENCE.
- Death of the last length-14 word: my bisection gives 112.49887 vs the
  record's 112.4989. Match.

### Claim 5 — family law W(a,b), death angles. **CONFIRMED (incl. new predictions)**

**Independent method:** `skeptic_family.py death` — my own float corridor
(different unfolding path), own y-for-γ solve, own 400-point-per-arc
bisection, 26 iterations.

| word | my death | record | law prediction |
|---|---|---|---|
| W(2,1) len 14 | 112.49887 | 112.4989 | 112.5 |
| W(2,2) len 18 | 119.99910 | 119.9991 | 120 |
| W(3,2) len 22 | 129.98806 | 129.9881 | 130 |
| W(3,3) len 26 | 134.99412 | 134.9941 | 135 |
| **W(5,4) len 38** | **147.59736** | *not measured by the attempt* | **147.6** |
| **W(7,7) len 58** | **157.49764** | *not measured by the attempt* | **157.5** |

The two bold rows are out-of-sample predictions of the SPECULATION laws
D(a,a) = 180a/(a+1)°, D(a,a−1) = 180(1−(2a−1)/(2a²))°, measured here for
the first time: both land within 0.003° of the prediction, always slightly
below — exactly the finite-sampling under-estimate direction the record
predicts. The laws remain unproven (correctly labelled SPECULATION) but they
survived a genuine refutation attempt.

### Claim 6 — coverage table and aggregation. **CONFIRMED (one figure corrected)**

**Independent method:** `skeptic_aggregate.py` recomputes everything from
the raw `map_g32.jsonl`/`anchor_g32.jsonl` with its own bookkeeping.

- Cells: 739 interior + 94 boundary; 264 CAND, 475 NOCAND, 260 anchored ✓.
- Interior area 0.1804, boundary fringe 0.0159 of π/16 ≈ 0.1963 ✓.
- Coverage rows (cand%/anchored%): 14: 7.6/7.2 · 18: 15.7/14.9 ·
  20: 22.7/21.9 · 22: 29.9/29.0 · 24: 30.7/30.0 · 26: 35.7/35.2 — all ✓;
  "no candidate ≤ 26 at midpoint" = 64.3% of interior ✓.
- γ bands: 92–95°: 54.3 · 95–100°: 75.7 · 100–105°: 60.6 ·
  105–112.3°: 59.1 · 112.3–120°: 60.0 · 120–135°: 38.6 · **135–180°: 0.0** —
  all match the record's 54/76/61/59/60/39/0 ✓.
- Staircase from the arcsweep files: ℓ=20 on 90.05–91.5°, 14 on 92–112.45°,
  18 on 112.5–119.5°, 22 on 120–129.5°, 26 on 130–134.5°, nothing ≤ 26 at
  135–150° ✓.
- Spot-check: 8 seeded-random cell midpoints, full 17,527-word universe with
  my own float corridor → my minimal candidate length equals theirs in all
  8 (including 4 NOCAND cells) ✓.
- **Discrepancy:** "certified area from 376 anchors totals ~2.2·10⁻⁶"
  (Why-it-failed section). 2.2456·10⁻⁶ is the **260 map anchors only** (it
  is `exactly_certified_area_anchors` in `coverage_map.json`); the sum over
  all 376 anchors is **≈ 4.3·10⁻⁶** (map 2.25 + arcs 1.93 + family
  0.12, ·10⁻⁶). The record understates its own certified area by ~2×; the
  qualitative point (microscopic certified area) is unchanged.

## Refutations found

No load-bearing claim was refuted. One numeric figure and three wording
scopings need correction, detailed here.

### Wording / scoping of the Outcome section

Mostly exemplary: SCOPE LIMIT on unstable orbits stated up front and repeated
under NOT-claimed; negative claims are sample-bounded; SPECULATION is
labelled inline; VERIFIED is applied to concrete certified boxes, not the
conjecture. Four flags, all minor:

1. **CORRECTED:** the ~2.2·10⁻⁶ certified-area figure (see Claim 6).
2. "**64.3% of the obtuse region** has no stable-orbit word of length ≤ 26
   at its cell midpoint" — strictly 64.3% **of the resolved interior**
   (91.9% of the region); the 0.0159 boundary fringe is excluded. The
   NOT-claimed section says so, but the bold sentence itself should say
   "of the interior cells".
3. The VERIFIED bullet says each anchor was "independently re-derived by
   exact rational simulation" — true **at the anchor centre**; the full
   9-sample box attack ran on every 5th (map) / 7th (arc) record only. The
   bullet should carry the "at the anchor centre" qualifier that step 4 of
   the What-was-done section does carry.
4. Kill-condition verdict "Measured answer: **no** — it grows
   hyperbolically" — measured only to γ ≈ 160° (a ≤ 8) and along the family
   curve; beyond that the hyperbolic form rests on the SPECULATION death
   law. Suggest "no super-exponential growth observed up to γ = 160°;
   hyperbolic form beyond is the (SPECULATION) law's extrapolation".

## Claims that survive

### Overall verdict

**CONFIRMED, with the one numeric correction above.** The headline —
coverage staircase, the 112.5°/120°/130°/135° frontier structure, the
W(a,b) family with its death-angle law, and the exact certificate base
(including the length-66 orbit at γ ≈ 159.25°) — survived an adversarial
pass in which every layer was re-implemented from scratch: 36/36 sampled
certificates verified at 180 exact apexes, the 17,527-word universe
reproduced word-for-word, the canonicalization argument re-derived and
exactly invariance-tested, the frontiers reproduced at 50× sampling
density, and the death law confirmed on two prediction rows the attempt
never measured. No headline claim was refuted.

**What I sampled rather than exhausted:** 36/376 certificates (5 apexes
each); 16/80 apexes at 112.6° and 6+10/60 at 120.3°/119.7° (exact); 8/739
map cells for the full-universe spot-check; family rows (2,1), (2,2),
(3,2), (3,3) re-measured plus new rows (5,4), (7,7). Total compute for this
review: well under 5 minutes.

## Corrections for the index

- `verifies`: ATTEMPT.md's claims 1–6 as scoped above.
- `refutes`: nothing.
- Correction to carry: certified-area total over all 376 anchors is
  ≈ 4.3·10⁻⁶, not ~2.2·10⁻⁶ (which is the map-anchor subtotal); plus the
  three wording tightenings listed above.

## Files created by this review

- `problems/billiards-triangles/explore/skeptic_orbit.py`
- `problems/billiards-triangles/explore/skeptic_enum.c` (+ compiled binary)
- `problems/billiards-triangles/explore/skeptic_family.py`
- `problems/billiards-triangles/explore/skeptic_aggregate.py`
- `problems/billiards-triangles/data/skeptic/check_{map,arc,family}_anchors.jsonl`
- `problems/billiards-triangles/data/skeptic/enum/words_L*.txt`
- `problems/billiards-triangles/data/skeptic/death_W{21,22,32,33,54,77}.json`
- `problems/billiards-triangles/data/skeptic/aggregate_check.txt`

## References

- `problems/billiards-triangles/attempts/001-word-census-coverage-map.md` — the record under review (filed from the blind working copy, where it was `ATTEMPT.md`).
- `harness/billiards-triangles/unfold.py`, `harness/billiards-triangles/verify_cover.py` — tier-0 harness, used only as cross-reference as stated in Tools.
- No external papers were consulted; this review ran entirely inside the tier-0 blind copy.

---

*Filing note (orchestrator): this review was written in the blind working copy as `REVIEW.md`; section organization was adjusted to the repo review shape on filing. Content, verdicts and numbers are unchanged.*
