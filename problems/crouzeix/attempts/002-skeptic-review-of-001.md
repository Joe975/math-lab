# 002 — Skeptic review of 001 (dim-3 local-maxima census): adversarial verification

- **Problem:** Crouzeix's conjecture, `problems/crouzeix/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed
- **Type:** adversarial verification of `001-blind-dim3-local-maxima-census.md`
  (default stance: refute). Certificates re-verified with an independent
  implementation sharing no algorithm with the tier-0 harness; stall verdicts
  re-probed with a different scheme, RNG, and consistent resolution; citations
  checked against the papers.
- **Outcome in one line:** 001's EVIDENCE stands — 30/30 stratified
  re-certifications consistent, both anchors reproduced to the digit, all 4
  stall verdicts confirmed, both attacked ratio-1 maxima survive — but the
  record miscounts its near-2 class (94, not 97), its index entry overstates
  "all 286 endpoints certified below 2" (the certificates prove *0 certified
  above 2*; 182/286 enclosure upper ends exceed 2), and the Calcolo 2021 paper
  it leans on is by Overton alone, not Greenbaum–Overton.
- **Tools:** `explore/skeptic_recert.py` (independent certifier: Sylvester
  minors by cofactor expansion + exact Rayleigh lower bounds + integer-isqrt
  root bounds + rounded-trig rational directions + Bernstein-coefficient
  segment bounds; no hull-ordering code; deterministic; ~2.8 s/endpoint),
  `explore/skeptic_probe.py` (sparse-perturbation escape probe, RNG 20260804,
  512-vs-512-sample confirmation, second float evaluator via complex Jacobi
  rotations); checkpoints `data/skeptic_recert.jsonl`, `data/skeptic_probe.jsonl`,
  logs alongside. Reproduce with:
  `python3 explore/skeptic_recert.py --selftest | --anchors | --sample` and
  `python3 explore/skeptic_probe.py stalls 300010 604002 406004 504010`,
  `... maxima 403004 603009`, `... stallcheck 400000 403000 405000`.
- **Sources:** 001 and its data/tools; `harness/crouzeix/ratio.py` /
  `verify_ratio.py` (read to audit soundness direction and understand data
  formats — the re-certifier transcribes nothing from them);
  arXiv:2105.14176 (Overton, Calcolo 2021) via the ar5iv machine rendering
  [T]; Greenbaum–Overton preprint optimization-online 5703 (2016), PDF
  fetched and text-extracted.

## Claims attacked

1. **Certificate validity** — the 286 den-256 rational endpoint enclosures
   (tol 1e-6, 32 half-planes), the two calibration anchors, the 1.969725
   record lower bound, and the 1.974 certification ceiling at the exact 2×2
   Jordan pair.
2. **Direction of the W(A) discretization** — does the 32-half-plane outer
   polygon inflate or deflate the denominator, and does the certified
   inequality survive it?
3. **Classification arithmetic** — the 13/97/175/4 accounting against
   `classified.jsonl` (13+97+175+4 = 289 ≠ 286), the 61% stall rate, and the
   per-family start counts against the raw census files.
4. **The stall verdicts** — the 4 deep-probed intermediate candidates
   (1.96764, 1.96026, 1.95884, 1.78117): does ascent really exist? The deep
   probe compared 128-sample escape values to a 512-sample baseline — a
   mismatch whose sign can manufacture apparent ascent.
5. **The ratio-1 maxima** — probe two of the 13 hard (ascent found there
   would refute the census's only "genuine local maximum" class), and check
   the claimed structure on all 13.
6. **Citations and transcriptions** — the ice-cream-cone attribution, the
   intermediate values ~1.433/~1.185 ([L]-transcribed in 001), the "~half a
   million heavy-tailed starts", and the claimed design differences from
   Greenbaum–Overton.
7. **Mode and scope** — the blind label, tier-1 contamination, and an
   overclaim audit of the Outcome section and the prior-art index entry.

## Refutations found

### R1. The near-2 class has 94 members, not 97

`classified.jsonl`: near2-double 48 + near2-other 44 + near2-triple 2 = **94**
(and 175 + 94 + 13 + 4 = 286, as it must). 001's own table says 48/44/2, but
its prose bullet says "The 97 near-2 endpoints", and the index `one_line`
repeats "97 Jordan-type endpoints near ratio 2". Corrected statement: **94
near-2 endpoints** (48 double-coalescence, 44 separated-eigenvalue, 2
triple-coalescence). The impossible sum 13+97+175+4 = 289 is exactly the
trace of this error.

### R2. "All 286 endpoints certified below 2" (index entry) is not what was proved

The certificates establish (i) a certified lower bound below 2 at every
endpoint and (ii) failure of the exact refutation test `n2_lo > 4·outer2` —
i.e. **no endpoint is certified above 2**. They do *not* certify any endpoint
*below* 2: in `certified.jsonl`, **182 of 286** enclosure upper ends
(`ratio_hi_f`) are ≥ 2 (up to 2.38), as expected for near-extremal points
enclosed with a 32-gon outer polygon. The attempt record's own Outcome text
states it correctly ("0 refutations; the largest certified lower bound … is
1.969725"); the `one_line` in `prior-art.json` overstates it. Corrected
statement: *all 286 endpoints certified — 0 refutations (none certified
above 2), max certified lower bound 1.969725*.

### R3. Calcolo 2021 / arXiv:2105.14176 is by Overton alone

001's References list "A. Greenbaum, M. L. Overton, *Local minimizers of the
Crouzeix ratio: a nonsmooth optimization case study*, Calcolo 58 (2021)", and
the prose attributes that paper's findings to "Greenbaum-Overton" throughout.
The paper is single-authored: **M. L. Overton** (verified against
arXiv:2105.14176). The 2018 LAA paper is correctly Greenbaum–Overton. Every
substantive claim 001 takes from the 2021 paper checks out (S6); only the
attribution is wrong.

### R4. Minor: the deep-probe ascent figures, and a resolution mismatch in the probe protocol

`deep_probe.log` shows ascents +0.00459, +0.00535, +0.00450, +0.00536 — i.e.
**0.0045–0.0054**, not "0.005-0.006" as the record says. More substantively,
`deep_probe.py` compares escape candidates evaluated at **128** boundary
samples against a baseline at **512**; coarser boundary sampling
under-estimates max_{∂W}|p| and so over-estimates the ratio, which can
manufacture ascent. The verdicts survive anyway (S4: consistent-resolution
re-probing finds genuine ascent at all four), but the protocol flaw is real,
and the same mismatch exists in `classify.py` (probe at 96 vs baseline at
512). Quantified: only 20/175 stalls show any self-inflation
(`objective(x, 96) > ratio512 + 5e-4`), and in all 20 the probe found ascent
beyond the point's own 96-sample value — so the 61% stall figure is robust;
future probe protocols should still compare like with like.

### R5. Minor: "single |p|-peak" is not clean on 2 of the 13 ratio-1 endpoints

Re-running the peak diagnostic at 512 samples: 11/13 endpoints show one peak;
seeds 500009 and 100003 show ~10 near-coincident peaks at the 0.999·max
threshold (|p| is near-flat along a boundary arc). In all 13 the *dominant*
peak sits at the boundary eigenvalue, so the configuration identification
stands; the sentence "the single |p|-peak at it" is an overstatement for
those two.

### R6. Framing: the census design was assigned, not blind-convergent

The blind label is accurate for what the worker read (tier-0 copy plus GO
papers via web). But the design brief — this exact slice (n=3, deg ≤ 3),
cost-measurement-before-sizing, and "state how the design differs from GO
before running it" — is the lab's own tier-1 queue item (STATUS.md queue
item 10, derived from `PRIOR-ART.md`'s editorial view), delivered through
the task assignment. No crouzeix tier-1 *findings* existed before 001, so
the results are uncontaminated; but the *choice of approach* must not be
counted as independent blind evidence about what an unprompted agent would
try. Nothing in the record's text presupposes non-public knowledge; this is
a caveat on how the mode field should be read, not a mislabel.

## Claims that survive (and what was thrown at them)

### S1. The certificates — independent re-implementation, 30/30 consistent

`skeptic_recert.py` was written from the PROBLEM.md verification contract
with disjoint algorithms end to end: eigenvalue upper bounds by bisection on
**Sylvester's criterion** (leading principal minors by cofactor expansion —
not LDL, not Sturm/Faddeev–LeVerrier), lower bounds by **exact Rayleigh
quotients** at rationalized power-iteration vectors, square-root bounds by
**integer `math.isqrt` scaling** (no bisection), W(A) outer polygons from
**40 rounded-trig rational directions** with support offsets rounded *up*
(vs the harness's 32 Gaussian-integer coprime pairs), and max|p|² bounded by
**Bernstein coefficients over all pairwise vertex segments** — no
convex-hull-ordering code exists in it, the one place a subtle geometry bug
could hide. Its selftest brackets the literature values (2 for (J2,z), √2
for (J3,z), 1 for a normal matrix).

Stratified sample of 30/286 (all 4 below2-max, 5 of 13 one-exact, 13 near-2
including the record holder, 8 stalls; `data/skeptic_recert.jsonl`): **all
30 consistent** — 001's certified lower bound never exceeds my certified
upper bound, my enclosure always overlaps or sits inside their slack, the
refutation test agrees (false) at every point, and every one of my lower
bounds is < 2. Highlights: record holder seed 406009: their lo 1.969725 vs
my enclosure [1.988100, 1.995295] — their bound is valid and conservative;
the four killed intermediates re-certify at [1.9395, 1.9519],
[1.9411, 1.9508], [1.9419, 1.9567], [1.7479, 1.7751]; the five sampled
ratio-1 endpoints re-certify at [1.000000, 1.000000] to 12 digits. Anchors:
001's `(J3, z²) → [1.980557335, 2.25]` and `(J2⊕[1/4], z) → [1.990266652,
2.0]` reproduce **to the last printed digit** (their pipeline, dir_bound 5,
tol 1e-9), and my independent certifier at 80 directions gives
[1.996454, 2.0] and [1.998226, 2.0] — both contain the true value 2, both
consistent with theirs. The certification ceiling reproduces too: their
pipeline at tol 1e-6 / 32 half-planes gives lo = 1.97417 at the exact
(J2, z) pair, so 1.969725 at rounded census endpoints is indeed "as close to
2 as the 32-gon slack allows".

### S2. The discretization direction is sound

Code audit of `ratio.py`: the outer polygon uses the *upper* ends of the
per-direction eigenvalue enclosures, the vertex-filter slack only *enlarges*
the polygon, and the interval-Horner edge bound is an upper bound — so the
denominator is **over**-estimated, hence `ratio_lo` and the refutation test
are conservative: the safe direction for "no refutation" claims (an
under-estimated denominator would have been fatal, per the contract, and
was the thing this review most expected to find). Empirically, on the
semi-analytic control A = diag(1, i, −1), p = z (true max_{W}|p| = 1, at the
eigenvalue hull): their `outer_max2` = 1.00000024 ≥ 1, `inner_max2` = 1.0
≤ 1, enclosure [0.99999940, 1.0] brackets the true ratio 1, no refutation.
The upper enclosure end correctly uses certified-inner points only. No
direction error found.

### S3. The counts — reconciled exactly

Raw census files: 79+83+79 (interleaved i4/i5/i6) + 19+13+13 (truncated
w1/w2/w3, all family "random", confirming the "family-biased first launch"
note) = 286. Per-family: random 81, jordan3 36, jordan2 34, nearnormal 34,
companion 34, peaks 34, real 33 — matches the record and the index `range`.
Class partition: 175 stall / 48+44+2 near-2 / 13 one-exact / 4 below2-max
= 286 (a true partition; the only error is R1's prose "97"). 175/286 =
61.2% ✓. `certified.jsonl`: 286 rows, 0 degenerate, 0 refutes, max
`ratio_lo` = 1.969724783 at seed 406009 ✓; max stall lower bound 1.9389
("up to 1.94" ✓). `certify.log`: 286 certified in 214 s ✓; top-8
independent verification all "AGREE", zero MISMATCH lines ✓.

### S4. The four stall verdicts — confirmed under a fair protocol

`skeptic_probe.py`: different scheme (900 sparse coordinate-subset Gaussian
perturbations, log-uniform scales 10^-3.5 to 10^-0.7, 3 jittered NM
restarts, fine compass), different RNG (20260804), and **every candidate
confirmed at 512-vs-512 samples** before counting, then cross-checked with
a second float evaluator built on complex Jacobi rotations (agrees with the
census evaluator to 1e-15 on controls). All four killed intermediates
escape genuinely: +0.0040 (seed 300010), +0.0010 (604002), +0.0100
(406004), +0.0219 (504010) — the two largest ascents *exceed* what 001's
deep probe found, so the stall verdicts are if anything understated. Three
randomly spot-checked ordinary stalls also escape at consistent resolution
(+0.0526, +0.2337, +0.0059). No hidden intermediate local maximum surfaced.

### S5. The ratio-1 maxima — survive attack, structure verified

Seeds 403004 and 603009 were probed with the full S4 arsenal: **zero
ascent** (best512 − v512 = +0.000000 for both). All 13 one-exact endpoints
re-checked structurally at 512 samples: in every case one eigenvalue sits
on ∂W(A) (margin −0.0 relative to diam W) and the dominant |p|-peak is
exactly at it, matching the published characterization of the ratio-1 local
minimizers of the reciprocal — A ≈ U diag(λ, B)U* with a dominant outside
scalar block (GO 2018, Theorem 2 and the surrounding discussion), for which
"W(A) has a vertex at λ, and often has the appearance of an ice cream cone"
(Overton 2021). The rediscovery-with-citation framing is correct (mod R3's
authorship fix), and the den-256 certificates pin all five sampled ratio-1
endpoints to [1.000000, 1.000000].

### S6. The transcribed GO facts — all check out against the papers

From the ar5iv rendering of arXiv:2105.14176 [T] and the text-extracted
2016 GO preprint: the intermediate locally minimal reciprocals include
**0.698 at (n=3, m=3)** → ratio 1.4327 ≈ "~1.433" ✓ and **0.844 at (n=2,
m=3) and (n=3, m=3)** → 1.1848 ≈ "~1.185" ✓; initialization is heavy-tailed
(entries x·e^{αx²}, α = 2, x normal), 10,000 starts per configuration,
**nearly 500,000 runs** in total → "~half a million heavy-tailed random
starts" ✓; GO 2018 uses **Chebfun** for bd W(A) plus **BFGS and Gradient
Sampling** from normally distributed random starts ✓. So 001's four claimed
design differences (fixed slice; structured start families; derivative-free
NM + compass with escape probing; certified endpoints) are all real
differences from the published methodology.

### S7. Scope and labels — no overclaim in the record body

The Outcome opens with `EVIDENCE` scoped by dimension, degree, and search
design, and explicitly disclaims any statement about the conjecture; the
"every surviving local maximum sits at already-known structure" verdict is
scoped to this design in the same paragraph and honestly hedged for the
near2-other subclass ("consistent with … coalescence only in the limit" —
44 of the 94 near-2 points have *not* collapsed their eigenvalues, and the
record says so rather than claiming they are Jordan pairs). The SPECULATION
label sits inline on the {1,2}-only basin claim, as required. The declared
gap (GO's intermediate maxima never reached; their basins untested) is
accurate and is the honest headline limitation. The index `gaps` entry
matches. Only the index `one_line` fails this audit (R1, R2).

## Verdict

| # | Claim | Verdict |
|---|-------|---------|
| a | 286 certified enclosures; anchors; 1.969725 record; 1.974 ceiling | **CONFIRMED** (30/30 stratified sample by independent implementation; anchors to the digit) |
| b | Discretization direction sound for "no refutation" | **CONFIRMED** (code audit + semi-analytic normal control) |
| c | Counts 13/97/175/4 and per-family | **CORRECTED** — near-2 is 94 (48/44/2); everything else reconciles exactly |
| d | 4 intermediate candidates are stalls | **CONFIRMED** (fair-protocol re-probe, ascents up to +0.022), protocol flaw noted and quantified harmless |
| e | Ratio-1 ice-cream-cone rediscovery; GO transcriptions | **CONFIRMED**, except **CORRECTED** authorship (Calcolo 2021 is Overton solo); "single peak" soft on 2/13 |
| f | Blind mode; design independence from GO | **CONFIRMED** as labelled, with the design-was-assigned caveat (R6) |
| g | Scope honesty; index fidelity | Record body **CONFIRMED**; index `one_line` **CORRECTED** (R1, R2) |

**Net assessment: 001's EVIDENCE stands.** The certificates are valid and
conservative, the landscape accounting is real, and the two headline
structures (ratio-1 cones, near-2 Jordan-type climbs) are what the record
says they are. The corrections are bookkeeping and attribution, not
substance — but two of them (94 vs 97, "certified below 2") live in the
index line future agents read first, which is exactly where errors compound.

## Residual risk

- The re-certification covered 30 endpoints + 2 anchors of 286; the
  remaining 256 rest on the stratified sample plus 001's own top-8
  cross-route check. The sample covers every class and every family.
- Both float evaluators share the boundary-point construction *concept*
  (support points from Hermitian eigenvectors); an error in that shared
  concept would evade both. The exact-rational certificates bound this
  risk: they bypass boundary sampling entirely.
- The Overton/GO paper details were read through machine renderings and
  text extraction ([T]); the specific numbers re-checked here (0.698,
  0.844, run counts, methods) are as reliable as those renderings.
- 001's Lead 1 (attack the published intermediate maxima directly) is
  untouched by this review; nothing here tests those basins either.

## References

- `problems/crouzeix/attempts/001-blind-dim3-local-maxima-census.md` (the
  record under review) and its data in `problems/crouzeix/data/`.
- This review's tools and checkpoints: `explore/skeptic_recert.py`,
  `explore/skeptic_probe.py`, `data/skeptic_recert.jsonl`,
  `data/skeptic_probe.jsonl`, `data/skeptic_recert.log`,
  `data/skeptic_probe.log`.
- M. L. Overton, *Local minimizers of the Crouzeix ratio: a nonsmooth
  optimization case study*, Calcolo 58 (2021); arXiv:2105.14176 — read via
  the ar5iv HTML rendering [T]; abstract page fetched directly.
- A. Greenbaum, M. L. Overton, *Numerical investigation of Crouzeix's
  conjecture*, Linear Algebra Appl. 542 (2018); preprint
  optimization-online 5703 (2016), PDF fetched and text-extracted for the
  methods check.
- Tier-0 harness: `harness/crouzeix/ratio.py`,
  `harness/crouzeix/verify_ratio.py` (audited, not reused, by this review's
  certifier).
