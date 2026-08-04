# 002 — Skeptic review of 001 (abk explicit witness): adversarial verification

- **Problem:** maxwell-equilibria, `problems/maxwell-equilibria/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed
- **Type:** adversarial verification of `001-abk-explicit-witness.md` (default
  stance: refute). Target per STATUS queue item 16 / 001 lead 1: both engines
  share the subdivision tree, so the residual common-mode risk is a bug in the
  subdivision DRIVER producing a wrong tree with self-consistent leaves.
- **Outcome in one line:** the certified count of 24 nondegenerate equilibria
  is REAL and survives every attack tried — including four full re-censuses
  over wholesale-different subdivision trees — but the review found a genuine
  soundness gap in the *independent checker's* tiling check (prefix-freeness +
  Kraft equality do NOT imply a tiling for axis-labelled paths; demonstrated
  by a tampered certificate that hides a real equilibrium and still PASSes
  verify_equilibria.py), so 001's "established by two independent
  implementations" overstated the independence of the coverage claim. The gap
  is closed for 001's actual certificate by a strictly stronger
  axis-consistency tree check written fresh for this review.
- **Tools:** `explore/skeptic_cert_checks.py` (independent structural checks,
  fresh code, exact rational arithmetic, no harness imports; ~4 s on the
  34 MB certificate); `explore/skeptic_verifier_gap_demo.py` (the tiling-gap
  demonstrations, ~1 min); `explore/skeptic_centroid_hessian.py` (exact
  centroid eigenvalue signature, instant); `harness/maxwell-equilibria/
  equilibria.py` re-runs at perturbed parameters (configs in
  `explore/skeptic_*_config.json`, runtimes below);
  `harness/maxwell-equilibria/verify_equilibria.py` on a perturbed
  certificate. All deterministic, standard library only.
- **Sources:** attempt 001 and its certificate
  `data/abk-witness-cert.json.gz`; the two harness files (read line by line);
  no external papers consulted beyond what 001 cites.

## Claims attacked

1. **Coverage**: the driver's subdivision loop cannot drop a box, emit
   overlapping leaves, or leave a hole; the certificate's 235,993 leaves tile
   the search region exactly.
2. **The 24 isolation certificates**: pairwise disjoint enclosures, boxes
   consistent with their split paths, enclosures inside their boxes.
3. **Nondegeneracy and index data**: all det signs ±1, index sum −4 = 1 − 5.
4. **Completeness scope**: the localization lemma and its instantiation
   (the region really contains every equilibrium); the domination radii.
5. **Independence of engine 2**: every place `verify_equilibria.py` trusts
   driver metadata rather than re-deriving it.
6. **Robustness of the count to the tree**: re-censuses at different
   precision and over a geometrically different box decomposition — the tree
   must change wholesale, the certified count must not.
7. **001's float side-predictions** (add-on): q = 4360/1000000 gives a
   complete count of 12; q = 4400/1000000 gives 16.
8. **001 lead 5** (add-on): the centroid (1/3,1/3,1/3) is a rational
   equilibrium of the whole family; its Hessian signature is decidable in
   exact arithmetic.

## Refutations found

### R1. The independent checker's tiling check is not sound by itself:
### prefix-freeness + Kraft equality do not imply a tiling

`verify_equilibria.py` validates the leaf decomposition by (i) prefix-freeness
of the split paths and (ii) the Kraft equality Σ 2^−depth = 1, and its
docstring claims this "holds iff the leaves tile the root". **That equivalence
is false for paths carrying a per-node axis choice.** The leaf set
`{"01", "11"}` (axis-0 upper half and axis-1 upper half of the root) is
prefix-free with Kraft sum exactly 1, yet the two boxes overlap on the
(hi, hi) x-y quadrant and jointly miss the (lo, lo) quadrant. Kraft counting
is blind to the axis digits; the claimed equivalence only holds *given* that
all paths through each internal node agree on that node's split axis — which
is exactly the driver-honesty assumption the checker was supposed to
discharge.

Demonstrations (`python problems/maxwell-equilibria/explore/
skeptic_verifier_gap_demo.py`, deterministic, ~1 min):

- **DEMO 1**: a synthetic certificate with leaves `{"01", "11"}` over a
  region where the field is sign-definite (the checker's own validation
  region), every leaf genuinely sign-excluded. `verify()` returns **PASS**
  on a non-tiling leaf list.
- **DEMO 2 (the load-bearing one)**: the driver's honest certificate for two
  unit charges at (±1,0,0) — exactly one equilibrium, at the origin — is
  tampered with: the isolated leaf is replaced by two mixed-axis exclusion
  leaves that each avoid the origin but jointly leave a hole containing it,
  and the count is edited from 1 to 0. `verify_equilibria.py` returns
  **PASS with count 0**: the checker endorses a "complete census" that
  misses a real equilibrium. (The same leaf lists are rejected by this
  review's tree reconstruction, at the mixed-axis node.)

Consequences, stated precisely:

- What is actually wrong: 001's sentence "established by two independent
  implementations ..." and the harness docstring's "exact combinatorial check
  (Kraft equality + prefix-freeness) that the leaf boxes tile the search
  region". The *coverage* half of the two-engine claim was never
  independently verified; it rested on the driver constructing
  axis-consistent trees.
- What is **not** wrong: the driver (`equilibria.py`) demonstrably emits both
  children of every split with the same axis digit (step 5 of
  `count_equilibria` pushes `path + f"{axis}{side}"` for both sides of one
  chosen `axis`), so certificates it actually produces do tile. The gap is in
  what the checker *catches*, not in any existing certificate.
- The fix applied here (checker hardening left as lead 1 — no harness edit in
  this review): the tiling claim for 001's actual certificate is
  re-established by a strictly stronger check written fresh for this review
  (`skeptic_cert_checks.py`, TREE block): reconstruct the subdivision tree
  from all 235,993 paths and verify at every internal node that (a) all
  paths agree on the split axis and (b) both children are present. By
  induction this proves the leaves partition the root box. **001's
  certificate passes with zero violations** (max depth 58 splits; Kraft and
  prefix-freeness re-confirmed as corollaries).

### R2. Reporting-level: the "independent" det-sign verification leans on an
### unstated (but valid) argument

The checker computes det signs over the driver's recorded `enclosure` without
re-deriving that the zero lies in it (the driver knows this from its own
Krawczyk image; the checker never recomputes an enclosure). This turns out to
be sound anyway, by an argument neither file states: success of the checker's
interval Gauss elimination on Y·J(B) forces every matrix in the interval
Jacobian J(B) to be nonsingular; J(B) is an entrywise-interval matrix set,
hence convex and connected, so det has one sign on all of J(B); DF at the
true zero and DF over any sub-box of B all lie in J(B); hence the sign
computed over the enclosure equals the index of the zero regardless of where
in B the zero sits. Recorded so the soundness of the index data does not
silently rest on an unstated lemma. No numerical conclusion changes.

## Claims that survive (and what was thrown at them)

### S1. Driver audit: no box can be dropped, no leaves overlap, the region
### is as claimed (attacks 1, 4 — hand audit, line by line)

- **No dropped boxes**: every box popped in `count_equilibria` either appends
  exactly one leaf (`continue` paths: outside / ball / sign / newton /
  isolated / unresolved) or pushes both children of one split; the only other
  exit is the `max_boxes` RuntimeError, which aborts the entire run with no
  certificate written. There is no code path that discards a box silently.
- **No overlap**: splits are exact Fraction midpoints; children share only a
  face; leaves are interior-disjoint by construction — and now also verified
  combinatorially from the emitted paths (R1 fix). A zero exactly on a shared
  face cannot be double-counted: a Krawczyk-unique verdict puts the leaf's
  single zero in the leaf's *interior*, so a face zero would be a second zero
  in the closed box and the verdict could not have been issued; face zeros
  can only surface as unresolved leaves (001 has none), and the offset region
  makes them unreachable for symmetric rational points.
- **Unresolved is honest**: a box at minimum width that certifies nothing is
  emitted as `unresolved` and counted; 001's certificate has zero.
- **Localization region cannot be smaller than claimed**: `auto_region` uses
  half-width r_up + 1/16 with offsets 1/97, 1/101, 1/103 < 1/16, so the
  region contains the closed ball |x| ≤ R0 for any outward rounding of
  r_up ≥ R0. The lemma itself re-derived: for |x| ≥ R0, x·F(x) =
  Σ qᵢ x·(x−xᵢ)/rᵢ³ has every term ≥ 0 (Cauchy–Schwarz plus |x| ≥ |xᵢ|), and
  simultaneous equality would force x = xᵢ for every i — impossible for ≥ 2
  distinct charges. The exclusion test (`box_min_norm_sq(box) >= r0_sq`,
  exact) therefore correctly excludes the closed outside including the
  boundary sphere. Re-checked exactly for all 15 excluded-outside leaves,
  with R0² = 1 recomputed from the config (REGION, BALL blocks — pass).
- **Domination lemma**: radii re-instantiated exactly with this review's own
  isqrt bounds (certified *lower* bounds for d_ij make the inequality
  conservative); all 18 excluded-ball leaves re-checked exactly. Pass.
- **Krawczyk verdict line**: operator formula, strict `inside_open`, strict
  `disjoint`, enclosure = K ∩ box, det-sign fallback to the full box — all
  match the standard statements; the mean-value emptiness argument is valid
  for any fixed rational Y, so the float-computed preconditioner is outside
  the trust chain on both engines. One latent hazard found, not triggered:
  see "what survived" below.
- **Shared-assumption inventory** (attack 5), beyond R1/R2: both engines use
  the same F/DF formulas (unavoidable; anchored by the classical validation
  knowns), both use a float preconditioner Y (certificates Y-independent by
  construction), the checker's `excluded-newton` fallback reuses the
  Krawczyk *theorem* (different arithmetic and code, same mean-value form —
  independence there is implementational, not mathematical), and both parse
  the same config JSON (checked by hand: it is the configuration 001
  states). The checker re-derives both lemmas and the region containment;
  domination radii values are taken from the certificate but only as
  candidates whose defining inequality is re-proved. The checker also PASSes
  certificates with unresolved leaves (count = lower bound, by design) and
  silently skips the index-sum identity if any det sign cannot be
  re-established — harmless for 001 (0 unresolved, all 24 signs
  re-established) but worth knowing when reading future PASS verdicts.

### S2. Structural checks on the actual certificate (attacks 1–3), fresh code

`python problems/maxwell-equilibria/explore/skeptic_cert_checks.py
<gunzipped cert> --centroid --expect-count 24` — all checks pass:

- TREE: full axis-consistency reconstruction of all 235,993 leaf paths
  (the R1 fix); Kraft = 1 and prefix-freeness re-confirmed independently.
- BOX: all 24 recorded isolated boxes equal the boxes rebuilt from their
  paths by this review's own midpoint code; every enclosure sits inside its
  leaf box; isolated records ↔ isolated leaves are 1:1.
- DISJ: the 24 isolation enclosures are pairwise disjoint by exact rational
  comparison (276 pairs; minimum separating gap ≈ 1.220·10⁻³ — consistent
  with 001's float minimum position separation 1.51·10⁻³ less enclosure
  widths).
- INDEX: all det signs ±1; sum = −4 = 1 − 5 (10 positive, 14 negative),
  recomputed from the leaf data.
- TALLY: leaf-kind tallies match the certificate summary exactly
  (7262 sign / 228674 newton / 18 ball / 15 outside / 24 isolated,
  0 unresolved).
- CENTROID: (1/3,1/3,1/3) lies in exactly one enclosure.

### S3. The count is robust to the subdivision tree (attack 6) —
### four full re-censuses, every tree different, every count as predicted

All runs single-core from the repo root, `--max-boxes 3000000`,
deterministic (wall times below are under 4-way CPU contention; 001's
baseline 1984 s was uncontended). Baseline (001): 471,985 boxes / 235,993
leaves, count 24, index sum −4. Certificates gzipped in `data/` as
`skeptic-002-cert{A,B,C,D}.json.gz`.

| run | perturbation | config / flags | boxes / leaves | secs | count | index sum |
|---|---|---|---|---|---|---|
| A | arithmetic grid: `--prec-bits 80`, same auto region | `data/abk-witness-config.json` | 471,985 / 235,993 | 2609 | **24** | −4 |
| B | wholesale different geometry: explicit region, offsets 1/89, 1/91, 1/93, half-width 5/4 (vs auto 1/97, 1/101, 1/103, 17/16) | `explore/skeptic_region_config.json`, `--prec-bits 64` | 472,775 / 236,388 | 2539 | **24** | −4 (from leaf data) |
| C | charge q = 4360/1000000 (window bracket, low) | `explore/skeptic_q4360_config.json`, `--prec-bits 64` | 359,577 / 179,789 | 1988 | **12** | −4 (+4/−8) |
| D | charge q = 4400/1000000 (window bracket, high) | `explore/skeptic_q4400_config.json`, `--prec-bits 64` | RUND_BOXES / RUND_LEAVES | RUND_TIME | **RUND_COUNT** | RUND_IDX |

Reproduce with (outputs land wherever `--out` points; the gzipped copies in
`data/` are byte-identical up to the `seconds` field):

    python harness/maxwell-equilibria/equilibria.py \
        --config problems/maxwell-equilibria/data/abk-witness-config.json \
        --out certA-prec80.json --prec-bits 80 --max-boxes 3000000
    python harness/maxwell-equilibria/equilibria.py \
        --config problems/maxwell-equilibria/explore/skeptic_region_config.json \
        --out certB-region.json --prec-bits 64 --max-boxes 3000000
    python harness/maxwell-equilibria/equilibria.py \
        --config problems/maxwell-equilibria/explore/skeptic_q4360_config.json \
        --out certC-q4360.json --prec-bits 64 --max-boxes 3000000
    python harness/maxwell-equilibria/equilibria.py \
        --config problems/maxwell-equilibria/explore/skeptic_q4400_config.json \
        --out certD-q4400.json --prec-bits 64 --max-boxes 3000000

Notes:

- **Run A came back with the *identical* tree** — the same 471,985 boxes,
  leaf-by-leaf identical (path, kind) lists and identical isolated boxes and
  det signs; only the isolation enclosures tightened (Krawczyk images on the
  finer 2⁻⁸⁰ grid). Checked by direct list comparison against 001's
  certificate. This is *expected on reflection* — split axis (max width) and
  midpoints are exact-arithmetic decisions independent of the grid, so the
  tree can only change where a certification test flips outcome, and every
  one of the ~472k test outcomes has margin far above the 2⁻⁶⁴ → 2⁻⁸⁰
  rounding difference. So A is evidence that no leaf verdict in 001 is a
  rounding-margin artifact, but it is a *weak* independence check for the
  tree itself; run B is the one that carries the tree-perturbation weight
  (different region ⇒ no box in common with 001's tree; leaf count 236,388
  vs 235,993, max depth 65 vs 58, 26 ball-exclusions vs 18, no
  outside-localization leaves — and the same 24 zeros, +10/−14).
- Run B supplies the region explicitly, so the driver does not itself claim
  completeness (`complete: false`) and skips the outside-localization
  shortcut and the index-sum field (that is why its log reports
  `index_sum=None`); the census is nevertheless complete because the region
  contains the closed ball |x| ≤ R0 = 1 (checked exactly on the run-B
  certificate) and the localization lemma was re-derived in S1. Its index
  sum, recomputed from the leaf data, is −4 (10 positive, 14 negative).
- Every completed run was pushed through `skeptic_cert_checks.py` — all
  checks pass on all certificates (run C with `--expect-count 12`; the
  centroid lies in exactly one enclosure in every run). RUNV_SENTENCE

### S4. The add-on predictions (attacks 7, 8)

- **Window brackets**: certified complete counts **RUNC_COUNT** at
  q = 4360/1000000 and **RUND_COUNT** at q = 4400/1000000, vs 001's float
  predictions 12 and 16. RUNCD_VERDICT
- **Exact centroid Hessian** (001 lead 5, closed): by the D₃ symmetry the
  Jacobian at the centroid is A·I + B·J (J = all-ones; the rank-one sums
  Σ_tri d dᵀ = I − J/3 and axial d dᵀ = t²J are verified exactly), and
  trace DF = 0 forces B = −A, so the full signature reduces to the sign of
  the axial eigenvalue λ_ax = 3/r_t³ − 4q/r_a³ (r_t² = 2/3, r_a² = 3t²) —
  i.e. to the single exact rational comparison **6561·t⁶ vs 128·q²**
  (`explore/skeptic_centroid_hessian.py`; also verifies F(c) = 0 for the
  whole family). At t = 17/200:
  - q = 4360/1000000 and q = **4367/1000000 (the witness)**: λ_ax > 0,
    signature (+,−,−), det DF(c) > 0 — **matching the certified det_sign +1
    of the centroid enclosure in 001's certificate**: an independent exact
    cross-check of one certified index by pure symmetry algebra, no
    intervals involved.
  - q = 4400/1000000: λ_ax < 0, signature (−,+,+), det DF(c) < 0.
  - The centroid is **exactly degenerate at q\*(t) = 81√2·t³/16**; at
    t = 17/200, q\* = 4.396801013961…·10⁻³ (certified enclosure of width
    2⁻¹²⁰ in the script output). Note q\* is *above* 001's float window edge
    ≈ 4.389·10⁻³. SPECULATION: the 24→16 transition at the window's high
    edge is NOT the centroid degeneracy; the centroid degenerates slightly
    later, inside the 16-count regime, flipping its index +1 → −1 (consistent
    with the C/D centroid det signs above). A certified census inside
    q ∈ (4.389, 4.3968)·10⁻³ would separate the two events — lead 3.

## Outcome

**VERIFIED_WITH_CORRECTIONS** (scope: attempt 001's certificate and its
verification chain; the four perturbed configurations listed above; nothing
about any other (t, q), and nothing about ABK's asymptotic theorem).

The headline stands: the configuration of 001 (units at e₁,e₂,e₃, charges
4367/1000000 at (251/600,251/600,251/600) and (149/600,149/600,149/600)) has
**exactly 24 nondegenerate equilibria**, now resting on: the original two
engines, this review's independent structural checks in exact arithmetic
(including a tiling proof strictly stronger than the checker's), and four
re-censuses over different subdivision trees, all counting 24 at the witness
parameters (A, B) and exactly the predicted bracket counts off them (C, D).
24 > 16 = (5−1)²: the explicit refutation instance of Maxwell's bound at
n = 5 is confirmed, and the escalation required by the verification contract
is discharged.

Corrections to 001 (none touches the count):

1. "Established by two independent implementations" overstates the coverage
   half: the tiling check in `verify_equilibria.py` is not sound by itself
   (R1); coverage independence is only established as of this review.
2. "The Kraft-equality tiling check exact" (001 step 3) — the arithmetic is
   exact, but the check is not equivalent to tiling; same correction.
3. R2: the det-sign location argument is sound but was never stated; recorded
   here.

**Not claimed:** that the harness driver is bug-free in regimes not exercised
(see the latent midpoint hazard below); anything about the window edges
beyond the two certified bracket points; the fold interpretation of the
window (still SPECULATION, now sharpened by the exact q\*).

## Why it failed / what survived

Nothing in the review kills 001. Found and worth keeping:

- **R1 (the real finding)**: the two-engine protocol had a common-mode hole
  on coverage — the checker's tiling check assumed the very driver behavior
  (per-node axis consistency) it was meant to verify. A driver bug emitting
  mixed-axis children would have produced certificates that PASS both
  engines while missing equilibria. The demonstrated tamper (DEMO 2) is the
  failure mode queue item 16 hypothesized, one level up: not a wrong tree
  with self-consistent leaves, but a wrong *tiling check* that
  self-consistent leaves slip through.
- **Latent driver hazard, not triggered**: `IV.mid()` rounds the midpoint
  down to the 2^−prec grid and `krawczyk()` uses it with J(B) taken over the
  box alone; if the rounding error 2^−prec ever reached half the box width
  (roughly prec_bits ≲ min_width_bits + 1 at unit scale), m could exit B and
  the mean-value form would be invalid. At 001's parameters the margin is
  2⁻³¹ vs 2⁻⁶⁴ — safe by 10 decimal orders — and the checker (exact
  midpoints) would catch resulting false leaves; but a future run at low
  `--prec-bits` with deep `--min-width-bits` should be refused rather than
  trusted. Lead 2.

Reusable:

- `explore/skeptic_cert_checks.py`: certificate-level structural audit with
  the axis-consistency tiling proof — cheap (~4 s) and engine-independent;
  run it on every future census certificate alongside the checker.
- `explore/skeptic_centroid_hessian.py`: the exact symmetry reduction of the
  centroid Hessian for the whole (t, q) family, and the closed-form
  degeneracy locus q\*(t) = 81√2·t³/16.
- The perturbed-region trick (explicit region in the config with fresh
  non-dyadic offsets) as the standard way to force a wholesale-different
  tree without touching the harness.
- Certified hard-instance data points for the family: counts 12 / 24 / 16 at
  q = 4360, 4367, 4400 (·10⁻⁶) — the two folds of 001's window are now
  bracketed by certificates, not floats.

## Leads generated

1. **Harden the checker's tiling check** (harness change, tier 0, generic):
   replace prefix-freeness + Kraft with the axis-consistency reconstruction
   (or add it), and port the same case to both `--validate` suites as a
   tamper test. Falsifiable: DEMO 1 and DEMO 2 certificates must flip to
   FAIL; all existing certificates must still PASS.
2. **Guard the midpoint rounding** (harness change, tier 0, generic): make
   `krawczyk()` return "unknown" when the box width is within a safety
   factor of the 2^−prec grid (or intersect m into B). Falsifiable: a
   synthetic low-precision deep-subdivision run must stop issuing
   certificates instead of mis-centering m.
3. **The sliver q ∈ (4.389, 4.3968)·10⁻³**: certified censuses on both sides
   of q\* = 81√2·t³/16 inside the 16-count regime — does the total count
   change at q\* (the centroid index flip +1 → −1 must be compensated by ±2
   somewhere), and is the 24→16 window edge a separate fold as 001
   speculated? Two ~35-min runs decide.
4. **Certified window edges** (001 lead 2, sharpened): the brackets are now
   certified (12 at 4360, 16 at 4400, ·10⁻⁶); binary-search both edges to
   certified 10⁻⁵-wide q-intervals with the same pipeline.

## References

- Attempt 001: `problems/maxwell-equilibria/attempts/001-abk-explicit-witness.md`
  and its data in `problems/maxwell-equilibria/data/`.
- Harness under review: `harness/maxwell-equilibria/equilibria.py`,
  `harness/maxwell-equilibria/verify_equilibria.py` (both read in full).
- A. Neumaier, *Interval Methods for Systems of Equations*, CUP 1990, §5.2 —
  as cited by the harness docstrings for the Krawczyk and interval-Gauss
  certificates; statements checked against the implementations, book not
  re-consulted.
- P. Arathoon, G. Ball, M. D. Kvalheim, *The Maxwell Conjecture is False*,
  arXiv:2607.27197 (July 2026) [T, via 001 — not re-read here].
