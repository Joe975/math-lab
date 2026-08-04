# 001 — Explicit certified witness for the five-charge Maxwell counterexample

- **Problem:** maxwell-equilibria, `problems/maxwell-equilibria/PROBLEM.md`
- **Date:** 2026-08-04
- **Mode:** informed
- **Type:** computational search + certification
- **Tools:** `problems/maxwell-equilibria/explore/abk_witness_scan.py` (float
  candidate census; deterministic, seeded RNG, ~10 s per parameter point);
  `harness/maxwell-equilibria/equilibria.py` (certified complete census;
  exact/outward-rounded arithmetic at 64 dyadic bits, run at
  `--max-boxes 3000000`); `harness/maxwell-equilibria/verify_equilibria.py`
  (independent re-verification). Runtimes in the body.
- **Sources:** arXiv:2607.27197 (Arathoon–Ball–Kvalheim, *The Maxwell
  Conjecture is False*, July 2026) [T — abstract and HTML body via web
  fetch, not the PDF]; arXiv:2607.28785 (three-charge bound 12 → 6) [T];
  both consumed as labelled inputs per STATUS queue item 16.

## Approach

The ABK refutation proves that five positive charges can have ≥ 24
nondegenerate equilibria — but only asymptotically: their configuration
family is valid "for all sufficiently small ε" with no explicit ε named, and
their supporting computations are computer-algebra, not certified. Queue
item 16: produce an *explicit rational* configuration and certify its exact
equilibrium count, making the counterexample a concrete checkable object.

Why this rather than the obvious alternative (certifying at their literal
coordinates): their triangle has vertices at irrational coordinates
(circumradius 1 in the xy-plane), while the equilateral triangle embeds
rationally in ℝ³ as e₁, e₂, e₃ (pairwise distance √2, symmetry axis along
(1,1,1)). Since equilibrium counts are similarity-invariant and charge
values are scale-invariant (F_λ(λx) = λ⁻²F(x)), the family maps to

    unit charges at e₁, e₂, e₃;  charges q at (1/3,1/3,1/3) ± t·(1,1,1)

with everything in ℚ once t, q are rational. The dictionary to ABK's frame:
ε = 3t/√2, and their charge law is q_ε = (3/4)ε³ − (5/32)ε⁵.

## What was done

**Statement under test.** The configuration

    q₁ = q₂ = q₃ = 1 at (1,0,0), (0,1,0), (0,0,1);
    q₄ = q₅ = 4367/1000000 at (251/600, 251/600, 251/600)
                        and (149/600, 149/600, 149/600)

(t = 17/200, i.e. ε ≈ 0.18031; q within 2·10⁻⁸ of the ABK law value)
has ≥ 17 isolated equilibria — refuting Maxwell's (5−1)² = 16 — and in fact
exactly 24, all nondegenerate, matching ABK's count.

**Step 1 — parameter hunt (floats; locate, never certify).**
`explore/abk_witness_scan.py --scan` sweeps (t, λ) with q = λ·q_ε. Findings
of record, all at float precision with a seeded global Newton census
(cluster seeds around the centroid + global seeds in the localization
ball), reproducible via the flags in the file:

- The full 24-count window exists only for t ≤ 0.085 (ε ≲ 0.18): at
  t = 0.0875 the maximum observed count is 18, and at t ≥ 0.1 the counts
  run 6 → 18 → 10 as q crosses the ABK law with no 24 anywhere on the
  sampled grid. "Sufficiently small ε" is doing real work — at ε ≈ 0.19
  the bifurcation has not yet released all 21 equilibria.
- At t = 17/200 the 24-window in q is roughly (4.366, 4.389)·10⁻³ — about
  half a percent wide, λ ∈ (1.000, 1.005) — with conditioning best at the
  low edge and count falling 24 → 16 past the high edge. Below the window
  the count is 12.
- Chosen witness: q = 4367/1000000, the low-edge sweet spot, which happens
  to sit within 2·10⁻⁸ of the exact ABK law value (3/4)ε³ − (5/32)ε⁵ at
  ε = (3/√2)(17/200). Float diagnostics there: 24 zeros, minimum pairwise
  separation 1.51·10⁻³ (centroid to each of the three nearest cluster
  points), minimum |det DF| ≈ 1.3·10⁻⁵, float index sum −4 = 1 − n.
- Structure of the 24 (float positions, D₃ₕ-symmetric): the exact centroid
  (1/3,1/3,1/3) — a rational equilibrium by symmetry, at all (t, q) — plus
  17 more within 0.0152 of it (orbits of sizes 3+3+6+2+3), plus two outer
  orbits of 3 (|z−c| ≈ 0.110 and 0.147), the latter being the persisted
  edge equilibria of the bare triangle.

**Step 2 — certified complete census.** Command:

    python harness/maxwell-equilibria/equilibria.py \
        --config problems/maxwell-equilibria/data/abk-witness-config.json \
        --out problems/maxwell-equilibria/data/abk-witness-cert.json \
        --prec-bits 64 --max-boxes 3000000 --progress-every 20000

All five charges are positive, so the localization lemma applies and the
census is complete (search region ⊇ ball of radius max|xᵢ| = 1): the
certificate covers every equilibrium of the configuration, not a sampled
region.

Result: **exactly 24 isolated equilibria, zero unresolved leaves**, in
471,985 boxes / 235,993 leaves / 1984 s (single core, 64 dyadic bits,
minimum width 2⁻³⁰ never reached — deepest leaf at 58 total splits ≈ 19
per axis). Leaf breakdown: 7,262 sign-exclusions, 228,674 Krawczyk-empty
exclusions, 18 charge-ball exclusions, 15 outside-localization exclusions,
24 Krawczyk isolation certificates. Every one of the 24 has a certified
det-sign (10 positive, 14 negative — all nondegenerate); enclosure widths
run 4.7·10⁻⁷ to 1.9·10⁻⁴. Certified positions agree with the float census
to ~9 digits. Artifacts in `problems/maxwell-equilibria/data/`:
`abk-witness-cert.json.gz` (full 34 MB certificate, 680 KB gzipped —
gunzip before re-verifying), `abk-witness-summary.json` (counts, lemma
data and the 24 enclosures, human-readable), and the run is
deterministic, so the certificate regenerates exactly from
`abk-witness-config.json` with the command above.

**Step 3 — independent re-verification.** Command:

    python harness/maxwell-equilibria/verify_equilibria.py \
        --cert problems/maxwell-equilibria/data/abk-witness-cert.json

Different arithmetic (directed-rounding decimal vs dyadic rational),
different isolation certificate (preconditioned interval Newton / interval
Gauss vs Krawczyk), independently re-derived lemmas, exact tiling check.

Result: **PASS, zero failures, 924 s.** All 235,993 leaves re-established:
the 24 isolation certificates by preconditioned interval Newton, all
235,969 exclusions, both lemmas re-derived from the configuration, the
Kraft-equality tiling check exact, det-signs matching, and the index-sum
identity re-confirmed.

**Cross-checks.**

- Index-sum identity: a complete nondegenerate census of 5 positive charges
  must have Σ sign det DF = 1 − 5 = −4. Holds: 10 positive and 14 negative
  certified det-signs, in both engines.
- The float census (independent implementation, different algorithm) found
  the same count and positions to ~9 digits before any certification ran.
- Harness validation suite (classical two-charge, collinear, equilateral
  knowns) passes before and after the run.

## Outcome

**VERIFIED** (scope: the single explicit configuration above, and nothing
else). The configuration

    q₁ = q₂ = q₃ = 1 at e₁, e₂, e₃;
    q₄ = q₅ = 4367/1000000 at (251/600,251/600,251/600), (149/600,149/600,149/600)

has **exactly 24 isolated equilibrium points, all nondegenerate** — a
complete certified census over a region provably containing every
equilibrium, established by two independent implementations using different
arithmetic and different fixed-point theorems. Since 24 > 16 = (5−1)², this
is an explicit, machine-checkable refutation of Maxwell's conjecture for
n = 5, and an independent confirmation of the count in arXiv:2607.27197 at
a concrete witness the preprint itself does not supply (its theorem is
asymptotic in ε with no ε₀ named).

Per the verification contract this is a record-adjacent count: it passed
both routes and is hereby **reported as requiring escalation** — a separate
skeptic-review attempt (fresh eyes on the harness itself, not just the
certificate) is queued before the ledger treats the refutation as settled
library fact.

**Not claimed:** that 24 is the maximum for five charges; anything about
any other (t, q), including ABK's asymptotic theorem itself (one instance
is certified, not the family); the float-level window observations in
step 1 (those are `EVIDENCE` at float precision, sampled grids only); and
nothing about the n = 3 gap, which stays open at 4-vs-6. The ABK preprint
is six days old and not peer-reviewed; this attempt confirms its headline
count at one explicit point, which is evidence for the paper, not a review
of its proof.

## Why it failed / what survived

Nothing failed. What remains unproven, labelled:

- The parameter-window observations (24-count window exists only for
  t ≤ 0.085, window ≈ half a percent wide in q) are float-level
  `EVIDENCE`; certified window boundaries are lead 2.
- `SPECULATION`: the 24-count window's low edge is a fold bifurcation
  (12 → 24, six zeros born in three symmetric pairs) and the high edge
  another (24 → 16). Consistent with the observed count sequence and
  separations; not certified.

Reusable:

- The witness configuration itself: five rational charges realizing 24
  equilibria — the hard-instance family for any future counting work.
- The rational equilateral embedding trick (e₁,e₂,e₃ + diagonal axis)
  turns any D₃-symmetric planar-triangle construction rational; it will
  serve every future triangle-based configuration here.
- Cost calibration for the harness at a near-degenerate instance:
  ~472k boxes / 33 min to certify a cluster with min separation 1.5·10⁻³
  and min |det DF| ≈ 1.3·10⁻⁵ at 64 bits; verification ~15 min. Deepest
  isolation at ~19 splits per axis. This prices lead 2 and queue item 17.
- The exploration scanner (`explore/abk_witness_scan.py`) as the
  candidate-hunting front end for any configuration family.

## Leads generated

1. **Skeptic review of this attempt** (required by the contract's
   escalation clause): independently audit the harness's Krawczyk and
   covering logic — the two engines share the subdivision tree, so a bug
   in the *driver* (not the certificates) is the residual common-mode
   risk. Concrete check: re-run the census at different precision,
   min-width, and box offsets (the tree changes wholesale; the count must
   not), and verify the 24 enclosures pairwise-disjoint by exact
   comparison. Outcome either way is recordable.
2. **Certify the window boundaries.** Run the census at q = 4360/1000000
   (expect 12) and q = 4400/1000000 (expect 16): certified counts on both
   sides bracket the two folds. Cost: ~2 × 35 min by the calibration
   above. Turns the float window into certified bifurcation brackets.
3. **Certified ε₀ probe.** Binary-search t ∈ (0.085, 0.0875) with
   certified censuses: the largest t with count 24 gives a certified
   lower bound on how far ABK's "sufficiently small" actually reaches —
   a quantitative datum the preprint leaves open.
4. **n = 4 variant.** Whether a triangle plus ONE axial charge (n = 4,
   Maxwell bound 9) can beat 9 by the same bifurcation mechanism.
   Falsifiable by the same pipeline: scan, then certify. If it can, the
   refutation extends to n = 4; if not, the mechanism's n-dependence is
   itself informative.
5. **Exact centroid equilibrium.** (1/3,1/3,1/3) is an equilibrium of
   every member of the family by symmetry — a rational point where the
   full Hessian is computable in ℚ exactly. The fold structure of lead 2
   should be decidable in exact arithmetic at that point alone
   (eigenvalue sign changes of a rational matrix), no census needed.

## References

## References

- P. Arathoon, G. Ball, M. D. Kvalheim, *The Maxwell Conjecture is False*,
  arXiv:2607.27197 (July 2026). [T]
- *From 12 to 6: Sharpening the Three-Charge Bound in Maxwell's Problem*,
  arXiv:2607.28785 (July 2026). [T]
- A. Gabrielov, D. Novikov, B. Shapiro, *Mystery of point charges*, Proc.
  London Math. Soc. 95 (2007). [T — abstract level only]
- Ya-Lun Tsai, *Maxwell's conjecture on three point charges with equal
  magnitudes*, Physica D 309 (2015). [T — abstract level only]
- This repo: `problems/maxwell-equilibria/PROBLEM.md` (verification
  contract), STATUS queue item 16 (route framing).
