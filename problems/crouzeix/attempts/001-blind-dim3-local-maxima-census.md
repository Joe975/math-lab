# Local-maxima census of the Crouzeix ratio, n = 3, deg p <= 3

```
problem: crouzeix
date: 2026-08-04
mode: blind
source-commit: a0bfd6ed660fb4ab8ced4158ba491c87a4667cd8
type: computational (optimization-landscape census; certified enclosures)
tools: explore/census.py, explore/classify.py, explore/survivors_detail.py,
       explore/certify_endpoints.py, explore/measure_cost.py (all pure
       stdlib), harness/crouzeix/ratio.py, harness/crouzeix/verify_ratio.py
sources: PROBLEM.md; published literature via web search (Greenbaum-Overton
         papers only; nothing about this lab)
```

## Approach

Question: starting local maximization of the Crouzeix ratio
`f(A,p) = ||p(A)|| / max_{z in W(A)} |p(z)|` from a designed sample of
(matrix, polynomial) configurations at the fixed slice **n = 3, deg p <= 3**,
do all local maxima terminate at already-known extremal structure, or does
anything new appear?

Why this: the conjecture is exactly the claim that the global maximum of f is
2. A census of *local* maxima at the smallest open dimension is the cheapest
way to (a) look for unexpected near-extremal structure and (b) produce
certified, reproducible landscape evidence, which the published numerical
work does not carry.

**How this differs from the published Greenbaum–Overton experiments** (their
LAA 2018 paper and the Calcolo 2021 follow-up; see References — stated up
front because reproducing their finding is not a result):

1. **Fixed slice.** GO sweep many dimensions and degrees; this census fixes
   n = 3, deg <= 3 and spends the budget on basin coverage of that one slice.
2. **Structured start families.** GO initialize from random (2021:
   heavy-tailed) matrices and polynomials. Here 7 named start families are
   used, five of them deliberately placed at or near non-generic structure:
   perturbed J3 with p near z^2; perturbed J2 (+) [c] with p near z;
   near-normal A with p-roots engineered between eigenvalues; companion-like
   strongly coupled triangular A with lacunary p = z^3 + noise; elliptic-W
   configurations with p-roots pushed just outside the boundary at spread
   angles (multiple near-coincident |p| peaks); plus a real slice and a
   random baseline.
3. **Optimizer.** Derivative-free (Nelder–Mead restarts + compass polish on a
   20-real-parameter Schur-form chart) instead of Chebfun + BFGS /
   gradient-sampling. Different method, different stall behavior — an
   endpoint here is only called a candidate local maximum after an explicit
   *escape probe* (240 random perturbations at 4 scales + one fresh
   Nelder–Mead restart) fails to ascend.
4. **Certified endpoints.** Every census endpoint, rounded entrywise to
   denominator 256, gets a rational-arithmetic certified enclosure of its
   ratio from the tier-0 harness (`ratio.py`), and the top endpoints are
   re-verified by the independent route (`verify_ratio.py`,
   Faddeev–LeVerrier + Sturm). GO's computations are floating-point
   throughout.

Parametrization note: the ratio is invariant under unitary similarity, so A
is kept upper triangular (Schur form, WLOG) — 6 complex entries; p has 4
complex coefficients, normalized to unit norm inside the objective (pure
gauge). 20 real parameters total.

## What was done

All commands run from the blind working copy root.

1. **Harness self-tests** (both pass):

   ```
   python3 harness/crouzeix/ratio.py --selftest          # 1.8 s
   python3 harness/crouzeix/verify_ratio.py --selftest   # 2.6 s
   ```

2. **Cost measurement before sizing** (`explore/measure_cost.py`): one
   certified enclosure of a dense 3x3, denominator-32 configuration with a
   degree-3 polynomial costs **0.35 s** (tol 1e-6, 16 half-planes) to
   **0.82 s** (tol 1e-9, 32 half-planes); denominator-256 endpoints cost
   about 1 s. So certifying every endpoint of a few-hundred-point census is
   affordable; the float exploration dominates the budget instead
   (~19 s per optimization start at ~4000 objective evaluations, three
   workers running concurrently beside two sibling agents).
   Census sized accordingly: target ~300 starts. RECORDED NARROWING: the
   original design wanted ~600+ starts; concurrency made the per-start cost
   ~2x the pilot estimate, so the interleaved run was cut at ~240 records
   (balanced across families by interleaving) plus 39 records from a first
   truncated launch (family-ordered, so family-biased; kept, labeled).

3. **Float validation.** The float ratio evaluator reproduces the literature
   values at this slice: (J3, z) -> 1.41421 (sqrt 2), (J3, z^2) -> 2.0000,
   (J2 (+) [0.2], z) -> 2.0000, normal matrix -> 1.0.

4. **Known-structure anchors, certified by both routes** (calibration for
   what "certified near 2" can look like at 80 half-planes):

   ```
   (J3, p=z^2):        ratio in [1.980557335, 2.25],  both routes AGREE
   (J2 (+) [1/4], p=z): ratio in [1.990266652, 2.0],  both routes AGREE
   ```

5. **Census.**

   ```
   python3 explore/census.py --run 25 --seed 4 --out data/census_i4.jsonl  # x3 workers, seeds 4,5,6
   # earlier truncated launch: --run 30 --seed {1,2,3} -> data/census_w{1,2,3}.jsonl
   ```

   286 starts total (241 interleaved + 45 from the truncated first launch;
   per-family: random 81, jordan2 34, jordan3 36, nearnormal 34, companion
   34, peaks 34, real 33 — random is overweighted by the truncated launch).
   Each start = NM(step .2) -> NM(step .05) -> compass polish at 192
   boundary samples -> diagnostics at 512 boundary samples.

6. **Classification + escape probe.**

   ```
   python3 explore/classify.py data/census_*.jsonl   # -> data/classified.jsonl
   ```

7. **Certification of every endpoint** (den 256, tol 1e-6, 32 half-planes)
   and independent verification of the top endpoints:

   ```
   python3 explore/certify_endpoints.py data/classified.jsonl --den 256 \
       --out data/certified.jsonl --top-verify 8
   ```

8. **Structural detail of probe survivors.**

   ```
   python3 explore/survivors_detail.py data/classified.jsonl below2-max
   ```

## Outcome

**Status: EVIDENCE** — scoped to dimension 3, polynomial degree <= 3, and
the search design above (286 starts across the 7 named families, escape
probing as the local-maximality test). Nothing here is a claim about the
conjecture.

**Census verdict: every surviving local maximum sits at already-known
extremal structure. Nothing new appeared.**

Endpoint classes (286 starts, after escape probing; full data in
`data/classified.jsonl`, per-endpoint certified enclosures in
`data/certified.jsonl`):

| class | n | meaning |
|---|---|---|
| stall | 175 (61%) | escape probe found ascent: an optimizer artifact, not a local maximum |
| near-2 (double eig) | 48 | ratio > 1.97, two eigenvalues coalescing: J2 (+) [c]-type basin |
| near-2 (other) | 44 | ratio > 1.97, eigenvalues still separated, W(A) near-disc (median disc residual 0.03) |
| near-2 (triple eig) | 2 | ratio > 1.97, all eigenvalues coalescing: J3-type basin |
| ratio exactly 1 | 13 | genuine local maximum at value 1 |
| strictly between 1 and 1.97 | 4 -> 0 | all four fell to deep probing (see below) |

- **The 13 ratio-1 endpoints** all show one eigenvalue of A lying *on* the
  boundary of W(A) with the single |p|-peak at it (boundary margins < 1e-3
  of diam W in `explore/survivors_detail.py` output). This is the
  ice-cream-cone configuration Greenbaum-Overton describe for their
  locally-minimal value 1 — a **rediscovery, with citation** (Calcolo 2021).
- **The 97 near-2 endpoints** have near-circular W(A) (median disc residual
  0.002-0.03 across the three subclasses) and climb toward the known
  Crouzeix pairs; the eigenvalue-coalescence classes match (J2 (+) [c], z)
  and (J3, z^2). The "other" subclass has not yet collapsed its eigenvalues
  — consistent with the extremal manifold being approached along
  non-defective matrices, with coalescence only in the limit.
- **The 4 intermediate candidates** (float ratios 1.96764, 1.96026, 1.95884,
  1.78117) all **escaped under deep probing** (`explore/deep_probe.py`:
  512-sample objective, 1200 random perturbations at 5 scales, 3 jittered
  NM restarts, fine compass) with ascents of 0.005-0.006. Verdict: stalls.
- **Certified enclosures: all 286 endpoints** (rounded to denominator 256)
  were certified by `ratio.py` (tol 1e-6, 32 half-planes, ~0.75 s/point,
  214 s total): **0 refutations**; the largest certified lower bound on any
  endpoint ratio is **1.969725** (compare the same pipeline's certified
  lower bound of 1.974 at the exact 2x2 Jordan pair — i.e. the near-2
  endpoints are as close to 2 as the 32-gon slack allows one to certify).
  The **top 8 endpoints by float ratio were re-verified by the independent
  route** (`verify_ratio.py`: Faddeev-LeVerrier + Sturm): all 8 AGREE, and
  both routes agree the refutation test fails.
- **Not seen:** the published intermediate locally-maximal values at this
  slice — reciprocals 0.698 and 0.844, i.e. ratios ~1.433 and ~1.185
  (values `[L]`-transcribed from a machine summary of arXiv:2105.14176; not
  read from the paper itself) — never appeared as endpoints of any of the
  286 basins. No start family here reproduced them.

## Why it failed / what survived

This was a census, and the dead-end is the finding:

- **All basins from this design terminate at known structure** (ratio-2
  Jordan-type pairs and ratio-1 ice-cream-cone points). No new extremal
  structure at n = 3, deg <= 3 from these seven families.
- **The design did not reach the published intermediate maxima.** GO's
  intermediate values came from ~half a million heavy-tailed random starts;
  286 structured starts found none of those basins. So this census can say
  nothing about them beyond: their basins, if genuinely locally maximal,
  are small relative to this start distribution. That is the main
  limitation of the design, and it was forced by the measured cost
  (recorded narrowing in "What was done", item 2).
- **What survived methodologically: the stall rate.** 61% of derivative-free
  maximization runs terminate at points that *look* like interior local
  maxima (ratios spread over [1.0, 1.97], certified lower bounds up to
  1.94) but escape under deeper probing. A census that skipped the escape
  probe would have reported dozens of fake "new local maxima". Any future
  landscape claim on this problem should be required to carry an
  escape-probe pass; the 4 deep-probed examples here are the concrete
  cautionary cases.
- SPECULATION: at n = 3, deg <= 3, the only locally maximal values of the
  Crouzeix ratio reachable with non-negligible basin measure from generic
  or J-structured starts are exactly 1 and 2. (The published intermediate
  values would then be either heavy-tail-accessible only, or nonsmooth
  stationary values that a stronger escape test might also dissolve —
  untested here.)

## Leads generated

1. **Attack the published intermediate maxima directly** (falsifiable):
   build starts matching GO's description of their intermediate minimizers
   — eigenvalues well separated, the three roots of p nearly coincident,
   several active |p|-peaks — and run the deep-probe test on the resulting
   endpoints. If they escape, there is a discrepancy with the published
   local-minimality claim worth escalating (first re-checking the paper
   itself, not the `[L]` transcription; informed-side task).
2. **Tighter certification of a single near-2 endpoint**: den 1024,
   dir_bound 7+, tol 1e-12 should certify an endpoint ratio inside
   [1.99, 2.001]; cost measured ~5-10 s/point, entirely affordable.
3. **Ridge geometry of the stalls**: the 175 stalls concentrate where the
   active set of max_{dW}|p| changes cardinality (3-4 near-coincident
   peaks). Characterizing that ridge could explain which nonsmooth
   stationary values exist between 1 and 2 — and whether GO's intermediate
   values live on it.
4. **Harness lead**: an escape-probe utility (random multi-scale + fresh-
   simplex restart) is generic for any optimization-landscape claim in this
   repo and could move into `harness/` as tier-0 tooling.

## References

- M. Crouzeix, C. Palencia, *The numerical range is a (1+sqrt 2)-spectral
  set*, SIAM J. Matrix Anal. Appl. 38 (2017).
- A. Greenbaum, M. L. Overton, *Numerical investigation of Crouzeix's
  conjecture*, Linear Algebra Appl. (2018); preprint
  optimization-online 5703 (2016).
- A. Greenbaum, M. L. Overton, *Local minimizers of the Crouzeix ratio: a
  nonsmooth optimization case study*, Calcolo 58 (2021); arXiv:2105.14176.
  Key published finding this census must be compared against: minimizing the
  reciprocal ratio ||p||_W / ||p(A)|| from heavy-tailed random starts, they
  found (i) 1 is a frequently occurring locally minimal value, (ii) many
  locally minimal values strictly between 0.5 and 1 exist, often the same
  over real and over complex parameters, and (iii) the smallest locally
  minimal value found is always 0.5 — i.e. in this record's convention,
  local maxima of the ratio at 1, at many values strictly between 1 and 2,
  and at 2, never above.
- Tier-0 harness: `harness/crouzeix/ratio.py`, `harness/crouzeix/verify_ratio.py`.
