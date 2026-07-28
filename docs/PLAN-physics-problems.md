# Plan: onboard three physics-flavored problems

Status: **plan only — nothing below is implemented.** Written 2026-07-28 as a
handoff for a fresh session. Branch: `claude/physics-math-problems-pgg8vx`.

Scoping decision (already made, do not re-litigate): add three problems whose
open frontiers are attackable with this lab's methodology — exact-arithmetic
censuses, falsifiable first steps, honest `EVIDENCE` framing. Candidates that
lost: hot-spots conjecture (certified PDE eigenfunction numerics in stdlib
Python is a bigger project than the problem work), Lieb–Thirring sharp
constants (open frontier is analytic; computation mostly reproduces knowns),
Saari's conjecture n≥4 (polynomial systems beyond our exact-algebra reach).

This file is **not tier 0** (see `tiers.json`): it contains this lab's framing
and route ideas, so it must never be copied into a blind checkout or cited in
`PROBLEM.md` prose.

## The three problems

### 1. `billiards-triangles` — periodic orbits in triangular billiards

- **Conjecture.** Every triangle admits a periodic billiard orbit. Settled for
  acute, right, rational-angled, and isoceles triangles; Schwartz's
  computer-assisted work covers all obtuse triangles with angles < 100°.
  Post-2008 extensions of that frontier are reported in the literature —
  **verify the exact current frontier during PROBLEM.md drafting** (search:
  R. Schwartz, obtuse triangular billiards; McBilliards) and cite precisely.
- **Physics lens.** Classical mechanics / dynamical billiards (new `dynamics`
  attack surface for the lab; the field tag exists in `mechanisms.json` but
  is thin on physical systems).
- **Key structure.** Unfolding: an orbit with bounce word `w` exists iff an
  exact linear-algebra condition on the angles holds, and one word certifies
  an open *region* of the (α, β) parameter triangle. The conjecture becomes a
  covering problem: tile the obtuse region with word-certificates. Same shape
  as the lonely-runner census: word ↔ tuple, region-certificate ↔ exact ML.
- **Harness** (`harness/billiards-triangles/`):
  - *unfold.py* — word → unfolding transformation → exact certificate that a
    word yields a periodic orbit on an angle region. Exact arithmetic over
    ℚ(angles) or exact rational sampling + interval bounds; standard library
    only (`fractions.Fraction`).
  - *verify_cover.py* — independent checker: takes a claimed (word, region)
    list and re-verifies by a *different* method (direct orbit simulation
    with interval arithmetic at region sample points plus boundary checks).
  - Verification contract in PROBLEM.md: a coverage claim states the angle
    region, the certificate list, and the arithmetic used; floating point is
    never sufficient for a certificate; partial coverage is `EVIDENCE` about
    the covered region only.
- **First attempts (queue lines):**
  1. Self-test: re-derive coverage of the right/isoceles/acute cases with our
     own implementation (this validates the harness; expected `VERIFIED`
     with scope = the re-derived classes).
  2. Word census by length L: which short words certify which obtuse
     regions; map where coverage stalls approaching and past 100°. The
     geometry of the uncovered set is a `MAP`-grade deliverable.
  - Kill condition to record at queue time: a region needing certificate
    words whose length grows super-exponentially as angle → frontier kills
    the finite-census route (that is itself worth recording).

### 2. `mahler-4d` — Mahler's volume-product conjecture in ℝ⁴

- **Conjecture.** Among centrally symmetric convex bodies K ⊂ ℝⁿ, the volume
  product vol(K)·vol(K°) is minimized by Hanner polytopes, value 4ⁿ/n!.
  Proved n ≤ 3 (n=3: Iriyeh–Shibata 2019). Open n ≥ 4. Known: Hanner
  polytopes are local minima (Kim; Kim–Reisner); Bourgain–Milman gives the
  lower bound cⁿ·4ⁿ/n!.
- **Physics lens.** Symplectic geometry / Hamiltonian mechanics: Mahler is
  implied by Viterbo's symplectic capacity conjecture, and Viterbo was
  **refuted** by Haim-Kislev–Ostrover (2024). Mahler survived its physics
  motivation — fresh landscape. When onboarded, run `/ripple` on the Viterbo
  refutation as an early informed-side exercise.
- **Key structure.** For rational polytopes everything is exact: polarity,
  triangulation, volume all live in ℚ. `fractions.Fraction` suffices; C
  kernel (`cc -O2`, repo convention) for the census inner loop. This is the
  best exact-arithmetic fit of the three.
- **Harness** (`harness/mahler-4d/`):
  - *polytope.py* — exact rational convex hull, polarity, volume by
    triangulation in ℝ⁴ (dimension-generic where cheap).
  - *verify_product.py* — independent re-computation of a claimed volume
    product by a different triangulation/ordering (adversarial re-run tool).
  - Verification contract: volume products exact in ℚ or not claimed; any
    census states its generator universe (e.g. vertex coordinates in
    {0,±1}, ≤ V vertices) and counts; a candidate below 4⁴/4! = 32/3 must be
    re-verified by *verify_product.py* before *any* ledger mention.
- **First attempts (queue lines):**
  1. Self-test: reproduce vol·vol° for cube, cross-polytope, ball
     approximants, and the ℝ³ knowns with our tooling.
  2. Exact census over small symmetric 0/±1-vertex polytopes in ℝ⁴: anything
     within ε of 32/3 that is not Hanner? (`EVIDENCE`, scope = universe.)
  3. Rational-perturbation check of local minimality at the two Hanner types
     in ℝ⁴ (independent confirmation of Kim's result = harness validation).
  - Kill/win condition: any symmetric polytope with product < 32/3 refutes
    Mahler outright — every search is kept honest by that.

### 3. `crouzeix` — Crouzeix's conjecture

- **Conjecture.** For every square complex matrix A and polynomial p,
  ‖p(A)‖ ≤ 2·max{|p(z)| : z ∈ W(A)} (W = numerical range / field of values).
  Known with constant 1+√2 (Crouzeix–Palencia 2017); constant 2 is attained
  in the limit (Crabb/Choi-type examples), proved for 2×2, and for various
  structured classes. Open in general.
- **Physics lens.** Operator/spectral theory; W(A) is the set of quantum
  expectation values; non-normal operator behavior. Adds an operator-theory
  attack surface to the lab.
- **Key structure.** Fully finite-dimensional optimization + structure
  census: local maxima of the Crouzeix ratio over (A, p) in dims 3–4,
  degrees ≤ 3. A ratio > 2 anywhere refutes the conjecture. Known
  experimental work (Greenbaum et al.) suggests all maxima flow to the known
  near-extremizers — our census asks that question independently.
- **Risk item.** Weakest exact-arithmetic fit: ‖p(A)‖ and W(A) need
  eigen/singular values. Contract must be built on *certified enclosures*:
  ratio claims as intervals via standard-library interval arithmetic
  (Gershgorin/residual bounds around floating-point candidates). Scope the
  harness spike (phase 2 below) before committing the census design; if
  certification in stdlib proves too heavy, narrow attempt 2 to structured
  families where ‖p(A)‖ has closed form, and record that narrowing in the
  attempt rather than silently shrinking scope.
- **Harness** (`harness/crouzeix/`):
  - *ratio.py* — compute the Crouzeix ratio with certified error bounds
    (Jacobi eigensolver for Hermitian parts, power iteration + residual
    certificates for norms; no numpy — repo is stdlib-only).
  - *verify_ratio.py* — independent re-computation (different algorithm,
    e.g. characteristic-polynomial root isolation with rational interval
    arithmetic) for any near-extremal claim.
  - Verification contract: any reported ratio is an interval with a proved
    enclosure; a refutation claim (ratio > 2) requires the interval's lower
    bound > 2 under *verify_ratio.py*; optimization landscapes are
    `EVIDENCE` scoped by dimension/degree/search design.
- **First attempts (queue lines):**
  1. Self-test: reproduce the 2×2 theorem numerically and the known
     near-extremal families' ratios → validates certification machinery.
  2. Local-maxima census, dim 3, degree ≤ 3: do all basins terminate at
     known extremal structure? (`EVIDENCE`.)
  - Kill condition: if certified enclosures cost more than ~10× the float
    computation, the wide census dies; fall back to structured families.

## Implementation phases (for the executing session)

Work informed-mode for the infrastructure; the *first attempts* on each new
problem should run **blind by design** — with no prior art yet, blind is
nearly free and seeds the anchoring dataset from attempt 001.

**Phase 0 — literature check (no repo changes).** Verify the published status
claims above before writing them into tier 0: billiards frontier post-2008;
Mahler ℝ⁴ status and local-minimality citations; Crouzeix known classes and
best constant. Anything in `PROBLEM.md` must be literature, nothing of ours
(tier-0 rule; CI-enforced leak scan).

**Phase 1 — scaffolding (one commit per problem).** For each problem slug:
- `problems/<slug>/PROBLEM.md` — statement, published status, and a hard
  **verification contract** section (model: `problems/lonely-runner/PROBLEM.md`).
- `problems/<slug>/prior-art.json` — empty-attempts index:
  `{"problem": "<slug>", "route_status": null, "route_summary": "no attempts yet", "attempts": []}`
  — validate the exact empty-index shape against `docs/prior-art.schema.json`
  and `tests/test_prior_art_schema.py` before committing.
- `problems/<slug>/PRIOR-ART.md` — stub ("No attempts yet.").
- `problems/<slug>/attempts/` `data/` — empty (add `.gitkeep` only if other
  problems do; match convention).
- No `tiers.json` change needed: existing `problems/*/PROBLEM.md` and
  `harness/**/*` globs cover the new slugs automatically. Confirm
  `scripts/blind.sh <slug> /tmp/x` works for each.

**Phase 2 — harness skeletons (one commit per problem).** The files listed
per problem above, each with a working self-test entry point but no census
runs yet. Check `harness/lonely-runner/` for structure conventions
(reference implementation + independent analyzer pair). Standard library
only; C kernels optional and deferred until a census actually needs speed.
Harness docstrings must not reference lab findings
(`test_harness_has_no_narrative_backreferences`).

**Phase 3 — registry + ledger (single commit).**
- `mechanisms.json`: add field lenses only if needed — `dynamics` and
  `geometry-topology` already exist and cover billiards/Mahler; consider
  adding `operator-theory` ("Spectral theory, numerical ranges, functional
  calculus, non-normal operator behavior.") for Crouzeix. New mechanism
  tags get added by attempts, not preemptively (check
  `tests/test_mechanisms.py` for registry invariants).
- `STATUS.md`: three new rows in the problem table (status "onboarded, no
  attempts"; budgets: billiards **high**, mahler **high**, crouzeix
  **medium** given the certification risk) and queue entries for the first
  attempts listed above, each with its kill condition inline.
- Check `tests/test_site.py` / `site/` for anything that enumerates
  problems and needs the new slugs.

**Phase 4 — verify + push.** `python -m pytest tests/ -q` must pass; run
`scripts/blind.sh` smoke test per problem; push branch. First actual
attempts (self-tests, then censuses) are separate sessions via
`scripts/new_attempt.py` — do **not** fold attempt work into onboarding.

## Non-goals for the onboarding session

- No attempt records, no census runs, no results — onboarding only.
- No edits to existing problems' files or existing attempt records (hard
  repo rule).
- No lab findings in tier-0 files, including "we plan to try X" framing.
- Do not renumber or reorder the existing STATUS.md queue; append.
