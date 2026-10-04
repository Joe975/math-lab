# 001 — Three routes toward an eight-body classification

- **Problem:** central-configurations-8
- **Date:** 2026-09-06
- **Mode:** informed
- **Type:** computational search, restricted-family proof, coverage barrier
- **Tools:** direct-force and independent polynomial certificate checkers; NumPy 2.5.2 and SciPy 1.18.1 for discovery; exact Fraction arithmetic for verification.
- **Sources:** publisher HTML and the authors' numerical data, cited below. No PDF transcription was used. This follows proposal 4 in docs/PLAN-physics-shortlist-2026-09.md.
- **Target shape:** mixed, with full classification a proof target.
- **Push rounds:** 3: fixed-mass continuation after seed search; pseudo-arclength after five failed paths; catalogue replication and global-coverage pilot after no new endpoints.

## Approach

Attack the classification through symmetry, continuation, and geometric exclusion.
First build independent local-root checkers so numerical roots cannot masquerade
as exact configurations. The key distinction is existence versus completeness:
this attempt establishes the former for known objects and exposes the latter's
remaining gap. No new central configuration is claimed.

## What was done

### Exact certificate contract

Work with eight unit masses, lambda=1, and y0=0. The 15-variable direct system
uses f_i=q_i-sum_j(q_i-q_j)/r_ij^3, omitting f0y. The independent 43-variable
system replaces inverse distances by 28 positive variables s_ij and imposes
s_ij^2*r_ij^2=1, using s_ij^3 in the force.

Every saved certificate contains an exact rational center, radius 1/100000000,
and a proposed rational inverse C. The checkers recompute the residual and
Jacobian enclosures. For T(x)=x-Cf(x), they verify strict self-map inclusion and
an infinity-norm derivative bound less than one. Banach gives a unique fixed
point. The same derivative bound implies C is nonsingular, hence f=0 there.
The torque identity sum_i q_i cross f_i=0 recovers f0y when x0 excludes zero;
summing all force equations gives center of mass zero. These are local gauge
certificates, not global uniqueness or stability certificates.

The independent verifier shares neither interval operations nor force/Jacobian
code with the direct verifier. It was written and attacked by a separate skeptic
agent; see 002 and data/skeptic-notes.md.

### A. Symmetry and discovery

With RNG seed 20260906, 24 structured starts and 400 random starts produced 19
certified pairwise-inequivalent classes. Structured starts included the octagon,
heptagon plus center, a collinear seed, aligned/twisted square pairs, and other
ring partitions. Starts are not an exhaustive enumeration of dihedral orbits.
The discovery search did not use the published coordinates as seeds.

Numerical isometry matching uses a consistent vertex permutation, not just a
sorted distance list. Its correctness test includes two noncongruent homometric
point sets. All 19 endpoints numerically match entries in the published catalogue.
The catalogue's twentieth entry was not recovered by this seed design.

Separately normalized, polished, and certified all 20 published candidates from
Doicu–Zhao–Doicu. Both checkers accept all 20. Exact squared-distance order
statistics separate every one of the 190 unordered pairs. The independent audit
re-establishes those separations from the polynomial boxes. This is a rigorous
lower bound of 20 similarity classes including the collinear class, and a
replication of known configurations; it is not evidence of novelty.

Source coordinate file SHA-256:
16f64ce75476006ecad654a5f94a1992770551cdfee1f796708cde0a93dc0b23.

### A. Complete aligned-square subfamily

Put the inner square at (+/-1,0),(0,+/-1), and the outer at (+/-t,0),(0,+/-t),
t>1. Let c4=(1+2sqrt(2))/4. Tangential forces cancel. Subtraction of the two
radial coefficients gives

    D(t)=c4(1-t^-3)-(1+1/t)/(t-1)^2+(1-1/t)/(t+1)^2.

The perpendicular cross terms cancel. D=0 is equivalent to

    c4=H(t)=2t^2(3t^2+1)/[(t^3-1)(t^2-1)^2].

For numerator N and denominator P,

    N'/N=2/t+6t/(3t^2+1)<4/t,
    P'/P=3t^2/(t^3-1)+4t/(t^2-1)>7/t.

Therefore H is strictly decreasing from infinity to zero on (1,infinity), so
there is exactly one ratio. The outer radial coefficient is positive, allowing
rescaling to lambda=1. Exchanging radii handles t<1; t=1 collides and is excluded.
This completes precisely the aligned, distinct concentric-square family modulo
similarity and permutation. It does not cover all D4-invariant shapes.

Independent skeptical re-derivation found no gap. Exact rational endpoint tests
in data/aligned-square-bracket.json enclose t near 2.24927674247929 in a bracket
of width 2^-70. The scalar reduction is also tested against the full force sum.
This is a restricted rediscovery of the nested-polygon result, not a new theorem;
Moeckel–Simo (1995) is the prior-art reference.

### B. Continuation, then a second push through folds

Vary the anchor mass toward 1/2 or 2 and back, keeping its label and the gauge
fixed. Among 38 attempted paths, 33 completed 20 prescribed steps and 5 failed.
The 99 perturbed equal-mass endpoint searches added no class. Failure of a
corrector is not evidence of a branch ending or a bifurcation.

Pseudo-arclength continuation then followed both directions from all 19 seed
classes, using an SVD tangent and a predictor/corrector step. There were 38 paths,
3984 accepted steps and 11 sampled crossings of the equal-mass slice. All returned
endpoints numerically matched existing classes. Twenty-six paths stopped at the
predeclared mass/gauge limits and 12 at the 140-step budget. These are EVIDENCE
about sampled curves; neither path continuity nor absence of other crossings
was rigorously certified. No branch-completeness claim is made.

### C. Published compactness bounds and coverage barrier

The proposed need to invent collision bounds was retired by the literature
check. In the unit-mass lambda=1 normalization, the published n-1 radius bound
rescales to R<=14. From I=U<=8*14^2, every pair has r_ij>1/1568. The equal-mass
radius lower bound gives R^3>=7/4, so a farthest-body gauge has x0>1. Every class
therefore has a representative in [1,14] x [-14,14]^14 with y0=0, subject to the
farthest-anchor constraint. These bounds are published consequences, not new.

The exact coverage pilot bisected the longest coordinate to depth 12, retaining
actual axes and rational cuts. It attempted exclusion using center of mass,
farthest-anchor, and minimum-separation constraints. All 4096 leaves remain
UNRESOLVED: zero boxes were excluded. The certificate is a complete partition
of the rectangular outer domain with honest unresolved labels, not a root cover.
No global volume or completeness claim can be extracted from it.

The concrete proof obstruction is excluding the complement of the local root
neighborhoods. Broad boxes straddle collision and cancellation regions, so these
necessary geometric predicates have no purchase at this resolution. SPECULATION:
stronger force/cluster predicates combined with permutation reduction could
improve this; the pilot supplies no quantitative evidence that they will scale.

### Reproduction

All paths below are absolute. Use the repo virtual environment; discovery
requirements are pinned in explore/requirements.txt. To preserve this record,
write reproduction output beneath C:/Repos/math-lab/out/cc8-reproduction and
substitute those output paths for input dependencies of later commands.

```powershell
& C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/search.py --output C:/Repos/math-lab/out/cc8-reproduction/run-001 --random-starts 400 --seed 20260906
& C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/continuation.py --source C:/Repos/math-lab/out/cc8-reproduction/run-001 --output C:/Repos/math-lab/out/cc8-reproduction/arclength-001 --steps 140
& C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/catalogue.py --source C:/Repos/math-lab/problems/central-configurations-8/data/published-CCNMAS8.dat --discovery C:/Repos/math-lab/out/cc8-reproduction/run-001 --output C:/Repos/math-lab/out/cc8-reproduction/catalogue-001
& C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/symmetric.py --output C:/Repos/math-lab/out/cc8-reproduction/aligned-square-bracket.json
& C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/coverage.py --depth 12 --output C:/Repos/math-lab/out/cc8-reproduction/coverage-depth12.json
& C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/audit_certificates.py
```

The last command rechecks the committed evidence and refreshes its verification
summary. All research outputs used here are under data/. Exact certificate
acceptance is deterministic; floating discovery trajectories can depend on
platform/BLAS versions even with a fixed seed. No total-runtime benchmark is
claimed; execution was interleaved with verification and writing.

### Repository validation

Final full suite: **507 passed, 7 skipped**, including blind-checkout isolation,
record/index consistency, site generation and all new mathematical tests.
`git diff --check` passed. The seven skips are pre-existing legacy-record
reference exemptions. NumPy/SciPy were installed only in the repo virtual
environment. To make Git Bash's utilities available in this Windows environment,
the full run was launched with its PATH explicitly set inside Bash:

```powershell
& 'C:/Program Files/Git/bin/bash.exe' -c 'export PATH=/usr/bin:/bin:$PATH; export GIT_CONFIG_COUNT=1 GIT_CONFIG_KEY_0=safe.directory GIT_CONFIG_VALUE_0=C:/Repos/math-lab; /c/Repos/math-lab/.venv/Scripts/python.exe -m pytest C:/Repos/math-lab/tests/ -q -p no:cacheprovider --basetemp=C:/Repos/math-lab/out/pytest-cc8-final-20260906 --tb=short'
```

The environment changes are process-local. No machine-wide Git or PATH settings
were changed. Choose a fresh basetemp path for a future run.

## Outcome

**VERIFIED:** 20 known, pairwise-inequivalent eight-body central configurations,
locally isolated by independent exact certificates; the complete restricted
aligned-square classification, independently re-derived and identified as a
rediscovery. **EVIDENCE:** 19 recovered classes from 424 starting points, finite
mass-continuation and pseudo-arclength probes. **MAP:** the unresolved global
coverage barrier. Full classification is not achieved; no new shape, stability
result, or new general theorem is claimed.

## Why it failed / what survived

The headline proof needs exhaustive exclusion, which this pilot does not supply.
The independent seed search misses even one known catalogue entry, exposing a
concrete failure of numerical completeness. The fold-following push does not
remove that failure. Pure coarse geometric subdivision excludes no boxes at
the measured depth. These are limits of the tested methods and ranges, not
refutations of stronger versions or proof that the catalogue is complete.

Reusable outputs are two exact local checkers, an independently checked catalogue,
consistent-permutation numerical matching, exact inequivalence witnesses, a small
complete family proof, continuation traces, and an explicit unresolved partition.
The target shape is mixed because local objects are checkable but global exclusion
is a proof step. This attempt stops at that named obstruction; budget-ended
continuation paths are not labeled dead routes.

## Leads generated

1. Start from the certified twentieth catalogue shape and measure its basin under
   controlled perturbations; compare against the same-radius perturbations of
   the other 19. Outcome: whether a narrow basin explains this search's omission.
2. On a fixed set of depth-12 unresolved boxes, add exact pair-force/cluster
   exclusion and compare the number certified empty. Zero additional exclusions
   rejects that predicate at this resolution; positive counts justify an adaptive
   pilot, not a complete n=8 run.
3. Prove or refute a coverage lemma that combines permutation ordering and
   farthest-body charts without losing asymmetric solutions. Test it against all
   20 certificates before using it for any global count.
4. SPECULATION: the twisted two-square family offers a second complete tractable
   slice. Compare the exact target against existing twisted-polygon classifications
   before doing algebra; a rediscovery is acceptable only when labeled.

## References

- Moczurad–Zgliczynski (2019), https://doi.org/10.1007/s10569-019-9920-6 .
- Doicu–Zhao–Doicu (2022), https://doi.org/10.1007/s10569-022-10075-7 .
- Coordinates: https://github.com/AlexandruDoicu/Balanced-and-Central-Configurations/blob/main/CentralConfig/CCNMAS8.dat . Retrieved 2026-09-06; exact file hash above. Its original mass and inertia normalization is explicitly replaced before certification.
- Moeckel–Simo (1995), https://doi.org/10.1137/S0036141093248414 . Nested-polygon prior art, also described in Corbera–Valls (2021), https://doi.org/10.1007/s00332-021-09743-z .
- Independent review: attempts/002-independent-certificate-audit.md and data/skeptic-notes.md.
