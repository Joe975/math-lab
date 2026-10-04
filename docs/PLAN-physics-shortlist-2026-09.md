# Five physics problems: selection and ideation handoff

Date: 2026-09-06. Status: **MAP — proposals only**. Mode: **informed**.
No problems have been onboarded, no approaches executed, and no agents dispatched.
This is a single-author preparatory sweep, adapting docs/IDEATE.md for candidate
problems that have no local attempt index. It is not an independently reviewed
attempt record. All prospects of improvement below are **SPECULATION**.
The approach families are established methods; novelty belongs to a future
specific result, not to renaming a method here.

## Project review and selection

Reviewed AGENTS.md, README.md, CONTRIBUTING.md, docs/IDEATE.md, mechanisms.json,
tiers.json, the current STATUS.md, previous physics/E&M onboarding plans, the
problem/harness inventory, and the almost-mathieu and three-phase-conductivity
indexes. This is a research-workflow review, not a correctness audit of all proofs.

- There are 16 problem directories; README's ten-problem overview is stale.
  STATUS and per-problem indexes are better starting points for assignments.
- Existing strengths: exact finite witnesses, independent certificate checking,
  symmetry reduction, algebraic arithmetic, and explicit records of failed routes.
- Existing lessons to carry over: inspect parameter boundaries, distinguish local
  from global optimization, verify coverage geometrically, and check the literature
  before treating a recovered construction as new.
- Existing physics coverage includes triangular billiards, Mahler, Crouzeix,
  Maxwell equilibria, almost Mathieu, and three-phase conductivity. Thomson was
  already proposed and deferred, so it is not counted as a new suggestion.
- The shortlist favors mathematical physics with checkable milestones. It is
  deliberately quantum-heavy; it does not promise discoveries of new physical laws.
  Unlike the earlier rejection of general n-body finiteness, the celestial-mechanics
  item below targets one fixed eight-body system and bounded certificate tasks.
- Literature checks found current primary research, including 2026 work. They are
  selection-level checks, not an exhaustive priority search. Pin full-text versions,
  supplements, exact benchmarks, and hypotheses before a research run begins.

| Priority | Proposed slug | Open target | Target shape | First useful deliverable |
|---|---|---|---|---|
| 1 | bell-grothendieck-3 | Sharpen the order-three quantum/classical correlation ratio | object / proof | Small exact Bell witness and independently checked classical bound |
| 2 | quantum-backflow | Sharpen rigorous bounds on maximum single-interval backflow | mixed | Certified trial-state bound, including integration error |
| 3 | mub-6 | Determine whether four mutually unbiased bases exist in C^6 | object / proof | Exact checker and a precisely scoped exclusion or construction |
| 4 | central-configurations-8 | Complete classification of eight equal planar gravitating masses | mixed | Certified configurations and explicit uncovered regions |
| 5 | sic-povm | Establish SIC existence in every finite dimension | object / proof | Exact certificate for a catalogue gap or a restricted algebraic lemma |

These priorities measure fit and likely useful output, not probability of solving
the full problem. A complete solution of any headline remains a long shot.

## 1. Bell nonlocality: Grothendieck constant of order three

**Physics.** How much can spin-singlet correlations outperform classical local
correlations? For a nonzero real rectangular matrix M define

    C(M) = max_{a_i,b_j in {-1,1}} |sum_ij M_ij a_i b_j|,
    Q3(M) = max_{u_i,v_j in S^2} |sum_ij M_ij (u_i dot v_j)|.

The open constant is K_G(3) = sup_M Q3(M)/C(M). A single certified matrix and set
of vectors can improve a lower bound; solving the quantum maximum is unnecessary
for that purpose. The unrestricted real K_G is a different target.
Designolle–Vértesi–Pokutta's [February 2026 paper](https://doi.org/10.1103/m5mn-c5dx)
improves finite-order bounds with rectangular instances and identifies exact
classical optimization as the certification bottleneck. Use its certified d=3
bound, checked against successors, as the launch benchmark; do not use an old
1.4367 figure as the current record.

**Checker first.** Integer/rational M, exact unit vectors (rational
stereographic parametrization suffices), and a certificate for C(M). For small
row count enumerate a, evaluating sum_j |sum_i M_ij a_i| exactly. Validate on
CHSH, whose ratio is sqrt(2), and deliberately corrupted witnesses.

| Approach / lens | First experiment | Success and route rejection |
|---|---|---|
| A. Symmetric spherical configurations / algebra | Form measurement sets from small polyhedral group orbits; optimize integer coefficients constant on pair orbits. Start with at most 16 settings on the smaller side so classical enumeration is affordable. | SPECULATION: symmetry compresses a competitive witness. Reject a fixed orbit family if a certified bound on its best ratio misses the benchmark; search failure alone is not exclusion. |
| B. Cutting planes on the local correlation polytope / convex geometry | Start from 8–16 unit vectors per party. Separate their correlation matrix from deterministic sign vertices, alternate vector improvement and separating-functional optimization, then rationalize M. | SPECULATION: asymmetric sparse instances certify more cheaply. Reject candidate gains that disappear under an exact classical upper bound; a heuristic classical maximum underestimates C and inflates the ratio. |
| C. Randomized rounding / probability | Seek a distribution of threshold rounding rules with E[a(u)b(v)] = c(u dot v), or an explicitly controlled correction permitting the same inequality after summation. Reproduce a known rounding inequality before varying thresholds. | SPECULATION: an order-three-specific rule improves an upper bound. Kill a proposed rule at a certified violating angle; pointwise approximate proportionality alone is insufficient for arbitrary signed M. |

Transfer candidates: branch-and-bound, assumed-symmetry-construction,
exact-rational-arithmetic. These are mechanism transfers, not an existing Bell
checker. A smaller witness matching a known bound is an engineering result only.

## 2. Quantum backflow

**Physics.** A free particle with momentum supported on p > 0 can temporarily
increase its probability of lying on x < 0. Define c_BM as the supremum of that
increase over normalized states and one time interval. Improve rigorous lower
or upper bounds; do not assume the supremum is an attained eigenvalue.

[Penz et al.](https://arxiv.org/abs/quant-ph/0511109) establish that the relevant
operator is bounded and self-adjoint but not compact. [Fewster–Kirk-Karakaya](https://arxiv.org/abs/2505.13184)
give an updated numerical estimate around 0.0384506, rather than the older
0.0384517. The [2026 thesis](https://etheses.whiterose.ac.uk/id/eprint/38105/)
reports a new rigorous lower bound. The numerical estimate is not a certified
interval. Read that bound before selecting an improvement target. Keep the
single-interval, positive-momentum problem separate from generalized backflow.

**Checker first.** An explicit positive-momentum wavefunction, certified norm,
and rigorously enclosed probability difference. Recompute using an independently
derived current integral or position-probability formulation, with justified
interchanges of integrals and explicit quadrature/tail errors. A truncated
matrix eigenvalue is not an upper bound on the full operator.

| Approach / lens | First experiment | Success and route rejection |
|---|---|---|
| A. Explicit trial functions / analysis | Optimize 2, 4, 8, then 16 compact-support piecewise-polynomial momentum functions; round coefficients to rationals and rigorously integrate the backflow kernel. | SPECULATION: a compact certificate approaches or improves the launch lower bound. Reject an apparent gain if normalization or integration error swallows it. A certified weaker bound remains a benchmark, not new progress. |
| B. Full-operator domination / operator theory | Split momentum into a finite interval and its complement; derive bounds for the complement and cross blocks, and combine these with a certified finite-block bound. | SPECULATION: oscillatory structure makes the resulting upper bound useful. First test whether the tail estimate improves with cutoff; noncompactness means norm convergence cannot be presumed. Stop at the precise non-decaying term if it defeats the method. |
| C. Wavepacket interference design / dynamics | Construct separated momentum packets with independently varied phases and widths; derive a tractable flux formula and certify the best small-packet state. Compare its cost/quality against A at equal parameter count. | SPECULATION: interference structure offers a more efficient trial state. Reject claimed positive-momentum examples with leakage into p < 0, and reject the finite packet ansatz as competitive if a certified ansatz bound is too low. |

Transfer: certified-enclosure and the lab's independent-checker discipline.
New work is required for validated oscillatory integration; this is not already
provided by harness/common.

## 3. Four mutually unbiased bases in dimension six

**Physics.** Measurements are mutually unbiased when an eigenstate of one gives
uniform outcome probabilities in another. Seek four 6-by-6 unitary matrices U_a
such that |(U_a^* U_b)_ij|^2 = 1/6 for a != b. Three bases are known; the
dimension-six problem remains open. The narrower four-basis target gives a
finite witness that would itself be significant, without demanding seven.

[Mortimer's 2025 primary paper](https://doi.org/10.1103/m8qd-x1js) studies polynomial
optimization, symmetry reduction, and smaller necessary constellations, explicitly
warning that the complete computation remains intractable. The [2026 literature
review](https://doi.org/10.22331/q-2026-04-01-2051) is an additional launch reading
list, not a substitute for checking the original exclusion theorems.

**Checker first.** Exact algebraic entries and exact orthogonality/unbiasedness,
or a rigorous existence certificate for the full real polynomial system. Fix
U_0 = I by unitary equivalence; distinguish this valid gauge choice from imposing
a special Hadamard family, which is a restriction. Validate on known triples and
complete prime-power examples. Residuals near zero certify nothing by themselves.

| Approach / lens | First experiment | Success and route rejection |
|---|---|---|
| A. Hadamard-family extension / algebra | Pick a published six-dimensional Hadamard family, audit existing extension exclusions, and reduce the equations for two additional bases on one uncovered parameter region. | SPECULATION: a not-yet-excluded family has tractable elimination. Exact inconsistency kills that family/region only; if literature already excludes it, discard before computing. |
| B. Finite phase alphabets / combinatorics | Use roots of unity of orders 6 and 12 as checker controls, then enumerate a predeclared unspent alphabet. Encode orthogonality and unbiasedness in its cyclotomic field, with symmetry-breaking and exact verification. | SPECULATION: a discrete structured witness exists. Exhaustive failure excludes that alphabet only; arbitrary complex solutions need not use roots of unity. Broaden to a structurally different route after exhaustion. |
| C. Necessary constellations / polynomial optimization | Extract one smaller necessary vector constellation from the 2025 formulation; reproduce a solved lower-dimensional certificate, then attempt rational sum-of-squares/branch certificates for the six-dimensional case. | SPECULATION: the smaller system admits a practical infeasibility certificate. A feasible constellation or zero relaxation bound defeats this relaxation, not the MUB problem. Every four-basis solution must be shown to contain the selected constellation. |

Transfer: bitmask/backtracking ideas, exact algebraic-field arithmetic from
almost-mathieu (inspect compatibility first), and certified-box-coverage.

## 4. Eight equal masses: planar central configurations

**Physics.** These shapes rotate rigidly or undergo homothetic gravitational
motion. For eight unit masses, solve

    sum_{j != i} (q_j - q_i)/|q_j-q_i|^3 = -lambda q_i,
    sum_i q_i = 0,    sum_i |q_i|^2 = 1,

with distinct q_i in R^2, modulo rotations, reflections, and permutations.
The ambitious target is a complete classification; the first target is certified
coverage of a stated compact, collision-free part of the normalized domain.

[Moczurad–Zgliczynski](https://doi.org/10.1007/s10569-019-9920-6) already classify
the equal-mass planar cases n=5,6,7 and certify asymmetric examples for n=8,9,10.
Thus finding any asymmetric eight-body example is not new. The literature search
did not locate a complete n=8 classification; confirm this against the authors'
current supplements before launch. [January 2026 work](https://arxiv.org/abs/2601.01165)
also treats selected exceptional unequal five-mass cases, which are not substitutes
for this fixed equal-mass eight-body task.

**Checker first.** Remove continuous symmetry with justified coordinate charts;
certify roots by interval methods in a nonsingular reduced system and verify all
original equations. Compare configurations through distance data plus consistent
permutations, not a sorted distance list alone. Count completeness requires a
domain cover and collision exclusion, not just isolated roots.

| Approach / lens | First experiment | Success and route rejection |
|---|---|---|
| A. Symmetry-class enumeration / algebra | Enumerate partitions into dihedral orbits for eight points, derive radial/angle equations, and certify roots in one bounded parameter box per class. Compare with published configurations. | SPECULATION: an unclassified symmetry class admits a complete small census. Failure outside these classes is invisible; reject any global-classification claim based on symmetric ansatz completeness. |
| B. Symmetry-breaking continuation / dynamics | Vary one mass in a known eight-body branch, monitor the reduced Jacobian, follow bifurcating branches, and return to the equal-mass slice. Independently isolate endpoints. | SPECULATION: asymmetric branches escape an orbit search. A returned endpoint equivalent to a known configuration kills its novelty claim; failure to return to equal masses invalidates it for this target. |
| C. Collision compactification and exclusion / geometry-analysis | Use moment-of-inertia normalization for an outer bound. First try to prove an explicit minimum-separation bound; otherwise certify only r_ij >= delta, starting with delta=1/10, and expose every uncovered collision stratum. | SPECULATION: cluster force estimates close the excluded region and enable a finite cover. If they do not, stop at the unresolved cluster inequality and keep the count explicitly partial. |

Transfer: maxwell's root-isolation/coverage patterns, canonicalization, and
assumed-symmetry-construction. The gravitational equations and gauge need a new
checker; electrostatic uniqueness certificates cannot simply be relabeled.

## 5. Symmetric informationally complete quantum measurements

**Physics.** SICs are highly symmetric measurements for quantum-state tomography.
Seek d^2 unit vectors psi_i in C^d with

    |<psi_i,psi_j>|^2 = 1/(d+1),  i != j.

Their projectors divided by d form the POVM. Existence for all finite d is the
headline; exact certification of a genuinely unfilled dimension/orbit or a
proved construction lemma is the achievable first target.

[Bengtsson–Grassl–McConnell](https://arxiv.org/abs/2403.02872) give constructions
in the family d=n^2+3=4p, but not a general proof of their recipe.
[June 2026 work](https://arxiv.org/abs/2606.23535) refines the proposed relation
between overlaps and Stark units, including non-minimal cases. Therefore do not
pick an arbitrary moderately large dimension and call it unexplored. First build
a cited ledger distinguishing numerical existence, exact existence, and orbit
classification. No particular unfilled dimension is asserted by this proposal.

**Checker first.** Exact algebraic-field arithmetic and independent Gram-matrix
verification. For Weyl–Heisenberg constructions verify normalization and every
nonidentity displacement overlap, then independently check the generated frame.
Use d=2 and d=3 as controls; handle the continuous family in d=3 rather than
assuming every fiducial is an isolated polynomial root.

| Approach / lens | First experiment | Success and route rejection |
|---|---|---|
| A. Fiducials in symmetry eigenspaces / representation theory | Reproduce a small known fiducial in an order-three Clifford eigenspace; benchmark exact solving before selecting the smallest documented catalogue gap compatible with that symmetry. | SPECULATION: symmetry makes the next gap affordable. Failure in the eigenspace does not refute general SIC existence or arbitrary covariance. |
| B. Algebraic reconstruction / number theory | Take a published numerical-only candidate with reproducible provenance; increase precision, propose minimal polynomials and a common field, then verify all overlaps exactly. | SPECULATION: field structure keeps certificate size manageable. Reject polynomial fits that fail at new precision or exact substitution; successful reconstruction may still be a rediscovery after catalogue comparison. |
| C. Gram-matrix completion / geometry | Search Hermitian d^2-by-d^2 matrices with diagonal 1, off-diagonal squared modulus 1/(d+1), and G^2=dG. Test low-dimensional block-phase patterns without imposing Weyl–Heisenberg covariance. | SPECULATION: another block structure yields an exact construction or exclusion lemma. Reject a proposed pattern by exact inconsistency; do not promote its failure to nonexistence in dimension d. Hermiticity and the polynomial identity give eigenvalues 0,d, with trace fixing rank d. |

Transfer: exact algebraic arithmetic, assumed symmetry, formal polynomial
identities. Large ray class fields may require Sage/PARI beyond the lab's usual
stdlib tooling; measure that cost before committing to this item.

## Later dispatch contract

1. Onboard selected problems with published-background-only PROBLEM.md files,
   empty indexes, site explainers, and checker tests. Keep this file tier 1.
2. Freeze a cited target/benchmark before optimization. For SIC and n=8 central
   configurations the precise catalogue gap remains a launch prerequisite.
3. Assign one route at a time with its stated search universe, first experiment,
   verifier, and failure criterion. Agents reading these routes are informed.
   A genuinely blind comparison must receive only the neutral statement/harness.
4. Build and test the checker before object search. Known examples must pass;
   malformed objects, wrong normalizations, missing coverage leaves, and incorrect
   symmetry identifications must fail where relevant.
5. Record explored ranges, unresolved regions, dependencies, seeds, commands,
   target_shape, and push_rounds. Separate an exhausted ansatz from a refuted
   conjecture. A budget stop is not a mathematical kill condition.
6. Use independent skeptical verification before any VERIFIED claim. Finite
   censuses remain EVIDENCE about their stated universe; proof steps stop at a
   precisely named obstruction. Register new mechanism vocabulary only when
   actual attempts justify it, following docs/IDEATE.md.

Recommended first pair: Bell route A/B for a cheap exact-object workflow, and
backflow route A for a physically direct certification task. MUB supplies the
clearest finite-object moonshot. Reserve broad classification and large SIC
fields until their pilot costs and catalogue checks are known.

## Validation of this planning change

The handoff is explicitly registered in tiers.json as tier 1. No executable
research logic changed. `git diff --check` passed. The tier-isolation suite
completed with 9 passes and 17 failures after using a repository-local temporary
directory and a process-local Git ownership exception. All 17 remaining failures
are blind-checkout invocations reporting `dirname`/`mkdir` unavailable in the
Git Bash environment; adding Git's usr/bin to the parent PATH did not resolve
them. The full blind-checkout validation therefore remains unverified in this
session. No machine-wide configuration was changed.

## Execution addendum: problem 4 (2026-09-06)

The user selected central-configurations-8. It is now onboarded, with the first
three-route attempt and an independent skeptical audit recorded in
`problems/central-configurations-8/prior-art.json`. Twenty known configurations
are independently certified; full classification remains open. The scoped
results, failures and reproduction commands are in attempt 001.

The subsequent full suite passed: 507 tests, 7 legacy skips. The earlier
blind-checkout environment failure was resolved for that run by exporting
PATH=/usr/bin:/bin:$PATH inside Git Bash before launching the repo's Python;
the process-local command is recorded in the attempt. No global settings changed.
