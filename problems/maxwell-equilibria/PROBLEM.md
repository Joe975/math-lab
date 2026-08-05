# Maxwell's problem on points of equilibrium

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Setting.** Fix n distinct points x₁, …, xₙ ∈ ℝ³ carrying nonzero charges
q₁, …, qₙ. The electrostatic potential is

  V(x) = Σᵢ qᵢ / |x − xᵢ|,

harmonic away from the charges. A **point of equilibrium** is a critical point
of V — equivalently a zero of the field E = −∇V — away from the charges. By
Earnshaw's theorem (harmonicity), V has no local extrema, so every
nondegenerate equilibrium is a saddle.

**Maxwell's claim** (*A Treatise on Electricity and Magnetism*, 1873, §113
and footnote). The number of equilibrium points of n point charges, when they
are finite in number, is at most (n−1)². Maxwell gave a Morse-theoretic
sketch, not a proof. The statement splits into the questions that are actually
studied:

1. Is the number of equilibria always finite (say, for generic
   configurations)?
2. What is the sharp upper bound on the number of isolated (or nondegenerate)
   equilibria, as a function of n?
3. What is the sharp bound for small n — in particular n = 3?

## Published status

- **Finiteness is open in general.** Gabrielov, Novikov and Shapiro (*Mystery
  of point charges*, Proc. London Math. Soc. 95, 2007) — the modern reference
  — prove via Khovanskii's fewnomial theory that the number of *isolated*
  equilibria admits a bound depending only on n (not on the dimension), but
  whether the equilibrium set of a positive-charge configuration in ℝ³ is
  always finite remains unproved. Their bound for three charges is 12.
- **Three charges.** Conjectured sharp bound: (3−1)² = 4. For three charges of
  **equal magnitude**, the sharp bound 4 is a theorem (Ya-Lun Tsai,
  *Maxwell's conjecture on three point charges with equal magnitudes*,
  Physica D 309, 2015): the possible counts of isolated equilibria are 0, 2,
  3, 4, and unit charges at the vertices of an equilateral triangle attain
  exactly 4 (one central, three edge-adjacent). For three **positive** charges
  of arbitrary magnitudes, a July 2026 preprint (arXiv:2607.28785, *From 12
  to 6: sharpening the three-charge bound in Maxwell's problem*) claims an
  unconditional improvement of the Gabrielov–Novikov–Shapiro bound from 12 to
  6 nondegenerate equilibria. **Whether 4 is the true maximum for general
  three positive charges is open**: the gap is between 4 and 6.
- **The general conjecture is claimed false** (preprint). Arathoon, Ball and
  Kvalheim, *The Maxwell Conjecture is False* (arXiv:2607.27197, July 2026),
  construct five positive charges whose potential has at least 24
  nondegenerate critical points, exceeding (5−1)² = 16. The configuration:
  three unit charges at the vertices of an equilateral triangle, plus two
  charges of magnitude q_ε = (3/4)ε³ − (5/32)ε⁵ at height ±ε on the
  triangle's symmetry axis. The three edge equilibria of the triangle persist
  and the central equilibrium bifurcates into 21 equilibria. The proof is a
  rigorous perturbation/bifurcation argument valid "for all sufficiently
  small ε > 0"; **no explicit value of ε is certified**, and the supporting
  computations are computer-algebra, not certified numerics. As a preprint of
  days ago at the time of writing, it has not been peer-reviewed.
- **Lower-bound constructions.** Maxwell's count is attained for small n in
  the plane: four point charges in a plane can produce nine equilibrium
  points (Physica D / Appl. Math. Letters literature, 2022). Which counts
  between the known constructions and the upper bounds are realizable is
  open.
- **Restricted geometries are classical.** For n collinear positive charges,
  all equilibria lie on the line and there is exactly one in each of the n−1
  open segments between consecutive charges (the on-axis field is strictly
  decreasing between consecutive charges and blows up with opposite signs at
  the endpoints); the count n−1 is exact and elementary. For coplanar
  positive charges all equilibria lie in the plane of the charges (the normal
  field component is z·Σqᵢ/rᵢ³, which vanishes only at z = 0).
- **Degree/index identity** (classical, Poincaré–Hopf for the field E = −∇V
  on a large ball minus small balls around the charges). If all charges are
  positive and all equilibria are nondegenerate, the indices (sign det DE)
  of the equilibria sum to 1 − n. This is an unconditional consistency
  constraint on any claimed complete census.

**What is open, concretely.** The sharp maximum for three positive charges
(4, 5, or 6); certified explicit witnesses for the five-charge refutation
(the preprint's construction is asymptotic in ε); the true growth of the
maximum in n now that (n−1)² is claimed dead; finiteness.

## Verification contract

Any claim recorded against this problem must meet the bar in
`CONTRIBUTING.md`. The equilibrium equations are not polynomial, but for
rational charge data the system becomes semi-algebraic after adjoining the
distances rᵢ = |x − xᵢ|, so every quantity in a certificate is decidable over
ℚ; the contract exploits that:

- **A count claim is a certificate, not a number.** It states the
  configuration exactly (rational charges and positions), the search region,
  and produces: (i) a list of pairwise interior-disjoint **isolating boxes**,
  each certified to contain exactly one equilibrium by a stated
  fixed-point/interval-Newton criterion evaluated in outward-rounded rational
  interval arithmetic; (ii) **exclusion certificates** covering the entire
  remainder of the region (sign-definite field component, charge-neighborhood
  domination, or an a-priori localization lemma, each with its inequality
  instantiated in exact rational arithmetic).
- **A complete count** (a claim about *all* equilibria of a configuration,
  not just those in a stated region) additionally requires a proved a-priori
  localization: for all-positive charges, no equilibrium lies at radius ≥
  max |xᵢ| from the origin (radial-component lemma); mixed-sign
  configurations need their own recorded localization argument before any
  completeness language is used.
- **Floating point locates candidates; it never certifies.** Approximate
  Newton iterates and float scans may steer the subdivision; every counted
  equilibrium and every excluded region must carry a rigorous certificate.
- **Nondegeneracy and index claims** require a certified sign of det DE over
  the isolating box. A complete census of an all-positive configuration must
  check the index-sum identity (Σ indices = 1 − n) and report the check.
- **Unresolved regions are reported, never dropped.** A run that cannot
  certify a subregion (near-degenerate equilibria, insufficient precision)
  states the unresolved boxes; its count is then a lower bound plus an
  explicit gap, not a count.
- **Record-adjacent counts escalate.** A configuration certified with a
  count that contradicts a published bound, or that would set a record
  (e.g., ≥ 5 isolated equilibria for three positive charges), must be
  re-verified by the independent checker (different arithmetic, different
  fixed-point criterion) before it is recorded anywhere outside the attempt
  record, and reported as requiring escalation rather than accepted.
- **Censuses are `EVIDENCE`**, scoped by the sampled family, grid, region
  and precision, all of which must be stated. Perturbation arguments valid
  "for small ε" without a certified ε are `SPECULATION` at the instantiation
  step no matter how rigorous the asymptotics.

## Harness (tier 0)

- `harness/maxwell-equilibria/equilibria.py` — the reference implementation.
  Certified equilibrium counting for rational configurations: adaptive
  bisection with exclusion tests in outward-rounded rational interval
  arithmetic (dyadic endpoints, exact integer square-root enclosures),
  Krawczyk-operator isolation certificates, the localization and
  charge-neighborhood lemmas instantiated in exact arithmetic, certified
  det-sign indices, and the index-sum consistency check. Emits a JSON
  certificate with the full leaf decomposition (split paths), so the covering
  is independently re-checkable. Standard library only.
- `harness/maxwell-equilibria/verify_equilibria.py` — the independent
  checker. Re-verifies a claimed certificate with different machinery:
  directed-rounding floating-point interval arithmetic (with an exact
  rational fallback), interval Newton with interval Gaussian elimination
  instead of the Krawczyk operator, its own derivations of the localization
  and domination lemmas, and an exact combinatorial check (Kraft equality +
  prefix-freeness) that the leaf boxes tile the search region. Run it on any
  certificate you intend to cite; a count confirmed by both routes at a
  record-adjacent value is reported as requiring escalation rather than
  accepted.
