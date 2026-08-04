# Plan: onboard electricity-and-magnetism problems

Status: **plan only — nothing below is implemented.** Written 2026-08-04 as a
handoff for an onboarding session. Branch:
`claude/electrical-magnetism-exploration-o4nztc`. Companion to
`docs/PLAN-physics-problems.md` (2026-07-28), which onboarded the first three
physics-flavored problems; the same phase structure applies and is not
repeated in full here.

This file is **not tier 0** (see `tiers.json`): it contains this lab's framing
and route ideas, so it must never be copied into a blind checkout or cited in
`PROBLEM.md` prose.

## Scoping decision

The request: open mathematics around electric fields and magnetism, attacked
in the hope of something physics-relevant. Same selection filter as the July
round — the open frontier must be attackable with this lab's methodology:
exact or certified arithmetic in stdlib Python (+ optional C kernels),
census-shaped falsifiable first steps, honest `EVIDENCE` framing.

**Expectation-setting, stated up front so the ledger stays honest.** None of
these problems will produce "new physics" in the discovery sense; they are
open *mathematics* whose statements are physical. The realistic deliverables
are the usual ones — approach library, dead ends with reasons, `EVIDENCE`
over stated ranges — plus one genuine long-shot upside: several of the
candidates below are refutable by a single certified finite object (a charge
configuration, a point configuration), and a refutation of a 150-year-old
electrostatics conjecture would be a real result with physical content
(trapping geometries, Earnshaw-adjacent structure). That is the honest
version of "physics-impacting": possible, census-shaped, and not the
expected outcome.

## Recommended: one problem now, one second

### 1. `maxwell-equilibria` — Maxwell's conjecture on points of equilibrium

- **Conjecture** (J.C. Maxwell, *Treatise on Electricity and Magnetism*,
  1873, §113 footnote). The electrostatic field of N fixed point charges in
  ℝ³ in general position has at most (N−1)² isolated equilibrium points
  (zeros of the field). **Open even for N = 3**, where the conjectured
  maximum is 4. Best published bound for three charges is 12, via Khovanskii
  fewnomial theory (Gabrielov–Novikov–Shapiro, "Mystery of point charges",
  Proc. LMS 2007). Sharper still: whether the equilibrium set is always
  *finite* is itself open in general. **Verify all of this precisely during
  PROBLEM.md drafting** — GNS 2007 is the anchor citation; search for
  post-2007 improvements (planar/collinear special cases, two-charge and
  symmetric-configuration results) and cite the exact current frontier.
- **Physics lens.** This *is* electrostatics: equilibria of Coulomb fields,
  Earnshaw's theorem (every such equilibrium is a saddle — no stable trap),
  Morse theory of the potential. Ion-trap adjacent. Fills the E&M gap in the
  lab's physics coverage with the most literal candidate available.
- **Key structure.** Equilibria are critical points of V(x) = Σ qᵢ/|x−xᵢ|.
  The system is not polynomial, but adjoining rᵢ = |x−xᵢ| as variables with
  rᵢ² = |x−xᵢ|², rᵢ > 0 makes it semi-algebraic — so exact-arithmetic and
  certified-enclosure machinery applies. For rational configurations
  everything lives over ℚ. A complete certified count over a bounded region
  is the shape of the work: enclose each equilibrium (interval Newton /
  Krawczyk), certify the rest of the region equilibrium-free (exclusion
  boxes), and derive per-configuration a-priori bounds that confine all
  equilibria to a computable ball (for all-positive charges they lie in the
  convex hull of the charges; the mixed-sign bound needs its own small lemma
  — a harness design item, not an afterthought). The affine-forms lesson
  from the billiards harness (STATUS insights, 2026-07-28) applies verbatim:
  iterated/recombined interval geometry must start from first-order forms.
- **Harness** (`harness/maxwell-equilibria/`):
  - *equilibria.py* — reference implementation. Input: rational charges and
    positions. Output: a certified equilibrium count — disjoint isolating
    boxes each certified to contain exactly one nondegenerate zero, plus a
    covering family of exclusion boxes for the remainder of the a-priori
    ball, plus the certified a-priori ball derivation itself. Stdlib only.
  - *verify_equilibria.py* — independent checker by a *different* method:
    re-certify isolating boxes via topological degree (sign analysis on box
    faces) rather than Krawczyk, and re-derive the count for coplanar
    configurations by an independent reduction. No shared code with
    equilibria.py beyond `fractions.Fraction`.
  - Verification contract for PROBLEM.md: a count claim states the
    configuration exactly (rational data), the a-priori ball and its proof,
    the isolation and exclusion certificates, and the arithmetic used.
    Floating point locates candidates; it never certifies. Degenerate
    (non-isolated or near-degenerate) configurations are reported as
    UNRESOLVED at stated tolerance, never silently skipped.
- **First attempts (queue lines):**
  1. Self-test: three collinear charges (the classical, fully analyzable
     case) and symmetric triangle configurations — reproduce known counts
     with certificates end to end. Expected `VERIFIED`, scope = those
     families. Run blind by design (nearly free before prior art exists).
  2. Census: 3-charge configurations, quotiented by similarity (2 shape
     parameters + 2 charge ratios, rational grid). Distribution of
     certified equilibrium counts; map the strata where the count changes.
     **Win condition, kept in view by every run: any configuration
     certified with ≥ 5 isolated equilibria refutes Maxwell's N=3 count
     outright.** `EVIDENCE` scoped by grid and tolerance.
  3. Later: N = 4 (conjectured max 9), guided by whatever the N=3 strata
     map shows about where counts jump.
  - Kill conditions to record at queue time: (i) if certified counting near
    the degenerate strata (where equilibria merge) forces resolution costs
    that grow without bound, the full-census route dies — record the
    reachable margin instead; (ii) if the mixed-sign a-priori ball lemma
    resists elementary proof, scope the census to all-positive charges and
    say so.
- **Budget: high.** Best physics-lens fit, decidable-over-ℚ certificates,
  century-old open statement with a one-object refutation path.

### 2. `thomson-sphere` — minimal Coulomb energy on the sphere

- **Conjecture / problem** (J.J. Thomson, 1904). N unit charges on S²
  minimizing Coulomb energy Σ 1/|xᵢ−xⱼ|. Global minimizers are **proven**
  only for N = 2, 3, 4, 6, 12 (linear-programming / universal-optimality
  methods) and N = 5 (R. Schwartz 2013, computer-assisted interval
  arithmetic). **N = 7 is open** — conjectured minimizer the pentagonal
  bipyramid. Related but distinct: Smale's 7th problem (logarithmic energy,
  algorithmic formulation) — keep the two separated in PROBLEM.md. **Verify
  the exact proven-N list and the N=7 status during drafting** (the N=5
  proof is the methodological anchor: it shows a stdlib-style
  interval-arithmetic optimality proof is a real, if heavy, object).
- **Physics lens.** Electron shells on a sphere; realized physically in
  multi-electron bubbles in liquid helium and (as a family of energies)
  colloidosomes/capsid structure. Electrostatics of confined charge.
- **Key structure.** Finite-dimensional certified global optimization.
  Energies are algebraic; comparisons between configurations reduce to
  certified enclosures. Symmetry reduction cuts the domain; the honest
  deliverable is a *landscape census*: all local minima for N = 7 (and
  neighbors) under a stated search design, certified energy enclosures, and
  the gap between the best two — `EVIDENCE` about the design, plus the
  building blocks a Schwartz-style N=7 proof would need.
- **Harness sketch** (`harness/thomson-sphere/`): *energy.py* (certified
  energy enclosures for exact/near-exact configurations; certified local
  minimality via interval Hessian) and *verify_energy.py* (independent
  re-computation, different parametrization). Contract: every reported
  energy is an interval; "global minimum" is never claimed, only "least
  found under design D" plus certified local minimality.
- **First attempts:** self-test on N ≤ 6, 12 knowns; then the N=7 landscape
  census. Kill condition: if enclosure widths near the putative minimizer
  cannot be driven below the energy gap to the second-best basin at
  affordable cost, the certification route dies at N=7 — measure and record
  the cost curve first, as with Crouzeix.
- **Budget: medium.** Heavier certification load than maxwell-equilibria;
  the full-optimality moonshot is real but expensive.

## Watch list (needs Phase-0 literature verification before any commitment)

- **Sendov's conjecture / Smale's mean value conjecture.** Both are exactly
  2D electrostatics of unit charges (critical points of the logarithmic
  potential of polynomial roots). By far the best exact-arithmetic fit on
  this page — pure polynomial algebra, censuses over rational-coefficient
  families are trivial to make exact. Sendov: proved for degree ≤ 8 and for
  sufficiently large degree (Tao 2020); open in between. Smale MVC: factor 4
  known since 1981, conjectured (n−1)/n. Reason not recommended now: the
  physics lens is the thinnest here (the electrostatic reading is an
  interpretation, not the subject), and the request was E&M. If a later
  round wants a near-guaranteed-tooling win, start here.
- **Almost Mathieu / Hofstadter spectral questions ("dry ten martini",
  critical coupling).** The one candidate whose physics is genuinely
  *magnetic* (Bloch electrons in a magnetic field). Transfer-matrix spectra
  admit certified enclosures in stdlib. But the open/solved boundary has
  moved fast (Ten Martini: Avila–Jitomirskaya 2009; dry-version and
  critical-case claims in 2023–2025 preprints) and this plan's author could
  not pin the current frontier from memory with tier-0-grade confidence.
  Do not onboard without a careful literature pass; if the dry problem
  stands in a census-able form, it becomes a strong third problem.

## Candidates that lost (recorded so the next sweep doesn't re-tread)

- **Yang–Mills mass gap.** No finite falsifiable step exists at any budget;
  nothing census-shaped to do. Out.
- **Pólya–Szegő plate-capacity conjecture** (the disk minimizes capacity
  among convex plates of given area, 1951 — open). Certified capacity needs
  certified PDE/integral-equation numerics in stdlib Python; loses for
  exactly the reason hot-spots lost in the July round.
- **Magnetic relaxation / helicity energy bounds (Arnold; Freedman–He
  asymptotic crossing number).** Frontier is geometric analysis; anything
  we could compute reproduces knowns. Out.
- **Anisotropic Calderón problem (EIT uniqueness).** Purely analytic
  frontier; no computational census touches it. Out.
- **2D crystallization / Abrikosov lattice / universal optimality in d=2**
  (triangular lattice as Coulomb-gas minimizer; open — solved in d = 8, 24
  by Fourier-interpolation machinery). The missing piece is sharp analytic
  machinery, not evidence; a census would reproduce what everyone already
  believes. Out, with regret — it is the best "magnetism" story of the
  losers (superconductor vortex lattices).

## Implementation phases

Same as `docs/PLAN-physics-problems.md` phases 0–4, with these deltas:

- **Phase 0 (literature check, no repo changes).** The load-bearing items:
  exact statement and best bounds for Maxwell's conjecture (GNS 2007 +
  successors; the finiteness question's status); the proven-N list and N=7
  status for Thomson; the two watch-list frontiers if pursued. Nothing
  enters PROBLEM.md that is not literature.
- **Phase 1–2 (scaffolding + harness).** Slugs `maxwell-equilibria` (now)
  and `thomson-sphere` (second, or deferred if the first round's harness
  work uncovers shared certified-optimization infrastructure worth
  extracting into `harness/common/` first — decide there, record the
  decision). Existing tiers.json globs cover new slugs automatically.
- **Phase 3 (registry + ledger).** No new field lens needed: `analysis`
  (potential theory), `geometry-topology`, and `computation` cover both
  problems; mechanism tags arrive with attempts, not preemptively. STATUS
  rows: maxwell-equilibria budget high, thomson-sphere medium. Queue
  entries as listed above, kill conditions inline, appended — never
  renumber the existing queue.
- **Phase 4.** `python -m pytest tests/ -q`, `scripts/blind.sh` smoke test
  per new slug, push.

## Execution addendum (2026-08-04, same branch)

Phases 0–4 for `maxwell-equilibria` were executed the day after this plan
was written, and Phase 0 materially changed the framing: **the general
Maxwell conjecture was refuted between planning and execution** —
arXiv:2607.27197 (July 29, 2026) constructs five charges with ≥ 24
nondegenerate equilibria (> 16), and arXiv:2607.28785 sharpens the
three-positive-charge bound to 6 unconditionally. The problem was onboarded
*reframed* around the surviving open questions (n = 3 max between 4 and 6;
explicit certified witness for the refutation, which is asymptotic-only;
growth of the max; finiteness) rather than the dead statement. See
PROBLEM.md and STATUS queue items 15–17.

Decision recorded per the Phase 1–2 delta above: `thomson-sphere` is
**deferred** to a follow-up onboarding session — the maxwell harness came
out problem-specific (its certificates lean on the charge-field structure),
so no shared certified-optimization infrastructure fell out to extract, and
the refutation shifted the value toward starting maxwell attempts sooner.

## Non-goals for the onboarding session

Identical to the July plan: onboarding only — no attempt records, no census
runs, no edits to existing problems or attempt records, no lab findings in
tier-0 files. First attempts run **blind by design** via
`scripts/new_attempt.py` in separate sessions.
