# Plan: onboard two conductance-flavored problems

Status: **implemented in the same session that wrote this** (2026-08-04,
branch `claude/eddy-currents-nanoscale-circuits-bh4fms`). Kept as the scoping
record: why these two problems, what was distilled away, and which candidates
lost. Read it before proposing a third problem in this direction.

Origin: a user question about eddy currents, nanoscale circuit geometry, and
alternative routes to superconductance. Two research directions from that
conversation were chosen for onboarding: (4) periodically driven / cavity
routes to superconductivity, and (6) metamaterial dielectric engineering of
pairing. Neither is itself a mathematical conjecture, so each was distilled to
the open mathematics at its core. That distillation is a judgment call and is
recorded here precisely so it can be challenged.

This file is **not tier 0** (see `tiers.json`): it contains this lab's framing
and route ideas, so it must never be copied into a blind checkout or cited in
`PROBLEM.md` prose.

## Distillation rationale

**Path 4 (Floquet/cavity dynamic routes).** The physics program — transient
superconducting-like states under optical driving (Mitrano et al., Nature 530,
2016, K₃C₆₀; verify citation before use in an attempt), Floquet engineering,
cavity-QED modification of pairing — has no census-attackable open conjecture
at the pairing level. Its rigorous mathematical core is the spectral theory of
field-threaded lattice models: the **critical almost Mathieu operator**
(Harper/Hofstadter model), which is simultaneously (a) the quasi-energy
operator of kicked/driven quantum systems, and (b) by de Gennes–Alexander
network theory, the operator whose spectral edge sets the critical-temperature
curve T_c(B) of a superconducting wire network (Pannetier–Chaussy–Rammal–
Villégier, PRL 53, 1984). So the "field-threaded nanoscale circuit" question
has a literal open-conjecture form. Onboarded as `almost-mathieu`.

**Path 6 (metamaterial dielectric engineering).** The physics program —
engineering the dielectric environment so the effective electron–electron
interaction favors pairing (Smolyaninov–Smolyaninova, "Is there a metamaterial
route to high temperature superconductivity?", 2014; ELC-structured Al
experiments; verify citations before use) — asks, mathematically: *which
effective response tensors can a geometric mixture of materials realize?* That
is G-closure / optimal-bounds theory for composites. The two-phase problem is
closed; the **three-phase problem in 2D conductivity is open**: Hashin–
Shtrikman bounds are not always attainable (Nesi 1995), improved bounds and
new optimal structures exist (Cherkaev 2009; Cherkaev–Zhang wheel
assemblages), and the exact attainable set in the intermediate regime is
unknown. Onboarded as `three-phase-conductivity`.

## Candidates that lost

- **Floquet heating / prethermalization bounds** (Abanin–De Roeck–Huveneers):
  genuinely open constants, but the frontier is analytic estimates with no
  exact-arithmetic census surface.
- **Rigorous BCS gap equation with engineered kernels** (Hainzl–Seiringer
  school): the existence/criticality theory is largely settled; what remains
  open is asymptotic analysis, again with no falsifiable finite computation
  for stdlib tooling.
- **Cavity-QED materials models**: no agreed-on mathematical model, so no
  conjecture to attack honestly.
- **Electron-hydrodynamics shape optimization** (superballistic constrictions):
  attractive, but certified Navier–Stokes-like PDE numerics in stdlib Python
  is a bigger project than the problem work — same reason the hot-spots
  conjecture lost in `PLAN-physics-problems.md`.

## The two problems

### 1. `almost-mathieu` — gaps of the critical Harper operator

- **Headline open problem.** Dry Ten Martini for the critical almost Mathieu
  operator: for λ = 1 and every irrational α, are *all* gaps predicted by the
  gap-labelling theorem open? Solved for λ ≠ ±1 (Avila–You–Zhou,
  arXiv:2306.16254); critical case open. Secondary open frontiers: exact
  Hausdorff dimension of the critical spectrum (≤ 1/2 by
  Jitomirskaya–Krasovsky; lower bound open, dimension is frequency-dependent)
  and the Thouless bandwidth constant (q·|σ| → 32C/π, C Catalan's constant —
  unproven in general).
- **Exact-arithmetic fit.** At rational flux p/q, Chambers' relation reduces
  the spectrum to {E : |P(E)| ≤ 4} with P ∈ ℤ[E] monic of degree q. P is
  computable exactly (transfer-matrix trace over ℤ[ζ_q], reduced mod the
  cyclotomic polynomial), and band edges / gap widths are then certified by
  Sturm-chain root isolation over ℚ. Everything about rational flux is
  `VERIFIED`-able; everything about irrational α is `EVIDENCE` via
  convergents, and the records must say so.
- **Harness**: `harness/almost-mathieu/chambers.py` (exact Chambers
  polynomial + certified bands/gaps), `harness/almost-mathieu/verify_bands.py`
  (independent route: multipoint evaluation + interpolation for P, separate
  root-isolation code). Both stdlib-only.
- **First attempts (queue lines):**
  1. Self-test: rational-flux gap census for all p/q, q ≤ 30 — re-derive the
     published rational-flux facts (all gaps open except the central touching
     at even q, van Mouche 1989 / Choi–Elliott–Yui 1990) with our own exact
     tooling. Expected `VERIFIED`, scope = the range of q run.
  2. Golden-mean convergents: exact minimal-gap-width and q·|σ| tables along
     Fibonacci p/q — scaling evidence toward critical Dry Ten Martini and the
     Thouless constant. `EVIDENCE`, scoped by the convergents reached.
  - Kill condition to record at queue time: if exact P computation past
    q ≈ 100 is out of reach even with a C kernel, the census route stalls at
    small denominators and says nothing asymptotic — record the wall.

### 2. `three-phase-conductivity` — optimal three-phase 2D composites

- **Headline open problem.** Characterize the attainable isotropic effective
  conductivity of a 2D composite of three isotropic phases in prescribed
  volume fractions — equivalently, close the gap between the best proved
  bounds (Hashin–Shtrikman 1962; Nesi 1995; Cherkaev 2009) and the best
  constructions (hierarchical laminates, wheel assemblages). Known: HS is
  attainable for two phases; for three phases it is *not* attainable when the
  best-conductor fraction is small; the exact attainability frontier is open.
- **Exact-arithmetic fit.** Effective tensors of hierarchical laminates with
  rational directions, fractions and conductivities are rational — the
  lamination formula is a rational function, so `fractions.Fraction` computes
  them exactly, and the classical bounds are rational formulas too.
  Attainment claims are exact equalities in ℚ; census claims are scoped by
  the structure class (lamination rank) and grid.
- **Harness**: `harness/three-phase-conductivity/laminate.py` (exact
  effective tensors of hierarchical laminates, volume-fraction accounting,
  Wiener and HS bounds, Keller–Dykhne duality check),
  `harness/three-phase-conductivity/verify_laminate.py` (independent
  re-computation via the dual/resistivity route plus independently coded
  bounds).
- **First attempts (queue lines):**
  1. Self-test: two-phase ground truth — rank-2 laminates attaining the 2D
     two-phase HS bound exactly in ℚ, Keller duality on nontrivial trees,
     series/parallel closed forms. Expected `VERIFIED`, scope = the checks
     run.
  2. Three-phase bounded-rank laminate census over a rational grid of volume
     fractions at fixed conductivities: map where laminates reach HS, where
     they stall against Nesi/Cherkaev-type behavior — a `MAP`-grade
     attainability landscape. Implement the Nesi/Cherkaev improved bounds as
     attempt tooling (they need careful transcription from the papers — mark
     transcription per CONTRIBUTING rule 6).
  - Kill condition: if bounded-rank laminate optima plateau strictly inside
    the bounds on the whole grid with no structure in the gap, the census
    route yields only a negative map — still recordable, but cap the budget.

## What onboarding does NOT claim

Neither problem being solved would produce a superconductor. The honest chain
is: these are the open mathematical cores of two physics programs, chosen
because progress on them is checkable to this lab's standard. The physics
motivation lives in the explainers and this file; `PROBLEM.md` carries only
published mathematics.
