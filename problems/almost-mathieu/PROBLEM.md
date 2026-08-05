# Gaps of the critical almost Mathieu operator

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** The almost Mathieu operator on ℓ²(ℤ) is

  (H_{λ,α,θ} u)_n = u_{n+1} + u_{n−1} + 2λ cos(2π(αn + θ)) u_n,

with coupling λ > 0, frequency α, and phase θ. At the **critical coupling
λ = 1** it is the Harper model: the Landau-gauge tight-binding model of an
electron on ℤ² threaded by magnetic flux α per unit cell, whose spectrum drawn
against α is the Hofstadter butterfly. For irrational α the spectrum is a
θ-independent compact set σ(λ, α).

The **gap-labelling theorem** says that on any spectral gap of H_{λ,α,·} with
α irrational, the integrated density of states takes a value {kα} (the
fractional part of kα) for some integer k ≠ 0; the gap with label k either is
open or has collapsed to a point. The **Dry Ten Martini problem** asks:

> For every irrational α, is every gap allowed by gap labelling actually
> *open* in the spectrum of the critical operator H_{1,α,·}?

For the non-critical operator (λ ≠ 1) the answer is known to be yes. The
critical case is open, and is the problem here.

## Published status

- **Ten Martini (Cantor spectrum).** σ(λ, α) is a Cantor set for every
  irrational α and every λ ≠ 0 — A. Avila, S. Jitomirskaya, *The Ten Martini
  Problem*, Ann. of Math. 170 (2009). This does not decide whether individual
  labelled gaps are open.
- **Dry Ten Martini, non-critical case.** All labelled gaps are open for all
  irrational α whenever λ ≠ 1 — A. Avila, J. You, Q. Zhou, *Dry Ten Martini
  Problem in the non-critical case*, arXiv:2306.16254. Earlier partial
  results: Choi–Elliott–Yui (Liouville α), Avila–Jitomirskaya (Diophantine α,
  via almost reducibility), Liu–Yuan. **The critical case λ = 1 is open for
  every single irrational α.**
- **Measure of the critical spectrum.** |σ(1, α)| = 0 for every irrational α
  (Avila–Krikorian, following Last's a.e. result). For rational α = p/q in
  lowest terms the spectrum is a union of q closed bands, and the total
  bandwidth obeys |σ(1, p/q)| > 0 with |σ(1, p/q)| → 0 as q → ∞.
- **Thouless bandwidth conjecture.** Thouless (1990) argued, with strong
  numerical support, that q·|σ(1, p/q)| → 32C/π as q → ∞, where
  C = Σ (−1)^k/(2k+1)² ≈ 0.9160 is Catalan's constant. Rigorous upper bounds
  of the conjectured order exist (Jitomirskaya–Konstantinov–Krasovsky, who
  also derive a new Chambers-type formula in the rational case,
  arXiv:2007.01005); the limit itself is unproven. The n-th-moment
  generalization is studied in Ouvry–Wu (arXiv:1703.09634).
- **Hausdorff dimension.** dim_H σ(1, α) ≤ 1/2 for every irrational α —
  S. Jitomirskaya, I. Krasovsky, arXiv:1909.04429. There is a dense G_δ of
  frequencies where the dimension is 0 (Last–Shamis), and frequency sets where
  it is positive; the exact value for a.e. α (Thouless's heuristic, tied to
  the bandwidth conjecture, predicts 1/2) is open.
- **Rational flux.** For α = p/q the analysis is exact. Chambers' relation:
  the discriminant of the period-q potential v_j = 2λcos(2πpj/q + ν)
  satisfies Δ_ν(E) = P(E) ± 2λ^q cos(qν) where P does not depend on ν, so
  σ(λ, p/q) = { E : |P(E)| ≤ 2 + 2λ^q }. At λ = 1, P is monic of degree q
  with coefficients in ℤ[2cos(2π/q)] — the ring of integers of the real
  cyclotomic field ℚ(ζ_q)⁺, rational precisely when φ(q) ≤ 2, i.e.
  q ∈ {1, 2, 3, 4, 6} (already at flux 1/5,
  P(E) = E⁵ − 10E³ + (15 − 10cos72°)E). The band edges are algebraic
  numbers — roots of P(E) = ±4. All q − 1 gaps of σ(1, p/q) are open
  **except** the central one: for q even the two middle bands touch at E = 0
  and nowhere else (P. van Mouche, Comm. Math. Phys. 122 (1989);
  Choi–Elliott–Yui, *Gauss polynomials and the rotation algebra*, Invent.
  Math. 99 (1990)).
- **Physical background.** The spectral edge of exactly this operator family
  sets the critical-temperature curve T_c(B) of a superconducting wire
  network via the de Gennes–Alexander network equations, an effect measured
  in aluminum networks (Pannetier, Chaussy, Rammal, Villégier, Phys. Rev.
  Lett. 53 (1984)); the same operator is the quasi-energy operator of the
  kicked Harper model in quantum chaos, and its spectrum has been imaged in
  cold-atom, photonic, moiré-graphene and superconducting-qubit platforms.

## Verification contract

Any claim recorded against this problem must meet the bar in `CONTRIBUTING.md`.
The split that keeps this problem honest is rational versus irrational flux:

- **Rational flux is exact or it is nothing.** A claim about σ(1, p/q) rests
  on the Chambers polynomial P computed in exact arithmetic over the real
  cyclotomic field ℚ(2cos(2π/q)), with the computation stating how the field
  elements were represented and how signs were decided. Band edges are
  reported as isolating rational intervals with certified
  one-root-per-interval counts (for a degree-q polynomial, q disjoint
  certified sign-change brackets is a complete count); a gap is "open" only
  with a certified positive rational lower bound on its width, and "closed"
  only with an exact multiple-root certificate. Floating point may screen,
  never certify, and every use of it must be labelled as a screen.
- **Every claim states its range**: the exact set of (p, q) computed. A
  census over q ≤ Q is `EVIDENCE` about those denominators and nothing more.
- **Irrational α is never reached by computation.** Band or gap data along a
  sequence of convergents p_n/q_n is `EVIDENCE` about a trend; the spectrum
  at irrational α is not the limit-in-any-naive-sense of rational spectra
  (gap edges converge, but band structure reorganizes at every level), and
  records must not conflate the two. Box-counting or bandwidth scaling at
  rational approximants does not bound dim_H σ(1, α); it is heuristic
  support only and is labelled `SPECULATION` wherever it carries weight.
- **Independent verification.** Any P(E), band edge, or gap bound intended
  for the ledger is re-derived by the second implementation
  (`verify_bands.py`), which must reach P by a different route and count
  roots by a different algorithm, before it is recorded anywhere outside the
  attempt record.

## Harness (tier 0)

- `harness/almost-mathieu/chambers.py` — reference implementation, critical
  coupling only. Computes P exactly in K[E], K = ℚ[y]/ψ_q(y) with
  y = 2cos(2π/q) and ψ_q extracted from the cyclotomic polynomial Φ_q, via
  the transfer-matrix product at phase ν = 0 (the potential values are then
  Chebyshev-style integer polynomials in y). The ν-structure of Chambers'
  relation is not assumed but pinned by two exact cross-phase identities:
  the same product recomputed in ℤ[x]/(x^{4q} − 1) at ν = π/(2q) and in
  ℤ[x]/(x^{2q} − 1) at ν = π/q must differ from the embedded K-result by the
  constants ±2 and ±4 with one consistent sign, which certifies P and the
  threshold |P| ≤ 4 under either sign convention. Band edges: exact integer
  deflation for rational roots (including the even-q touching at E = 0),
  then exact sign-change brackets found by bisection — with signs of
  K-elements decided exactly through a refinable rational enclosure of y —
  and refined to width 2⁻ᵇⁱᵗˢ. Output is the certified band list, per-gap
  open/touching certificates with rational width bounds, and a rational
  enclosure of q·|σ|.
- `harness/almost-mathieu/verify_bands.py` — independent re-computation. P is
  reached by a different route entirely: exact scalar evaluation of the
  discriminant at q + 1 integer energies in ℤ[x]/(x^{4q} − 1) at the phase
  ν = π/(2q), followed by exact Lagrange interpolation, compared against the
  claimed P through the subfield embedding y ↦ x⁴ + x^{4q−4}. Certificates
  are then re-checked with an independent sign oracle whose enclosure of
  2cos(2π/q) comes from certified rational series bounds (Machin's formula
  for π and the alternating cosine series) rather than from ψ_q. Run it on
  any spectral data you intend to claim; disagreement between the two tools
  is reported as requiring escalation rather than resolved by picking one.
