# Almost Mathieu (critical) — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None.** This problem was onboarded (2026-08-04) and has not been worked.
There is no prior art to be informed by, so `blind` and `informed` mode are
currently equivalent here — which makes the first attempts worth running
blind, since blind costs nothing while the record is empty.

## Editorial view of the attack surface

Spectral theory of quasi-periodic operators is a new attack surface for the
lab, arrived at from a physics direction (superconducting wire networks and
Floquet-driven lattices; see `docs/PLAN-conductance-problems.md` for the
scoping decision). The headline conjecture — Dry Ten Martini at critical
coupling — is a famous long shot: Avila, You and Zhou closed the non-critical
case, and what remains is exactly the hardest point of the parameter space.
Treat it like Collatz: the deliverables are exact rational-flux data, scaling
maps, and tooling, not a proof.

What makes it worth having anyway: rational flux is *fully* exact here. The
harness computes Chambers polynomials over the real cyclotomic field with no
floating point in any certificate, so every gap of every σ(1, p/q) in reach
is a `VERIFIED`-grade fact, and the open questions (critical gaps, Thouless
constant, dimension 1/2) all have falsifiable rational-approximant shadows.

A hard-won onboarding lesson, recorded so nobody re-learns it: the Chambers
polynomial does **not** have integer coefficients once φ(q) > 2 — flux 1/5
already needs ℤ[2cos72°]. The first harness draft assumed ℤ[E] and its own
integrality assertion caught the error. Any tool built here must live in the
real cyclotomic field or prove why it can avoid it.

Concrete lines, if you want them:

- Self-test census: all p/q with q ≤ 30 — certify every gap open except the
  even-q central touching, re-deriving van Mouche / Choi–Elliott–Yui in
  range. Expected `VERIFIED`, scope = the q range run. Watch the smallest
  certified gap width as a function of q: that trend line is the empirical
  shadow of critical Dry Ten Martini.
- Golden-mean convergents: exact minimal-gap and q·|σ| tables along
  Fibonacci p/q as far as the tooling reaches (the harness does q = 34 in
  seconds; a C kernel to the same contract should reach q ≈ 100+). Onboarding
  smoke runs gave q·|σ| = 9.2509 (q=13), 9.3199 (21), 9.3608 (34) against
  the conjectured 32C/π ≈ 9.3299 — oscillating convergence worth pinning
  down with an error model. `EVIDENCE`, scoped by the convergents reached.
- Box-counting the rational approximant spectra against the dimension-1/2
  prediction — strictly heuristic for irrational α (say so inline), but the
  scaling exponent is a falsifiable output.

Kill condition, inline: if exact arithmetic past q ≈ 100 is out of reach even
with a C kernel, the census route stalls at small denominators and says
nothing asymptotic — record the wall and stop rather than extrapolating.
