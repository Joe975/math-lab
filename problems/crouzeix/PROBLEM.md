# Crouzeix's Conjecture

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** For every square complex matrix A and every polynomial p,

  ‖p(A)‖ ≤ 2 · max{ |p(z)| : z ∈ W(A) },

where ‖·‖ is the spectral norm and W(A) = {x*Ax : x ∈ ℂⁿ, ‖x‖ = 1} is the
numerical range (field of values) of A — equivalently, that W(A) is a
2-spectral set for A. Conjectured by M. Crouzeix in 2004. The quantity being
bounded, ‖p(A)‖ / max_{W(A)} |p|, is the **Crouzeix ratio**.

W(A) is a compact convex subset of ℂ containing the eigenvalues
(Toeplitz–Hausdorff), and it is the intersection of the half-planes
{z : Re(e^{−iθ}z) ≤ λ_max(H_θ)} over θ ∈ [0, 2π), where H_θ is the Hermitian
part of e^{−iθ}A.

## Published status

Open in general. The constant 2 is known to be the smallest possible: for
A = [[0, 1], [0, 0]] and p(z) = z, ‖p(A)‖ = 1 while W(A) is the disc of radius
1/2, so the ratio is exactly 2. The conjecture is therefore a statement that
this 2×2 example is extremal.

- **Best proved constant.** 1 + √2 ≈ 2.414 — M. Crouzeix, C. Palencia, *The
  numerical range is a (1+√2)-spectral set*, SIAM J. Matrix Anal. Appl. 38
  (2017). Crouzeix's own earlier bound (2007) was 11.08.
- **Proved cases.** 2×2 matrices; normal matrices; matrices whose numerical
  range is a disc (via a theorem of Okubo–Ando); matrices that are nearly
  Jordan blocks; tridiagonal 3×3 matrices whose numerical range is an ellipse
  centred at an eigenvalue (Glader–Kurula–Lindström, 2018); and, for the
  complete-2-spectral-set form, weighted shift matrices.
- **Dimension-dependent improvements.** Malman–Mashreghi–O'Loughlin–Ransford
  (2024): for each fixed dimension N there is a constant C_N < 1 + √2 valid for
  all N×N matrices. A 2025 follow-up connects the bound to configuration
  constants for the Neumann–Poincaré operator, giving domain-specific
  improvements.
- **Reformulations.** arXiv:2306.12183 links the general conjecture to
  cyclicity — it holds in full generality if and only if it holds for the
  differentiation operator on a class of analytic functions — and shows that in
  the 3×3 symmetric case it is equivalent to the statement for analytic
  truncated Toeplitz operators.
- **Numerical work.** Greenbaum and Overton investigated the Crouzeix ratio by
  nonsmooth optimization over (A, p); no counterexample has been found, and the
  computed local maxima are reported to sit at known extremal structures.

## Verification contract

Any claim recorded against this problem must meet the bar in `CONTRIBUTING.md`.
This problem has the weakest exact-arithmetic fit of the problems here — both
‖p(A)‖ and max_{W(A)}|p| are limits of eigenvalue/singular-value problems, not
rational quantities — so the contract is built on certified enclosures rather
than on exactness:

- **Every reported ratio is an interval**, with the method that proved the
  enclosure stated (for example: a residual or Gershgorin bound around a
  floating-point eigenpair; interval arithmetic on a characteristic polynomial
  with isolated roots). A bare floating-point ratio is a candidate, never a
  claim.
- **Upper-bounding max_{W(A)}|p| requires an outer enclosure of W(A).**
  Sampling θ finitely and using *upper* bounds on λ_max(H_θ) yields a finite
  intersection of half-planes that provably contains W(A); the number of
  samples and the eigenvalue bound used must both be stated. An inner
  polygonal approximation, which is what connecting boundary points gives, is
  not valid for this direction.
- **A refutation claim (ratio > 2) requires the interval's lower bound to
  exceed 2** — that is, a certified lower bound on ‖p(A)‖ together with a
  certified outer enclosure of W(A) — and must be reproduced by an independent
  implementation using a *different* algorithm before it is recorded anywhere
  outside the attempt record.
- **Optimization landscapes are `EVIDENCE`**, scoped by matrix dimension,
  polynomial degree, and search design. A claim of the form "all local maxima
  terminate at known structure" is a claim about the sampling design, which
  must be stated in full; it is never a claim about the conjecture.
- If certified enclosures prove too expensive to run at the intended scale, the
  correct response is to narrow the claim to structured families where ‖p(A)‖
  has a closed form, and to record that narrowing explicitly. Silently
  substituting floating-point results for certified ones fails this contract.

## Harness (tier 0)

- `harness/crouzeix/ratio.py` — the reference implementation. Certified
  enclosures of the Crouzeix ratio with no floating point anywhere: eigenvalue
  bounds come from bisecting an exact rational LDL definiteness test, and
  W(A) is enclosed by half-planes in Gaussian-rational directions, which is
  possible because both sides of the support inequality scale with the
  direction and so it need not be a unit vector. The refutation test is
  decided in ℚ as ‖p(A)‖² > 4·(max|p|)², taking no square roots at all.
- `harness/crouzeix/verify_ratio.py` — independent re-computation. It factors
  nothing: characteristic polynomials by the Faddeev–LeVerrier recursion,
  largest roots by Sturm chains, and the characteristic polynomial itself
  checked against a determinant computed by elimination before being trusted.
  Run it on any near-extremal ratio you intend to claim; a ratio certified
  above 2 by both routes is reported as requiring escalation rather than
  accepted.
