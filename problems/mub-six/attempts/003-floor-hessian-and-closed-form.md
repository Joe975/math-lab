# 003 — The four-basis floor: second-order check and exact overlap values

- **Problem:** mub-six, `problems/mub-six/PROBLEM.md`
- **Date:** 2026-10-02
- **Mode:** blind
  (Same session as 001 and 002, so it is informed by those two records. It is
  blind in the lab's sense: no prior art existed beyond them. One published
  detail was read *after* the closed forms had been found from the numerics:
  Raynal–Lü–Englert's Eq. 20–21, used as an independent check. That order is
  recorded under "What was done".)
- **Type:** numerical second-order analysis + integer-relation detection +
  exact polynomial algebra
- **Tools:** `problems/mub-six/explore/hessian.py` (float Newton polish,
  finite-difference Hessian, explicit gauge subspace, Jacobi eigensolver;
  ~30 s); `explore/hiprec_polish.py` (50-digit Decimal polish, ~8 s);
  `explore/intrel.py` (LLL integer relations, exact integers);
  `explore/relations_check.py`; `explore/pattern_family.py` (exact Fraction
  algebra). All standard-library Python, deterministic. Outputs in
  `problems/mub-six/data/{hessian_d6k4,hiprec_d6k4,relations_check,pattern_family}.txt`
  and `data/polished_d6k4.json`.
- **Sources:** Raynal–Lü–Englert arXiv:1103.1025, Eq. 19–22 [T] (machine
  summary of the ar5iv HTML; read after the cubic was found).

## Approach

001/002 leads 1 and 2: is the floor L* = 0.0512492189962838 a genuine local
minimum of L on the full configuration space (not only inside the published
family where it was found), and do the three overlap values at the floor have
closed forms?

The second-order check comes first. A floor reached by 70% of random starts
could still be a saddle that first-order methods stall on, or the bottom of a
flat valley. Either would change how a negative search should be read. The
closed-form question is cheap once the point is polished to many digits,
because LLL gives exact-integer candidates that can then be tested at much
higher precision than was used to find them.

## What was done

**1. Second-order check** (`hessian.py data/best_d6k4_20k.json`).
- Coordinates: U_b → U_b exp(A_b) for the three moving bases, with A_b
  anti-Hermitian on the standard real basis of u(6), giving N = 108.
- Gauge subspace, built explicitly as tangent vectors: the column phases of
  each moving basis (3 × 6 = 18) and a common left diagonal unitary
  D·U_b (6, of which the global phase coincides with column phases). These
  leave every |⟨a_i|b_j⟩| fixed. Gram–Schmidt gives rank **23**.
- One Newton step on the complement of the gauge subspace (pseudo-inverse)
  takes ‖grad‖ from 7.3e-9 to 1.3e-15.
- Hessian by central differences of the analytic gradient (h = 1e-5),
  symmetrised. At a critical point the Hessian does not depend on the chart.
- Results:
  - The full Hessian has exactly **23** eigenvalues with |λ| < 1e-6 (all
    below 1.5e-10 in magnitude), equal to the gauge rank.
  - ‖H v‖ ≤ 2.4e-10 for every gauge vector v.
  - On the 85-dimensional complement, every eigenvalue is positive:
    **min 0.21488**, max 10.263. The full sorted list is in
    `data/hessian_d6k4.txt`.
  - The spectrum has heavy multiplicities (many 2- and 4-fold eigenvalues),
    consistent with a configuration that has a large discrete symmetry.

**2. High-precision polish** (`hiprec_polish.py data/polished_d6k4.json 50`).
- Decimal arithmetic at 60 working digits. The float Hessian is a fixed
  preconditioner; gradient, retraction and re-orthonormalisation are done in
  Decimal.
- ‖grad‖∞ goes 3.9e-16 → 1.0e-25 → 3.2e-35 → 1.1e-44 → 3.8e-54, linear
  convergence with ratio ~1e-10, as expected from a float preconditioner.
- At 50 digits:

      L* = 0.051249218996283827782882716219040904801850493547422
      p1 = 0.12439794566570321461037700050123282717085012795674  (x18)
      p2 = 0.15140197747747829881064201194373798309990083907013  (x18)
      p3 = 0.18105001921420462164474524688875729743231225824328  (x72)
      1/6                                                         (x108)

  The 108 entries equal to 1/6 are the three unbiased pairs. The 108 others
  are the three defective pairs, 36 entries each.

**3. Integer relations** (`intrel.py`, `relations_check.py`).
- LLL at **13 digits** proposed these cubics, with coefficients at most 784
  against ~10^3.25 expected for a chance relation:

      784 p1^3 - 552 p1^2 + 153 p1 - 12 = 0
      784 p2^3 +  24 p2^2 -  15 p2 -  1 = 0
      392 p3^3 - 228 p3^2 +  45 p3 -  3 = 0

- All three then hold at **50 digits**, with residuals 4.5e-50, 1.0e-49 and
  3.2e-51.
- A degree-4 candidate for L found at 13 digits (10L² + 19L − 1 = 0) **failed**
  at the 9th digit. It is recorded here as a reminder that the precision-limit
  test matters.
- At 45 digits, L* satisfies 5488 L³ − 7848 L² + 6183 L − 297 = 0
  (residual 2e-48).
- Searching for each value in the basis {1, s, s²}, where s = sin²θ_opt as
  defined by RLE's ASD formula, gives:

      p1 = 4c²/3,   p2 = (2c − 1)²,   p3 = c(3 − 4c)/3,   with c = 1 − s = cos²θ,

  with residuals ≤ 1e-50.

**4. Exact algebra on the pattern** (`pattern_family.py`, Fractions).
- p1 + p2 + 4p3 = 1 holds identically, so each row is a probability vector.
- With one defective pair's defect f = 6[(p1 − 1/6)² + (p2 − 1/6)² +
  4(p3 − 1/6)²] and L = 3f:

      L(c) = 448c⁴ − 768c³ + 504c² − 144c + 15
      L'(c) = 16 · (112c³ − 144c² + 63c − 9)        (exact identity)

- L' has exactly one sign change on [0, 1], at c* = 0.30544796487990783.
  L''(c*) = 102.07 > 0, and L(c*) reproduces L*.
- In s = 1 − c the cubic is 112s³ − 192s² + 111s − 22 = 0.

**5. The same value with three bases (002's observation; lead 3 run here).**
- `data/floor_d6k3_seed21.json` is a k = 3 endpoint at L*, found by running
  single starts of `mub_search 6 3 1 <seed>` for seeds 1, 2, … until one
  landed at L*; seed 21 was the first.
- `hessian.py` on it gives gauge rank **17** (12 column phases + 5) but
  **18** near-zero eigenvalues: one non-gauge null direction. The other 54
  complement eigenvalues lie in [0.250, 7.89].
- `flat_direction.py` steps along that null vector by t and re-minimises
  over everything except gauge and that vector:
  - raw L − L* grows like t⁴ (4.6e-14 at t = 1e-3, 4.6e-10 at 1e-2);
  - relaxed L − L* stays at roundoff (≤ 7e-17) for t up to 0.6;
  - the overlap moduli are unchanged along the way (to the ~1e-9 relaxation
    tolerance).
- `flat_direction.py … --bargmann` compares the gauge-invariant triple
  products ⟨a_i|b_j⟩⟨b_j|c_k⟩⟨c_k|a_i⟩ (as sorted multisets) between the
  start and the relaxed points. They change by up to 3.6e-3 (t = 0.1),
  1.8e-2 (t = 0.3) and 2.5e-2 (t = 0.6). The relaxed points are therefore
  genuinely different configurations, not gauge copies.

**6. Comparison with the literature, done after steps 1–5.**
- RLE Eq. 20 [T] is 112p⁶ − 192p⁴ + 111p² = 22 with p² = sin²θ_opt, the same
  cubic.
- Their radical solution [T], r = (21√3 − 36)^{1/3},
  sin²θ_opt = (3 + 16r − r²)/(28r), agrees with our s to 1.2e-52.
- Not circular: the polished overlaps p1, p2, p3 (from our numerics alone)
  match 4c²/3, (2c − 1)², c(3 − 4c)/3 at RLE's radical c to ≤ 3.3e-51.
- The quartic L(c) at that c matches the polished L* to 3.8e-52.

## Outcome

`EVIDENCE` for the second-order claim; exact algebra for the pattern
identities.

1. **Second order (EVIDENCE, float).** At the floor configuration of 001 the
   Hessian is positive definite on the 85-dimensional complement of the
   23-dimensional gauge orbit, with smallest eigenvalue 0.215 in the
   coordinates above, and the zero modes are exactly the gauge directions.
   So the floor is, numerically, a **nondegenerate local minimum of L on the
   full configuration space**, isolated up to symmetry. It is not a saddle,
   and it is not a point on a flat valley. Scope: one configuration; float
   finite differences with error ~1e-10 against a spectral gap of 0.2.
2. **Overlap pattern (EVIDENCE, 50 digits).** The floor's overlap values are
   4c²/3, (2c − 1)² and c(3 − 4c)/3 at c = cos²θ_opt, the root of
   112c³ − 144c² + 63c − 9 in [0, 1], with ≥ 49 digits of agreement.
3. **Three bases (EVIDENCE, float).** With k = 3, the same L* is attained
   along a **one-parameter family** of gauge-inequivalent configurations.
   They share the overlap moduli and differ in triple-product phases. So
   the three-basis floor is a flat valley, and adding the fourth "hub" basis
   (unbiased to the other three) cuts it to an isolated point. Scope: one
   valley, followed to t = 0.6 from one start.
4. **Pattern algebra (VERIFIED, exact arithmetic).** Given that pattern,
   L(c) = 448c⁴ − 768c³ + 504c² − 144c + 15, and its critical-point equation
   is exactly RLE's Eq. 20.

**Rediscovery, not new:** the cubic (RLE Eq. 20) and the value of the floor
(RLE Eq. 22) are published. The cubic was found here blind from
random-start numerics, then matched. **Possibly new, not claimed as new:**
- the closed forms of the three transition probabilities (the machine summary
  of RLE [T] says the paper does not state them; the PDF was not checked);
- the observation that the floor is the minimum of a single quartic along a
  fixed overlap pattern;
- the full-space Hessian spectrum (RLE optimise inside their family);
- the k = 3 valley. It may well be RLE's own free parameter seen from the
  three-basis side; that identification was not checked.

**Not claimed:**
- That L* is the global minimum of L.
- That the Hessian claim is a certificate. It is floats; see lead 1.
- That the pattern family is realisable for every c. Only c* is
  realised here, by the numerical configuration.

## Why it failed / what survived

Nothing failed. What survives:
- The floor is a strict local minimum modulo gauge, numerically. So a search
  that lands there is genuinely trapped, not stalled on a saddle.
- The mechanism has a geometric reading: a valley of frustrated three-basis
  configurations, from which the fourth basis selects a point. A route to
  four MU bases would have to leave the valley, not move along it.
- 001's "dominant funnel" now has an exact description. The funnel floor is
  the bottom of a quartic along a one-parameter overlap pattern.

**What this does not touch.** Whether a *different* basin reaches L = 0. A
strict local minimum says nothing about distant regions. The d = 7 caveat
from 001 stands unchanged.

**Reusable:**
- `hessian.py`: gauge-aware Hessian for any (d, k) configuration file.
- `hiprec_polish.py`: 50+ digit polish of any nondegenerate critical point.
- `intrel.py`: LLL relation finder, with the two-precision acceptance rule
  documented in its docstring.
- The method generalises: polish, LLL at low precision, confirm at high
  precision, then do exact algebra on the detected pattern. It turns a
  numerical floor into an exact algebraic statement.

## Leads generated

1. **Certificate.** Build the floor configuration exactly in the number field
   Q(c*) (cubic) using RLE's explicit family. Verify unitarity and the
   overlap pattern symbolically. Then bound the Hessian's smallest
   non-gauge eigenvalue below 0 via interval arithmetic. That would upgrade
   claim 1 to a proof that L* is a strict local minimum.
2. **Second floor.** Apply the same polish → LLL → algebra pipeline to the
   next local minima of 001's census (0.17887, 0.19565, 0.23401). If they
   also have closed forms on overlap patterns, the landscape near L = 0 may
   be classifiable pattern by pattern.
3. **Identify the k = 3 valley.** Parametrise it by one triple-product
   phase and check whether it is RLE's Fourier-transposed family with the
   constraint of their Eq. 7 [T]. Check also whether the point the hub basis
   selects is the unique point of the valley that admits a common unbiased
   basis. Both are definite yes/no computations.

## References

- P. Raynal, X. Lü, B.-G. Englert, *Mutually unbiased bases in dimension 6:
  The four most distant bases*, Phys. Rev. A 83 (2011), arXiv:1103.1025,
  Eq. 19–22 [T].
- A. K. Lenstra, H. W. Lenstra, L. Lovász, *Factoring polynomials with
  rational coefficients*, Math. Ann. 261 (1982) (LLL).
- `problems/mub-six/attempts/001-multistart-defect-landscape.md`,
  `002-skeptic-reimplementation-of-001.md`.
