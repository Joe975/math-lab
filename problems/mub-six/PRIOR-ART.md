# MUBs in dimension six — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`. Full records: `attempts/`.

## Editorial view of the attack surface

- Numerical searches for four MU bases are well trodden in the literature;
  a new negative search adds little unless it is calibrated or certified.
- The productive directions are certification (prove the numerical floor is a
  true minimum, ideally a global one within a family) and structured searches
  seeded from the known families of 6×6 complex Hadamard matrices.

## Attempts

### 001 — Multi-start census of the four-basis defect landscape · `EVIDENCE`

22,000 random starts of projected-gradient descent on four bases of C^6: no
start reaches L = 0. 70.7% of the 20,000-start run converge to
L = 0.0512492189962838 (ASD 0.998291692700), which reproduces the published
Raynal–Lü–Englert optimum (rediscovery, not new). At the floor all three
non-standard bases are Hadamard (unbiased to the standard basis). The
residual defect is shared equally by the three Hadamard–Hadamard pairs, and
each row of their overlap matrices is a permutation of {0.124398, 0.151402,
0.18105 ×4}.

Controls under identical settings: d = 2 sits at the proved floor 2/3; known
MU sets are found in d = 3, 4, 5, 7 and for three bases in d = 6. The d = 7
hit rate is only 0.7%, which is the main caveat on reading the d = 6 null.
Leads: certify the floor; closed forms for the three overlap values;
Hadamard-family seeds; equal-sample controls.

### 002 — Skeptic re-implementation of 001 · `VERIFIED` (001's claims, range only)

Independent code (exponential map + L-BFGS, different RNG and stopping rule;
written without seeing 001) reproduces the floor 0.0512492189962838 in
1431/2000 and 717/1000 starts, the secondary minima, the hub-basis structure,
the three-value overlap pattern, and all controls. No corrections. One
precision: the hub is *some* basis, not necessarily the fixed one. New: the
same value is also a 3-basis local minimum in d = 6 (about 10% of k = 3
starts, in both implementations).

### 003 — Floor Hessian and closed forms · `EVIDENCE` (+ exact algebra)

- **Second order.** The four-basis floor is a nondegenerate local minimum
  modulo the 23-dimensional gauge orbit (float Hessian, 108 coordinates).
- **Closed forms.** Polished to 50 digits, the three overlap values are
  4c²/3, (2c − 1)² and c(3 − 4c)/3, with c = cos²θ the root of
  112c³ − 144c² + 63c − 9 in [0, 1].
  - LLL found these at 13 digits; they were confirmed at 50.
  - Along this pattern, L(c) = 448c⁴ − 768c³ + 504c² − 144c + 15 exactly,
    and L′ is 16 × that cubic.
  - The cubic is RLE's Eq. 20, rediscovered blind and matched to their radical
    to about 1e-51.
- **Three bases.** The same L* is a one-parameter valley of
  gauge-inequivalent configurations: equal overlap moduli, varying
  triple-product phases. The fourth "hub" basis isolates one point of it.

### 004 — Skeptic Hessian review of 003 · `VERIFIED_WITH_CORRECTIONS`

- An exact automatic-differentiation Hessian in a Cayley chart (long double)
  confirms the signature: 23 gauge zeros and 85 positive eigenvalues. It also
  confirms the overlap values to 18 digits.
- Correction: 003's eigenvalue magnitudes depend on the generator
  normalisation and on which basis is fixed. Two correct implementations
  disagreed until both were pinned down; only the signature is invariant.
- Not reviewed: the k = 3 valley, and the closed forms (which are exact
  arithmetic).
