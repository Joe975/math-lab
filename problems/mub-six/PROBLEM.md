# Mutually Unbiased Bases in Dimension Six

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** Two orthonormal bases {a_i} and {b_j} of C^d are *mutually
unbiased* (MU) if |⟨a_i|b_j⟩|² = 1/d for every i, j: measuring a state
prepared in one basis in the other gives a perfectly uniform outcome. At most
d + 1 bases of C^d can be pairwise MU. The question: **how many MU bases exist
in C^6?** The conjecture (Zauner 1999) is that the answer is **3**, so in
particular that no 4 MU bases exist.

## Why qubits care

For d = 2 (one qubit) the three Pauli eigenbases (X, Y, Z) are MU; on the
Bloch sphere, a basis is an antipodal pair of points and two bases are MU
exactly when their axes are perpendicular, so 3 = d + 1 is the maximum.
Complete sets of d + 1 MU bases are optimal for quantum state tomography and
underlie several quantum-key-distribution protocols. d = 6 = 2 × 3 (a qubit
and a qutrit together) is the smallest dimension that is not a prime power,
and it is where the known constructions stop working.

## Published status

- **Prime powers.** For d = p^n a complete set of d + 1 MU bases exists
  (Ivanović 1981; Wootters–Fields 1989).
- **Composite d.** Taking tensor products of prime-power constructions gives at
  least p_min^{n} + 1 bases, where p^n is the smallest prime-power factor of d.
  For d = 6 this gives 3, and that is all that is known to exist.
- **Numerical evidence for d = 6.** Butterley–Hall (quant-ph/0701122)
  minimised a non-negative defect function over sets of bases and found no 4
  MU bases. Brierley–Weigert (arXiv:0808.1614, 0901.4051) searched for MU
  *constellations* (partial bases) and found only 18 of the 35 that a complete
  set would imply. Raynal–Lü–Englert (arXiv:1103.1025) found the maximal
  average squared distance of four bases in d = 6 to be ≈ 0.9983 < 1, attained
  inside the Fourier-transposed family.
- **Partial no-go theorems.** Jaming–Matolcsi–Móra–Szöllősi–Weiner
  (arXiv:0902.0882) proved that the standard basis together with any member
  of the Fourier family F(a, b) cannot be extended to four MU bases; other
  families of 6 × 6 complex Hadamard matrices have been excluded case by
  case. The full classification of 6 × 6
  complex Hadamard matrices is itself open.

## Standard reduction

WLOG the first basis is the standard basis, and then every other basis of an
MU set is (1/√6)·H for a 6 × 6 *complex Hadamard matrix* H (unimodular
entries, H H† = 6 I). Four MU bases ⇔ three complex Hadamard matrices H₁, H₂,
H₃ with every H_a† H_b / 6 also (1/√6)·Hadamard.

A standard defect function, for k bases B_1..B_k of C^d:

    L = Σ_{a<b} Σ_{i,j} ( |⟨a_i|b_j⟩|² − 1/d )²

L ≥ 0, with L = 0 iff the bases are MU. The Bengtsson distance between two
bases is D²_ab = 1 − (pair term of L)/(d − 1), and the average squared
distance of k bases is ASD = 1 − L / (C(k,2)·(d − 1)).

## What would refute the conjecture, and what checks it

A **refuting object** is 4 orthonormal bases of C^6 that are pairwise MU — a
finite, checkable object (three 6 × 6 matrices). Numerically it must have
L at roundoff level *and* survive an exact or certified check: either exact
algebraic entries (cyclotomic or other number field) verified symbolically,
or interval arithmetic proving a genuine solution of the polynomial system
lies within the stated tolerance. A float configuration with small L is not a
refutation; the defect function has near-misses.

A **proof** of the conjecture is a proof step, not an object, and no finite
computation is known to suffice, because the space of 6 × 6 complex Hadamard
matrices is not classified.

## Verification contract

- Any configuration claim ships its bases as JSON and is re-checked with
  `harness/mub-six/mub_check.py`, which recomputes orthonormality, L, the
  per-pair defects and the ASD independently of whatever produced the file.
- A multi-start search result is `EVIDENCE` about the starts actually run:
  state the number of starts, the seed, the stopping rule, and the
  distribution of local minima, not just the best value.
- A "minimum" found numerically is a local-search value, never a certified
  global minimum, unless accompanied by a certificate.

## Harness (tier 0)

- `harness/mub-six/mub_check.py` — standard-library checker for a set of
  bases stored as JSON (`{"d": .., "bases": [[[re, im], ...], ...]}`, where
  `bases[b][c]` is vector c of basis b).
