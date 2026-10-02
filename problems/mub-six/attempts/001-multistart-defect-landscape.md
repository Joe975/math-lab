# 001 — Multi-start census of the four-basis defect landscape in d = 6

- **Problem:** mub-six, `problems/mub-six/PROBLEM.md`
- **Date:** 2026-10-02
- **Mode:** blind
  (The problem was created in this same session and had no lab prior art, so
  there was nothing tier-1 to read. The only sources read were the published
  papers listed below, all of which are tier-0 background.)
- **Type:** computational search (multi-start local optimisation) with
  positive controls and an independent re-implementation
- **Tools:** `problems/mub-six/explore/mub_search.c` (C, deterministic given
  seed; ~9 ms per d = 6 start at 50,000 max iterations);
  `harness/mub-six/mub_check.py` (stdlib checker for saved configurations);
  the independent re-implementation of the search is a separate review
  record (002), so this record states only what one implementation found.
- **Sources:** Butterley–Hall quant-ph/0701122 (abstract only);
  Raynal–Lü–Englert arXiv:1103.1025 [T]; Brierley–Weigert arXiv:0808.1614
  (abstract only); Jaming–Matolcsi–Móra–Szöllősi–Weiner arXiv:0902.0882
  (abstract only). [T] = details read through a machine summary of the ar5iv
  HTML rendering, not the PDF: the definition of D², the ASD normalisation,
  the value 0.9983, the closed form (71 − 12 cos⁴θ)/70, and sin²θ_opt =
  0.6946.

## Approach

Minimise the standard defect

    L = Σ_{a<b} Σ_{i,j} ( |⟨a_i|b_j⟩|² − 1/6 )²

over four orthonormal bases of C^6 (first fixed to the standard basis), from
many random Haar-like starts, and **classify every endpoint** rather than
report only the best value. The point was not to expect a counterexample. It
was to produce a reproducible, *calibrated* negative: the same code and
settings run on cases where the answer is known (d = 2 impossible with a
provable floor; d = 3, 4, 5, 7 where MU sets exist; d = 6 with 3 bases), so
the d = 6, k = 4 negative can be read against measured hit rates.

Why not a structured search (seeding from known complex Hadamard families)?
That is the obvious stronger route and is listed as a lead. Random starts
come first because they carry no assumption about which family a fourth basis
would come from, and because the published optima were found inside one
specific family. A random search reaching the same optimum tests whether that
family's optimum is also the global landscape's attractor.

## What was done

**Optimiser.** Each non-fixed basis is a unitary U_b. The Euclidean gradient
of L with respect to conj(U_b) is assembled pairwise (for G = U_a† U_b,
∂L/∂conj(G_ij) = 2(|G_ij|² − 1/d) G_ij). It is projected to the tangent space
(P = E − U·herm(U†E)), a step U − ηP is taken, and the result is retracted by
modified Gram–Schmidt. The step is accepted only if L decreases (η × 1.3),
otherwise η is halved. A start stops at ‖P‖ < 1e-14, L < 1e-28, η < 1e-18, or
50,000 iterations. Starts are i.i.d. complex-Gaussian matrices
orthonormalised by Gram–Schmidt (splitmix64 RNG, Box–Muller).

**Commands** (from repo root; outputs committed in `problems/mub-six/data/`):

```bash
cc -O2 -o problems/mub-six/explore/mub_search problems/mub-six/explore/mub_search.c -lm
S=problems/mub-six/explore/mub_search; D=problems/mub-six/data
for dk in "2 4 200" "3 4 200" "4 4 200" "5 4 200" "7 4 200" "6 3 200" "4 5 200" "6 4 2000"; do
  set -- $dk; $S $1 $2 $3 20261002 50000 $D/best_d$1k$2.json > $D/run_d$1k$2.txt; done
$S 6 4 20000 7 50000 $D/best_d6k4_20k.json > $D/run_d6k4_20k.txt
$S 7 4 2000 11 50000 $D/best_d7k4_2k.json > $D/run_d7k4_2k.txt
python harness/mub-six/mub_check.py $D/best_d6k4_20k.json
```

**Analytic control (d = 2, k = 4).** For qubit bases with Bloch axes n_a,
|⟨a_i|b_j⟩|² = (1 ± n_a·n_b)/2, so each pair contributes 4·(cos θ_ab / 2)² =
cos² θ_ab and L = Σ_{a<b} (n_a·n_b)². The frame-potential bound for N unit
vectors in R³, Σ_{a,b}(n_a·n_b)² ≥ N²/3, gives with N = 4:
L ≥ (16/3 − 4)/2 = **2/3**, attained by the four cube diagonals
(cos² θ = 1/9 for every pair). This is a proof, and it is the code's
correctness check. (For d = 6, k = 4 the same Welch-type bound gives a
negative number, i.e. nothing. That is one way to see why the problem is hard
below k = d + 2.)

**Results** (hit = final L < 1e-20):

| case | starts | seed | hits | best L | note |
|---|---|---|---|---|---|
| d=2, k=4 | 200 | 20261002 | 0 | 0.66666666666666641 | all 200 at 2/3 to 1e-15 |
| d=3, k=4 | 200 | 20261002 | 200 | 1.1e-30 | complete set |
| d=4, k=4 | 200 | 20261002 | 200 | 3.1e-29 | |
| d=4, k=5 | 200 | 20261002 | 200 | 2.4e-29 | complete set |
| d=5, k=4 | 200 | 20261002 | 126 | 2.7e-29 | |
| d=7, k=4 | 200 | 20261002 | 2 | 6.5e-28 | |
| d=7, k=4 | 2000 | 11 | 14 | — | hit rate 0.7% |
| d=6, k=3 | 200 | 20261002 | 153 | 7.0e-29 | |
| d=6, k=4 | 2000 | 20261002 | 0 | 0.051249218996283839 | |
| **d=6, k=4** | **20000** | **7** | **0** | **0.051249218996283818** | 14,138 starts (70.7%) at the floor |

**Endpoint census, d = 6, k = 4, 20,000 starts.** Rounded to 5 decimals,
endpoints fall on 19 distinct values. The most populated: 0.05125 (14,138),
0.23401 (2,074), 0.19565 (1,897), 0.17887 (881), 0.26195 (255), 0.27713
(203), 0.21223 (199), 0.29415 (150), 0.26768 (133). The rest are rarer, up to
L ≈ 0.43.

**Structure of the best configuration** (`best_d6k4_20k.json`, rechecked with
`mub_check.py`; orthonormality error at roundoff):
- the pairs (standard, B_b) have defect ~1e-18 for b = 2, 3, 4, so all three
  are (1/√6)·complex Hadamard to roundoff;
- the three pairs (B_a, B_b), a, b ≥ 2, each carry defect 0.017083 (L/3);
- in each of those three pairs, **every row** of the 6×6 matrix of
  |⟨a_i|b_j⟩|² is a permutation of
  {0.124398, 0.151402, 0.18105, 0.18105, 0.18105, 0.18105}.

**Comparison with the literature.** ASD = 1 − L/(C(4,2)·5) = 1 − L/30 =
**0.998291692700124**, against Raynal–Lü–Englert's ≈ 0.9983 [T]. Inverting
their closed form ASD = (71 − 12 cos⁴θ)/70 [T] at our L gives
sin²θ = 0.694552, against their 0.6946 [T]. Agreement holds to every digit
they publish.

## Outcome

`EVIDENCE`, scoped to the runs above: in 20,000 + 2,000 random starts of this
local search with the stated stopping rule, no four-basis configuration in
C^6 had L below 0.0512492189962838, and the most frequent endpoint (70.7% of
starts) equals the published optimum of Raynal–Lü–Englert to the 4–5 digits
they quote. The d = 2 floor of exactly 2/3 is `VERIFIED` (proved above;
reproduced numerically to 1e-15).

**Not claimed:**
- that four MU bases do not exist in C^6 (local search proves nothing about
  a narrow basin it never entered — see the d = 7 control);
- that 0.0512492189962838 is the *global* minimum of L (not certified, see
  leads);
- any novelty for the floor value or for the negative search outcome. Both
  are **rediscoveries** of Butterley–Hall (2007) and Raynal–Lü–Englert
  (2011). The three-value row pattern in the overlap matrices may well be
  implicit in RLE's explicit family. That was not checked, so it is not
  claimed as new.

## Why it failed / what survived

The search "failed" in the expected way: no zero. The specific obstruction
the landscape shows is a **dominant funnel to a non-zero floor**. Every start
that gets within ~0.18 of zero ends at the same configuration type. There,
three Hadamard bases are each unbiased to the standard basis, and the whole
residual defect is shared equally among the three Hadamard–Hadamard pairs.
So the failure is entirely in making three Hadamard matrices mutually
unbiased *with each other*, which matches the reduction in PROBLEM.md.

`SPECULATION` (labelled, not a bound): the d = 6 landscape is more funnelled
than d = 7 (70.7% of starts reach one floor vs. 0.7% hit rate at d = 7). If a
four-basis MU set in d = 6 had a basin even one-tenth the relative size of
d = 7's, about 14 hits would be expected in 20,000 starts. Basin measures are
not comparable across dimensions, so this is a heuristic only.

Reusable: the C search kernel (any d ≤ 8, k ≤ 9); the stdlib checker in the
harness; the calibrated control table, which turns any future "no hit in N
starts" claim into something interpretable; the d = 2 analytic floor as a
regression test.

## Leads generated

1. **Certify the floor as a local minimum.** Compute the Riemannian Hessian
   of L at `best_d6k4_20k.json` (expect a zero eigenspace from the symmetry
   group of phases and permutations). Check that the rest is positive
   definite, with interval arithmetic. Either outcome is definite.
2. **Closed forms for the three overlap values.** Test whether 0.124398,
   0.151402 and 0.18105 are algebraic in RLE's θ_opt (sin²θ_opt ≈ 0.694552),
   e.g. by PSLQ at 30 digits after Newton-polishing the configuration.
3. **Structured starts.** Seed with every known family of 6×6 complex
   Hadamard matrices (Fourier F(a,b) and its transpose, Diță D_6, Björck B_6,
   Tao S_6, Karlsson's K_6 family). Report whether any start goes below
   0.0512492. Outcome either way is a definite statement about those seeds.
4. **Tighten the control calibration.** Measure d = 7, k = 4 and d = 8, k = 4
   hit rates at 20,000 starts, so the d = 6 null can be compared at equal
   sample size.

## References

- P. Butterley, W. Hall, *Numerical evidence for the maximum number of
  mutually unbiased bases in dimension six*, Phys. Lett. A 369 (2007),
  arXiv:quant-ph/0701122.
- S. Brierley, S. Weigert, *Maximal sets of mutually unbiased quantum states
  in dimension six*, Phys. Rev. A 78 (2008), arXiv:0808.1614.
- P. Raynal, X. Lü, B.-G. Englert, *Mutually unbiased bases in dimension 6:
  The four most distant bases*, Phys. Rev. A 83 (2011), arXiv:1103.1025 [T].
- P. Jaming, M. Matolcsi, P. Móra, F. Szöllősi, M. Weiner, *A generalized
  Pauli problem and an infinite family of MUB-triplets in dimension 6*,
  J. Phys. A 42 (2009), arXiv:0902.0882.
- Welch bound / frame potential: L. Welch, IEEE Trans. Inf. Theory 20 (1974);
  J. Benedetto, M. Fickus, Adv. Comput. Math. 18 (2003).
