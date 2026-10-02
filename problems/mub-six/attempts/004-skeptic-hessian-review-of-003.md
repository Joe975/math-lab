# 004 — Skeptic review of 003's second-order claim

- **Problem:** mub-six, `problems/mub-six/PROBLEM.md`
- **Date:** 2026-10-02
- **Mode:** blind
  (The skeptic was a separate Claude subagent. It was told not to read
  `problems/mub-six/explore/` or `problems/mub-six/attempts/`, and reports
  that it did not. It worked from its own 002 configuration and its own
  code. The director wrote this record and ran the cross-comparison marked
  "director".)
- **Type:** skeptic review (independent re-implementation of a computation)
- **Reviews:** `attempts/003-floor-hessian-and-closed-form.md`, claim 1
  (second order) and the overlap values of claim 2
- **Tools:** `problems/mub-six/explore/skeptic2/hess.c` (the skeptic's code:
  C, `long double`, hyperdual forward-mode automatic differentiation for an
  exact Hessian in a Cayley chart, plus a Richardson-extrapolated
  second-difference Hessian of L as an internal check; ~5.6 s). Output in
  `data/skeptic2/hess_out.txt`, polished point in
  `data/skeptic2/polished_d6k4.json` (wrapped with `d` for the harness).
  Director: `explore/hessian.py` functions, output in
  `data/hessian_d6k4_frobenius.txt`.
- **Sources:** none beyond 003's.

## Claims attacked

1. C1. At the four-basis floor the gradient can be driven to roundoff, so
   the point is a critical point.
2. C2. The Hessian has exactly 23 zero modes, and they are exactly the gauge
   directions (column phases + common left diagonal unitary).
3. C3. On the complement the Hessian is positive definite, so the point is a
   nondegenerate local minimum modulo gauge.
4. C4. The quoted spectrum: min 0.21488, max 10.263 "in the coordinates
   above".
5. C5. The distinct overlap values and multiplicities at the floor:
   0.1243979456657… ×18, 0.1514019774774… ×18, 1/6 ×108, 0.1810500192142… ×72.

## Refutations found

### C4 is coordinate-dependent in a way 003 does not spell out

003's eigenvalues are correct for 003's coordinates, but those coordinates
depend on two choices that 003 does not flag as affecting the numbers:

- **Generator normalisation.** 003's off-diagonal generators have Frobenius
  norm √2. That is a non-uniform rescaling, not a constant factor.
- **Which basis is fixed.** 003's configuration fixes the *hub* basis, the one
  unbiased to all others. The skeptic's configuration fixes a non-hub basis.
  Fixing a different basis changes the induced metric on the quotient, so it
  changes the eigenvalues.

Director's check (`data/hessian_d6k4_frobenius.txt`):
- 003's Hessian rescaled to Frobenius-orthonormal generators, on 003's point
  (hub fixed): min **0.10744**, max **5.1316**.
- The same code on the skeptic's point (non-hub fixed): min **0.0916047112**,
  max **4.722382207**. This matches the skeptic's exact AD values
  (9.1604711327e-02, 4.7223822083) to ~9 digits, and the full 85-value list
  matches entry by entry at 4 significant figures.

**Corrected statement of C4:** the *signature* of the Hessian (23 zero, 85
positive, 0 negative) is invariant. The eigenvalue *magnitudes* depend on the
coordinate normalisation and on which basis is held fixed. 003's
"min 0.215" is one representative. Under Frobenius normalisation the minimum
is 0.107 with the hub fixed and 0.0916 with a non-hub basis fixed. The
multiplicity pattern differs between the two gauge fixings too: 4-folds with
the hub fixed, doublets otherwise. So multiplicities are not invariants of
the minimum either. Nothing in 003's conclusions depends on the magnitudes.

## Claims that survive

- **C1 survives.** The skeptic polished in `long double` by Newton with
  pseudo-inverse and the true Cayley retraction: |g| went 2.0e-9 → 1.0e-18
  in one step, at L = 0.05124921899628382781. That matches 003's 50-digit
  L* = 0.0512492189962838277828… to all 20 digits shown. The attack: a
  different chart (Cayley, not exponential), different precision, and a
  different starting configuration (the skeptic's own 002 point, never
  aligned with 003's).
- **C2 survives.** The skeptic derived the symmetry count independently
  (18 + 6 − 1 = 23), built the 23 generators, and found them rank 23 with
  max |Hs|/|s| = 3.7e-19. Exactly 23 eigenvalues have |λ| < 1e-8 (largest
  2.5e-19), and the next is 0.0916, a gap of about 17 orders of magnitude.
  Every null eigenvector projects onto the symmetry span with squared norm
  1.000000000000000. The attack: an exact (AD) Hessian rather than
  finite differences, so a near-zero-but-nonzero curvature could not hide
  in FD noise.
- **C3 survives.** With the symmetry span projected out, all 85 remaining
  eigenvalues are positive. The skeptic's Richardson-extrapolated
  second-difference Hessian, a third method, agrees with its AD Hessian to
  3.0e-11 in eigenvalues.
- **C5 survives.** The skeptic's values at 18 digits:
  0.124397945665703215 ×18, 0.151401977477478299 ×18, 1/6 ×108,
  0.181050019214204622 ×72. They agree with 003's 50-digit values in every
  digit shown. The skeptic did not attempt closed forms, so 003's closed
  forms are checked only by 003's 50-digit test and RLE's published radical.

**Not addressed by this review:** 003 claim 3 (the k = 3 valley) and the
closed-form algebra (exact Fraction arithmetic, re-runnable by anyone).
Neither was independently re-implemented. The valley claim remains 003's
alone.

**Methodological lesson** for any Hessian-at-a-quotient claim in this lab:
report the signature, which is invariant, as the claim, and quote magnitudes
only together with the gauge fixing and the generator normalisation. Two
correct implementations here disagreed in every eigenvalue until both were
pinned down.

## References

- `problems/mub-six/attempts/003-floor-hessian-and-closed-form.md`
- Skeptic reproduction: `cd problems/mub-six/explore/skeptic2 && cc -O2 -o
  hess hess.c -lm && ./hess ../../data/skeptic/best_d6k4.json polished.json`
  (the input is the skeptic's 002 configuration in its original bare-list
  format; the copy in `data/skeptic/` is wrapped with `d` and may need
  unwrapping for `hess`).
