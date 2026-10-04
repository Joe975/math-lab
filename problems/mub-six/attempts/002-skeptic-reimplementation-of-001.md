# 002 — Skeptic re-implementation of 001's defect-landscape census

- **Problem:** mub-six, `problems/mub-six/PROBLEM.md`
- **Date:** 2026-10-02
- **Mode:** blind
  (The skeptic was a separate Claude subagent in the same session. It was
  given only the definition of L and the list of cases, and told not to read
  anything under `problems/mub-six/` or `harness/mub-six/`. It reports that it
  did not. It never saw 001's code, its numbers, or the literature values.
  The director, who did see 001, wrote this record and ran the cross-checks
  marked "director" below.)
- **Type:** skeptic review (independent re-implementation of a computation)
- **Reviews:** `attempts/001-multistart-defect-landscape.md`
- **Tools:** `problems/mub-six/explore/skeptic/mubexp.c` (C99; written from
  scratch by the skeptic), `problems/mub-six/explore/skeptic/verify.py`
  (the skeptic's stdlib checker), `harness/mub-six/mub_check.py` (director
  cross-check). The skeptic worked in a scratch directory. Its files were
  copied unchanged, so paths inside the logs refer to that scratch location.
- **Sources:** none beyond 001's. The skeptic used no literature.

## Claims attacked

From 001:

1. C1. In d = 6, k = 4, random-start local search finds no configuration with
   L below 0.0512492189962838.
2. C2. The dominant endpoint (about 70% of starts) is that floor value.
3. C3. ASD at the floor = 0.998291692700124.
4. C4. At the floor, one basis is unbiased to the other three. The residual is
   shared equally, L/3 each, by the other three pairs.
5. C5. In each of those three pairs, every row of |⟨a_i|b_j⟩|² is a
   permutation of {0.124398, 0.151402, 0.18105 ×4}.
6. C6. Controls: d = 2, k = 4 floor is exactly 2/3. Known MU sets are found
   in d = 3 and 5 (k = 4) and in d = 6 (k = 3).

## Refutations found

None. No claim of 001 failed under re-implementation.

One wording precision, which is not a correction. 001 states C4 with the
*standard* basis as the hub, because 001 fixes the standard basis. In the
skeptic's best configuration, which also fixes the standard basis as basis 1,
the hub is basis 2. The structure is the same up to relabelling which basis
is the hub. So C4 should be read as "some basis is unbiased to the other
three", not "the fixed basis is".

## Claims that survive

How the skeptic's method differs from 001's, so that agreement is not shared
error:

| | 001 | 002 (skeptic) |
|---|---|---|
| manifold step | projected gradient + Gram–Schmidt retraction | body-frame exponential map U ← U·exp(A), A anti-Hermitian (Taylor-20 with scaling and squaring) |
| optimiser | gradient descent, adaptive step | L-BFGS (memory 12) + Armijo backtracking |
| gradient | Euclidean, then projected | Lie-algebra analytic; checked against central differences, relative error ~1e-8 |
| stopping | ‖P‖ < 1e-14, η < 1e-18, or 50,000 iterations | ‖g‖ < 1e-12, L < 1e-28, or 30 stalled steps; 20,000 cap never hit |
| RNG | splitmix64, seeds 7 and 20261002 | xorshift64*, seeds 1–7 |

- **C1 and C2 survive.** Skeptic: d = 6, k = 4, 2000 starts (seed 6) and
  1000 starts (seed 7). Best L = 0.051249218996284 in both. 1431/2000 (71.6%)
  and 717/1000 (71.7%) of starts land within 1e-8 of it, and nothing goes
  below it. 001 had 70.7% of 20,000. The secondary minima also agree:
  0.23400568, 0.19565212, 0.17886734, 0.26194871, 0.21223226 appear in
  both censuses at similar relative frequencies. The attack this rules out:
  an optimiser-specific artefact, e.g. a retraction that cannot leave a
  manifold region, or a stopping rule that halts early. A different step
  map, optimiser and stopping rule land on the same values.
- **C3 survives.** The skeptic's configuration, rechecked by the director with
  `mub_check.py`, gives L = 0.0512492189962838 and ASD = 0.998291692700124,
  identical to 001 at printed precision.
- **C4 survives** (with the wording precision above). Director, on
  `data/skeptic/best_d6k4.json`: pairs (0,1), (1,2), (1,3) have defects
  4e-20 to 1e-19. Pairs (0,2), (0,3), (2,3) each have 1.708307e-02 = L/3.
- **C5 survives.** Director, same file: in all three defective pairs, every
  row's sorted overlaps are exactly [0.124398, 0.151402, 0.18105, 0.18105,
  0.18105, 0.18105] at 6 digits. The two configurations came from different
  code and were not aligned by hand, so the pattern is a property of the
  minimum, not of one run.
- **C6 survives.** Skeptic controls: d = 2 gives 200/200 at exactly 2/3;
  d = 3, k = 4 gives 200/200 to ~1e-26; d = 5, k = 4 gives 114/200 (001:
  126/200), with the same 0.353617013 secondary minimum; d = 6, k = 3 gives
  140/200 (001: 153/200).

**New observation from the skeptic, cross-checked by the director.**
0.0512492189963 is also a *local* minimum for d = 6 with only **three**
bases. The skeptic found it in 23 of 200 starts. The director recounted
001's own d = 6, k = 3 log: 20 of 200 starts at 0.05124922. Combined with C4,
this suggests a reading of the k = 4 optimum: a frustrated 3-basis local
minimum whose three bases all happen to be Hadamard with respect to one
common basis, so that adding that common basis costs nothing.
`SPECULATION`: that every frustrated 3-basis minimum at this value admits a
common unbiased basis was not tested; only the two k = 4 optima were
inspected.

**Not claimed:** this review verifies that two independent implementations
agree. Both are local searches from Haar-random starts. Neither certifies a
global minimum, and agreement between them does not address a narrow basin
that both samplers miss (001's d = 7 caveat stands).

## References

- `problems/mub-six/attempts/001-multistart-defect-landscape.md`
- Skeptic reproduction commands (run inside `problems/mub-six/explore/skeptic/`
  after `cc -O2 -o mubexp mubexp.c -lm`): `./mubexp 6 4 2000 6 20000
  best_d6k4.json`; `./mubexp 6 4 1000 7 20000`; controls `./mubexp 2 4 200 1
  20000`, `3 4 200 2`, `5 4 200 3`, `6 3 200 4`. Logs in
  `problems/mub-six/data/skeptic/`.
- Raynal–Lü–Englert, arXiv:1103.1025 (for context only; not used by the
  skeptic).
