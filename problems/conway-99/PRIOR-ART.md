# Conway's 99-graph — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

- **001** (blind, 2026-07-27, `MAP`) — exact feasibility (corrects PROBLEM.md's
  multiplicities to 54/44), the forced pair model and partner-regular split,
  lemmas L1–L4 on automorphisms (unreviewed), and the finding that the
  propagating search fails its own positive control.
- **002** (blind, same cycle, `MAP`) — orbit matrices under assumed symmetry:
  no *semiregular* order-33 automorphism (`VERIFIED`, two derivations); L5
  (order-3 automorphisms fix 0 or 3 points, unreviewed); four construction
  engines, all topping out below n = 45 on known graphs.
- **003** (informed, 2026-08-25, `LIVE`) — the order-7 profile (one fixed
  vertex, fourteen 7-orbits) admits **no orbit matrix**: exhaustive in 26 s
  once the enumerator prunes by row inner products and breaks symmetry
  correctly; the fixed-point-free order-11 profile likewise (0.2 s). Hence
  no automorphism of order 7 or 11, resting on 001's L3/L2 and awaiting an
  independent re-run. Order 9 has exactly one orbit matrix, R = 3I + J. Also: the character-block multiplicity
  lemma (54 = a + 6r) and a calibrated local-model worker family.

## Editorial view of the attack surface

Why this problem was added: a single sharply-posed finite existence question
where every standard obstruction has already been checked and passed, so the
problem is genuinely open rather than merely unexamined. It brings machinery
the library does not have — eigenvalue interlacing, orbit matrices,
automorphism-order elimination — while reusing the isomorph-rejection
discipline built for Erdős–Gyárfás.

**Calibrate the expectation.** Full exhaustion is astronomically infeasible,
and the published SAT attack (Keramatipour 2026) reports that solvers do not
settle it. The realistic output of a cycle here is a `MAP` or a partial
elimination with a stated case split, plus reusable tooling. An attempt
reporting that it has found the graph, or proved it does not exist, has almost
certainly made an error.

This problem should carry a **low budget**, like Collatz, for that reason.

Concrete lines, if you want them:

- Build the strongly-regular checker and exact feasibility calculator, and
  validate them by recovering published verdicts across a battery of parameter
  sets — the calibration step that makes any later claim meaningful.
- Re-derive the forced local structure (every neighbourhood is 7K₂) and
  enumerate how far a canonical extension search gets before it explodes.
  *Where* it explodes, stated precisely, is a publishable obstruction.
- The automorphism results (Makhnev–Minakova, Behbahani–Lam, Cesarz–Woldar)
  narrow Aut to a handful of tiny groups. Independently re-deriving one of
  those eliminations would be a real verification contribution, since almost
  nobody has checked them.
