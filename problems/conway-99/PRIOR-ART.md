# Conway's 99-graph — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None yet.** Added 2026-07-27; the first attempts are in flight. Until one
lands there is no prior art to be informed by, so `blind` and `informed` mode
are equivalent here.

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
