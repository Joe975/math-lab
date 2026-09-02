# Adjudication: our laminate vs the transcribed Cherkaev 2009 B2 bound

**STATUS: IN PROGRESS.** This file is a placeholder so the pointer from
`STATUS.md` resolves; it will be overwritten by the adjudication report.

## The question

At sigma = (1,2,5), m2 = 1/4 (so sqrt(m2) = 1/2 is rational and everything is
exact in Q), laminates built in this repo sit strictly BELOW the transcribed
Cherkaev 2009 Theorem 7.1 "B2" lower bound at m1 just under m11 = 1/8:

| m1 | our structure | HS_lo | transcribed B2 | ours - B2 |
|---|---|---|---|---|
| 31/250 | 3.0053470065 | 3.0053404539 | 3.0053523724 | -5.4e-06 |
| 3/25 | 3.0270098722 | 3.0268456376 | 3.0271504085 | -1.4e-04 |
| 11/100 | 3.0831692377 | 3.0816326531 | 3.0845383760 | -1.4e-03 |

The structures carry exact rational volume fractions hitting the target, an
isotropy residual below 1e-61, agreement between both harness algebra routes,
a passing Keller-Dykhne identity, and a passing `verify_laminate.py` run on
the saved record `attained/below-m11-f1-0.12.json`.

A structure below a valid lower bound is impossible, so something is wrong.
Per the verification contract in `PROBLEM.md`, the leading hypothesis is a
TRANSCRIPTION ERROR: the bound is marked [T] throughout and has never been
re-derived in this repo.

Complication worth stating up front: the transcribed B2 passes its most
natural internal consistency test. At the region boundary m1 = m11 it equals
HS_lo exactly (difference 0 in Q), so the piecewise bound is continuous there,
which a badly garbled formula would usually fail.

## What the adjudication must settle

1. Independent re-transcription of Theorem 7.1 and eqs 7.5-7.24.
2. Re-derivation of the underlying one-parameter translation family, maximized
   numerically at the three points independently of the closed-form split.
3. A from-scratch audit of our own structures, by an implementation sharing no
   code with the harness.
4. Whether a hypothesis of the bound (phase ordering, well-ordering, k3 = inf,
   exact isotropy, normalization) excludes our instance.
5. A verdict: transcription error / our structures invalid / scope mismatch /
   genuine discrepancy, with the evidence required for the last.

**Until this is settled, nothing is claimed below m11 in either direction.**
