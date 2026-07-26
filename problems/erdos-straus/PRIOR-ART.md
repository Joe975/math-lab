# Erdős–Straus — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`. Full records: `attempts/`.
Route-specific tooling: `explore/`.

## Editorial view of the attack surface

- Extend/complete the table of polynomial identities by residue class; map
  exactly which classes remain uncovered and why.
- Study the representation count f(n) statistically; look for structure in
  primes with few representations.
- Find primes with minimal representation counts and mine them for common
  structure.

## Attempts

### 001 — Residue-class identity coverage map mod 840 · `VERIFIED` + `EVIDENCE`

Twelve polynomial identity families, each machine-verified exactly over its
whole qualifying range up to 2·10^5 plus 500 random qualifying n up to 10^15.
They cover 834 of 840 residue classes; the uncovered set is exactly
{1, 121, 169, 289, 361, 529} mod 840 — the squares of units mod 840,
reproducing Mordell's classical coverage bound with machine-verified
identities.

`EVIDENCE`: every prime p < 10^5 has a solution (0 failures), with exact f(p)
computed for each. f(p) ≥ 9 for p > 1000.

**Observation that did not survive** (see 002): low-f primes concentrate in
QR-related classes {1, 49, 73, 97} mod 120, but non-QR class **601 mod 840**
held 6 of the bottom 50 — flagged as unexplained.

### 002 — Is the class-601 anomaly real? · `REFUTED` (the anomaly) + `VERIFIED` (the real signal)

Exact f(p) for all 9,732 primes ≡ 1 mod 24 up to 10^6, kernel triple-validated.

**The 601 anomaly is noise.** Under size-normalized selection it holds 0 of
the bottom 49 (expected 2.07) — if anything under-represented. 001's raw
selection was dominated by the smallest primes; the direction replicates but
at p = 0.13 across 24 classes tested, which is unremarkable.

**Methodological lesson, generalizable.** f(p) grows with p, so a "bottom by
raw f" set re-measures small primes. Rank within dyadic bands of p instead.
This flaw produced 001's false positive.

**The genuine signal.** The size-normalized bottom 2% is 98% inside Mordell's
six QR classes mod 840 (which hold 24.4% share) — p ≈ 10^−110. Mechanism:
low-f primes are *identity-poor*, with class mean f monotone in the number of
covering Type I families (r = 0.92); smooth p−1 depresses f further.
Band-minima of f grow like (log p)³.

## Open lines

- Prove the identity-poverty mechanism: why does QR-class membership mod 840
  force fewer Type I covering congruences? Target a theorem
  `f(p) ≥ g(N_typeI(p))`, or a disproof.
