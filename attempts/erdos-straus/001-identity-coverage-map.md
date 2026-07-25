# 001 — Residue-class identity coverage map mod 840

**Problem:** Erdős–Straus conjecture (4/n = 1/x + 1/y + 1/z for all n ≥ 2).
**Date:** 2026-07-25.
**Tool:** `tools/es_coverage.py` (re-runnable; see its docstring).

## Approach

Build the classical identity-side baseline: implement the polynomial identity
families that solve 4/n on residue classes, machine-verify each one, compute
exactly which classes mod 840 = 8·3·5·7 the set leaves uncovered, and test
Mordell's quadratic-residue obstruction against that set. Independently
cross-check with a bounded brute-force solver on every prime p < 10^5, and
mine the representation counts f(p) for structure in low-count primes.

## What was done

### 1. Twelve identity families, all machine-verified

Each family below was verified **exactly** (integer arithmetic, assert
`4xyz = n(xy+yz+zx)`) for **every** qualifying n ≤ 2·10^5 plus 500 random
qualifying n up to 10^15 per family (`--verify`, seeded, reproducible).
Since each identity after clearing denominators is a polynomial identity of
degree ≤ 3 in the class parameter, agreement at > 4 points already proves the
rational identity on the whole class; integrality of every denominator is
enforced by exact `divmod` (asserted at each evaluation, so the congruence
conditions were themselves machine-checked at every tested n).

| class | identity (n in class) | small-n checks |
|---|---|---|
| n ≡ 0 (2) | n=2m: 1/m + 1/(m+1) + 1/(m(m+1)) | 100000 |
| n ≡ 0 (3) | n=3m: 1/m + 1/(4m) + 1/(12m) | 66666 |
| n ≡ 2 (3) | n=3k+2: 1/(k+1) + 1/n + 1/(n(k+1)) | 66667 |
| n ≡ 3 (4) | n=4k+3: 1/(k+1) + 1/(2(k+1)n) + 1/(2(k+1)n) | 50000 |
| n ≡ 5 (8) | n=8k+5: 1/(2(k+1)) + 1/((k+1)n) + 1/(2(k+1)n) | 25000 |
| n ≡ 0 (5) | n=5m: 1/(2m) + 1/(4m) + 1/(20m) | 40000 |
| n ≡ 0 (7) | n=7m: 1/(2m) + 1/(28m) + 1/(28m) | 28571 |
| n ≡ 17 (20) | n=20d−3: 1/(5d) + 1/(2dn) + 1/(10dn) | 10000 |
| n ≡ 33 (40) | n=40d−7: 1/(10d) + 1/(2dn) + 1/(5dn) | 5000 |
| n ≡ 41 (56) | n=56t−15: 1/(14t) + 1/(14tn) + 1/(nt) | 3571 |
| n ≡ 145 (168) | n=168t−23: 1/(42t) + 1/(21tn) + 1/(2tn) | 1190 |
| n ≡ 73 (84) | n=84d−11: 1/(21d) + 1/(2dn) + 1/(42dn) | 2381 |

The last five are Mordell-style: they handle the residues that are quadratic
**non**-residues mod 5 (classes 2, 3 mod 5) and mod 7 (classes 3, 5, 6 mod 7)
inside n ≡ 1 mod 4 / mod 8. All twelve are instances of the Elsholtz–Tao
"Type I" shape 4/n = 1/(abd) + 1/(acdn) + 1/(bcdn) with 4abcd = cn + a + b;
e.g. the n ≡ 73 (84) family is (a,b,c) = (21,1,2), found by the systematic
sweep described below (this is the family that catches classes 241 and 409
mod 840, which the textbook mod-20/40/56/168 families miss).

### 2. Coverage map mod 840 (`--coverage`)

Coverage progression (classes mod 840 still uncovered after adding each
family): 420 → 280 → 140 → 70 → 35 → 28 → 24 → 18 → 12 → 10 → 8 → **6**.

**Uncovered residue classes mod 840: {1, 121, 169, 289, 361, 529}.**

Mordell obstruction check — precisely what was checked:
- Computed S = {u² mod 840 : 0 ≤ u < 840, gcd(u,840) = 1} exhaustively.
  |S| = 6 and S = {1, 121, 169, 289, 361, 529} = {1², 11², 13², 17², 19², 23²}
  mod 840.
- Verified **uncovered set == S exactly** (both inclusions), and that every
  uncovered class r has gcd(r, 840) = 1. So the classes our identity set
  misses are precisely the squares of units mod 840, as Mordell's analysis
  predicts. (Every non-unit class mod 840 is fully covered via the
  divisibility families, since a common factor with 840 persists across the
  whole class.)
- Independent sweep: enumerated all Type I families (a ≤ b ≤ 210 with
  ab | 210 — necessary for a family to cover *full* classes mod 840 — and
  all c | a+b; class covered iff r ≡ −(a+b)/c mod 4ab). The sweep's total
  coverage equals our 12 families' coverage exactly: it covers nothing we
  miss, and **no Type I family covers any of the 6 square classes** even at
  these exhaustive parameter bounds. Computational confirmation of the
  known theorem that finitely many polynomial identities of this divisor
  type cannot finish the problem.

### 3. Brute-force cross-check on primes (`--existence`, `--counts`)

Bounded standard enumeration (no unbounded search): for each prime p, x runs
over (p/4, 3p/4]; with A = 4x−p, B = px, solutions of 1/y + 1/z = A/B
correspond to divisors d | B², d ≤ B, A | B+d (y = (B+d)/A, z = (B+B²/d)/A;
y ≥ x enforced, so each unordered triple x ≤ y ≤ z is counted once; since
p ∤ x, divisors of B² are p^i·d₀ with d₀ | x², i ≤ 1 after the d ≤ B cut).

- **Existence:** every one of the 9592 primes p < 10^5 has a solution
  (0 failures; each found solution re-verified exactly). Pure Python,
  independent code path from the counter.
- **Counts:** exact f(p) = #{(x,y,z) : x ≤ y ≤ z} for all 9592 primes via a
  C helper compiled on the fly (~90 s); cross-checked against an independent
  pure-Python counter for all p < 3000 — **all match**. No prime has
  f(p) = 0, consistent with the existence pass. Raw counts reproducible via
  `--dump`.

### 4. Representation-count statistics

Bottom 20 primes by f(p) among p > 1000 (all ≡ 1 mod 24):

| f | p | p%120 | p%840 | QR class? |
|---|---|---|---|---|
| 9 | 2521 | 1 | 1 | yes |
| 11 | 1201 | 1 | 361 | yes |
| 15 | 3361 | 1 | 1 | yes |
| 19 | 1009 | 49 | 169 | yes |
| 19 | 1129 | 49 | 289 | yes |
| 19 | 1249 | 49 | 409 | no |
| 19 | 3169 | 49 | 649 | no |
| 20 | 4201 | 1 | 1 | yes |
| 21 | 1873 | 73 | 193 | no |
| 22 | 1153 | 73 | 313 | no |
| 22 | 1321 | 1 | 481 | no |
| 22 | 1489 | 49 | 649 | no |
| 23 | 1801 | 1 | 121 | yes |
| 23 | 2689 | 49 | 169 | yes |
| 23 | 9601 | 1 | 361 | yes |
| 23 | 20521 | 1 | 361 | yes |
| 24 | 3529 | 49 | 169 | yes |
| 25 | 2161 | 1 | 481 | no |
| 25 | 3049 | 49 | 529 | yes |
| 25 | 5281 | 1 | 241 | no |

(The unrestricted bottom-20 is the same list diluted by tiny primes 2, 3, 5,
7, ..., plus p = 193, 241, 73, 97 — again all the non-tiny ones ≡ 1 mod 24.)

Structure found:
- **The bottom 24 low-f primes > 1000 are all ≡ 1 mod 24** (first exception
  at rank 25: p = 1409 ≡ 17 mod 24, f = 27). In the bottom 50: 47/50 are
  ≡ 1 mod 24; mod 120 they concentrate on 1 (×23), 49 (×15), 73 (×6),
  97 (×3).
- **QR-mod-840 membership enriches but does not characterize:** 22/50 of the
  bottom-50 lie in the six uncovered classes vs a baseline of 273/9592
  (2.8%) of all primes — a ~15× enrichment — but classes like 601, 649, 409,
  193, 313, 481, 241 mod 840 (≡ 1 mod 24, **not** squares mod 840) also
  produce very low counts; class 601 alone holds 6 of the bottom 50.
- Class averages over primes in [10^4, 10^5): mean f = 90.1 for p ≡ 1 mod 24
  (min 23) rising monotonically-ish to mean 392.9 for p ≡ 23 mod 24
  (min 119). Max observed f = 1621 at p = 87359 ≡ 23 mod 24.
- The 273 primes < 10^5 in the six obstructed classes all have f ≥ 1
  (lowest: 2521 → 9, 1201 → 11, 3361 → 15, 1009 → 19, 1129 → 19).
- Side observation from the existence pass: primes whose *first* solution
  (smallest x) sits farthest above p/4 are all in the obstructed classes
  (99961 ≡ 1, 87481 ≡ 121, 67369 ≡ 169, 61681 ≡ 361 mod 840), with the
  forced z astronomically large (up to ~7·10^13) — "no easy small-x
  solution" is another face of the same obstruction.

## Outcome

- Re-runnable tool `tools/es_coverage.py`; all four phases pass.
- **Result (computational, exact):** the 12-family identity set covers 834
  of 840 residue classes; the uncovered classes are exactly
  {1, 121, 169, 289, 361, 529} mod 840 = the squares of units mod 840,
  replicating Mordell's classical coverage bound with machine-verified
  identities.
- **Evidence (not proof):** every prime p < 10^5 has a solution; f(p) ≥ 9
  for p > 1000; low-f primes are governed first by p ≡ 1 mod 24, then by
  finer quadratic-residue structure mod 120 / 840.

## Why it failed / what survived

This was a baseline/mapping attempt, not a proof attempt. The known hard
wall was confirmed computationally rather than evaded: no Type I polynomial
identity (exhaustive sweep, a·b ≤ 210, all valid c) covers any square class,
so finitely many identities of this type provably cannot finish the
conjecture — any future identity-flavored attack in this lab must state how
it escapes the quadratic-residue obstruction (e.g. identities whose
applicability depends on the factorization of n or on representability by
quadratic forms, not on n mod m alone).

Survived: the verified identity table, the exact coverage map, the C-backed
exact counter f(p) for all primes < 10^5 (90 s), and the low-f prime
dataset.

## Leads generated

1. **f(p) is stratified by more than the obstruction classes.** All of the
   lowest-f primes are ≡ 1 mod 24, but non-square classes mod 840 (601, 649,
   409, 313, 193, 481, 241) appear prominently among them. Hypothesis worth
   testing: f(p) correlates with *how many* identity families / small Type I
   parameter triples apply to p's class (classes covered by exactly one thin
   family, like 601 ≡ 41 mod 56, sit low). Quantify: regress f(p) against
   the count of (a,b,c) with 4abc | cp+a+b for small a,b,c.
2. **Class 601 mod 840 cluster** (6 of bottom 50): not a square mod 840
   (601 ≡ 6 mod 7), so not obstruction-explained. Deserves its own look —
   possibly a mod-7 secondary effect (6 is −1 mod 7).
3. **Minimal-f primes as test cases:** 2521, 1201, 3361, 4201, 9601, 20521
   (f = 9–23) are natural stress tests for any proposed constructive scheme;
   their full solution lists are cheap to enumerate with the existing tool.
4. **Extend the counter** to p up to 10^6–10^7 (the C kernel is O(p·d(x²))
   per prime; ~hours at 10^6) to test whether min f(p) among p ≡ 1 mod 24
   grows, and at what rate — Elsholtz–Tao average bounds give (log p)³-type
   growth on average; the *minimum's* growth is the interesting open
   statistic, and a slow-growing minimum sequence would pinpoint the
   genuinely hard prime families.
5. **First-solution offset statistic** (smallest x minus p/4): its extremes
   coincide with obstructed classes; as a cheap proxy for "hardness" it
   could screen much larger p than full counting can reach.

## Reproduction

```
python3 tools/es_coverage.py --all                      # ~2.5 min total
python3 tools/es_coverage.py --verify --coverage        # identities + map
python3 tools/es_coverage.py --existence --counts \
    --pmax 100000 --pycheck 3000 [--dump counts.txt]    # ~95 s (needs cc)
```
