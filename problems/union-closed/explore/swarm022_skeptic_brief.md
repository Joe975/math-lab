# Swarm brief — attempt 022 skeptic verifier (cross-family re-implementation)

You are an expert numerical/exact-arithmetic programmer. Write ONE
complete, standalone Python 3 program (standard library only — `fractions`,
`json`, `math`, `sys`) that CERTIFIES signs of certain census functionals
in exact rational arithmetic with directed rounding. Your program is a
skeptic re-implementation: it must be written from THIS SPEC ALONE, and
its whole value is independence, so implement everything yourself (your
own certified log2, your own root enclosure). Output the complete program
in one Python code block, nothing else after it.

## Input

argv[1] is a path to a JSON file: {"witnesses": [W, ...]} where each W is
  {"name": str, "kind": "mm_abs" | "aggregate" | "mm_abs_rplus",
   "n": int, "t": [num, den],
   "u": {bitstring: [num, den], ...}}
n <= 5. t = num/den > 1 is the "tilt" (a rational). u maps an atom
(bitstring of length n, leftmost char = bit n-1, i.e. Python
int(s, 2)) to a positive rational weight num/den.

## The mathematical objects

Work with exact Fractions throughout.

Coupling: for atoms A, B (integers, bit j of A = coordinate j+1), define
  pi(A, B) = u(A) * u(B) * t^{popcount(A & B)}      (unnormalized)
and Z = sum of pi over all pairs. The normalized coupling is pi/Z. The
measure is mu(A) = sum_B pi(A, B) / Z (exact rational).

Element marginals: marg_j = sum_{A : bit j set} mu(A), j = 0..n-1.
"in_regime" means every marg_j < 38271/100000 (strict, exact).

Census: for each coordinate i = 1..n-1 (EXCLUDING i = n), prefix mask
m_i = (1 << (i-1)) - 1, coordinate bit c_i = 1 << (i-1). Group the
normalized coupling by history pairs (a, b) = (A & m_i, B & m_i): each
group is a 2x2 table
  p00 = sum of pi/Z with A bit_i = 0, B bit_i = 0
  p01 = A bit 0, B bit 1;  p10 = A bit 1, B bit 0;  p11 = both 1.
A table is NONDEGENERATE iff all four cells > 0. For each nondegenerate
table define (all exact rationals):
  mass = p00+p01+p10+p11
  x = (p00+p01)/mass, y = (p00+p10)/mass     (zero-margins)
  OR = (p00*p11)/(p01*p10)
  dev = log2(OR / t)                          (needs certified log2)

Plackett root: z_t(x, y) = the unique root of
  F(z) = z*(1-x-y+z) - t*(x-z)*(y-z)
in the open interval (x*y, min(x,y)) (this bracket is valid since t > 1;
verify F changes sign on it and fail loudly if not). Enclose it by
bisection in exact rationals to width <= 2^-90.

Weight (up to a global positive constant that cannot affect any sign):
  w(x, y) = | hp(z) * dz |   evaluated at z = z_t(x,y), where
  hp(z) = log2((1-z)/z)                       (needs certified log2)
  dz    = (x-z)*(y-z) / ( (1-x-y+2z) + t*(x+y-2z) )
Certify an interval [w_lo, w_hi] for w from the z-interval: hp is
strictly decreasing in z; the numerator of dz is decreasing in z on
z < min(x,y); the denominator is linear in z (slope 2-2t < 0) — check it
stays > 0 on the z-interval and fail loudly otherwise. Use outward
(directed) rounding at every step: lower bounds round down, upper bounds
round up. For the absolute value: if the signed interval straddles 0,
use [0, max(|lo|, |hi|)].

## Certified log2

Implement your own: for rational r > 0, first shift by the unique
integer k with r' = r / 2^k in [1, 2), then extract >= 100 binary
fraction digits of log2(r') by repeated squaring (square r', if >= 2
emit digit 1 and halve, else digit 0). CRITICAL: exact squaring 100
times makes numerators astronomically large — instead run the squaring
loop TWICE, once on a lower track and once on an upper track, where
after every squaring you round the value to a dyadic rational with a
200-bit mantissa (floor(v * 2^200)/2^200 on the lower track,
ceil(...)/2^200 on the upper track). Rounding down can only delay digit
emission and rounding up can only advance it, so the lower track's digit
stream gives a valid lower bound for log2(r') and the upper track's
stream plus 2^-D a valid upper bound (D = digits extracted). Return
[k + lo_value, k + up_value + 2^-D] as exact Fractions. All operations
on exact Fractions/integers, no floats anywhere.

## What to certify, per witness kind

- kind "aggregate": S = sum over nondegenerate tables of mass * dev.
  Report certified [S_lo, S_hi].
- kind "mm_abs": N = sum over nondegenerate tables of
  mass * w(x,y) * dev, with interval products done with correct sign
  handling (four-products min/max). Report certified [N_lo, N_hi]. Also
  report the certified aggregate [S_lo, S_hi] from the same tables.
- kind "mm_abs_rplus": as "mm_abs", PLUS certify for EVERY nondegenerate
  table that z_t(x, y) < 1/2 strictly (upper end of the z-enclosure
  < 1/2); report a boolean all_rplus.

Always also report: in_regime (exact), max_marginal (as a float for
display), n_tables (count of nondegenerate tables), and dich_ok = True
iff every table with positive mass is either fully positive or has a
zero row/column pattern consistent with a degenerate conditional margin
(simply: report the count of tables having some but not all cells zero
— call it mixed_zero_tables; it should be 0... actually a degenerate
table has a full zero row or column; count tables violating that).

## Self-tests (the program runs these first; abort if any fails)

1. Product anchor: n = 3, u(A) = (3/10)^{popcount(A)} * (7/10)^{3-popcount(A)},
   t = 3/2: every nondegenerate table must have OR = t EXACTLY (rational
   equality), hence every dev-enclosure must contain 0, and the certified
   aggregate and mm_abs numerator intervals must both contain 0.
2. log2 anchor: certified log2(2) = [1,1+eps], log2(8/1) ~ 3,
   log2(1/2) ~ -1, each within 2^-90.
3. Plackett anchor: x = y = 1/2, t = 9 gives z_t = 3/8 exactly
   (check: F(3/8) = 3/8*3/8 - 9*(1/8)^2 = 0). The bisection enclosure
   must contain 3/8 with width <= 2^-90.

## Output

For each witness print ONE line of JSON:
  {"name": ..., "num_lo": str(Fraction), "num_hi": str(Fraction),
   "agg_lo": ..., "agg_hi": ..., "in_regime": bool, "max_marginal": float,
   "all_rplus": bool or null, "n_tables": int, "mixed_zero_tables": int}
(for kind "aggregate", num_* = null). Use str() of the Fraction (exact),
plus a float rendering for readability. Print self-test results first,
one line each, prefixed "SELFTEST".

Correctness over speed; but avoid quadratic blowups in the log2 digit
extraction (keep the squared value's numerator/denominator bounded by
reducing... Fractions auto-reduce; 100 digits of squaring on a ratio r'
in [1,2) is fine). n <= 5 and <= 16 atoms per witness: brute force over
all pairs is fine.
