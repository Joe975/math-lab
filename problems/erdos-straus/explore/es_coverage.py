#!/usr/bin/env python3
"""Erdos-Straus residue-class identity coverage map.

Conjecture: for every n >= 2 there are positive integers x, y, z with
    4/n = 1/x + 1/y + 1/z.

This tool:

  1. (--verify)   Implements 12 classical polynomial identity families that
                  solve 4/n on explicit residue classes (n even, n ≡ 0/2 mod 3,
                  n ≡ 3 mod 4, n ≡ 5 mod 8, n ≡ 0 mod 5/7, and mod-20/40/56/
                  84/168 families in Mordell's style).  Every family is
                  machine-verified exactly (integer arithmetic, no floats)
                  for ALL qualifying n up to --nmax and for random huge n up
                  to 10^15.  Since each identity, after clearing denominators,
                  is a polynomial identity in the class parameter of degree
                  <= 3, agreement at thousands of points proves the rational
                  identity itself; integrality of the denominators is forced
                  by the exact divmod construction (asserted at every use).

  2. (--coverage) Computes exactly which residue classes mod 840 = 8*3*5*7
                  remain uncovered by the implemented identity set, shows the
                  coverage progression as families are added, and checks
                  Mordell's quadratic-residue obstruction: the uncovered
                  classes are compared against {u^2 mod 840 : gcd(u,840)=1}.
                  Also runs an independent systematic sweep over Elsholtz-Tao
                  "Type I" parametrized identities
                      4/n = 1/(abd) + 1/(acdn) + 1/(bcdn),  4abcd = cn + a + b
                  (which cover the class n ≡ -(a+b)/c mod 4ab whenever
                  c | a+b and 4ab | 840) to confirm the hand-picked family
                  list attains everything this identity type can attain.

  3. (--existence) Brute-force cross-check: verifies every prime p < --pmax
                  (default 10^5) actually has a solution, via the standard
                  bounded enumeration: x runs over (p/4, 3p/4]; for each x,
                  with A = 4x - p and B = px, solutions of 1/y + 1/z =
                  A/B correspond to divisors d | B^2, d <= B, A | B + d
                  (then y = (B+d)/A, z = (B + B^2/d)/A).  No unbounded search.

  4. (--counts)   Exact representation counts f(p) = #{x<=y<=z} for all
                  primes p < --pmax via a small C helper (compiled on the fly;
                  pure-Python fallback/cross-check with --pycheck).  Reports
                  the bottom-20 primes by f(p) and their residue structure
                  mod 24 / 120 / 840.

Usage:
  python3 es_coverage.py --all              # everything, default bounds
  python3 es_coverage.py --verify --coverage
  python3 es_coverage.py --counts --pmax 100000 --pycheck 3000 [--dump f.txt]

Everything is exact integer arithmetic; randomness only chooses extra test
points for identity verification (seeded, reproducible).
"""

import argparse
import math
import os
import random
import subprocess
import sys
import tempfile

MOD = 840  # 8 * 3 * 5 * 7


# --------------------------------------------------------------------------
# 1. Identity families
# --------------------------------------------------------------------------

def exact_div(a, b):
    q, r = divmod(a, b)
    assert r == 0, f"exact_div: {b} does not divide {a}"
    return q


def s_even(n):       # n = 2m:  4/(2m) = 2/m = 1/m + 1/(m+1) + 1/(m(m+1))
    m = exact_div(n, 2)
    return (m, m + 1, m * (m + 1))


def s_div3(n):       # n = 3m:  4/(3m) = 1/m + 1/(4m) + 1/(12m)
    m = exact_div(n, 3)
    return (m, 4 * m, 12 * m)


def s_div5(n):       # n = 5m:  4/(5m) = 1/(2m) + 1/(4m) + 1/(20m)
    m = exact_div(n, 5)
    return (2 * m, 4 * m, 20 * m)


def s_div7(n):       # n = 7m:  4/(7m) = 1/(2m) + 1/(28m) + 1/(28m)
    m = exact_div(n, 7)
    return (2 * m, 28 * m, 28 * m)


def s_2mod3(n):      # n = 3k+2: 4/n = 1/(k+1) + 1/n + 1/(n(k+1))
    k = exact_div(n - 2, 3)
    return (k + 1, n, n * (k + 1))


def s_3mod4(n):      # n = 4k+3: 4/n = 1/(k+1) + 1/(2(k+1)n) + 1/(2(k+1)n)
    k = exact_div(n - 3, 4)
    return (k + 1, 2 * (k + 1) * n, 2 * (k + 1) * n)


def s_5mod8(n):      # n = 8k+5: 4/n = 1/(2(k+1)) + 1/((k+1)n) + 1/(2(k+1)n)
    k = exact_div(n - 5, 8)
    return (2 * (k + 1), (k + 1) * n, 2 * (k + 1) * n)


def s_17mod20(n):    # n = 20d-3: 4/n = 1/(5d) + 1/(2dn) + 1/(10dn)
    d = exact_div(n + 3, 20)
    return (5 * d, 2 * d * n, 10 * d * n)


def s_33mod40(n):    # n = 40d-7: 4/n = 1/(10d) + 1/(2dn) + 1/(5dn)
    d = exact_div(n + 7, 40)
    return (10 * d, 2 * d * n, 5 * d * n)


def s_41mod56(n):    # n = 56t-15: 4/n = 1/(14t) + 1/(14tn) + 1/(nt)
    t = exact_div(n + 15, 56)
    return (14 * t, 14 * t * n, n * t)


def s_145mod168(n):  # n = 168t-23: 4/n = 1/(42t) + 1/(21tn) + 1/(2tn)
    t = exact_div(n + 23, 168)
    return (42 * t, 21 * t * n, 2 * t * n)


def s_73mod84(n):    # n = 84d-11: 4/n = 1/(21d) + 1/(2dn) + 1/(42dn)
    d = exact_div(n + 11, 84)
    return (21 * d, 2 * d * n, 42 * d * n)


# (name, modulus m (divides 840), allowed residues mod m, solver, identity)
FAMILIES = [
    ("even",      2,  {0},   s_even,
     "n=2m: 1/m + 1/(m+1) + 1/(m(m+1))"),
    ("0 mod 3",   3,  {0},   s_div3,
     "n=3m: 1/m + 1/(4m) + 1/(12m)"),
    ("2 mod 3",   3,  {2},   s_2mod3,
     "n=3k+2: 1/(k+1) + 1/n + 1/(n(k+1))"),
    ("3 mod 4",   4,  {3},   s_3mod4,
     "n=4k+3: 1/(k+1) + 1/(2(k+1)n) + 1/(2(k+1)n)"),
    ("5 mod 8",   8,  {5},   s_5mod8,
     "n=8k+5: 1/(2(k+1)) + 1/((k+1)n) + 1/(2(k+1)n)"),
    ("0 mod 5",   5,  {0},   s_div5,
     "n=5m: 1/(2m) + 1/(4m) + 1/(20m)"),
    ("0 mod 7",   7,  {0},   s_div7,
     "n=7m: 1/(2m) + 1/(28m) + 1/(28m)"),
    ("17 mod 20", 20, {17},  s_17mod20,
     "n=20d-3: 1/(5d) + 1/(2dn) + 1/(10dn)   [n≡2 mod 5, n≡1 mod 4]"),
    ("33 mod 40", 40, {33},  s_33mod40,
     "n=40d-7: 1/(10d) + 1/(2dn) + 1/(5dn)   [n≡3 mod 5, n≡1 mod 8]"),
    ("41 mod 56", 56, {41},  s_41mod56,
     "n=56t-15: 1/(14t) + 1/(14tn) + 1/(nt)  [n≡6 mod 7, n≡1 mod 8]"),
    ("145 mod 168", 168, {145}, s_145mod168,
     "n=168t-23: 1/(42t) + 1/(21tn) + 1/(2tn) [n≡5 mod 7, 1 mod 8, 1 mod 3]"),
    ("73 mod 84", 84, {73}, s_73mod84,
     "n=84d-11: 1/(21d) + 1/(2dn) + 1/(42dn) [n≡3 mod 7, 1 mod 4, 1 mod 3]"),
]


def family_applies(fam, n):
    _, m, residues, _, _ = fam
    return n % m in residues


def check_solution(n, sol):
    x, y, z = sol
    assert x > 0 and y > 0 and z > 0
    assert isinstance(x, int) and isinstance(y, int) and isinstance(z, int)
    # exact: 4/n == 1/x+1/y+1/z  <=>  4xyz == n(xy+yz+zx)
    assert 4 * x * y * z == n * (x * y + y * z + z * x), (n, sol)


def verify_families(nmax, nrandom=500, seed=42):
    print(f"== Identity verification (exact, all qualifying n in [2, {nmax}]"
          f" + {nrandom} random n up to 10^15 per family) ==")
    rng = random.Random(seed)
    for fam in FAMILIES:
        name, m, residues, solver, ident = fam
        cnt = 0
        for n in range(2, nmax + 1):
            if n % m in residues:
                check_solution(n, solver(n))
                cnt += 1
        big = 0
        for _ in range(nrandom):
            r = rng.choice(sorted(residues))
            n = r + m * rng.randrange(1, 10**15 // m)
            if n >= 2:
                check_solution(n, solver(n))
                big += 1
        print(f"  [OK] {name:12s} verified on {cnt} small + {big} large n : "
              f"{ident}")
    print("  All identities machine-verified exactly (integer arithmetic).")
    print("  Note: each cleared-denominator identity is a polynomial identity"
          " of degree <= 3 in the class parameter, so agreement at >4 points"
          " already proves it for the whole class; integrality of every"
          " denominator is asserted via exact divmod at each evaluation.")


# --------------------------------------------------------------------------
# 2. Coverage map mod 840 and Mordell's obstruction
# --------------------------------------------------------------------------

def covered_residues(families):
    cov = set()
    for r in range(MOD):
        if any(family_applies(f, r if r >= 2 else r + MOD) for f in families):
            cov.add(r)
    return cov


def coverage_report():
    print("== Coverage map mod 840 ==")
    # progression
    prog = []
    for i in range(1, len(FAMILIES) + 1):
        cov = covered_residues(FAMILIES[:i])
        prog.append((FAMILIES[i - 1][0], MOD - len(cov)))
    print("  Coverage progression (family added -> classes still uncovered):")
    for name, left in prog:
        print(f"    + {name:12s} -> {left:3d} uncovered")

    cov = covered_residues(FAMILIES)
    unc = sorted(set(range(MOD)) - cov)
    print(f"  Final uncovered residue classes mod 840 ({len(unc)}): {unc}")

    # Mordell obstruction check
    units = [u for u in range(MOD) if math.gcd(u, MOD) == 1]
    sq_units = sorted({u * u % MOD for u in units})
    print(f"  Squares of units mod 840 ({len(sq_units)}): {sq_units}")
    print(f"  Uncovered == squares-of-units mod 840 : {unc == sq_units}")
    small_sq = sorted({u * u % MOD for u in [1, 11, 13, 17, 19, 23]})
    print(f"  (= {{1,11^2,13^2,17^2,19^2,23^2}} mod 840 : {unc == small_sq})")
    for r in unc:
        assert math.gcd(r, MOD) == 1
    print("  Every uncovered class is a unit class (gcd(r,840)=1): True")

    # Independent sweep over Type I identities 4abcd = cn + a + b
    # x=abd, y=acdn, z=bcdn.  Family (a,b,c) covers the full class
    # r mod 840 iff ab | 210 (so 4ab | 840), c | a+b, and
    # r == -(a+b)/c mod 4ab.
    sweep_cov = set()
    witnesses = {}
    for a in range(1, 211):
        for b in range(a, 211):
            if 210 % (a * b) != 0:
                continue
            m4 = 4 * a * b
            s = a + b
            for c in range(1, s + 1):
                if s % c:
                    continue
                r0 = (-(s // c)) % m4
                for r in range(r0 % m4, MOD, m4):
                    if r not in sweep_cov:
                        sweep_cov.add(r)
                        witnesses[r] = (a, b, c)
    # non-unit classes are also covered by divisibility scaling families
    for r in range(MOD):
        if math.gcd(r, MOD) > 1:
            sweep_cov.add(r)
    sweep_unc = sorted(set(range(MOD)) - sweep_cov)
    print(f"  Type I sweep (a,b<=210, ab|210, c|a+b) + divisibility families"
          f" leaves uncovered: {sweep_unc}")
    print(f"  Sweep uncovered == our uncovered : {sweep_unc == unc}")
    extra = sorted(sweep_cov - cov)
    print(f"  Classes the sweep covers that our 12 families miss: {extra}")
    return unc


# --------------------------------------------------------------------------
# 3. Brute-force existence for primes (bounded enumeration)
# --------------------------------------------------------------------------

def sieve_primes(n):
    comp = bytearray(n + 1)
    ps = []
    for i in range(2, n + 1):
        if not comp[i]:
            ps.append(i)
            for j in range(i * i, n + 1, i):
                comp[j] = 1
    return ps


def spf_sieve(n):
    spf = list(range(n + 1))
    for i in range(2, int(n ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, n + 1, i):
                if spf[j] == j:
                    spf[j] = i
    return spf


_div_cache = {}


def divisors_of_square(x, spf):
    """Sorted divisors of x^2."""
    got = _div_cache.get(x)
    if got is not None:
        return got
    fac = {}
    t = x
    while t > 1:
        q = spf[t]
        fac[q] = fac.get(q, 0) + 1
        t //= q
    divs = [1]
    for q, e in fac.items():
        divs = [d * q ** i for d in divs for i in range(2 * e + 1)]
    divs.sort()
    _div_cache[x] = divs
    return divs


def first_solution(p, spf):
    """First (x,y,z), x<=y<=z, with 4/p = 1/x+1/y+1/z; bounded enumeration."""
    for x in range(p // 4 + 1, 3 * p // 4 + 1):
        A = 4 * x - p
        B = p * x
        dmin = 2 * x * (2 * x - p) if 2 * x > p else 0
        d0s = divisors_of_square(x, spf)
        cands = list(d0s) + [d0 * p for d0 in d0s if d0 <= x]
        for d in cands:
            if d >= dmin and (B + d) % A == 0:
                y = (B + d) // A
                z = (B + (B * B) // d) // A
                if (B * B) % d == 0 and y <= z:
                    return (x, y, z)
    return None


def existence_check(pmax):
    print(f"== Brute-force existence check: all primes p < {pmax} ==")
    primes = sieve_primes(pmax - 1)
    spf = spf_sieve(3 * pmax // 4 + 2)
    hardest = []
    fails = []
    for p in primes:
        sol = first_solution(p, spf)
        if sol is None:
            fails.append(p)
            continue
        check_solution(p, sol)
        hardest.append((sol[0] - p // 4, p, sol))
    print(f"  primes checked: {len(primes)}, failures: {len(fails)} {fails}")
    assert not fails
    hardest.sort(reverse=True)
    print("  primes where the smallest denominator x sits farthest above p/4:")
    for off, p, sol in hardest[:5]:
        print(f"    p={p} (p mod 840 = {p % 840}): first solution {sol}")
    return len(primes)


# --------------------------------------------------------------------------
# 4. Representation counts f(p) via C helper
# --------------------------------------------------------------------------

C_SOURCE = r"""
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

int main(int argc, char **argv) {
    int pmax = argc > 1 ? atoi(argv[1]) : 100000;
    int xmax = 3 * pmax / 4 + 1;
    int *spf = calloc(xmax + 1, sizeof(int));
    for (int i = 2; i <= xmax; i++)
        if (!spf[i]) for (long j = i; j <= xmax; j += i)
            if (!spf[j]) spf[j] = i;
    uint64_t **divs = malloc((xmax + 1) * sizeof(uint64_t *));
    int *ndivs = calloc(xmax + 1, sizeof(int));
    for (int x = 1; x <= xmax; x++) {
        int pr[16], ex[16], np = 0, t = x;
        while (t > 1) {
            int q = spf[t], e = 0;
            while (t % q == 0) { t /= q; e++; }
            pr[np] = q; ex[np] = 2 * e; np++;
        }
        int cnt = 1;
        for (int i = 0; i < np; i++) cnt *= ex[i] + 1;
        uint64_t *D = malloc(cnt * sizeof(uint64_t));
        D[0] = 1; int m = 1;
        for (int i = 0; i < np; i++) {
            int m0 = m; uint64_t pe = 1;
            for (int e = 1; e <= ex[i]; e++) {
                pe *= (uint64_t)pr[i];
                for (int j = 0; j < m0; j++) D[m++] = D[j] * pe;
            }
        }
        divs[x] = D; ndivs[x] = cnt;
    }
    char *comp = calloc(pmax + 1, 1);
    for (long i = 2; i * i <= pmax; i++)
        if (!comp[i]) for (long j = i * i; j <= pmax; j += i) comp[j] = 1;
    for (int p = 2; p < pmax; p++) {
        if (p > 2 && comp[p]) continue;
        long count = 0;
        int xlo = p / 4 + 1, xhi = (3 * p) / 4;
        for (int x = xlo; x <= xhi; x++) {
            uint64_t A = (uint64_t)(4L * x - p);
            uint64_t B = (uint64_t)p * (uint64_t)x;
            long t2 = 2L * x - p;
            uint64_t dmin = t2 > 0 ? (uint64_t)(2L * x * t2) : 0;
            uint64_t *D = divs[x];
            int nd = ndivs[x];
            for (int i = 0; i < nd; i++) {
                uint64_t d0 = D[i];
                /* d = d0 (<= x^2 < px = B always since x < p) */
                if (d0 >= dmin && (B + d0) % A == 0) count++;
                if (d0 <= (uint64_t)x) {  /* d = d0*p <= B */
                    uint64_t d = d0 * (uint64_t)p;
                    if (d >= dmin && (B + d) % A == 0) count++;
                }
            }
        }
        printf("%d %ld\n", p, count);
    }
    return 0;
}
"""


def count_reps_python(p, spf):
    """Independent pure-Python f(p) (same math, separate code path)."""
    count = 0
    for x in range(p // 4 + 1, 3 * p // 4 + 1):
        A = 4 * x - p
        B = p * x
        dmin = 2 * x * (2 * x - p) if 2 * x > p else 0
        for d0 in divisors_of_square(x, spf):
            if d0 >= dmin and (B + d0) % A == 0:
                count += 1
            if d0 <= x:
                d = d0 * p
                if d >= dmin and (B + d) % A == 0:
                    count += 1
    return count


def run_counts(pmax, pycheck, dump):
    print(f"== Representation counts f(p) = #{{x<=y<=z}} for primes p < {pmax} ==")
    tmpd = tempfile.mkdtemp(prefix="es_cov_")
    src = os.path.join(tmpd, "es_count.c")
    exe = os.path.join(tmpd, "es_count")
    with open(src, "w") as f:
        f.write(C_SOURCE)
    cc = None
    for cand in ("cc", "gcc", "clang"):
        try:
            subprocess.run([cand, "--version"], capture_output=True, check=True)
            cc = cand
            break
        except Exception:
            continue
    if cc is None:
        print("  no C compiler found; falling back to (slow) pure Python")
        primes = sieve_primes(pmax - 1)
        spf = spf_sieve(3 * pmax // 4 + 2)
        counts = {p: count_reps_python(p, spf) for p in primes}
    else:
        subprocess.run([cc, "-O2", "-o", exe, src], check=True)
        out = subprocess.run([exe, str(pmax)], capture_output=True,
                             text=True, check=True).stdout
        counts = {}
        for line in out.splitlines():
            p, c = line.split()
            counts[int(p)] = int(c)
    print(f"  computed f(p) for {len(counts)} primes")

    if pycheck:
        spf = spf_sieve(3 * pycheck // 4 + 2)
        bad = []
        for p in sieve_primes(pycheck - 1):
            ref = count_reps_python(p, spf)
            if ref != counts[p]:
                bad.append((p, ref, counts[p]))
        print(f"  Python cross-check for p < {pycheck}: "
              f"{'ALL MATCH' if not bad else f'MISMATCHES {bad[:10]}'}")
        assert not bad

    zero = [p for p, c in counts.items() if c == 0]
    print(f"  primes with f(p) = 0: {zero if zero else 'none'}")
    assert not zero

    if dump:
        with open(dump, "w") as f:
            for p in sorted(counts):
                f.write(f"{p} {counts[p]}\n")
        print(f"  raw counts dumped to {dump}")

    # bottom-20 lists
    def show_bottom(items, label):
        print(f"  Bottom 20 primes by f(p) {label}:")
        print("    f(p)   p       p%24  p%120  p%840  QR840?")
        sq = {u * u % MOD for u in range(MOD) if math.gcd(u, MOD) == 1}
        for p, c in items:
            print(f"    {c:5d}  {p:6d}  {p % 24:4d}  {p % 120:5d}  "
                  f"{p % MOD:5d}  {p % MOD in sq}")

    order = sorted(counts.items(), key=lambda t: (t[1], t[0]))
    show_bottom([(p, c) for p, c in order[:20]], "(all primes)")
    big = [(p, c) for p, c in order if p > 1000][:20]
    show_bottom(big, "(primes > 1000)")

    # aggregate structure: f by residue class mod 24 for large p
    by24 = {}
    for p, c in counts.items():
        if p >= 10**4:
            by24.setdefault(p % 24, []).append(c)
    print("  mean/min f(p) by p mod 24 (primes in [10^4, pmax)):")
    for r in sorted(by24):
        v = by24[r]
        print(f"    p=={r:2d} mod 24: n={len(v):4d}  mean={sum(v)/len(v):9.1f}"
              f"  min={min(v):5d}  max={max(v):7d}")
    # uncovered-class primes
    sqset = {u * u % MOD for u in range(MOD) if math.gcd(u, MOD) == 1}
    hard = [(p, c) for p, c in counts.items() if p % MOD in sqset]
    hard_min = sorted(hard, key=lambda t: t[1])[:5]
    print(f"  primes in the 6 uncovered (QR) classes mod 840: {len(hard)}; "
          f"all have f(p) >= 1: {all(c > 0 for _, c in hard)}")
    print(f"  lowest among them: {hard_min}")
    return counts


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--coverage", action="store_true")
    ap.add_argument("--existence", action="store_true")
    ap.add_argument("--counts", action="store_true")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--nmax", type=int, default=200000,
                    help="identity verification bound (default 2*10^5)")
    ap.add_argument("--pmax", type=int, default=100000,
                    help="prime bound for existence/counts (default 10^5)")
    ap.add_argument("--pycheck", type=int, default=2000,
                    help="cross-check C counts with Python below this bound")
    ap.add_argument("--dump", type=str, default=None,
                    help="dump raw 'p f(p)' lines to this file")
    args = ap.parse_args()
    if args.all:
        args.verify = args.coverage = args.existence = args.counts = True
    if not (args.verify or args.coverage or args.existence or args.counts):
        ap.print_help()
        return
    if args.verify:
        verify_families(args.nmax)
    if args.coverage:
        coverage_report()
    if args.existence:
        existence_check(args.pmax)
    if args.counts:
        run_counts(args.pmax, args.pycheck, args.dump)


if __name__ == "__main__":
    main()
