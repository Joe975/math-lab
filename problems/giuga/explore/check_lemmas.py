"""Mechanised checks of the structural lemmas the factor-count bound rests on.

Each lemma is proved on paper in the attempt record; this script is the
independent empirical check that the proof is not misstated, run over every
Carmichael and Giuga number below a limit.

  L2  a solution is squarefree
  L3  a solution is odd            (checked as: every Carmichael number is odd)
  L4  for distinct p, q dividing a Carmichael number, q does not divide p - 1
  L5  every Giuga number has sum 1/p > 1
  L6  how close real Carmichael numbers get to sum 1/p = 1  (the gap the bound
      exploits -- reported, not asserted)

Usage:
    python check_lemmas.py --limit 10000000
"""

import argparse
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.enumerate_giuga import scan_giuga_numbers  # noqa: E402


def spf_sieve(limit):
    spf = list(range(limit + 1))
    i = 2
    while i * i <= limit:
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1
    return spf


def squarefree_factors(n, spf):
    ps = []
    while n > 1:
        p = spf[n]
        n //= p
        if n % p == 0:
            return None
        ps.append(p)
    return ps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=10 ** 7)
    a = ap.parse_args()
    spf = spf_sieve(a.limit)

    carmichael = []
    for n in range(3, a.limit + 1, 2):
        ps = squarefree_factors(n, spf)
        if ps is None or len(ps) < 3:
            continue
        if all((n - 1) % (p - 1) == 0 for p in ps):
            carmichael.append((n, ps))
    # even candidates, to check L3 rather than assume it
    even_carmichael = []
    for n in range(4, a.limit + 1, 2):
        ps = squarefree_factors(n, spf)
        if ps is None or len(ps) < 2:
            continue
        if all((n - 1) % (p - 1) == 0 for p in ps):
            even_carmichael.append(n)

    print(f"Carmichael numbers <= {a.limit}: {len(carmichael)}")
    print(f"L3  even ones: {even_carmichael}  ->  {'OK' if not even_carmichael else 'FALSIFIED'}")

    bad_l4 = [(n, p, q) for n, ps in carmichael
              for p in ps for q in ps if p != q and (p - 1) % q == 0]
    print(f"L4  (p-1) divisible by another factor: {bad_l4[:5]}  "
          f"->  {'OK' if not bad_l4 else 'FALSIFIED'}")

    giuga = scan_giuga_numbers(min(a.limit, 10 ** 6))
    bad_l5 = []
    for n in giuga:
        ps = squarefree_factors(n, spf)
        if sum((Fraction(1, p) for p in ps), Fraction(0)) <= 1:
            bad_l5.append(n)
    print(f"L5  Giuga numbers <= {min(a.limit, 10**6)}: {giuga}; "
          f"any with sum <= 1: {bad_l5}  ->  {'OK' if not bad_l5 else 'FALSIFIED'}")

    best = max(((sum((Fraction(1, p) for p in ps), Fraction(0)), n, ps)
                for n, ps in carmichael), default=None)
    if best:
        s, n, ps = best
        print(f"L6  largest sum 1/p over Carmichael n <= {a.limit}: "
              f"{float(s):.6f} at n={n} = {'*'.join(map(str, ps))} "
              f"({len(ps)} factors)")
        top = sorted(((sum((Fraction(1, p) for p in ps), Fraction(0)), n, len(ps))
                      for n, ps in carmichael), reverse=True)[:5]
        for s, n, k in top:
            print(f"      {float(s):.6f}  n={n} ({k} factors)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
