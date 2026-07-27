"""A second, deliberately naive implementation of the factor-count bound.

Shares nothing with harness/giuga/admissible.py: no bitmasks, no scaled
integers, no recursion, no pruning.  It enumerates *every* subset of the odd
primes <= X with itertools, throws away the ones that violate the pairwise rule,
and evaluates the leaf bound for each survivor with exact Fractions.

Exponential in pi(X), so only usable for small X -- which is the point: it is
an oracle for the pruned search, not a competitor to it.

Usage:
    python independent_bound.py --x 43 --pool 200000
"""

import argparse
import itertools
import sys
from fractions import Fraction


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    i = 2
    while i * i <= n:
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        i += 1
    return [i for i in range(n + 1) if sieve[i]]


def pairwise_ok(T):
    return all(not ((p - 1) % q == 0) for p in T for q in T if p != q)


def leaf_bound(T, X, tail_primes, cap):
    """|T| + (fewest primes above X, T-admissible, that can close the gap to 1)."""
    s = sum((Fraction(1, p) for p in T), Fraction(0))
    if s > 1:
        return len(T)
    t = 0
    for q in tail_primes:
        if any((q - 1) % p == 0 for p in T):
            continue
        s += Fraction(1, q)
        t += 1
        if s > 1:
            return len(T) + t
        if t > cap:
            return None
    # Pool exhausted: the t in-pool primes are the largest reciprocals available
    # above X, and together they still do not reach 1, so any completion needs
    # strictly more than t of them.  len(T) + t + 1 is therefore still a valid
    # lower bound (just a weak one).
    return len(T) + t + 1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x", type=int, default=43)
    ap.add_argument("--pool", type=int, default=200000)
    ap.add_argument("--cap", type=int, default=20000)
    a = ap.parse_args()

    ps = primes_upto(a.pool)
    small = [p for p in ps if 2 < p <= a.x]
    tail = [p for p in ps if p > a.x]
    print(f"X={a.x}: {len(small)} branchable primes, {2 ** len(small)} raw subsets")

    best, witness, kept = None, None, 0
    for r in range(len(small) + 1):
        for T in itertools.combinations(small, r):
            if not pairwise_ok(T):
                continue
            kept += 1
            b = leaf_bound(list(T), a.x, tail, a.cap)
            if b is not None and (best is None or b < best):
                best, witness = b, list(T)
    print(f"admissible subsets examined: {kept}")
    print(f"L({a.x}) >= {best}   witness T={witness}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
