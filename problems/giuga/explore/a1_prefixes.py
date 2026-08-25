"""The a = 1 prefixes, by recursion over prime sets rather than by sieving.

An a = 1 prefix is a squarefree N with sum_{p|N} 1/p + 1/N = 1 -- the case where
the last-two rectangle degenerates to (p-N)(q-N) = N^2-1 and completions are
plentiful.  These were first found by scanning every integer below a limit; this
enumerates them the other way round, by building prime sets, so the two share no
arithmetic and either would expose a bug in the other.

Pruning (sound): a completion has N_final >= N0, and its reciprocal sum is
1 - 1/N_final >= 1 - 1/N0.  The extra primes are all > last and their product
must keep N <= limit, which caps how many there can be; if the prefix sum plus
the largest sum those extras could contribute is still below 1 - 1/N0, no
completion exists.

Usage:
    python problems/giuga/explore/a1_prefixes.py --limit 20000000
"""

import argparse
import os
import sys
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.admissible import primes_upto  # noqa: E402


def search(limit):
    primes = primes_upto(limit)
    out = []

    def rec(last_idx, N, s):
        # record: a = 1 means s + 1/N == 1 exactly
        if N > 1 and s + Fraction(1, N) == 1:
            out.append(N)
        # how many more primes can fit, and the most they could contribute
        best, prod, i, extra = Fraction(0), N, last_idx + 1, []
        while i < len(primes):
            if prod * primes[i] > limit:
                break
            prod *= primes[i]
            best += Fraction(1, primes[i])
            extra.append(primes[i])
            i += 1
        if not extra:
            return
        if s + best < 1 - Fraction(1, N):
            return
        for j in range(last_idx + 1, len(primes)):
            p = primes[j]
            if N * p > limit:
                break
            rec(j, N * p, s + Fraction(1, p))

    rec(-1, 1, Fraction(0))
    return sorted(set(out))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=20_000_000)
    a = ap.parse_args()
    sols = search(a.limit)
    print(f"squarefree N <= {a.limit} with sum_{{p|N}} 1/p + 1/N = 1, by prime-set recursion:")
    for n in sols:
        fs = []
        m = n
        for p in primes_upto(int(n ** 0.5) + 1):
            while m % p == 0:
                fs.append(p)
                m //= p
        if m > 1:
            fs.append(m)
        print(f"   {n} = {' * '.join(map(str, fs))}   {'ODD' if n % 2 else 'even'}")
    print(f"odd ones: {[n for n in sols if n % 2]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
