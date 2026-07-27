"""How large can the admissible-set bound ever get?

The bound is  min { |S| : S admissible }.  Any *explicit* admissible set with
sum 1/p > 1 is therefore an upper bound on it -- and so a ceiling on everything
this relaxation can ever prove, no matter how far the branch limit X is pushed.

This script constructs such sets greedily (take the smallest prime that is still
allowed) from a few different starting points and reports their size.

Usage:
    python method_ceiling.py --pool 20000000 --skips "", "3", "3,5", "3,5,7"
"""

import argparse
import json
import os
import sys
from fractions import Fraction

sys.set_int_max_str_digits(2_000_000)

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.admissible import primes_upto  # noqa: E402


def greedy(primes, limit, is_prime_flag, skip):
    """Greedy admissible set, with the forbidden residue classes sieved.

    Testing each candidate against every chosen prime is O(pool * |S|) and does
    not finish.  Instead, when q joins S, mark every prime == 1 (mod q) as
    forbidden in one pass: sum over q of pool/q is about pool, not pool*|S|.
    """
    forbidden = bytearray(limit + 1)
    S, s = [], Fraction(0)
    for p in primes:
        if p in skip or forbidden[p]:
            continue
        S.append(p)
        s += Fraction(1, p)
        for m in range(1 + p, limit + 1, p):      # the class 1 mod p
            if is_prime_flag[m]:
                forbidden[m] = 1
        if s > 1:
            return S, s, True
    return S, s, False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pool", type=int, default=20_000_000)
    ap.add_argument("--skips", type=str, default="|3|3,5|3,5,7|3,5,7,11")
    ap.add_argument("--dump", type=str, default=None,
                    help="write the skip-{3} witness set to this JSON file")
    a = ap.parse_args()

    allp = primes_upto(a.pool)
    ps = [p for p in allp if p > 2]
    flag = bytearray(a.pool + 1)
    for p in allp:
        flag[p] = 1
    print(f"pool: {len(ps)} odd primes below {a.pool}")
    for spec in a.skips.split("|"):
        skip = {int(x) for x in spec.split(",") if x}
        S, s, ok = greedy(ps, a.pool, flag, skip)
        tag = f"skip {sorted(skip) if skip else 'nothing'}"
        if ok:
            prod = 1
            for p in S:
                prod *= p
            digits = len(str(prod))
            print(f"{tag:>22}: reached sum > 1 with |S| = {len(S)}, "
                  f"largest prime {S[-1]}, product has {digits} digits")
            if a.dump and skip == {3}:
                with open(a.dump, "w") as f:
                    json.dump({"skip": sorted(skip), "pool": a.pool,
                               "size": len(S), "digits": digits,
                               "primes": S}, f)
                print(f"    wrote {a.dump}")
        else:
            print(f"{tag:>22}: FAILED -- only reached {float(s):.6f} "
                  f"using all {len(S)} admissible primes in the pool")
    return 0


if __name__ == "__main__":
    sys.exit(main())
