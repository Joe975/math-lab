"""Recover the small Giuga numbers two ways and check the answers agree.

Usage:
    python calibrate_giuga_numbers.py [--scan-limit N] [--max-factors M]
"""

import argparse
import os
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.enumerate_giuga import (  # noqa: E402
    enumerate_by_prime_sets, k_upper_bound, product, scan_giuga_numbers, verify_set,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-limit", type=int, default=10 ** 6)
    ap.add_argument("--max-factors", type=int, default=6)
    ap.add_argument("--min-prime", type=int, default=2)
    a = ap.parse_args()

    t0 = time.time()
    scanned = scan_giuga_numbers(a.scan_limit)
    t_scan = time.time() - t0
    print(f"sieve scan n <= {a.scan_limit}: {scanned}   ({t_scan:.1f}s)")

    found = {}
    nodes = [0]
    for m in range(3, a.max_factors + 1):
        t0 = time.time()
        nodes[0] = 0
        sets = sorted(enumerate_by_prime_sets(m, min_prime=a.min_prime,
                                              on_node=lambda: nodes.__setitem__(0, nodes[0] + 1)))
        dt = time.time() - t0
        assert all(verify_set(s) for s in sets), "prime-set enumerator returned a non-Giuga set"
        found[m] = [product(s) for s in sets]
        print(f"m={m}  k<= {k_upper_bound(m, a.min_prime)}  nodes={nodes[0]:>9}  "
              f"{dt:8.2f}s  count={len(sets)}")
        for s in sets:
            print(f"        {product(s)} = {' * '.join(map(str, s))}")

    # cross-check: every Giuga number found by the sieve with <= max_factors
    # prime factors must appear in the prime-set enumeration, and vice versa
    # for those below the scan limit.
    from_sets = sorted(n for lst in found.values() for n in lst)
    small_from_sets = [n for n in from_sets if n <= a.scan_limit]
    print()
    print(f"sieve  -> {scanned}")
    print(f"sets(<= scan limit, <= {a.max_factors} factors) -> {small_from_sets}")
    print("AGREE" if scanned == small_from_sets else "*** DISAGREE ***")
    return 0 if scanned == small_from_sets else 1


if __name__ == "__main__":
    sys.exit(main())
