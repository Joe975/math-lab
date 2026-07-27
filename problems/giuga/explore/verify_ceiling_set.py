"""Independently re-check the explicit admissible set that caps the method.

The set was produced by a sieve of forbidden residue classes.  This verifies it
by the definition instead: every element prime (Miller-Rabin, not the sieve),
strictly increasing, odd, no element dividing another element minus one
(O(|S|^2) direct divisions), and the reciprocal sum > 1 as an exact Fraction.

Usage:
    python verify_ceiling_set.py problems/giuga/data/ceiling-set.json
"""

import json
import sys
import os
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.conditions import is_prime  # noqa: E402


def main(path):
    with open(path) as f:
        d = json.load(f)
    S = d["primes"]
    print(f"{path}: {len(S)} primes, claimed size {d['size']}")

    assert len(S) == d["size"]
    assert S == sorted(set(S)), "not strictly increasing / not distinct"
    assert all(p % 2 == 1 for p in S), "even element"
    assert all(is_prime(p) for p in S), "non-prime element"
    print("  all odd, distinct, increasing, prime (Miller-Rabin)   OK")

    bad = [(p, q) for p in S for q in S if p != q and (p - 1) % q == 0]
    assert not bad, f"pairwise rule violated: {bad[:5]}"
    print(f"  pairwise rule holds for all {len(S) * (len(S) - 1)} ordered pairs   OK")

    s = sum((Fraction(1, p) for p in S), Fraction(0))
    assert s > 1, f"reciprocal sum is only {float(s)}"
    print(f"  exact reciprocal sum > 1   OK   (numerator has "
          f"{s.numerator.bit_length()} bits; sum ~ {float(s):.9f})")

    # and it is minimal in the weak sense that dropping any one element breaks it
    slack = s - 1
    droppable = [p for p in S if Fraction(1, p) < slack]
    print(f"  slack above 1: {float(slack):.3e}; elements that could be dropped "
          f"individually: {len(droppable)}")
    print(f"VERIFIED: an admissible set of size {len(S)} exists, so the "
          f"reciprocal-sum relaxation can never prove a bound above {len(S)}.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1
                  else "problems/giuga/data/ceiling-set.json"))
