#!/usr/bin/env python3
"""Check a set of bases of C^d for mutual unbiasedness.

    python harness/mub-six/mub_check.py bases.json [--tol 1e-9]

Input JSON: {"d": d, "bases": [...]}, where bases[b][c] is the c-th vector of
basis b as a list of [re, im] pairs. Reports, independently of whatever
produced the file:

- the worst orthonormality error of each basis,
- the defect L = sum_{a<b} sum_{i,j} (|<a_i|b_j>|^2 - 1/d)^2 and each pair's
  share of it,
- the average squared (Bengtsson) distance ASD = 1 - L / (C(k,2)(d-1)).

Exit status 0 iff every basis is orthonormal to --tol and L <= --tol. That is
a floating-point check only: it is evidence, not a certificate, that a
configuration is MU (see PROBLEM.md). Standard library only.
"""

from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations


def inner(u: list[complex], v: list[complex]) -> complex:
    return sum(x.conjugate() * y for x, y in zip(u, v))


def load(path: str) -> tuple[int, list[list[list[complex]]]]:
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    d = int(data["d"])
    bases = [[[complex(re, im) for re, im in vec] for vec in basis]
             for basis in data["bases"]]
    for b, basis in enumerate(bases):
        if len(basis) != d or any(len(v) != d for v in basis):
            raise ValueError(f"basis {b} is not {d} vectors of length {d}")
    return d, bases


def orthonormality_error(basis: list[list[complex]]) -> float:
    d = len(basis)
    return max(abs(inner(basis[i], basis[j]) - (1 if i == j else 0))
               for i in range(d) for j in range(d))


def pair_defect(a: list[list[complex]], b: list[list[complex]], d: int) -> tuple[float, float]:
    """(sum of squared deviations, max |deviation|) of |<a_i|b_j>|^2 from 1/d."""
    total = worst = 0.0
    for u in a:
        for w in b:
            dev = abs(inner(u, w)) ** 2 - 1 / d
            total += dev * dev
            worst = max(worst, abs(dev))
    return total, worst


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("path")
    ap.add_argument("--tol", type=float, default=1e-9)
    args = ap.parse_args()

    d, bases = load(args.path)
    k = len(bases)
    print(f"d = {d}, k = {k}")
    ok = True
    for b, basis in enumerate(bases):
        err = orthonormality_error(basis)
        ok &= err <= args.tol
        print(f"basis {b}: orthonormality error {err:.3e}")
    L = 0.0
    for (i, a), (j, b) in combinations(enumerate(bases), 2):
        s, worst = pair_defect(a, b, d)
        L += s
        print(f"pair ({i},{j}): defect {s:.6e}  max |dev| {worst:.6f}  D^2 {1 - s / (d - 1):.8f}")
    pairs = k * (k - 1) // 2
    asd = 1 - L / (pairs * (d - 1)) if pairs else 1.0
    print(f"L = {L:.15g}")
    print(f"ASD = {asd:.15g}")
    ok &= L <= args.tol
    print("MU (to tolerance): " + ("yes" if ok else "no"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
