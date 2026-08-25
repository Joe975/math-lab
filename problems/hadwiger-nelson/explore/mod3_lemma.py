"""The lattice invariant that tracks the census split.

The canonical 3-colouring of the triangular lattice is reduction modulo the
Eisenstein prime lambda = 1 + zeta (norm 3): the colour of a + b*zeta is
(a - b) mod 3, since zeta = -1 mod lambda.

Lemma.  N(a, b) = a^2 + ab + b^2 is divisible by 3 exactly when a = b (mod 3),
i.e. exactly when the point has colour 0.

So "3 divides n" says precisely: every lattice point of norm n carries colour 0
in the canonical 3-colouring.  The rotation rho_n joins each norm-n point to its
own image, so this is the statement that all of rho_n's built-in unit edges
start from one colour class.

This script checks the lemma exhaustively over a range and reports the colour
multiset of each norm class.

Run:  python explore/mod3_lemma.py [range]
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def _harness_dir():
    """Locate harness/hadwiger-nelson by walking up from this file."""
    d = HERE
    while True:
        cand = os.path.join(d, "harness", "hadwiger-nelson")
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(d)
        if parent == d:
            raise RuntimeError("could not locate harness/hadwiger-nelson")
        d = parent


sys.path.insert(0, _harness_dir())
sys.path.insert(0, HERE)

from unit_distance import eisenstein_norm  # noqa: E402


def main():
    rng = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    bad = []
    checked = 0
    for a in range(-rng, rng + 1):
        for b in range(-rng, rng + 1):
            checked += 1
            n = eisenstein_norm(a, b)
            if (n % 3 == 0) != ((a - b) % 3 == 0):
                bad.append((a, b, n))
    print(f"lemma checked over |a|,|b| <= {rng}: {checked} points, "
          f"{len(bad)} counterexamples")

    print("\nnorm : colours of its lattice points (canonical 3-colouring)")
    for n in range(1, 40):
        pts = [(a, b) for a in range(-8, 9) for b in range(-8, 9)
               if eisenstein_norm(a, b) == n]
        if not pts:
            continue
        cols = sorted({(a - b) % 3 for a, b in pts})
        print(f"{n:>4} : {cols}   ({len(pts)} points, 3|n = {n % 3 == 0})")
    return 0 if not bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
