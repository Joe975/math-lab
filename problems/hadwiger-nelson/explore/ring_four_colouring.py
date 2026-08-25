"""An explicit 4-colouring of the whole ring Z[zeta][rho], and its verification.

Setting.  zeta = exp(i pi/3), E = Z[zeta], lambda = 1 + zeta (norm 3), and rho
the spindle rotation with cos = 5/6.  Since 3 rho^2 = 5 rho - 3, the ring
Z[zeta][rho] is the increasing union of H_k = 3^-k (E + rho E), and a point of
H_k is written (Z + rho W)/3^k with Z, W in E.

THE COLOURING.   c(p) = (lambda Z + W)  mod 2E,  four values, since E/2E has
four elements.

Well defined.  The inclusion H_k -> H_{k+1} is (Z, W) -> (3Z, 3W), and
3 = 1 mod 2, so the value is unchanged.  One rule therefore colours every level
at once, hence the whole ring.

Proper.  Two lemmas.

  L1.  If |Z + rho W|^2 is rational then W is a rational multiple of Z (or one
       of them is zero).  Writing conj(Z) W = (p + q i sqrt3)/2, the squared
       length is N(Z) + N(W) + (5p - q sqrt33)/6, and sqrt33 is irrational, so
       q = 0; conj(Z) W real means W is parallel to Z.

  L2.  Every unit vector of H_k has |Z + rho W|^2 = 9^k, which is odd.

  Now suppose a unit vector had lambda Z + W = 0 mod 2E.  N(lambda) = 3 is odd,
  so lambda is invertible in E/2E.  By L1, bW = aZ in E for coprime integers
  a, b (the degenerate cases Z = 0 and W = 0 are immediate).  Reducing mod 2E
  and splitting on the parities of a and b -- they are not both even -- each
  case forces Z = W = 0 in E/2E.  Then 4 divides |Z + rho W|^2 = 9^k, which is
  odd.  Contradiction, so lambda Z + W is never 0 mod 2E and adjacent points
  always differ.  QED

CONSEQUENCE.  Every finite unit-distance graph whose points lie in Z[zeta][rho]
is 4-colourable.  The Moser spindle lies in H_0, so the chromatic number of the
ring's unit-distance graph is exactly 4, and no construction inside this ring
can ever witness chi(R^2) >= 5.

This script re-derives none of that; it *checks* it, on real balls built by the
independent group machinery, edge by edge.

Run:  python problems/hadwiger-nelson/explore/ring_four_colouring.py
"""

from __future__ import annotations

import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def _harness_dir():
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

from homomorphism_colouring import eis_mul  # noqa: E402
from rho_group import (  # noqa: E402
    RhoGroup,
    ball,
    edges_from_unit_vectors,
    unit_vectors_structured,
)
from unit_distance import eisenstein_norm  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))

LAMBDA = (1, 1)  # 1 + zeta, norm 3


def colour(Z, W):
    """c(Z, W) = (lambda Z + W) mod 2E, encoded as one of four values."""
    lz = eis_mul(LAMBDA, Z)
    a = (lz[0] + W[0]) % 2
    b = (lz[1] + W[1]) % 2
    return a * 2 + b


def check_level(group, level, depth):
    gens = unit_vectors_structured(9 ** level)
    bad_gen = [g for g in gens if colour(g[0], g[1]) == 0]
    V = ball(group, gens, depth)
    E = edges_from_unit_vectors(V, gens)
    cols = [colour(c[0], c[1]) for c in V]
    mono = [(a, b) for a, b in E if cols[a] == cols[b]]
    return {
        "level": level,
        "unit_vectors": len(gens),
        "unit_vectors_with_colour_zero": len(bad_gen),
        "depth": depth,
        "vertices": len(V),
        "edges": len(E),
        "monochromatic_edges": len(mono),
        "colours_used": len(set(cols)),
        "proper": not mono and not bad_gen,
    }


def check_generator_lemma(group, max_level):
    """The whole proof reduces to: no unit vector has lambda Z + W = 0 mod 2E."""
    rows = []
    for level in range(max_level + 1):
        gens = unit_vectors_structured(9 ** level)
        zeros = [g for g in gens if colour(g[0], g[1]) == 0]
        parallel_violations = [
            g for g in gens
            if g[0] != (0, 0) and g[1] != (0, 0)
            and g[1][0] * g[0][1] - g[0][0] * g[1][1] != 0
        ]
        rows.append({
            "level": level,
            "unit_vectors": len(gens),
            "with_lambda_Z_plus_W_zero_mod_2": len(zeros),
            "non_parallel_unit_vectors": len(parallel_violations),
            "norms_seen": sorted({(eisenstein_norm(*g[0]), eisenstein_norm(*g[1]))
                                  for g in gens})[:8],
        })
    return rows


def main():
    group = RhoGroup(2)
    print("Lemma check: every unit vector has parallel components and is "
          "non-zero under the colour map\n")
    lem = check_generator_lemma(group, 8)
    for r in lem:
        print(f"  level {r['level']}: {r['unit_vectors']:>4} unit vectors, "
              f"{r['non_parallel_unit_vectors']} non-parallel, "
              f"{r['with_lambda_Z_plus_W_zero_mod_2']} killed by the colour map")

    print("\nBall checks (the colouring applied to real unit-distance graphs)\n")
    plan = [(0, 4), (1, 3), (2, 2), (3, 2), (4, 2), (5, 2)]
    rows = []
    for level, depth in plan:
        t0 = time.time()
        r = check_level(group, level, depth)
        rows.append(r)
        print(f"  level {r['level']} depth {r['depth']}: "
              f"|V|={r['vertices']:>6} |E|={r['edges']:>7} "
              f"mono={r['monochromatic_edges']} colours={r['colours_used']} "
              f"PROPER={r['proper']}  ({time.time() - t0:.1f}s)", flush=True)

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "ring_four_colouring.json")
    with open(path, "w") as fh:
        json.dump({"rule": "c(Z, W) = (lambda Z + W) mod 2E, lambda = 1 + zeta",
                   "lemma": lem, "balls": rows}, fh, indent=2)
    print(f"\nwrote {path}")
    ok = all(r["proper"] for r in rows) and all(
        r["with_lambda_Z_plus_W_zero_mod_2"] == 0 for r in lem)
    print("ALL PROPER" if ok else "A CHECK FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
