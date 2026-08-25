"""How many unit distances are there between a lattice patch and its rotation?

For a rotation rho_n about the origin, the cross-copy edges of
P(m) union rho_n(P(m)) are the pairs (p, rho_n q) of lattice points at distance
exactly 1.  This script measures that count as the patch grows, which is the
quantity that decides whether enlarging the patch can buy anything.

It also tests whether the union collapses into a single scaled triangular
lattice, by checking whether every pairwise squared distance in the union is a
rational number with denominator dividing n (a necessary condition for the
union to sit inside a lattice commensurable with the Eisenstein lattice).

Run:  python explore/interface_growth.py [--ns ...] [--patches ...]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction

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

from lattice import (  # noqa: E402
    eisenstein_patch,
    norms_up_to,
    rotated_union,
    rotation_field,
    rotation_for_norm,
)
from unit_distance import build_graph, sq_dist  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def commensurable(n):
    """True iff 4n-1 = 3m^2, i.e. n = 3k^2 + 3k + 1 for an integer k >= 0.

    Exactly the n for which sin(theta_n) lies in Q(sqrt3), so that rho_n has all
    matrix entries in Q(sqrt3) and maps the Eisenstein lattice into a lattice
    commensurable with it.
    """
    v = 4 * n - 1
    if v % 3:
        return False
    q = v // 3
    r = int(q ** 0.5)
    for cand in (r - 1, r, r + 1):
        if cand >= 0 and cand * cand == q:
            return True
    return False


def rationality_profile(points):
    """Denominators of the pairwise squared distances, when they are rational."""
    dens = set()
    irrational = 0
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            d2 = sq_dist(points[i], points[j])
            if d2.is_rational():
                dens.add(d2.as_fraction().denominator)
            else:
                irrational += 1
    return sorted(dens), irrational


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default="")
    ap.add_argument("--patches", default="1,3,4,7,9,12,13,16,19,21,25,28")
    ap.add_argument("--out", default=os.path.join(OUT, "interface_growth.json"))
    args = ap.parse_args()

    ns = ([int(x) for x in args.ns.split(",")] if args.ns
          else norms_up_to(21))
    patches = [int(x) for x in args.patches.split(",")]

    rows = []
    hdr = (f"{'n':>4} {'comm':>5} {'patch':>6} {'|P|':>5} {'|V|':>5} {'|E|':>6} "
           f"{'cross':>6} {'d2 denominators':>28} {'irr':>5}")
    print(hdr)
    print("-" * len(hdr))
    for n in ns:
        field = rotation_field(n)
        rot = rotation_for_norm(n, field)
        for m in patches:
            patch = eisenstein_patch(m)
            points, tags = rotated_union(field, patch, rot, copies=2)
            edges, stats = build_graph(points)
            cross = sum(1 for a, b in edges if tags[a][0] != tags[b][0])
            dens, irr = rationality_profile(points)
            row = {
                "n": n,
                "commensurable": commensurable(n),
                "patch_norm": m,
                "patch_size": len(patch),
                "vertices": len(points),
                "edges": len(edges),
                "cross_copy_edges": cross,
                "sqdist_denominators": dens,
                "irrational_sqdists": irr,
            }
            rows.append(row)
            print(f"{n:>4} {str(row['commensurable']):>5} {m:>6} {len(patch):>5} "
                  f"{len(points):>5} {len(edges):>6} {cross:>6} "
                  f"{str(dens)[:28]:>28} {irr:>5}")
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(rows, fh, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
