"""The ladder H_k = 3^-k (E + rho E), and the chromatic numbers of its boxes.

rho satisfies 3 rho^2 = 5 rho - 3, so rho is not an algebraic integer and the
ring Z[zeta][rho] is not a finitely generated Z[zeta]-module.  It is instead the
increasing union

    H_0 subset H_1 subset H_2 subset ...,    H_k = 3^-k (E + rho E),

and every unit-distance graph on points of the ring lies in some H_k.  Scaling
by 3^k is a similarity, so

    the unit-distance graph on H_k is the Cayley graph on E + rho E
    (a rank-4 group) whose generators are the vectors of squared length 9^k.

That is the whole ladder in one line: same vertex group at every level, a
strictly richer generator set as k grows.  Counts: 18 generators at k = 0,
42 at k = 1, 66 at k = 2.

Vertex sets are coefficient boxes {(z, w) : N(z) <= r, N(w) <= r}.  Edges come
from the generator table; a pair-by-pair exact recomputation on a prefix
cross-checks it.

Run:
  python problems/hadwiger-nelson/explore/rho_ladder.py --level 1 \
      --boxes "1 3 4 7 9 12" --verify-pairs 400
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from fractions import Fraction

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

import colouring_fast as F  # noqa: E402
from lattice import eisenstein_patch  # noqa: E402
from rho_group import RhoGroup, edges_from_unit_vectors  # noqa: E402
from unit_distance import Point, build_graph, certify, eisenstein_norm  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def scaled_point(group, coeffs, level):
    p = group.point(coeffs)
    s = group.field.rational(Fraction(1, 3 ** level))
    return Point(p.x * s, p.y * s)


def box_vertices(r):
    """Coefficient box {(z, w) : N(z) <= r, N(w) <= r} as coefficient tuples."""
    patch = eisenstein_patch(r)
    return [(z, w) for z in patch for w in patch]


def generator_span(gens):
    return max(max(eisenstein_norm(*ab) for ab in g) for g in gens)


def analyse(group, level, r, gens, verify_pairs, node_limit):
    t0 = time.time()
    coeffs = box_vertices(r)
    edges = edges_from_unit_vectors(coeffs, gens)
    n = len(coeffs)
    t_build = time.time() - t0

    xc = {"checked": 0, "agree": None}
    if verify_pairs > 0:
        k = min(verify_pairs, n)
        pts = [scaled_point(group, c, level) for c in coeffs[:k]]
        slow, _ = build_graph(pts)
        rep = certify(pts, slow)
        fast = sorted((a, b) for a, b in edges if a < k and b < k)
        xc = {
            "checked": k,
            "pairs": k * (k - 1) // 2,
            "certifier_ok": rep["ok"],
            "agree": sorted(map(tuple, slow)) == fast,
            "slow_edges": len(slow),
            "fast_edges": len(fast),
        }

    row = {
        "level": level, "box_norm": r,
        "generators": len(gens),
        "vertices": n, "edges": len(edges),
        "mean_degree": round(2 * len(edges) / n, 2) if n else 0,
        "cross_check": xc,
        "build_seconds": round(t_build, 2),
    }

    kept, sub, idx, shed = F.k_core(n, edges, 4)
    row["four_core_vertices"] = len(kept)
    row["four_core_edges"] = len(sub)

    t1 = time.time()
    try:
        col = F.k_colourable(n, edges, 4, node_limit=node_limit)
        row["four_colourable"] = col is not None
        if col is not None:
            row["four_colouring_verified"] = F.verify(n, edges, col, 4)
        else:
            row["chi_at_least"] = 5
    except TimeoutError:
        row["four_colourable"] = "timeout"
    row["four_seconds"] = round(time.time() - t1, 2)
    return row, coeffs, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--n", type=int, default=3)
    ap.add_argument("--boxes", default="1 3 4 7 9 12 13")
    ap.add_argument("--verify-pairs", type=int, default=300)
    ap.add_argument("--node-limit", type=int, default=0)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    group = RhoGroup(2, n=args.n)
    target = 9 ** args.level if args.n == 3 else None
    if target is None:
        raise SystemExit("the 3 rho^2 = 5 rho - 3 ladder is specific to n = 3")

    t0 = time.time()
    gens = group.vectors_of_norm(target)
    print(f"level {args.level}: |v|^2 = {target} gives {len(gens)} generators "
          f"(max coefficient norm {generator_span(gens)}), "
          f"{time.time() - t0:.1f}s")

    out = args.out or os.path.join(OUT, f"rho_ladder_L{args.level}.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    rows = []
    hdr = (f"{'box':>5} {'|V|':>7} {'|E|':>9} {'deg':>7} {'4core':>7} "
           f"{'4col':>8} {'xchk':>5} {'build':>7} {'4col_s':>8}")
    print(hdr)
    print("-" * len(hdr))
    for r in [int(x) for x in args.boxes.split()]:
        row, coeffs, edges = analyse(
            group, args.level, r, gens, args.verify_pairs, args.node_limit)
        rows.append(row)
        print(f"{r:>5} {row['vertices']:>7} {row['edges']:>9} "
              f"{row['mean_degree']:>7} {row['four_core_vertices']:>7} "
              f"{str(row['four_colourable']):>8} "
              f"{str(row['cross_check']['agree']):>5} "
              f"{row['build_seconds']:>7} {row['four_seconds']:>8}")
        with open(out, "w") as fh:
            json.dump(rows, fh, indent=2)
        if row["four_colourable"] is False:
            print("\n*** NOT 4-COLOURABLE -- extract and verify a core ***")
            break
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
