"""Chromatic numbers of finite pieces of G_m = E + rho E + ... + rho^(m-1) E.

Vertex sets are coefficient boxes: {sum_j rho^j z_j : N(z_j) <= r_j}.  Edges
come from the exhaustively enumerated unit vectors of the group, which is
exact; a cross-check against the pair-by-pair exact certifier is run on the
small cases so the fast path is not trusted on its own.

Run:
  python problems/hadwiger-nelson/explore/rho_group_census.py --m 2 \
      --boxes "3,3 3,7 7,7 7,13 13,13" --verify-pairs 400
"""

from __future__ import annotations

import argparse
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

import colouring_fast as F  # noqa: E402
from rho_group import (  # noqa: E402
    RhoGroup,
    build_vertices,
    describe_unit_vectors,
    edges_from_unit_vectors,
)
from unit_distance import build_graph, certify  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def cross_check_edges(group, points, coeffs, edges, limit):
    """Rebuild the graph pair-by-pair on a prefix and compare edge sets.

    The group route and the O(V^2) certifier route share no logic: one looks up
    coefficient differences in a table of unit vectors, the other computes every
    squared distance in the field.
    """
    if limit <= 0 or len(points) <= 1:
        return {"checked": 0, "agree": None}
    k = min(limit, len(points))
    sub_pts = points[:k]
    slow_edges, _ = build_graph(sub_pts)
    rep = certify(sub_pts, slow_edges)
    fast_edges = sorted((a, b) for a, b in edges if a < k and b < k)
    return {
        "checked": k,
        "pairs": k * (k - 1) // 2,
        "certifier_ok": rep["ok"],
        "agree": sorted(map(tuple, slow_edges)) == fast_edges,
        "slow_edges": len(slow_edges),
        "fast_edges": len(fast_edges),
    }


def analyse(group, patch_norms, verify_pairs, node_limit, want_core):
    t0 = time.time()
    uv = group.unit_vectors(group.coeff_bound(patch_norms))
    points, coeffs = build_vertices(group, patch_norms)
    edges = edges_from_unit_vectors(coeffs, uv)
    t_build = time.time() - t0

    xc = cross_check_edges(group, points, coeffs, edges, verify_pairs)

    n = len(points)
    row = {
        "m": group.m,
        "patch_norms": list(patch_norms),
        "unit_vectors": len(uv),
        "unit_vector_profiles": {str(k): v
                                 for k, v in describe_unit_vectors(group, uv).items()},
        "vertices": n,
        "edges": len(edges),
        "mean_degree": round(2 * len(edges) / n, 3) if n else 0,
        "cross_check": xc,
        "build_seconds": round(t_build, 2),
    }

    t1 = time.time()
    kept, sub, idx, shed = F.k_core(n, edges, 4)
    row["four_core_vertices"] = len(kept)
    row["four_core_edges"] = len(sub)

    try:
        col4 = F.k_colourable(n, edges, 4, node_limit=node_limit)
        row["four_colourable"] = col4 is not None
        if col4 is not None:
            row["four_colouring_verified"] = F.verify(n, edges, col4, 4)
        row["chi_seconds"] = round(time.time() - t1, 2)
    except TimeoutError:
        row["four_colourable"] = "timeout"
        row["chi_seconds"] = round(time.time() - t1, 2)

    if row.get("four_colourable") is False:
        row["chi_at_least"] = 5
    elif row.get("four_colourable") is True:
        col3 = F.k_colourable(n, edges, 3)
        row["three_colourable"] = col3 is not None
        row["chi"] = 3 if col3 is not None else 4
    return row, points, coeffs, edges


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--m", type=int, default=2)
    ap.add_argument("--n", type=int, default=3, help="rotation parameter for rho")
    ap.add_argument("--boxes", default="3,3 7,7 13,13")
    ap.add_argument("--verify-pairs", type=int, default=300)
    ap.add_argument("--node-limit", type=int, default=0)
    ap.add_argument("--cores", action="store_true")
    ap.add_argument("--out", default=os.path.join(OUT, "rho_group_census.json"))
    args = ap.parse_args()

    group = RhoGroup(args.m, n=args.n)
    boxes = [tuple(int(x) for x in b.split(",")) for b in args.boxes.split()]
    for b in boxes:
        if len(b) != args.m:
            raise SystemExit(f"box {b} does not have m={args.m} entries")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    rows = []
    hdr = (f"{'box':>14} {'|U|':>4} {'|V|':>6} {'|E|':>7} {'deg':>6} "
           f"{'4core':>6} {'4col':>6} {'chi':>4} {'xchk':>5} {'s':>8}")
    print(hdr)
    print("-" * len(hdr))
    for b in boxes:
        row, points, coeffs, edges = analyse(
            group, b, args.verify_pairs, args.node_limit, args.cores)
        rows.append(row)
        print(f"{str(b):>14} {row['unit_vectors']:>4} {row['vertices']:>6} "
              f"{row['edges']:>7} {row['mean_degree']:>6} "
              f"{row['four_core_vertices']:>6} {str(row['four_colourable']):>6} "
              f"{str(row.get('chi', '>=5')):>4} "
              f"{str(row['cross_check']['agree']):>5} "
              f"{row['build_seconds'] + row['chi_seconds']:>8.1f}")
        with open(args.out, "w") as fh:
            json.dump(rows, fh, indent=2)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
