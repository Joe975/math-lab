"""Census: how much chromatic number does one lattice rotation buy?

For each Eisenstein norm n we form the rotation rho_n about the origin with
cos(theta) = 1 - 1/(2n) -- the unique angle (up to sign) that sends a lattice
point of norm n to unit distance from itself -- and build the unit-distance
graph on

    P(m) union rho_n(P(m)) union ... union rho_n^(c-1)(P(m)),

where P(m) is the Eisenstein patch of norm <= m.  Every edge and non-edge is
decided by exact arithmetic in Q(sqrt3, sqrt(4n-1)).

Reports |V|, |E|, chi, the cross-copy edge count, and (optionally) the size of
a vertex-critical non-3-colourable core.

Run:  python explore/rotation_census.py [--max-n N] [--patches 3,4,7] [--copies 2,3]
      [--cores] [--out FILE]
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

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

import colouring as C  # noqa: E402
from exact_field import squarefree_part  # noqa: E402
from lattice import (  # noqa: E402
    eisenstein_patch,
    norms_up_to,
    rotated_union,
    rotation_field,
    rotation_for_norm,
)
from unit_distance import build_graph, certify  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def run_one(n, patch_norm, copies, want_core=False, certify_full=True):
    t0 = time.time()
    field = rotation_field(n)
    rot = rotation_for_norm(n, field)
    if not rot.is_orthogonal():
        raise AssertionError(f"rotation for n={n} is not orthogonal")
    patch = eisenstein_patch(patch_norm)
    points, tags = rotated_union(field, patch, rot, copies=copies)
    edges, stats = build_graph(points)
    if certify_full:
        rep = certify(points, edges)
        if not rep["ok"]:
            raise AssertionError(f"certificate failed for n={n}: {rep}")
    cross = sum(1 for a, b in edges if tags[a][0] != tags[b][0])
    chi = C.chromatic_number(len(points), edges, lo=1, hi=6)
    col = C.k_colouring(len(points), edges, chi)
    proper = C.verify_colouring(len(points), edges, col)
    row = {
        "n": n,
        "radicand_4n_minus_1": 4 * n - 1,
        "squarefree_part": squarefree_part(4 * n - 1)[1],
        "field": repr(field),
        "field_degree": field.dim,
        "patch_norm": patch_norm,
        "patch_size": len(patch),
        "copies": copies,
        "vertices": len(points),
        "edges": len(edges),
        "cross_copy_edges": cross,
        "chi": chi,
        "colouring_verified": proper,
        "seconds": round(time.time() - t0, 2),
    }
    if want_core and chi >= 4:
        core = C.vertex_critical_subgraph(len(points), edges, 3)
        kept, sub = core
        row["core_vertices"] = len(kept)
        row["core_edges"] = len(sub)
        row["core_tags"] = [list(tags[i]) for i in kept]
        row["core_degree_sequence"] = sorted(
            len(a) for a in C.adjacency(len(kept), sub))
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=49)
    ap.add_argument("--patches", default="3,4,7")
    ap.add_argument("--copies", default="2")
    ap.add_argument("--cores", action="store_true")
    ap.add_argument("--out", default=os.path.join(OUT, "rotation_census.json"))
    args = ap.parse_args()

    patches = [int(x) for x in args.patches.split(",")]
    copies = [int(x) for x in args.copies.split(",")]
    ns = [n for n in norms_up_to(args.max_n) if n >= 1]

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    rows = []
    header = (f"{'n':>4} {'4n-1':>6} {'sf':>5} {'field':>18} {'patch':>6} "
              f"{'cp':>3} {'|V|':>5} {'|E|':>6} {'cross':>6} {'chi':>4} {'core':>5} {'s':>7}")
    print(header)
    print("-" * len(header))
    for m in patches:
        for c in copies:
            for n in ns:
                row = run_one(n, m, c, want_core=args.cores)
                rows.append(row)
                print(f"{row['n']:>4} {row['radicand_4n_minus_1']:>6} "
                      f"{row['squarefree_part']:>5} {row['field']:>18} "
                      f"{row['patch_norm']:>6} {row['copies']:>3} "
                      f"{row['vertices']:>5} {row['edges']:>6} "
                      f"{row['cross_copy_edges']:>6} {row['chi']:>4} "
                      f"{row.get('core_vertices', '-'):>5} {row['seconds']:>7}")
                with open(args.out, "w") as fh:
                    json.dump(rows, fh, indent=2)
    print(f"\nwrote {args.out}  ({len(rows)} rows)")

    best = max(rows, key=lambda r: (r["chi"], -r["vertices"]))
    print(f"max chi seen: {best['chi']}  at n={best['n']} patch={best['patch_norm']} "
          f"copies={best['copies']} |V|={best['vertices']}")


if __name__ == "__main__":
    main()
