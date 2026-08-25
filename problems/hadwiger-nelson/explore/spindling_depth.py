"""How far does repeated rotation of one lattice patch go?

Sweeps patch size and number of rotated copies for a fixed rotation rho_n, and
also for unions driven by two different rotations.  Records chi at each step,
which is the quantity of interest: does stacking more rotated copies of the
triangular lattice ever buy a fifth colour, and if not, how does the graph grow
while chi stays put?

Every graph is built and certified exactly; chi is computed by the DSATUR
search and, whenever a non-3-colourable core small enough to fit is extracted,
re-decided by the independent inclusion-exclusion oracle.

Run:  python explore/spindling_depth.py [--ns 3,9,12,21] [--patches 3,4,7,13,21]
                                        [--copies 2,3,4,5]
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
from lattice import (  # noqa: E402
    eisenstein_patch,
    rotated_union,
    rotation_field,
    rotation_for_norm,
)
from unit_distance import build_graph, certify  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def analyse(points, tags, label, core_limit=22, want_core=True):
    t0 = time.time()
    edges, stats = build_graph(points)
    rep = certify(points, edges)
    if not rep["ok"]:
        raise AssertionError(f"certificate failed for {label}: {rep}")
    chi = C.chromatic_number(len(points), edges, lo=1, hi=6)
    col = C.k_colouring(len(points), edges, chi)
    row = {
        "label": label,
        "vertices": len(points),
        "edges": len(edges),
        "cross_copy_edges": sum(1 for a, b in edges if tags[a][0] != tags[b][0]),
        "chi": chi,
        "colouring_verified": C.verify_colouring(len(points), edges, col),
        "certified": True,
        "seconds": round(time.time() - t0, 2),
    }
    if want_core and chi >= 4:
        kept, sub = C.vertex_critical_subgraph(len(points), edges, chi - 1)
        row["core_vertices"] = len(kept)
        row["core_edges"] = len(sub)
        row["core_tags"] = [list(tags[i]) for i in kept]
        row["core_edge_list"] = [list(e) for e in sub]
        if len(kept) <= core_limit:
            row["core_ie_3colourable"] = C.is_k_colourable_ie(len(kept), sub, chi - 1)
            row["core_ie_kcolourable"] = C.is_k_colourable_ie(len(kept), sub, chi)
    row["seconds"] = round(time.time() - t0, 2)
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", default="3,9,12,21")
    ap.add_argument("--patches", default="3,4,7,13,21")
    ap.add_argument("--copies", default="2,3,4,5")
    ap.add_argument("--out", default=os.path.join(OUT, "spindling_depth.json"))
    args = ap.parse_args()

    ns = [int(x) for x in args.ns.split(",")]
    patches = [int(x) for x in args.patches.split(",")]
    copies = [int(x) for x in args.copies.split(",")]

    rows = []
    hdr = f"{'label':>22} {'|V|':>5} {'|E|':>6} {'cross':>6} {'chi':>4} {'core':>5} {'s':>8}"
    print(hdr)
    print("-" * len(hdr))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    for n in ns:
        field = rotation_field(n)
        rot = rotation_for_norm(n, field)
        for m in patches:
            patch = eisenstein_patch(m)
            for c in copies:
                points, tags = rotated_union(field, patch, rot, copies=c)
                label = f"n={n} patch<={m} x{c}"
                row = analyse(points, tags, label)
                row.update({"n": n, "patch_norm": m, "copies": c})
                rows.append(row)
                print(f"{label:>22} {row['vertices']:>5} {row['edges']:>6} "
                      f"{row['cross_copy_edges']:>6} {row['chi']:>4} "
                      f"{row.get('core_vertices', '-'):>5} {row['seconds']:>8}")
                with open(args.out, "w") as fh:
                    json.dump(rows, fh, indent=2)
    print(f"\nwrote {args.out}")
    print(f"max chi seen: {max(r['chi'] for r in rows)}")


if __name__ == "__main__":
    main()
