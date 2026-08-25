"""Calibration: rediscover a 4-chromatic unit-distance graph from scratch.

Nothing here is copied from a published coordinate list.  We take the
triangular lattice, rotate a patch of it by the angle that carries a
distance-sqrt(3) lattice point to unit distance from itself, and let the
certifier find the edges.  Then we ask what the smallest non-3-colourable
subgraph of the result is.

Run:  python explore/calibrate_spindle.py
"""

from __future__ import annotations

import json
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

import colouring as C  # noqa: E402
from exact_field import field_for  # noqa: E402
from lattice import (  # noqa: E402
    eisenstein_patch,
    rotated_union,
    rotation_for_norm,
)
from unit_distance import build_graph, certify, numeric_separation  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def main():
    n = 3
    field = field_for(3, 4 * n - 1)
    rot = rotation_for_norm(n, field)
    print(f"field            : {field}")
    print(f"rotation         : cos = {rot.c}, sin = {rot.s}")
    print(f"orthogonal (exact): {rot.is_orthogonal()}")

    patch = eisenstein_patch(3)
    print(f"lattice patch    : {len(patch)} points, norms <= 3")
    points, tags = rotated_union(field, patch, rot, copies=2)
    print(f"union of 2 copies: {len(points)} distinct points")

    edges, stats = build_graph(points)
    print(f"graph            : {stats}")

    rep = certify(points, edges)
    print(f"certificate      : ok={rep['ok']} pairs_checked={rep['pairs_checked']}")
    assert rep["ok"], rep

    sep = numeric_separation(points, edges, prec=80)
    print(f"numeric check    : max |d^2-1| over edges  = {sep['max_edge_residual'][0]}")
    print(f"                   min |d^2-1| over non-edges = "
          f"{sep['min_nonedge_gap'][0]}  at {sep['min_nonedge_gap'][1]}")

    chi = C.chromatic_number(len(points), edges)
    print(f"chi(union)       : {chi}")
    col = C.k_colouring(len(points), edges, chi)
    print(f"colouring proper : {C.verify_colouring(len(points), edges, col)}")

    core = C.vertex_critical_subgraph(len(points), edges, 3)
    assert core is not None, "union is 3-colourable"
    kept, sub = core
    print(f"4-critical core  : {len(kept)} vertices, {len(sub)} edges")
    print(f"core vertices    : {[tags[i] for i in kept]}")

    # three independent verdicts on the core
    chi_bt = C.chromatic_number(len(kept), sub)
    chi_bf = next(k for k in range(len(kept) + 1)
                  if C.brute_force_k(len(kept), sub, k) is not None)
    chi_poly = C.chromatic_number_via_polynomial(len(kept), sub)
    poly = C.chromatic_polynomial(len(kept), sub)
    print(f"core chi         : dsatur={chi_bt} bruteforce={chi_bf} polynomial={chi_poly}")
    print(f"core P(G,k)      : {poly}")
    print(f"core P(G,3)      : {C.poly_eval(poly, 3)}   P(G,4) = {C.poly_eval(poly, 4)}")

    degs = sorted(len(a) for a in C.adjacency(len(kept), sub))
    print(f"core degrees     : {degs}")

    # the core as an independently certifiable unit-distance graph
    core_pts = [points[i] for i in kept]
    core_edges, core_stats = build_graph(core_pts)
    core_rep = certify(core_pts, core_edges)
    print(f"core certificate : ok={core_rep['ok']} edges={core_stats['edges']}")
    assert set(map(tuple, map(sorted, core_edges))) == set(map(tuple, map(sorted, sub)))

    os.makedirs(OUT, exist_ok=True)
    record = {
        "construction": "Eisenstein patch (norm <= 3) union its image under "
                        "rotation with cos = 5/6, sin = sqrt11/6",
        "field": repr(field),
        "union_vertices": len(points),
        "union_edges": len(edges),
        "union_chi": chi,
        "core_vertices": len(kept),
        "core_edges": len(sub),
        "core_chi": {"dsatur": chi_bt, "brute_force": chi_bf, "polynomial": chi_poly},
        "core_chromatic_polynomial": poly,
        "core_degree_sequence": degs,
        "core_lattice_tags": [list(map(list, [tags[i]])) [0] for i in kept],
        "core_coordinates_exact": [
            {"x": [str(c) for c in points[i].x.coords()],
             "y": [str(c) for c in points[i].y.coords()]}
            for i in kept
        ],
        "coordinate_basis": ["1", "sqrt3", "sqrt11", "sqrt33"],
        "core_edge_list": [list(e) for e in sub],
    }
    path = os.path.join(OUT, "calibration_spindle.json")
    with open(path, "w") as fh:
        json.dump(record, fh, indent=2, default=str)
    print(f"wrote            : {path}")


if __name__ == "__main__":
    main()
