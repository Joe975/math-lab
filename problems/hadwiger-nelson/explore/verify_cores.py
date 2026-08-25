"""Independent verification of every 4-chromatic core the census produced.

For each rotation norm n whose two-copy union is 4-chromatic, this extracts a
vertex-critical non-3-colourable core, then subjects it to every check the
verification contract asks for:

  * geometry: every edge has squared distance exactly 1 and every non-edge does
    not, recomputed from the exact coordinates by the certifier;
  * a numeric separation report at 80 digits, so the exact verdict is visibly
    not a knife-edge;
  * non-3-colourability by the DSATUR search AND by the inclusion-exclusion
    oracle (a different algorithm that never builds a colouring) AND, when the
    core is small enough, by unpruned brute force over 3^v and by the
    chromatic polynomial;
  * 4-colourability by shipping an explicit colouring, checked by the linear
    checker.

Writes an evidence file with the exact coordinates of every core, so a reader
can redo all of this without rerunning the search.

Run:  python explore/verify_cores.py
"""

from __future__ import annotations

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
    norms_up_to,
    rotated_union,
    rotation_field,
    rotation_for_norm,
)
from unit_distance import build_graph, certify, numeric_separation  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))

BRUTE_LIMIT = 15      # 3^15 = 14M, a couple of minutes at worst
POLY_EDGE_LIMIT = 26  # deletion/contraction is exponential in the edge count
IE_LIMIT = 22         # 2^22 subsets


def verify_core(points, edges, tags, kept, sub, n, patch_norm):
    v = len(kept)
    core_pts = [points[i] for i in kept]

    rebuilt, stats = build_graph(core_pts)
    rep = certify(core_pts, rebuilt)
    same = sorted(map(tuple, rebuilt)) == sorted(tuple(sorted(e)) for e in sub)
    sep = numeric_separation(core_pts, rebuilt, prec=80)

    verdicts = {}
    t0 = time.time()
    verdicts["dsatur_3colourable"] = C.k_colouring(v, sub, 3) is not None
    verdicts["dsatur_seconds"] = round(time.time() - t0, 3)
    if v <= IE_LIMIT:
        t0 = time.time()
        verdicts["ie_3colourable"] = C.is_k_colourable_ie(v, sub, 3)
        verdicts["ie_4colourable"] = C.is_k_colourable_ie(v, sub, 4)
        verdicts["ie_seconds"] = round(time.time() - t0, 3)
    if v <= BRUTE_LIMIT:
        t0 = time.time()
        verdicts["brute_3colourable"] = C.brute_force_k(v, sub, 3) is not None
        verdicts["brute_seconds"] = round(time.time() - t0, 3)
    if len(sub) <= POLY_EDGE_LIMIT:
        poly = C.chromatic_polynomial(v, sub)
        verdicts["chromatic_polynomial"] = poly
        verdicts["P_at_3"] = C.poly_eval(poly, 3)
        verdicts["P_at_4"] = C.poly_eval(poly, 4)

    # Vertex-criticality, decided by the inclusion-exclusion oracle: deleting
    # any single vertex must leave a 3-colourable graph.  This is what rules out
    # the core containing a smaller 4-chromatic unit-distance graph (a Moser
    # spindle, say) as a subgraph.
    crit = []
    for drop in range(v):
        keep = [u for u in range(v) if u != drop]
        idx = {u: i for i, u in enumerate(keep)}
        e2 = [(idx[a], idx[b]) for a, b in sub if a != drop and b != drop]
        if v - 1 <= IE_LIMIT:
            ok = C.is_k_colourable_ie(v - 1, e2, 3)
        else:
            ok = C.k_colouring(v - 1, e2, 3) is not None
        crit.append(ok)
    verdicts["vertex_critical"] = all(crit)
    verdicts["criticality_method"] = "inclusion-exclusion" if v - 1 <= IE_LIMIT else "dsatur"

    col4 = C.k_colouring(v, sub, 4)
    return {
        "n": n,
        "patch_norm": patch_norm,
        "core_vertices": v,
        "core_edges": len(sub),
        "degree_sequence": sorted(len(a) for a in C.adjacency(v, sub)),
        "geometry_certificate_ok": rep["ok"],
        "rebuilt_edge_set_matches": same,
        "pairs_checked": rep["pairs_checked"],
        "max_edge_residual_80dp": str(sep["max_edge_residual"][0]),
        "min_nonedge_gap_80dp": str(sep["min_nonedge_gap"][0]),
        "verdicts": verdicts,
        "four_colouring": col4,
        "four_colouring_proper": C.verify_colouring(v, sub, col4),
        "lattice_tags": [list(tags[i]) for i in kept],
        "edge_list": [list(e) for e in sub],
        "coordinate_basis": None,   # filled in by caller
        "coordinates": [
            {"x": [str(c) for c in points[i].x.coords()],
             "y": [str(c) for c in points[i].y.coords()]}
            for i in kept
        ],
    }


def main():
    max_n = int(sys.argv[1]) if len(sys.argv) > 1 else 21
    patch_norm = int(sys.argv[2]) if len(sys.argv) > 2 else 21
    results = []
    for n in norms_up_to(max_n):
        field = rotation_field(n)
        rot = rotation_for_norm(n, field)
        patch = eisenstein_patch(patch_norm)
        points, tags = rotated_union(field, patch, rot, copies=2)
        edges, _ = build_graph(points)
        if C.k_colouring(len(points), edges, 3) is not None:
            print(f"n={n:>3}: union is 3-colourable, no core")
            continue
        kept, sub = C.vertex_critical_subgraph(len(points), edges, 3)
        rec = verify_core(points, edges, tags, kept, sub, n, patch_norm)
        basis = ["1"] + [f"sqrt{field._mask_prod[m]}" for m in range(1, field.dim)]
        rec["coordinate_basis"] = basis
        rec["field"] = repr(field)
        results.append(rec)
        v = rec["verdicts"]
        print(f"n={n:>3}: core {rec['core_vertices']}v/{rec['core_edges']}e "
              f"geom_ok={rec['geometry_certificate_ok']} "
              f"dsatur3={v['dsatur_3colourable']} "
              f"ie3={v.get('ie_3colourable')} "
              f"brute3={v.get('brute_3colourable')} "
              f"P(3)={v.get('P_at_3')} "
              f"critical={v.get('vertex_critical')} "
              f"4col={rec['four_colouring_proper']} "
              f"min_nonedge_gap={rec['min_nonedge_gap_80dp'][:12]}")

    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, "verified_cores.json")
    with open(path, "w") as fh:
        json.dump({"patch_norm": patch_norm, "max_n": max_n, "cores": results},
                  fh, indent=2)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
