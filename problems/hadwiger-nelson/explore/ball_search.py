"""Is any Cayley ball of the ring Z[zeta][rho] non-4-colourable?

V_L(k) is the set of points of H_k = 3^-k (E + rho E) reachable from the origin
in at most L unit steps.  It is the *maximal* candidate at that radius: any
connected unit-distance graph on ring points, of diameter at most L and
containing the origin, is a subgraph of it.  So a 4-colouring of V_L(k) rules
out every such construction at once, and a failure to 4-colour it is a
5-chromatic unit-distance graph.

The two directions are handled asymmetrically on purpose.  A colouring is
self-certifying, so the YES direction is answered by heuristics and then
verified in linear time.  Only the NO direction runs the complete search, and
only on the 4-core.

Checkpoints after every row, so a long run can be resumed or inspected.

Run:
  python problems/hadwiger-nelson/explore/ball_search.py \
      --levels 0,1,2 --depths 1,2,3 --tries 40
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
from rho_group import RhoGroup, ball, edges_from_unit_vectors  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def run(level, depth, group, tries, steps, node_limit, max_vertices):
    gens = group.vectors_of_norm(9 ** level)
    t0 = time.time()
    V = ball(group, gens, depth)
    if len(V) > max_vertices:
        return {"level": level, "depth": depth, "vertices": len(V),
                "skipped": "too large"}
    E = edges_from_unit_vectors(V, gens)
    t_build = time.time() - t0
    kept, sub, idx, shed = F.k_core(len(V), E, 4)

    t1 = time.time()
    verdict, col, how = F.decide_k_colourable(
        len(V), E, 4, tries=tries, steps=steps, node_limit=node_limit)
    row = {
        "level": level, "depth": depth,
        "generators": len(gens),
        "vertices": len(V), "edges": len(E),
        "mean_degree": round(2 * len(E) / len(V), 2),
        "four_core_vertices": len(kept), "four_core_edges": len(sub),
        "four_colourable": verdict, "decided_by": how,
        "colouring_verified": (F.verify(len(V), E, col, 4)
                               if verdict else None),
        "build_seconds": round(t_build, 2),
        "decide_seconds": round(time.time() - t1, 2),
    }
    return row


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--levels", default="0,1,2")
    ap.add_argument("--depths", default="1,2,3")
    ap.add_argument("--tries", type=int, default=25)
    ap.add_argument("--steps", type=int, default=400000)
    ap.add_argument("--node-limit", type=int, default=3000000)
    ap.add_argument("--max-vertices", type=int, default=60000)
    ap.add_argument("--out", default=os.path.join(OUT, "ball_search.json"))
    args = ap.parse_args()

    group = RhoGroup(2)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    rows = []
    hdr = (f"{'lvl':>4} {'L':>3} {'|V|':>8} {'|E|':>9} {'deg':>6} {'4core':>8} "
           f"{'4col':>7} {'how':>10} {'s':>8}")
    print(hdr)
    print("-" * len(hdr))
    for level in [int(x) for x in args.levels.split(",")]:
        for depth in [int(x) for x in args.depths.split(",")]:
            row = run(level, depth, group, args.tries, args.steps,
                      args.node_limit, args.max_vertices)
            rows.append(row)
            if "skipped" in row:
                print(f"{level:>4} {depth:>3} {row['vertices']:>8} "
                      f"{'-':>9} {'-':>6} {'-':>8} {'skip':>7} {'-':>10} {'-':>8}")
            else:
                print(f"{level:>4} {depth:>3} {row['vertices']:>8} "
                      f"{row['edges']:>9} {row['mean_degree']:>6} "
                      f"{row['four_core_vertices']:>8} "
                      f"{str(row['four_colourable']):>7} {row['decided_by']:>10} "
                      f"{row['build_seconds'] + row['decide_seconds']:>8.1f}",
                      flush=True)
            with open(args.out, "w") as fh:
                json.dump(rows, fh, indent=2)
            if row.get("four_colourable") is False:
                print("\n*** NOT 4-COLOURABLE at level "
                      f"{level}, depth {depth} -- extract and certify a core ***")
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
