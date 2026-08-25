"""Search for a pair of points that every 4-colouring must colour alike.

This is the mechanism that took 3 colours to 4, run one level up.  The Moser
rhombus has two vertices at distance sqrt3 that every *3*-colouring makes
equal; placing two rhombi so those vertices land at distance 1 contradicts
that, which is the spindle.  The same argument at the next level needs a
unit-distance graph H with vertices u, v such that every *4*-colouring of H has
c(u) = c(v).  Given one, rotate H about u until v and its image are at distance
1: both images force the same colour on two adjacent points, so the union is
not 4-colourable and chi(R^2) >= 5.

Two decidable questions per pair, both settled by a colouring solver:

* **forced equal**  -- H + uv is not 4-colourable.  Then every 4-colouring of H
  has c(u) = c(v).  This is the one that spindles.
* **forced unequal** -- H with u, v identified is not 4-colourable.  Then every
  4-colouring has c(u) != c(v): a "virtual edge" at distance |u - v|, which can
  be used as a building block at a distance the plane does not supply.

Both are searched here.  The YES direction is a colouring (self-certifying);
only the NO direction runs the complete search, under a node limit, and an
unresolved pair is reported as such rather than guessed.

Run:
  python problems/hadwiger-nelson/explore/forcing_search.py --level 1 --depth 2
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


def sq_distance(group, level, ca, cb):
    """Exact squared distance between two coefficient tuples at this level."""
    from fractions import Fraction

    diff = tuple((ca[j][0] - cb[j][0], ca[j][1] - cb[j][1]) for j in range(group.m))
    d2 = group.sq_norm(diff)
    return d2 * group.field.rational(Fraction(1, 9 ** level))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--depth", type=int, default=2)
    ap.add_argument("--tries", type=int, default=12)
    ap.add_argument("--steps", type=int, default=120000)
    ap.add_argument("--node-limit", type=int, default=400000)
    ap.add_argument("--max-pairs", type=int, default=0,
                    help="0 = every vertex paired with the origin")
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    group = RhoGroup(2)
    gens = group.vectors_of_norm(9 ** args.level)
    V = ball(group, gens, args.depth)
    E = edges_from_unit_vectors(V, gens)
    n = len(V)
    zero = tuple((0, 0) for _ in range(group.m))
    u = V.index(zero)
    adj = {a for a, b in E if b == u} | {b for a, b in E if a == u}
    print(f"level {args.level} depth {args.depth}: |V|={n} |E|={len(E)} "
          f"generators={len(gens)}; origin at index {u}, degree {len(adj)}")

    out = args.out or os.path.join(
        OUT, f"forcing_L{args.level}_d{args.depth}.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)

    targets = [v for v in range(n) if v != u and v not in adj]
    if args.max_pairs:
        targets = targets[:args.max_pairs]
    print(f"testing {len(targets)} non-adjacent pairs against the origin")

    rows = []
    forced_equal, forced_unequal, unresolved = [], [], []
    t_start = time.time()
    for i, v in enumerate(targets):
        # forced equal?  add the edge and ask for 4 colours
        e_plus = E + [(min(u, v), max(u, v))]
        vp, _, howp = F.decide_k_colourable(
            n, e_plus, 4, tries=args.tries, steps=args.steps,
            node_limit=args.node_limit)
        # forced unequal?  identify the two vertices and ask again
        cn, ce = F.contract(n, E, u, v)
        vc, _, howc = F.decide_k_colourable(
            cn, ce, 4, tries=args.tries, steps=args.steps,
            node_limit=args.node_limit)

        d2 = sq_distance(group, args.level, V[u], V[v])
        rec = {
            "v": v, "coeffs": [list(map(list, V[v]))],
            "sq_distance": str(d2),
            "plus_edge_4colourable": vp, "plus_edge_how": howp,
            "contracted_4colourable": vc, "contracted_how": howc,
        }
        if vp is False:
            rec["FORCED_EQUAL"] = True
            forced_equal.append(rec)
            print(f"  *** FORCED EQUAL at v={v}, d^2={d2} ***", flush=True)
        if vc is False:
            rec["FORCED_UNEQUAL"] = True
            forced_unequal.append(rec)
            print(f"  *** FORCED UNEQUAL (virtual edge) at v={v}, d^2={d2} ***",
                  flush=True)
        if vp is None or vc is None:
            unresolved.append(v)
        rows.append(rec)
        if (i + 1) % 25 == 0:
            print(f"  {i + 1}/{len(targets)} pairs, "
                  f"{len(forced_equal)} forced-equal, "
                  f"{len(forced_unequal)} forced-unequal, "
                  f"{len(unresolved)} unresolved, "
                  f"{time.time() - t_start:.0f}s", flush=True)
            with open(out, "w") as fh:
                json.dump({"level": args.level, "depth": args.depth,
                           "vertices": n, "edges": len(E),
                           "rows": rows}, fh, indent=2)
    with open(out, "w") as fh:
        json.dump({"level": args.level, "depth": args.depth,
                   "vertices": n, "edges": len(E), "rows": rows}, fh, indent=2)
    print(f"\nforced-equal pairs:   {len(forced_equal)}")
    print(f"forced-unequal pairs: {len(forced_unequal)}")
    print(f"unresolved:           {len(unresolved)}")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
