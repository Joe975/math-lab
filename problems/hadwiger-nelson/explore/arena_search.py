"""Search the arenas the mod-2 capping argument does not close.

``mod2_cap_sweep.py`` finds, for each group of rotated lattices, whether a
mod-2 functional 4-colours the whole group.  Where one exists the group is
closed and searching inside it is provably futile.  Where none exists the group
is a *candidate arena*, and this script does the actual object hunt there:
build Cayley balls, and ask a complete solver for a 4-colouring.

It also re-tests each candidate against larger quotients (E/NE)^s for several
N, since a capping colouring at some other modulus would close the arena
without any search at all.  That check is run first, because it is cheap and
can save the expensive one.

Run:
  python problems/hadwiger-nelson/explore/arena_search.py \
      --rotations "4;3,4;4,9" --depths 2,3,4
"""

from __future__ import annotations

import argparse
import itertools
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
from mod2_cap_sweep import capping_functional, check_direct_sum  # noqa: E402
from multi_rot_group import (  # noqa: E402
    MultiRotGroup,
    ball,
    edges_from_unit_vectors,
)

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def quotient_cap(uv, s, N, k=4, max_size=40000):
    """Is the Cayley graph on (E/NE)^s k-colourable?  If so the group is closed."""
    size = N ** (2 * s)
    if size > max_size:
        return None, "quotient too large"

    def enc(c):
        idx = 0
        for a, b in c:
            idx = (idx * N + a % N) * N + b % N
        return idx

    edges = set()
    for flat in itertools.product(*([range(N)] * (2 * s))):
        c = tuple((flat[2 * j], flat[2 * j + 1]) for j in range(s))
        i = enc(c)
        for u in uv:
            d = tuple((c[j][0] + u[j][0], c[j][1] + u[j][1]) for j in range(s))
            j2 = enc(d)
            if i == j2:
                return None, "self-loop"
            edges.add((min(i, j2), max(i, j2)))
    edges = sorted(edges)
    # A quotient that is *not* k-colourable proves nothing (the quotient has
    # extra edges), so an exhaustive refutation here would be wasted work.
    # Bound it and report the timeout honestly.
    verdict, col, how = F.decide_k_colourable(
        size, edges, k, tries=6, steps=60000, node_limit=200000)
    if verdict and F.verify(size, edges, col, k):
        return (size, len(edges)), "CAPPED"
    if verdict is None:
        return None, "undecided under the node limit"
    return None, f"not {k}-colourable"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rotations", default="4;3,4;4,9",
                    help="semicolon-separated rotation sets")
    ap.add_argument("--depths", default="2,3,4")
    ap.add_argument("--patch", type=int, default=9)
    ap.add_argument("--moduli", default="2,3,4,5")
    ap.add_argument("--max-vertices", type=int, default=60000)
    ap.add_argument("--node-limit", type=int, default=4000000)
    ap.add_argument("--out", default=os.path.join(OUT, "arena_search.json"))
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    rows = []
    for spec in args.rotations.split(";"):
        norms = [int(x) for x in spec.split(",")]
        words = [tuple([0] * len(norms))]
        for i in range(len(norms)):
            w = [0] * len(norms)
            w[i] = 1
            words.append(tuple(w))
        group = MultiRotGroup(norms, words)
        s = group.s
        coll, distinct = check_direct_sum(group, 3)
        uv = group.unit_vectors(args.patch)
        mu, nres = capping_functional(uv, s)
        rec = {
            "rotations": norms, "s": s, "field": repr(group.field),
            "unit_vectors": len(uv), "mod2_residues": nres,
            "direct_sum_collisions": coll,
            "mod2_capping_mu": [list(m) for m in mu] if mu else None,
        }
        print(f"\n=== rotations {norms} | field {group.field} | "
              f"{len(uv)} unit vectors | direct-sum collisions {coll} ===",
              flush=True)
        if coll:
            print("  the coefficient map is not injective: this group has "
                  "relations, so a coefficient-wise colouring says nothing "
                  "about the point set.  Skipping.")
            rec["skipped"] = "not a direct sum"
            rows.append(rec)
            continue

        rec["quotient_caps"] = {}
        capped = False
        for N in [int(x) for x in args.moduli.split(",")]:
            t0 = time.time()
            r, msg = quotient_cap(uv, s, N)
            rec["quotient_caps"][N] = msg
            print(f"  quotient (E/{N}E)^{s}: {msg}"
                  + (f"  {r[0]}v/{r[1]}e" if r else "")
                  + f"  ({time.time() - t0:.0f}s)", flush=True)
            if msg == "CAPPED":
                capped = True
                break
        rec["capped_by_quotient"] = capped
        if capped:
            print("  -> arena closed by a periodic colouring; no search needed")
            rows.append(rec)
            with open(args.out, "w") as fh:
                json.dump(rows, fh, indent=2)
            continue

        rec["balls"] = []
        for L in [int(x) for x in args.depths.split(",")]:
            t0 = time.time()
            V = ball(group, uv, L)
            if len(V) > args.max_vertices:
                print(f"  L={L}: |V|={len(V)} exceeds the cap, stopping")
                rec["balls"].append({"depth": L, "vertices": len(V),
                                     "skipped": True})
                break
            E = edges_from_unit_vectors(V, uv)
            kept, sub, idx, shed = F.k_core(len(V), E, 4)
            verdict, col, how = F.decide_k_colourable(
                len(V), E, 4, tries=25, steps=300000,
                node_limit=args.node_limit)
            row = {
                "depth": L, "vertices": len(V), "edges": len(E),
                "mean_degree": round(2 * len(E) / len(V), 2),
                "four_core_vertices": len(kept),
                "four_colourable": verdict, "decided_by": how,
                "colouring_verified": (F.verify(len(V), E, col, 4)
                                       if verdict else None),
                "seconds": round(time.time() - t0, 1),
            }
            rec["balls"].append(row)
            print(f"  L={L}: |V|={len(V):>6} |E|={len(E):>7} "
                  f"deg={row['mean_degree']:>6} 4core={len(kept):>6} "
                  f"4-colourable={verdict} ({how}) "
                  f"({row['seconds']}s)", flush=True)
            with open(args.out, "w") as fh:
                json.dump(rows + [rec], fh, indent=2)
            if verdict is False:
                print("\n  *** NOT 4-COLOURABLE -- extract and certify a core ***")
                break
        rows.append(rec)
        with open(args.out, "w") as fh:
            json.dump(rows, fh, indent=2)

    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
