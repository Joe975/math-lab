"""Which colour patterns can a gadget impose on a handful of its own points?

A unit-distance graph H with a small set T of marked vertices is a *constraint*
on the colours of T: the patterns realizable by proper k-colourings of H.  If
some pattern is missing, H behaves like an extra constraint that the plane does
not otherwise supply, and copies of H can be composed to build something the
raw distance graph cannot express.

Two missing patterns matter most:

* every pattern with c(u) != c(v) missing  =>  u and v are *forced equal*;
  rotating H about u until v meets its own image at distance 1 contradicts
  that, so the union needs k+1 colours.
* every pattern with c(u) == c(v) missing  =>  a *virtual edge* at |u - v|.

Patterns are counted up to permutation of colours, since a proper colouring
composed with a permutation is still proper; the canonical form of a pattern is
its restricted-growth string.

Run:
  python problems/hadwiger-nelson/explore/pattern_profile.py --level 1 --depth 2 \
      --terminals "norm3"
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
from rho_group import RhoGroup, ball, edges_from_unit_vectors  # noqa: E402
from unit_distance import eisenstein_norm  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))


def canonical(pattern):
    """Restricted-growth form: the pattern's class under colour permutation."""
    seen, out = {}, []
    for c in pattern:
        if c not in seen:
            seen[c] = len(seen)
        out.append(seen[c])
    return tuple(out)


def realizable_patterns(n, edges, terminals, k, node_limit=400000,
                        tries=6, steps=60000):
    """Which canonical patterns on ``terminals`` extend to a proper k-colouring.

    Each candidate pattern is imposed by identifying equal terminals and adding
    edges between unequal ones, then asking the solver.  Both moves are exact
    graph surgery, so a YES ships a colouring and a NO is a completed search.
    """
    results = {}
    t = len(terminals)
    for pattern in itertools.product(range(k), repeat=t):
        canon = canonical(pattern)
        if canon != pattern or max(pattern, default=-1) >= k:
            continue
        # build the constrained graph
        cur_n, cur_edges = n, list(edges)
        rep = {v: v for v in range(n)}

        def find(x):
            while rep[x] != x:
                rep[x] = rep[rep[x]]
                x = rep[x]
            return x

        groups = {}
        for term, c in zip(terminals, pattern):
            groups.setdefault(c, []).append(term)
        merged_edges = set((min(a, b), max(a, b)) for a, b in cur_edges)
        # union equal terminals
        for c, members in groups.items():
            for other in members[1:]:
                rep[find(other)] = find(members[0])
        # add edges between differently-coloured terminal groups
        keys = sorted(groups)
        for c1, c2 in itertools.combinations(keys, 2):
            a, b = find(groups[c1][0]), find(groups[c2][0])
            if a != b:
                merged_edges.add((min(a, b), max(a, b)))
        # relabel
        roots = sorted({find(v) for v in range(n)})
        idx = {r: i for i, r in enumerate(roots)}
        final = set()
        self_loop = False
        for a, b in merged_edges:
            aa, bb = idx[find(a)], idx[find(b)]
            if aa == bb:
                self_loop = True
                break
            final.add((min(aa, bb), max(aa, bb)))
        if self_loop:
            results[canon] = {"realizable": False, "how": "self-loop"}
            continue
        verdict, col, how = F.decide_k_colourable(
            len(roots), sorted(final), k, tries=tries, steps=steps,
            node_limit=node_limit)
        results[canon] = {"realizable": verdict, "how": how}
    return results


def pick_terminals(group, V, spec):
    zero = tuple((0, 0) for _ in range(group.m))
    idx = {c: i for i, c in enumerate(V)}
    if spec == "norm3":
        # the origin plus the two lattice points at distance sqrt3 that make
        # the classical rhombus, expressed as coefficient tuples
        picks = [zero]
        for ab in [(1, 1), (-1, 2)]:
            c = (ab, (0, 0))
            if c in idx:
                picks.append(c)
        return [idx[c] for c in picks]
    if spec == "origin-ring":
        picks = [zero]
        for ab in [(1, 1), (-1, 2), (-2, 1)]:
            c = (ab, (0, 0))
            if c in idx:
                picks.append(c)
        return [idx[c] for c in picks]
    raise SystemExit(f"unknown terminal spec {spec}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--level", type=int, default=1)
    ap.add_argument("--depth", type=int, default=2)
    ap.add_argument("--k", type=int, default=4)
    ap.add_argument("--terminals", default="norm3")
    ap.add_argument("--node-limit", type=int, default=600000)
    ap.add_argument("--out", default="")
    args = ap.parse_args()

    group = RhoGroup(2)
    gens = group.vectors_of_norm(9 ** args.level)
    V = ball(group, gens, args.depth)
    E = edges_from_unit_vectors(V, gens)
    terminals = pick_terminals(group, V, args.terminals)
    print(f"level {args.level} depth {args.depth}: |V|={len(V)} |E|={len(E)}; "
          f"terminals {terminals} ({args.terminals})")

    t0 = time.time()
    res = realizable_patterns(len(V), E, terminals, args.k,
                              node_limit=args.node_limit)
    missing = [p for p, r in res.items() if r["realizable"] is False]
    unknown = [p for p, r in res.items() if r["realizable"] is None]
    for p, r in sorted(res.items()):
        print(f"  pattern {p}: realizable={r['realizable']} ({r['how']})")
    print(f"\nmissing patterns:   {missing}")
    print(f"unresolved:         {unknown}")
    print(f"{time.time() - t0:.1f}s")

    out = args.out or os.path.join(
        OUT, f"pattern_L{args.level}_d{args.depth}_{args.terminals}.json")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out, "w") as fh:
        json.dump({"level": args.level, "depth": args.depth, "k": args.k,
                   "terminals": args.terminals, "vertices": len(V),
                   "edges": len(E),
                   "patterns": {str(p): r for p, r in res.items()}},
                  fh, indent=2)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
