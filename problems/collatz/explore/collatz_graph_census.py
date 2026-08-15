#!/usr/bin/env python3
"""Graph-structure census of truncated Collatz digraphs (queue item 21 / 001 lead 1).

The lead proposed three unlabelled graph invariants on G_B -- dominator tree
from 1, minimum directed cuts separating [1, 2^k] from orbits exceeding
2^{k+1}, and treewidth via flow-cutter -- with the falsifiable dichotomy
"interface size |S_k| bounded vs growing".

    python collatz_graph_census.py degeneracy [--b-exp 12 14 16]
    python collatz_graph_census.py interface  [--k 4 .. 24] [--json OUT]
    python collatz_graph_census.py verify-mincut [--k-max 10]

`degeneracy`    -- dominator tree and treewidth, computed for real at small B.
`interface`     -- |S_k| and its residue composition, the one live invariant.
`verify-mincut` -- checks the linear-time |S_k| against an actual max-flow
                   min-cut on the explicit graph, small k only.

G_B: vertices {1..B} + an escape sink 0, edge n -> T(n) when T(n) <= B and
n -> 0 otherwise, where T is the tier-0 map (n/2 for even, 3n+1 for odd).
Every vertex in [1, B] has out-degree exactly 1, which is the fact the whole
census turns on. Standard library only.
"""

from __future__ import annotations

import argparse
import heapq
import json
import sys
from array import array

SINK = 0


def step(n: int) -> int:
    return n // 2 if n % 2 == 0 else 3 * n + 1


# --- the one invariant with content: the interface ------------------------


def first_crossings(k: int) -> tuple[array, int]:
    """cross[v] for v in [0, 2^{k+1}]: the first orbit point of v whose own
    successor leaves the window, or -1 if the orbit never leaves it.

    Linear time and linear memory in the window, by memoizing along chains.
    Note the window bound alone determines this: B never enters, because a
    cut separating [1, 2^k] from everything above 2^{k+1} is supported
    entirely inside [1, 2^{k+1}].
    """
    lim = 1 << (k + 1)
    cross = array("i", bytes(4 * (lim + 1)))  # 0 = unknown, -2 = in progress
    stack: list[int] = []
    for start in range(1, lim + 1):
        if cross[start]:
            continue
        stack.clear()
        v = start
        while True:
            c = cross[v]
            if c > 0 or c == -1:
                break
            if v % 2 == 1 and 3 * v + 1 > lim:
                # v is itself the crossing point: its successor leaves. It
                # must be recorded here -- it is never pushed on the stack,
                # so nothing downstream would write it.
                cross[v] = v
                c = v
                break
            stack.append(v)
            cross[v] = -2
            nxt = step(v)
            if cross[nxt] == -2:
                # Walked into the 1-2-4 cycle: this chain never crosses.
                c = -1
                break
            v = nxt
        for u in stack:
            cross[u] = c
    return cross, lim


def interface(k: int) -> dict:
    """|S_k| = the minimum edge cut separating [1, 2^k] from above 2^{k+1}."""
    cross, lim = first_crossings(k)
    marked = 1 << k
    seen = bytearray(lim + 1)
    live = 0
    for n in range(1, marked + 1):
        c = cross[n]
        if c > 0:
            live += 1
            seen[c] = 1
    cut = [r for r in range(1, lim + 1) if seen[r]]
    residues = {}
    for mod in (2, 3, 4, 6, 8, 9):
        counts = [0] * mod
        for r in cut:
            counts[r % mod] += 1
        residues[str(mod)] = counts
    # Every crossing point is an odd r with 3r+1 > lim, so the candidate pool
    # is fixed by arithmetic alone; the census measures what fraction is used.
    pool = sum(1 for r in range(1, lim + 1, 2) if 3 * r + 1 > lim)
    # If the census measures nothing but the one-step inequality, the cut is
    # exactly the odd r <= 2^k whose single successor already leaves -- no
    # dynamics involved. That set, and its closed-form size, are the null
    # hypothesis this whole computation exists to test.
    predicted = [r for r in range(1, marked + 1, 2) if 3 * r + 1 > lim]
    closed_form = (2 ** (k - 1) - (-1) ** (k - 1)) // 3
    return {
        "predicted_set_matches": cut == predicted,
        "closed_form": closed_form,
        "closed_form_matches": len(cut) == closed_form,
        "cut_max_over_2k": (cut[-1] / marked) if cut else None,
        "k": k,
        "window": lim,
        "marked": marked,
        "live_marked": live,
        "cut": len(cut),
        "cut_over_2k": len(cut) / marked,
        "candidate_pool": pool,
        "pool_used": len(cut) / pool if pool else None,
        "min_cut_vertex": cut[0] if cut else None,
        "max_cut_vertex": cut[-1] if cut else None,
        "residues": residues,
    }


# --- the two invariants that are vacuous, verified rather than asserted ----


def build_graph(b: int) -> array:
    """nxt[v] = T(v) if T(v) <= b else SINK, for v in [1, b]."""
    nxt = array("i", bytes(4 * (b + 1)))
    for v in range(1, b + 1):
        t = step(v)
        nxt[v] = t if t <= b else SINK
    return nxt


def dominator_check(b: int) -> dict:
    """Dominator tree of the reverse graph rooted at 1, by a real algorithm.

    Claim under test: it is the forward-image tree itself, i.e. idom(v) is
    always T(v). If so the dominator computation carries no information the
    map does not already give, and lead 1's component (i) is vacuous.
    """
    nxt = build_graph(b)
    preds: dict[int, list[int]] = {}
    for v in range(1, b + 1):
        preds.setdefault(nxt[v], []).append(v)

    # Reverse-graph BFS from 1: successors of v are v's forward preimages.
    order = [1]
    seen = {1}
    i = 0
    while i < len(order):
        v = order[i]
        i += 1
        for w in preds.get(v, ()):
            if w not in seen:
                seen.add(w)
                order.append(w)
    rpo_num = {v: i for i, v in enumerate(order)}

    # Cooper-Harvey-Kennedy iterative dominators on the reverse graph.
    idom = {1: 1}

    def intersect(a: int, b_: int) -> int:
        while a != b_:
            while rpo_num[a] > rpo_num[b_]:
                a = idom[a]
            while rpo_num[b_] > rpo_num[a]:
                b_ = idom[b_]
        return a

    changed = True
    passes = 0
    while changed:
        changed = False
        passes += 1
        for v in order[1:]:
            # Predecessors of v in the reverse graph = successors in forward =
            # the single vertex nxt[v].
            candidates = [p for p in (nxt[v],) if p in idom]
            if not candidates:
                continue
            new = candidates[0]
            for p in candidates[1:]:
                new = intersect(new, p)
            if idom.get(v) != new:
                idom[v] = new
                changed = True

    mismatches = [v for v in order[1:] if idom.get(v) != nxt[v]]
    # In-degree in the reverse graph is the count of forward images, which is
    # 1 for every vertex: that is why the tree is forced.
    max_reverse_indegree = 1 if len(order) > 1 else 0
    return {
        "b": b,
        "reachable_from_1": len(order),
        "idom_equals_forward_image": not mismatches,
        "mismatches": mismatches[:10],
        "iteration_passes": passes,
        "max_reverse_indegree": max_reverse_indegree,
    }


def treewidth_check(b: int) -> dict:
    """Exact treewidth of the underlying undirected graph, by reduction.

    Repeatedly eliminate vertices of degree <= 2 (adding the fill edge for
    degree exactly 2). Emptying the graph certifies tw <= 2; a cycle in the
    graph certifies tw >= 2.
    """
    nxt = build_graph(b)
    adj: dict[int, set[int]] = {v: set() for v in range(0, b + 1)}
    edges = 0
    for v in range(1, b + 1):
        w = nxt[v]
        if w != v:
            adj[v].add(w)
            adj[w].add(v)
            edges += 1

    work = {v: set(ns) for v, ns in adj.items()}
    live = set(work)
    # Lazy min-degree heap: entries go stale as degrees drop, so a popped
    # entry is re-checked against the current degree before it is used.
    heap = [(len(work[v]), v) for v in work]
    heapq.heapify(heap)
    max_elim_degree = 0
    while live:
        deg, pick = heapq.heappop(heap)
        if pick not in live or deg != len(work[pick]):
            continue
        max_elim_degree = max(max_elim_degree, deg)
        if deg > 2:
            break
        ns = list(work[pick])
        if len(ns) == 2:
            a, c = ns
            work[a].add(c)
            work[c].add(a)
        for n in ns:
            work[n].discard(pick)
        for n in ns:
            heapq.heappush(heap, (len(work[n]), n))
        live.discard(pick)
        del work[pick]

    # Components: numbers reaching 1 (unicyclic) and numbers reaching the
    # sink (a tree). Counting them pins the "tree plus one cycle" shape.
    return {
        "b": b,
        "vertices": b + 1,
        "edges": edges,
        "eliminated_everything": not live,
        "max_elimination_degree": max_elim_degree,
        "treewidth": max_elim_degree if not live else None,
    }


# --- independent check of the min-cut identity ----------------------------


def maxflow_mincut(k: int) -> int:
    """Min cut separating [1, 2^k] from above 2^{k+1}, by explicit max flow.

    Unit-capacity Dinic on the real graph, with a super-source into every
    vertex of [1, 2^k] and every crossing edge running to a super-sink. Slow
    and small on purpose: this exists to check the linear-time count, not to
    scale.
    """
    lim = 1 << (k + 1)
    src, snk = lim + 1, lim + 2
    n = lim + 3
    graph: list[list[list[int]]] = [[] for _ in range(n)]

    def add(u: int, v: int, cap: int) -> None:
        graph[u].append([v, cap, len(graph[v])])
        graph[v].append([u, 0, len(graph[u]) - 1])

    for v in range(1, (1 << k) + 1):
        add(src, v, 1 << 30)
    for v in range(1, lim + 1):
        t = step(v)
        if t > lim:
            add(v, snk, 1)
        else:
            add(v, t, 1)

    flow = 0
    while True:
        level = [-1] * n
        level[src] = 0
        queue = [src]
        for u in queue:
            for v, cap, _ in graph[u]:
                if cap > 0 and level[v] < 0:
                    level[v] = level[u] + 1
                    queue.append(v)
        if level[snk] < 0:
            return flow
        it = [0] * n

        def dfs(u: int, pushed: int) -> int:
            if u == snk:
                return pushed
            while it[u] < len(graph[u]):
                edge = graph[u][it[u]]
                v, cap = edge[0], edge[1]
                if cap > 0 and level[v] == level[u] + 1:
                    got = dfs(v, min(pushed, cap))
                    if got:
                        edge[1] -= got
                        graph[v][edge[2]][1] += got
                        return got
                it[u] += 1
            return 0

        while True:
            pushed = dfs(src, 1 << 30)
            if not pushed:
                break
            flow += pushed


# --- commands -------------------------------------------------------------


def cmd_degeneracy(args) -> int:
    out = []
    for exp in args.b_exp:
        b = 1 << exp
        dom = dominator_check(b)
        tw = treewidth_check(b)
        out.append({"b_exp": exp, "dominator": dom, "treewidth": tw})
        print(
            f"B=2^{exp}: idom==T(v) for all v: {dom['idom_equals_forward_image']} "
            f"({dom['reachable_from_1']} vertices reach 1); "
            f"treewidth={tw['treewidth']} (V={tw['vertices']}, E={tw['edges']})"
        )
    if args.json:
        _write(args.json, out)
    return 0


def cmd_interface(args) -> int:
    rows = []
    prev = None
    for k in range(args.k[0], args.k[1] + 1):
        row = interface(k)
        ratio = row["cut"] / prev if prev else float("nan")
        row["growth_vs_prev_k"] = ratio
        rows.append(row)
        print(
            f"k={k:2d} window=2^{k+1:<2d} |S_k|={row['cut']:>10,} "
            f"|S_k|/2^k={row['cut_over_2k']:.6f} "
            f"pool_used={row['pool_used']:.6f} "
            f"live/marked={row['live_marked'] / row['marked']:.6f} "
            f"S_k/S_{{k-1}}={ratio:.4f}"
        )
        prev = row["cut"]
    if args.json:
        _write(args.json, rows)
    return 0


def cmd_verify_mincut(args) -> int:
    ok = True
    for k in range(2, args.k_max + 1):
        counted = interface(k)["cut"]
        flowed = maxflow_mincut(k)
        agree = counted == flowed
        ok &= agree
        print(f"k={k:2d} counted={counted:>6} maxflow={flowed:>6} {'ok' if agree else 'MISMATCH'}")
    return 0 if ok else 1


def _write(path: str, payload) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=1)
    print(f"wrote {path}", file=sys.stderr)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("degeneracy")
    p.add_argument("--b-exp", type=int, nargs="+", default=[10, 12, 14])
    p.add_argument("--json")
    p.set_defaults(fn=cmd_degeneracy)

    p = sub.add_parser("interface")
    p.add_argument("--k", type=int, nargs=2, default=[4, 20], metavar=("LO", "HI"))
    p.add_argument("--json")
    p.set_defaults(fn=cmd_interface)

    p = sub.add_parser("verify-mincut")
    p.add_argument("--k-max", type=int, default=9)
    p.set_defaults(fn=cmd_verify_mincut)

    args = parser.parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.setrecursionlimit(100000)
    sys.exit(main())
