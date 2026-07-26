#!/usr/bin/env python3
"""graceful.py — exact graceful-labeling counting statistics for all trees on n vertices.

Graceful labeling of a tree T on n vertices: injection f: V -> {0..n-1} such
that {|f(u)-f(v)| : uv in E} = {1..n-1}.

COUNTING CONVENTION: for each tree we take one fixed vertex-labeled
representative (as produced by networkx.nonisomorphic_trees) and count
distinct graceful injections f on it, divided by 2 for the complementation
symmetry f -> (n-1)-f (which always yields a distinct graceful labeling).
We do NOT quotient by tree automorphisms; counts of automorphism-rich trees
(paths, stars) therefore include automorphic images. A star's count is
exactly 1 under this convention (center must get label 0 or n-1).

Subcommands:
  validate          run the full validation suite (tree counts vs A000055,
                    brute force vs Python backtracking vs C core, stars)
  run  [--budget-min M] [--nmin A] [--nmax B] [--jobs J] [--outdir D]
                    count graceful labelings for ALL trees on n = A..B,
                    checkpointing a JSON per n; stops early (honestly) when
                    the projected time for the next n exceeds the remaining
                    budget.
  analyze [--outdir D]
                    aggregate the per-n JSON files: count distributions,
                    minimizer/maximizer structure, growth fit.

The heavy counting is done by a small C program (graceful_core.c, compiled
on demand with gcc); a pure-Python counter with identical semantics is used
for validation and as a fallback.
"""

import argparse
import itertools
import json
import math
import os
import statistics
import subprocess
import sys
import time

import networkx as nx
from networkx.generators.nonisomorphic_trees import nonisomorphic_trees

TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
CORE_SRC = os.path.join(TOOLS_DIR, "graceful_core.c")
CORE_BIN = os.path.join(TOOLS_DIR, "graceful_core")
DEFAULT_OUTDIR = os.path.normpath(
    os.path.join(TOOLS_DIR, "..", "attempts", "graceful-trees", "data"))

# A000055: number of trees on n nodes, n = 1, 2, 3, ...
A000055 = [1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159, 7741,
           19320, 48629, 123867, 317955, 823065]


# ---------------------------------------------------------------- tree prep

def bfs_parent_order(G):
    """Return parent array for vertices in a BFS order rooted at a max-degree
    vertex, children visited in degree-descending order (good pruning order).
    par[0] is unused; vertex i>=1 has parent par[i] < i in the new order."""
    deg = dict(G.degree())
    root = max(G.nodes(), key=lambda v: deg[v])
    order = [root]
    pos = {root: 0}
    par = [0] * G.number_of_nodes()
    head = 0
    while head < len(order):
        u = order[head]
        head += 1
        for w in sorted(G.neighbors(u), key=lambda v: -deg[v]):
            if w not in pos:
                pos[w] = len(order)
                par[pos[w]] = pos[u]
                order.append(w)
    assert len(order) == G.number_of_nodes()
    return par


def tree_features(G):
    """Structural features of tree G."""
    n = G.number_of_nodes()
    degs = sorted((d for _, d in G.degree()), reverse=True)
    maxdeg = degs[0] if degs else 0
    leaves = sum(1 for d in degs if d == 1)
    diam = nx.diameter(G) if n > 1 else 0

    def is_path_graph(H):
        if H.number_of_nodes() == 0:
            return True
        return (nx.is_connected(H)
                and max((d for _, d in H.degree()), default=0) <= 2)

    # caterpillar: removing all leaves yields a (possibly empty) path
    H1 = G.subgraph([v for v, d in G.degree() if d > 1]).copy()
    caterpillar = is_path_graph(H1)
    # lobster: removing all leaves yields a caterpillar
    if caterpillar:
        lobster = True
    else:
        H2 = H1.subgraph([v for v, d in H1.degree() if d > 1]).copy()
        lobster = is_path_graph(H2)
    # spider: at most one vertex of degree >= 3
    spider = sum(1 for d in degs if d >= 3) <= 1
    return {
        "maxdeg": maxdeg,
        "diam": diam,
        "leaves": leaves,
        "cat": caterpillar,
        "lob": lobster,
        "spider": spider,
    }


def gen_trees(n, check=True):
    """All non-isomorphic trees on n vertices (networkx), count-checked
    against OEIS A000055."""
    trees = list(nonisomorphic_trees(n))
    if check and n <= len(A000055):
        expected = A000055[n - 1]
        if len(trees) != expected:
            raise RuntimeError(
                f"tree generation mismatch at n={n}: got {len(trees)}, "
                f"A000055 says {expected}")
    return trees


# ------------------------------------------------------------- counters

def count_graceful_python(G):
    """Count graceful labelings of tree G up to complementation.
    Pure-Python backtracking; same algorithm as the C core."""
    n = G.number_of_nodes()
    if n == 1:
        return 1
    m = n - 1
    par = bfs_parent_order(G)
    lab = [0] * n

    def rec(i, usedL, usedD):
        if i == n:
            return 1
        total = 0
        pl = lab[par[i]]
        for a in range(m + 1):
            if usedL >> a & 1:
                continue
            d = a - pl if a > pl else pl - a
            if usedD >> d & 1:
                continue
            lab[i] = a
            total += rec(i + 1, usedL | 1 << a, usedD | 1 << d)
        return total

    half = 0
    for a in range((m + 1) // 2):     # root labels a with 2a < m
        lab[0] = a
        half += rec(1, 1 << a, 0)
    if m % 2 == 0:
        lab[0] = m // 2
        mid = rec(1, 1 << (m // 2), 0)
        assert mid % 2 == 0, "parity violation"
        half += mid // 2
    return half


def count_graceful_bruteforce(G):
    """Count graceful labelings up to complementation by checking all n!
    injections. Only for n <= 9. Independent of the backtracking code path."""
    n = G.number_of_nodes()
    if n == 1:
        return 1
    nodes = list(G.nodes())
    edges = list(G.edges())
    target = set(range(1, n))
    raw = 0
    for perm in itertools.permutations(range(n)):
        f = dict(zip(nodes, perm))
        if {abs(f[u] - f[v]) for u, v in edges} == target:
            raw += 1
    assert raw % 2 == 0
    return raw // 2


def ensure_core():
    """Compile the C core if needed; return path or None if gcc unavailable."""
    if (os.path.exists(CORE_BIN)
            and os.path.getmtime(CORE_BIN) >= os.path.getmtime(CORE_SRC)):
        return CORE_BIN
    try:
        subprocess.run(
            ["gcc", "-O3", "-march=native", "-o", CORE_BIN, CORE_SRC],
            check=True, capture_output=True)
        return CORE_BIN
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"warning: cannot build C core ({e}); using Python fallback",
              file=sys.stderr)
        return None


def count_batch_core(par_arrays):
    """Count a batch of trees (list of parent arrays) with the C core."""
    inp = "".join(
        f"{len(par)} " + " ".join(map(str, par[1:])) + "\n"
        for par in par_arrays)
    out = subprocess.run([CORE_BIN], input=inp, capture_output=True,
                         text=True, check=True)
    counts = [int(x) for x in out.stdout.split()]
    assert len(counts) == len(par_arrays)
    return counts


def count_all_trees(trees, jobs=1):
    """Counts for a list of trees, using the C core (parallel over `jobs`
    subprocesses) or the Python fallback."""
    pars = [bfs_parent_order(G) for G in trees]
    if ensure_core() is None:
        return [count_graceful_python(G) for G in trees]
    if jobs <= 1 or len(pars) < 4 * jobs:
        return count_batch_core(pars)
    # interleaved chunks so hard trees spread across workers
    chunks = [pars[k::jobs] for k in range(jobs)]
    procs = []
    for chunk in chunks:
        inp = "".join(
            f"{len(par)} " + " ".join(map(str, par[1:])) + "\n"
            for par in chunk)
        p = subprocess.Popen([CORE_BIN], stdin=subprocess.PIPE,
                             stdout=subprocess.PIPE, text=True)
        p.stdin.write(inp)
        p.stdin.close()
        procs.append(p)
    results = [None] * len(pars)
    for k, p in enumerate(procs):
        vals = [int(x) for x in p.stdout.read().split()]
        if p.wait() != 0:
            raise RuntimeError("C core failed")
        assert len(vals) == len(chunks[k])
        for j, v in enumerate(vals):
            results[k + j * jobs] = v
    return results


# ------------------------------------------------------------- validation

def validate():
    print("== validation suite ==")
    ok = True

    print("[1] tree generation counts vs OEIS A000055 (n=1..14):")
    for n in range(1, 15):
        got = len(list(nonisomorphic_trees(n))) if n >= 2 else 1
        exp = A000055[n - 1]
        status = "ok" if got == exp else "MISMATCH"
        if got != exp:
            ok = False
        print(f"    n={n:2d}: {got} (expected {exp}) {status}")

    print("[2] brute force vs Python backtracking vs C core, all trees n=4..8:")
    has_core = ensure_core() is not None
    for n in range(4, 9):
        trees = gen_trees(n)
        bf = [count_graceful_bruteforce(G) for G in trees]
        py = [count_graceful_python(G) for G in trees]
        agree = bf == py
        line = f"    n={n}: {len(trees)} trees, brute==python: {agree}"
        if has_core:
            cc = count_batch_core([bfs_parent_order(G) for G in trees])
            agree_c = bf == cc
            line += f", brute==C: {agree_c}"
            agree = agree and agree_c
        print(line)
        if not agree:
            ok = False
        if min(bf) < 1:
            print(f"    n={n}: TREE WITH 0 LABELINGS (bug!)")
            ok = False

    print("[3] stars K_{1,n-1} must have exactly 1 labeling up to complement:")
    for n in [5, 9, 12, 15, 18]:
        star = nx.star_graph(n - 1)
        c = (count_batch_core([bfs_parent_order(star)])[0] if has_core
             else count_graceful_python(star))
        status = "ok" if c == 1 else "FAIL"
        if c != 1:
            ok = False
        print(f"    n={n:2d}: count={c} {status}")

    print("[4] paths (sanity, counts reported for the record):")
    for n in range(4, 13):
        path = nx.path_graph(n)
        c = (count_batch_core([bfs_parent_order(path)])[0] if has_core
             else count_graceful_python(path))
        print(f"    n={n:2d}: path count={c}")

    print("VALIDATION", "PASSED" if ok else "FAILED")
    return ok


# ------------------------------------------------------------- production run

def run(nmin, nmax, budget_min, jobs, outdir):
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    budget = budget_min * 60.0
    last_elapsed = None
    growth = None
    for n in range(nmin, nmax + 1):
        spent = time.time() - t0
        if last_elapsed is not None:
            # projected cost of this n from observed per-n growth ratio
            proj = last_elapsed * (growth if growth else 8.0)
            if spent + proj > budget:
                print(f"[budget] stopping before n={n}: spent {spent:.0f}s, "
                      f"projected {proj:.0f}s for n={n} exceeds "
                      f"{budget:.0f}s budget")
                break
        tstart = time.time()
        trees = gen_trees(n)
        counts = count_all_trees(trees, jobs=jobs)
        elapsed = time.time() - tstart
        if min(counts) < 1:
            bad = trees[counts.index(0)]
            print(f"!!! n={n}: tree with ZERO graceful labelings: "
                  f"{nx.to_graph6_bytes(bad, header=False).decode().strip()} "
                  f"(almost certainly a bug — conjecture verified to n>=35)")
        recs = []
        for G, c in zip(trees, counts):
            g6 = nx.to_graph6_bytes(G, header=False).decode().strip()
            ft = tree_features(G)
            recs.append([g6, c, ft["maxdeg"], ft["diam"], ft["leaves"],
                         int(ft["cat"]), int(ft["lob"]), int(ft["spider"])])
        out = {
            "n": n,
            "convention": "labelings of fixed representative, up to "
                          "complementation f->(n-1)-f; not quotiented by "
                          "tree automorphisms",
            "fields": ["graph6", "count", "maxdeg", "diam", "leaves",
                       "caterpillar", "lobster", "spider"],
            "n_trees": len(trees),
            "elapsed_sec": round(elapsed, 2),
            "trees": recs,
        }
        path = os.path.join(outdir, f"graceful_counts_n{n:02d}.json")
        with open(path, "w") as fh:
            json.dump(out, fh)
        cs = sorted(counts)
        print(f"[n={n:2d}] {len(trees)} trees in {elapsed:.1f}s  "
              f"min={cs[0]} med={cs[len(cs)//2]} max={cs[-1]}  -> {path}",
              flush=True)
        if last_elapsed is not None and last_elapsed > 0.5:
            growth = elapsed / last_elapsed
        last_elapsed = elapsed


# ------------------------------------------------------------- analysis

def load_all(outdir):
    data = {}
    for fn in sorted(os.listdir(outdir)):
        if fn.startswith("graceful_counts_n") and fn.endswith(".json"):
            with open(os.path.join(outdir, fn)) as fh:
                d = json.load(fh)
            data[d["n"]] = d
    return data


def analyze(outdir):
    data = load_all(outdir)
    if not data:
        print("no data files found in", outdir)
        return

    print("== per-n count distribution (counts up to complementation) ==")
    print(f"{'n':>3} {'#trees':>7} {'min':>10} {'median':>12} {'max':>14} "
          f"{'mean':>14} {'time_s':>8}")
    for n in sorted(data):
        counts = [r[1] for r in data[n]["trees"]]
        print(f"{n:>3} {len(counts):>7} {min(counts):>10} "
              f"{int(statistics.median(counts)):>12} {max(counts):>14} "
              f"{statistics.mean(counts):>14.1f} "
              f"{data[n]['elapsed_sec']:>8.1f}")

    print("\n== minimizers and maximizers per n ==")
    for n in sorted(data):
        if n < 6:
            continue
        recs = data[n]["trees"]
        recs_sorted = sorted(recs, key=lambda r: r[1])
        lo = recs_sorted[:3]
        hi = recs_sorted[-3:]
        print(f"n={n}:")
        for tag, group in [("min", lo), ("max", hi)]:
            for r in group:
                g6, c, md, dm, lv, cat, lob, sp = r
                cls = ("star" if md == len(recs[0]) else
                       "spider" if sp and md >= 3 else
                       "caterpillar" if cat else
                       "lobster" if lob else "other")
                if md == n - 1:
                    cls = "star"
                elif dm == n - 1:
                    cls = "path"
                print(f"  {tag}  count={c:>12}  maxdeg={md:>2} diam={dm:>2} "
                      f"leaves={lv:>2} cat={cat} lob={lob} spider={sp} "
                      f"[{cls}]  {g6}")

    print("\n== structural class averages (largest n available) ==")
    nmax = max(data)
    recs = data[nmax]["trees"]
    classes = {
        "caterpillar": [r[1] for r in recs if r[5]],
        "lobster-not-cat": [r[1] for r in recs if r[6] and not r[5]],
        "spider(maxdeg>=3)": [r[1] for r in recs if r[7] and r[2] >= 3],
        "other": [r[1] for r in recs if not r[6]],
    }
    for name, cs in classes.items():
        if cs:
            print(f"  {name:>18}: {len(cs):>6} trees, "
                  f"min={min(cs)}, median={int(statistics.median(cs))}, "
                  f"max={max(cs)}")

    print("\n== growth of extremes ==")
    ns = sorted(data)
    mins = [min(r[1] for r in data[n]["trees"]) for n in ns]
    meds = [statistics.median(r[1] for r in data[n]["trees"]) for n in ns]
    maxs = [max(r[1] for r in data[n]["trees"]) for n in ns]
    print("  n:      ", ns)
    print("  min:    ", mins)
    print("  median: ", [int(x) for x in meds])
    print("  max:    ", maxs)
    if len(ns) >= 4:
        # fit log(min) ~ a + b*n on the last points where min > 0
        pts = [(n, math.log(v)) for n, v in zip(ns, mins) if v > 0][-6:]
        nbar = sum(p[0] for p in pts) / len(pts)
        ybar = sum(p[1] for p in pts) / len(pts)
        b = (sum((p[0] - nbar) * (p[1] - ybar) for p in pts)
             / sum((p[0] - nbar) ** 2 for p in pts))
        print(f"  log-linear fit of min over last {len(pts)} points: "
              f"min ~ C * {math.exp(b):.3f}^n  (b={b:.3f})")
        ratios = [mins[i + 1] / mins[i] for i in range(len(mins) - 1)
                  if mins[i] > 0]
        print(f"  successive min ratios: "
              f"{[round(r, 2) for r in ratios]}")


# ------------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("validate")
    r = sub.add_parser("run")
    r.add_argument("--nmin", type=int, default=4)
    r.add_argument("--nmax", type=int, default=18)
    r.add_argument("--budget-min", type=float, default=35.0)
    r.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    r.add_argument("--outdir", default=DEFAULT_OUTDIR)
    a = sub.add_parser("analyze")
    a.add_argument("--outdir", default=DEFAULT_OUTDIR)
    args = ap.parse_args()

    if args.cmd == "validate":
        sys.exit(0 if validate() else 1)
    elif args.cmd == "run":
        run(args.nmin, args.nmax, args.budget_min, args.jobs, args.outdir)
    elif args.cmd == "analyze":
        analyze(args.outdir)


if __name__ == "__main__":
    main()
