#!/usr/bin/env python3
"""Cycle-spectrum tooling for the Erdos-Gyarfas conjecture.

The conjecture: every graph with minimum degree 3 contains a cycle whose
length is a power of 2.  Cubic graphs are the hard case, so this tool

  1. generates all connected cubic graphs on n vertices (via nauty's geng,
     installed as `nauty-geng` from the Debian/Ubuntu `nauty` package;
     `sudo apt-get install nauty`),
  2. computes each graph's cycle spectrum (the set of cycle lengths present)
     by exact enumeration of all simple cycles, and
  3. checks whether the spectrum meets {4, 8, 16, ...} (powers of 2 that fit
     in an n-vertex graph).

Validation (`--validate`) checks the spectrum code against K4, K3,3 and the
Petersen graph, and cross-checks graph generation two ways: against a fully
independent pure-Python brute-force generator for n <= 8, and against the
known counts of connected cubic graphs (OEIS A002851) for larger n.

Usage:
  python3 cycle_spectrum.py --validate
  python3 cycle_spectrum.py --search --min-n 4 --max-n 16 [--json out.json]

Everything is exact; no randomness anywhere.
"""

import argparse
import json
import shutil
import subprocess
import sys
from itertools import combinations

GENG_CANDIDATES = ["nauty-geng", "geng"]

# OEIS A002851: number of connected cubic graphs on 2n vertices.
KNOWN_CUBIC_COUNTS = {4: 1, 6: 2, 8: 5, 10: 19, 12: 85, 14: 509,
                      16: 4060, 18: 41301, 20: 510489}


# ---------------------------------------------------------------------------
# graph6 parsing (n <= 62 is all we need; geng emits graph6 by default)
# ---------------------------------------------------------------------------

def parse_graph6(line):
    """Parse a graph6 string into an adjacency list. Supports n <= 62."""
    line = line.strip()
    if line.startswith(">>graph6<<"):
        line = line[10:]
    data = [ord(c) - 63 for c in line]
    if any(b < 0 or b > 63 for b in data):
        raise ValueError("invalid graph6 character in %r" % line)
    n = data[0]
    if n == 63:
        raise ValueError("graph6 n >= 63 not supported by this parser")
    bits = []
    for byte in data[1:]:
        for shift in range(5, -1, -1):
            bits.append((byte >> shift) & 1)
    need = n * (n - 1) // 2
    if len(bits) < need:
        raise ValueError("graph6 string too short for n=%d" % n)
    adj = [[] for _ in range(n)]
    k = 0
    for j in range(1, n):        # column-major upper triangle
        for i in range(j):
            if bits[k]:
                adj[i].append(j)
                adj[j].append(i)
            k += 1
    return adj


# ---------------------------------------------------------------------------
# Cycle spectrum by exact enumeration of simple cycles
# ---------------------------------------------------------------------------

def cycle_spectrum(adj):
    """Return the set of cycle lengths present in the graph.

    Enumerates every simple cycle exactly once: a cycle is rooted at its
    minimum vertex, and of the two traversal directions we keep the one whose
    second vertex is smaller than its last vertex.  DFS only visits vertices
    larger than the root.  Fast enough for cubic graphs up to ~n=20.
    """
    n = len(adj)
    spec = set()
    full = frozenset(range(3, n + 1))
    for root in range(n):
        # neighbors of root larger than root
        starts = [v for v in adj[root] if v > root]
        for first in starts:
            # iterative DFS; stack entries: (vertex, visited_mask, depth)
            # depth = number of edges on the path root..vertex
            stack = [(first, (1 << root) | (1 << first), 1)]
            while stack:
                u, visited, depth = stack.pop()
                for v in adj[u]:
                    if v == root:
                        if depth >= 2 and first < u:
                            spec.add(depth + 1)
                    elif v > root and not (visited >> v) & 1:
                        stack.append((v, visited | (1 << v), depth + 1))
        if spec == full:
            return spec
    return spec


def powers_of_two_upto(n):
    p, out = 4, set()          # cycles have length >= 3, so 4 is the first
    while p <= n:
        out.add(p)
        p *= 2
    return out


# ---------------------------------------------------------------------------
# Generation of connected cubic graphs via nauty geng
# ---------------------------------------------------------------------------

def find_geng():
    for name in GENG_CANDIDATES:
        path = shutil.which(name)
        if path:
            return path
    sys.exit("error: nauty's geng not found. Install it, e.g. "
             "`sudo apt-get install nauty` (provides nauty-geng).")


def generate_cubic_g6(n, geng=None, split=None):
    """Yield graph6 strings of all connected cubic graphs on n vertices.

    `split` is an optional "res/mod" string handed to geng to generate only
    the res-th of mod disjoint slices of the output (for parallel runs).
    """
    geng = geng or find_geng()
    args = [geng, "-c", "-d3", "-D3", "-q", str(n)]
    if split:
        args.append(split)
    proc = subprocess.Popen(args,
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        line = line.strip()
        if line:
            yield line
    if proc.wait() != 0:
        raise RuntimeError("geng failed for n=%d" % n)


# ---------------------------------------------------------------------------
# Independent brute-force generator (for cross-checking geng, n <= 8)
# ---------------------------------------------------------------------------

def _brute_labeled_cubic(n):
    """All labeled graphs on [n] with every degree exactly 3 (backtracking)."""
    results = []
    edges = []          # current edge set
    deg = [0] * n

    def place(v):
        if v == n:
            results.append(frozenset(edges))
            return
        need = 3 - deg[v]
        if need < 0:
            return
        if need == 0:
            place(v + 1)
            return
        cands = [u for u in range(v + 1, n) if deg[u] < 3]
        for combo in combinations(cands, need):
            for u in combo:
                deg[u] += 1
                edges.append((v, u))
            deg[v] += need
            place(v + 1)
            deg[v] -= need
            for u in combo:
                deg[u] -= 1
                edges.pop()

    place(0)
    return results


def _is_connected_edges(n, edgeset):
    adj = [[] for _ in range(n)]
    for a, b in edgeset:
        adj[a].append(b)
        adj[b].append(a)
    seen = {0}
    stack = [0]
    while stack:
        u = stack.pop()
        for v in adj[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
    return len(seen) == n


def _edges_to_adj(n, edgeset):
    adj = [[] for _ in range(n)]
    for a, b in edgeset:
        adj[a].append(b)
        adj[b].append(a)
    return adj


def _is_isomorphic(adj1, adj2):
    """Backtracking isomorphism test (all degrees equal; fine for n <= 8)."""
    n = len(adj1)
    if n != len(adj2):
        return False
    nb1 = [set(a) for a in adj1]
    nb2 = [set(a) for a in adj2]
    mapping = [-1] * n
    used = [False] * n

    def extend(v):
        if v == n:
            return True
        for w in range(n):
            if used[w]:
                continue
            ok = True
            for u in nb1[v]:
                if u < v and mapping[u] not in nb2[w]:
                    ok = False
                    break
            if ok:
                # also check non-adjacency of already-mapped vertices
                for u in range(v):
                    if u not in nb1[v] and mapping[u] in nb2[w]:
                        ok = False
                        break
            if ok:
                mapping[v] = w
                used[w] = True
                if extend(v + 1):
                    return True
                used[w] = False
                mapping[v] = -1
        return False

    return extend(0)


def brute_force_connected_cubic_count(n):
    """Number of connected cubic graphs on n vertices, up to isomorphism.

    Fully independent of nauty; usable for n <= 8.
    """
    reps = []
    for edgeset in _brute_labeled_cubic(n):
        if not _is_connected_edges(n, edgeset):
            continue
        adj = _edges_to_adj(n, edgeset)
        spec = tuple(sorted(cycle_spectrum(adj)))
        if not any(spec == s and _is_isomorphic(adj, r)
                   for r, s in reps):
            reps.append((adj, spec))
    return len(reps)


# ---------------------------------------------------------------------------
# Named graphs for validation
# ---------------------------------------------------------------------------

def k4():
    return [[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]]


def k33():
    return [[3, 4, 5], [3, 4, 5], [3, 4, 5], [0, 1, 2], [0, 1, 2], [0, 1, 2]]


def petersen():
    """Outer 5-cycle 0..4, inner pentagram 5..9, spokes i -- i+5."""
    adj = [[] for _ in range(10)]

    def add(a, b):
        adj[a].append(b)
        adj[b].append(a)

    for i in range(5):
        add(i, (i + 1) % 5)          # outer cycle
        add(5 + i, 5 + (i + 2) % 5)  # inner pentagram
        add(i, 5 + i)                # spokes
    return adj


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def run_validate():
    failures = []

    def check(name, got, want):
        status = "ok" if got == want else "FAIL"
        print("  [%s] %s: got %s, expected %s" % (status, name, got, want))
        if got != want:
            failures.append(name)

    print("Spectrum validation on named graphs:")
    check("K4 spectrum", cycle_spectrum(k4()), {3, 4})
    check("K3,3 spectrum", cycle_spectrum(k33()), {4, 6})
    check("Petersen spectrum", cycle_spectrum(petersen()), {5, 6, 8, 9})

    print("Generator cross-check (geng vs independent brute force):")
    geng = find_geng()
    for n in (4, 6, 8):
        geng_count = sum(1 for _ in generate_cubic_g6(n, geng))
        brute = brute_force_connected_cubic_count(n)
        check("n=%d counts (geng=%d, brute=%d)" % (n, geng_count, brute),
              (geng_count, brute), (KNOWN_CUBIC_COUNTS[n],) * 2)

    print("Generator count check vs OEIS A002851:")
    for n in (10, 12, 14):
        geng_count = sum(1 for _ in generate_cubic_g6(n, geng))
        check("n=%d count" % n, geng_count, KNOWN_CUBIC_COUNTS[n])

    print("Spot-check: geng output spectra for n=4 and n=6:")
    g6_n4 = list(generate_cubic_g6(4, geng))
    check("n=4 unique graph is K4",
          cycle_spectrum(parse_graph6(g6_n4[0])), {3, 4})
    specs6 = sorted(tuple(sorted(cycle_spectrum(parse_graph6(g))))
                    for g in generate_cubic_g6(6, geng))
    # The two cubic graphs on 6 vertices are K3,3 {4,6} and K4-prism/Mobius?
    # They are K3,3 and the 3-prism K3xK2 (spectrum {3,4,5,6}).
    check("n=6 spectra", specs6, [(3, 4, 5, 6), (4, 6)])

    if failures:
        print("\nVALIDATION FAILED: %s" % ", ".join(failures))
        return 1
    print("\nAll validations passed.")
    return 0


def run_search(min_n, max_n, json_out=None, near_miss_limit=20, split=None):
    geng = find_geng()
    grand_total = 0
    grand_bad = 0
    report = {}
    for n in range(min_n, max_n + 1, 2):
        pows = powers_of_two_upto(n)
        total = 0
        bad = []                    # counterexamples: no power-of-2 cycle
        inter_hist = {}             # spectrum-intersection pattern -> count
        near_misses = []            # graphs whose only pow2 cycle length is one value
        girth_hist = {}
        for g6 in generate_cubic_g6(n, geng, split=split):
            adj = parse_graph6(g6)
            spec = cycle_spectrum(adj)
            total += 1
            girth = min(spec)
            girth_hist[girth] = girth_hist.get(girth, 0) + 1
            inter = tuple(sorted(spec & pows))
            inter_hist[inter] = inter_hist.get(inter, 0) + 1
            if not inter:
                bad.append((g6, tuple(sorted(spec))))
            elif len(inter) == 1 and inter[0] != 4:
                # survives with a single power-of-2 length, and that length
                # is not the "easy" 4 -- the interesting near misses
                if len(near_misses) < near_miss_limit:
                    near_misses.append(
                        {"g6": g6, "spectrum": sorted(spec),
                         "pow2": list(inter)})
        singles = {k: v for k, v in inter_hist.items() if len(k) == 1}
        report[n] = {
            "graphs": total,
            "counterexamples": [{"g6": g, "spectrum": list(s)}
                                for g, s in bad],
            "powers_considered": sorted(pows),
            "intersection_histogram": {"+".join(map(str, k)) or "NONE": v
                                       for k, v in sorted(inter_hist.items())},
            "girth_histogram": girth_hist,
            "single_pow2_counts": {str(k[0]): v
                                   for k, v in sorted(singles.items())},
            "near_misses_non4": near_misses,
        }
        grand_total += total
        grand_bad += len(bad)
        print("n=%2d: %6d graphs, %d counterexamples; girths %s; "
              "pow2-intersection histogram %s"
              % (n, total, len(bad), girth_hist,
                 report[n]["intersection_histogram"]))
        for m in near_misses:
            print("       near miss (only pow2 length %s): %s spectrum=%s"
                  % (m["pow2"], m["g6"], m["spectrum"]))
    print("\nTOTAL: %d connected cubic graphs checked for n in [%d, %d]; "
          "%d counterexamples."
          % (grand_total, min_n, max_n, grand_bad))
    if grand_bad == 0:
        print("Every graph contains a cycle of length in {4, 8, 16, ...}: "
              "Erdos-Gyarfas holds for all connected cubic graphs in range.")
    if json_out:
        with open(json_out, "w") as f:
            json.dump({"min_n": min_n, "max_n": max_n,
                       "total": grand_total,
                       "counterexample_total": grand_bad,
                       "per_n": report}, f, indent=1)
        print("JSON report written to %s" % json_out)
    return 0 if grand_bad == 0 else 2


def run_block_search(n_list):
    """Search near-cubic 'bridge blocks' avoiding C4 and C8.

    A bridge in a cubic graph splits it into two sides, each a connected
    graph with exactly one degree-2 vertex and the rest degree 3 (n odd,
    (3n-1)/2 edges).  Cycles cannot cross a bridge, so if some side with
    <= 15 vertices avoided both C4 and C8 (C16 impossible at that size),
    two such sides joined by a bridge would be a cubic counterexample.
    This enumerates all C4-free such sides (geng -f) and C8-tests them.
    """
    geng = find_geng()
    any_found = False
    for n in n_list:
        if n % 2 == 0:
            print("n=%d skipped (near-cubic side needs odd n)" % n)
            continue
        e = (3 * n - 1) // 2
        proc = subprocess.run(
            [geng, "-c", "-d2", "-D3", "-f", "-q", str(n), "%d:%d" % (e, e)],
            capture_output=True, text=True, check=True)
        lines = [l for l in proc.stdout.split() if l]
        survivors = []
        for g6 in lines:
            adj = parse_graph6(g6)
            degs = sorted(len(a) for a in adj)
            assert degs.count(2) == 1 and degs.count(3) == n - 1
            spec = cycle_spectrum(adj)
            assert 4 not in spec, "geng -f violated"
            if 8 not in spec:
                survivors.append((g6, sorted(spec)))
        print("n=%2d: %6d C4-free near-cubic sides; C8-free among them: %d"
              % (n, len(lines), len(survivors)))
        for g6, spec in survivors:
            any_found = True
            print("   BLOCK FOUND (no C4, no C8): %s spectrum=%s" % (g6, spec))
    if not any_found:
        print("No C4+C8-free near-cubic side found: any bridged cubic "
              "counterexample needs both sides larger than the range tested.")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--validate", action="store_true",
                    help="run correctness validations and exit")
    ap.add_argument("--search", action="store_true",
                    help="search cubic graphs for power-of-2 cycles")
    ap.add_argument("--blocks", metavar="N", type=int, nargs="+",
                    help="search near-cubic bridge blocks (odd N) that avoid "
                         "C4 and C8")
    ap.add_argument("--min-n", type=int, default=4)
    ap.add_argument("--max-n", type=int, default=14)
    ap.add_argument("--json", metavar="FILE", default=None,
                    help="write JSON report of the search here")
    ap.add_argument("--split", metavar="RES/MOD", default=None,
                    help="process only geng's RES/MOD slice (parallel runs); "
                         "counts then cover only that slice")
    args = ap.parse_args()
    if args.validate:
        sys.exit(run_validate())
    if args.blocks:
        sys.exit(run_block_search(args.blocks))
    if args.search:
        sys.exit(run_search(args.min_n, args.max_n, args.json,
                            split=args.split))
    ap.print_help()
    sys.exit(1)


if __name__ == "__main__":
    main()
