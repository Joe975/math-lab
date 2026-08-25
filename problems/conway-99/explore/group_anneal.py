"""Search for a strongly regular graph invariant under an ASSUMED group.

Given permutations generating a group G on n points, G acts on unordered
vertex pairs.  A G-invariant graph is exactly a union of G-orbits of pairs, so
the search space collapses from C(n,2) binary variables to (number of pair
orbits) -- for n = 99 and |G| = 11 that is 4851 -> 441.

This supersedes both orbit_matrix.py + orbit_lift.py (which need a complete
orbit-matrix enumeration first, and whose search tree is intractable for
t >= 9) and cyclic_anneal.py (which assumes the action is semiregular).  It
handles fixed points, non-cyclic groups, and mixed orbit sizes uniformly,
because it never mentions orbit sizes at all.

State      the set of selected pair-orbits, held as an adjacency matrix.
Move       toggle one pair-orbit (all its edges at once).
Energy     sum over pairs of (|N(u) & N(v)| - target)^2, target = lam if
           adjacent else mu, PLUS sum over vertices of (deg(v) - k)^2 weighted,
           since toggling orbits does not preserve regularity.

Energy 0 iff the graph is SRG(n,k,lam,mu).  Witnesses are checked by
harness/conway-99/srg.py::check_srg, which shares no code with the energy.
A null result is NOT evidence of nonexistence.
"""

import argparse
import json
import math
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg


def popcount(x):
    return bin(x).count("1")


def pair_orbits(gens, n):
    """Orbits of <gens> acting on unordered pairs of {0..n-1}."""
    seen = [[False] * n for _ in range(n)]
    orbits = []
    for u in range(n):
        for v in range(u + 1, n):
            if seen[u][v]:
                continue
            orb = []
            stack = [(u, v)]
            seen[u][v] = seen[v][u] = True
            while stack:
                a, b = stack.pop()
                orb.append((a, b))
                for g in gens:
                    c, d = g[a], g[b]
                    if c > d:
                        c, d = d, c
                    if c != d and not seen[c][d]:
                        seen[c][d] = seen[d][c] = True
                        stack.append((c, d))
            orbits.append(orb)
    return orbits


class GroupSearch:
    def __init__(self, n, k, lam, mu, orbits, rng, degw=4):
        self.n, self.k, self.lam, self.mu = n, k, lam, mu
        self.orbits = orbits
        self.rng = rng
        self.degw = degw
        self.adj = [0] * n
        self.on = [False] * len(orbits)
        # random start: switch on orbits until roughly k-regular
        target_edges = n * k // 2
        order = list(range(len(orbits)))
        rng.shuffle(order)
        cur = 0
        for i in order:
            if cur >= target_edges:
                break
            self._set_orbit(i, True)
            cur += len(orbits[i])
        self.energy = self._full_energy()

    def _set_orbit(self, i, val):
        if self.on[i] == val:
            return
        for (a, b) in self.orbits[i]:
            self.adj[a] ^= 1 << b
            self.adj[b] ^= 1 << a
        self.on[i] = val

    def _full_energy(self):
        e = 0
        for u in range(self.n):
            au = self.adj[u]
            for v in range(u + 1, self.n):
                tgt = self.lam if (au >> v & 1) else self.mu
                d = popcount(au & self.adj[v]) - tgt
                e += d * d
        for v in range(self.n):
            d = popcount(self.adj[v]) - self.k
            e += self.degw * d * d
        return e

    def toggle_orbit(self, i):
        self._set_orbit(i, not self.on[i])
        self.energy = self._full_energy()

    def solved(self):
        return self.energy == 0


def anneal(n, k, lam, mu, gens, seed=0, iters=40000, restarts=8, t1=0.05,
           deadline=None, verbose=False, degw=4):
    rng = random.Random(seed)
    orbits = pair_orbits(gens, n)
    best = None
    for rs in range(restarts):
        if deadline and time.time() > deadline:
            break
        S = GroupSearch(n, k, lam, mu, orbits, rng, degw)
        e = S.energy
        t0_ = max(2.0, e / (4.0 * n))
        cur = e
        for it in range(iters):
            if S.energy == 0:
                break
            if deadline and (it & 255) == 0 and time.time() > deadline:
                break
            temp = t0_ * ((t1 / t0_) ** (it / iters))
            i = rng.randrange(len(orbits))
            e0 = S.energy
            S.toggle_orbit(i)
            de = S.energy - e0
            if de <= 0 or rng.random() < math.exp(-de / temp):
                cur = min(cur, S.energy)
            else:
                S.toggle_orbit(i)
        if best is None or cur < best:
            best = cur
        if verbose:
            print(f"    restart {rs}: best {cur} (overall {best})", flush=True)
        if S.energy == 0:
            return {"result": "solved", "adj": [x for x in S.adj],
                    "n_orbits": len(orbits)}
    return {"result": "no-witness-found", "best_energy": best,
            "n_orbits": len(orbits)}


# ---------------------------------------------------------------------------
# standard generators on 99 points
# ---------------------------------------------------------------------------

def cyclic_gen(n, cycles):
    """A permutation of {0..n-1} built from a list of cycle lengths.

    Points are consumed in order; a cycle of length 1 is a fixed point.
    """
    g = list(range(n))
    p = 0
    for L in cycles:
        for i in range(L):
            g[p + i] = p + (i + 1) % L
        p += L
    assert p == n, (p, n)
    return g


PROFILES_99 = {
    # semiregular, no fixed points
    "z3": [3] * 33,
    "z9": [9] * 11,
    "z11": [11] * 9,
    "z33": [33] * 3,
    # order 7 with exactly one fixed point (forced by blind lemma L3)
    "z7": [1] + [7] * 14,
    # involutions: an odd number of fixed points
    "z2_f1": [1] + [2] * 49,
    "z2_f3": [1] * 3 + [2] * 48,
    "z2_f9": [1] * 9 + [2] * 45,
    "z2_f11": [1] * 11 + [2] * 44,
    # order 3 with fixed points
    "z3_f3": [1] * 3 + [3] * 32,
    "z3_f9": [1] * 9 + [3] * 30,
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", required=True)
    ap.add_argument("-n", type=int, default=99)
    ap.add_argument("-k", type=int, default=14)
    ap.add_argument("--lam", type=int, default=1)
    ap.add_argument("--mu", type=int, default=2)
    ap.add_argument("--iters", type=int, default=40000)
    ap.add_argument("--restarts", type=int, default=8)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--minutes", type=float, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    cycles = PROFILES_99[args.profile]
    assert sum(cycles) == args.n
    g = cyclic_gen(args.n, cycles)
    orbits = pair_orbits([g], args.n)
    print(f"# profile {args.profile}: cycle type {sorted(set(cycles))} "
          f"({len(cycles)} cycles) on {args.n} points")
    print(f"# pair orbits: {len(orbits)} (from C({args.n},2) = "
          f"{args.n*(args.n-1)//2} pairs)")
    deadline = time.time() + args.minutes * 60 if args.minutes else None
    res = anneal(args.n, args.k, args.lam, args.mu, [g], seed=args.seed,
                 iters=args.iters, restarts=args.restarts, deadline=deadline,
                 verbose=True)
    print(f"result: {res['result']}")
    if res["result"] == "solved":
        ok, got = srg.check_srg(res["adj"])
        print(f"  INDEPENDENT CHECK (harness check_srg): {ok} {got}")
        if ok and args.out:
            json.dump({"profile": args.profile, "cycles": cycles,
                       "adjacency_bitmasks": res["adj"],
                       "verified": [bool(ok), list(got)]},
                      open(args.out, "w"), indent=1)
            print(f"  wrote {args.out}")
    else:
        print(f"  best energy: {res['best_energy']}")


if __name__ == "__main__":
    main()
