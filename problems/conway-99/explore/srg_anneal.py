"""Direct annealing search for a strongly regular graph, no symmetry assumed.

A structurally different attack from the orbit-matrix route: instead of
assuming an automorphism and searching the quotient, search the space of
k-regular graphs on n vertices directly.

State      a k-regular graph, held as adjacency bitmasks.
Move       a 2-swap: delete edges (a,b),(c,d) and insert (a,c),(b,d).  This is
           the standard degree-preserving move, so the search never leaves the
           k-regular manifold and degree never enters the objective.
Energy     sum over unordered pairs of (|N(u) & N(v)| - target)^2, where the
           target is lam for adjacent pairs and mu for non-adjacent ones.
           Energy 0 iff the graph is SRG(n,k,lam,mu).

The energy is maintained incrementally: toggling one edge (a,b) can only change
|N(u) & N(v)| for pairs (a,w) with w in N(b) and (b,w) with w in N(a), plus the
pair (a,b) itself whose target flips.  That is O(k) work per toggle instead of
O(n^2).

A witness is verified by harness/conway-99/srg.py::check_srg, an independently
written checker that shares no code with the energy function.  A failure to
find a witness is NOT evidence of nonexistence and is never reported as such.
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


class Search:
    def __init__(self, n, k, lam, mu, rng):
        self.n, self.k, self.lam, self.mu = n, k, lam, mu
        self.rng = rng
        self.adj = [0] * n
        self.edges = []
        self._random_regular()
        self.energy = self._full_energy()

    # -- construction -----------------------------------------------------
    def _random_regular(self):
        """A random k-regular graph.

        Built as a circulant -- connection set {+-1,...,+-(k//2)}, plus n/2 when
        k is odd -- which is k-regular by construction, then randomised by a
        long run of degree-preserving 2-swaps.  The configuration model with
        rejection is useless here: at n=99, k=14 a collision-free pairing
        essentially never occurs, so that approach never returns.
        """
        n, k = self.n, self.k
        if k % 2 and n % 2:
            raise ValueError(f"no k-regular graph with n={n}, k={k} both odd")
        adj = [0] * n
        edges = set()

        def link(a, b):
            if a == b or (adj[a] >> b & 1):
                return
            adj[a] |= 1 << b
            adj[b] |= 1 << a
            edges.add((min(a, b), max(a, b)))

        for v in range(n):
            for d in range(1, k // 2 + 1):
                link(v, (v + d) % n)
        if k % 2:
            for v in range(n):
                link(v, (v + n // 2) % n)
        assert all(popcount(adj[v]) == k for v in range(n)), \
            [popcount(x) for x in adj]
        self.adj = adj
        self.edges = sorted(edges)
        self.energy = 0                     # placeholder; set by caller
        # randomise while preserving regularity
        m = len(self.edges)
        rng = self.rng
        for _ in range(20 * m):
            ia, ic = rng.randrange(m), rng.randrange(m)
            a, b = self.edges[ia]
            c, d = self.edges[ic]
            if len({a, b, c, d}) != 4:
                continue
            if rng.random() < 0.5:
                c, d = d, c
            if (adj[a] >> c & 1) or (adj[b] >> d & 1):
                continue
            self._flip(a, b)
            self._flip(c, d)
            self._flip(a, c)
            self._flip(b, d)
            self.edges[ia] = (min(a, c), max(a, c))
            self.edges[ic] = (min(b, d), max(b, d))

    # -- energy -----------------------------------------------------------
    def _target(self, u, v):
        return self.lam if (self.adj[u] >> v & 1) else self.mu

    def _full_energy(self):
        e = 0
        for u in range(self.n):
            au = self.adj[u]
            for v in range(u + 1, self.n):
                d = popcount(au & self.adj[v]) - self._target(u, v)
                e += d * d
        return e

    def _pair_dev(self, u, v):
        d = popcount(self.adj[u] & self.adj[v]) - self._target(u, v)
        return d * d

    def _affected(self, a, b):
        """Pairs whose deviation can change when edge (a,b) is toggled."""
        s = {(min(a, b), max(a, b))}
        nb = self.adj[b]
        na = self.adj[a]
        for w in range(self.n):
            if w != a and (nb >> w & 1):
                s.add((min(a, w), max(a, w)))
            if w != b and (na >> w & 1):
                s.add((min(b, w), max(b, w)))
        return s

    def _flip(self, a, b):
        self.adj[a] ^= 1 << b
        self.adj[b] ^= 1 << a

    def toggle(self, a, b):
        """Flip edge (a,b), maintaining self.energy incrementally.

        The affected set is taken as the union of the affected sets before and
        after the flip, since N(a) and N(b) themselves change.
        """
        aff = self._affected(a, b)
        self._flip(a, b)
        aff |= self._affected(a, b)
        self._flip(a, b)                      # back to the original state
        before = sum(self._pair_dev(u, v) for (u, v) in aff)
        self._flip(a, b)                      # apply for real
        after = sum(self._pair_dev(u, v) for (u, v) in aff)
        self.energy += after - before

    def swap(self, ia, ic):
        """2-swap on edge indices ia, ic.  Returns an undo token or None."""
        a, b = self.edges[ia]
        c, d = self.edges[ic]
        if len({a, b, c, d}) != 4:
            return None
        if self.rng.random() < 0.5:
            c, d = d, c
        if (self.adj[a] >> c & 1) or (self.adj[b] >> d & 1):
            return None
        self.toggle(a, b)
        self.toggle(c, d)
        self.toggle(a, c)
        self.toggle(b, d)
        self.edges[ia] = (min(a, c), max(a, c))
        self.edges[ic] = (min(b, d), max(b, d))
        return (ia, ic, (a, b), (c, d))

    def undo(self, tok):
        ia, ic, (a, b), (c, d) = tok
        self.toggle(a, c)
        self.toggle(b, d)
        self.toggle(a, b)
        self.toggle(c, d)
        self.edges[ia] = (min(a, b), max(a, b))
        self.edges[ic] = (min(c, d), max(c, d))


def anneal(n, k, lam, mu, seed=0, iters=200000, restarts=5, t0=None, t1=0.05,
           deadline=None, verbose=False):
    rng = random.Random(seed)
    best = None
    for rs in range(restarts):
        if deadline and time.time() > deadline:
            break
        S = Search(n, k, lam, mu, rng)
        if t0 is None:
            t0_ = max(2.0, S.energy / (4.0 * n))
        else:
            t0_ = t0
        cur_best = S.energy
        m = len(S.edges)
        for it in range(iters):
            if S.energy == 0:
                break
            if deadline and (it & 2047) == 0 and time.time() > deadline:
                break
            temp = t0_ * ((t1 / t0_) ** (it / iters))
            ia, ic = rng.randrange(m), rng.randrange(m)
            e0 = S.energy
            tok = S.swap(ia, ic)
            if tok is None:
                continue
            de = S.energy - e0
            if de <= 0 or rng.random() < math.exp(-de / temp):
                if S.energy < cur_best:
                    cur_best = S.energy
            else:
                S.undo(tok)
        if best is None or cur_best < best[0]:
            best = (cur_best, [x for x in S.adj] if S.energy == 0 else None)
        if verbose:
            print(f"    restart {rs}: best {cur_best} (overall {best[0]})",
                  flush=True)
        if S.energy == 0:
            return {"result": "solved", "energy": 0, "adj": [x for x in S.adj]}
    return {"result": "no-witness-found", "best_energy": best[0]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, required=True)
    ap.add_argument("-k", type=int, required=True)
    ap.add_argument("--lam", type=int, required=True)
    ap.add_argument("--mu", type=int, required=True)
    ap.add_argument("--iters", type=int, default=200000)
    ap.add_argument("--restarts", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--minutes", type=float, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    deadline = time.time() + args.minutes * 60 if args.minutes else None
    t0 = time.time()
    res = anneal(args.n, args.k, args.lam, args.mu, seed=args.seed,
                 iters=args.iters, restarts=args.restarts, deadline=deadline,
                 verbose=True)
    dt = time.time() - t0
    print(f"SRG({args.n},{args.k},{args.lam},{args.mu}): {res['result']} "
          f"in {dt:.1f}s")
    if res["result"] == "solved":
        ok, got = srg.check_srg(res["adj"])
        print(f"  INDEPENDENT CHECK (harness check_srg): {ok} {got}")
        if ok and args.out:
            json.dump({"n": args.n, "k": args.k, "lam": args.lam,
                       "mu": args.mu,
                       "adjacency_bitmasks": res["adj"],
                       "verified": [bool(ok), list(got)]},
                      open(args.out, "w"), indent=1)
            print(f"  wrote {args.out}")
    else:
        print(f"  best energy reached: {res['best_energy']}")


if __name__ == "__main__":
    main()
