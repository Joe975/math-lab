"""Annealing inside the forced pair model.

Combines the two attacks: the pair model fixes everything about the graph
except the adjacency within the distance-2 layer D, and annealing searches only
that.  Compared with srg_anneal.py this
  * starts from a state where v0, N(v0) and every P-D edge are already exactly
    right (they are forced, not guessed), and
  * searches C(|D|,2) variables instead of C(n,2).

The moves are 2-swaps restricted to D x D, so every D-vertex keeps exactly
k-2 neighbours inside D, which is its forced D-degree.

Calibration matters more here than anywhere else: the same code is run on
k=4 (n=9) and k=22 (n=243), where a graph provably exists.  The propagating
search in pair_search.py fails that control; if this one passes it, results on
k=14 mean something.

Witnesses are verified with harness/conway-99/srg.py::check_srg, which shares
no code with the energy function.
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
from pair_search import Model


def popcount(x):
    return bin(x).count("1")


class PairSearch:
    """The 99-vertex graph, with only D x D free."""

    def __init__(self, k, rng):
        self.M = Model(k)
        self.k, self.lam, self.mu = k, 1, 2
        self.n = self.M.n
        self.rng = rng
        self.nd = self.M.nd
        self.degD = self.M.deg_d
        self.offset = 1 + k                      # D vertices start here
        self._build_seed()
        self._random_D()
        self.energy = self._full_energy()

    def _build_seed(self):
        """v0, P with its perfect matching, and every P-D edge (all forced)."""
        n, k, M = self.n, self.k, self.M
        self.adj = [0] * n

        def link(a, b):
            self.adj[a] |= 1 << b
            self.adj[b] |= 1 << a

        for x in range(k):
            link(0, 1 + x)
        for x in range(0, k, 2):
            link(1 + x, 1 + (x ^ 1))
        for i, p in enumerate(M.pairs):
            for x in p:
                link(1 + x, self.offset + i)

    def _random_D(self):
        """A random (k-2)-regular graph on D, as a circulant then swaps."""
        nd, deg = self.nd, self.degD
        off = self.offset
        self.dedges = set()
        for v in range(nd):
            for d in range(1, deg // 2 + 1):
                a, b = v, (v + d) % nd
                if a == b:
                    continue
                self.adj[off + a] |= 1 << (off + b)
                self.adj[off + b] |= 1 << (off + a)
                self.dedges.add((min(a, b), max(a, b)))
        if deg % 2:
            for v in range(nd):
                a, b = v, (v + nd // 2) % nd
                self.adj[off + a] |= 1 << (off + b)
                self.adj[off + b] |= 1 << (off + a)
                self.dedges.add((min(a, b), max(a, b)))
        self.dedges = sorted(self.dedges)
        # randomise
        m = len(self.dedges)
        for _ in range(20 * m):
            self._try_swap(self.rng.randrange(m), self.rng.randrange(m),
                           track=False)

    # -- energy ----------------------------------------------------------
    def _target(self, u, v):
        return self.lam if (self.adj[u] >> v & 1) else self.mu

    def _pair_dev(self, u, v):
        d = popcount(self.adj[u] & self.adj[v]) - self._target(u, v)
        return d * d

    def _full_energy(self):
        e = 0
        for u in range(self.n):
            au = self.adj[u]
            for v in range(u + 1, self.n):
                d = popcount(au & self.adj[v]) - self._target(u, v)
                e += d * d
        return e

    def _flip(self, a, b):
        self.adj[a] ^= 1 << b
        self.adj[b] ^= 1 << a

    def _affected(self, a, b):
        s = {(min(a, b), max(a, b))}
        na, nb = self.adj[a], self.adj[b]
        for w in range(self.n):
            if w != a and (nb >> w & 1):
                s.add((min(a, w), max(a, w)))
            if w != b and (na >> w & 1):
                s.add((min(b, w), max(b, w)))
        return s

    def toggle(self, a, b):
        aff = self._affected(a, b)
        self._flip(a, b)
        aff |= self._affected(a, b)
        self._flip(a, b)
        before = sum(self._pair_dev(u, v) for (u, v) in aff)
        self._flip(a, b)
        after = sum(self._pair_dev(u, v) for (u, v) in aff)
        self.energy += after - before

    def _try_swap(self, ia, ic, track=True):
        a, b = self.dedges[ia]
        c, d = self.dedges[ic]
        if len({a, b, c, d}) != 4:
            return None
        if self.rng.random() < 0.5:
            c, d = d, c
        off = self.offset
        if (self.adj[off + a] >> (off + c) & 1) or (self.adj[off + b] >> (off + d) & 1):
            return None
        if track:
            self.toggle(off + a, off + b)
            self.toggle(off + c, off + d)
            self.toggle(off + a, off + c)
            self.toggle(off + b, off + d)
        else:
            self._flip(off + a, off + b)
            self._flip(off + c, off + d)
            self._flip(off + a, off + c)
            self._flip(off + b, off + d)
        self.dedges[ia] = (min(a, c), max(a, c))
        self.dedges[ic] = (min(b, d), max(b, d))
        return (ia, ic, (a, b), (c, d))

    def undo(self, tok):
        ia, ic, (a, b), (c, d) = tok
        off = self.offset
        self.toggle(off + a, off + c)
        self.toggle(off + b, off + d)
        self.toggle(off + a, off + b)
        self.toggle(off + c, off + d)
        self.dedges[ia] = (min(a, b), max(a, b))
        self.dedges[ic] = (min(c, d), max(c, d))


def anneal(k, seed=0, iters=200000, restarts=6, t1=0.05, deadline=None,
           verbose=False):
    rng = random.Random(seed)
    best = None
    for rs in range(restarts):
        if deadline and time.time() > deadline:
            break
        S = PairSearch(k, rng)
        t0_ = max(2.0, S.energy / (4.0 * S.n))
        cur = S.energy
        m = len(S.dedges)
        for it in range(iters):
            if S.energy == 0:
                break
            if deadline and (it & 2047) == 0 and time.time() > deadline:
                break
            temp = t0_ * ((t1 / t0_) ** (it / iters))
            e0 = S.energy
            tok = S._try_swap(rng.randrange(m), rng.randrange(m))
            if tok is None:
                continue
            de = S.energy - e0
            if de <= 0 or rng.random() < math.exp(-de / temp):
                cur = min(cur, S.energy)
            else:
                S.undo(tok)
        if best is None or cur < best:
            best = cur
        if verbose:
            print(f"    restart {rs}: best {cur} (overall {best})", flush=True)
        if S.energy == 0:
            return {"result": "solved", "adj": [x for x in S.adj]}
    return {"result": "no-witness-found", "best_energy": best}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-k", type=int, default=14)
    ap.add_argument("--iters", type=int, default=200000)
    ap.add_argument("--restarts", type=int, default=6)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--minutes", type=float, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    deadline = time.time() + args.minutes * 60 if args.minutes else None
    M = Model(args.k)
    print(f"# pair-model annealing: k={args.k} n={M.n} |D|={M.nd} "
          f"free variables={M.nd*(M.nd-1)//2}")
    t0 = time.time()
    res = anneal(args.k, seed=args.seed, iters=args.iters,
                 restarts=args.restarts, deadline=deadline, verbose=True)
    dt = time.time() - t0
    print(f"result: {res['result']} in {dt:.1f}s")
    if res["result"] == "solved":
        ok, got = srg.check_srg(res["adj"])
        print(f"  INDEPENDENT CHECK (harness check_srg): {ok} {got}")
        if ok and args.out:
            json.dump({"k": args.k, "n": M.n,
                       "adjacency_bitmasks": res["adj"],
                       "verified": [bool(ok), list(got)]},
                      open(args.out, "w"), indent=1)
            print(f"  wrote {args.out}")
    else:
        print(f"  best energy: {res['best_energy']}")


if __name__ == "__main__":
    main()
