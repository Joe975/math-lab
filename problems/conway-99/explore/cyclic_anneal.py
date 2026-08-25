"""Annealing over Z_m-invariant graphs, without enumerating orbit matrices.

Model.  V = Z_m x {0..t-1}, (x,i) ~ (y,j) iff (y-x) mod m is in S_ij, with
S_ji = -S_ij and S_ii = -S_ii, 0 not in S_ii.  Every such graph admits the
semiregular automorphism z -> z+1, and every Z_m-invariant graph on m*t
vertices has this form, so the search covers the whole assumed-symmetry class
at once.

Why this rather than orbit_matrix.py + orbit_lift.py: the orbit-matrix
enumeration is a *complete* method but its search tree is enormous for t >= 9
(m = 11 on 99 vertices does not even reach a first solution in minutes).  This
route gives up completeness in exchange for searching the same space directly,
which is the right trade when the goal is to FIND an object rather than to rule
one out.

Objective.  With
    c_ij(y) = sum_l #{(a,b) in S_il x S_jl : a - b = y}
the graph is SRG(n,k,lam,mu) exactly when, for all i,j,y,

    c_ij(y) = k     if i == j and y == 0     (this IS the degree condition)
            = lam   if y in S_ij
            = mu    otherwise

Note the (i,i,0) entry: c_ii(0) = sum_l |S_il| is the degree of a vertex in
orbit i, so regularity needs no separate constraint or penalty.  Energy is the
sum of squared deviations; 0 iff the graph is strongly regular.

Moves toggle a single element of some S_ij (mirrored into S_ji), or a whole
negation orbit {a,-a} on the diagonal.  Sizes are free, so the search is not
confined to one orbit matrix.

Witnesses are checked by harness/conway-99/srg.py::check_srg, written
independently of the energy above.  A null result is not evidence of
nonexistence.
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


class CyclicModel:
    def __init__(self, m, t, k, lam, mu, rng):
        self.m, self.t, self.k, self.lam, self.mu = m, t, k, lam, mu
        self.n = m * t
        self.rng = rng
        self.negorbs = self._negation_orbits()
        self.S = [[set() for _ in range(t)] for _ in range(t)]
        self._random_start()

    def _negation_orbits(self):
        seen, out = set(), []
        for a in range(1, self.m):
            if a in seen:
                continue
            b = (-a) % self.m
            seen.add(a)
            seen.add(b)
            out.append(tuple(sorted({a, b})))
        return out

    def _random_start(self):
        """Random sets with roughly the right total degree."""
        m, t, k = self.m, self.t, self.k
        rng = self.rng
        for i in range(t):
            for j in range(i, t):
                if i == j:
                    for o in self.negorbs:
                        if rng.random() < k / (m * t):
                            self.S[i][i] |= set(o)
                else:
                    s = {x for x in range(m) if rng.random() < k / (m * t)}
                    self.S[i][j] = s
                    self.S[j][i] = {(-x) % m for x in s}

    # -- energy -----------------------------------------------------------
    def corr(self, i, j):
        """c_ij as a length-m vector."""
        c = [0] * self.m
        m = self.m
        for l in range(self.t):
            A, B = self.S[i][l], self.S[j][l]
            if not A or not B:
                continue
            for a in A:
                for b in B:
                    c[(a - b) % m] += 1
        return c

    def energy(self):
        e = 0
        for i in range(self.t):
            for j in range(self.t):
                c = self.corr(i, j)
                Sij = self.S[i][j]
                for y in range(self.m):
                    if i == j and y == 0:
                        tgt = self.k
                    elif y in Sij:
                        tgt = self.lam
                    else:
                        tgt = self.mu
                    d = c[y] - tgt
                    e += d * d
        return e

    def mutate(self):
        """Toggle one element (off-diagonal) or one negation orbit (diagonal).
        Returns an undo token."""
        t, m, rng = self.t, self.m, self.rng
        i = rng.randrange(t)
        j = rng.randrange(t)
        if i == j:
            o = rng.choice(self.negorbs)
            if set(o) <= self.S[i][i]:
                self.S[i][i] -= set(o)
            else:
                self.S[i][i] |= set(o)
            return ("d", i, o)
        x = rng.randrange(m)
        if x in self.S[i][j]:
            self.S[i][j].discard(x)
            self.S[j][i].discard((-x) % m)
        else:
            self.S[i][j].add(x)
            self.S[j][i].add((-x) % m)
        return ("o", i, j, x)

    def undo(self, tok):
        if tok[0] == "d":
            _, i, o = tok
            if set(o) <= self.S[i][i]:
                self.S[i][i] -= set(o)
            else:
                self.S[i][i] |= set(o)
        else:
            _, i, j, x = tok
            if x in self.S[i][j]:
                self.S[i][j].discard(x)
                self.S[j][i].discard((-x) % self.m)
            else:
                self.S[i][j].add(x)
                self.S[j][i].add((-x) % self.m)

    def build_graph(self):
        m, t = self.m, self.t
        n = self.n
        adj = [0] * n
        idx = lambda x, i: i * m + x
        for i in range(t):
            for j in range(t):
                for x in range(m):
                    for d in self.S[i][j]:
                        u, v = idx(x, i), idx((x + d) % m, j)
                        if u != v:
                            adj[u] |= 1 << v
                            adj[v] |= 1 << u
        return adj


def anneal(m, t, k, lam, mu, seed=0, iters=200000, restarts=10, t1=0.05,
           deadline=None, verbose=False):
    rng = random.Random(seed)
    best = None
    for rs in range(restarts):
        if deadline and time.time() > deadline:
            break
        M = CyclicModel(m, t, k, lam, mu, rng)
        e = M.energy()
        t0_ = max(2.0, e / (4.0 * m * t))
        cur = e
        for it in range(iters):
            if e == 0:
                break
            if deadline and (it & 1023) == 0 and time.time() > deadline:
                break
            temp = t0_ * ((t1 / t0_) ** (it / iters))
            tok = M.mutate()
            e2 = M.energy()
            if e2 <= e or rng.random() < math.exp(-(e2 - e) / temp):
                e = e2
                cur = min(cur, e)
            else:
                M.undo(tok)
        if best is None or cur < best:
            best = cur
        if verbose:
            print(f"    restart {rs}: best {cur} (overall {best})", flush=True)
        if e == 0:
            return {"result": "solved", "S": [[sorted(x) for x in row]
                                              for row in M.S],
                    "adj": M.build_graph()}
    return {"result": "no-witness-found", "best_energy": best}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", type=int, required=True)
    ap.add_argument("-n", type=int, default=99)
    ap.add_argument("-k", type=int, default=14)
    ap.add_argument("--lam", type=int, default=1)
    ap.add_argument("--mu", type=int, default=2)
    ap.add_argument("--iters", type=int, default=200000)
    ap.add_argument("--restarts", type=int, default=10)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--minutes", type=float, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    assert args.n % args.m == 0
    t = args.n // args.m
    deadline = time.time() + args.minutes * 60 if args.minutes else None
    print(f"# Z_{args.m}-invariant search: n={args.n} t={t} k={args.k} "
          f"lam={args.lam} mu={args.mu}")
    res = anneal(args.m, t, args.k, args.lam, args.mu, seed=args.seed,
                 iters=args.iters, restarts=args.restarts, deadline=deadline,
                 verbose=True)
    print(f"result: {res['result']}")
    if res["result"] == "solved":
        ok, got = srg.check_srg(res["adj"])
        print(f"  INDEPENDENT CHECK (harness check_srg): {ok} {got}")
        if ok and args.out:
            json.dump({"m": args.m, "t": t, "S": res["S"],
                       "adjacency_bitmasks": res["adj"],
                       "verified": [bool(ok), list(got)]},
                      open(args.out, "w"), indent=1)
            print(f"  wrote {args.out}")
    else:
        print(f"  best energy: {res['best_energy']}")


if __name__ == "__main__":
    main()
