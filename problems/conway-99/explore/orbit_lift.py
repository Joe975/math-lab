"""Lift an orbit matrix to an actual graph, for a semiregular cyclic group.

Model.  V = Z_m x {0..t-1}, and (x,i) ~ (y,j) iff (y-x) mod m lies in a
connection set S_ij <= Z_m, with

    |S_ij| = r_ij,      S_ji = -S_ij,      S_ii = -S_ii,   0 not in S_ii.

Then z -> z+1 on the first coordinate is a semiregular automorphism of order m,
and R = (r_ij) is exactly the orbit matrix.

The strongly-regular conditions become difference conditions.  The common
neighbours of (0,i) and (y,j) are the (z,l) with z in S_il and z-y in S_jl, so

    c_ij(y) := sum_l #{(a,b) in S_il x S_jl : a - b = y}

and the requirement is, for every (i,j,y) except (i=j, y=0),

    c_ij(y) = lam   if y in S_ij,
    c_ij(y) = mu    otherwise.

(The excluded case is automatic: c_ii(0) = sum_l |S_il| = k.)

Search.  Energy = sum of squared deviations of c_ij(y) from its target; a
solution is exactly energy 0.  Because the target for y in S_ij differs from
the target outside it, the energy is not a simple correlation objective and
plain hill-climbing stalls, so this uses annealing with restarts.

IMPORTANT.  This searches for an EXISTENCE witness.  A witness is verified by
`harness/conway-99/srg.py::check_srg`, which is an independently written
checker working from the definition of a strongly regular graph and sharing no
code with the energy function here.  A *failure* to find a witness proves
nothing whatsoever and is never reported as nonexistence.
"""

import argparse
import json
import os
import random
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------

def neg(s, m):
    return frozenset((-x) % m for x in s)


def negation_orbits(m):
    """The orbits of x -> -x on Z_m \\ {0}: pairs {a,-a}, plus the singleton
    {m/2} when m is even."""
    seen, out = set(), []
    for a in range(1, m):
        if a in seen:
            continue
        b = (-a) % m
        seen.add(a)
        seen.add(b)
        out.append(tuple(sorted({a, b})))
    return out


def energy(S, R, m, t, lam, mu):
    """Sum of squared deviations of c_ij(y) from its target.  0 iff the
    construction is a strongly regular graph with these parameters."""
    e = 0
    for i in range(t):
        for j in range(t):
            c = [0] * m
            for l in range(t):
                A, B = S[i][l], S[j][l]
                for a in A:
                    for b in B:
                        c[(a - b) % m] += 1
            Sij = S[i][j]
            for y in range(m):
                if i == j and y == 0:
                    continue
                target = lam if y in Sij else mu
                d = c[y] - target
                e += d * d
    return e


def random_state(R, m, t, rng):
    """A random assignment respecting |S_ij| = r_ij and the symmetry."""
    S = [[None] * t for _ in range(t)]
    negorbs = negation_orbits(m)
    for i in range(t):
        # diagonal: a symmetric subset of size r_ii avoiding 0.  The available
        # building blocks are the negation orbits {a,-a} (size 2) and, when m
        # is even, the singleton {m/2}.  Hitting the size exactly is a small
        # subset-sum, so it is solved by randomised backtracking rather than
        # greedily -- greedy can miss (e.g. need 2 in Z_4 after taking {2}).
        need = R[i][i]
        pool = list(negorbs)
        rng.shuffle(pool)

        def pick(idx, remaining, acc):
            if remaining == 0:
                return acc
            if idx == len(pool):
                return None
            if len(pool[idx]) <= remaining:
                got = pick(idx + 1, remaining - len(pool[idx]), acc + [pool[idx]])
                if got is not None:
                    return got
            return pick(idx + 1, remaining, acc)

        sel = pick(0, need, [])
        if sel is None:
            raise ValueError(f"no symmetric subset of Z_{m} of size {need} avoiding 0")
        S[i][i] = frozenset(x for o in sel for x in o)
    for i in range(t):
        for j in range(i + 1, t):
            S[i][j] = frozenset(rng.sample(range(m), R[i][j]))
            S[j][i] = neg(S[i][j], m)
    return S


def mutate(S, R, m, t, rng):
    """One symmetry-preserving move; returns a new state."""
    i = rng.randrange(t)
    j = rng.randrange(t)
    T = [row[:] for row in S]
    if i == j:
        negorbs = negation_orbits(m)
        cur = set(S[i][i])
        inside = [o for o in negorbs if set(o) <= cur]
        outside = [o for o in negorbs if not (set(o) & cur)]
        if not inside or not outside:
            return None
        a = rng.choice(inside)
        b = rng.choice(outside)
        if len(a) != len(b):
            return None
        cur -= set(a)
        cur |= set(b)
        T[i][i] = frozenset(cur)
    else:
        cur = set(S[i][j])
        out = [x for x in range(m) if x not in cur]
        if not cur or not out:
            return None
        a = rng.choice(sorted(cur))
        b = rng.choice(out)
        cur.discard(a)
        cur.add(b)
        T[i][j] = frozenset(cur)
        T[j][i] = neg(T[i][j], m)
    return T


def build_graph(S, m, t):
    """The explicit adjacency matrix, as bitmasks, for the independent checker."""
    n = m * t
    idx = lambda x, i: i * m + x
    adj = [0] * n
    for i in range(t):
        for j in range(t):
            for x in range(m):
                for d in S[i][j]:
                    y = (x + d) % m
                    u, v = idx(x, i), idx(y, j)
                    if u != v:
                        adj[u] |= 1 << v
                        adj[v] |= 1 << u
    return adj


# ---------------------------------------------------------------------------
# annealing search
# ---------------------------------------------------------------------------

def anneal(R, m, t, k, lam, mu, seed=0, iters=200000, restarts=20,
           t0=4.0, t1=0.02, verbose=False, deadline=None):
    rng = random.Random(seed)
    best_overall = None
    for rs in range(restarts):
        if deadline and time.time() > deadline:
            break
        try:
            S = random_state(R, m, t, rng)
        except ValueError as e:
            return {"result": "impossible-state", "reason": str(e)}
        e = energy(S, R, m, t, lam, mu)
        best, bestS = e, S
        for it in range(iters):
            if e == 0:
                break
            if deadline and (it & 1023) == 0 and time.time() > deadline:
                break
            temp = t0 * ((t1 / t0) ** (it / iters))
            T = mutate(S, R, m, t, rng)
            if T is None:
                continue
            e2 = energy(T, R, m, t, lam, mu)
            if e2 <= e or rng.random() < pow(2.718281828, -(e2 - e) / temp):
                S, e = T, e2
                if e < best:
                    best, bestS = e, [row[:] for row in S]
        if best_overall is None or best < best_overall[0]:
            best_overall = (best, bestS)
        if verbose:
            print(f"    restart {rs}: best energy {best} "
                  f"(overall {best_overall[0]})", flush=True)
        if best == 0:
            break
    return {"result": "solved" if best_overall[0] == 0 else "no-witness-found",
            "best_energy": best_overall[0],
            "S": [[sorted(x) for x in row] for row in best_overall[1]]}


def verify_witness(S, m, t, n, k, lam, mu):
    """Independent verification: build the graph and run the definition-based
    checker from the harness, which shares no code with `energy`."""
    adj = build_graph([[frozenset(x) for x in row] for row in S], m, t)
    ok, got = srg.check_srg(adj)
    return ok, got, adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrices", type=str, required=True,
                    help="json file written by orbit_matrix.py")
    ap.add_argument("--iters", type=int, default=200000)
    ap.add_argument("--restarts", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--minutes", type=float, default=None)
    ap.add_argument("--out", type=str, default=None)
    args = ap.parse_args()

    data = json.load(open(args.matrices))
    m, t, k = data["m"], data["t"], data["k"]
    lam, mu, n = data["lam"], data["mu"], data["n"]
    print(f"# lifting {data['count']} orbit matrices for m={m} t={t} n={n}")
    deadline = time.time() + args.minutes * 60 if args.minutes else None

    results = []
    for idx, R in enumerate(data["matrices"]):
        if deadline and time.time() > deadline:
            print("  (deadline reached)")
            break
        res = anneal(R, m, t, k, lam, mu, seed=args.seed + idx,
                     iters=args.iters, restarts=args.restarts,
                     verbose=True, deadline=deadline)
        line = f"  matrix {idx}: {res['result']} best_energy={res.get('best_energy')}"
        if res["result"] == "solved":
            ok, got, adj = verify_witness(res["S"], m, t, n, k, lam, mu)
            line += f"  INDEPENDENT CHECK: {ok} {got}"
            res["verified"] = [bool(ok), list(got) if ok else str(got)]
            if ok:
                print(line)
                print("  *** WITNESS FOUND ***")
                if args.out:
                    json.dump({"m": m, "t": t, "R": R, "S": res["S"],
                               "verified": res["verified"]},
                              open(args.out, "w"), indent=1)
                results.append(res)
                break
        print(line, flush=True)
        results.append({"matrix": idx, **{a: b for a, b in res.items() if a != "S"}})
    print(f"# done; {sum(1 for r in results if r.get('result')=='solved')} witnesses")


if __name__ == "__main__":
    main()
