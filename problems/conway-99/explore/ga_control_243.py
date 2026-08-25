"""The positive control that decides whether a null result at 99 means anything.

BvLS(243,22,1,2) is a Cayley graph on GF(3)^5, so translation by a fixed vector
is an order-3 automorphism acting semiregularly with 81 orbits.  That is the
same shape of assumption used on 99, at the same lam=1, mu=2 parameters, on a
graph that provably exists.

If group_anneal can rediscover a graph here, a null result on 99 is evidence
about 99.  If it cannot, the null is evidence about the search.
"""

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg
import constructions as C
import group_anneal as GA


def bvls_order3_automorphism():
    adj, reps = C.bvls243()
    rnum = {r: i for i, r in enumerate(reps)}
    code = C.ternary_golay_code()
    idx = {}
    for r in reps:
        for c in code:
            idx[tuple((a + b) % 3 for a, b in zip(r, c))] = r
    shift = tuple([1] + [0] * 10)
    perm = [0] * 243
    for r in reps:
        img = idx[tuple((a + b) % 3 for a, b in zip(r, shift))]
        perm[rnum[r]] = rnum[img]
    return adj, perm


def main():
    minutes = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0
    adj, perm = bvls_order3_automorphism()
    n = 243
    ok, params = srg.check_srg(adj)
    print(f"BvLS: {ok} {params}")
    bad = any((adj[u] >> v & 1) != (adj[perm[u]] >> perm[v] & 1)
              for u in range(n) for v in range(n))
    print(f"supplied permutation is an automorphism: {not bad}")
    cyc = {}
    seen = [False] * n
    for v in range(n):
        if seen[v]:
            continue
        c, w = 0, v
        while not seen[w]:
            seen[w] = True
            w = perm[w]
            c += 1
        cyc[c] = cyc.get(c, 0) + 1
    print(f"cycle type: {cyc}")
    orbits = GA.pair_orbits([perm], n)
    print(f"pair orbits: {len(orbits)} (from {n*(n-1)//2} pairs)")
    t0 = time.time()
    res = GA.anneal(n, 22, 1, 2, [perm], seed=9, iters=25000, restarts=12,
                    deadline=time.time() + minutes * 60, verbose=True)
    dt = time.time() - t0
    print(f"result: {res['result']} in {dt:.0f}s")
    if res["result"] == "solved":
        good, got = srg.check_srg(res["adj"])
        print(f"  INDEPENDENT CHECK: {good} {got}")
        print("  POSITIVE CONTROL PASSED: the method rebuilds a graph that "
              "exists at n=243, lam=1, mu=2")
    else:
        print(f"  best energy {res['best_energy']}")
        print("  POSITIVE CONTROL FAILED: any null result at n=99 from this "
              "method is a statement about the search, not about 99")


if __name__ == "__main__":
    main()
