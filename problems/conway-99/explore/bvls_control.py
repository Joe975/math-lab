"""Positive control for the order-7 computation, at the same shape.

The Berlekamp-van Lint-Seidel graph SRG(243, 22, 1, 2) is the coset graph of
the (cyclic) ternary Golay code, so the cyclic shift of coordinates is an
automorphism of order 11 = k/2. By the argument of 001's L3 it has exactly
one fixed vertex (the zero coset), the other 242 vertices form 22 orbits of
11, N(v0) is two whole orbits joined by the matching, and the orbit matrix
carries the same forced seed as the order-7 case of the 99-graph:

    r_01 = r_02 = 11, r_10 = r_20 = 1, r_11 = r_22 = 0, r_12 = r_21 = 1.

This script builds that automorphism, reads the real orbit matrix off the
graph, checks it satisfies (S), (R), (Q) and the seed, and then asks
orbit7.py -- with the same prunes used at 99 -- to find it. If the
enumeration with symmetry breaking does not contain it, nothing orbit7.py
says about 99 can be trusted.

    python problems/conway-99/explore/bvls_control.py [--node-budget N]
"""

import argparse
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import constructions as C  # noqa: E402
import srg  # noqa: E402
import orbit7  # noqa: E402
import char_blocks  # noqa: E402


def build():
    adj, reps = C.bvls243()
    idx = {r: i for i, r in enumerate(reps)}
    # cyclic shift of coordinates preserves weight, so shifted reps are reps
    perm = [idx[tuple(r[-1:] + r[:-1])] for r in reps]
    assert char_blocks.is_automorphism(adj, perm), "shift is not an automorphism"
    ok, params = srg.check_srg(adj)
    assert ok and params == (243, 22, 1, 2), params
    orbs = char_blocks.orbits_of(perm, 243)
    sizes = [len(o) for o in orbs]
    assert sizes.count(1) == 1 and sizes.count(11) == 22, sorted(sizes)
    # put the fixed point first, then the two N(v0) orbits, then the rest
    v0 = next(o[0] for o in orbs if len(o) == 1)
    nb = [o for o in orbs if len(o) == 11 and adj[v0] >> o[0] & 1]
    rest = [o for o in orbs if len(o) == 11 and not adj[v0] >> o[0] & 1]
    assert len(nb) == 2 and len(rest) == 20
    orbs = [[v0]] + nb + rest
    R = char_blocks.orbit_matrix(adj, orbs)
    return adj, perm, orbs, R


def seed_for(t):
    seed = {(0, 0): 0}
    for j in range(1, t):
        seed[(0, j)] = 11 if j in (1, 2) else 0
    seed[(1, 1)] = 0
    seed[(1, 2)] = 1
    seed[(2, 2)] = 0
    return seed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--node-budget", type=int, default=None)
    ap.add_argument("--spectral", action="store_true",
                    help="also run with the projector prune at the real quotient's a")
    args = ap.parse_args()
    adj, perm, orbs, R = build()
    sizes = [len(o) for o in orbs]
    k, lam, mu = 22, 1, 2
    ok, why = orbit7.check_orbit_matrix(R, sizes, k, lam, mu)
    print(f"real orbit matrix satisfies (S)(R)(Q): {ok} {why}")
    seed = seed_for(len(sizes))
    bad = {ij: (R[ij[0]][ij[1]], v) for ij, v in seed.items() if R[ij[0]][ij[1]] != v}
    print(f"seed entries match the real matrix: {not bad} {bad}")
    print(f"diagonal of the real matrix: {[R[i][i] for i in range(len(sizes))]}")
    print(f"trace {sum(R[i][i] for i in range(len(sizes)))}, "
          f"max off-diagonal among 11-orbits "
          f"{max(R[i][j] for i in range(1, 23) for j in range(1, 23) if i != j)}")
    target = orbit7.canonical(R, sizes, fixed=[0, 1, 2])
    t0 = time.time()
    found, stats = orbit7.enumerate_orbit_matrices(
        sizes, k, lam, mu, seed=seed, symmetry=True, node_budget=args.node_budget,
        limit=None)
    classes = {orbit7.canonical(M, sizes, fixed=[0, 1, 2]) for M in found}
    hit = target in classes
    print(f"orbit7 with symmetry: {len(found)} solutions, {len(classes)} classes, "
          f"nodes {stats['nodes']}, exhausted {stats['exhausted']}, "
          f"rows {stats.get('row_completions')}, {time.time() - t0:.1f}s")
    print(f"REAL ORBIT MATRIX FOUND: {hit}")
    if args.spectral:
        disc = (lam - mu) ** 2 + 4 * (k - mu)
        r = (lam - mu + int(round(disc ** 0.5))) // 2
        t = len(sizes)
        a = t - char_blocks.rank_q([[R[i][j] - (r if i == j else 0) for j in range(t)]
                                    for i in range(t)])
        print(f"real quotient: multiplicity of {r} is a = {a}")
        found2, stats2 = orbit7.enumerate_orbit_matrices(
            sizes, k, lam, mu, seed=seed, symmetry=True, spectral=a,
            node_budget=args.node_budget)
        classes2 = {orbit7.canonical(M, sizes, fixed=[0, 1, 2]) for M in found2}
        print(f"orbit7 with spectral a={a}: {len(found2)} solutions, "
              f"{len(classes2)} classes, nodes {stats2['nodes']}, "
              f"exhausted {stats2['exhausted']}, found real: {target in classes2}")
    return 0 if hit else 1


if __name__ == "__main__":
    sys.exit(main())
