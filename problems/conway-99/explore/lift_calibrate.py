"""End-to-end calibration of the orbit-matrix lifter.

Takes realised graphs with a genuine semiregular cyclic automorphism, reads off
their true orbit matrix, and asks the lifter to reconstruct *a* graph from it.
If the lifter cannot rebuild graphs that exist, a null result from it on 99
means nothing.

Deliberately spans several parameter sets, not just lam=1 mu=2, so the search
is not being tuned to one target.
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
import orbit_matrix as OM
import orbit_lift as OL
from orbit_calibrate import orbits_of, true_orbit_matrix


def calibrate(name, adj, perm, m, params, iters=60000, restarts=12, seed=1):
    n, k, lam, mu = params
    ok, got = srg.check_srg(adj)
    assert ok and got == params, (got, params)
    for u in range(n):
        for v in range(n):
            if (adj[u] >> v & 1) != (adj[perm[u]] >> perm[v] & 1):
                return f"{name}: permutation is not an automorphism"
    orbs = orbits_of(perm, n)
    if sorted({len(o) for o in orbs}) != [m]:
        return f"{name}: not semiregular of order {m}"
    t = len(orbs)
    R = true_orbit_matrix(adj, orbs)
    good, why = OM.check_orbit_matrix(R, [m] * t, k, lam, mu)
    if not good:
        return f"{name}: real orbit matrix fails conditions: {why}"

    t0 = time.time()
    res = OL.anneal(R, m, t, k, lam, mu, seed=seed, iters=iters,
                    restarts=restarts)
    dt = time.time() - t0
    if res["result"] != "solved":
        return (f"{name}: SRG{params} m={m} t={t} -- lifter FAILED to rebuild a "
                f"graph that exists (best energy {res['best_energy']}, {dt:.1f}s)")
    vok, vgot, _ = OL.verify_witness(res["S"], m, t, n, k, lam, mu)
    return (f"{name}: SRG{params} m={m} t={t} -- lifter rebuilt a graph in "
            f"{dt:.1f}s; independent check_srg says {vok} {vgot}")


def main():
    print("=== lifter calibration: can it rebuild graphs that exist? ===")

    # Paley(9), semiregular Z_3
    p9 = C.paley9()
    print(" ", calibrate("paley9", p9, [(v + 3) % 9 for v in range(9)],
                         3, (9, 4, 1, 2)))

    # rook(3), semiregular Z_3 (row shift)
    r3 = C.rook(3)
    perm = [((v // 3 + 1) % 3) * 3 + v % 3 for v in range(9)]
    print(" ", calibrate("rook3", r3, perm, 3, (9, 4, 1, 2)))

    # Petersen SRG(10,3,0,1): the 5-cycle rotation on both pentagon/pentagram
    pet = C.petersen()
    from itertools import combinations
    vs = list(combinations(range(5), 2))
    idx = {v: i for i, v in enumerate(vs)}
    perm = [idx[tuple(sorted(((a + 1) % 5, (b + 1) % 5)))] for (a, b) in vs]
    print(" ", calibrate("petersen", pet, perm, 5, (10, 3, 0, 1)))

    # rook(4) SRG(16,6,2,2), semiregular Z_4 (row shift)
    r4 = C.rook(4)
    perm = [((v // 4 + 1) % 4) * 4 + v % 4 for v in range(16)]
    print(" ", calibrate("rook4", r4, perm, 4, (16, 6, 2, 2)))

    # Clebsch SRG(16,5,0,2) as a Cayley graph on GF(2)^4: no element of order
    # 4, so use the Shrikhande graph on Z4 x Z4 instead, semiregular Z_4
    shr = C.shrikhande()
    perm = [(((v // 4) + 1) % 4) * 4 + v % 4 for v in range(16)]
    print(" ", calibrate("shrikhande", shr, perm, 4, (16, 6, 2, 2)))


if __name__ == "__main__":
    main()
