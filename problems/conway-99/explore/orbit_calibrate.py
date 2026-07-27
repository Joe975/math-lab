"""Calibration for the orbit-matrix machinery.

An enumeration of orbit matrices is only useful if it provably contains the
orbit matrix of a graph that exists.  This script:

  1. finds a genuine semiregular automorphism of a realised graph,
  2. reads off its true orbit matrix directly from the adjacency matrix,
  3. checks that matrix satisfies the conditions in orbit_matrix.py, and
  4. checks it appears in the enumerated list.

If step 4 ever fails, the enumerator is over-pruning and nothing it says about
99 can be trusted.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg
import constructions as C
import orbit_matrix as OM


def orbits_of(perm, n):
    seen, orbs = [False] * n, []
    for v in range(n):
        if seen[v]:
            continue
        o, w = [], v
        while not seen[w]:
            seen[w] = True
            o.append(w)
            w = perm[w]
        orbs.append(o)
    return orbs


def true_orbit_matrix(adj, orbs):
    """r_ij = |N(u) & O_j| for u the first element of O_i; verify it is
    independent of the choice of u, which is what makes R well defined."""
    t = len(orbs)
    R = [[0] * t for _ in range(t)]
    for i, Oi in enumerate(orbs):
        for j, Oj in enumerate(orbs):
            vals = {sum(1 for w in Oj if adj[u] >> w & 1) for u in Oi}
            if len(vals) != 1:
                raise ValueError(f"r_{i}{j} not constant on the orbit: {vals}")
            R[i][j] = vals.pop()
    return R


def cayley_shift(n, step):
    """The permutation v -> v+step on Z_n, as a semiregular automorphism of a
    Cayley graph on Z_n (used for Paley 9 via its Z_3 x Z_3 model)."""
    return [(v + step) % n for v in range(n)]


def check(name, adj, perm, m, k, lam, mu, enumerate_too=True):
    n = len(adj)
    ok, params = srg.check_srg(adj)
    assert ok, params
    # confirm perm really is an automorphism
    for u in range(n):
        for v in range(n):
            if (adj[u] >> v & 1) != (adj[perm[u]] >> perm[v] & 1):
                return f"{name}: supplied permutation is NOT an automorphism"
    orbs = orbits_of(perm, n)
    sizes = sorted({len(o) for o in orbs})
    if sizes != [m]:
        return f"{name}: orbits have sizes {sizes}, expected all {m}"
    R = true_orbit_matrix(adj, orbs)
    good, why = OM.check_orbit_matrix(R, [m] * len(orbs), k, lam, mu)
    if not good:
        return f"{name}: the REAL graph's orbit matrix fails our conditions: {why}"
    msg = (f"{name}: n={n} m={m} t={len(orbs)} -- real orbit matrix satisfies "
           f"(S) and (Q); diag={[R[i][i] for i in range(len(orbs))]} "
           f"trace={sum(R[i][i] for i in range(len(orbs)))}")
    if enumerate_too:
        sols, stats = OM.enumerate_semiregular(m, k, lam, mu, n)
        # compare up to simultaneous row/column permutation
        found = any(same_up_to_relabel(R, S) for S in sols)
        msg += (f"\n     enumerated {len(sols)} orbit matrices in "
                f"{stats['nodes']} nodes; real one present: {found}")
        if not found:
            msg += "   <-- ENUMERATOR IS OVER-PRUNING"
    return msg


def same_up_to_relabel(R, S):
    """Is R a simultaneous row+column permutation of S?  t is small here."""
    from itertools import permutations
    t = len(R)
    if t > 9:
        return canon(R) == canon(S)
    for p in permutations(range(t)):
        if all(R[i][j] == S[p[i]][p[j]] for i in range(t) for j in range(t)):
            return True
    return False


def canon(R):
    """Cheap invariant for larger t: sorted rows of the sorted matrix."""
    return sorted(tuple(sorted(row)) for row in R)


def main():
    print("=== orbit-matrix calibration against realised graphs ===")

    # Paley(9): built as a Cayley graph on GF(9); v -> v + 3 is the shift by
    # the element (1,0), which is semiregular of order 3.
    p9 = C.paley9()
    perm = [(v + 3) % 9 for v in range(9)]
    print(" ", check("paley9 / Z3", p9, perm, 3, 4, 1, 2))

    # rook(3) = 3x3 lattice: shifting rows is semiregular of order 3
    r3 = C.rook(3)
    perm = [((v // 3 + 1) % 3) * 3 + v % 3 for v in range(9)]
    print(" ", check("rook3 / Z3", r3, perm, 3, 4, 1, 2))

    # BvLS(243): Cayley on GF(3)^11 / Golay; translation by a fixed coset
    # representative is semiregular of order 3.
    adj, reps = C.bvls243()
    rnum = {r: i for i, r in enumerate(reps)}
    # build the coset-rep lookup once more so we can translate
    code = C.ternary_golay_code()
    idx = {}
    for r in reps:
        for c in code:
            idx[tuple((a + b) % 3 for a, b in zip(r, c))] = r
    shift = tuple([1] + [0] * 10)          # a weight-1 vector
    perm = [0] * 243
    for r in reps:
        img = idx[tuple((a + b) % 3 for a, b in zip(r, shift))]
        perm[rnum[r]] = rnum[img]
    print(" ", check("bvls243 / Z3", adj, perm, 3, 22, 1, 2,
                     enumerate_too=False))

    # a larger semiregular subgroup of BvLS: order 27, giving 9 orbits
    shift27 = tuple([0] * 10 + [1])
    perm2 = [0] * 243
    for r in reps:
        img = idx[tuple((a + b) % 3 for a, b in zip(r, shift27))]
        perm2[rnum[r]] = rnum[img]
    # compose the two shifts to get a bigger cyclic group? both have order 3;
    # instead use a single shift of order 3 on a different coordinate
    print(" ", check("bvls243 / Z3 (other coordinate)", adj, perm2, 3, 22, 1, 2,
                     enumerate_too=False))


if __name__ == "__main__":
    main()
