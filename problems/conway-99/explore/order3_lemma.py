"""L5: the fixed-point structure of an order-3 automorphism, and what it kills.

Setting: Gamma is an SRG(99,14,1,2), tau in Aut(Gamma) has order 3, and
F = Fix(tau) with f = |F|.

Step 1.  F is closed under common neighbours.
    For u,w in F the set N(u) & N(w) has size 1 or 2 and is tau-invariant;
    tau has odd order, so it fixes a set of size <= 2 pointwise.  Hence every
    common neighbour of two fixed vertices is fixed.

Step 2.  Every degree in F is 2 or 8.
    For v in F, tau permutes the 7 matching edges of N(v) in orbits of size 1
    or 3, so the number of fixed edges is congruent to 7 = 1 (mod 3), i.e.
    1, 4 or 7.  A tau-fixed edge is fixed setwise, and tau has odd order, so
    both its endpoints are fixed.  7 fixed edges would fix N(v) pointwise and
    force tau = id by L1.  So |F & N(v)| is 2 or 8.

Step 3.  A counting identity.  For any induced subgraph F closed under common
    neighbours, and any u in F,

        sum_{x in N_F(u)} deg_F(x) = 2f - 2

    Proof: count pairs (x,w) with x in N_F(u), w in N_F(x), w != u.  By x that
    is sum (deg_F(x) - 1).  By w it is sum over w != u of |N_F(u) & N_F(w)|,
    which by closure equals the Gamma-value, 1 if u~w and 2 otherwise, giving
    deg_F(u) + 2(f - 1 - deg_F(u)).  Equate and rearrange.

Step 4.  The case analysis.
    * If some u has degree 2, then 2f - 2 <= 2*8 = 16, so f <= 9.
    * If some u has degree 8, then 2f - 2 >= 8*2 = 16, so f >= 9.
    * If both degrees occur then f = 9, every degree-2 vertex has both
      neighbours of degree 8, and every degree-8 vertex has all eight
      neighbours of degree 2.  But lambda = 1 puts every edge in a triangle,
      and the third vertex of a triangle on a (deg 8, deg 2) edge would have to
      have degree 2 and degree 8 at once.  So F is regular.
    * F regular of degree 2: 2*2 = 2f - 2 gives f = 3, a triangle.
    * F regular of degree 8: 8*8 = 2f - 2 gives f = 33, and F is then an
      SRG(33,8,1,2) -- a parameter set the harness rejects (its multiplicities
      are not integers).  Impossible.

    Hence f is 0 or 3.

Consequences (each pure arithmetic on top of L5):
    * No automorphism of order 27.  Orbits of sigma of order 27 have size
      1, 3, 9 or 27; every point in an orbit of size <= 9 is fixed by
      sigma^9, which has order 3, so those points number 0 or 3.  Neither
      99 - 0 nor 99 - 3 is divisible by 27.
    * An automorphism of order 9 is semiregular with 11 orbits of size 9.
      Points in orbits of size 1 or 3 are fixed by sigma^3 (order 3), so they
      number 0 or 3; 99 - 3 = 96 is not divisible by 9, so the count is 0.
    * No automorphism of order 33.  Orbits have size 11 or 33 (a size 1 or 3
      orbit would give an order-11 element with a fixed point, contradicting
      L2), so 99 = 33a + 11b.  Points in 11-orbits are fixed by sigma^11 of
      order 3, and number 11b, which must be 0 or 3, forcing b = 0, a = 3 --
      the semiregular case, eliminated exhaustively in m33_exhaustive.py.
"""

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg
import constructions as C
from autos import automorphisms, order_of

FAIL = []


def counting_identity_holds(adj, F, lam, mu):
    """Check the general identity, for every u in F:

        sum_{x in N_F(u)} deg_F(x) = deg_F(u) * (1 + lam - mu) + mu * (f - 1)

    For (lam, mu) = (1, 2) the deg_F(u) coefficient vanishes and this collapses
    to the 2f - 2 used in L5.  Written in general so the check is meaningful on
    graphs outside the family too.
    """
    Fs = set(F)
    degF = {u: sum(1 for w in Fs if adj[u] >> w & 1) for u in Fs}
    f = len(Fs)
    for u in Fs:
        s = sum(degF[x] for x in Fs if adj[u] >> x & 1)
        want = degF[u] * (1 + lam - mu) + mu * (f - 1)
        if s != want:
            return False, (u, s, want)
    return True, None


def closed_under_common_neighbours(adj, F, n):
    Fs = set(F)
    for u in Fs:
        for w in Fs:
            if u >= w:
                continue
            for x in range(n):
                if (adj[u] >> x & 1) and (adj[w] >> x & 1) and x not in Fs:
                    return False
    return True


def check_on_graph(name, adj, params):
    """For every order-3 automorphism, verify the structural claims."""
    n, k, lam, mu = params
    auts = automorphisms(adj)
    order3 = [p for p in auts if order_of(p) == 3]
    print(f"  {name}: |Aut| = {len(auts)}, order-3 elements: {len(order3)}")
    sizes = set()
    for p in order3:
        F = [v for v in range(n) if p[v] == v]
        sizes.add(len(F))
        if not F:
            continue
        if not closed_under_common_neighbours(adj, F, n):
            print(f"    FAIL: Fix is not closed under common neighbours")
            FAIL.append(name)
            continue
        ok, why = counting_identity_holds(adj, F, lam, mu)
        if not ok:
            print(f"    FAIL: counting identity broken at {why}")
            FAIL.append(name)
    print(f"    fixed-point counts observed: {sorted(sizes)}")
    if lam == 1 and mu == 2:
        print(f"    (lam=1,mu=2: L5's analogue for k={k} predicts "
              f"{'f = 0 only' if k == 4 else 'f in {0,3}'})")
    return sizes


def main():
    print("=== L5 verification ===\n")
    print("1. the counting identity and closure, on realised graphs")
    # lam=1, mu=2 members
    check_on_graph("paley9 (9,4,1,2)", C.paley9(), (9, 4, 1, 2))
    check_on_graph("rook3 (9,4,1,2)", C.rook(3), (9, 4, 1, 2))
    # a different parameter set, to show the identity is not special to one
    check_on_graph("petersen (10,3,0,1)", C.petersen(), (10, 3, 0, 1))
    check_on_graph("rook4 (16,6,2,2)", C.rook(4), (16, 6, 2, 2))

    print("\n2. the arithmetic the case analysis rests on")
    ok, rep = srg.feasibility(33, 8, 1, 2)
    bad = srg.failed_conditions(rep)
    print(f"   SRG(33,8,1,2) feasible: {ok}  (fails: {bad})")
    if ok:
        FAIL.append("(33,8,1,2) should be infeasible")
    # the two regular solutions of d*d = 2f-2 with d in {2,8}
    for d in (2, 8):
        f = (d * d + 2) // 2
        print(f"   F regular of degree {d} -> f = {f}"
              + ("  (a triangle)" if d == 2 else "  (needs SRG(33,8,1,2))"))
    print(f"   so f in {{0, 3}}")

    print("\n3. consequences, as exact divisibility checks")
    for base, name in ((27, "order 27"), (9, "order 9"), (33, "order 33")):
        pass
    print(f"   order 27: 99 - 0 = 99 divisible by 27? {99 % 27 == 0}; "
          f"99 - 3 = 96 divisible by 27? {96 % 27 == 0}  -> no order-27 element")
    print(f"   order  9: 99 - 3 = 96 divisible by 9? {96 % 9 == 0}; "
          f"99 - 0 = 99 divisible by 9? {99 % 9 == 0}  -> fixed-point free, "
          f"{99 // 9} orbits of 9")
    print(f"   order 33: 99 = 33a + 11b; points in 11-orbits number 11b which "
          f"must be 0 or 3 -> b = 0 (11b = 3 has no integer solution: "
          f"{3 % 11 == 0}), so a = 3, the semiregular case")

    print()
    if FAIL:
        print(f"VERIFICATION FAILURES: {FAIL}")
        sys.exit(1)
    print("L5 checks pass.")


if __name__ == "__main__":
    main()
