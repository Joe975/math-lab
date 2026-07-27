"""The partner-regular case for SRG(n, k, 1, 2).

Call a vertex v0 *partner-regular* if, for every x in N(v0), the perfect
matching M_x induced on A_x is the partner matching -- that is, the two
distance-2 vertices {x,y} and {x,y'} are adjacent for every y.

Both realised members of the family (n=9 and n=243) are partner-regular at
every vertex; this is checked in the harness, and is forced there by their
Cayley structure over an elementary abelian group.  So partner-regularity is
the hypothesis that the sought graph resembles the two graphs that exist.

This script:
  1. imposes partner-regularity at v0 for a given k,
  2. runs the exact-count propagator to a fixed point,
  3. reports how much of the search space that settles,
and does it for k=22 (where the answer is known to be satisfiable) as a
control alongside k=14.
"""

import os
import sys
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import pair_search as PS


def partner_matching(M, x):
    """M_x = the partner matching on P \\ {x, x'}."""
    return {y: y ^ 1 for y in range(M.k) if y not in (x, x ^ 1)}


def impose_partner_regular(M, xs=None):
    """Fresh state with M_x = partner matching for every x in xs (default all)."""
    st = PS.State(M)
    q = []
    for x in (range(M.k) if xs is None else xs):
        PS.apply_matching(st, x, partner_matching(M, x), q)
    PS.propagate(st, q)
    return st


def stats(M, st):
    tot = M.nd * (M.nd - 1) // 2
    unk = sum(PS.popcount(st.unknown(v)) for v in range(M.nd)) // 2
    yes = sum(PS.popcount(st.yes[v]) for v in range(M.nd)) // 2
    settled = sum(1 for v in range(M.nd) if not st.unknown(v))
    return {"total_pairs": tot, "decided": tot - unk, "open": unk,
            "yes": yes, "settled_vertices": settled}


def blocks_report(M, st):
    """Under partner-regularity, check the derived block structure:
    a D-vertex's out-of-block neighbours must lie in DISJOINT blocks only."""
    pos = lambda p: p >> 1
    block = [frozenset(pos(p) for p in M.pairs[i]) for i in range(M.nd)]
    bad = 0
    forced_no = 0
    for i in range(M.nd):
        for j in range(M.nd):
            if i == j:
                continue
            if block[i] == block[j]:
                continue
            if block[i] & block[j]:
                s = st.status(i, j)
                if s is True:
                    bad += 1
                elif s is False:
                    forced_no += 1
    return bad, forced_no


def main():
    for k in (4, 14, 22):
        M = PS.Model(k)
        try:
            st = impose_partner_regular(M)
        except PS.Contradiction:
            print(f"k={k:3d} n={M.n:5d}: CONTRADICTION -- partner-regularity is "
                  f"impossible at any vertex")
            continue
        s = stats(M, st)
        bad, fno = blocks_report(M, st)
        print(f"k={k:3d} n={M.n:5d} |D|={M.nd:4d}: "
              f"{s['decided']:6d}/{s['total_pairs']:6d} pairs decided "
              f"({100*s['decided']/s['total_pairs']:.1f}%), "
              f"{s['open']:6d} open, {s['yes']:5d} edges forced, "
              f"{s['settled_vertices']:4d}/{M.nd} vertices settled")
        print(f"          derived block rule: out-of-block neighbours in a "
              f"MEETING block -- asserted true: {bad} (must be 0); "
              f"propagator already forced no: {fno}")
        # remaining degree freedom
        need = [M.deg_d - PS.popcount(st.yes[v]) for v in range(M.nd)]
        cand = [PS.popcount(st.unknown(v)) for v in range(M.nd)]
        print(f"          per-vertex: still needs {sorted(set(need))} more "
              f"neighbours from {sorted(set(cand))} candidates")


if __name__ == "__main__":
    main()
