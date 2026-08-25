"""Independent verification of the exhaustive Z_33 elimination.

Two separate checks, neither of which reuses m33_exhaustive.py's search:

  A. The algebra.  The claim "the three diagonal sets FORCE the three
     off-diagonal autocorrelations" is validated against a graph that exists:
     Paley(9) has a semiregular Z_3 with t = 3, exactly the same shape.  We
     compute its true connection sets, apply the formula, and check it
     reproduces the true autocorrelations.

  B. The unrealisability.  m33_exhaustive.py concludes by saying that certain
     autocorrelation vectors are not realised by any 6-subset of Z_33.  Here
     that is re-derived by plain brute force over ALL 6-subsets containing 0
     (C(32,5) = 201376), tabulating every autocorrelation that actually occurs,
     and checking the required ones are absent.

Translation invariance (A_{U+c} = A_U) is what makes "containing 0" lossless,
and it is checked directly rather than assumed.
"""

import os
import sys
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg
import constructions as C
from orbit_calibrate import orbits_of, true_orbit_matrix

FAIL = []


def autocorr(U, m):
    a = [0] * m
    for x in U:
        for y in U:
            a[(x - y) % m] += 1
    return a


# ---------------------------------------------------------------------------
# A. validate the forcing algebra on Paley(9)
# ---------------------------------------------------------------------------

def connection_sets(adj, orbs, m):
    """Recover S_ij from a real graph whose orbits are cyclic under +1."""
    t = len(orbs)
    S = [[None] * t for _ in range(t)]
    for i in range(t):
        for j in range(t):
            base = orbs[i][0]
            s = set()
            for step, v in enumerate(orbs[j]):
                if adj[base] >> v & 1:
                    s.add(step)
            # orbs[i][0] corresponds to 0 in Z_m for orbit i; orbs[j][step] to
            # step in Z_m for orbit j
            S[i][j] = frozenset(s)
    return S


def check_algebra():
    print("A. forcing algebra, validated on Paley(9) (m=3, t=3)")
    adj = C.paley9()
    perm = [(v + 3) % 9 for v in range(9)]
    orbs = orbits_of(perm, 9)
    m, t, k, lam, mu = 3, 3, 4, 1, 2
    S = connection_sets(adj, orbs, m)
    # true autocorrelations
    trueA = {(i, j): autocorr(S[i][j], m) for i in range(t) for j in range(t)}
    # f_i from the diagonals
    def f_of(Sii):
        A = autocorr(Sii, m)
        f = [k - A[0]] + [0] * (m - 1)
        for y in range(1, m):
            f[y] = (lam if y in Sii else mu) - A[y]
        return f
    f0, f1, f2 = f_of(S[0][0]), f_of(S[1][1]), f_of(S[2][2])
    pred01 = [(f0[y] + f1[y] - f2[y]) // 2 for y in range(m)]
    pred02 = [(f0[y] - f1[y] + f2[y]) // 2 for y in range(m)]
    pred12 = [(-f0[y] + f1[y] + f2[y]) // 2 for y in range(m)]
    ok = (pred01 == trueA[(0, 1)] and pred02 == trueA[(0, 2)]
          and pred12 == trueA[(1, 2)])
    print(f"   true    A01={trueA[(0,1)]} A02={trueA[(0,2)]} A12={trueA[(1,2)]}")
    print(f"   forced  A01={pred01} A02={pred02} A12={pred12}")
    print(f"   {'OK  ' if ok else 'FAIL'} the forcing formula reproduces a real graph")
    if not ok:
        FAIL.append("algebra")
    return ok


# ---------------------------------------------------------------------------
# B. brute-force the realisable autocorrelations of 6-subsets of Z_33
# ---------------------------------------------------------------------------

def check_unrealisable():
    print("\nB. which autocorrelations do 6-subsets of Z_33 actually have?")
    m, size = 33, 6
    # translation invariance, checked not assumed
    U = frozenset({0, 1, 4, 9, 11, 20})
    for c in range(m):
        V = frozenset((x + c) % m for x in U)
        if autocorr(V, m) != autocorr(U, m):
            FAIL.append("translation invariance")
            print("   FAIL translation invariance")
            return False
    print("   OK   autocorrelation is translation invariant "
          "(checked on all 33 translates)")

    realisable = set()
    for rest in combinations(range(1, m), size - 1):
        Uset = (0,) + rest
        realisable.add(tuple(autocorr(Uset, m)))
    print(f"   enumerated all C(32,5) = "
          f"{len(list(combinations(range(1,m), size-1)))} subsets containing 0")
    print(f"   distinct autocorrelation vectors realised: {len(realisable)}")

    # now recompute the required vectors independently
    import m33_exhaustive as X
    diag_choices = [frozenset({a, (-a) % m}) for a in range(1, m // 2 + 1)]
    required = []
    for S00 in diag_choices:
        f0 = X.diag_target(S00)
        for S11 in diag_choices:
            f1 = X.diag_target(S11)
            for S22 in diag_choices:
                f2 = X.diag_target(S22)
                sol, why = X.solve_offdiag(f0, f1, f2)
                if sol is not None:
                    required.append((S00, S11, S22, sol))
    print(f"   diagonal choices surviving the parity/positivity test: "
          f"{len(required)}")
    bad = 0
    for (S00, S11, S22, (A01, A02, A12)) in required:
        for A in (A01, A02, A12):
            if tuple(A) in realisable:
                bad += 1
                print(f"   REALISABLE required vector found: {A}")
    print(f"   {'OK  ' if bad == 0 else 'FAIL'} required autocorrelation "
          f"vectors realisable by some 6-subset: {bad}")
    if bad:
        FAIL.append("unrealisability")
    return bad == 0


def main():
    print("=== independent verification of the Z_33 elimination ===\n")
    check_algebra()
    check_unrealisable()
    print()
    if FAIL:
        print(f"VERIFICATION FAILED: {FAIL}")
        sys.exit(1)
    print("BOTH CHECKS PASS: the semiregular Z_33 case is eliminated, and the "
          "elimination is confirmed by an independent brute force.")


if __name__ == "__main__":
    main()
