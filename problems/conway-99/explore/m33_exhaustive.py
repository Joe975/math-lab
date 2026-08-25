"""Exhaustive treatment of the semiregular Z_33 case for SRG(99,14,1,2).

The unique orbit matrix (orbit_matrix.py -m 33) is

    R = [[2,6,6],[6,2,6],[6,6,2]]

so the graph would be V = Z_33 x {0,1,2} with connection sets S_ij of sizes
|S_ii| = 2 and |S_ij| = 6, S_ji = -S_ij, S_ii = -S_ii, 0 not in S_ii.

Key reduction.  Write A_U(y) = #{(a,b) in U x U : a - b = y} for the
autocorrelation of U, which is symmetric: A_U(-y) = A_U(y), since swapping the
pair negates the difference.  The three *diagonal* strongly-regular equations
c_ii(y) = target read

    A_{S00} + A_{S01} + A_{S02} = target_0        (from orbit 0)
    A_{S01} + A_{S11} + A_{S12} = target_1        (using S_10 = -S_01)
    A_{S02} + A_{S12} + A_{S22} = target_2        (using S_20 = -S_02, S_21 = -S_12)

Moving the (known) diagonal autocorrelations to the right gives three linear
equations in the three unknown functions A_{S01}, A_{S02}, A_{S12}, with an
invertible coefficient matrix.  Hence

    A_{S01} = (f0 + f1 - f2)/2,  A_{S02} = (f0 - f1 + f2)/2,
    A_{S12} = (-f0 + f1 + f2)/2

is FORCED by the choice of the three diagonal sets.  Each diagonal set is
{a, -a} for some a in 1..16, so there are only 16^3 = 4096 cases, and each one
either yields a legal autocorrelation vector for all three off-diagonal sets or
is dead on the spot.

That turns the Z_33 case into a finite, cheap, auditable computation.
"""

import os
import sys
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg

M = 33
K, LAM, MU = 14, 1, 2


def autocorr(U, m=M):
    a = [0] * m
    for x in U:
        for y in U:
            a[(x - y) % m] += 1
    return a


def diag_target(S, m=M):
    """f(y) for a diagonal orbit with diagonal set S: the required value of
    the sum of the two off-diagonal autocorrelations at y."""
    A = autocorr(S, m)
    f = [0] * m
    f[0] = K - A[0]
    for y in range(1, m):
        tgt = LAM if y in S else MU
        f[y] = tgt - A[y]
    return f


def solve_offdiag(f0, f1, f2, m=M):
    """A_{S01}, A_{S02}, A_{S12}, or None if not a legal autocorrelation."""
    out = []
    for combo in ((f0, f1, f2), (f0, f2, f1), (f1, f2, f0)):
        # (x + y - z)/2 pattern, applied as A01,A02,A12 below
        pass
    A01 = [0] * m
    A02 = [0] * m
    A12 = [0] * m
    for y in range(m):
        s01 = f0[y] + f1[y] - f2[y]
        s02 = f0[y] - f1[y] + f2[y]
        s12 = -f0[y] + f1[y] + f2[y]
        for s in (s01, s02, s12):
            if s % 2:
                return None, "half-integer autocorrelation"
            if s < 0:
                return None, "negative autocorrelation"
        A01[y], A02[y], A12[y] = s01 // 2, s02 // 2, s12 // 2
    for name, A in (("A01", A01), ("A02", A02), ("A12", A12)):
        if A[0] != 6:
            return None, f"{name}(0) = {A[0]}, must equal |S| = 6"
        if sum(A) != 36:
            return None, f"sum {name} = {sum(A)}, must equal 36"
        for y in range(m):
            if A[y] != A[(-y) % m]:
                return None, f"{name} is not symmetric"
            if A[y] > 6:
                return None, f"{name}({y}) = {A[y]} exceeds |S| = 6"
    return (A01, A02, A12), "ok"


def subsets_with_autocorr(target, size=6, m=M, limit=None):
    """All size-subsets of Z_m with the given autocorrelation.

    Normalised by translation: every such set can be shifted so that it
    contains 0, and translation does not change the autocorrelation.  The
    unnormalised solutions are exactly the m translates of these.
    """
    found = []
    # build by choosing increasing elements starting from 0
    def rec(chosen, start):
        if limit and len(found) >= limit:
            return
        if len(chosen) == size:
            if autocorr(chosen, m) == target:
                found.append(tuple(chosen))
            return
        for nxt in range(start, m):
            cur = chosen + [nxt]
            # partial prune: no difference may exceed the target count
            a = autocorr(cur, m)
            if any(a[y] > target[y] for y in range(m)):
                continue
            rec(cur, nxt + 1)
    rec([0], 1)
    return found


def main():
    print("=== exhaustive Z_33 analysis for SRG(99,14,1,2) ===")
    print(f"# unique orbit matrix R = [[2,6,6],[6,2,6],[6,6,2]]")
    print(f"# diagonal sets are {{a,-a}}, a in 1..16 -> {16**3} cases\n")

    diag_choices = [frozenset({a, (-a) % M}) for a in range(1, M // 2 + 1)]
    print(f"# distinct diagonal sets: {len(diag_choices)}")

    alive = []
    reasons = {}
    for i0, S00 in enumerate(diag_choices):
        f0 = diag_target(S00)
        for i1, S11 in enumerate(diag_choices):
            f1 = diag_target(S11)
            for i2, S22 in enumerate(diag_choices):
                f2 = diag_target(S22)
                sol, why = solve_offdiag(f0, f1, f2)
                if sol is None:
                    reasons[why] = reasons.get(why, 0) + 1
                else:
                    alive.append((S00, S11, S22, sol))
    total = len(diag_choices) ** 3
    print(f"# cases killed by the forced autocorrelation: {total - len(alive)}"
          f" of {total}")
    for why, c in sorted(reasons.items(), key=lambda kv: -kv[1]):
        print(f"    {c:5d}  {why}")
    print(f"# cases surviving: {len(alive)}\n")

    if not alive:
        print("RESULT: no diagonal choice admits legal off-diagonal "
              "autocorrelations, so the Z_33 orbit matrix does not lift.")
        return

    # for survivors, try to realise the required autocorrelations
    realisable = []
    for idx, (S00, S11, S22, (A01, A02, A12)) in enumerate(alive):
        sets01 = subsets_with_autocorr(A01)
        if not sets01:
            continue
        sets02 = subsets_with_autocorr(A02)
        if not sets02:
            continue
        sets12 = subsets_with_autocorr(A12)
        if not sets12:
            continue
        realisable.append((S00, S11, S22, sets01, sets02, sets12))
        print(f"  survivor {idx}: diag {sorted(S00)},{sorted(S11)},{sorted(S22)} "
              f"-> {len(sets01)}x{len(sets02)}x{len(sets12)} normalised sets")
    print(f"\n# survivors whose autocorrelations are realisable by actual "
          f"subsets: {len(realisable)}")

    if not realisable:
        print("RESULT: every surviving diagonal choice demands an "
              "autocorrelation that no 6-subset of Z_33 has, so the Z_33 "
              "orbit matrix does not lift.")
        return

    # full check over the remaining, now-small, space
    print("\n# checking the remaining space exhaustively")
    import orbit_lift as OL
    checked = 0
    for (S00, S11, S22, s01, s02, s12) in realisable:
        for T01 in s01:
            for sh1 in range(M):
                U01 = frozenset((x + sh1) % M for x in T01)
                for T02 in s02:
                    for sh2 in range(M):
                        U02 = frozenset((x + sh2) % M for x in T02)
                        for T12 in s12:
                            for sh3 in range(M):
                                U12 = frozenset((x + sh3) % M for x in T12)
                                S = [[S00, U01, U02],
                                     [frozenset((-x) % M for x in U01), S11, U12],
                                     [frozenset((-x) % M for x in U02),
                                      frozenset((-x) % M for x in U12), S22]]
                                checked += 1
                                if OL.energy(S, None, M, 3, LAM, MU) == 0:
                                    adj = OL.build_graph(S, M, 3)
                                    ok, got = srg.check_srg(adj)
                                    print(f"  *** WITNESS: {ok} {got} ***")
                                    print(f"      S = {[[sorted(x) for x in row] for row in S]}")
                                    return
    print(f"# combinations checked: {checked}")
    print("RESULT: no lift of the Z_33 orbit matrix exists.")


if __name__ == "__main__":
    main()
