"""Orbit matrices for a group acting with orbits of DIFFERENT sizes.

Generalises orbit_matrix.py, which assumes a semiregular action.  The
conditions are the same two families,

    r_ij * n_i = r_ji * n_j                                            (S)
    sum_l r_il r_lj = k delta_ij + lam r_ij + mu (n_j - delta_ij - r_ij)  (Q)
    sum_j r_ij = k

but R is no longer symmetric, so entries below the diagonal are *forced* by (S)
from entries above it, and the parity rule applies per orbit: for an orbit of
odd size n_i under a cyclic group, r_ii is even.

The case this exists for is an automorphism of order 7 on 99 vertices.  Blind
lemma L3 (see autos.py) forces exactly one fixed point, so the orbits are
1 + 14 x 7, and two further facts are forced before any search starts:

  * the fixed vertex v has k = 14 neighbours and every r_0j is 0 or 7 by (S),
    so N(v) is exactly two full orbits;
  * N(v) induces a perfect matching, so a vertex of one of those two orbits
    has exactly one neighbour inside N(v); with r_ii even this forces
    r_11 = r_22 = 0 and r_12 = r_21 = 1.

Both are imposed as a seed rather than searched.
"""

import argparse
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

from orbit_matrix import rhs, check_orbit_matrix


def enumerate_general(sizes, k, lam, mu, seed=None, limit=None, odd_diag_even=True,
                      report=None):
    """All orbit matrices for the given orbit size profile.

    `seed` is an optional dict {(i,j): value} of forced entries.
    """
    t = len(sizes)
    R = [[None] * t for _ in range(t)]
    forced = dict(seed or {})
    out = []
    stats = {"nodes": 0}

    def set_entry(i, j, v):
        R[i][j] = v
        # (S) forces the transpose entry
        num = v * sizes[i]
        if num % sizes[j]:
            return False
        R[j][i] = num // sizes[j]
        return True

    def clear_entry(i, j):
        R[i][j] = None
        R[j][i] = None

    def eq_ok(upto):
        for i in range(upto + 1):
            for j in range(upto + 1):
                s = 0
                for l in range(t):
                    a = R[i][l]
                    b = R[l][j]
                    if b is None and R[j][l] is not None:
                        num = R[j][l] * sizes[j]
                        if num % sizes[l]:
                            return False
                        b = num // sizes[l]
                    if a is None or b is None:
                        return True
                    s += a * b
                const, coeff = rhs(i, j, sizes, k, lam, mu)
                if s != const + coeff * R[i][j]:
                    return False
        return True

    def row_target_sq(i):
        const, coeff = rhs(i, i, sizes, k, lam, mu)
        return const + coeff * R[i][i]

    def fill(i, j, remaining):
        stats["nodes"] += 1
        if limit and len(out) >= limit:
            return
        if j == t:
            if remaining != 0:
                return
            # diagonal equation (Q) at (i,i) in full
            s = 0
            for l in range(t):
                b = R[l][i]
                if b is None:
                    if R[i][l] is None:
                        return
                    num = R[i][l] * sizes[i]
                    if num % sizes[l]:
                        return
                    b = num // sizes[l]
                s += R[i][l] * b
            if s != row_target_sq(i):
                return
            if not eq_ok(i):
                return
            if report and i >= report:
                print(f"    row {i} ok (nodes={stats['nodes']}, "
                      f"solutions={len(out)})", flush=True)
            if i == t - 1:
                out.append([row[:] for row in R])
                return
            nxt = i + 1
            pre = sum(R[nxt][c] for c in range(nxt) if R[nxt][c] is not None)
            fill(nxt, nxt, k - pre)
            return
        if R[i][j] is not None:
            fill(i, j + 1, remaining - R[i][j])
            return
        if (i, j) in forced:
            v = forced[(i, j)]
            if v > remaining or not set_entry(i, j, v):
                clear_entry(i, j)
                return
            fill(i, j + 1, remaining - v)
            clear_entry(i, j)
            return
        hi = min(sizes[j], remaining)
        for v in range(hi + 1):
            if i == j:
                if sizes[i] == 1 and v != 0:
                    continue
                if odd_diag_even and sizes[i] % 2 == 1 and v % 2:
                    continue
            if not set_entry(i, j, v):
                clear_entry(i, j)
                continue
            fill(i, j + 1, remaining - v)
            clear_entry(i, j)

    # row 0
    pre = 0
    fill(0, 0, k - pre)
    return out, stats


def z7_profile():
    """Orbit sizes and the forced seed for an order-7 automorphism of a
    (99,14,1,2) graph: one fixed point and fourteen 7-orbits."""
    sizes = [1] + [7] * 14
    seed = {(0, 0): 0}
    # v is adjacent to exactly two whole orbits; label them 1 and 2
    for j in range(1, 15):
        seed[(0, j)] = 7 if j in (1, 2) else 0
    # N(v) is a perfect matching, so r_11 = r_22 = 0 and r_12 = 1
    seed[(1, 1)] = 0
    seed[(1, 2)] = 1
    seed[(2, 2)] = 0
    return sizes, seed


def z33_profiles():
    """Every possible orbit-size profile of an order-33 automorphism.

    Let sigma have order 33.  A fixed point, or an orbit of size 3, would give
    a point whose stabiliser contains the order-11 subgroup <sigma^3>; but an
    order-11 automorphism with a fixed point is the identity (blind lemma L2),
    which contradicts |sigma| = 33.  So every orbit has size 11 or 33, and

        99 = 33a + 11b,   a >= 1

    (a = 0 would make every stabiliser contain <sigma^11>, forcing sigma^11 to
    fix every point, i.e. |sigma| = 11).  That leaves exactly three profiles.
    """
    out = {}
    for a in range(1, 4):
        rem = 99 - 33 * a
        if rem % 11:
            continue
        b = rem // 11
        out[f"z33_{a}x33_{b}x11"] = [33] * a + [11] * b
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default="z7")
    ap.add_argument("--sizes", type=str, default=None,
                    help="comma-separated orbit sizes, overrides --profile")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--report", type=int, default=None)
    args = ap.parse_args()

    k, lam, mu = 14, 1, 2
    if args.sizes:
        sizes = [int(x) for x in args.sizes.split(",")]
        seed = {}
    elif args.profile in z33_profiles():
        sizes, seed = z33_profiles()[args.profile], {}
    else:
        sizes, seed = z7_profile()
    print(f"# orbit sizes: {sizes}")
    print(f"# forced seed: {len(seed)} entries "
          f"(fixed vertex sees two whole orbits; N(v) is a perfect matching)")
    sols, stats = enumerate_general(sizes, k, lam, mu, seed=seed,
                                    limit=args.limit, report=args.report)
    print(f"# nodes: {stats['nodes']}")
    print(f"# orbit matrices: {len(sols)}")
    good = 0
    for R in sols:
        ok, why = check_orbit_matrix(R, sizes, k, lam, mu)
        if ok:
            good += 1
        else:
            print(f"  VERIFICATION FAILURE: {why}")
    print(f"# independently verified: {good}/{len(sols)}")
    for R in sols[:5]:
        print("  diag:", [R[i][i] for i in range(len(sizes))])
    if args.out:
        json.dump({"sizes": sizes, "k": k, "lam": lam, "mu": mu,
                   "count": len(sols), "matrices": sols},
                  open(args.out, "w"), indent=1)
        print(f"# wrote {args.out}")


if __name__ == "__main__":
    main()
