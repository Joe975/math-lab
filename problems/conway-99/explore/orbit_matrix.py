"""Orbit matrices for SRG(n, k, lam, mu) under an assumed automorphism group.

If a group H <= Aut(Gamma) has orbits O_1..O_t with |O_i| = n_i, then
r_ij := |N(u) & O_j| is independent of the choice of u in O_i, and the matrix
R = (r_ij) satisfies two exact families of integer conditions.

Counting edges between O_i and O_j:

    r_ij * n_i = r_ji * n_j                                            (S)

Projecting A^2 = k I + lam A + mu (J - I - A) onto orbits.  For u in O_i,

    sum_{v in O_j} (A^2)_{uv} = sum_w A_uw |N(w) & O_j| = sum_l r_il r_lj

while the right-hand side contributes k delta_ij + lam r_ij
+ mu (n_j - delta_ij - r_ij), so

    sum_l r_il r_lj = k delta_ij + lam r_ij + mu (n_j - delta_ij - r_ij)   (Q)

Together with 0 <= r_ij <= n_j and sum_j r_ij = k, these pin R hard enough to
enumerate exhaustively for small t.  Everything here is integer arithmetic.

Parity note: for a cyclic H acting semiregularly with |H| = m odd, the subgraph
induced on an orbit is a circulant on Z_m whose connection set is closed under
negation and misses 0, so r_ii is even.  For m even the element m/2 is its own
negative, so r_ii may be odd.

This module ENUMERATES orbit matrices.  An orbit matrix is a necessary
condition, not a graph: lifting one to an actual graph is a separate search
(see orbit_lift.py), and many orbit matrices do not lift.
"""

import argparse
import itertools
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))


# ---------------------------------------------------------------------------
# the conditions, written out plainly so they can be checked directly
# ---------------------------------------------------------------------------

def rhs(i, j, sizes, k, lam, mu):
    """k delta_ij + lam r_ij + mu (n_j - delta_ij - r_ij), minus the r_ij part.

    Returned as (const, coeff) with the equation reading
        sum_l r_il r_lj = const + coeff * r_ij
    """
    d = 1 if i == j else 0
    const = k * d + mu * (sizes[j] - d)
    coeff = lam - mu
    return const, coeff


def check_orbit_matrix(R, sizes, k, lam, mu):
    """Independent verification of a candidate R, by direct evaluation.

    Deliberately written as plain nested loops rather than reusing anything the
    search does, so it is a second implementation of the conditions.
    """
    t = len(sizes)
    if len(R) != t or any(len(row) != t for row in R):
        return False, "shape"
    for i in range(t):
        if sum(R[i]) != k:
            return False, f"row {i} sums to {sum(R[i])}, not k={k}"
        for j in range(t):
            if not (0 <= R[i][j] <= sizes[j]):
                return False, f"r[{i}][{j}]={R[i][j]} out of range"
            if R[i][j] * sizes[i] != R[j][i] * sizes[j]:
                return False, f"(S) fails at {i},{j}"
    for i in range(t):
        for j in range(t):
            lhs = 0
            for l in range(t):
                lhs += R[i][l] * R[l][j]
            const, coeff = rhs(i, j, sizes, k, lam, mu)
            if lhs != const + coeff * R[i][j]:
                return False, (f"(Q) fails at {i},{j}: lhs={lhs} "
                               f"rhs={const + coeff * R[i][j]}")
    return True, "ok"


# ---------------------------------------------------------------------------
# enumeration for the semiregular case (all orbits the same size)
# ---------------------------------------------------------------------------

def enumerate_semiregular(m, k, lam, mu, n, diag_even=None, limit=None,
                          progress=None):
    """All orbit matrices for a semiregular cyclic group of order m.

    n/m orbits, all of size m, so R is symmetric.  Rows are built one at a
    time over the upper triangle; after row i is complete every equation (Q)
    with both indices <= i can be evaluated in full, because for l > i the
    entry r_lj equals r_jl and row j is already finished.  That is what makes
    the pruning strong.
    """
    assert n % m == 0
    t = n // m
    if diag_even is None:
        diag_even = (m % 2 == 1)

    R = [[None] * t for _ in range(t)]
    out = []
    stats = {"nodes": 0, "rows_completed": 0}

    # per-row necessary condition from the diagonal equation (Q) at (i,i):
    #   sum_l r_il^2 = const + coeff * r_ii
    def row_ok(i):
        s = sum(R[i])
        if s != k:
            return False
        sq = sum(x * x for x in R[i])
        const, coeff = rhs(i, i, [m] * t, k, lam, mu)
        return sq == const + coeff * R[i][i]

    def equations_ok(upto):
        sizes = [m] * t
        for i in range(upto + 1):
            for j in range(upto + 1):
                lhs = 0
                for l in range(t):
                    a, b = R[i][l], R[l][j]
                    if a is None or b is None:
                        return True   # cannot evaluate yet
                    lhs += a * b
                const, coeff = rhs(i, j, sizes, k, lam, mu)
                if lhs != const + coeff * R[i][j]:
                    return False
        return True

    const_ii, coeff_ii = rhs(0, 0, [m] * t, k, lam, mu)

    def fill_row(i, j, remaining, sq):
        """sq is the running sum of squares of row i, over all columns set so
        far (including those forced by symmetry).  Once r_ii is known the
        diagonal equation fixes the final value of sq exactly, which prunes
        hard: sum_l r_il^2 = const_ii + coeff_ii * r_ii."""
        stats["nodes"] += 1
        if limit and len(out) >= limit:
            return
        if R[i][i] is not None:
            target_sq = const_ii + coeff_ii * R[i][i]
            if sq > target_sq:
                return
            # the remaining `remaining` units spread over (t-j) cells add at
            # least ceil(remaining^2/(t-j)) and at most remaining^2 to sq
            cells = t - j
            if cells > 0:
                lo_add = (remaining * remaining + cells - 1) // cells
                hi_add = remaining * remaining
            else:
                lo_add = hi_add = 0
            if sq + lo_add > target_sq or sq + hi_add < target_sq:
                return
        if j == t:
            if remaining != 0:
                return
            if not row_ok(i):
                return
            if not equations_ok(i):
                return
            stats["rows_completed"] += 1
            if progress and i >= progress:
                print(f"    row {i} completed (nodes={stats['nodes']}, "
                      f"solutions so far={len(out)})", flush=True)
            if i == t - 1:
                out.append([row[:] for row in R])
                return
            nx = i + 1
            pre = R[nx][:nx]
            fill_row(nx, nx, k - sum(pre), sum(x * x for x in pre))
            return
        if R[i][j] is not None:          # forced by symmetry from an earlier row
            fill_row(i, j + 1, remaining - R[i][j], sq + R[i][j] ** 2)
            return
        hi = min(m, remaining)
        lo = 0
        for v in range(lo, hi + 1):
            if i == j and diag_even and v % 2:
                continue
            R[i][j] = v
            if i != j:
                R[j][i] = v
            fill_row(i, j + 1, remaining - v, sq + v * v)
            R[i][j] = None
            if i != j:
                R[j][i] = None

    fill_row(0, 0, k, 0)
    return out, stats


# ---------------------------------------------------------------------------
# spectral pre-filter: which traces are possible at all
# ---------------------------------------------------------------------------

def admissible_traces(m, k, lam, mu, n):
    """Possible values of tr(R), from the spectrum of R.

    On the all-ones vector R has eigenvalue k.  On its orthogonal complement
    (Q) reduces to theta^2 - (lam-mu) theta - (k-mu) = 0, i.e. theta is a
    restricted eigenvalue r or s of the parameter set.  So
    tr(R) = k + a*r + b*s with a + b = t - 1.
    """
    import srg
    t = n // m
    spec = srg.spectrum(n, k, lam, mu)
    if not spec["ok"] or spec.get("conference"):
        return None
    r, s = spec["r"], spec["s"]
    out = []
    for a in range(t):
        b = t - 1 - a
        tr = k + a * r + b * s
        if tr < 0:
            continue
        if m % 2 == 1 and tr % 2:
            continue                      # every r_ii is even
        if tr > t * (m if m % 2 == 0 else m - 1):
            continue
        out.append({"a_mult_of_r": a, "b_mult_of_s": b, "trace": tr})
    return {"r": r, "s": s, "t": t, "options": out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", type=int, required=True,
                    help="order of the assumed semiregular cyclic group")
    ap.add_argument("-n", type=int, default=99)
    ap.add_argument("-k", type=int, default=14)
    ap.add_argument("--lam", type=int, default=1)
    ap.add_argument("--mu", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--out", type=str, default=None)
    ap.add_argument("--progress", type=int, default=None)
    args = ap.parse_args()

    n, k, lam, mu, m = args.n, args.k, args.lam, args.mu, args.m
    if n % m:
        print(f"{m} does not divide {n}: no semiregular action")
        return
    t = n // m
    print(f"# assumed group: cyclic of order {m}, acting semiregularly")
    print(f"# {t} orbits of size {m} on {n} vertices; k={k} lam={lam} mu={mu}")
    tr = admissible_traces(m, k, lam, mu, n)
    if tr:
        print(f"# restricted eigenvalues r={tr['r']} s={tr['s']}; "
              f"admissible traces: {[o['trace'] for o in tr['options']]}")

    sols, stats = enumerate_semiregular(m, k, lam, mu, n, limit=args.limit,
                                        progress=args.progress)
    print(f"# search nodes: {stats['nodes']}")
    print(f"# orbit matrices found: {len(sols)}")

    verified = 0
    for R in sols:
        ok, why = check_orbit_matrix(R, [m] * t, k, lam, mu)
        if not ok:
            print(f"  VERIFICATION FAILURE: {why}")
        else:
            verified += 1
    print(f"# independently verified: {verified}/{len(sols)}")

    for idx, R in enumerate(sols[:10]):
        print(f"  solution {idx}: diag={[R[i][i] for i in range(t)]} "
              f"trace={sum(R[i][i] for i in range(t))}")
        for row in R:
            print("     ", row)
    if len(sols) > 10:
        print(f"  ... and {len(sols)-10} more")

    if args.out:
        with open(args.out, "w") as fh:
            json.dump({"m": m, "n": n, "k": k, "lam": lam, "mu": mu,
                       "t": t, "count": len(sols), "nodes": stats["nodes"],
                       "verified": verified,
                       "admissible_traces": tr,
                       "matrices": sols}, fh, indent=1)
        print(f"# wrote {args.out}")


if __name__ == "__main__":
    main()
