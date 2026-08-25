"""Orbit matrices with row-inner-product pruning, and the order-7 profile.

The conditions are the same as in orbit_matrix.py / orbit_general.py:

    (S)  r_ij * n_i = r_ji * n_j
    (R)  sum_j r_ij = k
    (Q)  sum_l r_il r_lj = k [i=j] + lam r_ij + mu (n_j - [i=j] - r_ij)

What is different is the search. Rows are filled in order, and (S) makes the
column below the diagonal a function of the rows above it, so once rows a and
i are both complete the (Q) equation at (a, i) is *fully determined*. While
row i is being filled, therefore, the (Q) equations with every earlier row are
linear constraints on row i's free entries with known coefficients, and the
diagonal (Q) equation at (i, i) is a fixed weighted sum of squares. Each is
bounded after every assignment. orbit_general.py only checked (Q) once a row
was complete, which is why it never reached a first solution at 15 orbits.

Two further prunes, both optional and both switched on for the Z_7 profile:

* Symmetry breaking. Orbits with identical seed rows and equal size are
  interchangeable, so among them the key (seeded column entries, diagonal)
  is required to be non-increasing. Every equivalence class has a
  representative in that order, so no class is lost; a class can still
  appear more than once, which is why solutions are afterwards reduced by
  an exact canonical form (`canonical`).

* Trace constraint. Decomposing the adjacency matrix over the characters of
  Z_7: the trivial character carries the 15x15 quotient R, and each of the 6
  non-trivial characters carries a 14x14 block whose eigenvalues are 3 and
  -4 with the same multiplicity r for all six (Galois conjugation preserves
  rank). With a = multiplicity of 3 in R this gives 54 = a + 6r, so a is a
  multiple of 6, and tr R = 14 + 3a - 4(14 - a) = 7a - 42 is then 0 or 42
  (84 exceeds the maximum 12 * 6). Passed as `--trace 0,42`; the (Q)
  conditions alone only force tr R to be a multiple of 7.

Calibration lives in tests/test_conway99_orbit7.py: the enumeration must
contain the real orbit matrix of Paley(9)/Z_3, Petersen/Z_3, Paley(9)/Z_2 and
rook(3)/Z_3, and reproduce 002's Z_33 counts.
"""

import argparse
import itertools
import json
import os
import sys
import time
from fractions import Fraction
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)


def lcm(a, b):
    return a * b // gcd(a, b)


def rhs(i, j, sizes, k, lam, mu):
    d = 1 if i == j else 0
    return k * d + mu * (sizes[j] - d), lam - mu


def check_orbit_matrix(R, sizes, k, lam, mu):
    """Direct check of (S), (R), (Q) by plain loops."""
    t = len(sizes)
    for i in range(t):
        if sum(R[i]) != k:
            return False, f"row {i} sums to {sum(R[i])}"
        for j in range(t):
            if not 0 <= R[i][j] <= sizes[j]:
                return False, f"r[{i}][{j}] out of range"
            if R[i][j] * sizes[i] != R[j][i] * sizes[j]:
                return False, f"(S) fails at {i},{j}"
    for i in range(t):
        for j in range(t):
            lhs = sum(R[i][l] * R[l][j] for l in range(t))
            c, e = rhs(i, j, sizes, k, lam, mu)
            if lhs != c + e * R[i][j]:
                return False, f"(Q) fails at {i},{j}"
    return True, "ok"


def enumerate_orbit_matrices(sizes, k, lam, mu, seed=None, odd_diag_even=True,
                             limit=None, node_budget=None, symmetry=True,
                             trace_allowed=None, max_entry=None, spectral=None,
                             symmetry_rows=None, dump_rows=None, progress=None):
    """All orbit matrices for `sizes`, up to the optional symmetry breaking.

    Returns (matrices, stats). stats["exhausted"] is False if `limit` or
    `node_budget` cut the search short.
    """
    t = len(sizes)
    seed = dict(seed or {})
    # (S) makes seeds symmetric-ish: a seeded r_ij fixes r_ji.
    for (i, j), v in list(seed.items()):
        num = v * sizes[i]
        if num % sizes[j]:
            return [], {"nodes": 0, "exhausted": True, "reason": "seed breaks (S)"}
        seed.setdefault((j, i), num // sizes[j])
    R = [[None] * t for _ in range(t)]
    out = []
    stats = {"nodes": 0, "exhausted": True, "rows_reached": 0}
    L = 1
    for n in sizes:
        L = lcm(L, n)

    # Interchangeable classes for symmetry breaking: same size, same seed row.
    def seed_signature(i):
        # the diagonal seed is a property of the orbit, not of its label
        return (sizes[i], tuple(sorted((-1 if j == i else j, v)
                                       for (a, j), v in seed.items() if a == i)))

    cls = {}
    for i in range(t):
        cls.setdefault(seed_signature(i), []).append(i)
    prev_same = [None] * t  # previous orbit in the same class, if any
    for members in cls.values():
        for a, b in zip(members, members[1:]):
            prev_same[b] = a
    distinguished = sorted(m[0] for m in cls.values() if len(m) == 1)

    swaps = [(a, b) for members in cls.values() for a, b in zip(members, members[1:])]

    def lex_ok(upto):
        """Lex-leader (maximum) under each adjacent same-class transposition.

        Relabelling orbits a <-> b maps R to P R P^T, another orbit matrix
        of the same class. Requiring R >= P R P^T in row-major order for
        every such adjacent pair keeps the lex-maximal representative of
        each class (and possibly more; `canonical` dedupes afterwards).
        Only positions decided so far are compared; the comparison stops
        undecided, without pruning, at the first unknown position."""
        def get(p, q):
            if R[p][q] is not None:
                return R[p][q]
            if R[q][p] is not None:
                return R[q][p] * sizes[q] // sizes[p]
            return None

        for a, b in swaps:
            def img(x):
                return b if x == a else a if x == b else x
            decided = True
            for p in range(t if symmetry_rows is None else min(t, symmetry_rows)):
                for q in range(t):
                    x, y = get(p, q), get(img(p), img(q))
                    if x is None or y is None:
                        decided = False
                        break
                    if x > y:
                        break
                    if x < y:
                        return False
                else:
                    continue
                break
            if not decided:
                continue
        return True

    stats["spectral_cuts"] = 0
    stats["spectral_checks"] = 0
    stats["max_rank_r"] = stats["max_rank_s"] = 0
    # trace pruning during the search: the diagonal is bounded entrywise
    diag_max = [0 if n == 1 else (n - 1 if (odd_diag_even and n % 2 == 1) else n)
                for n in sizes]
    for (i, j), v in seed.items():
        if i == j:
            diag_max[i] = v
    diag_sufmax = [0] * (t + 1)
    for i in range(t - 1, -1, -1):
        diag_sufmax[i] = diag_sufmax[i + 1] + diag_max[i]
    n_total = sum(sizes)
    scale = 1
    for x in (n_total, *sizes):
        scale = lcm(scale, x)
    scale *= 4  # covers the (r - s) denominator for any integral spectrum

    def spectral_ok(upto):
        """Leading blocks of both eigenvalue projectors must be PSD, and
        their ranks cannot exceed the multiplicities a (for r) and
        t - 1 - a (for s) that `spectral` asserts."""
        a = spectral
        stats["spectral_checks"] += 1
        Zr, Zs, _, _ = spectral_blocks(R, sizes, upto, k, lam, mu, scale)
        ok, rk = psd_rank(Zr)
        stats["max_rank_r"] = max(stats["max_rank_r"], rk)
        if not ok or rk > a:
            return False
        ok, rk = psd_rank(Zs)
        stats["max_rank_s"] = max(stats["max_rank_s"], rk)
        if not ok or rk > t - 1 - a:
            return False
        return True

    def fill_row(i):
        """Fill row i's free entries (j >= i), rows < i complete."""
        if node_budget is not None and stats["nodes"] >= node_budget:
            stats["exhausted"] = False
            return
        if limit is not None and len(out) >= limit:
            stats["exhausted"] = False
            return
        if i == t:
            if trace_allowed is not None and sum(R[a][a] for a in range(t)) not in trace_allowed:
                return
            M = [row[:] for row in R]
            ok, why = check_orbit_matrix(M, sizes, k, lam, mu)
            assert ok, why
            out.append(M)
            return
        stats["rows_reached"] = max(stats["rows_reached"], i)
        if dump_rows is not None and i == dump_rows:
            stats.setdefault("dumped", []).append([row[:] for row in R[:i]])
        stats.setdefault("row_completions", {})
        stats["row_completions"][i - 1] = stats["row_completions"].get(i - 1, 0) + 1
        if progress and i <= progress[0]:
            progress[1](i, stats["nodes"], len(out))
        # entries j < i are forced by (S) from rows above
        pre_sum = 0
        for j in range(i):
            R[i][j] = R[j][i] * sizes[j] // sizes[i]
            pre_sum += R[i][j]
        free = list(range(i, t))
        # linear constraints: for each complete row a < i,
        #   sum_l r_il * c_la = const_ia + coeff * r_ia,   c_la = r_al n_a / n_l
        # which after moving known terms is  sum_{j>=i} x_j * w_aj = target_a
        lin = []
        for a in range(i):
            c, e = rhs(i, a, sizes, k, lam, mu)
            target = c + e * R[i][a]
            for l in range(i):
                target -= R[i][l] * (R[a][l] * sizes[a] // sizes[l])
            w = [R[a][j] * sizes[a] // sizes[j] for j in free]
            lin.append((w, target))
        # diagonal (Q) scaled by L: sum_l r_il^2 (n_i L / n_l) = L(c + e r_ii)
        c_ii, e_ii = rhs(i, i, sizes, k, lam, mu)
        sq_known = sum(R[i][l] ** 2 * (sizes[i] * L // sizes[l]) for l in range(i))
        sq_w = [sizes[i] * L // sizes[j] for j in free]
        sum_target = k - pre_sum
        # per-entry upper bounds
        ub = []
        for j in free:
            if (i, j) in seed:
                ub.append(seed[(i, j)])
            elif max_entry is not None and j != i and sizes[j] == sizes[i]:
                ub.append(min(sizes[j], max_entry))
            else:
                ub.append(sizes[j])
        # suffix maxima for bound pruning
        nfree = len(free)
        lin_sufmax = []
        for w, _ in lin:
            s = [0] * (nfree + 1)
            for p in range(nfree - 1, -1, -1):
                s[p] = s[p + 1] + w[p] * ub[p]
            lin_sufmax.append(s)
        sq_sufmax = [0] * (nfree + 1)
        for p in range(nfree - 1, -1, -1):
            sq_sufmax[p] = sq_sufmax[p + 1] + sq_w[p] * ub[p] ** 2
        sum_sufmax = [0] * (nfree + 1)
        for p in range(nfree - 1, -1, -1):
            sum_sufmax[p] = sum_sufmax[p + 1] + ub[p]

        lin_acc = [0] * len(lin)

        def rec(p, rem_sum, sq_acc):
            stats["nodes"] += 1
            if node_budget is not None and stats["nodes"] >= node_budget:
                stats["exhausted"] = False
                return
            if limit is not None and len(out) >= limit:
                return
            if p == nfree:
                if rem_sum != 0:
                    return
                if sq_acc != L * (c_ii + e_ii * R[i][i]):
                    return
                for q, (w, target) in enumerate(lin):
                    if lin_acc[q] != target:
                        return
                if symmetry and not lex_ok(i):
                    return
                if trace_allowed is not None:
                    tr = sum(R[a][a] for a in range(i + 1))
                    if not any(tr <= x <= tr + diag_sufmax[i + 1] for x in trace_allowed):
                        return
                if spectral is not None and not spectral_ok(i):
                    stats["spectral_cuts"] += 1
                    return
                fill_row(i + 1)
                return
            j = free[p]
            hi = min(ub[p], rem_sum)
            if j == i and R[i][i] is not None:
                hi = min(hi, R[i][i])
            for v in range(hi + 1):
                if j == i:
                    if sizes[i] == 1 and v:
                        break
                    if odd_diag_even and sizes[i] % 2 == 1 and v % 2:
                        continue
                if (v * sizes[i]) % sizes[j]:
                    continue  # (S) needs an integer transpose entry
                if (i, j) in seed and seed[(i, j)] != v:
                    continue
                if rem_sum - v > sum_sufmax[p + 1]:
                    continue
                nsq = sq_acc + sq_w[p] * v * v
                # the sum-of-squares target is known once r_ii is set, which
                # happens at p == 0 since the diagonal is the first free entry
                ok = True
                if R[i][i] is not None or j == i:
                    dii = v if j == i else R[i][i]
                    tgt = L * (c_ii + e_ii * dii)
                    if nsq > tgt or nsq + sq_sufmax[p + 1] < tgt:
                        ok = False
                if not ok:
                    continue
                for q, (w, target) in enumerate(lin):
                    acc = lin_acc[q] + w[p] * v
                    if acc > target or acc + lin_sufmax[q][p + 1] < target:
                        ok = False
                        break
                if not ok:
                    continue
                R[i][j] = v
                for q, (w, _) in enumerate(lin):
                    lin_acc[q] += w[p] * v
                rec(p + 1, rem_sum - v, nsq)
                for q, (w, _) in enumerate(lin):
                    lin_acc[q] -= w[p] * v
                R[i][j] = None

        rec(0, sum_target, sq_known)
        for j in free:
            R[i][j] = None
        for j in range(i):
            R[i][j] = None

    fill_row(0)
    return out, stats


def canonical(R, sizes, fixed=()):
    """Exact canonical form under simultaneous relabelling of orbits of equal
    size, refined by (size, row-multiset) and then brute-forced within cells.
    `fixed` orbits are never moved. Feasible when the refinement cells are
    small, which they are for everything enumerated here."""
    t = len(sizes)
    inv = [(sizes[i], R[i][i], tuple(sorted(R[i])), tuple(sorted(R[j][i] for j in range(t))))
           for i in range(t)]
    # iterate refinement with neighbour multisets
    colour = [inv[i] for i in range(t)]
    for _ in range(t):
        new = [(colour[i], tuple(sorted((colour[j], R[i][j], R[j][i]) for j in range(t))))
               for i in range(t)]
        ranks = {c: n for n, c in enumerate(sorted(set(new)))}
        nc = [ranks[c] for c in new]
        if len(set(nc)) == len(set(colour)):
            colour = nc
            break
        colour = nc
    cells = {}
    for i in range(t):
        cells.setdefault(colour[i], []).append(i)
    fixed = set(fixed)
    cell_list = [sorted(v) for _, v in sorted(cells.items())]
    best = None
    for perms in itertools.product(*[
            [tuple(c)] if any(i in fixed for i in c) else list(itertools.permutations(c))
            for c in cell_list]):
        order = [i for p in perms for i in p]
        M = tuple(tuple(R[a][b] for b in order) for a in order)
        if best is None or M < best:
            best = M
    return best


def psd_rank(M):
    """(is_psd, rank) of a symmetric matrix of integers/Fractions, by exact
    LDL^T without pivoting. For a PSD matrix a zero pivot forces its whole
    row and column to be zero at that stage, which is checked; a negative
    pivot means not PSD."""
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    rank = 0
    for p in range(n):
        d = A[p][p]
        if d < 0:
            return False, rank
        if d == 0:
            if any(A[p][q] != 0 for q in range(p + 1, n)):
                return False, rank
            continue
        rank += 1
        for q in range(p + 1, n):
            f = A[q][p] / d
            if f == 0:
                continue
            for r in range(q, n):
                A[q][r] -= f * A[p][r]
            for r in range(q, n):
                A[r][q] = A[q][r]
    return True, rank


def spectral_blocks(R, sizes, upto, k, lam, mu, scale):
    """Integer matrices congruent to the eigenvalue projectors of the
    quotient, on orbits 0..upto. With eigenvalues k, r, s of the graph,
    B = D^{1/2} R D^{-1/2} = k ww^T + r P_r + s P_s and I = ww^T + P_r + P_s,
    so   P_r = (B - sI - (k - s) ww^T) / (r - s)
         P_s = ((r I - B) + (k - r) ww^T) / (r - s)
    and D^{-1/2} P D^{-1/2} (same signature and rank) has entries
         P_r: (r_ij/n_j - s d_ij/n_j - (k - s)/n) / (r - s)
         P_s: (r d_ij/n_j - r_ij/n_j + (k - r)/n) / (r - s)
    Multiplied by `scale` (a common denominator) these are integers."""
    n = sum(sizes)
    disc = (lam - mu) ** 2 + 4 * (k - mu)
    sq = int(round(disc ** 0.5))
    assert sq * sq == disc, "non-integral eigenvalues"
    r = Fraction(lam - mu + sq, 2)
    s_ = Fraction(lam - mu - sq, 2)
    Zr, Zs = [], []
    for i in range(upto + 1):
        rowr, rows = [], []
        for j in range(upto + 1):
            d = 1 if i == j else 0
            nj = sizes[j]
            vr = (Fraction(R[i][j], nj) - s_ * d / nj - Fraction(k - s_, n)) / (r - s_)
            vs = (r * d / nj - Fraction(R[i][j], nj) + Fraction(k - r, n)) / (r - s_)
            rowr.append(vr * scale)
            rows.append(vs * scale)
        Zr.append(rowr)
        Zs.append(rows)
    return Zr, Zs, r, s_


def z7_profile():
    sizes = [1] + [7] * 14
    seed = {(0, 0): 0}
    for j in range(1, 15):
        seed[(0, j)] = 7 if j in (1, 2) else 0
    seed[(1, 1)] = 0
    seed[(1, 2)] = 1
    seed[(2, 2)] = 0
    return sizes, seed


def parse_seed(text):
    seed = {}
    if not text:
        return seed
    for item in text.split(";"):
        ij, v = item.split("=")
        i, j = ij.split(",")
        seed[(int(i), int(j))] = int(v)
    return seed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--profile", default=None, help="z7 for the order-7 profile")
    ap.add_argument("--sizes", default=None)
    ap.add_argument("--seed", default="")
    ap.add_argument("--k", type=int, default=14)
    ap.add_argument("--lam", type=int, default=1)
    ap.add_argument("--mu", type=int, default=2)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--node-budget", type=int, default=None)
    ap.add_argument("--no-symmetry", action="store_true")
    ap.add_argument("--dump-rows", type=int, default=None,
                    help="record the completed rows 0..N-1 every time row N is "
                         "reached (the seeds for an independent re-run)")
    ap.add_argument("--symmetry-rows", type=int, default=None,
                    help="only break symmetry using the first N rows (a weaker, "
                         "trivially valid mode used as an internal cross-check)")
    ap.add_argument("--trace", default=None, help="comma list of allowed traces")
    ap.add_argument("--diag-zero", action="store_true",
                    help="seed every diagonal entry to 0 (see the docstring)")
    ap.add_argument("--max-entry", type=int, default=None,
                    help="cap on off-diagonal entries (see the docstring)")
    ap.add_argument("--spectral", type=int, default=None,
                    help="multiplicity a of the larger eigenvalue in the quotient; "
                         "enables the projector PSD/rank prune")
    ap.add_argument("--out", default=None)
    ap.add_argument("--progress-rows", type=int, default=-1)
    args = ap.parse_args()
    if args.profile == "z7":
        sizes, seed = z7_profile()
    else:
        sizes = [int(x) for x in args.sizes.split(",")]
        seed = parse_seed(args.seed)
    trace = None if args.trace is None else {int(x) for x in args.trace.split(",")}
    if args.diag_zero:
        for i in range(len(sizes)):
            seed[(i, i)] = 0
    t0 = time.time()

    def report(i, nodes, sols):
        print(f"  row {i} nodes={nodes} solutions={sols} {time.time() - t0:.0f}s", flush=True)

    mats, stats = enumerate_orbit_matrices(
        sizes, args.k, args.lam, args.mu, seed=seed, limit=args.limit,
        node_budget=args.node_budget, symmetry=not args.no_symmetry,
        trace_allowed=trace, max_entry=args.max_entry, spectral=args.spectral,
        symmetry_rows=args.symmetry_rows, dump_rows=args.dump_rows,
        progress=(args.progress_rows, report))
    classes = {}
    for M in mats:
        classes.setdefault(canonical(M, sizes, fixed=[i for (i, j) in seed]), M)
    print(json.dumps({"sizes": sizes, "solutions": len(mats), "classes": len(classes),
                      "nodes": stats["nodes"], "exhausted": stats["exhausted"],
                      "spectral_cuts": stats.get("spectral_cuts"),
                      "spectral_checks": stats.get("spectral_checks"),
                      "max_rank": [stats.get("max_rank_r"), stats.get("max_rank_s")],
                      "row_completions": stats.get("row_completions"),
                      "rows_reached": stats["rows_reached"],
                      "seconds": round(time.time() - t0, 1)}))
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump({"sizes": sizes, "seed": {f"{i},{j}": v for (i, j), v in seed.items()},
                       "dumped": stats.get("dumped"),
                       "k": args.k, "lam": args.lam, "mu": args.mu, "stats": stats,
                       "trace_allowed": sorted(trace) if trace else None,
                       "diag_zero": args.diag_zero, "max_entry": args.max_entry,
                       "spectral": args.spectral,
                       "symmetry": not args.no_symmetry,
                       "classes": [list(map(list, c)) for c in classes]}, f)
    else:
        for c in classes:
            print(json.dumps([list(r) for r in c]))


if __name__ == "__main__":
    main()
