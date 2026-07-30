#!/usr/bin/env python3
"""Identify census masks that are Hanner polytopes.

Two tests, in order of strength:
  1. B4 test: the mask's canonical form (minimum over the 384 signed
     coordinate permutations -- an INDEPENDENT Python reimplementation of
     census.c's canon) equals that of one of the harness-generated Hanner
     vertex sets in R^4.
  2. GL test: an explicit invertible rational matrix T with T(V) = W as
     point sets, for some harness Hanner vertex set W.  Found by fixing a
     linear basis v1..v4 inside V and trying all ordered 4-tuples of W.
     A found T is a certificate; failure for all W proves the body is not
     linearly equivalent to any Hanner polytope (the search is exhaustive
     over images of a basis, and T is determined by where a basis goes).

Usage:
  hanner_match.py list                 the harness Hanner masks + B4 orbits
  hanner_match.py match HEX [HEX...]   classify census masks
  hanner_match.py matchfile LIST.txt   classify first field of each line
"""
import itertools
import os
import sys
from fractions import Fraction

HARNESS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "..", "..", "harness", "mahler-4d")
sys.path.insert(0, HARNESS)
import polytope
from polytope import Q


def pair_reps():
    out = []
    for p in itertools.product((-1, 0, 1), repeat=4):
        if p == (0, 0, 0, 0):
            continue
        if next(x for x in p if x) == 1:
            out.append(p)
    return out


REPS = pair_reps()
IDX = {p: i for i, p in enumerate(REPS)}


def normalize(p):
    fz = next((x for x in p if x), 0)
    return tuple(-x for x in p) if fz < 0 else p


def pair_perms():
    perms = []
    for perm in itertools.permutations(range(4)):
        for signs in itertools.product((1, -1), repeat=4):
            perms.append(tuple(
                IDX[normalize(tuple(signs[c] * p[perm[c]] for c in range(4)))]
                for p in REPS))
    return perms


PERMS = pair_perms()


def canon(mask):
    bits = [i for i in range(40) if mask >> i & 1]
    best = mask
    for perm in PERMS:
        m = 0
        for i in bits:
            m |= 1 << perm[i]
        best = min(best, m)
    return best


def points_to_mask(points):
    mask = 0
    for p in points:
        ip = tuple(int(x) for x in p)
        assert all(abs(x) <= 1 for x in ip), "point outside {0,+-1}^4"
        mask |= 1 << IDX[normalize(ip)]
    return mask


def mask_points(mask):
    pts = []
    for i in range(40):
        if mask >> i & 1:
            p = tuple(Q(x) for x in REPS[i])
            pts.append(p)
            pts.append(tuple(-x for x in p))
    return pts


def hanner_masks():
    out = []
    for poly in polytope.hanner_polytopes(4):
        out.append((points_to_mask(list(poly)), len(poly)))
    return out


def solve_T(basis, images):
    """T with T basis_i = images_i, as a 4x4 Fraction matrix, or None."""
    # T M = N  =>  columns: T = N M^-1; solve via polytope.solve_square per row
    # M has basis vectors as columns.
    M = [[basis[j][i] for j in range(4)] for i in range(4)]
    rows = []
    for r in range(4):
        rhs = [images[j][r] for j in range(4)]
        # solve M^T? We need T row r: sum_c T[r][c] * basis_j[c] = images_j[r]
        sol = polytope.solve_square([list(b) for b in basis], rhs)
        if sol is None:
            return None
        rows.append(sol)
    del M
    return rows


def apply_T(T, p):
    return tuple(sum(T[r][c] * p[c] for c in range(4)) for r in range(4))


def gl_equivalent(V, W):
    """An invertible T with T(V) = W as sets, or None. |V| == |W| assumed.

    T is linear, so it maps antipodal pairs to antipodal pairs; the search
    runs over ordered 4-tuples of W's pair representatives with sign choices,
    the first sign fixed (+-T are equivalent since W = -W).  Exhaustive: any
    valid T sends the chosen basis pairs of V to 4 distinct pairs of W with
    some signs, and one of +-T has first sign +1."""
    vp = sorted({normalize(tuple(int(x) for x in p)) for p in V})
    wp = sorted({normalize(tuple(int(x) for x in p)) for p in W})
    # fix a linearly independent 4-subset of V's pair reps
    basis = None
    for combo in itertools.combinations(vp, 4):
        if polytope.rank([list(p) for p in combo]) == 4:
            basis = [tuple(Q(x) for x in p) for p in combo]
            break
    assert basis is not None
    Wset = set(W)
    for reps in itertools.permutations(wp, 4):
        for signs in itertools.product((1, -1), repeat=3):
            images = [tuple(Q(x) for x in reps[0])] + [
                tuple(Q(s * x) for x in r) for s, r in zip(signs, reps[1:])]
            T = solve_T(basis, images)
            if T is None:
                continue
            if polytope.rank([list(r) for r in T]) != 4:
                continue        # singular: would collapse V into W
            img = {apply_T(T, p) for p in V}
            if img == Wset:     # true bijection of vertex sets
                return T
    return None


_HM = None


def classify(mask):
    global _HM
    if _HM is None:
        _HM = hanner_masks()
    hm = _HM
    c = canon(mask)
    for m, nv in hm:
        if canon(m) == c:
            return f"B4-equivalent to harness Hanner ({nv} vertices, mask {m:x})"
    V = mask_points(mask)
    for m, nv in hm:
        if nv != len(V):
            continue
        T = gl_equivalent(V, mask_points(m))
        if T is not None:
            rows = "; ".join(",".join(str(x) for x in r) for r in T)
            return f"GL-equivalent to harness Hanner mask {m:x} via T=[{rows}]"
    return ("NOT-HANNER: not B4-equivalent and no linear map onto any "
            "harness Hanner vertex set exists (exhaustive basis-image search)")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    if sys.argv[1] == "list":
        seen = {}
        for m, nv in hanner_masks():
            seen.setdefault(canon(m), []).append((m, nv))
        print(f"{len(hanner_masks())} harness Hanner vertex sets, "
              f"{len(seen)} B4-orbits:")
        for c, members in sorted(seen.items()):
            print(f"  canon {c:x} ({members[0][1]} vertices): "
                  + " ".join(f"{m:x}" for m, _ in members))
        return 0
    if sys.argv[1] == "match":
        for h in sys.argv[2:]:
            print(f"{h}: {classify(int(h, 16))}")
        return 0
    if sys.argv[1] == "matchfile":
        with open(sys.argv[2]) as f:
            for line in f:
                fs = line.split()
                if fs and not fs[0].startswith("#"):
                    print(f"{fs[0]}: {classify(int(fs[0], 16))}")
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main())
