#!/usr/bin/env python3
"""Exact rational polytope arithmetic: facets, polarity, and volume, for the
volume product P(K) = vol(K) * vol(K-polar) of a centrally symmetric polytope.

Everything here is exact over Q (`fractions.Fraction`).  The question this
harness serves is a comparison of P(K) against a rational threshold, so a
floating-point volume cannot support a claim -- it can only propose a body
worth computing exactly.

Representations
---------------
A polytope is a list of points ("V-representation"), each point a tuple of
Fractions of the same length d.  The hull is conv(vertices); redundant points
inside the hull are allowed and are ignored by every routine here.  The
"H-representation" is a list of (a, b) meaning {x : <a,x> <= b}, with a an
integer-primitive OUTWARD normal.

Facet enumeration (brute force).  A hyperplane spanned by d affinely
independent vertices is a facet iff every other vertex lies weakly on one
side.  Enumerating all d-subsets is O(V^d) and is the honest choice at this
scale: it is exact, it has no degenerate-input branches, and for the vertex
counts of interest (a few dozen in R^4) it is fast enough that a subtler hull
algorithm would only add ways to be wrong.

Volume (divergence theorem, by recursion on dimension).  For a polytope P with
facets F in hyperplanes <a,x> = b, with a outward,

    vol_d(P) = (1/d) * sum_F  b_F * vol_{d-1}(F) / ||a_F||.

The norm is irrational, but the projection pi_j that deletes any coordinate j
with a_j != 0 is injective on the facet's hyperplane and scales (d-1)-volume by
exactly |a_j| / ||a||, so

    vol_d(P) = (1/d) * sum_F  b_F * vol_{d-1}(pi_j F) / |a_{F,j}|,

which is rational throughout.  b_F is signed, so the origin need not lie inside
P.  Recursion bottoms out at d = 1, where the volume is max - min.

Polarity.  For K = conv(V) with 0 in the interior, K-polar is {y : <v,y> <= 1
for all v in V}, i.e. one halfspace per vertex.  Its vertices are recovered by
solving every d-subset of those constraints with equality and keeping the
solutions that satisfy the rest -- again exact, again brute force.

Usage:
  polytope.py --selftest              run the validation suite
  polytope.py --product cube 4        volume product of a named body in R^n
                                      (cube, cross, 24cell), exact
  polytope.py --hanner 4              every Hanner polytope in R^n, with its
                                      vertex count and exact volume product
"""

import argparse
import itertools
import sys
from fractions import Fraction
from math import gcd, factorial

Q = Fraction


# ----------------------------------------------------------------------
# exact linear algebra
# ----------------------------------------------------------------------

def rref(rows):
    """Reduced row echelon form over Q.  Returns (matrix, pivot columns)."""
    m = [list(r) for r in rows]
    if not m:
        return m, []
    ncols = len(m[0])
    pivots = []
    r = 0
    for c in range(ncols):
        p = next((i for i in range(r, len(m)) if m[i][c] != 0), None)
        if p is None:
            continue
        m[r], m[p] = m[p], m[r]
        inv = Q(1) / m[r][c]
        m[r] = [x * inv for x in m[r]]
        for i in range(len(m)):
            if i != r and m[i][c] != 0:
                f = m[i][c]
                m[i] = [x - f * y for x, y in zip(m[i], m[r])]
        pivots.append(c)
        r += 1
        if r == len(m):
            break
    return m, pivots


def rank(rows):
    return len(rref(rows)[1])


def nullspace_line(rows, ncols):
    """A nonzero a with <row, a> = 0 for every row, when the null space is
    exactly one-dimensional.  Returns None otherwise (degenerate input)."""
    m, pivots = rref(rows)
    if ncols - len(pivots) != 1:
        return None
    free = next(c for c in range(ncols) if c not in pivots)
    a = [Q(0)] * ncols
    a[free] = Q(1)
    for i, c in enumerate(pivots):
        a[c] = -m[i][free]
    return tuple(a)


def solve_square(rows, rhs):
    """Unique solution of a square system, or None if singular."""
    n = len(rows)
    aug = [list(r) + [v] for r, v in zip(rows, rhs)]
    m, pivots = rref(aug)
    if pivots != list(range(n)):
        return None
    return tuple(m[i][n] for i in range(n))


def primitive(vec):
    """Scale a rational vector to integer entries with gcd 1, sign preserved."""
    den = 1
    for x in vec:
        den = den * x.denominator // gcd(den, x.denominator)
    ints = [int(x * den) for x in vec]
    g = 0
    for x in ints:
        g = gcd(g, abs(x))
    if g == 0:
        return tuple(ints), Q(den)
    return tuple(x // g for x in ints), Q(den, g)


# ----------------------------------------------------------------------
# facets
# ----------------------------------------------------------------------

def dimension(points):
    """Affine dimension of a point set."""
    if not points:
        return -1
    base = points[0]
    return rank([tuple(x - y for x, y in zip(p, base)) for p in points[1:]])


def facets(vertices):
    """All facets of conv(vertices), as (a, b, on) with a an outward
    integer-primitive normal, b the offset, and `on` the frozenset of indices
    of the given points lying in the hyperplane.

    Raises ValueError if the hull is not full-dimensional: every routine here
    assumes a body, and a silently-flat input would return volume 0 rather
    than an error.
    """
    d = len(vertices[0])
    if dimension(vertices) != d:
        raise ValueError(f"point set is not {d}-dimensional")
    seen = {}
    for combo in itertools.combinations(range(len(vertices)), d):
        base = vertices[combo[0]]
        diffs = [
            tuple(x - y for x, y in zip(vertices[i], base)) for i in combo[1:]
        ]
        a = nullspace_line(diffs, d)
        if a is None:
            continue
        b = sum(x * y for x, y in zip(a, base))
        pos = neg = False
        for v in vertices:
            s = sum(x * y for x, y in zip(a, v)) - b
            if s > 0:
                pos = True
            elif s < 0:
                neg = True
            if pos and neg:
                break
        if pos and neg:
            continue
        if pos:  # normal points inward; flip it
            a = tuple(-x for x in a)
            b = -b
        ai, scale = primitive(a)
        bi = b * scale
        if ai in seen:
            continue
        on = frozenset(
            i for i, v in enumerate(vertices)
            if sum(x * y for x, y in zip(ai, v)) == bi
        )
        seen[ai] = (ai, bi, on)
    return list(seen.values())


# ----------------------------------------------------------------------
# volume
# ----------------------------------------------------------------------

def volume(vertices):
    """Exact d-dimensional volume of conv(vertices)."""
    d = len(vertices[0])
    if d == 1:
        xs = [v[0] for v in vertices]
        return max(xs) - min(xs)
    total = Q(0)
    for a, b, on in facets(vertices):
        j = next(k for k in range(d) if a[k] != 0)
        face = [
            tuple(x for k, x in enumerate(vertices[i]) if k != j) for i in on
        ]
        total += b * volume(face) / abs(a[j])
    return total / d


# ----------------------------------------------------------------------
# polarity
# ----------------------------------------------------------------------

def polar_halfspaces(vertices):
    """H-representation of the polar body: <v, y> <= 1 for every vertex v."""
    return [(tuple(v), Q(1)) for v in vertices]


def polar(vertices):
    """Vertices of the polar body.  Requires 0 in the interior of the hull,
    which is what makes the polar bounded; violating it raises."""
    d = len(vertices[0])
    for a, b, _ in facets(vertices):
        if b <= 0:
            raise ValueError("origin is not interior; the polar is unbounded")
    out = set()
    for combo in itertools.combinations(range(len(vertices)), d):
        y = solve_square([vertices[i] for i in combo], [Q(1)] * d)
        if y is None:
            continue
        if all(sum(x * z for x, z in zip(v, y)) <= 1 for v in vertices):
            out.add(y)
    return sorted(out)


def volume_product(vertices):
    """P(K) = vol(K) * vol(K-polar), exact."""
    return volume(vertices) * volume(polar(vertices))


def is_symmetric(vertices):
    """K = -K, tested on the point set."""
    s = {tuple(v) for v in vertices}
    return all(tuple(-x for x in v) in s for v in s)


# ----------------------------------------------------------------------
# named bodies
# ----------------------------------------------------------------------

def segment():
    return [(Q(-1),), (Q(1),)]


def cube(n):
    return [tuple(Q(s) for s in signs)
            for signs in itertools.product((-1, 1), repeat=n)]


def cross(n):
    out = []
    for i in range(n):
        for s in (-1, 1):
            v = [Q(0)] * n
            v[i] = Q(s)
            out.append(tuple(v))
    return out


def cell24():
    """The 24-cell: all (+-1, +-1, 0, 0) up to coordinate position.  A
    self-dual symmetric 4-polytope that is not a Hanner polytope, so it is a
    useful non-trivial value to compute."""
    out = []
    for i, j in itertools.combinations(range(4), 2):
        for si in (-1, 1):
            for sj in (-1, 1):
                v = [Q(0)] * 4
                v[i], v[j] = Q(si), Q(sj)
                out.append(tuple(v))
    return out


def linf_sum(p, q):
    """K x L, the l-infinity (product) sum.  Polar-dual to the l1 sum."""
    return [tuple(list(u) + list(v)) for u in p for v in q]


def l1_sum(p, q):
    """conv(K x {0}, {0} x L), the l1 sum."""
    dp, dq = len(p[0]), len(q[0])
    return ([tuple(list(u) + [Q(0)] * dq) for u in p]
            + [tuple([Q(0)] * dp + list(v)) for v in q])


def hanner_polytopes(n):
    """Every Hanner polytope in R^n, as generated by the definition: segments
    combined by l1 and l-infinity sums in all ways.  Duplicates that coincide
    as point sets are removed; linearly equivalent copies are NOT identified,
    since deciding linear equivalence is a separate problem from generating
    the family."""
    if n == 1:
        return [tuple(segment())]
    out = set()
    for k in range(1, n):
        for left in hanner_polytopes(k):
            for right in hanner_polytopes(n - k):
                for combine in (l1_sum, linf_sum):
                    out.add(tuple(sorted(combine(list(left), list(right)))))
    return sorted(out, key=lambda p: (len(p), p))


BODIES = {"cube": cube, "cross": cross, "24cell": lambda n: cell24()}


# ----------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------

def selftest():
    ok = True

    def check(label, got, want, show=True):
        nonlocal ok
        good = got == want
        ok &= good
        detail = f"{got} (expect {want})" if show else f"{len(got)} points"
        print(f"[{label}] {detail} {'OK' if good else 'FAIL'}")
        return good

    # 1. Volumes of the two bodies whose volumes are elementary.
    for n in range(1, 5):
        check(f"vol cube R^{n}", volume(cube(n)), Q(2) ** n)
        check(f"vol cross R^{n}", volume(cross(n)), Q(2 ** n, factorial(n)))

    # 2. Polarity exchanges them, and is an involution.
    for n in range(1, 5):
        check(f"polar cube R^{n}", set(polar(cube(n))), set(cross(n)), False)
        check(f"polar cross R^{n}", set(polar(cross(n))), set(cube(n)), False)
        check(f"bipolar cube R^{n}", set(polar(polar(cube(n)))), set(cube(n)),
              False)

    # 3. Every Hanner polytope has volume product 4^n/n!.
    for n in range(1, 5):
        want = Q(4 ** n, factorial(n))
        family = hanner_polytopes(n)
        counts = sorted({len(p) for p in family})
        bad = [p for p in family if volume_product(list(p)) != want]
        ok &= not bad
        print(f"[hanner R^{n}] {len(family)} generated, vertex counts {counts}, "
              f"all P = {want} {'OK' if not bad else 'FAIL'}")

    # 4. The volume product is a linear invariant.  A unimodular map changes
    #    every intermediate quantity and must leave the product alone.
    t = [[1, 1, 0, 0], [0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1]]
    for name, body in (("cube", cube(4)), ("cross", cross(4)),
                       ("24cell", cell24())):
        mapped = [tuple(sum(Q(t[i][k]) * v[k] for k in range(4))
                        for i in range(4)) for v in body]
        check(f"invariance {name}", volume_product(mapped),
              volume_product(body))

    # 5. A symmetric body that is not a Hanner polytope must not come out
    #    below the conjectured bound; if it did, the bug is here.
    p24 = volume_product(cell24())
    bound = Q(4 ** 4, factorial(4))
    good = p24 >= bound and is_symmetric(cell24())
    ok &= good
    print(f"[24cell] P = {p24} = {float(p24):.6f}, bound {bound} = "
          f"{float(bound):.6f} {'OK' if good else 'FAIL'}")

    print("SELFTEST", "PASSED" if ok else "FAILED")
    return ok


# ----------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--product", nargs=2, metavar=("BODY", "N"),
                    help=f"one of {sorted(BODIES)}, and a dimension")
    ap.add_argument("--hanner", type=int, metavar="N")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)
    if args.product:
        name, n = args.product[0], int(args.product[1])
        body = BODIES[name](n)
        p = volume_product(body)
        print(f"{name} in R^{len(body[0])}: vol = {volume(body)}, "
              f"vol polar = {volume(polar(body))}, P = {p} = {float(p):.10f}")
        return
    if args.hanner:
        n = args.hanner
        for poly in hanner_polytopes(n):
            p = volume_product(list(poly))
            print(f"{len(poly):3d} vertices  P = {p} = {float(p):.10f}")
        print(f"conjectured minimum 4^{n}/{n}! = {Q(4 ** n, factorial(n))}")
        return
    ap.print_help()


if __name__ == "__main__":
    main()
