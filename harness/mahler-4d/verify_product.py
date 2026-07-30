#!/usr/bin/env python3
"""Adversarial re-computation of a claimed volume product.

`polytope.py` computes volume by the divergence theorem: a sum over facets of
projected lower-dimensional volumes.  Re-running it cannot catch an error in
that decomposition, so this tool computes the same number a different way and
compares.

The independent method here is Fubini slicing.  For a d-polytope P, let A(t) be
the (d-1)-volume of the slice P ∩ {x_d = t}.  Between consecutive vertex
heights A is a polynomial of degree at most d-1, because the slice is the image
of a fixed combinatorial cross-section whose vertices move affinely in t.  So
sampling A at d interior points of each height interval determines it, and the
exact integral of the interpolating polynomial gives the exact contribution of
that interval.  The slice's own volume is obtained by the same slicing
recursion, one dimension down, bottoming out at d = 1.

The slice vertices are the intersections of the edges of P with the hyperplane,
which is why edges are computed here: two vertices span an edge exactly when
the facets containing both have normals of rank d-1.

What is and is not independent, stated rather than implied.  The decomposition
is independent at every dimension -- no facet-projection sum is used anywhere
in this file's volume routine.  Both routines do share facet enumeration, so a
facet-enumeration bug would fool both; that shared component is instead checked
structurally, by the duality identities below, which any error in the facet set
would break:

  * the facets of the polar body are exactly the vertices of the original,
  * the polar of the polar is the original,
  * the volume product is unchanged by a unimodular map.

Usage:
  verify_product.py --selftest             run the validation suite
  verify_product.py --check CLAIM.json     re-verify a claim record
  verify_product.py --body cube 4          re-verify a named body

A claim record is JSON: {"vertices": [[...], ...], "volume": "a/b",
"volume_polar": "a/b", "product": "a/b"}, rationals as "numerator/denominator"
strings.  Any field present is checked; missing fields are computed and
reported.
"""

import argparse
import json
import os
import sys
from fractions import Fraction
from math import factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import polytope
from polytope import Q, facets, polar, primitive, rank, volume


# ----------------------------------------------------------------------
# exact polynomial interpolation and integration
# ----------------------------------------------------------------------

def poly_mul(p, q):
    out = [Q(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def poly_integral(p, lo, hi):
    """Definite integral of a polynomial given as ascending coefficients."""
    total = Q(0)
    for k, c in enumerate(p):
        total += c * (hi ** (k + 1) - lo ** (k + 1)) / (k + 1)
    return total


def lagrange_integral(nodes, values, lo, hi):
    """Exact integral over [lo, hi] of the polynomial through (nodes, values)."""
    total = Q(0)
    for i, xi in enumerate(nodes):
        basis, denom = [Q(1)], Q(1)
        for j, xj in enumerate(nodes):
            if j != i:
                basis = poly_mul(basis, [-xj, Q(1)])
                denom *= (xi - xj)
        total += values[i] * poly_integral(basis, lo, hi) / denom
    return total


# ----------------------------------------------------------------------
# combinatorics needed for slicing
# ----------------------------------------------------------------------

def vertices_and_edges(points):
    """Extreme points of conv(points), and the pairs spanning an edge.

    A point is a vertex when the normals of the facets through it span R^d;
    two vertices span an edge when the normals of their common facets have
    rank d-1, so the smallest face containing both is one-dimensional.
    """
    d = len(points[0])
    fs = facets(points)
    normals = [f[0] for f in fs]
    through = [
        {k for k, f in enumerate(fs) if i in f[2]} for i in range(len(points))
    ]
    verts = [i for i in range(len(points))
             if rank([normals[k] for k in through[i]]) == d]
    edges = []
    for a in range(len(verts)):
        for b in range(a + 1, len(verts)):
            i, j = verts[a], verts[b]
            common = through[i] & through[j]
            if common and rank([normals[k] for k in common]) == d - 1:
                edges.append((i, j))
    return verts, edges


# ----------------------------------------------------------------------
# volume by slicing
# ----------------------------------------------------------------------

def slice_points(points, edges, t):
    """Vertices of the slice {x_d = t}, projected to R^(d-1).

    t must be strictly between two consecutive vertex heights, so no vertex of
    the polytope lies in the hyperplane and every slice vertex comes from an
    edge crossing it.
    """
    out = set()
    for i, j in edges:
        u, v = points[i], points[j]
        hu, hv = u[-1], v[-1]
        if (hu < t < hv) or (hv < t < hu):
            s = (t - hu) / (hv - hu)
            out.add(tuple(a + s * (b - a) for a, b in zip(u[:-1], v[:-1])))
    return sorted(out)


def volume_by_slicing(points):
    """Exact d-volume of conv(points), by Fubini in the last coordinate."""
    d = len(points[0])
    if d == 1:
        xs = [p[0] for p in points]
        return max(xs) - min(xs)
    verts, edges = vertices_and_edges(points)
    heights = sorted({points[i][-1] for i in verts})
    total = Q(0)
    for lo, hi in zip(heights, heights[1:]):
        # d interior nodes: A is a polynomial of degree at most d-1 here.
        nodes = [lo + (hi - lo) * Q(k + 1, d + 1) for k in range(d)]
        values = [volume_by_slicing(slice_points(points, edges, t))
                  for t in nodes]
        total += lagrange_integral(nodes, values, lo, hi)
    return total


# ----------------------------------------------------------------------
# structural duality checks (they cover the shared facet enumeration)
# ----------------------------------------------------------------------

def facets_match_vertices(points):
    """The facets of the polar body are exactly the vertices of the original."""
    pol = polar(points)
    got = {(a, b) for a, b, _ in facets(pol)}
    # only the extreme points of the original correspond to facets
    verts, _ = vertices_and_edges(points)
    want = set()
    for i in verts:
        a, scale = primitive(points[i])
        want.add((a, Q(1) * scale))
    return got == want


def bipolar_is_identity(points):
    verts, _ = vertices_and_edges(points)
    return set(polar(polar(points))) == {points[i] for i in verts}


UNIMODULAR = {
    2: [[1, 1], [0, 1]],
    3: [[1, 1, 0], [0, 1, 1], [0, 0, 1]],
    4: [[1, 1, 0, 0], [0, 1, 0, 1], [0, 0, 1, 1], [0, 0, 0, 1]],
}


def apply_map(points, matrix):
    d = len(points[0])
    return [tuple(sum(Q(matrix[i][k]) * p[k] for k in range(d))
                  for i in range(d)) for p in points]


# ----------------------------------------------------------------------
# claim checking
# ----------------------------------------------------------------------

def report(points, claimed=None):
    """Recompute a body both ways and print a verdict.  Returns True if
    everything agreed and nothing fell below the conjectured bound."""
    d = len(points[0])
    ok = True
    v_div = volume(points)
    v_sli = volume_by_slicing(points)
    agree = v_div == v_sli
    ok &= agree
    print(f"  vol           divergence {v_div}  slicing {v_sli}  "
          f"{'AGREE' if agree else 'MISMATCH'}")

    pol = polar(points)
    p_div = volume(pol)
    p_sli = volume_by_slicing(pol)
    agree = p_div == p_sli
    ok &= agree
    print(f"  vol polar     divergence {p_div}  slicing {p_sli}  "
          f"{'AGREE' if agree else 'MISMATCH'}")

    for label, got in (("facets of polar == vertices",
                        facets_match_vertices(points)),
                       ("polar of polar == body", bipolar_is_identity(points))):
        ok &= got
        print(f"  {label}: {'OK' if got else 'FAIL'}")

    product = v_div * p_div
    mapped = apply_map(points, UNIMODULAR[d])
    inv = volume(mapped) * volume(polar(mapped))
    agree = inv == product
    ok &= agree
    print(f"  unimodular invariance: {inv} {'OK' if agree else 'FAIL'}")

    bound = Q(4 ** d, factorial(d))
    print(f"  P = {product} = {float(product):.10f}   bound 4^{d}/{d}! = "
          f"{bound} = {float(bound):.10f}")
    if product < bound:
        print("  *** BELOW THE CONJECTURED BOUND -- UNVERIFIED, ESCALATE. "
              "Check central symmetry and convex position before recording. ***")
        ok = False

    if claimed:
        for key, got in (("volume", v_div), ("volume_polar", p_div),
                         ("product", product)):
            if key in claimed:
                want = Fraction(claimed[key])
                agree = want == got
                ok &= agree
                print(f"  claimed {key} = {want}: "
                      f"{'CONFIRMED' if agree else 'CONTRADICTED, got ' + str(got)}")
    return ok


def check_claim(path):
    data = json.load(open(path))
    points = [tuple(Fraction(str(x)) for x in v) for v in data["vertices"]]
    print(f"{path}: {len(points)} points in R^{len(points[0])}, "
          f"symmetric={polytope.is_symmetric(points)}")
    return report(points, data)


# ----------------------------------------------------------------------

def selftest():
    ok = True
    cases = [
        ("cube R^2", polytope.cube(2)),
        ("cross R^3", polytope.cross(3)),
        ("cube R^3", polytope.cube(3)),
        ("cube R^4", polytope.cube(4)),
        ("cross R^4", polytope.cross(4)),
        ("24cell R^4", polytope.cell24()),
    ]
    for name, body in cases:
        print(f"[{name}]")
        ok &= report(body)

    # Every Hanner polytope in R^3 must survive both routines, not just the
    # two extreme ones.
    for poly in polytope.hanner_polytopes(3):
        body = list(poly)
        if volume(body) * volume(polar(body)) != volume_by_slicing(body) * \
                volume_by_slicing(polar(body)):
            ok = False
            print(f"[hanner R^3] MISMATCH on {len(body)}-vertex body")
    print(f"[hanner R^3] {len(polytope.hanner_polytopes(3))} bodies, both "
          f"routines agree {'OK' if ok else 'FAIL'}")

    # A deliberately wrong claim must be contradicted, or this tool would
    # confirm anything handed to it.
    print("[negative control] claim vol(cube R^3) = 7")
    got = report(polytope.cube(3), {"volume": "7"})
    ok &= not got
    print(f"[negative control] rejected as expected: {'OK' if not got else 'FAIL'}")

    print("SELFTEST", "PASSED" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--check", metavar="CLAIM.json")
    ap.add_argument("--body", nargs=2, metavar=("BODY", "N"))
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)
    if args.check:
        sys.exit(0 if check_claim(args.check) else 1)
    if args.body:
        body = polytope.BODIES[args.body[0]](int(args.body[1]))
        print(f"[{args.body[0]} R^{len(body[0])}]")
        sys.exit(0 if report(body) else 1)
    ap.print_help()


if __name__ == "__main__":
    main()
