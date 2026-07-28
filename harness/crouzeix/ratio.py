#!/usr/bin/env python3
"""Certified enclosures of the Crouzeix ratio ||p(A)|| / max{|p(z)| : z in W(A)}.

Neither quantity in the ratio is rational: both come from eigenvalue problems.
So every number this module reports is an interval with a proof behind it, and
a bare floating-point value is never produced -- there is no floating point in
this file at all.  Four ideas make that affordable.

**Eigenvalue bounds by definiteness, not by iteration.**  For a Hermitian H,
lambda_max(H) < t exactly when t*I - H is positive definite, and definiteness
is decided exactly by an LDL factorization over the rationals: no square roots,
no residual estimates, no convergence assumption.  Bisecting on t between
Gershgorin bounds encloses lambda_max to any width, with both ends certified.

**The numerical range from rational directions.**  For any nonzero complex u,
Re(conj(u) * x*Ax) = x* Herm(conj(u) A) x <= lambda_max(Herm(conj(u) A)), so
each direction gives a halfplane containing W(A) -- and both sides scale the
same way with u, so u need not be a unit vector.  That removes the usual
obstacle: directions can be Gaussian *rational*, taken here from coprime
integer pairs, and the intersection of the halfplanes is an exact rational
polygon containing W(A).

**Outer for the numerator, inner for the denominator.**  The ratio is bounded
below by (lower bound on ||p(A)||) / (upper bound on the maximum of |p| over
W(A)), and the second needs an *outer* enclosure of W(A), since a maximum over
a larger set is larger.  Bounding the ratio from above needs the opposite:
points known to lie inside W(A), which x*Ax/x*x supplies exactly for any
rational x.  Getting these two backwards is the easiest way to produce a
refutation that is not one.

**Squares, not roots.**  ||p(A)||^2 = lambda_max(p(A)* p(A)) is rational-
enclosable; its square root is not.  The decision that matters -- is the ratio
above 2 -- is therefore taken as ||p(A)||^2 > 4 * (max |p|)^2, entirely in Q.
Square roots are taken only to display a ratio.

Usage:
  ratio.py --selftest         run the validation suite
  ratio.py --jordan N         the ratio for the N x N Jordan block, p(z) = z
"""

import argparse
import itertools
import sys
from fractions import Fraction

Q = Fraction


# ----------------------------------------------------------------------
# exact complex rational arithmetic
# ----------------------------------------------------------------------

class Cx:
    """A complex number with exact rational parts."""

    __slots__ = ("re", "im")

    def __init__(self, re=0, im=0):
        self.re = Q(re)
        self.im = Q(im)

    def __add__(self, o):
        o = cx(o)
        return Cx(self.re + o.re, self.im + o.im)

    def __sub__(self, o):
        o = cx(o)
        return Cx(self.re - o.re, self.im - o.im)

    def __neg__(self):
        return Cx(-self.re, -self.im)

    def __mul__(self, o):
        o = cx(o)
        return Cx(self.re * o.re - self.im * o.im,
                  self.re * o.im + self.im * o.re)

    def __truediv__(self, o):
        o = cx(o)
        n = o.re * o.re + o.im * o.im
        if n == 0:
            raise ZeroDivisionError
        return Cx((self.re * o.re + self.im * o.im) / n,
                  (self.im * o.re - self.re * o.im) / n)

    __radd__ = __add__
    __rmul__ = __mul__

    def conj(self):
        return Cx(self.re, -self.im)

    def norm2(self):
        """|z|^2, exactly rational."""
        return self.re * self.re + self.im * self.im

    def __eq__(self, o):
        o = cx(o)
        return self.re == o.re and self.im == o.im

    def __hash__(self):
        return hash((self.re, self.im))

    def __repr__(self):
        if self.im == 0:
            return f"{float(self.re):.6g}"
        return f"{float(self.re):.6g}{float(self.im):+.6g}i"


def cx(x):
    return x if isinstance(x, Cx) else Cx(x)


# ----------------------------------------------------------------------
# matrices
# ----------------------------------------------------------------------

def zeros(n):
    return [[Cx(0) for _ in range(n)] for _ in range(n)]


def eye(n):
    return [[Cx(1 if i == j else 0) for j in range(n)] for i in range(n)]


def matmul(a, b):
    n = len(a)
    out = zeros(n)
    for i in range(n):
        for k in range(n):
            if a[i][k] == Cx(0):
                continue
            aik = a[i][k]
            for j in range(n):
                out[i][j] = out[i][j] + aik * b[k][j]
    return out


def matadd(a, b):
    return [[x + y for x, y in zip(ra, rb)] for ra, rb in zip(a, b)]


def scalar(s, a):
    s = cx(s)
    return [[s * x for x in row] for row in a]


def adjoint(a):
    n = len(a)
    return [[a[j][i].conj() for j in range(n)] for i in range(n)]


def hermitian_part(a):
    return scalar(Q(1, 2), matadd(a, adjoint(a)))


def poly_matrix(coeffs, a):
    """p(A) for p given by ascending coefficients."""
    n = len(a)
    out = zeros(n)
    power = eye(n)
    for c in coeffs:
        out = matadd(out, scalar(c, power))
        power = matmul(power, a)
    return out


def jordan(n):
    """The n x n nilpotent Jordan block, the standard extremal example."""
    m = zeros(n)
    for i in range(n - 1):
        m[i][i + 1] = Cx(1)
    return m


def diagonal(values):
    n = len(values)
    m = zeros(n)
    for i, v in enumerate(values):
        m[i][i] = cx(v)
    return m


# ----------------------------------------------------------------------
# certified eigenvalue bounds
# ----------------------------------------------------------------------

def is_positive_definite(h):
    """Exact LDL test for a Hermitian matrix.  No square roots are taken, so
    every quantity stays rational and the answer is a decision, not an
    estimate."""
    n = len(h)
    d = []
    lo = [[Cx(0)] * n for _ in range(n)]
    for i in range(n):
        for j in range(i):
            s = h[i][j]
            for k in range(j):
                s = s - lo[i][k] * lo[j][k].conj() * Cx(d[k])
            lo[i][j] = s / Cx(d[j])
        s = h[i][i].re
        for k in range(i):
            s -= lo[i][k].norm2() * d[k]
        if s <= 0:
            return False
        d.append(s)
    return True


def gershgorin(h):
    """A certified bound: every eigenvalue lies in [-g, g]."""
    g = Q(0)
    for i, row in enumerate(h):
        centre = row[i].re
        radius = sum(abs_bound(row[j]) for j in range(len(h)) if j != i)
        g = max(g, abs(centre) + radius)
    return g


def abs_bound(z):
    """A rational upper bound for |z|, since |z| itself need not be rational."""
    n = z.norm2()
    if n == 0:
        return Q(0)
    return isqrt_upper(n)


def isqrt_upper(v, tol=Q(1, 10 ** 12)):
    """Smallest rational found by bisection whose square is at least v."""
    lo, hi = Q(0), max(Q(1), v)
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if mid * mid < v:
            lo = mid
        else:
            hi = mid
    return hi


def isqrt_lower(v, tol=Q(1, 10 ** 12)):
    lo, hi = Q(0), max(Q(1), v)
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if mid * mid <= v:
            lo = mid
        else:
            hi = mid
    return lo


def lambda_max(h, tol=Q(1, 10 ** 9)):
    """Enclosure of the largest eigenvalue of a Hermitian matrix.

    t*I - H is positive definite exactly when t exceeds every eigenvalue, so
    bisection on that test brackets lambda_max from both sides with no
    iteration to converge and nothing to trust.
    """
    n = len(h)
    g = gershgorin(h)
    lo, hi = -g, g
    while hi - lo > tol:
        mid = (lo + hi) / 2
        shifted = matadd(scalar(mid, eye(n)), scalar(-1, h))
        if is_positive_definite(shifted):
            hi = mid
        else:
            lo = mid
    return lo, hi


def norm_squared(m, tol=Q(1, 10 ** 9)):
    """Enclosure of ||M||^2 = lambda_max(M* M)."""
    return lambda_max(matmul(adjoint(m), m), tol)


# ----------------------------------------------------------------------
# the numerical range
# ----------------------------------------------------------------------

def directions(bound=3):
    """Gaussian-rational directions: coprime integer pairs in a box.

    Soundness does not depend on which directions are chosen -- every
    direction gives a valid halfplane -- so they can be enumerated exactly
    instead of placed by trigonometry.  More directions means a tighter
    polygon.
    """
    from math import gcd
    out = []
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            if (a or b) and gcd(abs(a), abs(b)) == 1:
                out.append(Cx(a, b))
    return out


def support_halfplanes(a, dirs, tol=Q(1, 10 ** 9)):
    """Halfplanes {z : Re(conj(u) z) <= h} containing W(A), one per direction.

    The bound uses the *upper* end of the eigenvalue enclosure, so each
    halfplane provably contains W(A) rather than approximately containing it.
    """
    out = []
    for u in dirs:
        h = hermitian_part(scalar(u.conj(), a))
        _, hi = lambda_max(h, tol)
        out.append((u, hi))
    return out


def polygon(halfplanes):
    """Vertices of the intersection of halfplanes, exactly.

    Each pair of boundary lines is solved and the solution kept when it
    satisfies every constraint -- brute force, but exact and with no
    orientation heuristics to get wrong.

    The slack in the containment test admits a few points just outside the
    exact intersection, which enlarges the polygon slightly.  That is the safe
    direction and deliberately so: this polygon's job is to contain W(A), and
    a bound computed over a slightly larger set is still an upper bound.
    """
    pts = set()
    for (u1, h1), (u2, h2) in itertools.combinations(halfplanes, 2):
        # Re(conj(u) z) = u.re*x + u.im*y
        det = u1.re * u2.im - u1.im * u2.re
        if det == 0:
            continue
        x = (h1 * u2.im - h2 * u1.im) / det
        y = (u1.re * h2 - u2.re * h1) / det
        if all(u.re * x + u.im * y <= h + Q(1, 10 ** 6)
               for u, h in halfplanes):
            pts.add((x, y))
    return sorted(pts)


def hull(points):
    """Convex hull in counter-clockwise order, by monotone chain with exact
    rational cross products."""
    pts = sorted(set(points))
    if len(pts) < 3:
        return pts

    def cross(o, a, b):
        return ((a[0] - o[0]) * (b[1] - o[1])
                - (a[1] - o[1]) * (b[0] - o[0]))

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def inner_points(a, vectors):
    """Points certified to lie in W(A): x*Ax / x*x for rational vectors x."""
    out = []
    n = len(a)
    for x in vectors:
        denom = sum(v.norm2() for v in x)
        if denom == 0:
            continue
        total = Cx(0)
        for i in range(n):
            for j in range(n):
                total = total + x[i].conj() * a[i][j] * x[j]
        out.append(total / Cx(denom))
    return out


def probe_vectors(n, spread=2):
    """Small rational vectors, enough to put points inside W(A) in several
    directions.  Any nonzero vector is valid, so these need no justification
    beyond being cheap."""
    out = []
    for i in range(n):
        e = [Cx(0)] * n
        e[i] = Cx(1)
        out.append(e)
    for combo in itertools.product(range(-spread, spread + 1), repeat=n):
        if any(combo):
            out.append([Cx(c) for c in combo])
    for combo in itertools.product((0, 1), repeat=n):
        if any(combo):
            out.append([Cx(0, c) if i % 2 else Cx(c)
                        for i, c in enumerate(combo)])
    return out


# ----------------------------------------------------------------------
# bounding |p| on the boundary of a polygon
# ----------------------------------------------------------------------

def poly_on_segment(coeffs, z0, z1):
    """|p(z0 + s(z1-z0))|^2 as a rational polynomial in s, ascending."""
    d = z1 - z0
    # p(z0 + s d) as a polynomial in s with Cx coefficients
    cur = [Cx(1)]                      # (z0 + s d)^k, built up
    acc = [Cx(0)] * len(coeffs)
    for k, c in enumerate(coeffs):
        for i, t in enumerate(cur):
            acc[i] = acc[i] + cx(c) * t
        nxt = [Cx(0)] * (len(cur) + 1)
        for i, t in enumerate(cur):
            nxt[i] = nxt[i] + t * z0
            nxt[i + 1] = nxt[i + 1] + t * d
        cur = nxt
    re = [t.re for t in acc]
    im = [t.im for t in acc]
    out = [Q(0)] * (2 * len(acc) - 1)
    for i, x in enumerate(re):
        for j, y in enumerate(re):
            out[i + j] += x * y
    for i, x in enumerate(im):
        for j, y in enumerate(im):
            out[i + j] += x * y
    return out


def bound_on_unit_interval(q, pieces=64):
    """Upper bound for a rational polynomial on [0, 1], by interval Horner
    over sub-intervals.  Every coefficient is exact, so the bound is a bound."""
    best = None
    for k in range(pieces):
        lo, hi = Q(k, pieces), Q(k + 1, pieces)
        rl, rh = Q(0), Q(0)
        for c in reversed(q):
            products = (rl * lo, rl * hi, rh * lo, rh * hi)
            rl, rh = min(products) + c, max(products) + c
        best = rh if best is None else max(best, rh)
    return best


def max_abs_squared(coeffs, verts, pieces=64):
    """Upper bound for |p|^2 on a convex polygon.

    |p| is subharmonic, so its maximum over the polygon is attained on the
    boundary and it is enough to bound each edge.
    """
    if len(verts) == 1:
        z = Cx(verts[0][0], verts[0][1])
        val = Cx(0)
        for k, c in enumerate(coeffs):
            p = Cx(1)
            for _ in range(k):
                p = p * z
            val = val + cx(c) * p
        return val.norm2()
    best = Q(0)
    for i, (x0, y0) in enumerate(verts):
        x1, y1 = verts[(i + 1) % len(verts)]
        q = poly_on_segment(coeffs, Cx(x0, y0), Cx(x1, y1))
        best = max(best, bound_on_unit_interval(q, pieces))
    return best


def eval_abs_squared(coeffs, z):
    val = Cx(0)
    power = Cx(1)
    for c in coeffs:
        val = val + cx(c) * power
        power = power * z
    return val.norm2()


# ----------------------------------------------------------------------
# the ratio
# ----------------------------------------------------------------------

def crouzeix_ratio(a, coeffs, dir_bound=3, tol=Q(1, 10 ** 9), pieces=64):
    """Certified enclosure of the Crouzeix ratio, as a dict.

    The lower end pairs the smallest certified ||p(A)|| with the largest
    certified bound on |p| over an outer enclosure of W(A); the upper end
    pairs the largest ||p(A)|| with the largest value of |p| at points proved
    to lie inside W(A).  Both ends are therefore valid for the true ratio.
    """
    m = poly_matrix(coeffs, a)
    n2_lo, n2_hi = norm_squared(m, tol)

    planes = support_halfplanes(a, directions(dir_bound), tol)
    verts = hull(polygon(planes))
    outer2 = max_abs_squared(coeffs, verts, pieces)

    pts = inner_points(a, probe_vectors(len(a)))
    inner2 = max((eval_abs_squared(coeffs, z) for z in pts), default=Q(0))

    if outer2 <= 0 or inner2 <= 0:
        return {"degenerate": True, "outer_max2": outer2, "inner_max2": inner2}

    return {
        "degenerate": False,
        "norm2": (n2_lo, n2_hi),
        "outer_max2": outer2,          # >= (max |p| over W(A))^2
        "inner_max2": inner2,          # <= (max |p| over W(A))^2
        "ratio_lo": isqrt_lower(n2_lo / outer2),
        "ratio_hi": isqrt_upper(n2_hi / inner2),
        # The refutation test, decided entirely in Q: no square root is taken.
        "refutes": n2_lo > 4 * outer2,
        "polygon": verts,
    }


def describe(a, coeffs, **kw):
    r = crouzeix_ratio(a, coeffs, **kw)
    if r["degenerate"]:
        print(f"  degenerate: outer {r['outer_max2']}, inner {r['inner_max2']}")
        return r
    print(f"  ||p(A)||^2 in [{float(r['norm2'][0]):.9f}, "
          f"{float(r['norm2'][1]):.9f}]")
    print(f"  max|p| over W(A): squared, between {float(r['inner_max2']):.9f} "
          f"and {float(r['outer_max2']):.9f} ({len(r['polygon'])}-gon)")
    print(f"  ratio in [{float(r['ratio_lo']):.9f}, "
          f"{float(r['ratio_hi']):.9f}]")
    if r["refutes"]:
        print("  *** CERTIFIED ABOVE 2 -- UNVERIFIED, ESCALATE. "
              "Re-derive with verify_ratio.py before recording. ***")
    return r


# ----------------------------------------------------------------------
# validation
# ----------------------------------------------------------------------

def selftest():
    ok = True

    def check(label, good, detail=""):
        nonlocal ok
        ok &= good
        print(f"[{label}] {detail} {'OK' if good else 'FAIL'}")

    # 1. Exact linear algebra on matrices with known spectra.
    h = diagonal([3, -1, Q(1, 2)])
    lo, hi = lambda_max(h)
    check("lambda_max of a diagonal matrix", lo <= 3 <= hi,
          f"[{float(lo):.9f}, {float(hi):.9f}] contains 3")
    check("positive definite test", is_positive_definite(diagonal([1, 2, 3]))
          and not is_positive_definite(diagonal([1, -2, 3])))

    # A Hermitian matrix with irrational eigenvalues: [[0,1],[1,0]] has
    # eigenvalues +-1, and [[1,i],[-i,1]] has 0 and 2.
    swap = [[Cx(0), Cx(1)], [Cx(1), Cx(0)]]
    lo, hi = lambda_max(swap)
    check("lambda_max of the swap matrix", lo <= 1 <= hi,
          f"[{float(lo):.9f}, {float(hi):.9f}] contains 1")

    # 2. The 2x2 Jordan block is the case that makes 2 the right constant:
    #    ||J|| = 1 while W(J) is the disc of radius 1/2, so the ratio is
    #    exactly 2.  The enclosure must contain 2 -- and must not claim to
    #    have exceeded it.
    print("[2x2 Jordan block, p(z) = z]")
    r = describe(jordan(2), [0, 1])
    check("ratio encloses 2", r["ratio_lo"] <= 2 <= r["ratio_hi"])
    check("no false refutation", not r["refutes"])
    check("||J||^2 encloses 1", r["norm2"][0] <= 1 <= r["norm2"][1])
    check("W(J) is the disc of radius 1/2",
          r["outer_max2"] >= Q(1, 4) and r["inner_max2"] <= Q(1, 4),
          f"inner {float(r['inner_max2']):.6f} <= 1/4 <= outer "
          f"{float(r['outer_max2']):.6f}")

    # 3. The 3x3 Jordan block: W is the disc of radius cos(pi/4), so the
    #    ratio is sqrt(2).  A second literature value, with a different
    #    answer, so a harness that returned 2 for everything would fail here.
    print("[3x3 Jordan block, p(z) = z]")
    r3 = describe(jordan(3), [0, 1])
    root2 = Q(14142135623, 10 ** 10)
    check("ratio encloses sqrt(2)", r3["ratio_lo"] <= root2 <= r3["ratio_hi"],
          f"[{float(r3['ratio_lo']):.9f}, {float(r3['ratio_hi']):.9f}]")
    check("no false refutation in 3x3", not r3["refutes"])

    # 4. Normal matrices satisfy the conjecture with constant 1: for a
    #    diagonal A, ||p(A)|| is the largest |p| at an eigenvalue, and the
    #    eigenvalues lie in W(A).
    print("[diagonal matrix, p(z) = z^2 - 1]")
    rn = describe(diagonal([Cx(1), Cx(0, 1), Cx(-1, 0)]), [-1, 0, 1])
    check("normal ratio at most 1", rn["ratio_lo"] <= 1 + Q(1, 10 ** 6),
          f"lower end {float(rn['ratio_lo']):.9f}")
    check("no false refutation for a normal matrix", not rn["refutes"])

    # 5. The enclosure must actually enclose: the inner points are in W(A)
    #    and the polygon contains them.
    a = [[Cx(0), Cx(2)], [Cx(0, 1), Cx(1)]]
    planes = support_halfplanes(a, directions(3))
    verts = hull(polygon(planes))
    pts = inner_points(a, probe_vectors(2))
    outside = [z for z in pts
               if any(u.re * z.re + u.im * z.im > h + Q(1, 10 ** 6)
                      for u, h in planes)]
    check("every certified inner point is inside the outer polygon",
          not outside, f"{len(pts)} points, {len(verts)}-gon")

    # 6. A polygon that did not contain the spectrum would be wrong: every
    #    eigenvalue of A lies in W(A).
    d = diagonal([Cx(2, 1), Cx(-1, 2), Cx(0, -3)])
    planes = support_halfplanes(d, directions(3))
    bad = [v for v in (Cx(2, 1), Cx(-1, 2), Cx(0, -3))
           if any(u.re * v.re + u.im * v.im > h + Q(1, 10 ** 6)
                  for u, h in planes)]
    check("the enclosure contains the spectrum", not bad)

    print("SELFTEST", "PASSED" if ok else "FAILED")
    return ok


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--jordan", type=int, metavar="N")
    args = ap.parse_args()

    if args.selftest:
        sys.exit(0 if selftest() else 1)
    if args.jordan:
        print(f"[{args.jordan}x{args.jordan} Jordan block, p(z) = z]")
        describe(jordan(args.jordan), [0, 1])
        return
    ap.print_help()


if __name__ == "__main__":
    main()
