#!/usr/bin/env python3
"""Adversarial re-computation of a Crouzeix ratio enclosure, by root isolation.

`ratio.py` bounds every eigenvalue it needs by bisecting on an exact
definiteness test.  That is one idea used everywhere in it: the norm, each
support value, and therefore the numerical range all rest on the same LDL
factorization.  A sign error there would be invisible to a re-run.

So this tool never factors anything.  It forms the characteristic polynomial
of each Hermitian matrix by the Faddeev-LeVerrier recursion, which is a
sequence of matrix products and traces, and then isolates its largest real
root with a Sturm chain, which counts sign variations.  Both steps are exact
over Q and neither shares a line of reasoning with a definiteness test.  Two
routes to the same eigenvalue that agree are worth considerably more than one
route run twice.

The characteristic polynomial is itself checked before it is trusted, by
evaluating it at rational points and comparing against a determinant computed
by elimination -- otherwise an error in the recursion would propagate into
every root this tool isolates.

Usage:
  verify_ratio.py --selftest        run the validation suite
  verify_ratio.py --jordan N        re-verify the N x N Jordan block
"""

import argparse
import sys
import os
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ratio
from ratio import Cx, Q, adjoint, cx, diagonal, eye, hermitian_part, jordan, \
    matadd, matmul, scalar, zeros


# ----------------------------------------------------------------------
# rational polynomials, ascending coefficients
# ----------------------------------------------------------------------

def trim(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def padd(a, b):
    out = [Q(0)] * max(len(a), len(b))
    for i, c in enumerate(a):
        out[i] += c
    for i, c in enumerate(b):
        out[i] += c
    return trim(out)


def pneg(a):
    return [-c for c in a]


def pmul(a, b):
    if not a or not b:
        return []
    out = [Q(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return trim(out)


def pdivmod(a, b):
    a = list(a)
    q = [Q(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and a:
        shift = len(a) - len(b)
        factor = a[-1] / b[-1]
        q[shift] = factor
        for i, c in enumerate(b):
            a[i + shift] -= factor * c
        trim(a)
    return trim(q), a


def pderiv(a):
    return trim([c * i for i, c in enumerate(a)][1:])


def pgcd(a, b):
    a, b = list(a), list(b)
    while b:
        a, b = b, pdivmod(a, b)[1]
    return [c / a[-1] for c in a] if a else a


def peval(a, x):
    out = Q(0)
    for c in reversed(a):
        out = out * x + c
    return out


# ----------------------------------------------------------------------
# characteristic polynomial and its own check
# ----------------------------------------------------------------------

def charpoly(h):
    """det(lambda I - H) by Faddeev-LeVerrier, as ascending rational
    coefficients.  H is Hermitian, so the coefficients are real; that they
    come out real is asserted rather than assumed."""
    n = len(h)
    coeffs = [None] * (n + 1)
    coeffs[n] = Q(1)
    m = eye(n)
    for k in range(1, n + 1):
        if k > 1:
            m = matadd(matmul(h, m), scalar(coeffs[n - k + 1], eye(n)))
        prod = matmul(h, m)
        tr = Cx(0)
        for i in range(n):
            tr = tr + prod[i][i]
        if tr.im != 0:
            raise ValueError("non-real trace from a Hermitian matrix")
        coeffs[n - k] = -tr.re / k
    return coeffs


def determinant(m):
    """det of a complex-rational matrix, by elimination.  Used only to check
    the characteristic polynomial, so it must not share its method."""
    a = [row[:] for row in m]
    n = len(a)
    det = Cx(1)
    for c in range(n):
        p = next((r for r in range(c, n) if not (a[r][c] == Cx(0))), None)
        if p is None:
            return Cx(0)
        if p != c:
            a[c], a[p] = a[p], a[c]
            det = det * Cx(-1)
        det = det * a[c][c]
        inv = Cx(1) / a[c][c]
        a[c] = [x * inv for x in a[c]]
        for r in range(c + 1, n):
            if not (a[r][c] == Cx(0)):
                f = a[r][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return det


def charpoly_agrees(h, samples=(Q(0), Q(1), Q(-2), Q(7, 3))):
    """The polynomial must agree with an independently computed determinant."""
    p = charpoly(h)
    n = len(h)
    for x in samples:
        want = determinant(matadd(scalar(x, eye(n)), scalar(-1, h)))
        if want.im != 0 or want.re != peval(p, x):
            return False
    return True


# ----------------------------------------------------------------------
# Sturm chains
# ----------------------------------------------------------------------

def sturm_chain(p):
    """The Sturm chain of a rational polynomial, after making it square-free.

    Repeated roots would break the sign-variation count, so the polynomial is
    divided by gcd(p, p') first; that leaves the same root *set*, which is all
    a count of distinct roots needs.
    """
    g = pgcd(p, pderiv(p))
    if len(g) > 1:
        p = pdivmod(p, g)[0]
    chain = [p, pderiv(p)]
    while len(chain[-1]) > 1:
        rem = pdivmod(chain[-2], chain[-1])[1]
        if not rem:
            break
        chain.append(pneg(rem))
    return chain


def variations(chain, x):
    signs = []
    for p in chain:
        v = peval(p, x)
        if v != 0:
            signs.append(1 if v > 0 else -1)
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)


def roots_in(chain, a, b):
    """Number of distinct real roots in (a, b]."""
    return variations(chain, a) - variations(chain, b)


def cauchy_bound(p):
    """Every real root lies in [-B, B]."""
    lead = p[-1]
    return 1 + max((abs(c / lead) for c in p[:-1]), default=Q(0))


def largest_root(p, tol=Q(1, 10 ** 9)):
    """Enclosure of the largest real root, by bisection on a Sturm count."""
    chain = sturm_chain(p)
    b = cauchy_bound(p)
    lo, hi = -b, b
    if roots_in(chain, lo, hi) == 0:
        raise ValueError("no real roots, which a Hermitian matrix must have")
    while hi - lo > tol:
        mid = (lo + hi) / 2
        if roots_in(chain, mid, hi) > 0:
            lo = mid
        else:
            hi = mid
    return lo, hi


def lambda_max(h, tol=Q(1, 10 ** 9)):
    """Largest eigenvalue of a Hermitian matrix, via its characteristic
    polynomial.  Nothing here factors, bisects on definiteness, or iterates
    toward an eigenvector."""
    return largest_root(charpoly(h), tol)


# ----------------------------------------------------------------------
# re-verification
# ----------------------------------------------------------------------

def overlaps(a, b):
    return not (a[1] < b[0] or b[1] < a[0])


def verify(a, coeffs, dir_bound=3, tol=Q(1, 10 ** 9), pieces=64):
    """Recompute a ratio enclosure by the independent route and compare."""
    ok = True
    claim = ratio.crouzeix_ratio(a, coeffs, dir_bound, tol, pieces)
    if claim["degenerate"]:
        print("  reference reports a degenerate case; nothing to compare")
        return False

    m = ratio.poly_matrix(coeffs, a)
    gram = matmul(adjoint(m), m)
    good = charpoly_agrees(gram)
    ok &= good
    print(f"  characteristic polynomial agrees with the determinant: "
          f"{'OK' if good else 'FAIL'}")

    mine = lambda_max(gram, tol)
    theirs = claim["norm2"]
    good = overlaps(mine, theirs)
    ok &= good
    print(f"  ||p(A)||^2  definiteness [{float(theirs[0]):.9f}, "
          f"{float(theirs[1]):.9f}]  roots [{float(mine[0]):.9f}, "
          f"{float(mine[1]):.9f}]  {'AGREE' if good else 'MISMATCH'}")

    # Rebuild the outer enclosure of W(A) from independently computed support
    # values, then re-derive the bound on |p| over it.
    planes = []
    disagreements = 0
    for u in ratio.directions(dir_bound):
        h = hermitian_part(scalar(u.conj(), a))
        lo, hi = lambda_max(h, tol)
        _, ref_hi = ratio.lambda_max(h, tol)
        if not overlaps((lo, hi), (ref_hi - tol, ref_hi + tol)):
            disagreements += 1
        planes.append((u, hi))
    ok &= disagreements == 0
    print(f"  {len(planes)} support values, {disagreements} disagreeing with "
          f"the reference {'OK' if not disagreements else 'FAIL'}")

    verts = ratio.hull(ratio.polygon(planes))
    outer2 = ratio.max_abs_squared(coeffs, verts, pieces)
    good = outer2 > 0 and abs(outer2 - claim["outer_max2"]) <= \
        Q(1, 10 ** 6) * max(Q(1), outer2)
    ok &= good
    print(f"  outer bound on max|p|^2  reference {float(claim['outer_max2']):.9f}"
          f"  independent {float(outer2):.9f} "
          f"{'AGREE' if good else 'MISMATCH'}")

    lo = ratio.isqrt_lower(mine[0] / outer2)
    hi = ratio.isqrt_upper(mine[1] / claim["inner_max2"])
    good = overlaps((lo, hi), (claim["ratio_lo"], claim["ratio_hi"]))
    ok &= good
    print(f"  ratio  reference [{float(claim['ratio_lo']):.9f}, "
          f"{float(claim['ratio_hi']):.9f}]  independent [{float(lo):.9f}, "
          f"{float(hi):.9f}]  {'AGREE' if good else 'MISMATCH'}")

    refutes = mine[0] > 4 * outer2
    if refutes != claim["refutes"]:
        ok = False
        print("  MISMATCH on the refutation test itself")
    if refutes:
        print("  *** BOTH ROUTES CERTIFY A RATIO ABOVE 2. This would refute "
              "the conjecture. Check the matrix and polynomial are what you "
              "meant before recording anything. ***")
    return ok


# ----------------------------------------------------------------------

def selftest():
    ok = True

    def check(label, good, detail=""):
        nonlocal ok
        ok &= good
        print(f"[{label}] {detail} {'OK' if good else 'FAIL'}")

    # 1. The polynomial machinery, on answers that can be written down.
    p = pmul([-1, 1], [-2, 1])            # (x-1)(x-2) = x^2 - 3x + 2
    check("polynomial product", p == [Q(2), Q(-3), Q(1)], str(p))
    lo, hi = largest_root(p)
    check("largest root of (x-1)(x-2)", lo <= 2 <= hi,
          f"[{float(lo):.9f}, {float(hi):.9f}]")
    lo, hi = largest_root([-2, 0, 1])     # x^2 - 2, largest root sqrt(2)
    root2 = Q(14142135623, 10 ** 10)
    check("largest root of x^2 - 2", lo <= root2 <= hi,
          f"[{float(lo):.9f}, {float(hi):.9f}]")
    # A repeated root must not break the count.
    lo, hi = largest_root(pmul([-3, 1], pmul([-3, 1], [-1, 1])))
    check("repeated roots handled", lo <= 3 <= hi,
          f"[{float(lo):.9f}, {float(hi):.9f}]")

    # 2. Characteristic polynomials against an independent determinant.
    for name, h in (("diagonal", diagonal([3, -1, Q(1, 2)])),
                    ("swap", [[Cx(0), Cx(1)], [Cx(1), Cx(0)]]),
                    ("complex hermitian",
                     [[Cx(1), Cx(0, 1)], [Cx(0, -1), Cx(1)]])):
        check(f"charpoly of the {name} matrix", charpoly_agrees(h),
              str([str(c) for c in charpoly(h)]))

    # 3. The two eigenvalue routes must agree where both apply.
    for name, h in (("diagonal", diagonal([3, -1, Q(1, 2)])),
                    ("swap", [[Cx(0), Cx(1)], [Cx(1), Cx(0)]]),
                    ("complex hermitian",
                     [[Cx(1), Cx(0, 1)], [Cx(0, -1), Cx(1)]])):
        mine, theirs = lambda_max(h), ratio.lambda_max(h)
        check(f"two routes agree on the {name} matrix", overlaps(mine, theirs),
              f"roots [{float(mine[0]):.9f}, {float(mine[1]):.9f}] vs "
              f"definiteness [{float(theirs[0]):.9f}, {float(theirs[1]):.9f}]")

    # 4. Full re-verification of the two literature cases.
    print("[re-verify the 2x2 Jordan block]")
    ok &= verify(jordan(2), [0, 1])
    print("[re-verify the 3x3 Jordan block]")
    ok &= verify(jordan(3), [0, 1])
    print("[re-verify a normal matrix]")
    ok &= verify(diagonal([Cx(1), Cx(0, 1), Cx(-1)]), [-1, 0, 1])

    # 5. Negative control.  A verifier that cannot disagree is decoration, so
    #    hand the comparison a deliberately wrong eigenvalue and check it is
    #    caught.
    h = diagonal([3, -1])
    wrong = (Q(5), Q(5) + Q(1, 10 ** 9))
    check("a wrong eigenvalue is not accepted",
          not overlaps(lambda_max(h), wrong),
          f"claimed 5, actual {float(lambda_max(h)[0]):.6f}")

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
        print(f"[re-verify the {args.jordan}x{args.jordan} Jordan block]")
        sys.exit(0 if verify(jordan(args.jordan), [0, 1]) else 1)
    ap.print_help()


if __name__ == "__main__":
    main()
