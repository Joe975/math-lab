#!/usr/bin/env python3
"""Closed-form analysis of two screened rank-3/rank-4 laminate shapes for the
three-phase isotropic lower bound, sigma=(1,2,5), f=(f1,f2,f3).

Everything here is exact (fractions.Fraction) and every closed form is
cross-checked against harness/three-phase-conductivity/laminate.py:effective()
on concrete rational parameter values (--selftest).  Optimisation (locating an
isotropic point, minimising over a free parameter) is done by bisection /
golden-section on exact Fraction-valued objective functions -- "optimising the
exact closed form", not a float screen -- though a golden-section search itself
uses floating comparisons to pick step sizes; every point it lands on is
evaluated exactly, and the final answer is verified against effective() at a
concrete rational point.

RANK-3 SHAPE (as literally scored by explore/tp_search.py; see selftest for
the exact tree and the note on the initial mislabelling in the human summary
that this script corrects -- normals are e1/e1/e2, not e1/e1/e2-outer as first
described):

    L12 = rank-1 laminate of p1 (fraction m_L) and p2, normal e1
    L13 = rank-1 laminate of p1 (fraction m_R) and p3, normal e2
    result = laminate of L12 (fraction m_top) and L13, normal e1

RANK-4 SHAPE (screened tree, explore/tp_search.py --rank 4):

    X = rank-1 laminate of p3 (fraction b) and p1, normal e2
    Z = rank-1 laminate of p1 (fraction d) and p3, normal e1
    Y = rank-1 laminate of Z (fraction c) and p2, normal e2
    result = laminate of X (fraction a) and Y, normal e1

Usage:
    python tp_shapes.py --selftest
    python tp_shapes.py --report
"""
from __future__ import annotations

import argparse
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "three-phase-conductivity"))
import laminate as L  # noqa: E402

S1, S2, S3 = F(1), F(2), F(5)


# ---------------------------------------------------------------------------
# elementary weighted means, m weights the FIRST argument
# ---------------------------------------------------------------------------

def harm(m: F, p: F, q: F) -> F:
    return 1 / (m / p + (1 - m) / q)


def arith(m: F, p: F, q: F) -> F:
    return m * p + (1 - m) * q


# ---------------------------------------------------------------------------
# RANK 3
# ---------------------------------------------------------------------------

def rank3_tree(m_top: F, m_L: F, m_R: F) -> dict:
    L12 = {"fraction": str(m_L), "normal": [1, 0], "layers": [{"phase": "p1"}, {"phase": "p2"}]}
    L13 = {"fraction": str(m_R), "normal": [0, 1], "layers": [{"phase": "p1"}, {"phase": "p3"}]}
    return {"fraction": str(m_top), "normal": [1, 0], "layers": [L12, L13]}


def rank3_from_fractions(m_top: F, f2: F, f3: F, sig=(S1, S2, S3)):
    """m_L, m_R implied by target f2, f3 (see write-up for the derivation)."""
    s1, s2, s3 = sig
    m_L = 1 - f2 / m_top
    m_R = 1 - f3 / (1 - m_top)
    return m_L, m_R


def rank3_tensor(m_top: F, f2: F, f3: F, sig=(S1, S2, S3)):
    """Closed-form (sigma11, sigma22) of the rank-3 shape as a function of the
    single free parameter m_top, with m_L, m_R eliminated via the fraction
    constraints f2 = m_top(1-m_L), f3 = (1-m_top)(1-m_R)."""
    s1, s2, s3 = sig
    m_L, m_R = rank3_from_fractions(m_top, f2, f3, sig)
    if not (0 < m_L < 1) or not (0 < m_R < 1):
        return None
    h12 = harm(m_L, s1, s2)          # L12 normal e1: T11 harmonic
    a13 = arith(m_R, s1, s3)         # L13 normal e2: T11 arithmetic
    a12 = arith(m_L, s1, s2)         # L12 normal e1: T22 arithmetic
    h13 = harm(m_R, s1, s3)          # L13 normal e2: T22 harmonic
    sigma11 = h12 * a13 / ((1 - m_top) * h12 + m_top * a13)   # harm(m_top; h12, a13)
    sigma22 = m_top * a12 + (1 - m_top) * h13
    return sigma11, sigma22


def rank3_gap(m_top: F, f2: F, f3: F, sig=(S1, S2, S3)):
    t = rank3_tensor(m_top, f2, f3, sig)
    if t is None:
        return None
    return t[0] - t[1]


def rank3_domain(f2: F, f3: F):
    return f2, 1 - f3


def rank3_isotropic_point(f2: F, f3: F, sig=(S1, S2, S3), iters=200, eps=F(1, 10**12)):
    """Exact bisection for the isotropic m_top in the domain (f2, 1-f3), if
    the gap changes sign there; else None (no isotropic member of this exact
    shape at these fractions)."""
    lo, hi = rank3_domain(f2, f3)
    lo, hi = lo + eps, hi - eps
    if lo >= hi:
        return None
    glo = rank3_gap(lo, f2, f3, sig)
    ghi = rank3_gap(hi, f2, f3, sig)
    if glo is None or ghi is None or (glo > 0) == (ghi > 0):
        return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        g = rank3_gap(mid, f2, f3, sig)
        if (g > 0) == (glo > 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# integer numerator polynomial (low->high) of the isotropy equation at
# f=(1/8,1/8,6/8), sigma=(1,2,5), CLEARED of the spurious factor m_top (an
# extraneous root introduced by clearing denominators at m_top=0) -- see
# selftest for the derivation and the exact-match check.
RANK3_QUADRATIC_1_8 = (F(520), F(-3374), F(2917))  # 2917 x^2 - 3374 x + 520


def rank3_quadratic_root(iters=300):
    """Exact bisection root of the quadratic above (the low root, matching
    the isotropic m_top at f=(1/8,1/8,6/8))."""
    a, b, c = RANK3_QUADRATIC_1_8[2], RANK3_QUADRATIC_1_8[1], RANK3_QUADRATIC_1_8[0]

    def p(x):
        return a * x * x + b * x + c

    lo, hi = F(1, 8), F(1, 4)
    if (p(lo) > 0) == (p(hi) > 0):
        raise AssertionError("expected sign change")
    plo = p(lo)
    for _ in range(iters):
        mid = (lo + hi) / 2
        if (p(mid) > 0) == (plo > 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


# ---------------------------------------------------------------------------
# RANK 4
# ---------------------------------------------------------------------------

def rank4_tree(a: F, b: F, c: F, d: F) -> dict:
    X = {"fraction": str(b), "normal": [0, 1], "layers": [{"phase": "p3"}, {"phase": "p1"}]}
    Z = {"fraction": str(d), "normal": [1, 0], "layers": [{"phase": "p1"}, {"phase": "p3"}]}
    Y = {"fraction": str(c), "normal": [0, 1], "layers": [Z, {"phase": "p2"}]}
    return {"fraction": str(a), "normal": [1, 0], "layers": [X, Y]}


def rank4_c_of_a(a: F, f2: F) -> F:
    return 1 - f2 / (1 - a)


def rank4_b_of(a: F, d: F, c: F, f1: F) -> F:
    return 1 - (f1 - (1 - a) * c * d) / a


def rank4_tensor(a: F, d: F, f1: F, f2: F, sig=(S1, S2, S3)):
    s1, s2, s3 = sig
    c = rank4_c_of_a(a, f2)
    if not (0 < c < 1):
        return None
    b = rank4_b_of(a, d, c, f1)
    if not (0 < b < 1) or not (0 < d < 1):
        return None
    T11_X, T22_X = arith(b, s3, s1), harm(b, s3, s1)
    T11_Y = c * harm(d, s1, s3) + (1 - c) * s2
    T22_Y = 1 / (c / arith(d, s1, s3) + (1 - c) / s2)
    T11o = 1 / (a / T11_X + (1 - a) / T11_Y)
    T22o = a * T22_X + (1 - a) * T22_Y
    return T11o, T22o, b, c


def rank4_gap(a: F, d: F, f1: F, f2: F, sig=(S1, S2, S3)):
    t = rank4_tensor(a, d, f1, f2, sig)
    return None if t is None else t[0] - t[1]


def rank4_solve_d(a: F, f1: F, f2: F, sig=(S1, S2, S3), n=80, iters=70):
    lo, hi = F(1, 1000), F(999, 1000)
    prev = None
    prev_d = None
    bracket = None
    for k in range(n + 1):
        d = lo + (hi - lo) * F(k, n)
        g = rank4_gap(a, d, f1, f2, sig)
        if g is None:
            continue
        s = g > 0
        if prev is not None and s != prev:
            bracket = (prev_d, d)
            break
        prev, prev_d = s, d
    if bracket is None:
        return None
    plo, phi = bracket
    glo = rank4_gap(a, plo, f1, f2, sig)
    for _ in range(iters):
        mid = (plo + phi) / 2
        g = rank4_gap(a, mid, f1, f2, sig)
        if g is None:
            return None
        if (g > 0) == (glo > 0):
            plo = mid
        else:
            phi = mid
    return (plo + phi) / 2


def rank4_value_at(a: F, f1: F, f2: F, sig=(S1, S2, S3)):
    d = rank4_solve_d(a, f1, f2, sig)
    if d is None:
        return None, None
    t = rank4_tensor(a, d, f1, f2, sig)
    return t[0], d


def rank4_optimize(f1: F, f2: F, sig=(S1, S2, S3), lo=F(1, 100), hi=F(85, 100), iters=45, tol=F(1, 10**7)):
    """Golden-section search (float step sizes, exact Fraction evaluation at
    every point) for the minimum of sigma* over the free parameter a, along
    the isotropic curve d*(a)."""
    gr = (math.sqrt(5) - 1) / 2

    def val(a):
        v, _ = rank4_value_at(a, f1, f2, sig)
        return v if v is not None else F(10**9)

    c1 = hi - F(gr).limit_denominator(10**9) * (hi - lo)
    c2 = lo + F(gr).limit_denominator(10**9) * (hi - lo)
    f1v, f2v = val(c1), val(c2)
    for _ in range(iters):
        if f1v < f2v:
            hi, c2, f2v = c2, c1, f1v
            c1 = hi - F(gr).limit_denominator(10**9) * (hi - lo)
            f1v = val(c1)
        else:
            lo, c1, f1v = c1, c2, f2v
            c2 = lo + F(gr).limit_denominator(10**9) * (hi - lo)
            f2v = val(c2)
        if hi - lo < tol:
            break
    a_star = (lo + hi) / 2
    v_star, d_star = rank4_value_at(a_star, f1, f2, sig)
    return a_star, d_star, v_star


def hs_lo(f, sig=(S1, S2, S3)):
    lo, _ = L.hs_bounds(list(f), list(sig))
    return lo


# ---------------------------------------------------------------------------
# selftest: every closed form checked exactly against harness.effective()
# ---------------------------------------------------------------------------

def selftest() -> None:
    f1, f2, f3 = F(1, 8), F(1, 8), F(6, 8)
    phases = {"p1": S1, "p2": S2, "p3": S3}

    # rank-3 closed form vs harness, several m_top
    for m_top in (F(3, 20), F(1, 5), F(9, 40)):
        m_L, m_R = rank3_from_fractions(m_top, f2, f3)
        assert 0 < m_L < 1 and 0 < m_R < 1
        tree = rank3_tree(m_top, m_L, m_R)
        t = L.effective(tree, phases)
        assert t[1] == 0
        got = rank3_tensor(m_top, f2, f3)
        assert (t[0], t[2]) == got, (t, got, m_top)
        fr = L.fractions_of(tree)
        assert fr["p2"] == f2 and fr["p3"] == f3 and fr["p1"] == 1 - f2 - f3
        assert L.keller_check(tree, phases)

    # rank-3 isotropic point at f=(1/8,1/8,6/8): matches the screen (3.3848002...)
    m_star = rank3_isotropic_point(f2, f3)
    assert m_star is not None
    sig11, sig22 = rank3_tensor(m_star, f2, f3)
    assert abs(sig11 - sig22) < F(1, 10**20)
    val = float(sig11)
    assert abs(val - 3.3848001699) < 1e-6, val

    # the same point also solves the hand-derived quadratic
    # 2917 x^2 - 3374 x + 520 = 0 (isotropy numerator with the extraneous
    # factor of x removed -- see rank3_quadratic_root)
    a, b, c = F(2917), F(-3374), F(520)
    resid = a * m_star * m_star + b * m_star + c
    assert abs(resid) < F(1, 10**15), resid
    q_root = rank3_quadratic_root()
    assert abs(q_root - m_star) < F(1, 10**20)

    # discriminant is not a perfect square -> quadratic irreducible over Q
    disc = b * b - 4 * a * c
    assert disc == 5316516
    n = int(disc)
    r = math.isqrt(n)
    assert r * r != n, "discriminant must not be a perfect square"

    # rank-3: near team's rounded screen point, tree from harness matches too
    m_top0 = F(183106, 1000000)
    m_L0, m_R0 = F(317337, 1000000), F(81888, 1000000)
    tree0 = rank3_tree(m_top0, m_L0, m_R0)
    t0 = L.effective(tree0, phases)
    assert abs(float(t0[0]) - 3.3848017) < 1e-5
    assert abs(float(t0[2]) - 3.3847987) < 1e-5

    # rank-4 closed form vs harness: exact self-consistent point near the
    # screen (a, d free; b, c derived from the fraction constraints so the
    # tree hits f1=f2=1/8, f3=6/8 exactly, unlike the screen's rounded floats)
    a4, d4 = F(631443, 1000000), F(259276, 1000000)
    got4 = rank4_tensor(a4, d4, F(1, 8), f2)  # f1 = 1/8
    assert got4 is not None
    T11, T22, b4, c4 = got4
    tree4 = rank4_tree(a4, b4, c4, d4)
    t4 = L.effective(tree4, phases)
    assert t4[1] == 0
    assert (T11, T22) == (t4[0], t4[2]), (T11, T22, t4)
    fr4 = L.fractions_of(tree4)
    assert fr4["p1"] == f1 and fr4["p2"] == f2 and fr4["p3"] == f3
    assert abs(float(T11) - 3.364222) < 1e-4  # near the screen's 3.3642227
    assert L.keller_check(tree4, phases)

    # rank-4 optimum improves on (or matches) the screen and stays above HS_lo
    a_star, d_star, v_star = rank4_optimize(f1, f2)
    assert v_star is not None
    lo_bound = hs_lo((f1, f2, f3))
    assert v_star > lo_bound, "must not undercut the certified HS lower bound"
    assert v_star <= float(t4[0]) + F(1, 10**6), "exact optimum should not be worse than the screen"
    # verify the optimum against the harness at a nearby rational point
    c_star = rank4_c_of_a(a_star, f2)
    b_star = rank4_b_of(a_star, d_star, c_star, f1)

    def snap(x, den=10**9):
        return F(round(x * den), den)

    tree_star = rank4_tree(snap(a_star), snap(b_star), snap(c_star), snap(d_star))
    t_star = L.effective(tree_star, phases)
    assert abs(float(t_star[0]) - float(t_star[2])) < 1e-6
    assert abs(float(t_star[0]) - float(v_star)) < 1e-5

    print("selftest: all checks passed")


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def report() -> None:
    f1, f2, f3 = F(1, 8), F(1, 8), F(6, 8)
    lo_bound = hs_lo((f1, f2, f3))
    print(f"HS_lo at f=(1,1,6)/8, sigma=(1,2,5): {lo_bound} = {float(lo_bound):.10f}")

    m_star = rank3_isotropic_point(f2, f3)
    sig11, sig22 = rank3_tensor(m_star, f2, f3)
    print(f"\nRank-3 isotropic m_top* ~ {float(m_star):.12f}")
    print(f"  sigma* ~ {float(sig11):.10f}  (gap to HS_lo: {float(sig11 - lo_bound):.6e})")
    print("  minimal polynomial: 2917 x^2 - 3374 x + 520 = 0, discriminant 5316516 = 2^2 3^4 61 269 (not a square)")
    print(f"  x* = (1687 - 9*sqrt(16409))/2917 ~ {(F(1687) - 9 * F(math.sqrt(16409)).limit_denominator(10**12)) / F(2917)}")

    print("\nRank-3 leading-order-in-f1 sequence (f2:f3 fixed at 1:6):")
    for f1k in (F(1, 8), F(1, 16), F(1, 32)):
        f2k, f3k = (1 - f1k) / 7, 6 * (1 - f1k) / 7
        m = rank3_isotropic_point(f2k, f3k)
        if m is None:
            print(f"  f1={f1k}: NO isotropic point in this shape's family")
            continue
        s11, s22 = rank3_tensor(m, f2k, f3k)
        g = s11 - hs_lo((f1k, f2k, f3k))
        print(f"  f1={float(f1k):.6f}: sigma*={float(s11):.8f}  gap={float(g):.6e}  gap/f1={float(g/f1k):.6f}")
    for f1k in (F(1, 128), F(1, 512)):
        f2k, f3k = (1 - f1k) / 7, 6 * (1 - f1k) / 7
        m = rank3_isotropic_point(f2k, f3k)
        print(f"  f1={float(f1k):.6f}: {'isotropic point exists' if m else 'NO isotropic point (shape fails)'}")

    a_star, d_star, v_star = rank4_optimize(f1, f2)
    print(f"\nRank-4 optimum over the free parameter: a*~{float(a_star):.8f} d*~{float(d_star):.8f}")
    print(f"  sigma* ~ {float(v_star):.10f}  gap to HS_lo: {float(v_star - lo_bound):.6e}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--report", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
    if args.report:
        report()
    if not args.selftest and not args.report:
        ap.print_help()


if __name__ == "__main__":
    main()
