#!/usr/bin/env python3
"""Exact second-order expansion of the below-m11 gap, and a closed form for
015's constant c(sigma, m2).

015 measured, for the 006 family at m2 = 1/4,

    ratio(m1) = (family min - HS_lo) / (B2 - HS_lo)   ->   c(sigma)

as m1 -> m11 from below, with c between 0.53 and 0.69 on five triples, and left
"find a closed form for c, or refute that one exists" as its lead 1.

The observation that makes it exact.  At m1 = m11 all three quantities agree:
004's structure attains HS_lo there, and the transcribed B2 meets HS_lo there
(tp_cherkaev_bound selftest A).  Both differences therefore VANISH at m11, and
both vanish to SECOND order.  So the limit is a ratio of curvatures

    c  =  gamma / beta,
    beta  = the u^2 coefficient of B2 - HS_lo at m1 = m11 + u,
    gamma = the u^2 coefficient of the family's minimum minus HS_lo,

and gamma is a Schur complement.  Writing g(m1, a0) = value - HS_lo, which is
>= 0 everywhere and 0 at the base point (m11, 1-r), the gradient vanishes there,
so with g = A20 u^2 + A11 u v + A02 v^2 + O(3), u = m1 - m11, v = a0 - (1-r),

    min over v  =  (A20 - A11^2/(4 A02)) u^2,     gamma = A20 - A11^2/(4 A02).

Everything on the right is a rational number when sigma and r = sqrt(m2) are
rational, so c is EXACTLY rational -- computed here by truncated Taylor ("jet")
arithmetic over Fraction, with a5 obtained as a series by Newton on the isotropy
condition.

The laminate algebra is re-implemented here over the jet ring rather than
reusing harness/laminate.py (which coerces its inputs with Fraction()); the
constant terms are cross-checked against the harness on every run, so this
doubles as a second independent implementation of the family's tensor.

Usage:
    python tp_curvature.py --selftest
    python tp_curvature.py --table
    python tp_curvature.py --m2-scan
Standard library only.
"""

from __future__ import annotations

import argparse
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "three-phase-conductivity"))
import laminate as L  # noqa: E402

sys.path.insert(0, HERE)
import tp_cherkaev_bound as CB  # noqa: E402


# ---------------------------------------------------------------------------
# jets: truncated multivariate Taylor series in 2 variables over Q
# ---------------------------------------------------------------------------

class Jet:
    """Truncated Taylor series sum c_{ij} u^i v^j, i + j <= deg."""

    __slots__ = ("c", "deg")

    def __init__(self, c, deg):
        self.c = {k: v for k, v in c.items() if v != 0}
        self.deg = deg

    @staticmethod
    def const(x, deg):
        return Jet({(0, 0): Fr(x)}, deg)

    @staticmethod
    def var(i, deg):
        return Jet({(1, 0) if i == 0 else (0, 1): Fr(1)}, deg)

    def _co(self, other):
        return other if isinstance(other, Jet) else Jet.const(other, self.deg)

    def __add__(self, o):
        o = self._co(o)
        c = dict(self.c)
        for k, v in o.c.items():
            c[k] = c.get(k, Fr(0)) + v
        return Jet(c, self.deg)

    __radd__ = __add__

    def __neg__(self):
        return Jet({k: -v for k, v in self.c.items()}, self.deg)

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) + (-self)

    def __mul__(self, o):
        o = self._co(o)
        c = {}
        for (i, j), a in self.c.items():
            for (p, q), b in o.c.items():
                if i + p + j + q > self.deg:
                    continue
                k = (i + p, j + q)
                c[k] = c.get(k, Fr(0)) + a * b
        return Jet(c, self.deg)

    __rmul__ = __mul__

    def __pow__(self, n):
        assert isinstance(n, int) and n >= 0
        out = Jet.const(1, self.deg)
        for _ in range(n):
            out = out * self
        return out

    def inv(self):
        a0 = self.c.get((0, 0), Fr(0))
        if a0 == 0:
            raise ZeroDivisionError("jet with zero constant term")
        e = (self * (1 / a0)) - 1          # nilpotent past deg
        out = Jet.const(1, self.deg)
        term = Jet.const(1, self.deg)
        for _ in range(self.deg):
            term = term * (-e)
            out = out + term
        return out * (1 / a0)

    def __truediv__(self, o):
        return self * self._co(o).inv()

    def __rtruediv__(self, o):
        return self._co(o) * self.inv()

    def coeff(self, i, j):
        return self.c.get((i, j), Fr(0))

    def __repr__(self):
        return "Jet(" + ", ".join(f"{k}:{v}" for k, v in sorted(self.c.items())) + ")"


# ---------------------------------------------------------------------------
# the family's effective tensor, over any commutative ring of scalars
# ---------------------------------------------------------------------------

E1 = (1, 0)
E2 = (0, 1)


def lam(A, B, m, n):
    """Layering A (fraction m) with B (1-m), unnormalized integer normal n.

    Independent re-implementation of the lamination formula: no assertions on m
    (a jet here), no Fraction() coercion of the inputs.
    """
    u, v = n
    d0, d1, d2 = A[0] - B[0], A[1] - B[1], A[2] - B[2]
    w0 = d0 * u + d1 * v
    w1 = d1 * u + d2 * v
    c0 = A[0] * (1 - m) + B[0] * m
    c1 = A[1] * (1 - m) + B[1] * m
    c2 = A[2] * (1 - m) + B[2] * m
    denom = c0 * (u * u) + c1 * (2 * u * v) + c2 * (v * v)
    k = m * (1 - m) / denom
    return (A[0] * m + B[0] * (1 - m) - k * w0 * w0,
            A[1] * m + B[1] * (1 - m) - k * w0 * w1,
            A[2] * m + B[2] * (1 - m) - k * w1 * w1)


def family_tensor(sig, m2, m1, a0, a3, a5):
    """The 004/006 topology.

        A    = lam(p1 at a1, p3; e2)      D    = lam(p3 at a5, p1; e1)
        C    = lam(p2 at a4, D;  e2)      B    = lam(C  at a3, p2; e2)
        root = lam(A  at a0, B;  e1)

    a4 comes from the phase-2 volume constraint and a1 from phase-1's.
    """
    s1, s2, s3 = sig
    z = s1 * 0
    p1, p2, p3 = (s1, z, s1), (s2, z, s2), (s3, z, s3)
    Q = (1 - a0) - m2
    a4 = (m2 / (1 - a0) - (1 - a3)) / a3
    a1 = (m1 - Q * (1 - a5)) / a0
    A = lam(p1, p3, a1, E2)
    D = lam(p3, p1, a5, E1)
    C = lam(p2, D, a4, E2)
    B = lam(C, p2, a3, E2)
    return lam(A, B, a0, E1), a1, a4


def hs_lo_gen(sig, m1, m2):
    s1, s2, s3 = sig
    m3 = 1 - m1 - m2
    return 1 / (m1 / (s1 + s1) + m2 / (s2 + s1) + m3 / (s3 + s1)) - s1


def B2_gen(sig, m1, s):
    """Cherkaev's B2 [T], written so it accepts jets (cf. tp_cherkaev_bound)."""
    k1, k2, k3 = sig
    m2 = s * s
    m3 = 1 - m1 - m2
    Z5 = m1 * k1 ** 2 - m1 * k2 ** 2 + 2 * m3 * k1 * (k3 - k2)
    Z6 = (((1 - s) ** 2 + (1 - m1 - s) ** 2) * k1
          + m1 * (1 - s) ** 2 * k2 + m1 * m3 * k3)
    return k2 + (1 - s) ** 2 * Z5 / Z6


# ---------------------------------------------------------------------------
# the expansion
# ---------------------------------------------------------------------------

def base_point(sig, s):
    """004's attaining parameters at m1 = m11 (s = sqrt(m2) = r)."""
    k1, k2, k3 = sig
    Th = k1 * (k3 - k2) / ((k2 + k1) * (k3 - k1))
    return {"theta": Th, "r": s, "m2": s * s, "m11": 2 * Th * s * (1 - s),
            "a0": 1 - s, "a1": s * Th, "a5": 1 - Th}


def expand(sig, s, deg=3, a3=None):
    """Jets g(u, v) = value - HS_lo and d(u) = B2 - HS_lo about the base point,
    plus beta, gamma and c = gamma/beta."""
    bp = base_point(sig, s)
    m2, m11, a0s, a5s = bp["m2"], bp["m11"], bp["a0"], bp["a5"]
    if a3 is None:
        a3 = (1 + (1 - bp["r"])) / 2          # any value in (1-r, 1)

    U, V = Jet.var(0, deg), Jet.var(1, deg)
    m1 = Jet.const(m11, deg) + U
    a0 = Jet.const(a0s, deg) + V
    sigj = tuple(Jet.const(x, deg) for x in sig)
    a3j, m2j = Jet.const(a3, deg), Jet.const(m2, deg)

    # d(sigma11 - sigma22)/d(a5) at the base point, exactly
    W = Jet.var(0, 1)
    t, _, _ = family_tensor(tuple(Jet.const(x, 1) for x in sig),
                            Jet.const(m2, 1), Jet.const(m11, 1),
                            Jet.const(a0s, 1), Jet.const(a3, 1),
                            Jet.const(a5s, 1) + W)
    A = (t[0] - t[2]).coeff(1, 0)
    if A == 0:
        raise RuntimeError("isotropy condition is degenerate in a5")

    # chord-Newton for a5 as a series in (u, v); each pass gains one degree
    Wk = Jet.const(0, deg)
    for _ in range(deg + 3):
        t, a1, a4 = family_tensor(sigj, m2j, m1, a0, a3j, Jet.const(a5s, deg) + Wk)
        Wk = Wk - (t[0] - t[2]) / A
    t, a1, a4 = family_tensor(sigj, m2j, m1, a0, a3j, Jet.const(a5s, deg) + Wk)

    value = (t[0] + t[2]) / 2
    hs = hs_lo_gen(sigj, m1, m2j)
    g = value - hs
    d = B2_gen(sigj, m1, Jet.const(s, deg)) - hs

    A20, A11, A02 = g.coeff(2, 0), g.coeff(1, 1), g.coeff(0, 2)
    beta = d.coeff(2, 0)
    gamma = A20 - A11 * A11 / (4 * A02) if A02 != 0 else None
    c = gamma / beta if (gamma is not None and beta != 0) else None
    return {"bp": bp, "g": g, "d": d, "a5_series": Wk, "a1": a1, "a4": a4,
            "iso_residual": t[0] - t[2], "off_diag": t[1],
            "A20": A20, "A11": A11, "A02": A02,
            "beta": beta, "gamma": gamma, "c": c,
            "vslope": (-A11 / (2 * A02)) if A02 != 0 else None}


# ---------------------------------------------------------------------------
# closed forms
# ---------------------------------------------------------------------------
#
# beta, gamma and c are rational functions of (x, y, r) with x = k2/k1,
# y = k3/k1 (the problem is homogeneous in sigma, so only the ratios enter),
# reconstructed by exact rational-function interpolation and then verified as
# identities at thousands of rational points by --closed-form.  Writing
#
#     p  = (1 + x)/(y - x),          q  = (y - 1)/(x - 1),
#     A  = p (y - 1),                n1 = (1 + x)(y - 1)/((x - 1)(y + 1)),
#     K  = k1 p^3 (y - 1)^2 (y + 1)/4,
#     q2(r) = C2 + B2r r - r^2,
#     C2 = (1 + x)^3 y (y - 1) / (2 (y - x)^2 (x y - 1)),
#     B2r = (1 + x)[2x(y^3 - 1) - (x^2 + 2x + 3) y^2 + (3x^2 + 2x + 1) y]
#           / (2 (y - x)^2 (x y - 1)),
#
# the three quantities are
#
#     beta  = K (r + n1) / [(1 - r)(r + p)^2 (r + p q)],
#     gamma = K r (A - r) / [(1 - r)(r + p)^2 q2(r)],
#     c     = r (A - r)(r + p q) / [(r + n1) q2(r)],
#
# the K and the (1 - r)(r + p)^2 cancelling between the two curvatures.


def closed_forms(sig, r):
    """beta, gamma, c from the closed forms above (no jets)."""
    k1, k2, k3 = sig
    x, y = k2 / k1, k3 / k1
    p = (1 + x) / (y - x)
    q = (y - 1) / (x - 1)
    A = p * (y - 1)
    n1 = (1 + x) * (y - 1) / ((x - 1) * (y + 1))
    K = k1 * p ** 3 * (y - 1) ** 2 * (y + 1) / 4
    den = 2 * (y - x) ** 2 * (x * y - 1)
    C2 = (1 + x) ** 3 * y * (y - 1) / den
    B2r = ((1 + x) * (2 * x * (y ** 3 - 1) - (x * x + 2 * x + 3) * y * y
                      + (3 * x * x + 2 * x + 1) * y) / den)
    q2 = C2 + B2r * r - r * r
    beta = K * (r + n1) / ((1 - r) * (r + p) ** 2 * (r + p * q))
    gamma = K * r * (A - r) / ((1 - r) * (r + p) ** 2 * q2)
    c = r * (A - r) * (r + p * q) / ((r + n1) * q2)
    return {"beta": beta, "gamma": gamma, "c": c, "p": p, "q": q, "A": A,
            "n1": n1, "K": K, "q2": q2}


TRIPLES = [(1, 2, 5), (1, 3, 7), (1, 4, 9), (1, 2, 9), (2, 5, 11)]

# 015's measured limiting ratios, for the cross-check
MEASURED = {(1, 2, 5): 0.532369, (2, 5, 11): 0.532696, (1, 3, 7): 0.576304,
            (1, 4, 9): 0.598509, (1, 2, 9): 0.690125}


def _harness_tree(bp, a1, a3, a4):
    return {"fraction": str(bp["a0"]), "normal": [1, 0], "layers": [
        {"fraction": str(a1), "normal": [0, 1],
         "layers": [{"phase": "p1"}, {"phase": "p3"}]},
        {"fraction": str(a3), "normal": [0, 1], "layers": [
            {"fraction": str(a4), "normal": [0, 1], "layers": [
                {"phase": "p2"},
                {"fraction": str(bp["a5"]), "normal": [1, 0],
                 "layers": [{"phase": "p3"}, {"phase": "p1"}]}]},
            {"phase": "p2"}]}]}


def selftest():
    F = Fr
    for trip in TRIPLES:
        sig = tuple(F(x) for x in trip)
        s = F(1, 2)
        R = expand(sig, s, deg=3)
        bp = R["bp"]

        # 0. the jet ring's constant term must reproduce the harness exactly
        a3 = (1 + (1 - bp["r"])) / 2
        t0, a1c, a4c = family_tensor(sig, bp["m2"], bp["m11"], bp["a0"], a3, bp["a5"])
        assert a1c == bp["a1"], (a1c, bp["a1"])
        tree = _harness_tree(bp, a1c, a3, a4c)
        ph = {"p1": sig[0], "p2": sig[1], "p3": sig[2]}
        assert L.effective(tree, ph) == t0, "jet ring disagrees with the harness"
        assert L.keller_check(tree, ph), "Keller-Dykhne fails at the base point"
        fr = L.fractions_of(tree)
        assert fr["p1"] == bp["m11"] and fr["p2"] == bp["m2"], fr
        assert t0[0] == t0[2] and t0[1] == 0, "base point must be isotropic"
        assert t0[0] == CB.hs_lo(sig, bp["m11"], bp["m2"]), "004 must attain"

        # 1. the isotropy series is solved to the truncation order, and the
        #    structure stays isotropic off the base point
        assert R["iso_residual"].c == {}, R["iso_residual"]
        assert R["off_diag"].c == {}, R["off_diag"]

        # 2. gap and gradient vanish at the base point, EXACTLY in Q
        g, d = R["g"], R["d"]
        assert g.coeff(0, 0) == 0 and g.coeff(1, 0) == 0 and g.coeff(0, 1) == 0, g
        assert d.coeff(0, 0) == 0 and d.coeff(1, 0) == 0, d

        # 3. the quadratic form is PSD -- it must be, HS is a valid lower bound
        assert R["A20"] > 0 and R["A02"] > 0, (R["A20"], R["A02"])
        assert 4 * R["A20"] * R["A02"] >= R["A11"] ** 2, "gap Hessian not PSD"
        assert R["gamma"] > 0 and R["beta"] > 0

        # 4. a3 really is redundant: a different a3 gives the same jets
        R2 = expand(sig, s, deg=3, a3=(1 + 3 * (1 - bp["r"])) / 4)
        assert R2["g"].c == g.c, "a3 is not redundant in the expansion"
        assert R2["c"] == R["c"]

        # 5. the closed form reproduces 015's measured limit.  015's values are
        #    consistently ~4e-5 LOW: they were read off at a finite m1 rather
        #    than extrapolated, and the ratio increases as m1 -> m11.  The
        #    direct check that pins the limit is --verify, where the ratio at
        #    h = 1e-5 agrees with c to 6-8 digits.
        assert 0 <= float(R["c"]) - MEASURED[trip] < 1e-4, (trip, float(R["c"]))

        # 6. 0 < c < 1: the family's curvature is strictly under B2's
        assert 0 < R["c"] < 1, R["c"]

        # 7. the closed forms reproduce the jet computation exactly
        for r2 in (F(1, 3), F(1, 2), F(3, 4)):
            J = expand(sig, r2, deg=2)
            C = closed_forms(sig, r2)
            assert (J["beta"], J["gamma"], J["c"]) == (C["beta"], C["gamma"], C["c"]),                 (trip, r2)
    print("selftest: base point attains HS in Q on 5 triples; gap and gradient "
          "vanish exactly; Hessian PSD; a3 redundant; 0 < c < 1; and "
          "gamma/beta sits 0..1e-4 ABOVE 015's measured c everywhere")


def table():
    print(f"{'sigma':>12} {'m11':>10} {'beta':>15} {'gamma':>15} "
          f"{'c':>11} {'015 measured':>13}")
    for trip in TRIPLES:
        sig = tuple(Fr(x) for x in trip)
        R = expand(sig, Fr(1, 2), deg=2)
        print(f"{str(trip):>12} {str(R['bp']['m11']):>10} {float(R['beta']):15.8f} "
              f"{float(R['gamma']):15.8f} {float(R['c']):11.8f} "
              f"{MEASURED[trip]:13.6f}")
        print(f"{'':>12} beta = {R['beta']}")
        print(f"{'':>12} gamma = {R['gamma']}")
        print(f"{'':>12} c = {R['c']}")


def _iso_a5_exact(sig, m2, m1, a0, a3, lo, hi, iters=180):
    """Bisect a5 for sigma11 == sigma22, exactly in Q (no jets involved)."""
    def aniso(x):
        t, _, _ = family_tensor(sig, m2, m1, a0, a3, x)
        return t[0] - t[2]
    glo = aniso(lo)
    if (glo > 0) == (aniso(hi) > 0):
        return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        if (aniso(mid) > 0) == (glo > 0):
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def _family_value_float(sig, m2, m1, a0, a3):
    """Same family, in floats, with a5 found by float bisection.  The bracket
    is a5 > 1 - m1/Q, which is what keeps a1 positive."""
    def aniso(x):
        t, _, _ = family_tensor(sig, m2, m1, a0, a3, x)
        return t[0] - t[2]
    Q = (1 - a0) - m2
    if Q <= 0:
        return None
    lo, hi = max(0.0, 1 - m1 / Q) + 1e-12, 1 - 1e-12
    glo = aniso(lo)
    if (glo > 0) == (aniso(hi) > 0):
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (aniso(mid) > 0) == (glo > 0):
            lo = mid
        else:
            hi = mid
    a5 = 0.5 * (lo + hi)
    t, a1, a4 = family_tensor(sig, m2, m1, a0, a3, a5)
    if not (0 < a1 < 1 and 0 < a4 < 1 and 0 < a5 < 1):
        return None
    return 0.5 * (t[0] + t[2])


def verify(trips=None, hs=(2, 3, 4, 5)):
    """Independent numeric confirmation that gap_min(h) / h^2 -> gamma and
    ratio(h) -> c, using the harness algebra and a direct search over a0 rather
    than the jet expansion."""
    F = Fr
    for trip in (trips or TRIPLES):
        sig = tuple(F(x) for x in trip)
        s = F(1, 2)
        R = expand(sig, s, deg=2)
        bp = R["bp"]
        a3 = (1 + (1 - bp["r"])) / 2
        gamma, beta, c = R["gamma"], R["beta"], R["c"]
        print(f"sigma = {trip}   gamma = {float(gamma):.10f}   "
              f"beta = {float(beta):.10f}   c = {c} = {float(c):.10f}")

        sf = tuple(float(x) for x in sig)
        m2f, m11f, a0f, a3f = (float(bp["m2"]), float(bp["m11"]),
                               float(bp["a0"]), float(a3))
        for e in hs:
            h = 10.0 ** (-e)
            m1f = m11f - h
            # ternary search over a0, seeded at the jet-predicted optimum
            centre = a0f + float(R["vslope"]) * (-h)
            lo, hi = centre - 50 * h, centre + 50 * h
            for _ in range(200):
                x1 = lo + (hi - lo) / 3
                x2 = hi - (hi - lo) / 3
                v1 = _family_value_float(sf, m2f, m1f, x1, a3f)
                v2 = _family_value_float(sf, m2f, m1f, x2, a3f)
                if v1 is None or v2 is None:
                    break
                if v1 < v2:
                    hi = x2
                else:
                    lo = x1
            a0best = 0.5 * (lo + hi)
            val = _family_value_float(sf, m2f, m1f, a0best, a3f)
            hsv = float(hs_lo_gen(sig, F(m11f) - F(1, 10 ** e), bp["m2"]))
            b2v = float(B2_gen(sig, F(m11f) - F(1, 10 ** e), s))
            gap = val - hsv
            print(f"   h=1e-{e}: gap/h^2 = {gap / h / h:14.8f}  (gamma "
                  f"{float(gamma):14.8f})   ratio = {gap / (b2v - hsv):.8f}  "
                  f"(c {float(c):.8f})   a0-a0* = {a0best - a0f:+.3e}  "
                  f"(pred {float(R['vslope']) * -h:+.3e})")

        # one fully exact point on the predicted optimal path
        e = 3
        u = -F(1, 10 ** e)
        m1e = bp["m11"] + u
        a0e = bp["a0"] + R["vslope"] * u
        Qe = (1 - a0e) - bp["m2"]
        a5e = _iso_a5_exact(sig, bp["m2"], m1e, a0e, a3,
                            max(F(0), 1 - m1e / Qe) + F(1, 10 ** 12),
                            1 - F(1, 10 ** 12))
        t, a1e, a4e = family_tensor(sig, bp["m2"], m1e, a0e, a3, a5e)
        tree = _harness_tree({"a0": a0e, "a5": a5e}, a1e, a3, a4e)
        ph = {"p1": sig[0], "p2": sig[1], "p3": sig[2]}
        te = L.effective(tree, ph)
        assert te == t, "harness and local algebra disagree off the base point"
        fr = L.fractions_of(tree)
        assert fr["p1"] == m1e and fr["p2"] == bp["m2"], fr
        gape = (te[0] + te[2]) / 2 - hs_lo_gen(sig, m1e, bp["m2"])
        print(f"   exact at h=1e-{e} on the predicted path: gap/h^2 = "
              f"{float(gape / u / u):.8f}  (gamma {float(gamma):.8f}), "
              f"anisotropy {float(te[0] - te[2]):.2e}, verified by the harness")


def closed_form_check(nsig=None):
    """The closed forms against the jet computation, then a wide cheap sweep of
    0 < c < 1 using the closed form alone."""
    F = Fr
    tests = []
    for k1 in (1, 2, 3, 5):
        for k2 in range(k1 + 1, k1 + 9):
            for k3 in range(k2 + 1, k2 + 9):
                tests.append((F(k1), F(k2), F(k3)))
    tests += [(F(1), F(5, 2), F(11, 2)), (F(2), F(7, 3), F(9, 2)),
              (F(3, 2), F(4), F(17, 3)), (F(1, 3), F(4, 7), F(9, 5)),
              (F(7), F(23, 2), F(101, 4))]
    rs = [F(1, 7), F(1, 5), F(1, 3), F(2, 5), F(1, 2), F(3, 5), F(2, 3),
          F(3, 4), F(7, 9), F(9, 10)]
    n = bad = 0
    for S in tests[:nsig] if nsig else tests:
        for r in rs:
            R = expand(S, r, deg=2)
            C = closed_forms(S, r)
            n += 1
            if (R["beta"], R["gamma"], R["c"]) != (C["beta"], C["gamma"], C["c"]):
                bad += 1
                if bad < 4:
                    print("MISMATCH", S, r, R["c"], C["c"])
    print(f"closed forms vs the jet expansion: {n - bad}/{n} exact, {bad} bad")

    lo, hi, m = None, None, 0
    viol = []
    for k2n in range(3, 60):
        for k3n in range(k2n + 1, 90):
            x, y = F(k2n, 2), F(k3n, 2)
            if not (1 < x < y):
                continue
            for rn in range(1, 24):
                r = F(rn, 24)
                m2 = r * r
                th = (y - x) / ((x + 1) * (y - 1))
                m11 = 2 * th * r * (1 - r)
                if m11 <= 0 or 1 - m11 - m2 <= 0:
                    continue
                c = closed_forms((F(1), x, y), r)["c"]
                m += 1
                if not (0 < c < 1):
                    viol.append((x, y, r, c))
                if lo is None or c < lo[0]:
                    lo = (c, x, y, r)
                if hi is None or c > hi[0]:
                    hi = (c, x, y, r)
    print(f"closed-form sweep: {m} admissible (sigma, m2) points, "
          f"{len(viol)} with c outside (0, 1)")
    print(f"   min c = {float(lo[0]):.8f} at x={lo[1]} y={lo[2]} r={lo[3]}")
    print(f"   max c = {float(hi[0]):.8f} at x={hi[1]} y={hi[2]} r={hi[3]}")
    if viol:
        print("   VIOLATIONS:", viol[:5])


def m2_scan():
    """c at several rational r = sqrt(m2): does c depend on m2 as well?"""
    print(f"{'sigma':>12} {'r':>7} {'m2':>9} {'m11':>12} {'c':>12}   c exact")
    for trip in [(1, 2, 5), (1, 3, 7)]:
        sig = tuple(Fr(x) for x in trip)
        for r in [Fr(1, 4), Fr(1, 3), Fr(2, 5), Fr(1, 2), Fr(3, 5), Fr(2, 3), Fr(3, 4)]:
            R = expand(sig, r, deg=2)
            c = R["c"]
            print(f"{str(trip):>12} {str(r):>7} {str(r * r):>9} "
                  f"{str(R['bp']['m11']):>12} {float(c):12.8f}   "
                  f"{c.numerator}/{c.denominator}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--table", action="store_true")
    ap.add_argument("--m2-scan", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--closed-form", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    if a.table:
        table()
    if a.m2_scan:
        m2_scan()
    if a.verify:
        verify()
    if a.closed_form:
        closed_form_check()
