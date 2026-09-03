#!/usr/bin/env python3
"""Does a richer class beat the below-m11 family curve?  Answered locally, by
Schur complement, instead of by search.

015 lead 2 -- "push a richer class below the family's minimum" -- was attacked
in 016 by a two-parameter grid search over one altered topology.  The curvature
picture from 017 turns the same question into finite exact linear algebra.

Near m1 = m11 the family's gap is  gamma (m1 - m11)^2  with gamma the Schur
complement of the gap's Hessian at 004's attaining point (017).  Any class that
CONTAINS that point inherits the same base point, and enlarging the parameter
vector can only lower the Schur complement.  So "does a richer class beat the
curve?" becomes "does the Schur complement drop when the extra parameters are
switched on?" -- one Hessian, no search, and an answer that holds for every m1
just below m11 at once.

The enrichment used here is the largest one available at fixed topology: give
each of the five internal nodes a tilted normal,

    root (1, t0)   A (tA, 1)   D (1, tD)   C (tC, 1)   B (tB, 1),

all reducing to the axis-aligned base at t = 0.  Tilting breaks the automatic
vanishing of sigma12, so isotropy is now TWO conditions; a5 and t0 are solved
from them as series (the 2x2 Jacobian is invertible at the base), leaving
a0, a3, tA, tD, tC, tB free.  a3 is redundant in the axis-aligned family (005)
but not once normals tilt, so it is carried as a genuine parameter.

Usage:
    python tp_enrich.py --selftest
    python tp_enrich.py --report
Standard library only.
"""

from __future__ import annotations

import argparse
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import tp_curvature as TC  # noqa: E402


# ---------------------------------------------------------------------------
# degree-2 jets in N variables over Q
# ---------------------------------------------------------------------------

class J2:
    """Truncated Taylor series to total degree 2 in N variables."""

    __slots__ = ("c", "n")

    def __init__(self, c, n):
        self.c = {k: v for k, v in c.items() if v != 0}
        self.n = n

    @staticmethod
    def const(x, n):
        return J2({(): Fr(x)}, n)

    @staticmethod
    def var(i, n):
        return J2({(i,): Fr(1)}, n)

    def _co(self, o):
        return o if isinstance(o, J2) else J2.const(o, self.n)

    def __add__(self, o):
        o = self._co(o)
        c = dict(self.c)
        for k, v in o.c.items():
            c[k] = c.get(k, Fr(0)) + v
        return J2(c, self.n)

    __radd__ = __add__

    def __neg__(self):
        return J2({k: -v for k, v in self.c.items()}, self.n)

    def __sub__(self, o):
        return self + (-self._co(o))

    def __rsub__(self, o):
        return self._co(o) + (-self)

    def __mul__(self, o):
        o = self._co(o)
        c = {}
        for k1, a in self.c.items():
            for k2, b in o.c.items():
                if len(k1) + len(k2) > 2:
                    continue
                k = tuple(sorted(k1 + k2))
                c[k] = c.get(k, Fr(0)) + a * b
        return J2(c, self.n)

    __rmul__ = __mul__

    def __pow__(self, k):
        out = J2.const(1, self.n)
        for _ in range(k):
            out = out * self
        return out

    def inv(self):
        a0 = self.c.get((), Fr(0))
        if a0 == 0:
            raise ZeroDivisionError("jet with zero constant term")
        e = (self * (1 / a0)) - 1
        return (1 - e + e * e) * (1 / a0)

    def __truediv__(self, o):
        return self * self._co(o).inv()

    def __rtruediv__(self, o):
        return self._co(o) * self.inv()

    def const_term(self):
        return self.c.get((), Fr(0))

    def grad(self, i):
        return self.c.get((i,), Fr(0))

    def hess(self, i, j):
        """The Hessian entry d^2/dx_i dx_j (not the monomial coefficient)."""
        k = tuple(sorted((i, j)))
        return self.c.get(k, Fr(0)) * (2 if i == j else 1)


# ---------------------------------------------------------------------------
# exact rational linear algebra
# ---------------------------------------------------------------------------

def solve(M, b):
    """Solve M z = b exactly.  Returns (z, rank_deficient_directions).

    When M is singular the system is solved on the range and the free
    directions are reported; the caller checks b lies in the range.
    """
    n = len(M)
    A = [row[:] + [b[i]] for i, row in enumerate(M)]
    piv, r = [], 0
    for c in range(n):
        pr = None
        for i in range(r, n):
            if A[i][c] != 0:
                pr = i
                break
        if pr is None:
            continue
        A[r], A[pr] = A[pr], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c)
        r += 1
    for i in range(r, n):
        if A[i][n] != 0:
            raise ValueError("inconsistent system: b is outside the range of M")
    z = [Fr(0)] * n
    for i, c in enumerate(piv):
        z[c] = A[i][n]
    return z, [c for c in range(n) if c not in piv]


def solve2(J, rhs):
    """Exact 2x2 solve."""
    (a, b), (c, d) = J
    det = a * d - b * c
    if det == 0:
        raise RuntimeError("isotropy Jacobian is singular at the base point")
    r0, r1 = rhs
    return ((d * r0 - b * r1) / det, (a * r1 - c * r0) / det)


# ---------------------------------------------------------------------------
# the tilted family
# ---------------------------------------------------------------------------

def tilted_tensor(sig, m2, m1, a0, a3, a5, t0, tA, tD, tC, tB):
    """The 004 topology with every internal normal tilted; t = 0 is the base."""
    s1, s2, s3 = sig
    z = s1 * 0
    p1, p2, p3 = (s1, z, s1), (s2, z, s2), (s3, z, s3)
    Q = (1 - a0) - m2
    a4 = (m2 / (1 - a0) - (1 - a3)) / a3
    a1 = (m1 - Q * (1 - a5)) / a0
    one = (s1 * 0) + 1
    A = TC.lam(p1, p3, a1, (tA, one))
    D = TC.lam(p3, p1, a5, (one, tD))
    C = TC.lam(p2, D, a4, (tC, one))
    B = TC.lam(C, p2, a3, (tB, one))
    return TC.lam(A, B, a0, (one, t0)), a1, a4


# free directions, in order; index 0 is always u = m1 - m11
NAMES = ["u", "a0", "a3", "tA", "tD", "tC", "tB"]


def gap_jet(sig, s, active):
    """The gap g = value - HS_lo as a degree-2 jet in u plus the `active`
    parameters, with a5 and t0 eliminated by the two isotropy conditions."""
    bp = TC.base_point(sig, s)
    m2, m11, a0s, a5s = bp["m2"], bp["m11"], bp["a0"], bp["a5"]
    a3s = (1 + (1 - bp["r"])) / 2
    idx = {nm: i for i, nm in enumerate(["u"] + list(active))}
    n = len(idx)

    def build(vals, ring_n):
        """vals maps a name to a scalar/jet; missing names take base values."""
        g = lambda nm, d: vals.get(nm, d)      # noqa: E731
        return tilted_tensor(
            tuple(vals["sig"]), vals["m2"], vals["m1"], g("a0", vals["a0b"]),
            g("a3", vals["a3b"]), vals["a5"], vals["t0"],
            g("tA", vals["zero"]), g("tD", vals["zero"]),
            g("tC", vals["zero"]), g("tB", vals["zero"]))

    # --- the 2x2 isotropy Jacobian at the base, in (a5, t0)
    def base_vals(ring_n, extra):
        v = {"sig": [J2.const(x, ring_n) for x in sig],
             "m2": J2.const(m2, ring_n), "m1": J2.const(m11, ring_n),
             "a0b": J2.const(a0s, ring_n), "a3b": J2.const(a3s, ring_n),
             "a5": J2.const(a5s, ring_n), "t0": J2.const(0, ring_n),
             "zero": J2.const(0, ring_n)}
        v.update(extra)
        return v

    v = base_vals(2, {})
    v["a5"] = J2.const(a5s, 2) + J2.var(0, 2)
    v["t0"] = J2.const(0, 2) + J2.var(1, 2)
    t, _, _ = build(v, 2)
    Jac = [[(t[0] - t[2]).grad(0), (t[0] - t[2]).grad(1)],
           [t[1].grad(0), t[1].grad(1)]]

    # --- chord-Newton for (a5, t0) as series in the active variables
    W5 = J2.const(0, n)
    W0 = J2.const(0, n)
    for _ in range(5):
        v = base_vals(n, {})
        v["m1"] = J2.const(m11, n) + J2.var(idx["u"], n)
        for nm in active:
            base = {"a0": a0s, "a3": a3s}.get(nm, Fr(0))
            v[nm] = J2.const(base, n) + J2.var(idx[nm], n)
        v["a5"] = J2.const(a5s, n) + W5
        v["t0"] = W0
        t, a1, a4 = build(v, n)
        d5, d0 = solve2(Jac, [t[0] - t[2], t[1]])
        W5, W0 = W5 - d5, W0 - d0

    v = base_vals(n, {})
    v["m1"] = J2.const(m11, n) + J2.var(idx["u"], n)
    for nm in active:
        base = {"a0": a0s, "a3": a3s}.get(nm, Fr(0))
        v[nm] = J2.const(base, n) + J2.var(idx[nm], n)
    v["a5"] = J2.const(a5s, n) + W5
    v["t0"] = W0
    t, a1, a4 = build(v, n)

    value = (t[0] + t[2]) / 2
    hs = TC.hs_lo_gen([J2.const(x, n) for x in sig], v["m1"], J2.const(m2, n))
    return (value - hs), t, idx, bp


def schur(sig, s, active):
    """gamma over the class with `active` free parameters switched on."""
    g, t, idx, bp = gap_jet(sig, s, active)
    n = len(idx)
    assert (t[0] - t[2]).c == {} and t[1].c == {}, "isotropy series unsolved"
    assert g.const_term() == 0, g.const_term()
    for i in range(n):
        assert g.grad(i) == 0, ("gradient does not vanish", NAMES[i], g.grad(i))
    H = [[g.hess(i, j) for j in range(n)] for i in range(n)]
    Huu = H[0][0]
    if n == 1:
        return Huu / 2, [], bp
    M = [[H[i][j] for j in range(1, n)] for i in range(1, n)]
    b = [H[0][j] for j in range(1, n)]
    z, free = solve(M, b)
    gamma = (Huu - sum(bi * zi for bi, zi in zip(b, z))) / 2
    return gamma, free, bp


TRIPLES = [(1, 2, 5), (1, 3, 7), (1, 4, 9), (1, 2, 9), (2, 5, 11)]
ALL = ["a0", "a3", "tA", "tD", "tC", "tB"]


def _det(M):
    """Exact determinant by fraction-free elimination."""
    A = [row[:] for row in M]
    n = len(A)
    det = Fr(1)
    for c in range(n):
        pr = next((i for i in range(c, n) if A[i][c] != 0), None)
        if pr is None:
            return Fr(0)
        if pr != c:
            A[c], A[pr] = A[pr], A[c]
            det = -det
        det *= A[c][c]
        inv = 1 / A[c][c]
        A[c] = [x * inv for x in A[c]]
        for i in range(c + 1, n):
            if A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[c])]
    return det


def selftest():
    F = Fr
    for trip in TRIPLES:
        sig = tuple(F(x) for x in trip)
        s = F(1, 2)
        # 1. with only a0 active the tilted machinery must reproduce 017's gamma
        g1, free1, bp = schur(sig, s, ["a0"])
        ref = TC.closed_forms(sig, s)
        assert g1 == ref["gamma"], (trip, g1, ref["gamma"])
        # 2. a3 alone is inert in the axis-aligned family (005's collapse)
        g2, _, _ = schur(sig, s, ["a0", "a3"])
        assert g2 == g1, (trip, g2, g1)
        # 3. the tilts are REAL: a nonzero tilt breaks isotropy before the
        #    (a5, t0) re-solve, so they are not being silently ignored
        bp = TC.base_point(sig, s)
        a3b = (1 + (1 - bp["r"])) / 2
        t, _, _ = tilted_tensor(sig, bp["m2"], bp["m11"], bp["a0"], a3b,
                                bp["a5"], F(0), F(1, 10), F(0), F(0), F(0))
        assert t[1] != 0 and t[0] != t[2], "a tilt did nothing"

        # 4. the u-row of the Hessian vanishes on EVERY tilt direction, which
        #    is why no tilt can move gamma
        g, tt, idx, _ = gap_jet(sig, s, ALL)
        for nm in ("tA", "tD", "tC", "tB", "a3"):
            assert g.hess(0, idx[nm]) == 0, (trip, nm, g.hess(0, idx[nm]))
        assert g.hess(0, idx["a0"]) != 0, "a0 must couple to m1"

        # 5. the tilt block is positive semidefinite: tilting is uphill
        ti = [idx[nm] for nm in ("tA", "tD", "tC", "tB")]
        for k in range(1, len(ti) + 1):
            sub = [[g.hess(i, j) for j in ti[:k]] for i in ti[:k]]
            assert _det(sub) >= 0, (trip, "tilt block not PSD", k)

        # 6. the full tilt class: gamma is UNCHANGED, and c stays inside (0, 1)
        gf, freef, _ = schur(sig, s, ALL)
        assert gf == g1, (trip, gf, g1)
        assert 0 < gf / ref["beta"] < 1, gf / ref["beta"]
    print("selftest: tilts are real (they break isotropy before the re-solve) "
          "but the gap Hessian's m1-row vanishes on all five of them, the tilt "
          "block is PSD, and the six-parameter Schur complement is EQUAL to the "
          "one-parameter gamma on all 5 triples")


def report():
    F = Fr
    print("gamma over nested classes, m2 = 1/4.  c = gamma/beta with beta from 017.")
    for trip in TRIPLES:
        sig = tuple(F(x) for x in trip)
        s = F(1, 2)
        beta = TC.closed_forms(sig, s)["beta"]
        print(f"\nsigma = {trip}   beta = {float(beta):.8f}")
        rows = [["a0"], ["a0", "a3"], ["a0", "tA"], ["a0", "tD"], ["a0", "tC"],
                ["a0", "tB"], ["a0", "a3", "tB"], ALL]
        base = None
        for act in rows:
            gm, free, _ = schur(sig, s, act)
            if base is None:
                base = gm
            tag = "+".join(act)
            print(f"   {tag:<28} gamma = {float(gm):14.8f}  c = "
                  f"{float(gm / beta):.8f}  drop = {float(1 - gm / base):.3e}"
                  + (f"  [flat dirs {free}]" if free else ""))


# ---------------------------------------------------------------------------
# 016 re-checked
# ---------------------------------------------------------------------------

def _b_normal_tensor(sig, m2, m1, a0, a3, a5, nB):
    """016's altered topology: node B laminates C with p2 along nB."""
    s1, s2, s3 = sig
    z = s1 * 0
    p1, p2, p3 = (s1, z, s1), (s2, z, s2), (s3, z, s3)
    Q = (1 - a0) - m2
    a4 = (m2 / (1 - a0) - (1 - a3)) / a3
    a1 = (m1 - Q * (1 - a5)) / a0
    A = TC.lam(p1, p3, a1, TC.E2)
    D = TC.lam(p3, p1, a5, TC.E1)
    C = TC.lam(p2, D, a4, TC.E2)
    B = TC.lam(C, p2, a3, nB)
    return TC.lam(A, B, a0, TC.E1), a1, a4


def _iso_value_exact(sig, m2, m1, a0, a3, nB, iters=90):
    """Bisect a5 for isotropy in exact Q and return the effective value."""
    Q = (1 - a0) - m2
    if Q <= 0 or not (0 < a0 < 1 and 0 < a3 < 1):
        return None
    a4 = (m2 / (1 - a0) - (1 - a3)) / a3
    if not 0 < a4 < 1:
        return None
    lo = max(Fr(0), 1 - m1 / Q) + Fr(1, 10 ** 12)
    hi = 1 - Fr(1, 10 ** 12)

    def aniso(x):
        try:
            t, _, _ = _b_normal_tensor(sig, m2, m1, a0, a3, x, nB)
        except ZeroDivisionError:
            return None
        return t[0] - t[2]

    glo, ghi = aniso(lo), aniso(hi)
    if glo is None or ghi is None or (glo > 0) == (ghi > 0):
        return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        g = aniso(mid)
        if g is None:
            return None
        if (g > 0) == (glo > 0):
            lo = mid
        else:
            hi = mid
    a5 = (lo + hi) / 2
    t, a1, a4 = _b_normal_tensor(sig, m2, m1, a0, a3, a5, nB)
    if not (0 < a1 < 1 and 0 < a5 < 1):
        return None
    return (t[0] + t[2]) / 2


def recheck_016():
    """016 reported the e1-at-B class landing far ABOVE B2 (ratio 66 and 6.7).
    Its infimum is not there: it decreases monotonically as a3 -> 1, where the
    p2 layer at B has vanishing weight and B's normal stops mattering, and it
    converges to the ORIGINAL family's minimum."""
    F = Fr
    sig = (F(1), F(2), F(5))
    m2 = F(1, 4)
    print("sigma = (1,2,5), m2 = 1/4; node B relaid along e1 (016's change).")
    print("Exact arithmetic; a5 by exact isotropy bisection; a0 on a 1e-4 grid.")
    for m1, fam in ((F(3, 25), "3.0270075702"), (F(11, 100), "3.0831690087")):
        hs = TC.hs_lo_gen(sig, m1, m2)
        b2 = TC.B2_gen(sig, m1, F(1, 2))
        print()
        print(f" m1 = {m1}   HS_lo = {float(hs):.10f}   B2 = {float(b2):.10f}"
              f"   original family min ~ {fam}")
        for a3 in (F(1, 2), F(3, 5), F(9, 10), F(99, 100), F(999, 1000),
                   F(10 ** 5 - 1, 10 ** 5)):
            best, lo, hi, step = None, F(35, 100), F(74, 100), F(1, 500)
            for _ in range(5):                    # coarse grid, then refine
                a0 = lo
                while a0 <= hi:
                    v = _iso_value_exact(sig, m2, m1, a0, a3, TC.E1)
                    if v is not None and (best is None or v < best[0]):
                        best = (v, a0)
                    a0 += step
                if best is None:
                    break
                lo, hi = best[1] - step, best[1] + step
                step /= 10
            if best is None:
                print(f"   a3 = {a3}: no admissible member")
                continue
            v, a0 = best
            print(f"   a3 = {str(a3):>13}  min over a0 = {float(v):.10f}"
                  f"   ratio = {float((v - hs) / (b2 - hs)):9.4f}   a0 = {float(a0):.4f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--recheck-016", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    if a.report:
        report()
    if a.recheck_016:
        recheck_016()
