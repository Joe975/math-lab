#!/usr/bin/env python3
"""Cherkaev (2009) Section 8.1's HS-lower-bound-attaining laminate, built and
verified exactly.

Target: sigma = (1, 2, 5), f = (1/8, 1/8, 3/4).  HS_lo(f, sigma) = 37/11.
This point is BELOW Milton's classical threshold (2*Theta*(1-m2) = 0.4375,
see attempt 001) so the plain coated (Milton) assemblage cannot reach it; the
best laminate found by bounded-rank search in attempts 001/003 stalled
5.46e-4 above the bound at rank 4.

THE OBSTRUCTION THIS SCRIPT RESOLVES: the point m1 = 1/8 sits strictly above
Cherkaev's threshold m11 = 2*Theta*sqrt(m2)*(1-sqrt(m2)) = 0.114277 (with the
Theta of this sigma), so per Cherkaev 2009 Section 8.1 the attaining
structure is not the "plain" T^2-structure ("L13,2,13") -- that one attains
HS_lo only at the SINGLE point m1 = m11 -- but that structure with EXTRA
phase-1 coated on in two orthogonal equal-fraction layers ("L13,2,13,1,1"),
raising m1 from m11 up to the target 1/8.  Both the plain-T^2 attaining point
and the amount of extra coating needed are IRRATIONAL numbers: the T^2's own
m11-fixed isotropic member forces sqrt(m2_inner) to appear, and matching the
outer target ratio f2:f3 = 1:6 forces the T^2's OWN m2 (not the final m2!) to
solve a quadratic with an irrational root.  So the exact attaining tree has
irrational (algebraic, in Q(sqrt(105))) node fractions -- this is normal (see
001's rank-3 shape, whose isotropic point is irrational too) and does not
prevent an EXACT equality: the final effective tensor and the final volume
fractions come out exactly RATIONAL despite every intermediate quantity
living in Q(sqrt(105)).  That is proved below by exact field arithmetic (no
floating point anywhere in the proof), independently by TWO different
lamination algorithms (mirroring laminate.py's projection formula and
verify_laminate.py's rotate-and-average algebra), and cross-checked by a
rational sandwich (two purely-Fraction trees, obtained by snapping sqrt(105)
from below and above to a tiny width) that is fed through the REAL
harness/three-phase-conductivity/laminate.py and verify_laminate.py.

THE CLOSED FORM (new, found by direct algebraic solve rather than transcribed
from the paper -- see attempt record for the derivation and why the
transcription route (leads 1-2 of attempt 003) was abandoned): for this sigma
(Theta = 1/4 exactly) and ANY inner phase-2 fraction m2 in (0,1), the rank-4
axis-normal tree

    X = lam(p3, p1; frac b; normal e2)
    Z = lam(p1, p3; frac d; normal e1)
    Y = lam(Z,  p2; frac c; normal e2)
    T2 = lam(X, Y;  frac a; normal e1)          <- "L13,2,13" (T^2-structure)

with

    u = sqrt(m2),  a = c = 1 - u,  d = Theta = 1/4,  b = 1 - u/4

is isotropic and EXACTLY equal to HS_lo(m11(m2), m2, 1 - m11(m2) - m2) at
m1 = m11(m2) = 2*Theta*u*(1-u).  Verified symbolically over Q(sqrt(m2)) at
five different m2 (--selftest) and used here at the specific m2 that makes
the coated-up structure land exactly on the target fractions.

Standard library only (a small exact quadratic-field class, no sympy).

Usage:
    python tp_cherkaev.py --selftest
    python tp_cherkaev.py --build      builds+verifies the target structure,
                                        writes the record and enclosure
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "three-phase-conductivity"))
import laminate as L          # noqa: E402
import verify_laminate as V   # noqa: E402

DATA = os.path.join(HERE, "..", "data", "cherkaev")

E1, E2 = (1, 0), (0, 1)


# ---------------------------------------------------------------------------
# Exact quadratic field Q(sqrt(D)), elements p + q*sqrt(D), p, q in Fraction.
# Only +,-,*,/ and equality are needed; both are field operations so this is
# a genuine field (D is fixed at construction time and never mixed across
# instances with different D in the same computation -- guarded by assert).
# ---------------------------------------------------------------------------


class Quad:
    __slots__ = ("D", "p", "q")

    def __init__(self, D: int, p, q=0):
        self.D = D
        self.p = Fr(p)
        self.q = Fr(q)

    def _c(self, o):
        if isinstance(o, Quad):
            assert o.D == self.D, "field mismatch"
            return o
        return Quad(self.D, o, 0)

    def __add__(self, o):
        o = self._c(o)
        return Quad(self.D, self.p + o.p, self.q + o.q)

    __radd__ = __add__

    def __sub__(self, o):
        o = self._c(o)
        return Quad(self.D, self.p - o.p, self.q - o.q)

    def __rsub__(self, o):
        return self._c(o) - self

    def __neg__(self):
        return Quad(self.D, -self.p, -self.q)

    def __mul__(self, o):
        o = self._c(o)
        return Quad(self.D, self.p * o.p + self.D * self.q * o.q, self.p * o.q + self.q * o.p)

    __rmul__ = __mul__

    def __truediv__(self, o):
        o = self._c(o)
        den = o.p * o.p - self.D * o.q * o.q
        assert den != 0, "division by zero in Q(sqrt(D))"
        num = self * Quad(self.D, o.p, -o.q)
        return Quad(self.D, num.p / den, num.q / den)

    def __rtruediv__(self, o):
        return self._c(o) / self

    def __eq__(self, o):
        o = self._c(o)
        return self.p == o.p and self.q == o.q

    def __hash__(self):
        return hash((self.D, self.p, self.q))

    def __repr__(self):
        return f"({self.p}+{self.q}*sqrt{self.D})"

    def to_float(self) -> float:
        return float(self.p) + float(self.q) * (self.D ** 0.5)

    def is_rational(self) -> bool:
        return self.q == 0


def qi(D, x):
    """Lift a Fraction/int into Q(sqrt(D))."""
    return Quad(D, x, 0)


# ---------------------------------------------------------------------------
# Generic 2x2 symmetric-tensor lamination over ANY field supporting the
# overloaded operators -- two independent algorithms.
# ---------------------------------------------------------------------------


def laminate_proj(A, B, m, n):
    """laminate.py's frame-free projection formula, generic in the field."""
    u, v = n
    dA = (A[0] - B[0], A[1] - B[1], A[2] - B[2])
    w = (dA[0] * u + dA[1] * v, dA[1] * u + dA[2] * v)
    one_m = 1 - m
    comp = (A[0] * one_m + B[0] * m, A[1] * one_m + B[1] * m, A[2] * one_m + B[2] * m)
    denom = comp[0] * u * u + comp[1] * (2 * u * v) + comp[2] * v * v
    coeff = (m * one_m) / denom
    mean = (A[0] * m + B[0] * one_m, A[1] * m + B[1] * one_m, A[2] * m + B[2] * one_m)
    return (
        mean[0] - coeff * w[0] * w[0],
        mean[1] - coeff * w[0] * w[1],
        mean[2] - coeff * w[1] * w[1],
    )


def laminate_frame(A, B, m, n, D):
    """verify_laminate.py's rotate + Backus-average algebra, generic in the
    field.  u, v are plain Python ints (integer normal); division by
    u^2+v^2 is done after lifting it into the field."""
    u, v = n
    d_uv = qi(D, u * u + v * v)

    def conj(t, inverse):
        uu, vv = (u, -v) if inverse else (u, v)
        a, b, c = t
        m00 = uu * (a * uu + b * vv) + vv * (b * uu + c * vv)
        m01 = uu * (-a * vv + b * uu) + vv * (-b * vv + c * uu)
        m11 = -vv * (-a * vv + b * uu) + uu * (-b * vv + c * uu)
        return (m00 / d_uv, m01 / d_uv, m11 / d_uv)

    Af, Bf = conj(A, False), conj(B, False)
    wA, wB = m, 1 - m
    inv_a = wA / Af[0] + wB / Bf[0]
    s11 = 1 / inv_a
    b_over_a = wA * Af[1] / Af[0] + wB * Bf[1] / Bf[0]
    schur = wA * (Af[2] - Af[1] * Af[1] / Af[0]) + wB * (Bf[2] - Bf[1] * Bf[1] / Bf[0])
    s12 = s11 * b_over_a
    s22 = schur + s11 * b_over_a * b_over_a
    return conj((s11, s12, s22), True)


def iso_leaf(D, s):
    return (qi(D, s), qi(D, 0), qi(D, s))


# ---------------------------------------------------------------------------
# The T^2 ("L13,2,13") family, symbolic in u = sqrt(m2), sigma = (1,2,5).
# ---------------------------------------------------------------------------


SIGMA = (Fr(1), Fr(2), Fr(5))
THETA = Fr(1, 4)  # s1(s3-s2) / ((s2+s1)(s3-s1)) = 1*3/(3*4)


def t2_structure(D, u):
    """u = sqrt(m2) in Q(sqrt(D)).  Returns (a, b, c, d) node fractions."""
    a = 1 - u
    c = 1 - u
    d = qi(D, THETA)
    b = 1 - u / 4
    return a, b, c, d


def t2_tensor_proj(D, u, algo="proj"):
    s1, s2, s3 = (qi(D, s) for s in SIGMA)
    a, b, c, d = t2_structure(D, u)
    lam = laminate_proj if algo == "proj" else (lambda A, B, m, n: laminate_frame(A, B, m, n, D))
    X = lam(iso_leaf(D, SIGMA[2]), iso_leaf(D, SIGMA[0]), b, E2)
    Z = lam(iso_leaf(D, SIGMA[0]), iso_leaf(D, SIGMA[2]), d, E1)
    Y = lam(Z, iso_leaf(D, SIGMA[1]), c, E2)
    T2 = lam(X, Y, a, E1)
    return T2


def hs_lo_field(D, fracs, sigs):
    s1 = min(sigs, key=lambda x: x.to_float() if isinstance(x, Quad) else x)
    total = sum((f / (s + s1) for f, s in zip(fracs, sigs)), qi(D, 0))
    return 1 / total - s1


# ---------------------------------------------------------------------------
# selftest: the T^2 closed form attains HS_lo exactly at m1 = m11(m2), for
# several m2, by BOTH lamination algorithms.
# ---------------------------------------------------------------------------


def selftest() -> None:
    D = 2  # sqrt(m2) for these test points lives in Q(sqrt(2)) up to a
    # rational scale only when m2 is itself 1/(2k^2); use a generic quadratic
    # per point instead, D = numerator-cleared discriminant of m2.
    test_m2 = [Fr(1, 8), Fr(1, 16), Fr(1, 4), Fr(1, 3), Fr(1, 20)]
    for m2 in test_m2:
        # sqrt(m2): write m2 = a/b in lowest terms; sqrt(m2) = sqrt(ab)/b
        num, den = m2.numerator, m2.denominator
        Dloc = num * den  # sqrt(num/den) = sqrt(num*den)/den
        u = Quad(Dloc, 0, Fr(1, den))
        assert u * u == qi(Dloc, m2), (u * u, m2)
        m11 = 2 * qi(Dloc, THETA) * u * (1 - u)
        f1, f2 = m11, qi(Dloc, m2)
        f3 = 1 - f1 - f2
        for algo in ("proj", "frame"):
            T2 = t2_tensor_proj(Dloc, u, algo)
            assert T2[1] == qi(Dloc, 0), ("off-diagonal must vanish", algo, m2, T2)
            assert T2[0] == T2[2], ("not isotropic", algo, m2, T2)
            hs = hs_lo_field(Dloc, [f1, f2, f3], [qi(Dloc, s) for s in SIGMA])
            assert T2[0] == hs, ("does not attain HS_lo", algo, m2, T2[0], hs)
    print(f"selftest: T^2 structure attains HS_lo exactly at m1=m11(m2), "
          f"{len(test_m2)} values of m2, both lamination algorithms")

    # cross-check the two lamination algorithms agree on a fully generic
    # (non-diagonal, non-axis-normal) tree over Q(sqrt(3)), rational leaves
    D3 = 3
    A = (qi(D3, 2), qi(D3, Fr(1, 3)), qi(D3, 5))
    B = (qi(D3, 7), qi(D3, -1), qi(D3, 4))
    for n in [(1, 0), (0, 1), (2, 3), (1, -1)]:
        r1 = laminate_proj(A, B, qi(D3, Fr(2, 5)), n)
        r2 = laminate_frame(A, B, qi(D3, Fr(2, 5)), n, D3)
        assert r1 == r2, (n, r1, r2)
    print("selftest: laminate_proj and laminate_frame agree on a generic "
          "anisotropic tree, several normals")


# ---------------------------------------------------------------------------
# The actual target construction: sigma=(1,2,5), f=(1/8,1/8,3/4)
# ---------------------------------------------------------------------------


def build_target():
    D = 105
    x = Quad(D, Fr(-1, 26), Fr(1, 26))  # x = (-1+sqrt(105))/26, root of 13x^2+x-2=0
    assert 13 * x * x + x - 2 == qi(D, 0)

    g2 = x * x
    g1 = 1 - 7 * g2
    g3 = 6 * g2
    assert g1 + g2 + g3 == qi(D, 1)
    m11_check = 2 * qi(D, THETA) * x * (1 - x)
    assert m11_check == g1, "g1 must equal m11(g2)"

    a, b, c, d = t2_structure(D, x)
    for algo in ("proj", "frame"):
        T2 = t2_tensor_proj(D, x, algo)
        assert T2[1] == qi(D, 0) and T2[0] == T2[2]
        hs_inner = hs_lo_field(D, [g1, g2, g3], [qi(D, s) for s in SIGMA])
        assert T2[0] == hs_inner, (algo, T2[0], hs_inner)
    core = t2_tensor_proj(D, x, "proj")[0]  # isotropic scalar value

    Fcoat = qi(D, Fr(1, 8)) / g2       # core's fraction of the final whole
    c1 = (1 + Fcoat) / 2
    c2 = Fcoat / c1

    s1q = qi(D, SIGMA[0])
    T2tensor = (core, qi(D, 0), core)
    stage1 = laminate_proj(T2tensor, iso_leaf(D, SIGMA[0]), c1, E1)
    whole = laminate_proj(stage1, iso_leaf(D, SIGMA[0]), c2, E2)
    whole_frame = laminate_frame(T2tensor, iso_leaf(D, SIGMA[0]), c1, E1, D)
    whole_frame = laminate_frame(whole_frame, iso_leaf(D, SIGMA[0]), c2, E2, D)
    assert whole == whole_frame, "the two lamination algorithms disagree"

    target = qi(D, Fr(37, 11))
    assert whole[1] == qi(D, 0), "final tensor must be diagonal"
    assert whole[0] == whole[2] == target, (whole, target)

    # exact final volume fractions
    f2f = Fcoat * g2
    f3f = Fcoat * g3
    f1f = 1 - f2f - f3f
    assert (f1f, f2f, f3f) == (qi(D, Fr(1, 8)), qi(D, Fr(1, 8)), qi(D, Fr(3, 4)))

    return {
        "D": D, "x": x, "g": (g1, g2, g3), "abcd": (a, b, c, d),
        "core": core, "Fcoat": Fcoat, "c1": c1, "c2": c2, "value": whole[0],
    }


# ---------------------------------------------------------------------------
# Rational sandwich: snap sqrt(105) from below/above to a Fraction of huge
# denominator, rebuild the SAME tree with plain Fractions, and run it through
# the real harness (laminate.py + verify_laminate.py).  This does not add
# rigor to the exact proof above (which needs no floating point and no
# enclosure), but it lets the two REPO verifiers, which only understand
# Fraction trees, certify a numeric bracket around 37/11 as an independent
# sanity layer, per the contract's "every record passes verify_laminate.py".
# ---------------------------------------------------------------------------


def sqrt_bracket(n: int, digits_bits: int = 200):
    """Rational lo <= sqrt(n) <= hi via integer Newton on n * 4^k."""
    k = digits_bits
    scaled = n << (2 * k)
    r = isqrt(scaled)
    lo = Fr(r, 1 << k)
    hi = Fr(r + 1, 1 << k)
    return lo, hi


def isqrt(n: int) -> int:
    if n < 0:
        raise ValueError
    if n == 0:
        return 0
    x = 1 << ((n.bit_length() + 1) // 2)
    while True:
        y = (x + n // x) // 2
        if y >= x:
            return x
        x = y


def rational_tree(sqrt105: Fr):
    """The same tree as build_target(), with sqrt(105) replaced by a Fraction
    approximation -- gives an all-rational tree the real harness can run."""
    x = (Fr(-1, 26) + sqrt105 / 26)
    g2 = x * x
    g1 = 1 - 7 * g2
    g3 = 6 * g2
    u = x
    a = 1 - u
    c = 1 - u
    d = THETA
    b = 1 - u / 4

    def leaf(p):
        return {"phase": p}

    def node(frac, normal, layers):
        return {"fraction": str(frac), "normal": list(normal), "layers": layers}

    X = node(b, E2, [leaf("p3"), leaf("p1")])
    Z = node(d, E1, [leaf("p1"), leaf("p3")])
    Y = node(c, E2, [Z, leaf("p2")])
    T2 = node(a, E1, [X, Y])

    Fcoat = Fr(1, 8) / g2
    c1 = (1 + Fcoat) / 2
    c2 = Fcoat / c1
    stage1 = node(c1, E1, [T2, leaf("p1")])
    whole = node(c2, E2, [stage1, leaf("p1")])
    return whole, {"p1": SIGMA[0], "p2": SIGMA[1], "p3": SIGMA[2]}


def build_and_verify(out_dir=DATA):
    os.makedirs(out_dir, exist_ok=True)
    res = build_target()
    print("Exact symbolic construction (Q(sqrt(105))):")
    print(f"  inner T^2 fractions g = (g1, g2, g3) = "
          f"({res['g'][0]}, {res['g'][1]}, {res['g'][2]})")
    print(f"  inner T^2 node fractions a=c={res['abcd'][0]}, b={res['abcd'][1]}, d={res['abcd'][2]}")
    print(f"  outer coat fraction F = {res['Fcoat']}, c1={res['c1']}, c2={res['c2']}")
    print(f"  final tensor value = {res['value']}  (== 37/11 exactly: "
          f"{res['value'] == qi(105, Fr(37, 11))})")

    # rational sandwich, width ~ 2^-190
    lo, hi = sqrt_bracket(105, digits_bits=190)
    assert lo * lo < 105 < hi * hi  # as Fractions this compares exactly
    results = {}
    for tag, approx in (("lo", lo), ("hi", hi)):
        tree, phases = rational_tree(approx)
        t = L.effective(tree, phases)
        assert t[1] == 0
        fr = L.fractions_of(tree)
        tol = Fr(1, 10**50)
        assert abs(fr["p2"] - Fr(1, 8)) < tol, fr["p2"]
        assert abs(fr["p3"] - Fr(3, 4)) < tol, fr["p3"]
        assert abs(fr["p1"] - Fr(1, 8)) < tol, fr["p1"]
        assert L.keller_check(tree, phases)
        rec = {
            "tree": tree,
            "phases": {k: str(v) for k, v in phases.items()},
            "tensor": [str(x) for x in t],
            "fractions": {k: str(v) for k, v in fr.items()},
            "wiener": [str(x) for x in L.wiener_bounds([fr[n] for n in sorted(fr)],
                                                        [phases[n] for n in sorted(fr)])],
            "hs": [str(x) for x in L.hs_bounds([fr[n] for n in sorted(fr)],
                                                [phases[n] for n in sorted(fr)])],
        }
        path = os.path.join(out_dir, f"cherkaev_rational_{tag}.json")
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(rec, fh, indent=1)
        problems = V.verify_record(path)
        assert not problems, problems
        results[tag] = t[0]
        print(f"  rational snap ({tag}, sqrt105 approx width 2^-190): "
              f"tensor[0]-37/11 = {t[0]-Fr(37,11)}  ({float(t[0]-Fr(37,11)):.3e}); "
              f"verify_laminate.py: PASS")

    # NOTE: the map sqrt(105) -> tensor value need not be monotone through
    # this whole pipeline, so (lo, hi) is not guaranteed to bracket 37/11 in
    # that order; report both residuals as independent witnesses instead
    # (the exact symbolic proof above already settles the value with no
    # floating point or enclosure needed -- this is a numeric sanity layer,
    # run through the real harness, not additional rigor).
    span = max(abs(results["lo"] - Fr(37, 11)), abs(results["hi"] - Fr(37, 11)))
    print(f"  both rational-snap witnesses agree with 37/11 to within "
          f"{float(span):.3e} at sqrt(105) precision 2^-190")

    # write the symbolic record (irrational fields, so JSON stores algebraic
    # descriptions, not a Fraction-only tree the harness can parse directly)
    sym = {
        "problem": "three-phase-conductivity",
        "target": {"sigma": ["1", "2", "5"], "f": ["1/8", "1/8", "3/4"], "hs_lo": "37/11"},
        "field": "Q(sqrt(105))",
        "x": {"p": str(res["x"].p), "q": str(res["x"].q), "note": "x = (-1+sqrt(105))/26, root of 13x^2+x-2=0"},
        "inner_T2_fractions_g1_g2_g3": [str(v) for v in res["g"]],
        "inner_T2_node_fractions": {"a": [str(res["abcd"][0].p), str(res["abcd"][0].q)],
                                     "b": [str(res["abcd"][1].p), str(res["abcd"][1].q)],
                                     "c": [str(res["abcd"][2].p), str(res["abcd"][2].q)],
                                     "d": [str(res["abcd"][3].p), str(res["abcd"][3].q)]},
        "coating": {"F": [str(res["Fcoat"].p), str(res["Fcoat"].q)],
                    "c1": [str(res["c1"].p), str(res["c1"].q)],
                    "c2": [str(res["c2"].p), str(res["c2"].q)]},
        "final_value": [str(res["value"].p), str(res["value"].q)],
        "final_value_is_exactly_37_11": bool(res["value"] == qi(105, Fr(37, 11))),
        "rational_sandwich": {
            "sqrt105_lo": str(lo), "sqrt105_hi": str(hi),
            "tensor_lo": str(results["lo"]), "tensor_hi": str(results["hi"]),
            "max_residual_from_37_11": str(span),
        },
    }
    with open(os.path.join(out_dir, "cherkaev_target.json"), "w", encoding="utf-8") as fh:
        json.dump(sym, fh, indent=1)
    print(f"\nrecords written to {out_dir}")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--build", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    if a.build:
        build_and_verify()
    if not a.selftest and not a.build:
        ap.print_help()
