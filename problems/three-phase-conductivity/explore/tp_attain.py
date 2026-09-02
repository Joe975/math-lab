#!/usr/bin/env python3
"""A closed-form rank-5 laminate attaining the 2D three-phase Hashin-Shtrikman
LOWER bound at the threshold fraction f1 = m11, for arbitrary phases.

This is the object the census of attempts 001-003 was circling.  With

    Theta = s1 (s3 - s2) / ((s2 + s1)(s3 - s1)),        r = sqrt(m2),

set the volume fractions to

    f2 = r^2,      f1 = m11 = 2 Theta r (1 - r),      f3 = 1 - f1 - f2,

and build the tree

    A    = laminate(p1 at a1, p3),   normal e2
    D    = laminate(p3 at a5, p1),   normal e1
    C    = laminate(p2 at a4, D),    normal e2
    B    = laminate(C  at a3, p2),   normal e2
    root = laminate(A  at a0, B),    normal e1

with

    a0 = 1 - r
    a1 = r * Theta
    a3 free in (1 - r, 1)                      <-- a ONE-PARAMETER FAMILY
    a4 = (a3 - (1 - r)) / a3
    a5 = (s2 s3 - s1^2) / ((s2 + s1)(s3 - s1))

Then sigma* is isotropic and equals HS_lo EXACTLY, in Q.

Why this matters.  Milton's classical coated assemblage attains the bound only
for f1 >= 2 Theta (1 - m2); here f1 = 2 Theta r (1 - r) = 2 Theta sqrt(m2)
(1 - sqrt(m2)), which is strictly smaller for every m2 in (0,1).  Together with
the coating lemma of attempt 003 -- coating an HS-optimal structure with extra
comparison-medium phase preserves HS-optimality, so attainability is upward
closed in f1 along a fixed f2:f3 ray -- this construction gives attainment for
ALL f1 >= m11.  That is a constructive route to the statement Cherkaev (2009)
Theorem 8.1 makes for his structure L13,2,13,1,1.

Parameter sanity, proved rather than assumed (see selftest):
  * a5 in (0,1): s2 s3 > s1^2 gives a5 > 0, and s3 > s2 gives a5 < 1.
  * a1 in (0,1): a1 = r Theta with Theta in (0, 1) for s1 < s2 < s3.
  * a4 in (0,1) exactly when a3 > 1 - r, which is the stated range.
So the only genuine restriction is a3 > 1 - r, plus f3 > 0.

Usage:
    python tp_attain.py --selftest
    python tp_attain.py --sigma 1,2,5 --root 1/2 [--a3 4/5] [--json OUT]
Standard library only.
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
import laminate as L  # noqa: E402
import verify_laminate as V  # noqa: E402


def theta_of(sig):
    s1, s2, s3 = sig
    return s1 * (s3 - s2) / ((s2 + s1) * (s3 - s1))


def params(sig, r: Fr, a3: Fr):
    """The closed-form lamination fractions (a0, a1, a3, a4, a5)."""
    s1, s2, s3 = sig
    th = theta_of(sig)
    a0 = 1 - r
    a1 = r * th
    a5 = (s2 * s3 - s1 * s1) / ((s2 + s1) * (s3 - s1))
    a4 = (a3 - (1 - r)) / a3
    return a0, a1, a3, a4, a5


def fractions_for(sig, r: Fr):
    """(f1, f2, f3) this construction is built for: f1 = m11, f2 = r^2."""
    th = theta_of(sig)
    f1 = 2 * th * r * (1 - r)
    f2 = r * r
    return f1, f2, 1 - f1 - f2


def tree_for(sig, r: Fr, a3: Fr):
    a0, a1, a3, a4, a5 = params(sig, r, a3)
    if not all(0 < x < 1 for x in (a0, a1, a3, a4, a5)):
        return None
    A = {"fraction": str(a1), "normal": [0, 1],
         "layers": [{"phase": "p1"}, {"phase": "p3"}]}
    D = {"fraction": str(a5), "normal": [1, 0],
         "layers": [{"phase": "p3"}, {"phase": "p1"}]}
    C = {"fraction": str(a4), "normal": [0, 1], "layers": [{"phase": "p2"}, D]}
    B = {"fraction": str(a3), "normal": [0, 1], "layers": [C, {"phase": "p2"}]}
    return {"fraction": str(a0), "normal": [1, 0], "layers": [A, B]}


def attains(sig, r: Fr, a3: Fr):
    """True iff the tree is isotropic and equals HS_lo exactly, both routes."""
    tree = tree_for(sig, r, a3)
    if tree is None:
        return False, None, None
    f1, f2, f3 = fractions_for(sig, r)
    if f1 <= 0 or f3 <= 0:
        return False, None, None
    phases = {"p1": sig[0], "p2": sig[1], "p3": sig[2]}
    e = L.effective(tree, phases)
    fr = L.fractions_of(tree)
    hs, _ = L.hs_bounds([f1, f2, f3], list(sig))
    good = (e[0] - hs == 0 and e[1] == 0 and e[0] == e[2]
            and fr["p1"] == f1 and fr["p2"] == f2
            and e == V.effective(tree, phases)
            and L.keller_check(tree, phases))
    return good, tree, hs


def record_for(sig, r: Fr, a3: Fr):
    tree = tree_for(sig, r, a3)
    phases = {"p1": sig[0], "p2": sig[1], "p3": sig[2]}
    e = L.effective(tree, phases)
    fr = L.fractions_of(tree)
    names = sorted(fr)
    fracs = [fr[n] for n in names]
    sigs = [phases[n] for n in names]
    return {
        "tree": tree,
        "phases": {k: str(v) for k, v in phases.items()},
        "tensor": [str(x) for x in e],
        "fractions": {k: str(v) for k, v in fr.items()},
        "wiener": [str(x) for x in L.wiener_bounds(fracs, sigs)],
        "hs": [str(x) for x in L.hs_bounds(fracs, sigs)],
    }


def selftest() -> None:
    F = Fr
    # 1. the worked instance, fully exact
    sig = (F(1), F(2), F(5))
    ok, tree, hs = attains(sig, F(1, 2), F(4, 5))
    assert ok and hs == 3, (ok, hs)
    a0, a1, a3, a4, a5 = params(sig, F(1, 2), F(4, 5))
    assert (a0, a1, a5) == (F(1, 2), F(1, 8), F(3, 4)), (a0, a1, a5)
    assert fractions_for(sig, F(1, 2)) == (F(1, 8), F(1, 4), F(5, 8))

    # 2. a3 is genuinely free: a one-parameter family, all attaining
    for a3 in [F(3, 5), F(2, 3), F(4, 5), F(7, 8), F(9, 10)]:
        ok, _, _ = attains(sig, F(1, 2), a3)
        assert ok, a3

    # 3. f1 = m11 is strictly below Milton's threshold 2*Theta*(1-m2)
    for r in [F(1, 3), F(1, 2), F(2, 3), F(3, 4)]:
        th = theta_of(sig)
        f1, f2, _ = fractions_for(sig, r)
        assert f1 == 2 * th * r * (1 - r)
        assert f1 < 2 * th * (1 - f2), (f1, 2 * th * (1 - f2))

    # 4. many conductivity triples and roots, exact
    cases = [((F(1), F(3), F(7)), F(1, 2)), ((F(2), F(5), F(11)), F(1, 2)),
             ((F(1), F(2), F(100)), F(1, 4)), ((F(1), F(4), F(5)), F(1, 2)),
             ((F(1, 3), F(5, 2), F(9)), F(2, 5)), ((F(1), F(2), F(5)), F(1, 3)),
             ((F(3, 2), F(7, 2), F(20)), F(3, 5))]
    n = 0
    for s, r in cases:
        for a3 in [(1 - r + 1) / 2, F(9, 10)]:
            if a3 <= 1 - r or a3 >= 1:
                continue
            ok, _, _ = attains(s, r, a3)
            f1, f2, f3 = fractions_for(s, r)
            if f3 <= 0:
                continue
            assert ok, (s, r, a3)
            n += 1
    assert n >= 8, n

    # 5. the parameter-range claims in the docstring
    for s, r in cases:
        s1, s2, s3 = s
        a5 = (s2 * s3 - s1 * s1) / ((s2 + s1) * (s3 - s1))
        assert 0 < a5 < 1
        assert 0 < theta_of(s) < 1
    print("selftest: closed-form construction attains HS_lo exactly in every case")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--sigma", default="1,2,5")
    ap.add_argument("--root", default="1/2", help="r = sqrt(m2)")
    ap.add_argument("--a3", default=None, help="free parameter in (1-r, 1)")
    ap.add_argument("--json", metavar="PATH")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        return
    sig = tuple(Fr(x) for x in a.sigma.split(","))
    r = Fr(a.root)
    a3 = Fr(a.a3) if a.a3 else (1 - r + 1) / 2
    ok, tree, hs = attains(sig, r, a3)
    f1, f2, f3 = fractions_for(sig, r)
    print(f"sigma = {tuple(str(x) for x in sig)}, r = {r}, a3 = {a3}")
    print(f"fractions (f1=m11, f2, f3) = ({f1}, {f2}, {f3})")
    print(f"Theta = {theta_of(sig)}, Milton threshold 2*Theta*(1-m2) = {2*theta_of(sig)*(1-f2)}")
    print(f"params (a0,a1,a3,a4,a5) = {tuple(str(x) for x in params(sig, r, a3))}")
    print(f"HS_lo = {hs};  attains exactly: {ok}")
    if a.json and tree is not None:
        with open(a.json, "w", encoding="utf-8") as fh:
            json.dump(record_for(sig, r, a3), fh, indent=1)
        print(f"record written to {a.json}")


if __name__ == "__main__":
    main()
