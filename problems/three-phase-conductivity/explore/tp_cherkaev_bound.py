#!/usr/bin/env python3
"""The Cherkaev (2009) three-material lower bound, TRANSCRIBED [T], plus the
internal consistency tests it passes and the discrepancy it currently shows.

TRANSCRIPTION NOTICE.  Every formula below is [T] -- taken from a machine
transcription of Cherkaev 2009 (Mech. Mater. 41, 411-433; preprint
arXiv:1009.3060) Theorem 7.1 and eqs 7.5-7.18, recorded in
`data/literature-check.md`.  The PDF was NOT read directly by this repo.  Per
`PROBLEM.md`'s verification contract a transcribed bound must be re-derived or
cross-checked before anything is killed against it, so nothing here kills
anything: this module exists to make the transcription reproducible and to
record precisely what it does and does not survive.

The bound is piecewise in m1, with B1 equal to the Hashin-Shtrikman lower
bound:

    B(m1, m2) = B1  for m1 >= m11
                B2  for m12 <= m1 <= m11
                B3  for m1 <= m12

THE OPEN ESCALATION.  At sigma = (1,2,5), m2 = 1/4 (so sqrt(m2) = 1/2 is
rational and everything is exact in Q), laminates built in this repo sit
STRICTLY BELOW B2 just under m11 = 1/8 -- by 5.4e-06 at m1 = 31/250, 1.4e-04
at 3/25 and 1.4e-03 at 11/100.  Those structures carry exact volume fractions,
an off-diagonal of exactly 0, a diagonal difference below 1e-61, agreement
between both harness algebra routes, a passing Keller-Dykhne identity and a
passing verify_laminate.py run.  A structure below a valid lower bound is
impossible.

What makes this hard to dismiss as a transcription slip is that the transcribed
B2 passes three independent internal tests (all reproduced by --selftest):

  A. continuity at the upper boundary: B2(m11) == HS_lo(m11) EXACTLY in Q;
  B. continuity at the lower boundary: B2(m12) - B3(m12) ~ 1e-34, i.e. zero;
  C. B2 >= HS_lo throughout its own region, as a tighter bound must be.

An adjudication is running (`data/bound-adjudication.md`).  Until it returns,
NOTHING is claimed below m11 in either direction.

Usage:
    python tp_cherkaev_bound.py --selftest
Standard library only.
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "three-phase-conductivity"))
import laminate as L  # noqa: E402


def hs_lo(k, m1: Fr, m2: Fr) -> Fr:
    k1, k2, k3 = k
    return L.hs_bounds([m1, m2, 1 - m1 - m2], [k1, k2, k3])[0]


def theta(k):
    k1, k2, k3 = k
    return k1 * (k3 - k2) / ((k2 + k1) * (k3 - k1))


def m11_of(k, s: Fr) -> Fr:
    """[T] m11 = 2 sqrt(m2) (1 - sqrt(m2)) k1 (k3-k2) / ((k3-k1)(k1+k2))."""
    k1, k2, k3 = k
    return 2 * s * (1 - s) * k1 * (k3 - k2) / ((k3 - k1) * (k1 + k2))


def m12_of(k, s: Fr):
    """[T] eqs 7.5-7.18. Involves sqrt(Z2), so generally irrational."""
    k1, k2, k3 = k
    m2 = s * s
    Z3 = k2 * (k1 - k2) * (k1 - k3) * (k1 - k2 + 2 * k3)
    Z4 = (k1 - k2) ** 2 * (k1 ** 2 + 6 * k1 * k2 - 4 * k1 * k3
                           - 4 * k2 * k3 + 4 * k3 ** 2 + k2 ** 2)
    Z2 = 4 * k2 ** 2 * (k3 - k1) ** 2 + 4 * s * Z3 + m2 * Z4
    root = Fr(math.sqrt(float(Z2))).limit_denominator(10 ** 12)
    Z0 = 2 * k2 * (k3 - k1) + s * (k1 + k2) * (2 * k3 - k1 - k2) - root
    return (1 - s) / (4 * k2 * (k3 - k1)) * Z0, Z2


def B2(k, m1: Fr, s: Fr) -> Fr:
    """[T] B2 = k2 + (1-sqrt(m2))^2 Z5/Z6."""
    k1, k2, k3 = k
    m2 = s * s
    m3 = 1 - m1 - m2
    Z5 = m1 * k1 ** 2 - m1 * k2 ** 2 + 2 * m3 * k1 * (k3 - k2)
    Z6 = ((1 - s) ** 2 + (1 - m1 - s) ** 2) * k1 + m1 * (1 - s) ** 2 * k2 + m1 * m3 * k3
    return k2 + (1 - s) ** 2 * Z5 / Z6


def B3(k, m1: Fr, s: Fr) -> Fr:
    """[T] B3 = -k2 + 1/(m2/(2 k2) + Z7)."""
    k1, k2, k3 = k
    m2 = s * s
    m3 = 1 - m1 - m2
    Z7 = (((k1 - k2) * m1 ** 2 + (2 * k1 - k2 + k3) * m1 * m3 + 2 * k1 * m3 ** 2)
          / ((k1 ** 2 - k2 ** 2) * m1 + 2 * k1 * (k2 + k3) * m3))
    return -k2 + 1 / (m2 / (2 * k2) + Z7)


def bound(k, m1: Fr, s: Fr):
    """(region, value) per the piecewise statement."""
    m11 = m11_of(k, s)
    m12, _ = m12_of(k, s)
    if m1 >= m11:
        return "B1 (=HS)", hs_lo(k, m1, s * s)
    if m1 >= m12:
        return "B2", B2(k, m1, s)
    return "B3", B3(k, m1, s)


def selftest() -> None:
    F = Fr
    k = (F(1), F(2), F(5))
    s = F(1, 2)          # sqrt(m2) rational, so tests A and C are exact in Q
    m2 = s * s
    m11 = m11_of(k, s)
    assert m11 == F(1, 8), m11
    assert m11 == 2 * theta(k) * s * (1 - s)

    # A. continuity at the upper boundary, EXACT in Q
    assert B2(k, m11, s) == hs_lo(k, m11, m2), "B2 must meet HS at m11"

    # B. continuity at the lower boundary (m12 is irrational; float-level)
    m12, _ = m12_of(k, s)
    assert abs(float(B2(k, m12, s) - B3(k, m12, s))) < 1e-20, "B2 must meet B3 at m12"

    # C. B2 >= HS throughout its own region, exactly at every sampled point
    for i in range(1, 200):
        m1 = m12 + (m11 - m12) * F(i, 200)
        assert B2(k, m1, s) >= hs_lo(k, m1, m2), m1

    # D. the recorded discrepancy, so it cannot silently disappear
    ours = {F(31, 250): F("3.0053470065").limit_denominator(10 ** 10),
            F(3, 25): F("3.0270098722").limit_denominator(10 ** 10),
            F(11, 100): F("3.0831692377").limit_denominator(10 ** 10)}
    for m1, val in ours.items():
        region, b = bound(k, m1, s)
        assert region == "B2", (m1, region)
        assert val > hs_lo(k, m1, m2), "our structure must exceed HS"
        assert val < b, "the OPEN discrepancy: our structure is below B2"
    print("selftest: transcription passes continuity at m11 (exact) and m12, "
          "and B2 >= HS on its region; the below-B2 discrepancy is reproduced")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
