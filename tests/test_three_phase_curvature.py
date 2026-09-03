"""The below-m11 curve as a ratio of curvatures at m11 (attempt 017).

These fail if the second-order tangency at m11 breaks, if the closed forms for
beta, gamma and c stop reproducing the jet expansion, or if c leaves (0, 1).
"""

from __future__ import annotations

import sys
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "harness" / "three-phase-conductivity"))
sys.path.insert(0, str(ROOT / "problems" / "three-phase-conductivity" / "explore"))

import tp_curvature as TC  # noqa: E402

SIG = (Fr(1), Fr(2), Fr(5))
HALF = Fr(1, 2)


def test_second_order_tangency_at_m11():
    """B2 and the family both meet HS_lo at m11 to FIRST order, so the ratio
    is a ratio of second derivatives."""
    R = TC.expand(SIG, HALF, deg=2)
    assert R["g"].coeff(0, 0) == 0 and R["g"].coeff(1, 0) == 0
    assert R["g"].coeff(0, 1) == 0
    assert R["d"].coeff(0, 0) == 0 and R["d"].coeff(1, 0) == 0
    assert R["beta"] > 0 and R["gamma"] > 0


def test_c_exact_at_the_worked_instance():
    assert TC.expand(SIG, HALF, deg=2)["c"] == Fr(189, 355)
    assert TC.closed_forms(SIG, HALF)["c"] == Fr(189, 355)


def test_closed_forms_match_the_jet_expansion():
    for trip in [(1, 2, 5), (1, 3, 7), (1, 4, 9), (1, 2, 9), (2, 5, 11),
                 (2, 7, 9), (3, 5, 8)]:
        sig = tuple(Fr(x) for x in trip)
        for r in (Fr(1, 3), Fr(1, 2), Fr(3, 5), Fr(3, 4)):
            J = TC.expand(sig, r, deg=2)
            C = TC.closed_forms(sig, r)
            assert (J["beta"], J["gamma"], J["c"]) == (C["beta"], C["gamma"],
                                                       C["c"]), (trip, r)


def test_c_depends_on_m2_not_only_on_sigma():
    """015 recorded c as a function of sigma alone; it moves with m2 too."""
    cs = [TC.closed_forms(SIG, r)["c"] for r in (Fr(1, 4), Fr(1, 2), Fr(3, 4))]
    assert cs[0] < cs[1] < cs[2]
    assert float(cs[2]) - float(cs[0]) > 0.2


def test_c_strictly_between_zero_and_one():
    for k2n in range(3, 20):
        for k3n in range(k2n + 1, 40):
            x, y = Fr(k2n, 2), Fr(k3n, 2)
            for rn in (1, 3, 6, 9, 11):
                r = Fr(rn, 12)
                th = (y - x) / ((x + 1) * (y - 1))
                m11 = 2 * th * r * (1 - r)
                if m11 <= 0 or 1 - m11 - r * r <= 0:
                    continue
                c = TC.closed_forms((Fr(1), x, y), r)["c"]
                assert 0 < c < 1, (x, y, r, c)


def test_c_approaches_one_only_at_infinite_contrast():
    r = Fr(19, 20)
    cs = [TC.closed_forms((Fr(1), Fr(3, 2), Fr(y)), r)["c"]
          for y in (10, 100, 1000, 10 ** 5)]
    assert all(a < b for a, b in zip(cs, cs[1:]))
    assert cs[-1] < 1 and float(cs[-1]) > 0.999
