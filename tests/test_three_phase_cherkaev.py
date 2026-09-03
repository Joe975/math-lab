"""The coated-T2 laminate attaining HS_lo at the original target point.

Guards attempt 007: the construction lives in Q(sqrt(105)), but its two
rational-snap witnesses are plain-Fraction trees the harness can check, so the
attainment is testable without the quadratic-field machinery.
"""

from __future__ import annotations

import json
import sys
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "harness" / "three-phase-conductivity"))

import laminate as L  # noqa: E402
import verify_laminate as V  # noqa: E402

DATA = ROOT / "problems" / "three-phase-conductivity" / "data" / "cherkaev"
TARGET = Fr(37, 11)


def _records():
    for name in ("cherkaev_rational_lo.json", "cherkaev_rational_hi.json"):
        path = DATA / name
        if path.exists():
            yield name, json.loads(path.read_text(encoding="utf-8"))


def test_rational_snaps_exist():
    assert list(_records()), "expected the two rational-snap witnesses"


def test_snaps_are_isotropic_and_hit_the_bound():
    """Value equals 37/11 to the precision of the sqrt(105) approximation."""
    for name, rec in _records():
        phases = {k: Fr(v) for k, v in rec["phases"].items()}
        tree = rec["tree"]
        e = L.effective(tree, phases)
        assert e == V.effective(tree, phases), name
        assert L.keller_check(tree, phases), name
        # isotropic to the snap precision
        assert abs(e[1]) < Fr(1, 10 ** 40), name
        assert abs(e[0] - e[2]) < Fr(1, 10 ** 40), name
        # and on the bound to the same precision
        assert abs(e[0] - TARGET) < Fr(1, 10 ** 50), (name, float(e[0] - TARGET))


def test_snaps_carry_the_target_volume_fractions():
    for name, rec in _records():
        fr = {k: Fr(v) for k, v in rec["fractions"].items()}
        for phase, want in (("p1", Fr(1, 8)), ("p2", Fr(1, 8)), ("p3", Fr(3, 4))):
            assert abs(fr[phase] - want) < Fr(1, 10 ** 40), (name, phase)


def test_the_bound_at_the_target_is_37_over_11():
    lo, _ = L.hs_bounds([Fr(1, 8), Fr(1, 8), Fr(3, 4)], [Fr(1), Fr(2), Fr(5)])
    assert lo == TARGET
