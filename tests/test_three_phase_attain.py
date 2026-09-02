"""The closed-form rank-5 laminate that attains the HS lower bound at f1 = m11.

These fail if the construction stops attaining, if its parameter ranges break,
or if it stops beating Milton's classical threshold.
"""

from __future__ import annotations

import sys
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "harness" / "three-phase-conductivity"))
sys.path.insert(0, str(ROOT / "problems" / "three-phase-conductivity" / "explore"))

import laminate as L  # noqa: E402
import tp_attain as A  # noqa: E402
import tp_fields as TF  # noqa: E402

SIG = (Fr(1), Fr(2), Fr(5))
PHASES = {"p1": SIG[0], "p2": SIG[1], "p3": SIG[2]}


def test_worked_instance_attains_exactly():
    ok, tree, hs = A.attains(SIG, Fr(1, 2), Fr(4, 5))
    assert ok and hs == 3
    assert A.fractions_for(SIG, Fr(1, 2)) == (Fr(1, 8), Fr(1, 4), Fr(5, 8))
    assert A.params(SIG, Fr(1, 2), Fr(4, 5)) == (
        Fr(1, 2), Fr(1, 8), Fr(4, 5), Fr(3, 8), Fr(3, 4))
    assert L.effective(tree, PHASES) == (Fr(3), Fr(0), Fr(3))


def test_a3_is_a_free_parameter():
    """A one-parameter family, every member attaining."""
    for a3 in [Fr(3, 5), Fr(2, 3), Fr(4, 5), Fr(7, 8), Fr(9, 10)]:
        ok, _, _ = A.attains(SIG, Fr(1, 2), a3)
        assert ok, a3


def test_beats_miltons_threshold():
    """m11 = 2*Theta*r*(1-r) is strictly below 2*Theta*(1-m2) for r in (0,1)."""
    for sig in [SIG, (Fr(1), Fr(3), Fr(7)), (Fr(2), Fr(5), Fr(11))]:
        th = A.theta_of(sig)
        for r in [Fr(1, 4), Fr(1, 3), Fr(1, 2), Fr(2, 3), Fr(4, 5)]:
            f1, f2, _ = A.fractions_for(sig, r)
            assert f1 == 2 * th * r * (1 - r)
            assert f1 < 2 * th * (1 - f2)


def test_attains_across_conductivity_triples():
    cases = [((Fr(1), Fr(3), Fr(7)), Fr(1, 2)),
             ((Fr(2), Fr(5), Fr(11)), Fr(1, 2)),
             ((Fr(1), Fr(2), Fr(100)), Fr(1, 4)),
             ((Fr(1), Fr(4), Fr(5)), Fr(1, 2)),
             ((Fr(1, 3), Fr(5, 2), Fr(9)), Fr(2, 5)),
             ((Fr(1), Fr(2), Fr(5)), Fr(1, 3))]
    for sig, r in cases:
        a3 = (1 - r + 1) / 2
        f1, f2, f3 = A.fractions_for(sig, r)
        assert f3 > 0
        ok, _, _ = A.attains(sig, r, a3)
        assert ok, (sig, r)


def test_derived_parameter_identities():
    """a1 = r*Theta and a5 = (s2 s3 - s1^2)/((s2+s1)(s3-s1)), both in (0,1)."""
    for sig in [SIG, (Fr(1), Fr(3), Fr(7)), (Fr(2), Fr(5), Fr(11)),
                (Fr(1, 3), Fr(5, 2), Fr(9))]:
        s1, s2, s3 = sig
        a5_expected = (s2 * s3 - s1 * s1) / ((s2 + s1) * (s3 - s1))
        assert 0 < a5_expected < 1
        for r in [Fr(1, 3), Fr(1, 2), Fr(3, 5)]:
            a0, a1, a3, a4, a5 = A.params(sig, r, (1 - r + 1) / 2)
            assert a1 == r * A.theta_of(sig)
            assert a5 == a5_expected
            assert 0 < a0 < 1 and 0 < a1 < 1 and 0 < a4 < 1


def test_attaining_structure_satisfies_the_002_field_conditions():
    """Phases 2 and 3 uniform exactly on target; the gap decomposition is tight."""
    E0 = (Fr(1), Fr(1))
    _, tree, _ = A.attains(SIG, Fr(1, 2), Fr(4, 5))
    fr = L.fractions_of(tree)
    names = sorted(fr)
    _, targets = TF.hs_target_fields([fr[n] for n in names],
                                     [PHASES[n] for n in names], E0)
    got = TF.per_phase(tree, PHASES, E0)
    for i, nm in enumerate(names):
        if nm == "p1":
            continue
        vals = {E for _, E in got[nm]}
        assert len(vals) == 1, nm
        assert vals == {targets[i]}, nm
    d = TF.gap_decomposition(tree, PHASES, E0)
    assert d["gap"] == 0
    assert d["mean_term"] == 0
    assert d["var_term"] == d["var_required"]
