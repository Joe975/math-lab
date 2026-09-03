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


def test_family_attains_at_m11_and_not_below():
    """Within this rank-5 family m11 is exactly the threshold: the ISOTROPIC
    member attains at f1 = m11 and lies strictly above the bound just below it.
    EVIDENCE about this family and parametrisation only, not all microstructures.
    """
    r = Fr(1, 2)
    f2 = r * r
    m11 = 2 * A.theta_of(SIG) * r * (1 - r)
    assert m11 == Fr(1, 8)
    ok, _, _ = A.attains(SIG, r, Fr(4, 5))
    assert ok

    a0, _, a3, a4, _ = A.params(SIG, r, Fr(4, 5))
    Q = (1 - a0) - f2

    def build(f1, a5):
        a1 = (f1 - Q * (1 - a5)) / a0
        if not all(0 < x < 1 for x in (a0, a1, a3, a4, a5)):
            return None
        return {"fraction": str(a0), "normal": [1, 0], "layers": [
            {"fraction": str(a1), "normal": [0, 1],
             "layers": [{"phase": "p1"}, {"phase": "p3"}]},
            {"fraction": str(a3), "normal": [0, 1], "layers": [
                {"fraction": str(a4), "normal": [0, 1], "layers": [
                    {"phase": "p2"},
                    {"fraction": str(a5), "normal": [1, 0],
                     "layers": [{"phase": "p3"}, {"phase": "p1"}]}]},
                {"phase": "p2"}]}]}

    def isotropic_member(f1, iters=90):
        """Bisect a5 for sigma11 == sigma22; returns the exact tree or None."""
        def aniso(x):
            t = build(f1, x)
            if t is None:
                return None
            e = L.effective(t, PHASES)
            return e[0] - e[2]
        lo, hi = max(Fr(0), 1 - f1 / Q) + Fr(1, 10 ** 6), 1 - Fr(1, 10 ** 6)
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
        return build(f1, (lo + hi) / 2)

    for f1 in [Fr(3, 25), Fr(11, 100)]:      # 0.12 and 0.11, both below m11
        t = isotropic_member(f1)
        assert t is not None, f1
        fr = L.fractions_of(t)
        assert fr["p1"] == f1
        names = sorted(fr)
        lo_b, _ = L.hs_bounds([fr[n] for n in names], [PHASES[n] for n in names])
        e = L.effective(t, PHASES)
        # near-isotropic by construction; compare the quadratic form
        val = (e[0] + e[2]) / 2
        assert val > lo_b, (f1, val, lo_b)


def test_parameters_are_derived_from_the_forced_fields():
    """a5 and a1 follow from tangential continuity, not from pattern-matching.

    From attempt 010 the attainment fields are forced: phases 2 and 3 isotropic
    at rho_i = c/(s_i+s1), phase 1 carrying the mirror pair with y = s3*rho_3/s1.
    Then:
      In Y = lam(D at w, p2; e2) the e1 components are tangential and equal, so
      D's e1 component a5*rho_3 + (1-a5)*y must equal rho_2, giving
          a5 = (y - rho_2)/(y - rho_3).
      At the root lam(A at a0, Y; e1) the e2 components are tangential and equal,
      and the applied field is Z0 = I so that shared value is 1, giving
          a1 = (1 - rho_3)/(y - rho_3).
    """
    cases = [((Fr(1), Fr(2), Fr(5)), Fr(1, 2)),
             ((Fr(1), Fr(3), Fr(7)), Fr(1, 2)),
             ((Fr(2), Fr(5), Fr(11)), Fr(1, 2)),
             ((Fr(1), Fr(4), Fr(9)), Fr(1, 3)),
             ((Fr(1, 3), Fr(5, 2), Fr(9)), Fr(2, 5))]
    for sig, r in cases:
        a3 = (1 - r + 1) / 2
        ok, _, hs = A.attains(sig, r, a3)
        assert ok, (sig, r)
        s1, s2, s3 = sig
        c = hs + s1
        rho2, rho3 = c / (s2 + s1), c / (s3 + s1)
        y = s3 * rho3 / s1
        _, a1_built, _, _, a5_built = A.params(sig, r, a3)
        assert (y - rho2) / (y - rho3) == a5_built, (sig, r)
        assert (1 - rho3) / (y - rho3) == a1_built, (sig, r)


def test_a5_reduction_is_an_identity():
    """(y - rho_2)/(y - rho_3) reduces to (s2 s3 - s1^2)/((s2+s1)(s3-s1)) = 1-Theta."""
    for sig in [(Fr(1), Fr(2), Fr(5)), (Fr(1), Fr(3), Fr(7)),
                (Fr(2), Fr(5), Fr(11)), (Fr(1, 3), Fr(5, 2), Fr(9))]:
        s1, s2, s3 = sig
        for c in [Fr(3), Fr(9, 2), Fr(77, 10), Fr(1, 2)]:
            rho2, rho3 = c / (s2 + s1), c / (s3 + s1)
            y = s3 * rho3 / s1
            if y == rho3:
                continue
            assert ((y - rho2) / (y - rho3)
                    == (s2 * s3 - s1 * s1) / ((s2 + s1) * (s3 - s1))
                    == 1 - A.theta_of(sig))
