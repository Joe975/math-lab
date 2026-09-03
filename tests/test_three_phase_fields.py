"""Fast deterministic tests for the exact laminate field tool.

These fail if the field propagation, the Hashin-Shtrikman attainment-field
formula, or the Keller-Dykhne duality reduction is broken.
"""

from __future__ import annotations

import sys
from fractions import Fraction as Fr
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "harness" / "three-phase-conductivity"))
sys.path.insert(0, str(ROOT / "problems" / "three-phase-conductivity" / "explore"))

import laminate as L  # noqa: E402
import tp_coated as TC  # noqa: E402
import tp_fields as TF  # noqa: E402

SIGMA = [Fr(1), Fr(2), Fr(5)]
PHASES = {"p1": SIGMA[0], "p2": SIGMA[1], "p3": SIGMA[2]}
E0 = (Fr(1), Fr(1))

TREE = {
    "fraction": "1/3", "normal": [1, 0],
    "layers": [{"phase": "p1"},
               {"fraction": "3/7", "normal": [1, 2],
                "layers": [{"phase": "p2"}, {"phase": "p3"}]}],
}


def test_field_identities_exact():
    """Average field, average flux and energy all reproduce, in Q."""
    for e in [(Fr(1), Fr(0)), (Fr(0), Fr(1)), (Fr(1), Fr(1)), (Fr(3), Fr(-2))]:
        TF.check_identities(TREE, PHASES, e)


def test_milton_assemblage_realises_attainment_fields():
    """Phases 2 and 3 carry exactly E_i = (s_HS+s1)/(s_i+s1) E0."""
    f = [Fr(3, 8), Fr(2, 8), Fr(3, 8)]
    tree, ph, _ = TC.milton(f, SIGMA, 0)
    assert tree is not None
    hs_lo, targets = TF.hs_target_fields(f, SIGMA, E0)
    got = TF.per_phase(tree, ph, E0)
    for i, name in [(1, "p2"), (2, "p3")]:
        for _, E in got[name]:
            assert E == targets[i]
    # phase 1 is the comparison medium: free pointwise, but its volume
    # average is still pinned to the target.
    tot = sum(v for v, _ in got["p1"])
    avg = tuple(sum(v * E[k] for v, E in got["p1"]) / tot for k in (0, 1))
    assert avg == targets[0]


def test_middle_phase_field_condition_is_miltons_threshold():
    """E_2 <= E0 exactly when HS_lo <= sigma2, i.e. m1 >= 2*Theta*(1-m2)."""
    s1, s2, s3 = SIGMA
    theta = s1 * (s3 - s2) / ((s2 + s1) * (s3 - s1))
    for i in range(1, 8):
        for j in range(1, 8 - i):
            f = [Fr(i, 8), Fr(j, 8), Fr(8 - i - j, 8)]
            hs_lo, targets = TF.hs_target_fields(f, SIGMA, E0)
            ratio = targets[1][0] / E0[0]  # |E_2| / |E0|
            assert (ratio <= 1) == (hs_lo <= s2)
            assert (hs_lo <= s2) == (f[0] >= 2 * theta * (1 - f[1]))


def test_duality_reduces_upper_side_to_lower_side():
    """HS_hi(f, sigma) inverts to HS_lo(f, 1/sigma); isotropic trees follow."""
    f = [Fr(1, 8), Fr(1, 8), Fr(6, 8)]
    lo, hi = L.hs_bounds(f, SIGMA)
    dlo, dhi = L.hs_bounds(f, [1 / s for s in SIGMA])
    assert dhi == 1 / lo and dlo == 1 / hi
    tree, ph, _ = TC.milton([Fr(3, 8), Fr(2, 8), Fr(3, 8)], SIGMA, 0)
    sig = L.effective(tree, ph)
    dual = L.effective(tree, {k: 1 / v for k, v in ph.items()})
    assert sig[1] == 0 and sig[0] == sig[2]
    assert dual == (1 / sig[0], Fr(0), 1 / sig[0])


def test_required_phase1_variance_closed_form():
    """Var_1 = (c^2 S2 - 1)/f1 matches energy bookkeeping, and is positive."""
    n0 = Fr(2)
    cases = [
        ([Fr(3, 8), Fr(2, 8), Fr(3, 8)], [Fr(1), Fr(2), Fr(5)]),
        ([Fr(3, 4), Fr(1, 8), Fr(1, 8)], [Fr(1), Fr(2), Fr(5)]),
        ([Fr(2, 5), Fr(1, 5), Fr(2, 5)], [Fr(1), Fr(7, 2), Fr(11)]),
    ]
    for f, s in cases:
        hs_lo, tg = TF.hs_target_fields(f, s, E0)
        bracket = sum(f[i] * s[i] * (tg[i][0] ** 2 + tg[i][1] ** 2)
                      for i in range(3)) / n0
        assert TF.required_variance(f, s) == (hs_lo - bracket) / (f[0] * s[0])
        assert TF.required_variance(f, s) > 0


def test_gap_decomposition_is_exact():
    """The mean/variance split of the HS gap reproduces the gap in Q."""
    f = [Fr(3, 8), Fr(2, 8), Fr(3, 8)]
    tree, ph, _ = TC.milton(f, SIGMA, 0)
    d = TF.gap_decomposition(tree, ph, E0)
    assert d["gap"] == 0
    assert d["mean_term"] == 0
    assert d["var_term"] == d["var_required"]
    # a structure that does NOT attain the bound still decomposes exactly
    # (gap_decomposition asserts the identity internally)
    d2 = TF.gap_decomposition(TREE, PHASES, E0)
    assert d2["gap"] > 0


def test_coating_with_comparison_medium_preserves_hs_optimality():
    """Coating an HS-optimal core with extra phase 1 stays HS-optimal, exactly.

    Proof: the coating is a two-phase HS-optimal laminate of core value v0
    (fraction F) in matrix s1, so 1/(value+s1) = F/(v0+s1) + (1-F)/(2 s1).
    The core's own HS identity is 1/(v0+s1) = (1/F) sum_inner f_i/(s_i+s1),
    and the outer phase-1 contributes (1-F)/(2 s1), so the right side is the
    sum over ALL phases, i.e. 1/(HS_lo+s1) at the new fractions.

    Coating scales f2 and f3 equally, so their ratio is preserved while f1
    rises: attainability along a fixed f2:f3 ray is upward closed in f1.
    """
    for base in [[Fr(3, 8), Fr(2, 8), Fr(3, 8)], [Fr(5, 8), Fr(2, 8), Fr(1, 8)]]:
        core, phc, _ = TC.milton(base, SIGMA, 0)
        assert core is not None
        v0 = L.effective(core, phc)[0]
        for Fv in [Fr(9, 10), Fr(1, 2), Fr(1, 4)]:
            tree, rho = TC.iso_coat(core, L.iso(v0), "p1", SIGMA[0], Fv)
            assert tree is not None and rho == Fr(1, 2)
            fr = L.fractions_of(tree)
            names = sorted(fr)
            lo, _ = L.hs_bounds([fr[n] for n in names], [PHASES[n] for n in names])
            val = L.effective(tree, PHASES)
            assert val[1] == 0 and val[0] == val[2]
            assert val[0] == lo
            # the f2:f3 ratio survives the coating
            assert fr["p2"] / fr["p3"] == Fr(base[1], 1) / Fr(base[2], 1)


def test_equal_value_grains_compose_to_an_hs_optimal_mixture():
    """Grains each HS-optimal at their own fractions, sharing a common value v,
    mix to an HS-optimal composite at the aggregate fractions.

    Reciprocal-form proof: each grain has 1/(v+s1) = sum_i g_i/(s_i+s1); the
    aggregate fractions are f_i = sum_g w_g g_i^(g), so
    sum_i f_i/(s_i+s1) = sum_g w_g/(v+s1) = 1/(v+s1). This generalises both the
    Milton assemblage and the coating lemma.
    """
    s1, s2, s3 = SIGMA
    for v in [Fr(3, 2), Fr(7, 4), Fr(19, 10)]:
        def phi_for(si, v=v):
            return (1 / (v + s1) - 1 / (2 * s1)) / (1 / (si + s1) - 1 / (2 * s1))
        p2, p3 = phi_for(s2), phi_for(s3)
        assert 0 < p2 <= 1 and 0 < p3 <= 1
        g2, _ = TC.iso_coat(TC.leaf("p2"), L.iso(s2), "p1", s1, p2)
        g3, _ = TC.iso_coat(TC.leaf("p3"), L.iso(s3), "p1", s1, p3)
        assert L.effective(g2, PHASES) == L.iso(v)
        assert L.effective(g3, PHASES) == L.iso(v)
        for w in [Fr(1, 3), Fr(1, 2), Fr(3, 4)]:
            tree = TC.node(w, (1, 1), g2, g3)
            fr = L.fractions_of(tree)
            names = sorted(fr)
            lo, _ = L.hs_bounds([fr[n] for n in names], [PHASES[n] for n in names])
            val = L.effective(tree, PHASES)
            assert val[1] == 0 and val[0] == val[2]
            assert val[0] == lo


def test_no_phase1_free_structure_attains_the_bound():
    """With f1 = 0 the best {2,3} structure strictly exceeds the sigma1-
    comparison value, so the attainability threshold f1* is strictly positive.
    """
    s1, s2, s3 = SIGMA
    for f2 in [Fr(1, 4), Fr(1, 2), Fr(3, 4), Fr(9, 10)]:
        f3 = 1 - f2
        target = 1 / (f2 / (s2 + s1) + f3 / (s3 + s1)) - s1
        best23 = 1 / (f2 / (2 * s2) + f3 / (s3 + s2)) - s2
        assert best23 > target


def _field_matrix_invariants(tree, phases):
    """(S, D*, D**, V) per leaf, from the two orthogonal applied fields.

    Cherkaev's field variable is Z = grad u for a PAIR of potentials; row i is
    the field for potential i. Common factor 1/sqrt(2) dropped throughout.
    """
    a = TF.leaf_fields(tree, phases, (Fr(1), Fr(0)))
    b = TF.leaf_fields(tree, phases, (Fr(0), Fr(1)))
    out = {}
    for (na, va, Ea), (nb, vb, Eb) in zip(a, b):
        assert na == nb and va == vb
        Z11, Z12, Z21, Z22 = Ea[0], Ea[1], Eb[0], Eb[1]
        out.setdefault(na, []).append(
            (Z11 + Z22, Z11 - Z22, Z12 + Z21, Z12 - Z21))
    return out


def test_attaining_structure_meets_cherkaevs_optimality_conditions():
    """The structure that attains HS_lo satisfies (4.24) and (4.25) exactly:
    V vanishes everywhere, and the field in the MOST conducting phase is a
    single isotropic value. Corroborates the published framework."""
    import tp_attain as A
    _, tree, _ = A.attains(SIGMA, Fr(1, 2), Fr(4, 5))
    inv = _field_matrix_invariants(tree, PHASES)
    for name, rows in inv.items():
        for _, _, _, V in rows:
            assert V == 0, name              # (4.24)
    p3 = {(S, Ds, Dss) for S, Ds, Dss, _ in inv["p3"]}
    assert len(p3) == 1                      # (4.25) constant
    (S, Ds, Dss), = p3
    assert Ds == 0 and Dss == 0              # (4.25) isotropic


def test_below_m11_structure_violates_condition_4_25():
    """The candidate counterexample of attempt 006 satisfies V=0 but its most
    conducting phase carries two distinct, anisotropic field values. This is
    the named localization of the open escalation, not a claim either way."""
    import json
    path = (ROOT / "problems" / "three-phase-conductivity" / "data"
            / "attained" / "below-m11-m1-3_25.json")
    rec = json.loads(path.read_text(encoding="utf-8"))
    phases = {k: Fr(v) for k, v in rec["phases"].items()}
    inv = _field_matrix_invariants(rec["tree"], phases)
    for name, rows in inv.items():
        for _, _, _, V in rows:
            assert V == 0, name              # (4.24) still holds
    p3 = {(S, Ds, Dss) for S, Ds, Dss, _ in inv["p3"]}
    assert len(p3) > 1                       # (4.25) constancy fails
    assert any(Ds != 0 or Dss != 0 for S, Ds, Dss in p3)   # and isotropy fails


def _slack_4_26(tree, phases):
    """min over phases 1,2 of (S - varsigma_N)^2 - D^2. Negative violates (4.26)."""
    inv = _field_matrix_invariants(tree, phases)
    p3 = inv["p3"]
    # varsigma_N is the S value in the most conducting phase; volume-weight it
    # when it is not constant (which is itself a violation of (4.25)).
    a = TF.leaf_fields(tree, phases, (Fr(1), Fr(0)))
    vols = {}
    for name, vol, _ in a:
        vols.setdefault(name, []).append(vol)
    tot = sum(vols["p3"])
    sN = sum(v * row[0] for v, row in zip(vols["p3"], p3)) / tot
    return min((S - sN) ** 2 - (Ds * Ds + Dss * Dss)
               for nm in ("p1", "p2") for S, Ds, Dss, _ in inv[nm])


def test_attaining_structure_makes_4_26_exactly_active():
    """At the attaining structure the constraint is tight: slack exactly 0.

    That is the derivation's own equality case, where the minimiser takes
    D = +/- Sigma_1^(1/2)(S).
    """
    import tp_attain as A
    _, tree, _ = A.attains(SIGMA, Fr(1, 2), Fr(4, 5))
    assert _slack_4_26(tree, PHASES) == 0


def test_structure_below_the_bound_violates_4_26():
    """The candidate counterexample exceeds the cap the bound's estimate uses.

    Section 4.2's coefficient on the D^2 integral is negative, so exceeding the
    cap is exactly what lets a structure land under B2. This pins the mechanism
    of the escalation recorded in attempts 006/008/009.
    """
    import json
    path = (ROOT / "problems" / "three-phase-conductivity" / "data"
            / "attained" / "below-m11-m1-3_25.json")
    rec = json.loads(path.read_text(encoding="utf-8"))
    phases = {k: Fr(v) for k, v in rec["phases"].items()}
    assert _slack_4_26(rec["tree"], phases) < 0


def test_attainment_field_pattern_is_forced_and_makes_4_26_active():
    """Two algebraic identities behind the attaining construction.

    With c = sigma* + s1 and rho_i = c/(s_i+s1), normal-flux continuity across a
    p1|p3 interface forces phase 1's larger component y = s3*rho_3/s1. Then:

      (a) the phase-1 mean rho_1 = c/(2 s1) forces the smaller component
          2*rho_1 - y to equal rho_3 -- an identity, not a fitted value;
      (b) S_1 - varsigma_N = (y + rho_3) - 2*rho_3 = y - rho_3 = D_1, so
          D_1^2 = (S_1 - varsigma_N)^2 EXACTLY.

    (b) explains attempt 009's observation that the attaining structure makes
    Cherkaev's (4.26) exactly active: it is forced, not a coincidence.
    """
    cases = [((Fr(1), Fr(2), Fr(5)), Fr(4)),
             ((Fr(1), Fr(3), Fr(7)), Fr(9, 2)),
             ((Fr(2), Fr(5), Fr(11)), Fr(77, 10)),
             ((Fr(1, 3), Fr(5, 2), Fr(9)), Fr(3))]
    for (s1, s2, s3), c in cases:
        rho3 = c / (s3 + s1)
        rho1 = c / (2 * s1)
        y = s3 * rho3 / s1
        assert 2 * rho1 - y == rho3                      # (a)
        S1, D1, sN = y + rho3, y - rho3, 2 * rho3
        assert D1 * D1 == (S1 - sN) ** 2                 # (b)


def test_derived_field_pattern_matches_the_built_structure():
    """The forced pattern is what tp_attain actually produces."""
    import tp_attain as A
    for sig, r in [((Fr(1), Fr(2), Fr(5)), Fr(1, 2)),
                   ((Fr(1), Fr(3), Fr(7)), Fr(1, 2)),
                   ((Fr(1), Fr(4), Fr(9)), Fr(1, 3))]:
        ok, tree, hs = A.attains(sig, r, (1 - r + 1) / 2)
        assert ok
        phases = {"p1": sig[0], "p2": sig[1], "p3": sig[2]}
        c = hs + sig[0]
        rho2, rho3 = c / (sig[1] + sig[0]), c / (sig[2] + sig[0])
        y = sig[2] * rho3 / sig[0]
        a = TF.leaf_fields(tree, phases, (Fr(1), Fr(0)))
        b = TF.leaf_fields(tree, phases, (Fr(0), Fr(1)))
        seen = {}
        for (na, _, Ea), (nb, _, Eb) in zip(a, b):
            seen.setdefault(na, set()).add((Ea[0], Eb[1]))
        assert seen["p2"] == {(rho2, rho2)}
        assert seen["p3"] == {(rho3, rho3)}
        assert seen["p1"] == {(rho3, y), (y, rho3)}
