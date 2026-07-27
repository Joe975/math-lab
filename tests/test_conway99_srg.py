"""Calibration tests for the strongly-regular-graph harness.

These are the checks that make the harness evidence rather than assertion: the
feasibility calculator has to recover published verdicts in *both* directions,
and the local reduction has to hold on graphs that actually exist.

Fast and deterministic; no randomness, no network, standard library only.
"""

import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))

import constructions as C  # noqa: E402
import local_model as L  # noqa: E402
import srg  # noqa: E402


# --------------------------------------------------------------------------
# exact feasibility
# --------------------------------------------------------------------------

REALISABLE = [
    (9, 4, 1, 2),      # Paley(9)
    (10, 3, 0, 1),     # Petersen
    (16, 6, 2, 2),     # rook(4) and Shrikhande
    (16, 5, 0, 2),     # Clebsch
    (50, 7, 0, 1),     # Hoffman-Singleton
    (243, 22, 1, 2),   # Berlekamp-van Lint-Seidel
    (99, 14, 1, 2),    # open, and feasible -- that is why it is open
]

NOT_REALISABLE = [
    ((28, 9, 0, 4), "absolute-bound-g"),
    ((33, 8, 1, 2), "integral-spectrum"),
    ((21, 10, 4, 5), "conference-sum-two-squares"),
]


@pytest.mark.parametrize("params", REALISABLE)
def test_feasible_sets_are_admitted(params):
    ok, report = srg.feasibility(*params)
    assert ok, f"{params} wrongly rejected: {srg.failed_conditions(report)}"


@pytest.mark.parametrize("params,why", NOT_REALISABLE)
def test_infeasible_sets_are_rejected_for_the_right_reason(params, why):
    ok, report = srg.feasibility(*params)
    assert not ok, f"{params} wrongly admitted"
    assert why in srg.failed_conditions(report), (
        f"{params} rejected, but not by {why}: "
        f"{srg.failed_conditions(report)}")


def test_target_spectrum_is_exact_and_satisfies_the_trace_identities():
    """(99,14,1,2) has eigenvalues 3 and -4 with multiplicities 54 and 44.

    Pinned by tr(A) = 0 and tr(A^2) = nk, both checked here in integers.  A
    multiplicity pair that misses these is not the spectrum of any graph.
    """
    s = srg.spectrum(99, 14, 1, 2)
    assert s["ok"] and not s["conference"]
    assert (s["r"], s["s"]) == (3, -4)
    assert (s["f"], s["g"]) == (54, 44)
    n, k, r, sv, f, g = 99, 14, s["r"], s["s"], s["f"], s["g"]
    assert f + g == n - 1
    assert k + f * r + g * sv == 0             # tr(A) = 0
    assert k * k + f * r * r + g * sv * sv == n * k   # tr(A^2) = nk


def test_lambda1_mu2_family_is_exactly_five_sets_up_to_2e6():
    fam = [(n, k) for (n, k, _, _, _) in srg.enumerate_feasible(1, 2, 2_000_000)]
    assert fam == [(9, 4), (99, 14), (243, 22), (6273, 112), (494019, 994)]


def test_sum_of_two_squares_is_exact():
    assert srg.is_sum_of_two_squares(9)
    assert srg.is_sum_of_two_squares(25)
    assert srg.is_sum_of_two_squares(45)
    assert not srg.is_sum_of_two_squares(21)
    assert not srg.is_sum_of_two_squares(3)


# --------------------------------------------------------------------------
# explicit graphs
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name", sorted(C.CATALOGUE))
def test_constructions_have_their_claimed_parameters(name):
    build, params = C.CATALOGUE[name]
    ok, got = srg.check_srg(build(), *params)
    assert ok, f"{name}: {got}"


def test_shrikhande_and_rook4_are_not_isomorphic():
    """Same parameters, different graphs -- the standard sanity pair.

    Distinguished here without any isomorphism test, by a local invariant:
    the rook graph's neighbourhoods induce 2K_3, Shrikhande's induce C_6.
    """
    assert L_comp(C.rook(4)) != L_comp(C.shrikhande())


def L_comp(adj):
    return sorted(srg.neighbourhood_components(adj, v) for v in range(len(adj)))


def test_ternary_golay_code_is_perfect():
    code = C.ternary_golay_code()
    assert len(code) == 3 ** 6
    weights = sorted({sum(1 for c in w if c) for w in code})
    assert weights[0] == 0 and weights[1] == 5, weights


# --------------------------------------------------------------------------
# the local reduction
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name,builder,k", [
    ("paley9", C.paley9, 4),
    ("rook3", lambda: C.rook(3), 4),
])
def test_local_reduction_holds_at_every_vertex(name, builder, k):
    adj = builder()
    for v0 in range(len(adj)):
        ok, dec = L.verify_decomposition(adj, v0)
        assert ok, f"{name} at {v0}: {dec}"
        assert len(dec["far"]) == k * (k - 1) // 2 - k // 2
        assert len(dec["vert_of"]) == len(dec["far"])


def test_local_reduction_holds_on_bvls243():
    """The large realised member of the lam=1, mu=2 family.

    Checked at a few base vertices rather than all 243 to keep the suite fast;
    problems/conway-99/explore/make_evidence.py does all of them.
    """
    adj, _ = C.bvls243()
    ok, params = srg.check_srg(adj)
    assert ok and params == (243, 22, 1, 2)
    for v0 in (0, 1, 7, 100, 242):
        ok, dec = L.verify_decomposition(adj, v0)
        assert ok, f"bvls243 at {v0}: {dec}"
        assert len(dec["far"]) == 220
        assert len(dec["vert_of"]) == 220
        assert all(len(m) == 20 for m in dec["match"].values())


def test_lambda_one_forces_a_perfect_matching_neighbourhood():
    for builder in (C.paley9, lambda: C.rook(3)):
        adj = builder()
        for v in range(len(adj)):
            assert srg.neighbourhood_components(adj, v) == [2, 2]


def test_decompose_rejects_a_graph_outside_the_family():
    """Petersen has lambda = 0, so the reduction must refuse it rather than
    silently returning something."""
    with pytest.raises(ValueError):
        L.decompose(C.petersen(), 0)
