"""Tests for the assumed-symmetry construction machinery.

The orbit-matrix conditions and the lifter are only worth running on an open
case if they demonstrably recover known ones, so most of these are calibration
tests against graphs that exist.

Fast and deterministic: fixed seeds, standard library only.
"""

import os
import random
import sys
from itertools import combinations

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, os.path.join(ROOT, "problems", "conway-99", "explore"))

import constructions as C  # noqa: E402
import orbit_general as OG  # noqa: E402
import orbit_lift as OL  # noqa: E402
import orbit_matrix as OM  # noqa: E402
import srg  # noqa: E402
import srg_anneal as SA  # noqa: E402
from orbit_calibrate import orbits_of, true_orbit_matrix  # noqa: E402


# --------------------------------------------------------------------------
# the orbit-matrix conditions
# --------------------------------------------------------------------------

def test_real_graphs_satisfy_the_orbit_conditions():
    """A genuine semiregular automorphism must produce a legal orbit matrix."""
    cases = [
        ("paley9", C.paley9(), [(v + 3) % 9 for v in range(9)], 3, (9, 4, 1, 2)),
        ("rook3", C.rook(3),
         [((v // 3 + 1) % 3) * 3 + v % 3 for v in range(9)], 3, (9, 4, 1, 2)),
        ("rook4", C.rook(4),
         [((v // 4 + 1) % 4) * 4 + v % 4 for v in range(16)], 4, (16, 6, 2, 2)),
    ]
    for name, adj, perm, m, (n, k, lam, mu) in cases:
        orbs = orbits_of(perm, n)
        assert sorted({len(o) for o in orbs}) == [m], name
        R = true_orbit_matrix(adj, orbs)
        ok, why = OM.check_orbit_matrix(R, [m] * len(orbs), k, lam, mu)
        assert ok, f"{name}: {why}"


def test_enumeration_contains_the_real_orbit_matrix():
    adj = C.paley9()
    orbs = orbits_of([(v + 3) % 9 for v in range(9)], 9)
    R = true_orbit_matrix(adj, orbs)
    sols, _ = OM.enumerate_semiregular(3, 4, 1, 2, 9)
    assert any(OM.canonical_form(R) == OM.canonical_form(S) for S in sols)


def test_isomorph_rejection_preserves_every_class():
    """The sorted-diagonal filter is a canonical-form restriction, so it must
    not lose a single isomorphism class."""
    for (m, n, k, lam, mu) in [(3, 9, 4, 1, 2), (4, 16, 6, 2, 2),
                               (5, 10, 3, 0, 1), (33, 99, 14, 1, 2)]:
        raw, _ = OM.enumerate_semiregular(m, k, lam, mu, n)
        cut, _ = OM.enumerate_semiregular(m, k, lam, mu, n, sorted_diagonal=True)
        assert ({OM.canonical_form(R) for R in raw}
                == {OM.canonical_form(R) for R in cut}), (m, n)


def test_general_enumerator_matches_the_semiregular_one():
    for (m, n, k, lam, mu) in [(3, 9, 4, 1, 2), (5, 10, 3, 0, 1),
                               (33, 99, 14, 1, 2)]:
        t = n // m
        a, _ = OM.enumerate_semiregular(m, k, lam, mu, n)
        b, _ = OG.enumerate_general([m] * t, k, lam, mu)
        assert ({tuple(map(tuple, R)) for R in a}
                == {tuple(map(tuple, R)) for R in b}), (m, n)


# --------------------------------------------------------------------------
# the Z_33 case for (99,14,1,2)
# --------------------------------------------------------------------------

def test_z33_has_exactly_one_orbit_matrix_bruteforce():
    """Independent re-derivation of the m=33 enumeration.

    Brute force over every symmetric 3x3 candidate with even diagonal and row
    sums 14, checking (Q) directly, rather than reusing the backtracking
    search.
    """
    m, k, lam, mu = 33, 14, 1, 2
    sizes = [m] * 3
    found = []
    for r00 in range(0, m + 1, 2):
        for r01 in range(0, k + 1):
            r02 = k - r00 - r01
            if not (0 <= r02 <= m):
                continue
            for r11 in range(0, m + 1, 2):
                r12 = k - r01 - r11
                if not (0 <= r12 <= m):
                    continue
                r22 = k - r02 - r12
                if not (0 <= r22 <= m) or r22 % 2:
                    continue
                R = [[r00, r01, r02], [r01, r11, r12], [r02, r12, r22]]
                ok, _ = OM.check_orbit_matrix(R, sizes, k, lam, mu)
                if ok:
                    found.append(R)
    assert len(found) == 1, found
    assert found[0] == [[2, 6, 6], [6, 2, 6], [6, 6, 2]]

    sols, _ = OM.enumerate_semiregular(m, k, lam, mu, 99)
    assert sols == found


def test_z33_forced_autocorrelations_kill_all_but_46_cases():
    import m33_exhaustive as X
    diag = [frozenset({a, (-a) % 33}) for a in range(1, 17)]
    alive = 0
    for S00 in diag:
        f0 = X.diag_target(S00)
        for S11 in diag:
            f1 = X.diag_target(S11)
            for S22 in diag:
                f2 = X.diag_target(S22)
                sol, _ = X.solve_offdiag(f0, f1, f2)
                if sol is not None:
                    alive += 1
    assert alive == 46


def test_autocorrelation_is_translation_invariant():
    import m33_exhaustive as X
    U = [0, 1, 4, 9, 11, 20]
    base = X.autocorr(U)
    for c in range(33):
        assert X.autocorr([(x + c) % 33 for x in U]) == base


# --------------------------------------------------------------------------
# the lifter and the annealer must rebuild graphs that exist
# --------------------------------------------------------------------------

@pytest.mark.parametrize("name,m,params", [
    ("paley9", 3, (9, 4, 1, 2)),
    ("petersen", 5, (10, 3, 0, 1)),
])
def test_lifter_rebuilds_a_graph_that_exists(name, m, params):
    n, k, lam, mu = params
    if name == "paley9":
        adj, perm = C.paley9(), [(v + 3) % 9 for v in range(9)]
    else:
        adj = C.petersen()
        vs = list(combinations(range(5), 2))
        idx = {v: i for i, v in enumerate(vs)}
        perm = [idx[tuple(sorted(((a + 1) % 5, (b + 1) % 5)))] for (a, b) in vs]
    orbs = orbits_of(perm, n)
    R = true_orbit_matrix(adj, orbs)
    res = OL.anneal(R, m, len(orbs), k, lam, mu, seed=3, iters=20000, restarts=6)
    assert res["result"] == "solved", res.get("best_energy")
    ok, got, _ = OL.verify_witness(res["S"], m, len(orbs), n, k, lam, mu)
    assert ok and got == params


def test_annealer_incremental_energy_matches_full_recompute():
    rng = random.Random(11)
    S = SA.Search(20, 6, 1, 2, rng)
    assert S.energy == S._full_energy()
    for _ in range(120):
        tok = S.swap(rng.randrange(len(S.edges)), rng.randrange(len(S.edges)))
        if tok is None:
            continue
        assert S.energy == S._full_energy()


def test_annealer_starts_from_a_regular_graph():
    rng = random.Random(5)
    for (n, k) in [(9, 4), (16, 6), (99, 14), (27, 10)]:
        S = SA.Search(n, k, 1, 2, rng)
        assert all(bin(S.adj[v]).count("1") == k for v in range(n)), (n, k)


@pytest.mark.parametrize("params", [(9, 4, 1, 2), (15, 6, 1, 3)])
def test_annealer_finds_graphs_that_exist(params):
    n, k, lam, mu = params
    res = SA.anneal(n, k, lam, mu, seed=2, iters=30000, restarts=8)
    assert res["result"] == "solved", res.get("best_energy")
    ok, got = srg.check_srg(res["adj"])
    assert ok and got == params


# --------------------------------------------------------------------------
# L5: order-3 automorphisms, and the order-33 elimination
# --------------------------------------------------------------------------

def test_order3_counting_identity_on_real_graphs():
    """The identity sum_{x in N_F(u)} deg_F(x) = deg_F(u)(1+lam-mu) + mu(f-1)
    must hold for the fixed subgraph of every order-3 automorphism."""
    import order3_lemma as L5
    from autos import automorphisms, order_of
    cases = [(C.paley9(), (9, 4, 1, 2)), (C.petersen(), (10, 3, 0, 1)),
             (C.rook(4), (16, 6, 2, 2))]
    for adj, (n, k, lam, mu) in cases:
        for p in automorphisms(adj):
            if order_of(p) != 3:
                continue
            F = [v for v in range(n) if p[v] == v]
            if not F:
                continue
            assert L5.closed_under_common_neighbours(adj, F, n), (n, k)
            ok, why = L5.counting_identity_holds(adj, F, lam, mu)
            assert ok, (n, k, why)


def test_order3_fixed_points_are_zero_for_k4():
    """L5's analogue at k=4: with only two matching edges, an order-3 map must
    fix both, hence fix N(v) pointwise, hence be the identity.  So a genuine
    order-3 automorphism of a (9,4,1,2) graph has no fixed points at all."""
    from autos import automorphisms, order_of
    for adj in (C.paley9(), C.rook(3)):
        for p in automorphisms(adj):
            if order_of(p) == 3:
                assert not [v for v in range(9) if p[v] == v]


def test_srg_33_8_1_2_is_infeasible():
    """The step that kills the 8-regular branch of L5."""
    ok, rep = srg.feasibility(33, 8, 1, 2)
    assert not ok
    assert "integral-spectrum" in srg.failed_conditions(rep)


def test_order33_profiles_are_exactly_three_and_two_die_immediately():
    """99 = 33a + 11b with a >= 1; the b > 0 profiles are killed by L5 because
    the 11-orbit points are exactly Fix(sigma^11), which must number 0 or 3."""
    profs = OG.z33_profiles()
    assert len(profs) == 3
    for name, sizes in profs.items():
        eleven = sum(1 for s in sizes if s == 11) * 11
        if eleven:
            assert eleven not in (0, 3), (name, eleven)   # so L5 forbids it


def test_order33_1x33_6x11_has_no_orbit_matrix():
    """Independent of L5: this profile dies at the orbit-matrix level."""
    sizes = [33] + [11] * 6
    sols, _ = OG.enumerate_general(sizes, 14, 1, 2)
    assert sols == []


def test_mixed_size_enumerator_finds_real_orbit_matrices():
    """Calibration for the mixed-size path, without which a count of zero
    above would prove nothing."""
    from orbit_calibrate import orbits_of, true_orbit_matrix
    # Paley(9) under x -> -x : orbit sizes [1,2,2,2,2]
    p9 = C.paley9()
    perm = []
    for v in range(9):
        a, b = divmod(v, 3)
        perm.append(((-a) % 3) * 3 + ((-b) % 3))
    orbs = orbits_of(perm, 9)
    sizes = [len(o) for o in orbs]
    R = true_orbit_matrix(p9, orbs)
    ok, why = OM.check_orbit_matrix(R, sizes, 4, 1, 2)
    assert ok, why
    sols, _ = OG.enumerate_general(sizes, 4, 1, 2)
    assert any(all(R[i][j] == S[i][j] for i in range(len(sizes))
                   for j in range(len(sizes))) for S in sols)
