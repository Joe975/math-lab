"""Calibration for problems/conway-99/explore/orbit7.py (attempt 003).

The enumerator's claims about the order-7 profile of the 99-graph are only
worth anything if (a) it provably contains the real orbit matrices of graphs
that exist, including ones with a fixed point, (b) its counts agree with the
independent enumerator from attempt 002 (orbit_general.py) when symmetry
breaking is off, (c) its symmetry breaking loses no isomorphism class, and
(d) the character-block multiplicity lemma it uses holds on real graphs.
"""

import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "harness" / "conway-99"))
sys.path.insert(0, str(ROOT / "problems" / "conway-99" / "explore"))

import constructions as C  # noqa: E402
import orbit7  # noqa: E402
import orbit_general  # noqa: E402
import char_blocks  # noqa: E402


def real_orbit_matrix(adj, perm):
    orbs = char_blocks.orbits_of(perm, len(adj))
    assert char_blocks.is_automorphism(adj, perm)
    return char_blocks.orbit_matrix(adj, orbs), [len(o) for o in orbs]


REAL = {
    "paley9/Z3": (C.paley9(), [(v + 3) % 9 for v in range(9)], (4, 1, 2)),
    "paley9/Z2": (C.paley9(), [0, 2, 1, 6, 8, 7, 3, 5, 4], (4, 1, 2)),
    "rook4/Z3": (char_blocks.rook(4),
                 [{0: 1, 1: 2, 2: 0, 3: 3}[v // 4] * 4 + {0: 1, 1: 2, 2: 0, 3: 3}[v % 4]
                  for v in range(16)], (6, 2, 2)),
}
_adj, _vs, _idx = char_blocks.petersen()
REAL["petersen/Z3"] = (
    _adj,
    [_idx[tuple(sorted({0: 1, 1: 2, 2: 0}.get(x, x) for x in v))] for v in _vs],
    (3, 0, 1),
)


@pytest.mark.parametrize("name", sorted(REAL))
def test_enumeration_contains_the_real_orbit_matrix(name):
    adj, perm, (k, lam, mu) = REAL[name]
    R, sizes = real_orbit_matrix(adj, perm)
    ok, why = orbit7.check_orbit_matrix(R, sizes, k, lam, mu)
    assert ok, why
    # sizes must be in the enumerator's fixed-point-first order; permute
    order = sorted(range(len(sizes)), key=lambda i: sizes[i])
    sizes_o = [sizes[i] for i in order]
    R_o = [[R[a][b] for b in order] for a in order]
    target = orbit7.canonical(R_o, sizes_o)
    found, stats = orbit7.enumerate_orbit_matrices(sizes_o, k, lam, mu, symmetry=True)
    assert stats["exhausted"]
    assert target in {orbit7.canonical(M, sizes_o) for M in found}


@pytest.mark.parametrize("sizes,k,lam,mu", [
    ([3, 3, 3], 4, 1, 2),
    ([3, 3, 3, 1], 3, 0, 1),
    ([1, 3, 3, 3, 3, 3], 6, 2, 2),
    ([1, 2, 2, 2, 2], 4, 1, 2),
    ([33, 33, 11, 11, 11], 14, 1, 2),
])
def test_counts_match_the_independent_enumerator_from_002(sizes, k, lam, mu):
    a, _ = orbit7.enumerate_orbit_matrices(sizes, k, lam, mu, symmetry=False)
    b, _ = orbit_general.enumerate_general(sizes, k, lam, mu)
    assert len(a) == len(b)
    assert {orbit7.canonical(M, sizes) for M in a} == {orbit7.canonical(M, sizes) for M in b}


@pytest.mark.parametrize("sizes,k,lam,mu", [
    ([1, 3, 3, 3, 3, 3], 6, 2, 2),
    ([1, 2, 2, 2, 2], 4, 1, 2),
    ([33, 33, 11, 11, 11], 14, 1, 2),
    ([2, 2, 2, 2, 2, 2], 3, 0, 1),
])
def test_symmetry_breaking_loses_no_class(sizes, k, lam, mu):
    with_sym, s1 = orbit7.enumerate_orbit_matrices(sizes, k, lam, mu, symmetry=True)
    without, s2 = orbit7.enumerate_orbit_matrices(sizes, k, lam, mu, symmetry=False)
    assert s1["exhausted"] and s2["exhausted"]
    c1 = {orbit7.canonical(M, sizes) for M in with_sym}
    c2 = {orbit7.canonical(M, sizes) for M in without}
    assert c1 == c2
    assert len(with_sym) <= len(without)


def test_psd_rank():
    assert orbit7.psd_rank([[2, 1], [1, 2]]) == (True, 2)
    assert orbit7.psd_rank([[1, 1], [1, 1]]) == (True, 1)
    assert orbit7.psd_rank([[1, 2], [2, 1]]) == (False, 1)
    assert orbit7.psd_rank([[0, 0], [0, 3]]) == (True, 1)
    assert orbit7.psd_rank([[0, 1], [1, 3]]) == (False, 0)


@pytest.mark.parametrize("name", ["rook4/Z3", "petersen/Z3", "paley9/Z3"])
def test_spectral_prune_keeps_the_real_quotient(name):
    """With `spectral` set to the real quotient's multiplicity of the larger
    restricted eigenvalue, the real orbit matrix must survive."""
    adj, perm, (k, lam, mu) = REAL[name]
    R, sizes = real_orbit_matrix(adj, perm)
    order = sorted(range(len(sizes)), key=lambda i: sizes[i])
    sizes_o = [sizes[i] for i in order]
    R_o = [[R[a][b] for b in order] for a in order]
    disc = (lam - mu) ** 2 + 4 * (k - mu)
    r = (lam - mu + int(round(disc ** 0.5))) // 2
    t = len(sizes_o)
    a = t - char_blocks.rank_q([[R_o[i][j] - (r if i == j else 0) for j in range(t)]
                                for i in range(t)])
    found, stats = orbit7.enumerate_orbit_matrices(sizes_o, k, lam, mu, spectral=a)
    assert stats["exhausted"]
    assert orbit7.canonical(R_o, sizes_o) in {orbit7.canonical(M, sizes_o) for M in found}
    # and the wrong multiplicity must cut it
    found_wrong, _ = orbit7.enumerate_orbit_matrices(sizes_o, k, lam, mu, spectral=a + 1
                                                     if a + 1 < t - 1 else a - 1)
    assert orbit7.canonical(R_o, sizes_o) not in {orbit7.canonical(M, sizes_o)
                                                  for M in found_wrong}


def test_character_block_multiplicity_lemma_on_real_graphs(capsys):
    assert char_blocks.main() == 0
    assert "MISMATCH" not in capsys.readouterr().out


def test_bvls_order11_control_seed_and_deep_completion():
    """BvLS(243) under its order-11 shift has the same shape as the order-7
    case at 99 (one fixed point, N(v0) two orbits joined by the matching).
    Its real orbit matrix must carry the forced seed, and the enumerator
    seeded with its first 14 rows must complete to exactly that matrix."""
    import bvls_control  # noqa: E402  (slow import: builds the 243-graph)
    adj, perm, orbs, R = bvls_control.build()
    sizes = [len(o) for o in orbs]
    ok, why = orbit7.check_orbit_matrix(R, sizes, 22, 1, 2)
    assert ok, why
    for (i, j), v in bvls_control.seed_for(len(sizes)).items():
        assert R[i][j] == v, (i, j, R[i][j], v)
    seed = {(i, j): R[i][j] for i in range(14) for j in range(len(sizes))}
    found, stats = orbit7.enumerate_orbit_matrices(sizes, 22, 1, 2, seed=seed,
                                                   symmetry=False)
    assert stats["exhausted"]
    assert found == [R]


def test_order7_profile_has_no_orbit_matrix_within_budget():
    """The headline computation, budget-capped for the suite: no solution in
    the first 2M nodes and the row-completion profile is the recorded one.
    Set MATHLAB_SLOW=1 to run it to exhaustion (about 30 s)."""
    sizes, seed = orbit7.z7_profile()
    budget = None if os.environ.get("MATHLAB_SLOW") else 2_000_000
    found, stats = orbit7.enumerate_orbit_matrices(
        sizes, 14, 1, 2, seed=seed, node_budget=budget, dump_rows=3)
    assert found == []
    if budget is None:
        assert stats["exhausted"]
        assert stats["row_completions"][3] == 33
        assert stats["row_completions"][4] == 414
        # rows 1 and 2 (the two orbits in N(v0)) have exactly one completion
        # up to relabelling of orbits 3..14: (2,2,2,1^6,0^3) and its complement
        dumped = stats["dumped"]
        assert len(dumped) == 1
        assert dumped[0][1] == [1, 0, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 0, 0, 0]
        assert dumped[0][2] == [1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2]
