"""Tests for harness/giuga.

Fast and deterministic. The point of each test is that it fails if one of the
three formulations of Giuga's congruence, or one of the pruning rules in the
two enumerators, is broken -- those are the places where a wrong number could
be produced without leaving a trace in the output.
"""

import sys
from fractions import Fraction
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "harness"))

from giuga.admissible import Pool, leaf_bound_exact, search  # noqa: E402
from giuga.conditions import (  # noqa: E402
    agoh_holds, bernoulli, counterexample_conditions, factorize,
    giuga_sum_condition, is_carmichael, is_counterexample, is_giuga_number, is_prime,
    power_sum_mod, power_sum_mod_structural, satisfies_congruence,
    satisfies_congruence_structural,
)
from giuga.enumerate_giuga import (  # noqa: E402
    enumerate_by_prime_sets, product, scan_giuga_numbers, verify_set,
)

# The four Giuga numbers below 10^6 and the full lists by factor count that the
# prime-set recursion is expected to reproduce.
GIUGA_BY_FACTORS = {
    3: [(2, 3, 5)],
    4: [(2, 3, 7, 41), (2, 3, 11, 13)],
    5: [(2, 3, 11, 17, 59)],
    6: [(2, 3, 7, 43, 3041, 4447), (2, 3, 11, 23, 31, 47057)],
}


def test_power_sum_agrees_with_structural_residue():
    for n in range(2, 400):
        s = power_sum_mod_structural(n)
        if s is not None:
            assert power_sum_mod(n) == s, n


def test_congruence_characterises_primality():
    for n in range(2, 400):
        assert satisfies_congruence(n) is is_prime(n), n
        assert satisfies_congruence_structural(n) is is_prime(n), n


def test_agoh_agrees_with_primality():
    b = bernoulli(80)
    for n in range(2, 81):
        assert agoh_holds(n, b) is is_prime(n), n


def test_no_even_n_passes_agoh():
    b = bernoulli(80)
    assert [n for n in range(4, 81, 2) if agoh_holds(n, b)] == []


def test_counterexample_criterion_matches_the_congruence():
    for n in range(2, 400):
        composite = not is_prime(n)
        assert is_counterexample(n) is (composite and satisfies_congruence(n)), n


def test_sieve_finds_the_small_giuga_numbers():
    assert scan_giuga_numbers(100_000) == [30, 858, 1722, 66198]


def test_prime_set_enumeration_matches_the_sieve():
    for m, expected in GIUGA_BY_FACTORS.items():
        got = sorted(enumerate_by_prime_sets(m))
        assert got == sorted(expected), (m, got)
        for s in got:
            assert verify_set(s)                       # pointwise divisibility
            assert giuga_sum_condition(s)              # rational form
            assert is_giuga_number(product(s))


def test_disabling_a_pruning_rule_does_not_change_the_answer():
    # An unsound prune would silently delete solutions; toggling must only
    # change the cost.
    for m in (4, 5, 6):
        base = sorted(enumerate_by_prime_sets(m))
        assert sorted(enumerate_by_prime_sets(m, use_p4=False)) == base, m
        assert sorted(enumerate_by_prime_sets(m, use_p7=False)) == base, m


def test_known_carmichael_numbers():
    assert is_carmichael(561) and is_carmichael(1105) and is_carmichael(1729)
    assert not is_carmichael(560)
    # 561 is Carmichael but not Giuga, so it is not a counterexample
    assert not is_giuga_number(561)
    assert not is_counterexample(561)


def test_carmichael_numbers_are_odd_and_obey_the_pairwise_rule():
    # If p, q both divide a Carmichael number then q does not divide p-1, and
    # no Carmichael number is even. These two constraints are what the
    # admissible-set bound is built on, so a counterexample here would
    # invalidate it.
    found = 0
    for n in range(3, 100_000):
        if not is_carmichael(n):
            continue
        found += 1
        assert n % 2 == 1, n
        fs = sorted(factorize(n))
        for p in fs:
            for q in fs:
                if p != q:
                    assert (p - 1) % q != 0, (n, p, q)
    assert found == 16          # Carmichael numbers below 100000


def test_giuga_numbers_have_reciprocal_sum_above_one():
    for n in scan_giuga_numbers(100_000):
        ps = sorted(factorize(n))
        assert sum((Fraction(1, p) for p in ps), Fraction(0)) > 1, n


def test_counterexample_conditions_pointwise_form():
    # 30 is a Giuga number but not Carmichael, so it fails the joint condition.
    assert not counterexample_conditions([2, 3, 5])
    assert is_giuga_number(30)
    assert not is_carmichael(30)


def test_admissible_bound_is_stable_under_pruning():
    for X, expected in ((13, 202), (23, 554)):
        pool = Pool(200_000, X)
        pruned, witness, _ = search(pool, best=3000, prune=True)
        unpruned, _, _ = search(pool, best=3000, prune=False)
        assert pruned == unpruned == expected, X
        # third opinion: exact Fractions, no scaling, no bitmasks
        assert leaf_bound_exact(witness[0], X, 200_000) == expected, X
