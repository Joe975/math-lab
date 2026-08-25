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


def test_factorint_round_trips_and_reports_failure_honestly():
    from giuga.factorint import divisors, factorize_full
    for n in (1, 2, 97, 2 * 3 * 5 * 7, 10 ** 12 + 39, 2 ** 20 * 3 ** 7,
              999_999_000_001, 1234567891011121314151617):
        fac, cof = factorize_full(n)
        prod = cof
        for p, e in fac.items():
            assert is_prime(p), (n, p)
            prod *= p ** e
        assert prod == n, n
        if cof == 1:
            ds = sorted(divisors(fac))
            assert ds[0] == 1 and ds[-1] == n
            assert all(n % d == 0 for d in ds)
            assert len(ds) == len(set(ds))
    # a hard semiprime with a tiny budget must report the cofactor, not lie
    hard = 32416190071 * 32416189381
    fac, cof = factorize_full(hard, small_limit=1000, budget=50)
    assert cof != 1 and cof > 1
    prod = cof
    for p, e in fac.items():
        prod *= p ** e
    assert prod == hard


def test_p8_reproduces_p7_on_every_known_giuga_number():
    # The scan and the factorisation closed form are independently written and
    # must agree; disagreement means one of the two derivations is wrong.
    for m in GIUGA_BY_FACTORS:
        unresolved = []
        got = sorted(enumerate_by_prime_sets(m, use_p8=True, unresolved=unresolved))
        assert unresolved == []
        assert got == sorted(GIUGA_BY_FACTORS[m]), m


def test_rectangle_identity_holds_on_known_giuga_numbers():
    # (a*p - N)(a*q - N) == N^2 - a, with rem = a/N in lowest terms.
    for m, sets in GIUGA_BY_FACTORS.items():
        for s in sets:
            head, p, q = list(s[:-2]), s[-2], s[-1]
            N = product(head)
            rem = 1 - sum((Fraction(1, x) for x in head), Fraction(0))
            assert rem.denominator == N, (s, rem)   # the denominator lemma
            a = rem.numerator
            assert (a * p - N) * (a * q - N) == N * N - a, s


def test_giuga_sequences_contain_exactly_the_giuga_numbers():
    from giuga.sequences import (
        all_prime, enumerate_sequences, verify_sequence,
    )
    for m, expected in GIUGA_BY_FACTORS.items():
        unresolved = []
        seqs = sorted(enumerate_sequences(m, unresolved=unresolved))
        assert unresolved == []
        for s in seqs:
            assert verify_sequence(s), s
        assert sorted(s for s in seqs if all_prime(s)) == sorted(expected), m
    # length 5 has a genuinely non-prime member, so the wider search is not
    # just the prime search in disguise
    seqs5 = sorted(enumerate_sequences(5))
    assert (2, 3, 7, 83, 85) in seqs5
    assert not all_prime((2, 3, 7, 83, 85))


def test_no_short_odd_giuga_sequence():
    from giuga.enumerate_giuga import k_upper_bound as prime_k_upper_bound
    from giuga.sequences import enumerate_sequences
    # The eight smallest odd PRIMES have reciprocal sum 0.998956 < 1, so no odd
    # Giuga number has fewer than 9 prime factors -- no search needed.
    assert prime_k_upper_bound(8, 3) == 0
    assert prime_k_upper_bound(9, 3) == 1
    # Odd Giuga *sequences* may use composite odd entries, so the sum reaches 1
    # sooner (3,5,7,9,11,13,15 already exceeds it); those still have to be
    # searched, and there are none up to length 9.
    for m in range(3, 10):
        assert list(enumerate_sequences(m, min_entry=3, odd_only=True)) == [], m


def test_parity_relation_between_factor_count_and_k():
    from giuga.conditions import giuga_parity_ok
    # every known Giuga number has k = 1, checked directly from the definition
    for m, sets in GIUGA_BY_FACTORS.items():
        for s in sets:
            n = product(s)
            assert sum(n // p for p in s) - 1 == n, s      # k == 1
    # for odd n and k = 1 the factor count must be even
    assert giuga_parity_ok(10, 1, odd=True)
    assert not giuga_parity_ok(9, 1, odd=True)
    assert not giuga_parity_ok(11, 1, odd=True)
    assert giuga_parity_ok(9, 2, odd=True)
    # the relation makes no claim about even n
    for m in range(3, 9):
        assert giuga_parity_ok(m, 1, odd=False)
    # and it is consistent with the searches: the odd prime-set enumeration
    # returns nothing at odd m, which the parity relation predicts a priori
    for m in (9, 11):
        assert list(enumerate_by_prime_sets(m, min_prime=3)) == [], m


def test_parity_identity_on_random_odd_squarefree_integers():
    # The load-bearing arithmetic behind the parity law: n odd and squarefree
    # with m prime factors implies sum_{p|n} n/p == m (mod 2). The law's
    # conclusion cannot be checked against examples (no odd Giuga number is
    # known), but this step can.
    import random
    from giuga.enumerate_giuga import PRIMES
    rng = random.Random(7)
    odd_primes = [p for p in PRIMES.upto(500) if p > 2]
    for _ in range(500):
        m = rng.randint(2, 7)
        fs = rng.sample(odd_primes, m)
        n = product(fs)
        assert sum(n // p for p in fs) % 2 == m % 2, (n, fs)


def test_a1_prefixes_agree_between_sieve_and_prime_set_recursion():
    # a = 1 prefixes: squarefree N with sum_{p|N} 1/p + 1/N = 1.
    import importlib.util
    path = (Path(__file__).resolve().parent.parent
            / "problems" / "giuga" / "explore" / "a1_prefixes.py")
    spec = importlib.util.spec_from_file_location("a1_prefixes", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    limit = 100_000
    by_recursion = mod.search(limit)

    by_sieve = []
    for n in range(2, limit + 1):
        fs = sorted(factorize(n))
        if product(fs) != n:          # not squarefree
            continue
        if sum(n // p for p in fs) + 1 == n:
            by_sieve.append(n)
    assert by_recursion == by_sieve == [2, 6, 42, 1806, 47058]
    assert all(n % 2 == 0 for n in by_sieve)


def test_admissible_bound_is_stable_under_pruning():
    for X, expected in ((13, 202), (23, 554)):
        pool = Pool(200_000, X)
        pruned, witness, _ = search(pool, best=3000, prune=True)
        unpruned, _, _ = search(pool, best=3000, prune=False)
        assert pruned == unpruned == expected, X
        # third opinion: exact Fractions, no scaling, no bitmasks
        assert leaf_bound_exact(witness[0], X, 200_000) == expected, X
