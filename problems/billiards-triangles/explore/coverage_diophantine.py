#!/usr/bin/env python3
"""Exact checker for the coverage Diophantine lemma (006 lead 2 / queue 12).

LEMMA.  For every real t in (0, 1) there exist integers a >= 1, b >= 1 with

        a + b  <  t * a * (b+1)  <  a + b + 1.

Substituting q = b+1 >= 2 and dividing by a*q, the two inequalities say
exactly that t lies in the open interval

        I(a, q) = ( 1/q + (q-1)/(a*q) ,  1/q + 1/a ).

For fixed q the intervals I(a, q), a >= q, form an overlapping chain
(consecutive overlap iff a > q-1) whose union is (1/q, 2/q); the chain
TOUCHES at a = q-1 (upper endpoint of I(q, q-1)... see selftest: the
endpoints coincide exactly, so a >= q is sharp).  Since the (1/q, 2/q)
cover (0, 1) as q runs over the integers >= 2, the lemma follows, with the
constructive witness

        q = floor(1/t) + 1,          a = floor((q-1)/(q*t - 1)) + 1,
        b = q - 1.

This module verifies all of that in exact arithmetic, twice over:

  * the constructive witness path works in Fraction arithmetic on the
    interval form of the inequality;
  * an independent brute-force path scans a+b = 2, 3, ... and checks the
    ORIGINAL inequality in pure integer arithmetic
    ( r*(a+b) < p*a*(b+1) < r*(a+b+1)  for t = p/r ),
    sharing no formulas with the witness construction.

Commands (deterministic, stdlib only):

    python coverage_diophantine.py selftest            # full stress suite
    python coverage_diophantine.py witness P/R         # witness for t = P/R
    python coverage_diophantine.py gamma G_NUM/G_DEN   # witness for an
                                       obtuse angle gamma (deg); uses
                                       t = (180 - gamma)/90

For an obtuse angle gamma the witness (a, b) names the family member
W(a, b) whose SPECULATION window (006) strictly contains gamma; minimal
word length for the member is 4*(a+b) + 2.  The lemma itself is
unconditional; its billiards meaning is conditional on 006's unproven
birth law and exact-window claim.
"""

from __future__ import annotations

import random
import sys
from fractions import Fraction

# ---------------------------------------------------------------------------
# the two exact forms of the condition


def holds_fraction(t: Fraction, a: int, b: int) -> bool:
    """Interval form, Fraction arithmetic: t interior to I(a, b+1)."""
    q = b + 1
    lo = Fraction(1, q) + Fraction(q - 1, a * q)
    hi = Fraction(1, q) + Fraction(1, a)
    return lo < t < hi


def holds_integer(p: int, r: int, a: int, b: int) -> bool:
    """Original form for t = p/r, pure integer arithmetic (no Fraction)."""
    mid = p * a * (b + 1)
    return r * (a + b) < mid < r * (a + b + 1)


# ---------------------------------------------------------------------------
# constructive witness


def witness(t: Fraction) -> tuple[int, int]:
    """The (a, b) from the proof: q = floor(1/t)+1, a = floor((q-1)/(qt-1))+1.

    Raises AssertionError (never expected) if any derived bound fails; the
    selftest leans on these assertions.
    """
    if not (0 < t < 1):
        raise ValueError("t must be in (0, 1)")
    q = (1 / t).__floor__() + 1
    # q is the integer promised by |(1/t, 2/t)| = 1/t > 1:
    assert q >= 2 and Fraction(1, q) < t < Fraction(2, q)
    d = q * t - 1                       # > 0 since t > 1/q
    a = ((q - 1) / d).__floor__() + 1
    # a lies in the open integer window ((q-1)/(qt-1), q/(qt-1)), and the
    # lower bound already forces a >= q because qt - 1 < 1:
    assert Fraction(q - 1, 1) / d < a < Fraction(q, 1) / d
    assert a >= q
    b = q - 1
    assert holds_fraction(t, a, b)
    return a, b


def witness_bruteforce(p: int, r: int, cap: int = 100_000) -> tuple[int, int]:
    """Minimal-(a+b) witness by exhaustive integer scan; independent of the
    constructive formulas.  Guaranteed to terminate by the lemma; `cap` on
    a+b only guards against a checker bug."""
    if not (0 < p < r):
        raise ValueError("need 0 < p < r")
    for s in range(2, cap + 1):
        for a in range(1, s):
            if holds_integer(p, r, a, s - a):
                return a, s - a
    raise AssertionError(f"no witness with a+b <= {cap} for t = {p}/{r}")


# ---------------------------------------------------------------------------
# stress suite


def _check_both_paths(p: int, r: int, brute: bool) -> tuple[int, int]:
    """Constructive witness + exact check in BOTH arithmetic layers; optional
    brute-force cross-check.  Returns the constructive (a, b)."""
    t = Fraction(p, r)
    a, b = witness(t)
    assert holds_fraction(t, a, b), (p, r, a, b)
    assert holds_integer(p, r, a, b), (p, r, a, b)
    if brute:
        a2, b2 = witness_bruteforce(p, r)
        assert holds_fraction(t, a2, b2), (p, r, a2, b2)
        assert a2 + b2 <= a + b, (p, r, a, b, a2, b2)
    return a, b


def selftest() -> int:
    checks = 0

    # 1. Equivalence of the two forms (random exact spot check).
    rng = random.Random(20260815)
    for _ in range(20_000):
        a = rng.randint(1, 60)
        b = rng.randint(1, 60)
        p = rng.randint(1, 199)
        assert holds_fraction(Fraction(p, 200), a, b) == holds_integer(
            p, 200, a, b
        )
        checks += 1
    print(f"[1] interval form == original form at {checks} random points")

    # 2. Chain overlap is exactly a > q-1; the chain touches at a = q-1.
    n2 = 0
    for q in range(2, 121):
        for a in range(1, 2 * q + 60):
            overlap = Fraction(1, a + 1) > Fraction(q - 1, a * q)
            assert overlap == (a > q - 1), (a, q)
            n2 += 1
        # sharpness: at a = q-1 the would-be overlap is exact equality
        if q >= 3:
            assert Fraction(1, q) == Fraction(q - 1, (q - 1) * q)
    checks += n2
    print(f"[2] consecutive-overlap criterion a > q-1 exact at {n2} pairs; "
          "touching (equality) confirmed at a = q-1")

    # 3. Farey sweep: every rational t = p/r, r <= 300, both layers.
    n3 = 0
    worst = (0, (0, 0, 0))  # (a+b, (p, r, ...))
    from math import gcd
    for r in range(2, 301):
        for p in range(1, r):
            if gcd(p, r) != 1:
                continue
            a, b = _check_both_paths(p, r, brute=(r <= 60))
            if a + b > worst[0]:
                worst = (a + b, (p, r, a))
            n3 += 1
    checks += n3
    print(f"[3] Farey sweep all p/r, r <= 300: {n3} values, constructive "
          f"witness valid in both layers (brute-force cross-check for "
          f"r <= 60); largest witness a+b = {worst[0]} at t = "
          f"{worst[1][0]}/{worst[1][1]}")

    # 4. Boundary adversaries: window endpoints and 1/q, 2/q exactly.
    #    At its own endpoints I(a, q) fails by strictness; the lemma must
    #    recover via a DIFFERENT (a, q).
    n4 = 0
    for q in range(2, 41):
        for a in range(q, q + 31):
            for t in (
                Fraction(1, q) + Fraction(q - 1, a * q),   # lower endpoint
                Fraction(1, q) + Fraction(1, a),           # upper endpoint
            ):
                if not (0 < t < 1):
                    continue
                assert not holds_fraction(t, a, q - 1)     # strictness bites
                _check_both_paths(t.numerator, t.denominator, brute=False)
                n4 += 1
        for t in (Fraction(1, q), Fraction(2, q)):
            if 0 < t < 1:
                _check_both_paths(t.numerator, t.denominator, brute=False)
                n4 += 1
    checks += n4
    print(f"[4] endpoint adversaries: {n4} exact boundary rationals all "
          "recovered by a different member")

    # 5. Near the edges of (0, 1): tiny and huge denominators.
    n5 = 0
    edge_cases = []
    for k in range(1, 19):
        big = 10 ** k
        edge_cases += [(1, big), (1, big + 1), (1, big + 3), (2, big + 1),
                       (big - 1, big), (big - 3, big + 1) if big > 3 else (1, 2)]
    edge_cases += [(1, 10 ** 50), (10 ** 50 - 1, 10 ** 50),
                   (3, 10 ** 50 + 7), (10 ** 50 + 1, 2 * 10 ** 50 + 1)]
    for p, r in edge_cases:
        from math import gcd as _g
        g = _g(p, r)
        _check_both_paths(p // g, r // g, brute=False)
        n5 += 1
    checks += n5
    print(f"[5] edge stress near 0 and 1 (denominators to 2e50): {n5} values")

    # 6. Random large rationals.
    n6 = 0
    for _ in range(2_000):
        r = rng.randint(2, 10 ** 40)
        p = rng.randint(1, r - 1)
        from math import gcd as _g
        g = _g(p, r)
        _check_both_paths(p // g, r // g, brute=False)
        n6 += 1
    checks += n6
    print(f"[6] random rationals, denominators to 1e40: {n6} values")

    # 7. Informational: constructive vs minimal witness on the coverage-grid
    #    angles of 006 (gamma = 90.5 .. 165 step 0.5 -> t = (180-gamma)/90).
    ratios = []
    for g2 in range(181, 330):  # gamma*2
        t = Fraction(360 - g2, 180)
        a, b = witness(t)
        a2, b2 = witness_bruteforce(t.numerator, t.denominator)
        ratios.append((a + b) / (a2 + b2))
    print(f"[7] info: on 006's 149 half-degree arcs the constructive "
          f"witness's a+b is within x{max(ratios):.2f} of minimal "
          f"(mean x{sum(ratios)/len(ratios):.2f}); minimal word length "
          "4(a+b)+2 is governed by the brute-force column")

    print(f"\nselftest PASS: {checks} exact checks, all strict inequalities "
          "verified in two independent arithmetic layers")
    return 0


# ---------------------------------------------------------------------------


def _parse_fraction(s: str) -> Fraction:
    if "/" in s:
        num, den = s.split("/", 1)
        return Fraction(int(num), int(den))
    return Fraction(s)


def main(argv: list[str]) -> int:
    if len(argv) >= 1 and argv[0] == "selftest":
        return selftest()
    if len(argv) == 2 and argv[0] == "witness":
        t = _parse_fraction(argv[1])
        a, b = witness(t)
        ok = holds_integer(t.numerator, t.denominator, a, b)
        print(f"t = {t}: a = {a}, b = {b}  "
              f"(q = {b+1}, a+b = {a+b}, exact re-check: {ok})")
        return 0
    if len(argv) == 2 and argv[0] == "gamma":
        gamma = _parse_fraction(argv[1])
        if not (90 < gamma < 180):
            print("gamma must be strictly between 90 and 180", file=sys.stderr)
            return 2
        t = (180 - gamma) / 90
        a, b = witness(t)
        ok = holds_integer(t.numerator, t.denominator, a, b)
        print(f"gamma = {gamma} deg -> t = {t}: member W({a},{b}), "
              f"word length {4*(a+b)+2}, exact re-check: {ok}")
        print("(strict window containment assumes 006's SPECULATION "
              "birth law; the arithmetic itself is exact)")
        return 0
    print(__doc__)
    return 2


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
