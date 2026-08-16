"""Skeptic re-implementation for attempt 013 (the coverage Diophantine lemma).

Independence protocol (013 lead 1b): everything above the DIFF HARNESS marker
was written from the prose and formulas of
attempts/013-coverage-diophantine-lemma.md ALONE, before reading
explore/coverage_diophantine.py. Only the diff harness at the bottom imports
the committed module, and it was added after this file's selftest passed.

The lemma under review:

    For every real t in (0, 1) there exist integers a >= 1, b >= 1 with
        a + b < t * a * (b + 1) < a + b + 1.

Constructive witness claimed by 013 (q = b + 1):

    q = floor(1/t) + 1,   a = floor((q-1)/(q*t - 1)) + 1,   b = q - 1.

Usage:
    python cdsk_review.py selftest    # all independent checks, deterministic
    python cdsk_review.py diff        # witness diff vs the committed checker
    python cdsk_review.py gamma 144   # witness for an obtuse angle (degrees)

Stdlib only. Deterministic: RNG seeded.
"""

import random
import sys
from fractions import Fraction
from math import gcd


# ---------------------------------------------------------------------------
# The witness and the two arithmetic layers (from the record's formulas)
# ---------------------------------------------------------------------------

def witness(t):
    """Constructive witness (a, b) of 013 for rational t in (0, 1)."""
    t = Fraction(t)
    if not (0 < t < 1):
        raise ValueError("t must be in (0, 1)")
    p, r = t.numerator, t.denominator
    q = r // p + 1                    # floor(1/t) + 1
    # (q-1)/(q*t - 1) = (q-1)*r / (q*p - r); q*p - r > 0 since q > 1/t
    a = (q - 1) * r // (q * p - r) + 1
    b = q - 1
    return a, b


def ok_int(t, a, b):
    """Integer layer: r(a+b) < p*a*(b+1) < r(a+b+1) for t = p/r, pure ints."""
    t = Fraction(t)
    p, r = t.numerator, t.denominator
    mid = p * a * (b + 1)
    return r * (a + b) < mid < r * (a + b + 1)


def interval(a, q):
    """I(a, q) = (1/q + (q-1)/(a q), 1/q + 1/a) as exact Fractions."""
    lo = Fraction(1, q) + Fraction(q - 1, a * q)
    up = Fraction(1, q) + Fraction(1, a)
    return lo, up


def ok_frac(t, a, b):
    """Fraction layer: t strictly inside I(a, b+1)."""
    t = Fraction(t)
    lo, up = interval(a, b + 1)
    return lo < t < up


def brute_min(t, cap):
    """Minimal a+b witness by exhaustive scan of the integer layer.

    Shares no formulas with witness(): scans s = a+b ascending, b ascending.
    Returns (a, b) or None if nothing with a+b <= cap works.
    """
    t = Fraction(t)
    for s in range(2, cap + 1):
        for b in range(1, s):
            if ok_int(t, s - b, b):
                return s - b, b
    return None


# ---------------------------------------------------------------------------
# Independent checks
# ---------------------------------------------------------------------------

def check_form_equivalence(n=20000, seed=13):
    """Integer chain <=> t in I(a, q), on random (a, b, p, r)."""
    rng = random.Random(seed)
    bad = 0
    for _ in range(n):
        a = rng.randint(1, 500)
        b = rng.randint(1, 500)
        r = rng.randint(2, 10 ** 6)
        p = rng.randint(1, r - 1)
        t = Fraction(p, r)
        if ok_int(t, a, b) != ok_frac(t, a, b):
            bad += 1
    return n, bad


def check_overlap_criterion(qmax=120, extra=60):
    """Overlap of I(a+1, q) with I(a, q) iff a > q - 1; touching at a = q-1."""
    checked = bad = 0
    for q in range(2, qmax + 1):
        for a in range(1, 2 * q + extra):
            lo_a, _ = interval(a, q)
            _, up_a1 = interval(a + 1, q)
            checked += 1
            if (up_a1 > lo_a) != (a > q - 1):
                bad += 1
            if a == q - 1 and up_a1 != lo_a:
                bad += 1                       # the claimed exact touching
    return checked, bad


def check_union_identity(qmax=60):
    """Step 3: union over a >= q of I(a, q) is exactly (1/q, 2/q).

    Independent argument re-check: for t in (1/q, 2/q) take minimal a0 >= q
    with lower(a0) < t; then t < upper(a0). Tested on a dense exact grid,
    including every interval endpoint interior to the union.
    """
    checked = bad = 0
    for q in range(2, qmax + 1):
        pts = set()
        # dense rational grid strictly inside (1/q, 2/q)
        for k in range(1, 400):
            pts.add(Fraction(1, q) + Fraction(k, 400 * q))
        # every lower endpoint of I(a, q), a in [q, q+80] -- these are points
        # of the union that sit ON an endpoint and must be caught by a LATER
        # interval (013 lead 1c)
        for a in range(q, q + 81):
            lo, _ = interval(a, q)
            pts.add(lo)
        lo_bound, up_bound = Fraction(1, q), Fraction(2, q)
        for t in pts:
            if not (lo_bound < t < up_bound):
                continue
            checked += 1
            # find minimal a0 >= q with lower(a0) < t (must exist; walk up)
            a0 = q
            while interval(a0, q)[0] >= t:
                a0 += 1
                if a0 > 10 ** 7:
                    bad += 1
                    break
            else:
                lo, up = interval(a0, q)
                if not (lo < t < up):
                    bad += 1
    return checked, bad


def farey_sweep(rmax=300, brute_rmax=60):
    """Every reduced p/r, r <= rmax: witness passes BOTH layers, a >= q.

    For r <= brute_rmax additionally run the independent brute-force scan:
    it must find a witness with a+b <= the constructive a+b.
    """
    checked = bad = brute_checked = brute_bad = 0
    worst_ratio = Fraction(0)
    for r in range(2, rmax + 1):
        for p in range(1, r):
            if gcd(p, r) != 1:
                continue
            t = Fraction(p, r)
            a, b = witness(t)
            checked += 1
            if not (ok_int(t, a, b) and ok_frac(t, a, b) and a >= b + 1):
                bad += 1
            if r <= brute_rmax:
                brute_checked += 1
                mb = brute_min(t, a + b)
                if mb is None:
                    brute_bad += 1
                else:
                    ratio = Fraction(a + b, mb[0] + mb[1])
                    worst_ratio = max(worst_ratio, ratio)
    return checked, bad, brute_checked, brute_bad, worst_ratio


def endpoint_adversaries(qmax=40, arange=40):
    """Interval endpoints: member (a, q-1) must FAIL by strictness, yet the
    constructive witness for that t must succeed (with a different member)."""
    checked = bad = 0
    for q in range(2, qmax + 1):
        for a in range(q, q + arange):
            lo, up = interval(a, q)
            for t in (lo, up):
                if not (0 < t < 1):
                    continue
                checked += 1
                if ok_int(t, a, q - 1):
                    bad += 1                   # strict inequality must fail
                    continue
                wa, wb = witness(t)
                if (wa, wb) == (a, q - 1) or not ok_int(t, wa, wb):
                    bad += 1
        for t in (Fraction(1, q), Fraction(2, q)):
            if 0 < t < 1:
                checked += 1
                wa, wb = witness(t)
                if not (ok_int(t, wa, wb) and ok_frac(t, wa, wb)):
                    bad += 1
    return checked, bad


def edge_stress(seed=1301, n=2000):
    """Huge denominators near both ends of (0,1), plus random big rationals."""
    checked = bad = 0
    big = 2 * 10 ** 50
    ends = [Fraction(1, big), Fraction(1, big - 1), Fraction(2, big + 1),
            Fraction(big - 1, big), Fraction(big - 2, big - 1),
            Fraction(1, 2) + Fraction(1, big), Fraction(1, 2) - Fraction(1, big)]
    for t in ends:
        checked += 1
        a, b = witness(t)
        if not (ok_int(t, a, b) and ok_frac(t, a, b)):
            bad += 1
    rng = random.Random(seed)
    for _ in range(n):
        r = rng.randint(2, 10 ** 40)
        p = rng.randint(1, r - 1)
        t = Fraction(p, r)
        checked += 1
        a, b = witness(t)
        if not (ok_int(t, a, b) and ok_frac(t, a, b)):
            bad += 1
    return checked, bad


def check_a1_impossible(seed=7, n=5000):
    """013 lead 4: no t in (0,1) has an a = 1 witness.

    Algebra: a = 1 needs 1 + b < t(b+1), i.e. (b+1) < t(b+1), i.e. t > 1.
    Machine side: integer layer with a = 1 fails on random (t, b).
    """
    rng = random.Random(seed)
    bad = 0
    for _ in range(n):
        r = rng.randint(2, 10 ** 6)
        p = rng.randint(1, r - 1)
        b = rng.randint(1, 1000)
        if ok_int(Fraction(p, r), 1, b):
            bad += 1
    return n, bad


def halfdeg_arc_ratios():
    """013 check 6: constructive vs brute-minimal a+b on the 149 half-degree
    arcs gamma = 90.5, 91.0, ..., 164.5 (t = (180 - gamma)/90)."""
    ratios = []
    g = Fraction(181, 2)
    while g <= Fraction(329, 2):
        t = (180 - g) / 90
        a, b = witness(t)
        mb = brute_min(t, a + b)
        ratios.append((g, Fraction(a + b, mb[0] + mb[1]), (a, b), mb))
        g += Fraction(1, 2)
    worst = max(ratios, key=lambda x: x[1])
    mean = sum(x[1] for x in ratios) / len(ratios)
    return len(ratios), worst, mean


def spot_checks():
    """The record's consistency spot-checks, re-done."""
    out = []
    # gamma = 135  =>  t = 1/2; record says witness (5, 2)
    out.append(("gamma=135", witness(Fraction(1, 2))))
    # gamma = 144  =>  t = 2/5; record says "W(4, 3)" -- suspect
    t = Fraction(2, 5)
    out.append(("gamma=144 constructive", witness(t)))
    out.append(("gamma=144 (4,3) valid?", ok_int(t, 4, 3)))
    out.append(("gamma=144 brute minimal", brute_min(t, 50)))
    # the garbled Step-2 parenthetical: I(q-1, q) vs (1/q, 2/q)
    q = 5
    lo, up = interval(q - 1, q)
    out.append(("I(q-1,q) at q=5", (lo, up), "lower == 2/q?", lo == Fraction(2, q)))
    return out


def selftest():
    fails = 0

    n, bad = check_form_equivalence()
    print(f"form equivalence: {n} triples, {bad} mismatches")
    fails += bad

    checked, bad = check_overlap_criterion()
    print(f"overlap criterion (incl. a=q-1 touching): {checked} pairs, {bad} bad")
    fails += bad

    checked, bad = check_union_identity()
    print(f"union identity incl. interior endpoints: {checked} points, {bad} bad")
    fails += bad

    checked, bad, bchecked, bbad, worst = farey_sweep()
    print(f"farey sweep r<=300: {checked} witnesses, {bad} bad; "
          f"brute r<=60: {bchecked} scans, {bbad} bad, worst ratio {worst} "
          f"(~{float(worst):.2f})")
    fails += bad + bbad

    checked, bad = endpoint_adversaries()
    print(f"endpoint adversaries: {checked} boundary rationals, {bad} bad")
    fails += bad

    checked, bad = edge_stress()
    print(f"edge stress (denoms to 2e50): {checked} values, {bad} bad")
    fails += bad

    n, bad = check_a1_impossible()
    print(f"a=1 impossibility: {n} trials, {bad} spurious witnesses")
    fails += bad

    n, worst, mean = halfdeg_arc_ratios()
    print(f"half-degree arcs: {n} arcs, worst overshoot x{float(worst[1]):.2f} "
          f"at gamma={worst[0]} (constructive {worst[2]}, minimal {worst[3]}), "
          f"mean x{float(mean):.2f}")

    for row in spot_checks():
        print("spot:", row)

    print("SELFTEST", "FAIL" if fails else "PASS", f"({fails} failures)")
    return 1 if fails else 0


# ---------------------------------------------------------------------------
# DIFF HARNESS -- the only part allowed to touch the committed module.
# ---------------------------------------------------------------------------

def diff_vs_committed(rmax=300):
    """Compare this file's witness against the committed checker's on the
    full Farey sweep r <= rmax. Added after selftest passed."""
    import importlib.util
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    spec = importlib.util.spec_from_file_location(
        "committed", os.path.join(here, "coverage_diophantine.py"))
    committed = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(committed)

    checked = mismatches = 0
    for r in range(2, rmax + 1):
        for p in range(1, r):
            if gcd(p, r) != 1:
                continue
            t = Fraction(p, r)
            mine = witness(t)
            theirs = committed.witness(t)
            checked += 1
            if mine != tuple(theirs):
                mismatches += 1
                if mismatches <= 10:
                    print(f"MISMATCH t={t}: mine {mine}, committed {theirs}")
    print(f"diff vs committed: {checked} rationals, {mismatches} mismatches")
    return 1 if mismatches else 0


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "selftest"
    if cmd == "selftest":
        sys.exit(selftest())
    elif cmd == "diff":
        sys.exit(diff_vs_committed())
    elif cmd == "gamma":
        g = Fraction(sys.argv[2])
        t = (180 - g) / 90
        a, b = witness(t)
        print(f"gamma={g} t={t} witness W({a},{b}) "
              f"ok_int={ok_int(t, a, b)} ok_frac={ok_frac(t, a, b)}")
    else:
        raise SystemExit(f"unknown command {cmd}")
