"""Cross-checks for the harness: two field implementations, three colouring methods.

Run:  python harness/hadwiger-nelson/selftest.py
Deterministic apart from a fixed-seed random sample; prints PASS/FAIL per block.
"""

from __future__ import annotations

import os
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from itertools import product

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from exact_field import MQField, field_for, squarefree_part  # noqa: E402
from exact_field_poly import PolyQuad  # noqa: E402
import colouring as C  # noqa: E402


def _to_poly(el, mods):
    """Translate an exact_field element into the reference representation."""
    k = len(mods)
    terms = {}
    for mask, c in enumerate(el.coords()):
        e = tuple((mask >> i) & 1 for i in range(k))
        terms[e] = c
    return PolyQuad(mods, terms)


def check_fields(trials=400, seed=20260727):
    rng = random.Random(seed)
    failures = []
    for gens in [(3,), (3, 11), (3, 5), (3, 7), (2, 3, 5)]:
        F = MQField(gens)
        mods = gens
        for _ in range(trials):
            def rnd():
                num = [rng.randint(-9, 9) for _ in range(F.dim)]
                den = rng.randint(1, 7)
                return type(F.one())(F, tuple(num), den)

            a, b = rnd(), rnd()
            for name, fast, ref in [
                ("add", a + b, _to_poly(a, mods) + _to_poly(b, mods)),
                ("sub", a - b, _to_poly(a, mods) - _to_poly(b, mods)),
                ("mul", a * b, _to_poly(a, mods) * _to_poly(b, mods)),
                ("sq", a * a, _to_poly(a, mods) * _to_poly(a, mods)),
            ]:
                if _to_poly(fast, mods) != ref:
                    failures.append((gens, name, repr(a), repr(b)))
    return failures


def check_numeric(trials=120, seed=7, prec=140):
    """Third net: the exact value and a 140-digit numeric evaluation must agree.

    Catches a wrong multiplication table that both symbolic implementations
    could in principle share, since the numeric path knows nothing about the
    subset basis.
    """
    getcontext().prec = prec
    rng = random.Random(seed)
    worst = Decimal(0)
    for gens in [(3,), (3, 11), (3, 5), (2, 3, 5)]:
        F = MQField(gens)
        E = type(F.one())
        for _ in range(trials):
            a = E(F, tuple(rng.randint(-9, 9) for _ in range(F.dim)), rng.randint(1, 7))
            b = E(F, tuple(rng.randint(-9, 9) for _ in range(F.dim)), rng.randint(1, 7))
            lhs = (a * b).to_decimal(prec)
            rhs = a.to_decimal(prec) * b.to_decimal(prec)
            worst = max(worst, abs(lhs - rhs))
    return worst


def check_inverse():
    bad = []
    for gens in [(3,), (3, 11), (2, 3, 5)]:
        F = MQField(gens)
        E = type(F.one())
        for coords in product(range(-2, 3), repeat=F.dim):
            if not any(coords):
                continue
            a = E(F, coords, 1)
            if a * a.inverse() != F.one():
                bad.append((gens, coords))
        if len(gens) == 3:
            break
    return bad


def check_sqrt_int():
    F = field_for(3, 11)
    assert repr(F) == "Q(sqrt3, sqrt11)", repr(F)
    assert field_for(3, 27).gens == (3,), field_for(3, 27).gens
    assert field_for(3, 33).gens == (3, 11), field_for(3, 33).gens
    assert field_for(3, 75).gens == (3,)
    assert F.sqrt_int(99) * F.sqrt_int(99) == F.rational(99)
    assert F.sqrt_int(27) * F.sqrt_int(27) == F.rational(27)
    try:
        F.sqrt_int(5)
    except ValueError:
        return True
    return False


# --- colouring cross-checks --------------------------------------------------

def _random_graph(n, p, rng):
    return [(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < p]


def check_colouring_methods(trials=60, seed=99):
    rng = random.Random(seed)
    bad = []
    for _ in range(trials):
        n = rng.randint(3, 9)
        edges = _random_graph(n, rng.choice([0.3, 0.5, 0.7]), rng)
        chi_bt = C.chromatic_number(n, edges)
        chi_poly = C.chromatic_number_via_polynomial(n, edges)
        chi_ie = C.chromatic_number_ie(n, edges)
        chi_bf = None
        for k in range(0, n + 1):
            if C.brute_force_k(n, edges, k) is not None:
                chi_bf = k
                break
        if not (chi_bt == chi_poly == chi_bf == chi_ie):
            bad.append((n, edges, chi_bt, chi_poly, chi_bf, chi_ie))
            continue
        # the polynomial and the inclusion-exclusion oracle must agree on the
        # *verdict* for every k (they count different things, so not on values)
        poly = C.chromatic_polynomial(n, edges)
        for k in range(0, n + 1):
            if (C.poly_eval(poly, k) > 0) != C.is_k_colourable_ie(n, edges, k):
                bad.append(("verdict-mismatch", n, edges, k))
        col = C.k_colouring(n, edges, chi_bt)
        if not C.verify_colouring(n, edges, col):
            bad.append(("improper", n, edges, col))
    return bad


def check_fast_solver(trials=400, seed=4242):
    """The scaling solver must agree with the exhaustive ones, verdict and colouring.

    Random graphs are deliberately small so `colouring.py` can decide them
    outright: the fast solver's reductions (k-core, components, propagation)
    are exactly where a scaling solver goes wrong, and they fire on small
    sparse graphs too.
    """
    import colouring_fast as F

    rng = random.Random(seed)
    bad = []
    for _ in range(trials):
        n = rng.randint(1, 11)
        p = rng.choice([0.15, 0.3, 0.5, 0.75])
        edges = _random_graph(n, p, rng)
        for k in range(1, 6):
            slow = C.k_colouring(n, edges, k) is not None
            fast_col = F.k_colourable(n, edges, k)
            fast = fast_col is not None
            ie = C.is_k_colourable_ie(n, edges, k)
            if not (slow == fast == ie):
                bad.append(("verdict", n, edges, k, slow, fast, ie))
                continue
            if fast and not F.verify(n, edges, fast_col, k):
                bad.append(("improper colouring", n, edges, k, fast_col))
        # and the same with the reduction switched off, to isolate it
        for k in (3, 4):
            a = F.k_colourable(n, edges, k, reduce=True) is not None
            b = F.k_colourable(n, edges, k, reduce=False) is not None
            if a != b:
                bad.append(("k-core reduction changed the verdict", n, edges, k))
    # heuristic layer: a found colouring must be proper, and must never
    # contradict the complete search
    for _ in range(120):
        n = rng.randint(4, 14)
        edges = _random_graph(n, rng.choice([0.2, 0.4, 0.6]), rng)
        for k in (3, 4):
            verdict, col, how = F.decide_k_colourable(
                n, edges, k, tries=3, steps=4000, seed=rng.randrange(10 ** 6))
            truth = C.is_k_colourable_ie(n, edges, k)
            if verdict is not True and verdict is not False:
                continue
            if verdict != truth:
                bad.append(("decide verdict", n, edges, k, verdict, truth, how))
            if verdict and not F.verify(n, edges, col, k):
                bad.append(("decide colouring improper", n, edges, k, col))

    # contraction helper
    for _ in range(60):
        n = rng.randint(3, 9)
        edges = _random_graph(n, 0.4, rng)
        u, v = rng.sample(range(n), 2)
        if (min(u, v), max(u, v)) in set(edges):
            continue
        cn, ce = F.contract(n, edges, u, v)
        # contracting u,v is k-colourable iff some k-colouring has c(u) == c(v)
        for k in (3, 4):
            merged = F.k_colourable(cn, ce, k) is not None
            found = False
            for col in _all_colourings(n, edges, k):
                if col[u] == col[v]:
                    found = True
                    break
            if merged != found:
                bad.append(("contraction", n, edges, u, v, k, merged, found))
    return bad


def _all_colourings(n, edges, k):
    for tail in product(range(k), repeat=n - 1):
        col = (0,) + tail
        if all(col[a] != col[b] for a, b in edges):
            yield col


def main():
    ok = True

    f = check_fields()
    print(f"[field] fast vs reference implementation: {'PASS' if not f else 'FAIL'}"
          f"  ({len(f)} mismatches)")
    ok &= not f

    w = check_numeric()
    print(f"[field] exact vs 140-digit numeric: worst |lhs-rhs| = {w:.3e}"
          f"  {'PASS' if w < Decimal('1e-100') else 'FAIL'}")
    ok &= w < Decimal("1e-100")

    b = check_inverse()
    print(f"[field] inverse: {'PASS' if not b else 'FAIL'}  ({len(b)} failures)")
    ok &= not b

    s = check_sqrt_int()
    print(f"[field] sqrt_int / field_for: {'PASS' if s else 'FAIL'}")
    ok &= bool(s)

    c = check_colouring_methods()
    print(f"[colour] DSATUR vs brute force vs chromatic polynomial: "
          f"{'PASS' if not c else 'FAIL'}  ({len(c)} mismatches)")
    ok &= not c

    s = check_fast_solver()
    print(f"[colour] scaling solver vs exhaustive methods: "
          f"{'PASS' if not s else 'FAIL'}  ({len(s)} mismatches)")
    for row in s[:5]:
        print("         ", row)
    ok &= not s

    print("ALL PASS" if ok else "FAILURES PRESENT")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
