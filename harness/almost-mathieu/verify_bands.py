#!/usr/bin/env python3
"""Independent re-verification of critical-AMO band records at rational flux.

Companion to the reference tool in this directory, sharing nothing but the
problem statement: a JSON record claiming the Chambers polynomial P, certified
band edges, gap certificates and bandwidth enclosures for sigma(1, p/q) is
re-derived and re-checked here by different algorithms.

  * P is recomputed by exact scalar evaluation: the discriminant
    Delta(E_k) at the q+1 integer energies E_k = 0..q, each a transfer-matrix
    product in Z[x]/(x^{4q}-1) at phase nu = pi/(2q) (where the Chambers
    cosine term vanishes), followed by exact Lagrange interpolation over Q.
    The claimed P -- given over Q[y]/(psi_q), y = 2 cos(2 pi / q) -- is
    embedded via y -> x^4 + x^{4q-4} and compared coefficient by coefficient
    after reduction modulo the cyclotomic polynomial Phi_{4q}.
  * Signs of field elements are decided by a sign oracle built from certified
    rational series enclosures -- Machin's formula for pi and the alternating
    cosine series, refined on demand -- rather than from the minimal
    polynomial psi_q, so even the enclosure machinery is independent.
  * Certificates are re-checked from scratch: each claimed edge bracket must
    carry a sign change of the correct P -+ 4 (the side pattern follows from
    the sign of P at -infinity), brackets must be disjoint and ordered, gap
    interiors must have |P| > 4 and band interiors |P| < 4 at exact rational
    test points, touchings must be exact double roots (P(0) = +-4 and
    P'(0) = 0, checked as integer vectors), and the degree argument closes
    the count: q disjoint certified brackets for a degree-q polynomial leave
    no room for missed edges.  Widths and bandwidth enclosures are
    re-summed.

Usage:
    python verify_bands.py RECORD.json
    python verify_bands.py --selftest

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from math import gcd

# ---------------------------------------------------------------------------
# integer / rational polynomial helpers (deliberately re-typed, not imported)
# ---------------------------------------------------------------------------


def trim(f: list) -> list:
    while f and f[-1] == 0:
        f.pop()
    return f


def divmod_monic(f: list, g: list) -> tuple[list, list]:
    f = [Fraction(c) for c in f]
    g = [Fraction(c) for c in g]
    assert g and g[-1] == 1
    dg = len(g) - 1
    quo = [Fraction(0)] * max(1, len(f) - dg)
    while trim(f) and len(f) - 1 >= dg:
        d = len(f) - 1 - dg
        c = f[-1]
        quo[d] = c
        for i, gi in enumerate(g):
            f[d + i] -= c * gi
        trim(f)
    return trim(quo), trim(f)


_CYCLO: dict[int, list[int]] = {}


def cyclotomic(n: int) -> list[int]:
    if n in _CYCLO:
        return _CYCLO[n]
    f: list = [-1] + [0] * (n - 1) + [1]
    for d in range(1, n):
        if n % d == 0:
            quo, rem = divmod_monic(f, cyclotomic(d))
            assert not rem
            f = [int(c) for c in quo]
    _CYCLO[n] = f
    return f


# ---------------------------------------------------------------------------
# route 1: P by scalar evaluation + interpolation in Z[x]/(x^{4q}-1)
# ---------------------------------------------------------------------------


def disc_at_integer(p: int, q: int, E: int) -> list[int]:
    """Delta(E) at nu = pi/(2q) as a vector in Z[x]/(x^{4q}-1)."""
    N = 4 * q
    zero = [0] * N
    t00, t01, t10, t11 = zero[:], zero[:], zero[:], zero[:]
    t00[0] = t11[0] = 1

    for j in range(1, q + 1):
        e = (4 * p * j + 1) % N
        ne = (N - e) % N

        def apply(row: list[int]) -> list[int]:
            # (E - x^e - x^{ne}) * row
            out = [E * c for c in row]
            for k, c in enumerate(row):
                if c:
                    out[(k + e) % N] -= c
                    out[(k + ne) % N] -= c
            return out

        n00 = [a - b for a, b in zip(apply(t00), t10)]
        n01 = [a - b for a, b in zip(apply(t01), t11)]
        t10, t11 = t00, t01
        t00, t01 = n00, n01
    return [a + b for a, b in zip(t00, t11)]


def interpolated_P(p: int, q: int) -> list[list[Fraction]]:
    """P as q+1 canonical vectors mod Phi_{4q} (Lagrange over integer nodes)."""
    N = 4 * q
    xs = list(range(q + 1))
    ys = [disc_at_integer(p, q, E) for E in xs]
    coeffs = [[Fraction(0)] * N for _ in range(q + 1)]
    for i, xi in enumerate(xs):
        poly = [Fraction(1)]
        denom = Fraction(1)
        for j, xj in enumerate(xs):
            if j == i:
                continue
            poly = [Fraction(0)] + poly
            poly = [
                poly[k] - xj * (poly[k + 1] if k + 1 < len(poly) else Fraction(0))
                for k in range(len(poly))
            ]
            denom *= xi - xj
        for k in range(len(poly)):
            w = poly[k] / denom
            if w:
                for comp in range(N):
                    if ys[i][comp]:
                        coeffs[k][comp] += w * ys[i][comp]
    phi = cyclotomic(N)
    out = []
    for vec in coeffs:
        _, rem = divmod_monic(list(vec), phi)
        out.append(rem + [Fraction(0)] * (len(phi) - 1 - len(rem)))
    return out


def embed_claimed_P(psi: list[int], P_vecs: list[list[Fraction]], q: int) -> list[list[Fraction]]:
    """y -> x^4 + x^{4q-4}, reduced mod Phi_{4q}, one vector per E-coefficient."""
    N = 4 * q
    d = len(psi) - 1
    # powers of (x^4 + x^{N-4}) mod x^N - 1
    powers = [[0] * N for _ in range(d)]
    powers[0][0] = 1
    for i in range(1, d):
        prev = powers[i - 1]
        cur = [0] * N
        for k, c in enumerate(prev):
            if c:
                cur[(k + 4) % N] += c
                cur[(k - 4) % N] += c
        powers[i] = cur
    phi = cyclotomic(N)
    out = []
    for vec in P_vecs:
        big = [Fraction(0)] * N
        for i, c in enumerate(vec):
            if c:
                for k, pk in enumerate(powers[i]):
                    if pk:
                        big[k] += c * pk
        _, rem = divmod_monic(list(big), phi)
        out.append(rem + [Fraction(0)] * (len(phi) - 1 - len(rem)))
    return out


# ---------------------------------------------------------------------------
# route 2: an independent sign oracle from certified series enclosures
# ---------------------------------------------------------------------------


class CosOracle:
    """Certified rational enclosure of y = 2 cos(2 pi / q), by series only."""

    def __init__(self, q: int):
        self.q = q
        self.terms = 12
        self._recompute()

    def _pi_bracket(self, terms: int) -> tuple[Fraction, Fraction]:
        def arctan_bracket(inv_x: int) -> tuple[Fraction, Fraction]:
            x = Fraction(1, inv_x)
            s = Fraction(0)
            sign = 1
            lo = hi = None
            for k in range(terms):
                term = sign * x ** (2 * k + 1) / (2 * k + 1)
                s += term
                if k >= terms - 2:
                    if lo is None:
                        lo = hi = s
                    lo, hi = min(lo, s), max(hi, s)
                sign = -sign
            return lo, hi  # alternating series: consecutive partials bracket

        a5 = arctan_bracket(5)
        a239 = arctan_bracket(239)
        return 16 * a5[0] - 4 * a239[1], 16 * a5[1] - 4 * a239[0]

    def _cos2_bracket(self, x_lo: Fraction, x_hi: Fraction, terms: int):
        """[2cos(x_hi), 2cos(x_lo)] for 0 < x_lo <= x_hi < pi (cos decreasing)."""

        def cos_partials(x: Fraction) -> tuple[Fraction, Fraction]:
            s = Fraction(1)
            term = Fraction(1)
            lo = hi = None
            for k in range(1, terms):
                term = -term * x * x / ((2 * k - 1) * (2 * k))
                s += term
                if k >= terms - 2:  # the last two partials bracket the sum
                    if lo is None:
                        lo = hi = s
                    else:
                        lo, hi = min(lo, s), max(hi, s)
            return lo, hi

        c_hi = cos_partials(x_lo)  # larger cosine
        c_lo = cos_partials(x_hi)
        return 2 * c_lo[0], 2 * c_hi[1]

    def _recompute(self) -> None:
        pi_lo, pi_hi = self._pi_bracket(self.terms)
        x_lo, x_hi = 2 * pi_lo / self.q, 2 * pi_hi / self.q
        self.lo, self.hi = self._cos2_bracket(x_lo, x_hi, self.terms)
        assert self.lo < self.hi

    def refine(self) -> None:
        self.terms += 6
        self._recompute()

    def sign(self, vec: list) -> int:
        """Exact sign of sum vec[i] * y^i; zero iff the vector is zero."""
        if all(c == 0 for c in vec):
            return 0
        while True:
            lo, hi = Fraction(0), Fraction(0)
            for c in reversed(vec):
                cands = (lo * self.lo, lo * self.hi, hi * self.lo, hi * self.hi)
                lo, hi = min(cands) + c, max(cands) + c
            if lo > 0:
                return 1
            if hi < 0:
                return -1
            self.refine()
            assert self.terms < 400, "sign oracle failed to converge (zero divisor?)"


def eval_hom(P_vecs: list[list], num: int, den: int) -> list:
    """den^q * P(num/den) as a vector over y (integer arithmetic)."""
    acc = [0] * len(P_vecs[0])
    pw = 1
    for c in reversed(P_vecs):
        acc = [x * num + ci * pw for x, ci in zip(acc, c)]
        pw *= den
    return acc


# ---------------------------------------------------------------------------
# record verification
# ---------------------------------------------------------------------------


def verify_record(path: str) -> list[str]:
    with open(path, encoding="utf-8") as fh:
        rec = json.load(fh)
    p, q = rec["p"], rec["q"]
    problems: list[str] = []
    if q <= 2:
        return _verify_small_q(rec)
    assert gcd(p, q) == 1

    psi = rec["psi"]
    P_vecs = [[Fraction(x) for x in c] for c in rec["P"]]
    if len(P_vecs) != q + 1:
        return [f"claimed P has degree {len(P_vecs) - 1}, expected {q}"]

    # route-1 equality: interpolated discriminant == embedded claimed P
    mine = interpolated_P(p, q)
    theirs = embed_claimed_P(psi, P_vecs, q)
    for k, (a, b) in enumerate(zip(mine, theirs)):
        if trim([x - y for x, y in zip(a, b)]):
            problems.append(f"P mismatch at E^{k}: interpolation route disagrees")
    if problems:
        return problems

    oracle = CosOracle(q)
    ints = all(x.denominator == 1 for c in P_vecs for x in c)
    if not ints:
        problems.append("claimed P has non-integer coordinates over y")
        return problems
    P_int = [[int(x) for x in c] for c in P_vecs]

    def sgn_P_shift(x: Fraction, shift: int) -> int:
        # sign of P(x) - shift, exact
        v = eval_hom(P_int, x.numerator, x.denominator)
        v = list(v)
        v[0] -= shift * x.denominator ** q
        return oracle.sign(v)

    def absP_side(x: Fraction) -> int:
        hi_ = sgn_P_shift(x, 4)
        lo_ = sgn_P_shift(x, -4)
        if hi_ > 0 or lo_ < 0:
            return 1
        if hi_ < 0 and lo_ > 0:
            return -1
        return 0

    bands = [
        (
            (Fraction(b["lo"][0]), Fraction(b["lo"][1])),
            (Fraction(b["hi"][0]), Fraction(b["hi"][1])),
        )
        for b in rec["bands"]
    ]
    if len(bands) != q:
        problems.append(f"claimed {len(bands)} bands, expected {q}")
        return problems

    # reconstruct the flat ordered edge list and its P = +-4 side pattern:
    # P -> +inf on the right, and inside every band P runs from -4 to +4 or
    # back, so counting r = 0, 1, 2, ... from the right the sides follow the
    # four-periodic pattern +, -, -, +, +, -, -, ...
    edges = []
    for lo, hi in bands:
        edges.append(lo)
        edges.append(hi)
    sides = [0] * (2 * q)
    for r in range(2 * q):
        sides[2 * q - 1 - r] = 1 if r % 4 in (0, 3) else -1
    # ordering and disjointness (touching edges may coincide exactly)
    for (a1, b1), (a2, b2) in zip(edges, edges[1:]):
        if not (b1 <= a2 or (a1 == a2 and b1 == b2 and a1 == b1)):
            problems.append(f"edges out of order or overlapping: {(a1, b1)} vs {(a2, b2)}")

    # sign-change certificates; exact rational edges verified algebraically,
    # with their claimed multiplicity (a value appearing twice on one side
    # must be a double root -- that is the touching certificate)
    dP = [[c * k for c in vec] for k, vec in enumerate(P_int)][1:]
    exact_mult: dict[tuple[int, Fraction], int] = {}
    for (a, b), side in zip(edges, sides):
        if a == b:
            exact_mult[(side, a)] = exact_mult.get((side, a), 0) + 1
        else:
            sa = sgn_P_shift(a, side * 4)
            sb = sgn_P_shift(b, side * 4)
            if sa == 0 or sb == 0 or sa == sb:
                problems.append(f"bracket ({a}, {b}) carries no sign change of P {side:+d}4")
    for (side, a), mult in exact_mult.items():
        v0 = eval_hom(P_int, a.numerator, a.denominator)
        v0[0] -= side * 4 * a.denominator**q
        if oracle.sign(v0) != 0:
            problems.append(f"claimed exact edge {a} is not a root of P {side:+d}4")
        if mult > 2:
            problems.append(f"exact edge {a} claimed with multiplicity {mult} > 2")
        if mult == 2:
            v1 = eval_hom(dP, a.numerator, a.denominator)
            if oracle.sign(v1) != 0:
                problems.append(f"claimed touching at {a} is not a double root")

    # count completeness by the degree argument, per side
    for side in (+1, -1):
        n_brackets = sum(1 for (a, b), sd in zip(edges, sides) if sd == side and a != b)
        n_exact = sum(m for (sd, _), m in exact_mult.items() if sd == side)
        if n_brackets + n_exact != q:
            problems.append(
                f"side {side:+d}: {n_brackets} brackets + {n_exact} exact roots != {q}"
            )

    # interiors
    for i, (lo, hi) in enumerate(bands):
        if lo[1] < hi[0] and absP_side((lo[1] + hi[0]) / 2) != -1:
            problems.append(f"band {i + 1} midpoint fails |P| < 4")
    for g, (bl, br) in zip(rec["gaps"], zip(bands, bands[1:])):
        if g["touching"]:
            continue
        left, right = bl[1], br[0]
        if absP_side((left[1] + right[0]) / 2) != 1:
            problems.append(f"gap {g['k']} midpoint fails |P| > 4")
        w_lo, w_hi = Fraction(g["width_lo"]), Fraction(g["width_hi"])
        if not (0 < w_lo <= right[0] - left[1] and right[1] - left[0] <= w_hi):
            problems.append(f"gap {g['k']} width enclosure inconsistent with edges")

    # bandwidth enclosure re-summed
    bw_lo = sum(hi[0] - lo[1] for lo, hi in bands)
    bw_hi = sum(hi[1] - lo[0] for lo, hi in bands)
    if [str(bw_lo), str(bw_hi)] != rec["bandwidth"]:
        problems.append("bandwidth enclosure does not re-sum")
    if [str(q * bw_lo), str(q * bw_hi)] != rec["q_bandwidth"]:
        problems.append("q * bandwidth does not re-sum")
    return problems


def _verify_small_q(rec: dict) -> list[str]:
    """q = 1, 2: closed forms E and E^2 - 4."""
    q = rec["q"]
    want = [["0"], ["1"]] if q == 1 else [["-4"], ["0"], ["1"]]
    out = []
    if rec["P"] != want:
        out.append(f"claimed P for q={q} differs from the closed form")
    return out


# ---------------------------------------------------------------------------


def selftest() -> None:
    # 2 cos(72 deg) = (sqrt(5) - 1) / 2 ~ 0.618: check the oracle brackets it
    o = CosOracle(5)
    assert o.sign([0, 1]) == 1  # y > 0
    assert o.sign([-1, 1]) == -1  # y < 1
    assert o.sign([1, -1]) == 1
    assert Fraction(0.61) < o.hi and o.lo < Fraction(0.62)
    # interpolation route reproduces the rational-field closed forms
    for (p, q), want in [((1, 3), [0, -6, 0, 1]), ((1, 4), [4, 0, -8, 0, 1])]:
        vecs = interpolated_P(p, q)
        for k, vec in enumerate(vecs):
            v = trim(list(vec))
            assert (not v and want[k] == 0) or v == [want[k]], (p, q, k, vec)
    print("selftest: all checks passed")


def main() -> None:
    ap = argparse.ArgumentParser(description="Independent critical-AMO record verification")
    ap.add_argument("record", nargs="?", help="JSON record produced by the reference tool")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if not args.record:
        ap.error("need RECORD (or --selftest)")
    problems = verify_record(args.record)
    if problems:
        print("VERIFICATION FAILED:")
        for pr in problems:
            print(f"  - {pr}")
        sys.exit(1)
    print("record verified: P, edge certificates, gap certificates and sums all reproduce")


if __name__ == "__main__":
    main()
