#!/usr/bin/env python3
"""Certified band/gap data for the critical almost Mathieu operator at rational flux.

For flux alpha = p/q (lowest terms) and critical coupling lambda = 1, Chambers'
relation says the discriminant of the period-q potential
v_j = 2 cos(2 pi p j / q + nu) has the form

    Delta_nu(E) = P(E) - c * cos(q nu)          (single harmonic in nu),

with P and the constant c independent of nu, so the spectrum is

    sigma(1, p/q) = { E : |P(E)| <= 2 + |c| },   and  |c| = 2  here.

The coefficients of P are NOT rational integers in general: they lie in
Z[2 cos(2 pi p / q)], the ring of integers of the real cyclotomic field
Q(zeta_q)^+ -- rational exactly when phi(q) <= 2, i.e. q in {1, 2, 3, 4, 6}.
(Flux 1/5 already has P(E) = E^5 - 10 E^3 + (15 - 10 cos 72deg) E.)

Everything here is exact:

  * P is computed in K = Q[y]/(psi_q(y)), y = 2 cos(2 pi / q), psi_q the
    minimal polynomial of y (extracted exactly from the cyclotomic polynomial
    Phi_q).  The potential at nu = 0 is v_j = D_{pj mod q}(y) with D_k the
    Chebyshev-style recurrence D_0 = 2, D_1 = y, D_k = y D_{k-1} - D_{k-2},
    so the transfer-matrix product lives in K[E] with integer vectors.
  * The nu-structure is pinned by two exact cross-phase identities: the same
    transfer product is recomputed in Z[x]/(x^{4q}-1) at nu = pi/(2q) (where
    cos(q nu) = 0) and in Z[x]/(x^{2q}-1) at nu = pi/q (where cos(q nu) = -1),
    both are compared with the embedded K result, and the differences are
    asserted to be the constants s*2 and s*4 for one sign s.  That certifies
    P = Delta_0 + s*2 and threshold |P| <= 4 given Chambers' single-harmonic
    theorem, whatever the sign convention.
  * Band edges are the real roots of P -+ 4.  Rational candidates are removed
    by exact deflation (this includes the known touching at E = 0 for even
    q); the remaining simple roots are bracketed by exact sign changes found
    by bisection.  q disjoint sign-change brackets for a degree-q polynomial
    is a complete root count -- no root can hide.
  * Signs of elements of K are decided exactly: zero by the canonical
    representation, nonzero by interval arithmetic over a refinable rational
    enclosure of y (the largest root of psi_q, isolated by a rational Sturm
    chain).  No floating point participates in any certificate; floats appear
    only in the --selftest screen and display columns, both labelled.

Usage:
    python chambers.py P Q [--bits B] [--json PATH]
    python chambers.py --selftest

Standard library only.  Reference implementation: clarity over speed -- large
censuses want a C kernel built to the same contract.
"""

from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from math import gcd

# ---------------------------------------------------------------------------
# integer polynomials (dense, index = degree)
# ---------------------------------------------------------------------------


def ipoly_trim(f: list[int]) -> list[int]:
    while f and f[-1] == 0:
        f.pop()
    return f


def ipoly_divmod(f: list[int], g: list[int]) -> tuple[list[int], list[int]]:
    """Exact division of integer polynomials; g must be monic."""
    assert g and g[-1] == 1, "divisor must be monic"
    f = list(f)
    dg = len(g) - 1
    quo = [0] * max(1, len(f) - dg)
    while ipoly_trim(f) and len(f) - 1 >= dg:
        d = len(f) - 1 - dg
        c = f[-1]
        quo[d] = c
        for i, gi in enumerate(g):
            f[d + i] -= c * gi
        ipoly_trim(f)
    return ipoly_trim(quo), ipoly_trim(f)


_CYCLO_CACHE: dict[int, list[int]] = {}


def cyclotomic(n: int) -> list[int]:
    """Phi_n as an integer coefficient list, by the recursive quotient formula."""
    if n in _CYCLO_CACHE:
        return _CYCLO_CACHE[n]
    f = [-1] + [0] * (n - 1) + [1]  # x^n - 1
    for d in range(1, n):
        if n % d == 0:
            f, rem = ipoly_divmod(f, cyclotomic(d))
            assert not rem, f"cyclotomic division not exact at n={n}, d={d}"
    _CYCLO_CACHE[n] = f
    return f


def minpoly_2cos(q: int) -> list[int]:
    """Minimal polynomial psi_q of y = 2 cos(2 pi / q), for q >= 3.

    Phi_q is palindromic of degree 2d, d = phi(q)/2, and
    x^d psi_q(x + 1/x) = Phi_q(x); the coefficients are extracted by matching
    degrees downward against T_i(x) = x^{d-i} (x^2+1)^i.
    """
    assert q >= 3
    phi_q = cyclotomic(q)
    d = (len(phi_q) - 1) // 2
    rest = list(phi_q)
    psi = [0] * (d + 1)
    for i in range(d, -1, -1):
        c = rest[d + i] if d + i < len(rest) else 0
        psi[i] = c
        if c:
            # subtract c * x^{d-i} * (x^2+1)^i
            ti = [0] * (d - i) + _binomial_x2_plus_1(i)
            for k, tk in enumerate(ti):
                rest[k] -= c * tk
    assert not ipoly_trim(rest), f"psi extraction left a remainder for q={q}"
    assert psi[d] == 1
    return psi


def _binomial_x2_plus_1(i: int) -> list[int]:
    """(x^2 + 1)^i as a coefficient list."""
    out = [1]
    for _ in range(i):
        nxt = [0] * (len(out) + 2)
        for k, c in enumerate(out):
            nxt[k] += c
            nxt[k + 2] += c
        out = nxt
    return out


# ---------------------------------------------------------------------------
# rational polynomials + Sturm (used only to isolate the largest root of psi)
# ---------------------------------------------------------------------------


def qp_trim(f: list[Fraction]) -> list[Fraction]:
    while f and f[-1] == 0:
        f.pop()
    return f


def qp_eval(f: list[Fraction], x: Fraction) -> Fraction:
    acc = Fraction(0)
    for c in reversed(f):
        acc = acc * x + c
    return acc


def qp_divmod(f: list[Fraction], g: list[Fraction]) -> tuple[list[Fraction], list[Fraction]]:
    f = list(f)
    g = qp_trim(list(g))
    assert g, "division by zero polynomial"
    dg = len(g) - 1
    quo = [Fraction(0)] * max(1, len(f) - dg)
    while qp_trim(f) and len(f) - 1 >= dg:
        d = len(f) - 1 - dg
        c = f[-1] / g[-1]
        quo[d] = c
        for i, gi in enumerate(g):
            f[d + i] -= c * gi
        qp_trim(f)
    return qp_trim(quo), qp_trim(f)


def sturm_chain(f: list[Fraction]) -> list[list[Fraction]]:
    deriv = [c * k for k, c in enumerate(f)][1:] or [Fraction(0)]
    chain = [qp_trim(list(f)), qp_trim(deriv)]
    while len(chain[-1]) > 1:
        _, r = qp_divmod(chain[-2], chain[-1])
        if not r:
            break
        chain.append([-c for c in r])
    return chain


def sign_variations(chain: list[list[Fraction]], x: Fraction) -> int:
    signs = []
    for f in chain:
        v = qp_eval(f, x)
        if v != 0:
            signs.append(1 if v > 0 else -1)
    return sum(1 for a, b in zip(signs, signs[1:]) if a != b)


def largest_root_bracket(psi: list[int]) -> tuple[Fraction, Fraction]:
    """An isolating rational interval for the largest real root of psi.

    All roots of psi_q are 2 cos(2 pi k / q) in (-2, 2); the largest is
    y = 2 cos(2 pi / q).
    """
    f = [Fraction(c) for c in psi]
    if len(f) == 2:  # linear: exact rational root
        r = -f[0] / f[1]
        return r, r
    chain = sturm_chain(f)
    lo, hi = Fraction(-2), Fraction(2)
    # peel off the rightmost root by bisection on Sturm counts
    while True:
        mid = (lo + hi) / 2
        if sign_variations(chain, mid) - sign_variations(chain, hi) >= 1:
            lo = mid
        else:
            hi = mid
        if sign_variations(chain, lo) - sign_variations(chain, hi) == 1 and (
            qp_eval(f, lo) != 0 and qp_eval(f, hi) != 0
        ):
            return lo, hi


# ---------------------------------------------------------------------------
# K = Q[y]/(psi_q): vectors of length d, plus an exact sign oracle
# ---------------------------------------------------------------------------


class Field:
    """The real cyclotomic field Q(2 cos(2 pi / q)) with exact sign decisions."""

    def __init__(self, q: int):
        assert q >= 3
        self.q = q
        self.psi = minpoly_2cos(q)
        self.d = len(self.psi) - 1
        self.lo, self.hi = largest_root_bracket(self.psi)

    # -- ring operations (vectors of length d; int or Fraction entries) -----

    def zero(self):
        return [0] * self.d

    def one(self):
        e = [0] * self.d
        e[0] = 1
        return e

    def scalar(self, c):
        e = [0] * self.d
        e[0] = c
        return e

    def add(self, a, b):
        return [x + y for x, y in zip(a, b)]

    def sub(self, a, b):
        return [x - y for x, y in zip(a, b)]

    def scale(self, a, c):
        return [x * c for x in a]

    def mul(self, a, b):
        conv = [0] * (2 * self.d - 1)
        for i, ai in enumerate(a):
            if ai:
                for j, bj in enumerate(b):
                    if bj:
                        conv[i + j] += ai * bj
        # reduce mod psi (monic)
        for k in range(len(conv) - 1, self.d - 1, -1):
            c = conv[k]
            if c:
                conv[k] = 0
                for i, pi in enumerate(self.psi[:-1]):
                    conv[k - self.d + i] -= c * pi
        return conv[: self.d]

    def gen_power_2cos(self, k: int):
        """D_k(y) = 2 cos(2 pi k / q) as an element of K (Chebyshev recurrence)."""
        k %= self.q
        a, b = self.scalar(2), None  # D_0, then D_1 = y
        if k == 0:
            return a
        y = self.zero()
        if self.d >= 2:
            y[1] = 1
        else:
            y[0] = -self.psi[0]  # linear psi: y is rational
        b = y
        for _ in range(k - 1):
            a, b = b, self.sub(self.mul(y, b), a)
        return b

    # -- exact sign decisions ------------------------------------------------

    def refine_bracket(self) -> None:
        f = [Fraction(c) for c in self.psi]
        mid = (self.lo + self.hi) / 2
        vm = qp_eval(f, mid)
        assert vm != 0, "psi has a rational root above degree 1?"
        # keep the sub-interval with the sign change (largest root already isolated)
        if (qp_eval(f, self.lo) > 0) == (vm > 0):
            self.lo = mid
        else:
            self.hi = mid

    def sign(self, a) -> int:
        """Exact sign of a in K under y -> 2 cos(2 pi / q)."""
        if all(x == 0 for x in a):
            return 0
        if self.d == 1 or self.lo == self.hi:
            v = sum(Fraction(x) * self.lo**i for i, x in enumerate(a))
            assert v != 0, "nonzero element evaluated to zero at rational generator"
            return 1 if v > 0 else -1
        while True:
            lo, hi = Fraction(0), Fraction(0)
            for c in reversed(a):
                cands = (lo * self.lo, lo * self.hi, hi * self.lo, hi * self.hi)
                lo, hi = min(cands) + c, max(cands) + c
            if lo > 0:
                return 1
            if hi < 0:
                return -1
            self.refine_bracket()

    def approx(self, a) -> float:
        """Display only: float value of a."""
        y = 2 * math.cos(2 * math.pi / self.q)
        return float(sum(float(x) * y**i for i, x in enumerate(a)))


# ---------------------------------------------------------------------------
# transfer-matrix products in three exact models of the discriminant
# ---------------------------------------------------------------------------


def transfer_K(field: Field, p: int, q: int) -> list[list]:
    """Delta_0(E) in K[E]: potential v_j = D_{pj mod q}(y), phase nu = 0."""
    t00, t01 = [field.one()], [field.zero()]
    t10, t11 = [field.zero()], [field.one()]
    for j in range(1, q + 1):
        vj = field.gen_power_2cos(p * j)

        def step(top, bot):
            out = [field.zero()] + [list(c) for c in top]  # E * top
            for k, c in enumerate(top):
                out[k] = field.sub(out[k], field.mul(c, vj))
            for k, c in enumerate(bot):
                out[k] = field.sub(out[k], c)
            return out

        n00, n01 = step(t00, t10), step(t01, t11)
        t10, t11 = t00, t01
        t00, t01 = n00, n01
    deg = max(len(t00), len(t11)) - 1
    out = []
    for k in range(deg + 1):
        c = field.zero()
        if k < len(t00):
            c = field.add(c, t00[k])
        if k < len(t11):
            c = field.add(c, t11[k])
        out.append(c)
    return out


def transfer_cyclic(p: int, q: int, N: int, exps: list[int]) -> list[list[int]]:
    """Discriminant in Z[x]/(x^N - 1) with v_j = x^{exps[j]} + x^{-exps[j]}."""
    zero = [0] * N

    def mul_sparse(a, e):
        out = [0] * N
        for ee in (e % N, (N - e) % N):
            for k, ak in enumerate(a):
                if ak:
                    out[(k + ee) % N] += ak
        return out

    one = zero[:]
    one[0] = 1
    t00, t01 = [one[:]], [zero[:]]
    t10, t11 = [zero[:]], [one[:]]
    for j in range(1, q + 1):
        e = exps[j - 1]

        def step(top, bot):
            out = [zero[:]] + [c[:] for c in top]
            for k, c in enumerate(top):
                m = mul_sparse(c, e)
                out[k] = [x - y for x, y in zip(out[k], m)]
            for k, c in enumerate(bot):
                out[k] = [x - y for x, y in zip(out[k], c)]
            return out

        n00, n01 = step(t00, t10), step(t01, t11)
        t10, t11 = t00, t01
        t00, t01 = n00, n01
    deg = max(len(t00), len(t11)) - 1
    out = []
    for k in range(deg + 1):
        c = zero[:]
        if k < len(t00):
            c = [x + y for x, y in zip(c, t00[k])]
        if k < len(t11):
            c = [x + y for x, y in zip(c, t11[k])]
        out.append(c)
    return out


def embed_K_in_cyclic(field: Field, a: list, N: int) -> list:
    """y -> x^4 + x^{N-4} in Z[x]/(x^N - 1) (valid for N = 4q)."""
    ybig = [0] * N
    ybig[4 % N] += 1
    ybig[(N - 4) % N] += 1
    out = [0] * N
    power = [0] * N
    power[0] = 1  # y^0
    for i, c in enumerate(a):
        if c:
            out = [o + c * pw for o, pw in zip(out, power)]
        if i + 1 < len(a):
            nxt = [0] * N
            for k, pk in enumerate(power):
                if pk:
                    for m, ym in enumerate(ybig):
                        if ym:
                            nxt[(k + m) % N] += pk * ym
            power = nxt
    return out


def cyclic_reduce(a: list, N: int) -> list:
    """Canonical representative mod Phi_N (rational vector allowed)."""
    f = [Fraction(c) for c in a]
    g = [Fraction(c) for c in cyclotomic(N)]
    _, rem = qp_divmod(f, g)
    return rem


def chambers_poly_K(p: int, q: int) -> tuple[Field, list]:
    """P in K[E] with sigma(1, p/q) = {E : |P(E)| <= 4}, exactly certified.

    Given Chambers' single-harmonic theorem Delta_nu = P - c cos(q nu), the
    three exact evaluations below overdetermine (P, c) and are asserted
    consistent; the sign convention is detected, not assumed.
    """
    assert q >= 3 and 1 <= p < q and gcd(p, q) == 1
    field = Field(q)
    delta0 = transfer_K(field, p, q)  # P - c
    # nu = pi/(2q): cos(q nu) = 0 -> P itself, in Z[x]/(x^{4q}-1)
    d_mid = transfer_cyclic(p, q, 4 * q, [(4 * p * j + 1) % (4 * q) for j in range(1, q + 1)])
    # nu = pi/q: cos(q nu) = -1 -> P + c, in Z[x]/(x^{2q}-1)
    d_pi = transfer_cyclic(p, q, 2 * q, [(2 * p * j + 1) % (2 * q) for j in range(1, q + 1)])

    N = 4 * q
    assert len(delta0) == len(d_mid) == len(d_pi) == q + 1
    sign_s = None
    for k in range(q + 1):
        e0 = cyclic_reduce(embed_K_in_cyclic(field, delta0[k], N), N)
        em = cyclic_reduce(d_mid[k], N)
        # spread the 2q-ring vector into the 4q ring: x_{2q} = x_{4q}^2
        spread = [0] * N
        for i, c in enumerate(d_pi[k]):
            spread[(2 * i) % N] += c
        ep = cyclic_reduce(spread, N)
        dm = qp_trim([x - y for x, y in zip(em + [Fraction(0)] * N, e0 + [Fraction(0)] * N)])
        dp = qp_trim([x - y for x, y in zip(ep + [Fraction(0)] * N, e0 + [Fraction(0)] * N)])
        if k < q:
            want = Fraction(2) if k == 0 else None
        if k == 0:
            assert len(dm) == 1 and abs(dm[0]) == 2, f"mid-phase identity broke: {dm}"
            assert len(dp) == 1 and abs(dp[0]) == 4, f"pi-phase identity broke: {dp}"
            s_m = 1 if dm[0] > 0 else -1
            s_p = 1 if dp[0] > 0 else -1
            assert s_m == s_p, "phase identities disagree on the sign convention"
            sign_s = s_m
        else:
            assert not dm and not dp, f"phase identities have E^{k} residue: {dm} {dp}"

    # P = Delta_0 + s*2
    P = [list(c) for c in delta0]
    P[0] = field.add(P[0], field.scalar(2 * sign_s))
    assert field.sign(field.sub(P[-1], field.one())) == 0, "P is not monic"
    return field, P


# ---------------------------------------------------------------------------
# band / gap certification over K
# ---------------------------------------------------------------------------


def keval_hom(field: Field, f: list, num: int, den: int) -> list:
    """den^deg(f) * f(num/den) as an integer K-vector (den > 0: sign-preserving).

    All coefficient vectors of f must be integer vectors, which they are
    throughout: P is integral over Z[y] and deflation only happens at integer
    roots.  Pure integer arithmetic -- this is what keeps large q feasible.
    """
    acc = field.zero()
    pw = 1
    for c in reversed(f):
        acc = [x * num + ci * pw for x, ci in zip(acc, c)]
        pw *= den
    return acc


def sgn_at(field: Field, f: list, x: Fraction) -> int:
    return field.sign(keval_hom(field, f, x.numerator, x.denominator))


def kdeflate(field: Field, f: list, r: int) -> list:
    """Exact synthetic division of f in K[E] by (E - r), integer r."""
    partials = []
    acc = field.zero()
    for c in reversed(f):
        acc = [x * r + ci for x, ci in zip(acc, c)]
        partials.append(acc)
    assert field.sign(partials[-1]) == 0, "deflation at a non-root"
    return list(reversed(partials[:-1]))


def coeff_bound(field: Field, f: list) -> Fraction:
    """A rational M with all real roots of f in (-M, M) (Cauchy, monic f)."""
    worst = Fraction(0)
    for c in f[:-1]:
        # |c| <= sum |c_i| * 2^i since |y| < 2
        bound = sum(abs(Fraction(x)) * 2**i for i, x in enumerate(c))
        worst = max(worst, bound)
    return worst + 1


def _float_root_seeds(field: Field, f: list, M: Fraction) -> list[float]:
    """Screen only: approximate real roots of f by float sampling + bisection.

    Feeds the exact certification below with candidate separators; never
    certifies anything itself.
    """
    coeffs = [field.approx(c) for c in f]

    def ev(x: float) -> float:
        acc = 0.0
        for c in reversed(coeffs):
            acc = acc * x + c
        return acc

    lo, hi = float(-M), float(M)
    n = max(64, 40 * (len(f) - 1))
    xs = [lo + (hi - lo) * i / n for i in range(n + 1)]
    vals = [ev(x) for x in xs]
    roots = []
    for (x1, v1), (x2, v2) in zip(zip(xs, vals), zip(xs[1:], vals[1:])):
        if v1 == 0.0:
            roots.append(x1)
            continue
        if v1 * v2 < 0:
            a, b, va = x1, x2, v1
            for _ in range(60):
                m = (a + b) / 2
                vm = ev(m)
                if vm == 0.0:
                    a = b = m
                    break
                if (vm > 0) == (va > 0):
                    a, va = m, vm
                else:
                    b = m
            roots.append((a + b) / 2)
    return roots


def find_brackets(field: Field, f: list, count: int) -> list[tuple[Fraction, Fraction]]:
    """Exactly `count` disjoint sign-change brackets for the simple roots of f.

    Complete by the degree argument: deg f = count, so `count` certified
    disjoint sign changes leave no room for missed roots.  Roots are located
    by a float screen and certified by exact signs at rational separators; if
    the screen misses (near-degenerate roots), an exact adaptive subdivision
    takes over.
    """
    # every band edge lies in [-4, 4] (the spectrum sits inside it), so a
    # fixed M = 5 keeps the search range tight; a root outside would surface
    # as a stalled count below, never as a silent omission
    M = min(coeff_bound(field, f), Fraction(5))
    seeds = _float_root_seeds(field, f, M)
    seps: list[Fraction] = [-M]
    for r1, r2 in zip(seeds, seeds[1:]):
        seps.append(Fraction((r1 + r2) / 2).limit_denominator(10**15))
    seps.append(M)

    pts = []
    for x in seps:
        s = sgn_at(field, f, x)
        if s == 0:  # a separator hit a root exactly; nudge it
            x = x + Fraction(1, 10**9)
            s = sgn_at(field, f, x)
            assert s != 0
        pts.append((x, s))
    pts.sort()

    rounds = 0
    while True:
        brackets = [(a, b) for (a, sa), (b, sb) in zip(pts, pts[1:]) if sa != sb]
        if len(brackets) == count:
            return brackets
        assert len(brackets) < count, "more sign changes than roots?"
        rounds += 1
        assert rounds <= 60, (
            "bracket search stalled: unexpected multiple or near-degenerate roots"
        )
        # split every interval; the range is tight, so the depth needed to
        # separate the closest pair stays small
        nxt = [pts[0]]
        for (a, sa), (b, sb) in zip(pts, pts[1:]):
            mid = (a + b) / 2
            sm = sgn_at(field, f, mid)
            assert sm != 0, f"unexpected rational root at {mid}"
            nxt.append((mid, sm))
            nxt.append((b, sb))
        pts = nxt


def refine_bracket(field: Field, f: list, a: Fraction, b: Fraction, width: Fraction):
    sa = sgn_at(field, f, a)
    while b - a > width:
        mid = (a + b) / 2
        sm = sgn_at(field, f, mid)
        assert sm != 0, "unexpected rational root during refinement"
        if sm == sa:
            a = mid
        else:
            b = mid
    return a, b


class Edge:
    """A band edge: exact rational root (lo == hi) or a sign-change bracket."""

    __slots__ = ("lo", "hi", "mult", "poly")

    def __init__(self, lo, hi, mult, poly):
        self.lo, self.hi, self.mult, self.poly = lo, hi, mult, poly

    def shrink(self, field: Field) -> None:
        if self.poly is not None:
            self.lo, self.hi = refine_bracket(field, self.poly, self.lo, self.hi, (self.hi - self.lo) / 4)


def spectrum(p: int, q: int, bits: int = 30) -> dict:
    """Certified bands, gaps and bandwidth enclosure for sigma(1, p/q)."""
    if q <= 2:
        # closed forms: P = E (q=1), P = E^2 - 4 (q=2); handled via K on q=4's field
        return _spectrum_small_q(p, q, bits)
    field, P = chambers_poly_K(p, q)
    width = Fraction(1, 2**bits)
    edges: list[Edge] = []
    for side in (+1, -1):
        f = [list(c) for c in P]
        f[0] = field.sub(f[0], field.scalar(4 * side))
        count = 0
        for r in range(-4, 5):  # all edges lie in [-4, 4]; rational roots are integers
            mult = 0
            while len(f) > 1 and field.sign(keval_hom(field, f, r, 1)) == 0:
                f = kdeflate(field, f, r)
                mult += 1
            if mult:
                edges.append(Edge(Fraction(r), Fraction(r), mult, None))
                count += mult
        if len(f) > 1:
            rest = len(f) - 1
            for a, b in find_brackets(field, f, rest):
                lo, hi = refine_bracket(field, f, a, b, width)
                edges.append(Edge(lo, hi, 1, f))
                count += 1
        assert count == q, f"P {side:+d}4 edge count {count} != {q}"

    changed = True
    while changed:
        changed = False
        edges.sort(key=lambda e: (e.lo, e.hi))
        for u, v in zip(edges, edges[1:]):
            if v.lo < u.hi and not (u.lo == v.lo and u.hi == v.hi):
                u.shrink(field)
                v.shrink(field)
                changed = True

    flat: list[Edge] = []
    for e in sorted(edges, key=lambda e: (e.lo, e.hi)):
        flat.extend([e] * e.mult)
    assert len(flat) == 2 * q, f"expected {2 * q} band edges, got {len(flat)}"

    def absP_side(x: Fraction) -> int:
        """+1 if |P(x)| > 4, -1 if < 4, else 0 (exact)."""
        v = keval_hom(field, P, x.numerator, x.denominator)
        scale = 4 * x.denominator ** (len(P) - 1)
        hi_ = field.sign(field.sub(v, field.scalar(scale)))
        lo_ = field.sign(field.add(v, field.scalar(scale)))
        if hi_ > 0 or lo_ < 0:
            return 1
        if hi_ < 0 and lo_ > 0:
            return -1
        return 0

    bands = [(flat[2 * i], flat[2 * i + 1]) for i in range(q)]
    for i, (lo_e, hi_e) in enumerate(bands):
        if lo_e.hi < hi_e.lo:
            assert absP_side((lo_e.hi + hi_e.lo) / 2) == -1, f"band {i + 1} interior check failed"

    pinv = pow(p, -1, q)
    gaps = []
    for k in range(1, q):
        left, right = bands[k - 1][1], bands[k][0]
        if left is right:
            gaps.append({"k": k, "touching": True, "width_lo": Fraction(0), "width_hi": Fraction(0)})
            continue
        w_lo = right.lo - left.hi
        tries = 0
        while w_lo <= 0 and tries < 200:
            left.shrink(field)
            right.shrink(field)
            w_lo = right.lo - left.hi
            tries += 1
        assert w_lo > 0, f"gap {k} not certified open after refinement"
        assert absP_side((left.hi + right.lo) / 2) == 1, f"gap {k} interior check failed"
        ell = (k * pinv) % q
        if ell > q // 2:
            ell -= q
        gaps.append(
            {"k": k, "touching": False, "width_lo": w_lo, "width_hi": right.hi - left.lo, "label": ell}
        )

    bw_lo = sum(hi.lo - lo.hi for lo, hi in bands)
    bw_hi = sum(hi.hi - lo.lo for lo, hi in bands)
    return {
        "p": p,
        "q": q,
        "lambda": "1",
        "threshold": 4,
        "bits": bits,
        "psi": field.psi,
        "P": [[str(x) for x in c] for c in P],
        "P_float_display": [field.approx(c) for c in P],
        "bands": [
            {"lo": [str(lo.lo), str(lo.hi)], "hi": [str(hi.lo), str(hi.hi)]} for lo, hi in bands
        ],
        "gaps": [
            {
                "k": g["k"],
                "touching": g["touching"],
                "width_lo": str(g["width_lo"]),
                "width_hi": str(g["width_hi"]),
                **({"label": g["label"]} if "label" in g else {}),
            }
            for g in gaps
        ],
        "bandwidth": [str(bw_lo), str(bw_hi)],
        "q_bandwidth": [str(q * bw_lo), str(q * bw_hi)],
    }


def _spectrum_small_q(p: int, q: int, bits: int) -> dict:
    """q = 1, 2 exactly: P = E and P = E^2 - 4 (from the same transfer algebra)."""
    width = Fraction(1, 2**bits)
    if q == 1:
        bands = [(Fraction(-4), Fraction(-4), Fraction(4), Fraction(4))]
        P_repr = [["0"], ["1"]]
        Pf = [0.0, 1.0]
        gaps = []
    else:
        # edges: roots of E^2 - 8 (outer) and E^2 = 0 (touching)
        lo, hi = Fraction(2), Fraction(3)
        f = [Fraction(-8), Fraction(0), Fraction(1)]
        while hi - lo > width:
            mid = (lo + hi) / 2
            if qp_eval(f, mid) < 0:
                lo = mid
            else:
                hi = mid
        bands = [(-hi, -lo, Fraction(0), Fraction(0)), (Fraction(0), Fraction(0), lo, hi)]
        P_repr = [["-4"], ["0"], ["1"]]
        Pf = [-4.0, 0.0, 1.0]
        gaps = [{"k": 1, "touching": True, "width_lo": "0", "width_hi": "0"}]
    bw_lo = sum(b[2] - b[1] for b in bands)
    bw_hi = sum(b[3] - b[0] for b in bands)
    return {
        "p": p,
        "q": q,
        "lambda": "1",
        "threshold": 4,
        "bits": bits,
        "psi": [1],
        "P": P_repr,
        "P_float_display": Pf,
        "bands": [{"lo": [str(b[0]), str(b[1])], "hi": [str(b[2]), str(b[3])]} for b in bands],
        "gaps": gaps,
        "bandwidth": [str(bw_lo), str(bw_hi)],
        "q_bandwidth": [str(q * bw_lo), str(q * bw_hi)],
    }


# ---------------------------------------------------------------------------
# self-test (floats appear here ONLY, as a screen -- never as a certificate)
# ---------------------------------------------------------------------------


def _float_disc(p: int, q: int, nu: float, E: float) -> float:
    t = [[1.0, 0.0], [0.0, 1.0]]
    for j in range(1, q + 1):
        v = 2 * math.cos(2 * math.pi * p * j / q + nu)
        t = [
            [(E - v) * t[0][0] - t[1][0], (E - v) * t[0][1] - t[1][1]],
            [t[0][0], t[0][1]],
        ]
    return t[0][0] + t[1][1]


def selftest() -> None:
    # rational-field fluxes have integer Chambers polynomials
    for (p_, q_), want in [
        ((1, 3), [0, -6, 0, 1]),
        ((2, 3), [0, -6, 0, 1]),
        ((1, 4), [4, 0, -8, 0, 1]),
        ((1, 6), [-4, 0, 24, 0, -12, 0, 1]),
    ]:
        field, P = chambers_poly_K(p_, q_)
        got = [c[0] if field.d == 1 else None for c in P]
        assert field.d == 1 and [Fraction(w) for w in want] == [Fraction(g) for g in got], (
            f"flux {p_}/{q_}: expected {want}, got {got}"
        )
    # irrational-field fluxes: exact P against a float transfer-matrix screen
    for p_, q_ in [(1, 5), (2, 5), (1, 7), (3, 7), (1, 8), (3, 8)]:
        field, P = chambers_poly_K(p_, q_)
        for E in (-1.7, 0.3, 2.1):
            exact = sum(field.approx(c) * E**k for k, c in enumerate(P))
            screen = _float_disc(p_, q_, math.pi / (2 * q_), E)
            assert abs(exact - screen) < 1e-6, (
                f"flux {p_}/{q_} at E={E}: exact {exact} vs screen {screen}"
            )
    # flux 1/2 closed form: touching at 0, outer edges at 2*sqrt(2)
    s = spectrum(1, 2, bits=40)
    assert s["gaps"][0]["touching"]
    assert abs(float(Fraction(s["bands"][1]["hi"][1])) - 2 * math.sqrt(2)) < 1e-9
    # structural facts: q bands, all gaps open except central touching iff q even
    for q_ in range(3, 9):
        for p_ in [p for p in range(1, q_) if gcd(p, q_) == 1]:
            s = spectrum(p_, q_, bits=20)
            assert len(s["bands"]) == q_
            touches = [g["k"] for g in s["gaps"] if g["touching"]]
            assert touches == ([q_ // 2] if q_ % 2 == 0 else []), (
                f"unexpected touchings {touches} at flux {p_}/{q_}"
            )
            for g in s["gaps"]:
                assert g["touching"] or Fraction(g["width_lo"]) > 0
    print("selftest: all checks passed")


# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description="Exact critical-AMO band data at rational flux")
    ap.add_argument("p", type=int, nargs="?")
    ap.add_argument("q", type=int, nargs="?")
    ap.add_argument("--bits", type=int, default=30, help="edge enclosure width 2^-bits")
    ap.add_argument("--json", metavar="PATH", help="write the full certified record here")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return
    if args.p is None or args.q is None:
        ap.error("need flux numerator and denominator (or --selftest)")

    s = spectrum(args.p, args.q, bits=args.bits)
    catalan = sum((-1) ** k / (2 * k + 1) ** 2 for k in range(200000))
    bw = [Fraction(x) for x in s["q_bandwidth"]]
    print(f"flux {args.p}/{args.q}, lambda = 1 (critical), {args.q} bands")
    open_gaps = [g for g in s["gaps"] if not g["touching"]]
    touch = [g["k"] for g in s["gaps"] if g["touching"]]
    print(f"gaps: {len(open_gaps)} certified open" + (f", touching at k={touch}" if touch else ""))
    if open_gaps:
        wmin = min(Fraction(g["width_lo"]) for g in open_gaps)
        print(f"smallest certified gap width >= {float(wmin):.3e}")
    print(f"q * |sigma|  in  [{float(bw[0]):.9f}, {float(bw[1]):.9f}]")
    print(f"  (display only: Thouless constant 32C/pi = {32 * catalan / math.pi:.9f})")
    if args.json:
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(s, fh, indent=1)
        print(f"record written to {args.json}")


if __name__ == "__main__":
    main()
