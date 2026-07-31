#!/usr/bin/env python3
"""Attempt 004: skeptic review of 003's death-angle laws.  Default REFUTE.

Everything here is written from scratch for this review; nothing is imported
from deathlaw_{measure,exact,symbolic,prove}.py (the code under review).
skeptic_family (002's independent float corridor) is imported ONLY inside the
selftest as a cross-reference, never in a verdict path.

Independent implementation choices, deliberately different from 003's:

  * Unfolding by COMPOSED ISOMETRY MAPS applied to the base triangle's sides
    (gate k = M_{k-1}(base side s_k), M_k = M_{k-1} o R_{s_k}), instead of
    003's reflected-vertex chains.
  * The glide projection action p(u(z)) = 2m - p(z) and the axis offset
    m = Im(conj(delta) w)/2 are re-derived by hand (see attempt record) and
    then verified as ring identities on all base vertices, rather than taken
    from 003's decomposition.
  * The ring identities I1-I3/D1/D2 are checked twice: (a) in my own Laurent
    ring; (b) proof-grade for W(2,1) by exact rational evaluation OFF the
    unit torus, using a formal-conjugate pair representation
    (value, star-value) at a grid of rational points exceeding an a-priori
    degree bound (interpolation argument; see record), plus seeded random
    rational spot points for all 15 members.
  * The interval stack for Lemma C and the gamma certifications is my own:
    pi by Machin's formula with alternating-series remainder bounds, sin/cos
    by my own Taylor enclosures.  unfold.py's intervals are NOT used.
  * Death measurement uses my own adaptive arc sampler and my own float
    corridor (composed maps, complex arithmetic).

Subcommands (all deterministic; seeds fixed):
  selftest    base geometry + corridor cross-checks
  identities  ring proof of I1,I2,I3,D1,D2 + glide facts, my ring, 15 members
  gridproof   off-torus grid certification for one member (default W(2,1))
  reduction   numeric check: full corridor == half-corridor cap glide image;
              positive width => m strictly inside every half-gate projection
  casetree    adversarial search for a sign-pattern counterexample to the
              theorem's case analysis at theta <= theta_d
  alivecheck  re-verify 003's exact alive certificates (own exact corridor,
              own certified gamma enclosures)
  death       my own adaptive death bisection for one member
  deadhunt    search for alive apexes AT/ABOVE gamma_d (kill attempt)
  lemmac      my own certification of Lemma C for one member

Usage examples (repo root):
  python3 problems/billiards-triangles/explore/deathlaw_skeptic.py identities \
      --out problems/billiards-triangles/data/dlsk_identities.json
"""

import argparse
import cmath
import json
import math
import os
import random
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

SIDE_ENDS = {0: (1, 2), 1: (2, 0), 2: (0, 1)}
MEMBERS = [(2, 1), (2, 2), (3, 2), (3, 3), (4, 2), (4, 3), (4, 4),
           (5, 3), (5, 4), (5, 5), (6, 6), (7, 7), (8, 8),
           (10, 9), (12, 12)]


def family_word(a, b):
    return ((0,) + (1, 2) * a + (0, 2) * b) * 2


def theta_d(a, b):
    return Fr(90 * (a + b), a * (b + 1))


def gamma_d(a, b):
    return 180 - theta_d(a, b)


# =====================================================================
# Q(i): complex numbers with Fraction components, as (re, im) tuples
# =====================================================================

QZERO = (Fr(0), Fr(0))
QONE = (Fr(1), Fr(0))
QI = (Fr(0), Fr(1))


def qadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def qsub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def qmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def qneg(u):
    return (-u[0], -u[1])


def qconj(u):
    return (u[0], -u[1])


def qdiv_real(u, r):
    return (u[0] / r, u[1] / r)


# =====================================================================
# Backend 1: my Laurent ring Q(i)[X^, Y^], X = e^{i alpha/2}, Y = e^{i beta/2}
# elements: dict {(m, n): (re, im)}
# =====================================================================

class Ring:
    zero = {}
    one = {(0, 0): QONE}

    @staticmethod
    def E(m, n, c=QONE):
        return {(m, n): c} if c != QZERO else {}

    @staticmethod
    def add(p, q):
        r = dict(p)
        for k, c in q.items():
            s = qadd(r.get(k, QZERO), c)
            if s == QZERO:
                r.pop(k, None)
            else:
                r[k] = s
        return r

    @staticmethod
    def sub(p, q):
        return Ring.add(p, {k: qneg(c) for k, c in q.items()})

    @staticmethod
    def mul(p, q):
        r = {}
        for (m1, n1), c1 in p.items():
            for (m2, n2), c2 in q.items():
                k = (m1 + m2, n1 + n2)
                s = qadd(r.get(k, QZERO), qmul(c1, c2))
                if s == QZERO:
                    r.pop(k, None)
                else:
                    r[k] = s
        return r

    @staticmethod
    def star(p):
        return {(-m, -n): qconj(c) for (m, n), c in p.items()}

    @staticmethod
    def scale(p, c):
        return {k: qmul(v, c) for k, v in p.items() if qmul(v, c) != QZERO}

    @staticmethod
    def is_zero(p):
        return all(c == QZERO for c in p.values())


def r_sin(m, n):
    """sin((m alpha + n beta)/2) as a ring element."""
    return Ring.scale(Ring.sub(Ring.E(m, n), Ring.E(-m, -n)),
                      (Fr(0), Fr(-1, 2)))          # 1/(2i)


def r_cos(m, n):
    return Ring.scale(Ring.add(Ring.E(m, n), Ring.E(-m, -n)),
                      (Fr(1, 2), Fr(0)))


def r_im(p):
    """(p - star p)/(2i): the 'imaginary part' as a ring element."""
    return Ring.scale(Ring.sub(p, Ring.star(p)), (Fr(0), Fr(-1, 2)))


def ring_eval(p, alpha, beta):
    """Float value at angles (radians)."""
    z = 0j
    for (m, n), c in p.items():
        z += complex(c[0], c[1]) * cmath.exp(1j * (m * alpha + n * beta) / 2)
    return z


# =====================================================================
# Backend 2: pair numbers (value, star-value) at a fixed rational (u0, v0)
# off the unit torus.  star is a ring automorphism with star(X) = 1/X,
# star(Y) = 1/Y, star(i) = -i, so tracking the pair evaluates any expression
# built from +,-,*,star exactly, WITHOUT any polynomial representation.
# =====================================================================

class Pair:
    __slots__ = ("v", "s")

    def __init__(self, v, s):
        self.v, self.s = v, s

    def __add__(self, o):
        return Pair(qadd(self.v, o.v), qadd(self.s, o.s))

    def __sub__(self, o):
        return Pair(qsub(self.v, o.v), qsub(self.s, o.s))

    def __mul__(self, o):
        return Pair(qmul(self.v, o.v), qmul(self.s, o.s))

    def star(self):
        return Pair(self.s, self.v)

    def is_zero(self):
        return self.v == QZERO and self.s == QZERO


class PairOps:
    """Same interface as Ring, bound to one evaluation point (u0, v0)."""

    def __init__(self, u0, v0):
        self.u0, self.v0 = Fr(u0), Fr(v0)
        self.zero = Pair(QZERO, QZERO)
        self.one = Pair(QONE, QONE)

    def E(self, m, n, c=QONE):
        val = (self.u0 ** m * self.v0 ** n, Fr(0))
        sval = (self.u0 ** -m * self.v0 ** -n, Fr(0))
        return Pair(qmul(c, val), qmul(qconj(c), sval))

    @staticmethod
    def add(p, q):
        return p + q

    @staticmethod
    def sub(p, q):
        return p - q

    @staticmethod
    def mul(p, q):
        return p * q

    @staticmethod
    def star(p):
        return p.star()

    @staticmethod
    def scale(p, c):
        return p * Pair(c, qconj(c))

    @staticmethod
    def is_zero(p):
        return p.is_zero()


def ops_sin(ops, m, n):
    return ops.scale(ops.sub(ops.E(m, n), ops.E(-m, -n)), (Fr(0), Fr(-1, 2)))


def ops_cos(ops, m, n):
    return ops.scale(ops.add(ops.E(m, n), ops.E(-m, -n)), (Fr(1, 2), Fr(0)))


def ops_im(ops, p):
    return ops.scale(ops.sub(p, ops.star(p)), (Fr(0), Fr(-1, 2)))


# =====================================================================
# Symbolic unfolding by composed maps, generic over the backend.
#
# Circumdiameter-1 normalization (003's): A = 0, B = sin(alpha+beta),
# C = sin(beta) e^{i alpha}.  Mirror data (line through p_s, direction
# squared E_s):  side 0 = BC: p = B, E = e^{2i(pi-beta)} = Y^-4;
# side 1 = CA: p = A = 0, E = e^{2i alpha} = X^4; side 2 = AB: p = 0, E = 1.
# Reflection R_s(z) = E_s star(z) + c_s, c_s = p_s - E_s star(p_s).
# Composed map M_k = M_{k-1} o R_{s_k}; gate k = M_{k-1}(base side s_k).
# =====================================================================

def base_geometry(ops):
    if isinstance(ops, PairOps):
        A = ops.zero
        B = ops_sin(ops, 2, 2)
        C = ops.mul(ops_sin(ops, 0, 2), ops.E(2, 0))
        Es = [ops.E(0, -4), ops.E(4, 0), ops.E(0, 0)]
    else:
        A = Ring.zero
        B = r_sin(2, 2)
        C = Ring.mul(r_sin(0, 2), Ring.E(2, 0))
        Es = [Ring.E(0, -4), Ring.E(4, 0), Ring.E(0, 0)]
    Ps = [B, A, A]
    cs = [ops.sub(Ps[s], ops.mul(Es[s], ops.star(Ps[s]))) for s in range(3)]
    return [A, B, C], Es, cs


def compose2(ops, word):
    """Correct composition.  M(z) = L z + t (sigma=+1) or L star(z) + t
    (sigma=-1).  M o R_s with R_s(z) = E star(z) + c:
      sigma=+1:  z -> L(E star z + c) + t   = (L E) star z + (L c + t)
      sigma=-1:  z -> L star(E star z + c) + t
                                            = (L star E) z + (L star(c) + t)
    """
    verts, Es, cs = base_geometry(ops)
    sigma, L, t = 1, ops.one, ops.zero
    gates = []
    for s in word:
        i, j = SIDE_ENDS[s]
        gates.append((apply_map(ops, sigma, L, t, verts[i]),
                      apply_map(ops, sigma, L, t, verts[j])))
        if sigma == 1:
            L, t = ops.mul(L, Es[s]), ops.add(ops.mul(L, cs[s]), t)
            sigma = -1
        else:
            newL = ops.mul(L, ops.star(Es[s]))
            t = ops.add(ops.mul(L, ops.star(cs[s])), t)
            L = newL
            sigma = 1
    return (sigma, L, t), gates


def apply_map(ops, sigma, L, t, z):
    zz = z if sigma == 1 else ops.star(z)
    return ops.add(ops.mul(L, zz), t)


# ---------------------------------------------------------------------
# my own integer bookkeeping for the composed linear part:
# angles Lambda = p*alpha + q*beta + r*pi with linear part e^{i Lambda} z
# (rot) or e^{i Lambda} conj z (ref).  R_s has ref-angle 2*theta_s:
# side 0: 2(pi - beta) = (0,-2,2); side 1: 2 alpha = (2,0,0); side 2: 0.
# Compose M o r: rot(Lam) o ref(E) = ref(Lam + E);
#                ref(Lam) o ref(E) = rot(Lam - E).
# ---------------------------------------------------------------------

MIRROR_ANGLE = {0: (0, -2, 2), 1: (2, 0, 0), 2: (0, 0, 0)}


def linear_angle(word):
    kind, lam = "rot", (0, 0, 0)
    for s in word:
        e = MIRROR_ANGLE[s]
        if kind == "rot":
            lam = tuple(lam[i] + e[i] for i in range(3))
            kind = "ref"
        else:
            lam = tuple(lam[i] - e[i] for i in range(3))
            kind = "rot"
    return kind, lam


# =====================================================================
# identity targets (003's I1-I3, D1, D2), built in the given backend
# =====================================================================

def identity_residuals(ops, a, b):
    """Returns dict name -> residual (backend element) that must be zero,
    plus the glide data used.  Follows MY derivation:
      m = Im(star(delta) w)/2, p(z) = Im(star(delta) z),
      glide action p(u(z)) = 2m - p(z)."""
    word = family_word(a, b)
    half = word[:len(word) // 2]
    (sigma, L, w), gates = compose2(ops, half)
    assert sigma == -1, "half word must be orientation-reversing"

    kind, (lp, lq, lr) = linear_angle(half)
    assert kind == "ref"
    assert lp % 2 == 0 and lq % 2 == 0, "odd half-exponent in glide angle"
    # mu = e^{i Lambda} = (-1)^r X^{2p} Y^{2q}; delta = sqrt with delta^2 = mu
    if lr % 2 == 0:
        delta = ops.E(lp, lq)
    else:
        delta = ops.E(lp, lq, QI)
    # check L really is that monomial (residual must vanish)
    mu_res = ops.sub(L, ops.mul(delta, delta))

    dstar = ops.star(delta)
    proj = lambda z: ops_im(ops, ops.mul(dstar, z))
    m2 = proj(w)                       # 2*m
    m = ops.scale(m2, (Fr(1, 2), Fr(0)))
    tau = ops.add(w, ops.mul(L, ops.star(w)))

    verts, Es, cs = base_geometry(ops)
    A, B, C = verts
    A1 = cs[0]                          # R_0(A) = E_0 star(0) + c_0 = c_0
    # gate check: gate1 = (B, C); gate2 = (M_1(C), M_1(A)) = (C, A1)
    g1, g2 = gates[0], gates[1]
    gate1_res = (ops.sub(g1[0], B), ops.sub(g1[1], C))
    gate2_res = (ops.sub(g2[0], C), ops.sub(g2[1], A1))

    pB, pC, pA1 = proj(B), proj(C), proj(A1)

    def S(mm, nn):
        return ops_sin(ops, mm, nn)

    def Co(mm, nn):
        return ops_cos(ops, mm, nn)

    I1_rhs = ops.sub(ops.zero, ops.add(
        ops.mul(ops.mul(Co(2 * a, 0), S(2, 0)), S(0, 2 * b)),
        ops.mul(ops.mul(S(2 * a, 0), S(0, 2)), Co(2, 2 * (b + 1)))))
    I2_rhs = ops.mul(ops.mul(Co(2 * a, 0), S(0, 2 * (b + 1))), S(2, 2))
    I3_rhs = ops.mul(ops.mul(Co(0, 2 * (b + 1)), S(2 * a, 0)), S(2, 2))
    D1_rhs = ops.mul(S(2, 0), S(2 * a, -2 * b))
    D2_rhs = ops.sub(ops.zero,
                     ops.mul(S(0, 2), S(2 * (a - 1), -2 * (b + 1))))

    res = {
        "mu_monomial": mu_res,
        "tau_parallel": ops_im(ops, ops.mul(dstar, tau)),
        "gate1_u": gate1_res[0], "gate1_v": gate1_res[1],
        "gate2_u": gate2_res[0], "gate2_v": gate2_res[1],
        "I1": ops.sub(ops.sub(m, pC), I1_rhs),
        "I2": ops.sub(ops.sub(pA1, m), I2_rhs),
        "I3": ops.sub(ops.sub(pB, m), I3_rhs),
        "D1": ops.sub(ops.sub(pB, pC), D1_rhs),
        "D2": ops.sub(ops.sub(pA1, pC), D2_rhs),
    }
    # glide action p(u(P)) + p(P) - 2m = 0 for all base vertices and A1
    for name, P in (("gA", A), ("gB", B), ("gC", C), ("gA1", A1)):
        uP = apply_map(ops, -1, L, w, P)
        res["glide_action_" + name] = ops.sub(ops.add(proj(uP), proj(P)), m2)
    return res


# =====================================================================
# float + exact corridors in APEX coordinates (A=(0,0), B=(1,0), C=(x,y)),
# the harness normalization -- independent of the circumdiameter one.
# Composed-map unfolding (different from 003's and 002's vertex chains).
# =====================================================================

def corridor_apex_float(word, x, y):
    """(width, lo, hi, tau) via composed maps in complex floats."""
    V = (0j, 1 + 0j, complex(x, y))
    Ep, Pp = [], []
    for s in range(3):
        i, j = SIDE_ENDS[s]
        d = V[j] - V[i]
        Ep.append(d * d / (d.real * d.real + d.imag * d.imag))
        Pp.append(V[i])
    cvals = [Pp[s] - Ep[s] * Pp[s].conjugate() for s in range(3)]
    sigma, L, t = 1, 1 + 0j, 0j
    gates = []
    for s in word:
        i, j = SIDE_ENDS[s]
        if sigma == 1:
            gates.append((L * V[i] + t, L * V[j] + t))
            L, t = L * Ep[s], L * cvals[s] + t
            sigma = -1
        else:
            gates.append((L * V[i].conjugate() + t, L * V[j].conjugate() + t))
            newL = L * Ep[s].conjugate()
            t = L * cvals[s].conjugate() + t
            L = newL
            sigma = 1
    tau = t if sigma == 1 else None
    if sigma != 1:
        return None
    # translation check: L must be 1
    if abs(L - 1) > 1e-9:
        return None
    lo, hi = -math.inf, math.inf
    for (zu, zv) in gates:
        pu = (tau.conjugate() * zu).imag
        pv = (tau.conjugate() * zv).imag
        lo = max(lo, min(pu, pv))
        hi = min(hi, max(pu, pv))
    return (hi - lo, lo, hi, tau, gates, (sigma, L, t))


def width_apex_float(word, x, y):
    c = corridor_apex_float(word, x, y)
    return -math.inf if c is None else c[0]


def corridor_apex_exact(word, x, y):
    """Exact rational corridor (lo, hi, tau) at rational apex."""
    x, y = Fr(x), Fr(y)
    V = ((Fr(0), Fr(0)), (Fr(1), Fr(0)), (x, y))
    Ep, cvals = [], []
    for s in range(3):
        i, j = SIDE_ENDS[s]
        d = qsub(V[j], V[i])
        n2 = d[0] * d[0] + d[1] * d[1]
        E = qdiv_real(qmul(d, d), n2)
        Ep.append(E)
        cvals.append(qsub(V[i], qmul(E, qconj(V[i]))))
    sigma, L, t = 1, QONE, QZERO
    gates = []
    for s in word:
        i, j = SIDE_ENDS[s]
        if sigma == 1:
            gates.append((qadd(qmul(L, V[i]), t), qadd(qmul(L, V[j]), t)))
            L, t = qmul(L, Ep[s]), qadd(qmul(L, cvals[s]), t)
            sigma = -1
        else:
            gates.append((qadd(qmul(L, qconj(V[i])), t),
                          qadd(qmul(L, qconj(V[j])), t)))
            newL = qmul(L, qconj(Ep[s]))
            t = qadd(qmul(L, qconj(cvals[s])), t)
            L = newL
            sigma = 1
    if sigma != 1 or L != QONE:
        return None
    tau = t
    if tau == QZERO:
        return None
    lo = hi = None
    for (zu, zv) in gates:
        pu = qmul(qconj(tau), zu)[1]
        pv = qmul(qconj(tau), zv)[1]
        if pu > pv:
            pu, pv = pv, pu
        lo = pu if lo is None else max(lo, pu)
        hi = pv if hi is None else min(hi, pv)
    return lo, hi, tau


# =====================================================================
# my own interval stack (NOT unfold.py's)
# =====================================================================

class MyIv:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = Fr(lo)
        self.hi = Fr(lo if hi is None else hi)
        assert self.lo <= self.hi

    def __add__(self, o):
        o = miv(o)
        return MyIv(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = miv(o)
        return MyIv(self.lo - o.hi, self.hi - o.lo)

    def __neg__(self):
        return MyIv(-self.hi, -self.lo)

    def __mul__(self, o):
        o = miv(o)
        c = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return MyIv(min(c), max(c))

    def __truediv__(self, o):
        o = miv(o)
        assert not (o.lo <= 0 <= o.hi)
        c = (self.lo / o.lo, self.lo / o.hi, self.hi / o.lo, self.hi / o.hi)
        return MyIv(min(c), max(c))

    __radd__ = __add__
    __rmul__ = __mul__

    def __rsub__(self, o):
        return miv(o) - self


def miv(x):
    return x if isinstance(x, MyIv) else MyIv(x)


def _arctan_series(x, terms):
    """Enclosure of arctan(x) for rational 0 < x < 1: alternating series,
    remainder bounded by the first omitted term."""
    x = Fr(x)
    total = Fr(0)
    xx = x * x
    term = x
    sign = 1
    for k in range(terms):
        total += sign * term / (2 * k + 1)
        term *= xx
        sign = -sign
    r = term / (2 * terms + 1)
    return MyIv(total - r, total + r)


def my_pi():
    """pi = 16 arctan(1/5) - 4 arctan(1/239) (Machin)."""
    return _arctan_series(Fr(1, 5), 60) * 16 - _arctan_series(Fr(1, 239), 30) * 4


PI = my_pi()          # width < 1e-70


def _sin_pt(x, terms=30):
    x = Fr(x)
    total, term = Fr(0), x
    for k in range(terms):
        total += term
        term = -term * x * x / ((2 * k + 2) * (2 * k + 3))
    r = abs(term)
    return MyIv(total - r, total + r)


def _cos_pt(x, terms=30):
    x = Fr(x)
    total, term = Fr(0), Fr(1)
    for k in range(terms):
        total += term
        term = -term * x * x / ((2 * k + 1) * (2 * k + 2))
    r = abs(term)
    return MyIv(total - r, total + r)


def sin_rad(a):
    a = miv(a)
    h = (a.hi - a.lo) / 2
    v = _sin_pt((a.lo + a.hi) / 2)
    return MyIv(v.lo - h, v.hi + h)      # 1-Lipschitz


def cos_rad(a):
    a = miv(a)
    h = (a.hi - a.lo) / 2
    v = _cos_pt((a.lo + a.hi) / 2)
    return MyIv(v.lo - h, v.hi + h)


def sin_deg(x):
    return sin_rad(miv(x) * PI / 180)


def cos_deg(x):
    return cos_rad(miv(x) * PI / 180)


# =====================================================================
# subcommands
# =====================================================================

def cmd_selftest(args):
    ok = True
    rng = random.Random(20260731)

    # base geometry: side lengths in circumdiameter normalization
    ops = Ring
    verts, Es, cs = base_geometry(ops)
    al, be = 0.53, 0.29
    A0, B0, C0 = (ring_eval(v, al, be) for v in verts)
    good = (abs(abs(B0 - A0) - math.sin(al + be)) < 1e-12
            and abs(abs(C0 - A0) - math.sin(be)) < 1e-12
            and abs(abs(C0 - B0) - math.sin(al)) < 1e-12)
    ok &= good
    print("[circumdiameter base geometry]", "OK" if good else "FAIL")

    # ring reflections fix their own mirror endpoints: R_0(B)=B, R_0(C)=C
    R0 = lambda z: Ring.add(Ring.mul(Es[0], Ring.star(z)), cs[0])
    good = Ring.is_zero(Ring.sub(R0(verts[1]), verts[1])) \
        and Ring.is_zero(Ring.sub(R0(verts[2]), verts[2]))
    ok &= good
    print("[R_0 fixes B and C exactly in the ring]", "OK" if good else "FAIL")

    # my float corridor vs 002's independent skeptic_family corridor
    from skeptic_family import width as sk_width
    agree = True
    for (fa, fb) in [(2, 1), (3, 3), (5, 3)]:
        w = family_word(fa, fb)
        for _ in range(40):
            g = rng.uniform(91.0, 175.0)
            th = math.radians(180.0 - g)
            a1 = rng.uniform(0.05 * th, 0.95 * th)
            x, y = apex_from_angles(math.degrees(a1),
                                    math.degrees(th - a1))
            w1 = width_apex_float(w, x, y)
            w2 = sk_width(w, x, y)
            if abs(w1 - w2) > 1e-9 * max(1.0, abs(w1)):
                agree = False
    ok &= agree
    print("[my corridor == 002 skeptic corridor, 120 samples]",
          "OK" if agree else "FAIL")

    # exact corridor sign == harness certify at a point box (2 samples)
    sys.path.insert(0, os.path.join(HERE, "..", "..", "..",
                                    "harness", "billiards-triangles"))
    import unfold
    w = family_word(2, 1)
    for (x, y, expect_alive) in [(Fr(4899, 16384), Fr(9699, 32768), True),
                                 (Fr(1, 3), Fr(1, 5), False)]:
        c = corridor_apex_exact(w, x, y)
        alive_mine = c is not None and c[1] > c[0]
        verdict = unfold.certify(w, (unfold.Iv(x), unfold.Iv(y)))[0]
        good = (alive_mine == (verdict == "TRUE"))
        if expect_alive is not None:
            good &= (alive_mine == expect_alive)
        ok &= good
        print(f"[exact corridor vs harness certify at point ({x},{y})]",
              "OK" if good else f"FAIL (mine {alive_mine}, harness {verdict})")

    # my pi enclosure sanity: rational bracket + tightness (a float compare
    # would be blind here: the double nearest pi sits below PI.lo)
    good = (Fr(314159265, 100000000) < PI.lo <= PI.hi
            < Fr(31415926536, 10000000000)
            and PI.hi - PI.lo < Fr(1, 10 ** 60))
    ok &= good
    print("[my Machin pi enclosure]", "OK" if good else "FAIL")

    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def apex_from_angles(alpha_deg, beta_deg):
    ta = math.tan(math.radians(alpha_deg))
    tb = math.tan(math.radians(beta_deg))
    return tb / (ta + tb), ta * tb / (ta + tb)


def cmd_identities(args):
    results = []
    allok = True
    rng = random.Random(20260731)
    for (a, b) in MEMBERS:
        res = identity_residuals(Ring, a, b)
        member_ok = {k: Ring.is_zero(v) for k, v in res.items()}
        # backend 2 spot checks: 20 seeded random rational off-torus points
        spots_ok = True
        for _ in range(args.spots):
            u0 = Fr(rng.randint(2, 40), rng.randint(1, 7))
            v0 = Fr(rng.randint(2, 40), rng.randint(1, 7))
            if u0 == v0:
                v0 += 1
            pres = identity_residuals(PairOps(u0, v0), a, b)
            if not all(p.is_zero() for p in pres.values()):
                spots_ok = False
                break
        entry = {"a": a, "b": b, "ring": member_ok,
                 "pair_spot_points": args.spots, "pair_spots_ok": spots_ok,
                 "ok": all(member_ok.values()) and spots_ok}
        allok &= entry["ok"]
        results.append(entry)
        print(f"W({a},{b}): ring={'OK' if all(member_ok.values()) else member_ok} "
              f"spots={'OK' if spots_ok else 'FAIL'}")
    out = {"command": "identities", "members": results, "all_ok": allok}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print("IDENTITIES", "ALL CONFIRMED" if allok else "FAILURE")
    sys.exit(0 if allok else 1)


def cmd_gridproof(args):
    """Proof-grade certification for one member by off-torus grid evaluation.

    Degree bound (per variable, absolute Laurent exponent), derived in the
    attempt record: every residual checked has exponents <= D = 6*nh + 8,
    nh = half-word length.  A Laurent polynomial with exponents in [-D, D]
    in each variable vanishing on a (2D+1) x (2D+1) grid of distinct
    rational abscissae is identically zero (multiply by u^D v^D, then
    univariate interpolation in each variable in turn)."""
    a, b = args.a, args.b
    nh = 2 * a + 2 * b + 1
    D = 6 * nh + 8
    N = 2 * D + 1
    names = None
    for iu in range(N):
        u0 = Fr(iu + 2)
        for iv_ in range(N):
            v0 = Fr(iv_ + 2, 2 * D + 4)     # distinct values, != u-values
            res = identity_residuals(PairOps(u0, v0), a, b)
            if names is None:
                names = sorted(res)
            for k, p in res.items():
                if not p.is_zero():
                    print(f"NONZERO residual {k} at u={u0}, v={v0}")
                    sys.exit(1)
        if iu % 25 == 0:
            print(f"  grid row {iu + 1}/{N} done")
    out = {"command": "gridproof", "a": a, "b": b, "half_len": nh,
           "degree_bound": D, "grid": [N, N],
           "residuals_all_zero": names,
           "conclusion": (f"I1,I2,I3,D1,D2 + glide facts for W({a},{b}) are "
                          f"exact ring identities (grid interpolation proof)")}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print(json.dumps(out, indent=1))


def glide_data_float(word, x, y):
    """Half-word glide data in apex coords: (m, delta, half_gate_projs,
    full_corridor, half map)."""
    full = corridor_apex_float(word, x, y)
    if full is None:
        return None
    wfull, LO, HI, tau, gates, _ = full
    nh = len(word) // 2
    # recompute composed map for half word
    half = word[:nh]
    V = (0j, 1 + 0j, complex(x, y))
    Ep, Pp = [], []
    for s in range(3):
        i, j = SIDE_ENDS[s]
        d = V[j] - V[i]
        Ep.append(d * d / (d.real ** 2 + d.imag ** 2))
        Pp.append(V[i])
    cvals = [Pp[s] - Ep[s] * Pp[s].conjugate() for s in range(3)]
    sigma, L, t = 1, 1 + 0j, 0j
    for s in half:
        if sigma == 1:
            L, t = L * Ep[s], L * cvals[s] + t
            sigma = -1
        else:
            newL = L * Ep[s].conjugate()
            t = L * cvals[s].conjugate() + t
            L = newL
            sigma = 1
    assert sigma == -1
    delta = cmath.sqrt(L)
    m = (delta.conjugate() * t).imag / 2
    projs = []
    for (zu, zv) in gates[:nh]:
        pu = (delta.conjugate() * zu).imag
        pv = (delta.conjugate() * zv).imag
        projs.append((min(pu, pv), max(pu, pv)))
    return {"full_width": wfull, "full": (LO, HI), "tau": tau,
            "delta": delta, "m": m, "half_projs": projs, "halfmap": (L, t)}


def cmd_reduction(args):
    rng = random.Random(20260731)
    allok = True
    rows = []
    for (a, b) in [(2, 1), (3, 3), (4, 2), (5, 3), (5, 4), (10, 9)]:
        w = family_word(a, b)
        n_alive = n_checks = 0
        worst_gate_margin_fail = None
        eq_fail = 0
        for _ in range(args.samples):
            g = rng.uniform(91.0, float(gamma_d(a, b)) + 2.0)
            th = math.radians(180.0 - g)
            a1 = rng.uniform(0.03 * th, 0.97 * th)
            x, y = apex_from_angles(math.degrees(a1), math.degrees(th - a1))
            gd = glide_data_float(w, x, y)
            if gd is None:
                continue
            n_checks += 1
            # (1) full corridor == intersection over half gates of
            #     I_k cap (2m - I_k), up to |tau|/|delta| scaling: compare
            #     via the delta-projection directly.
            tau, delta, m = gd["tau"], gd["delta"], gd["m"]
            # tau parallel to delta?
            par = abs((delta.conjugate() * tau).imag) / abs(tau)
            if par > 1e-9:
                eq_fail += 1
                continue
            lo = max(p[0] for p in gd["half_projs"])
            hi = min(p[1] for p in gd["half_projs"])
            glo, ghi = max(lo, 2 * m - hi), min(hi, 2 * m - lo)
            w_red = ghi - glo                      # reduced corridor width
            scale = abs(tau)                       # full projs use tau-normal
            w_full = gd["full_width"] / scale      # rescaled to unit normal
            if abs(w_full - w_red) > 1e-8 * max(1.0, abs(w_red)):
                eq_fail += 1
            # (2) positive width => m strictly inside every half gate
            if gd["full_width"] > 1e-12:
                n_alive += 1
                margin = min(min(m - p[0], p[1] - m)
                             for p in gd["half_projs"])
                if margin <= 0:
                    worst_gate_margin_fail = (g, x, y, margin)
        # targeted near-death alive apexes (windows are too thin for the
        # random sampler on the larger members)
        for dg in (1e-3, 1e-5, 1e-7):
            g = float(gamma_d(a, b)) - dg
            bw, arg, _ = best_alive_arc(w, a, b, g, 800, True)
            if arg is None or bw <= 1e-13:
                continue
            x, y = apex_from_angles(*arg)
            gd = glide_data_float(w, x, y)
            if gd is None:
                continue
            n_checks += 1
            tau, m = gd["tau"], gd["m"]
            lo = max(p[0] for p in gd["half_projs"])
            hi = min(p[1] for p in gd["half_projs"])
            glo, ghi = max(lo, 2 * m - hi), min(hi, 2 * m - lo)
            w_full = gd["full_width"] / abs(tau)
            if abs(w_full - (ghi - glo)) > 1e-8 * max(1.0, abs(ghi - glo)):
                eq_fail += 1
            if gd["full_width"] > 1e-13:
                n_alive += 1
                margin = min(min(m - p[0], p[1] - m)
                             for p in gd["half_projs"])
                if margin <= 0:
                    worst_gate_margin_fail = (g, x, y, margin)
        okrow = (eq_fail == 0 and worst_gate_margin_fail is None)
        allok &= okrow
        rows.append({"a": a, "b": b, "checks": n_checks, "alive": n_alive,
                     "reduction_mismatches": eq_fail,
                     "m_outside_gate_at": worst_gate_margin_fail,
                     "ok": okrow})
        print(f"W({a},{b}): {n_checks} triangles, {n_alive} alive, "
              f"mismatches={eq_fail}, m-violation={worst_gate_margin_fail}")
    out = {"command": "reduction", "samples_per_member": args.samples,
           "seed": 20260731, "rows": rows, "all_ok": allok}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print("REDUCTION", "CONFIRMED" if allok else "FAILURE")
    sys.exit(0 if allok else 1)


def signs_from_identities(a, b, alpha, beta):
    """(N1, N2, N3) floats from the closed forms (degrees in, radians used)."""
    al, be = math.radians(alpha), math.radians(beta)
    G = (math.cos(a * al) * math.sin(al) * math.sin(b * be)
         + math.sin(a * al) * math.sin(be) * math.cos(al + (b + 1) * be))
    N1 = -G
    N2 = math.cos(a * al) * math.sin((b + 1) * be) * math.sin(al + be)
    N3 = math.cos((b + 1) * be) * math.sin(a * al) * math.sin(al + be)
    return N1, N2, N3


def cmd_casetree(args):
    """Search for (alpha, beta) with theta <= theta_d yet Case I or II
    satisfied -- a counterexample to 003's case analysis."""
    rng = random.Random(20260731)
    hits = []
    total = 0
    for (a, b) in MEMBERS:
        td = float(theta_d(a, b))
        cand_alphas = []
        # targeted: branch boundaries of a*alpha and (b+1)*beta
        for base in (90.0, 180.0, 270.0, 360.0):
            for eps in (1e-6, 1e-4, 1e-2, 1.0, 5.0):
                cand_alphas += [(base - eps) / a, (base + eps) / a]
        for _ in range(args.samples):
            mode = rng.random()
            if mode < 0.5:
                theta = td * rng.random()
            elif mode < 0.9:
                theta = td * (1 - 10 ** -rng.uniform(0, 12))
            else:
                theta = td          # boundary itself
            if rng.random() < 0.3 and cand_alphas:
                alpha = rng.choice(cand_alphas)
            elif rng.random() < 0.5:
                # near the death corner alpha = 90/a
                alpha = 90.0 / a - 10 ** -rng.uniform(0, 12)
            else:
                alpha = theta * rng.random()
            beta = theta - alpha
            if alpha <= 0 or beta <= 0:
                continue
            total += 1
            N1, N2, N3 = signs_from_identities(a, b, alpha, beta)
            tol = args.margin
            if (N1 > tol and N2 > tol and N3 > tol) or \
               (N1 < -tol and N2 < -tol and N3 < -tol):
                hits.append({"a": a, "b": b, "alpha": alpha, "beta": beta,
                             "theta": alpha + beta, "theta_d": td,
                             "N": [N1, N2, N3]})
        # positive control: just above theta_d at the corner, Case I fires
        ctl = False
        for k in range(2, 10):
            theta = td * (1 + 10 ** -k)
            alpha = 90.0 / a - td * 10 ** -k * (b / (b + 1.0)) * 0.5
            beta = theta - alpha
            N1, N2, N3 = signs_from_identities(a, b, alpha, beta)
            if N1 > 0 and N2 > 0 and N3 > 0:
                ctl = True
                break
        print(f"W({a},{b}): searched, hits so far {len(hits)}, "
              f"positive-control fires above theta_d: {ctl}")
    out = {"command": "casetree", "samples_per_member": args.samples,
           "seed": 20260731, "float_margin": args.margin,
           "total_points": total, "counterexamples": hits}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print(f"CASETREE: {len(hits)} counterexample candidates in {total} points")
    sys.exit(1 if hits else 0)


def certify_gamma_bounds(x, y, g_lo, g_hi):
    """My own certification: gamma(x,y) > g_lo and < g_hi (degrees,
    rational), for rational apex with obtuse angle at C."""
    x, y = Fr(x), Fr(y)
    dot = (-x) * (1 - x) + y * y
    if dot >= 0:
        return False, False
    r = dot * dot / ((x * x + y * y) * ((1 - x) ** 2 + y * y))
    c_lo = cos_deg(miv(Fr(g_lo)))
    c_hi = cos_deg(miv(Fr(g_hi)))
    # gamma in (90,180); cos^2 increasing there; gamma > g iff r > cos^2(g)
    sq_lo = MyIv(min(c_lo.lo * c_lo.lo, c_lo.hi * c_lo.hi)
                 if not (c_lo.lo <= 0 <= c_lo.hi) else Fr(0),
                 max(c_lo.lo * c_lo.lo, c_lo.hi * c_lo.hi))
    sq_hi = MyIv(min(c_hi.lo * c_hi.lo, c_hi.hi * c_hi.hi)
                 if not (c_hi.lo <= 0 <= c_hi.hi) else Fr(0),
                 max(c_hi.lo * c_hi.lo, c_hi.hi * c_hi.hi))
    return r > sq_lo.hi, r < sq_hi.lo


def cmd_alivecheck(args):
    data_dir = os.path.join(HERE, "..", "data")
    files = sorted(f for f in os.listdir(data_dir)
                   if f.startswith("deathlaw_exact_alive_W"))
    rows = []
    allok = True
    for f in files:
        d = json.load(open(os.path.join(data_dir, f)))
        a, b = d["a"], d["b"]
        x = Fr(d["apex"][0])
        y = Fr(d["apex"][1])
        w = family_word(a, b)
        c = corridor_apex_exact(w, x, y)
        alive = c is not None and c[1] > c[0]
        width_match = None
        if alive:
            width_match = (str(c[1] - c[0]) == d["exact_corridor_width"])
        claimed_lo = Fr(d["certified_gamma_lower_bound_deg"])
        gd = gamma_d(a, b)
        lo_ok, hi_ok = certify_gamma_bounds(x, y, claimed_lo, gd)
        okrow = alive and lo_ok and hi_ok
        allok &= okrow
        rows.append({"file": f, "a": a, "b": b,
                     "exact_corridor_positive": alive,
                     "width_matches_003_exactly": width_match,
                     "my_certified_gamma_gt": str(claimed_lo),
                     "my_certified_gamma_lt_gamma_d": hi_ok,
                     "ok": okrow})
        print(f"{f}: alive={alive} width_match={width_match} "
              f"gamma>({claimed_lo})={lo_ok} gamma<gamma_d={hi_ok}")
    out = {"command": "alivecheck", "rows": rows, "all_ok": allok}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print("ALIVECHECK", "CONFIRMED" if allok else "FAILURE")
    sys.exit(0 if allok else 1)


def best_alive_arc(word, a, b, g, uniform, full):
    """(best width, best (alpha, beta), best_by_half) on arc gamma = g.
    best_by_half = (best width alpha>=beta, best width alpha<beta)."""
    theta = 180.0 - g
    edge_a = 90.0 / a
    edge_b = 90.0 / (b + 1)
    lo_cut = 0.0 if full else theta / 2
    alphas = [lo_cut + (theta - lo_cut) * (k + 0.5) / uniform
              for k in range(uniform)]
    for k in range(60):
        d = 0.5 * 2.0 ** -k
        for e in (edge_a - d, edge_a + d,
                  theta - (edge_b - d), theta - (edge_b + d)):
            if lo_cut < e < theta:
                alphas.append(e)
    best, arg = -math.inf, None
    best_ge = best_lt = -math.inf
    for al in alphas:
        be = theta - al
        if al <= 0 or be <= 0:
            continue
        x, y = apex_from_angles(al, be)
        v = width_apex_float(word, x, y)
        if v > best:
            best, arg = v, (al, be)
        if al >= be:
            best_ge = max(best_ge, v)
        else:
            best_lt = max(best_lt, v)
    return best, arg, (best_ge, best_lt)


ALIVE_TOL = 1e-12


def cmd_death(args):
    a, b = args.a, args.b
    w = family_word(a, b)
    pred = float(gamma_d(a, b))
    wl = best_alive_arc(w, a, b, pred - 0.5, args.uniform, args.full)
    wh = best_alive_arc(w, a, b, pred + 0.5, args.uniform, args.full)
    if not (wl[0] > ALIVE_TOL and wh[0] <= ALIVE_TOL):
        print(f"# bracket failure: w(lo)={wl[0]:.3g} w(hi)={wh[0]:.3g}")
        sys.exit(1)
    g_lo, g_hi = pred - 0.5, pred + 0.5
    last = None
    for _ in range(args.iters):
        mid = 0.5 * (g_lo + g_hi)
        mw, ar, halves = best_alive_arc(w, a, b, mid, args.uniform, args.full)
        if mw > ALIVE_TOL:
            g_lo, last = mid, (mw, ar, halves)
        else:
            g_hi = mid
    res = {"command": "death", "a": a, "b": b, "len": len(w),
           "uniform": args.uniform, "iters": args.iters, "full": args.full,
           "death_between": [g_lo, g_hi],
           "predicted": pred,
           "measured_minus_predicted": 0.5 * (g_lo + g_hi) - pred}
    if last:
        mw, ar, halves = last
        res["last_alive"] = {"alpha": ar[0], "beta": ar[1], "width": mw,
                             "best_width_alpha_ge_beta": halves[0],
                             "best_width_alpha_lt_beta": halves[1]}
    if args.out:
        json.dump(res, open(args.out, "w"), indent=1)
    print(json.dumps(res))


def cmd_deadhunt(args):
    """Try to find an alive apex at gamma >= gamma_d.  Any float positive
    is exact-checked before being reported as a kill."""
    a, b = args.a, args.b
    w = family_word(a, b)
    gd = float(gamma_d(a, b))
    rows = []
    kills = []
    for dg in (0.0, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 0.5):
        g = gd + dg
        best, arg, halves = best_alive_arc(w, a, b, g, args.uniform, True)
        row = {"gamma_minus_gamma_d": dg, "best_float_width": best,
               "at": arg}
        if best > ALIVE_TOL:
            # exact confirmation attempt: snap and require certified gamma
            al, be = arg
            x, y = apex_from_angles(al, be)
            xq = Fr(round(x * 2 ** 45), 2 ** 45)
            yq = Fr(round(y * 2 ** 45), 2 ** 45)
            c = corridor_apex_exact(w, xq, yq)
            exact_alive = c is not None and c[1] > c[0]
            lo_ok, _ = certify_gamma_bounds(xq, yq, gamma_d(a, b), Fr(180))
            row["exact_alive"] = exact_alive
            row["gamma_certified_ge_gamma_d"] = lo_ok
            if exact_alive and lo_ok:
                kills.append(row)
        rows.append(row)
        print(row)
    out = {"command": "deadhunt", "a": a, "b": b, "uniform": args.uniform,
           "rows": rows, "kills": kills}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print("DEADHUNT:", "KILL FOUND" if kills else
          "no alive point at or above gamma_d")
    sys.exit(1 if kills else 0)


# ----------------------------------------------------------------- Lemma C

def h2_iv(v, a, b, v0):
    return (sin_deg(v) * sin_deg((miv(90) - v) * Fr(1, a)) * cos_deg(v0)
            - sin_deg(miv(v0)) * sin_deg(v * Fr(b, a)) * cos_deg(v))


def h2p_iv(v, a, b, v0):
    """dH2/dv up to the positive factor pi/180 (hand-derived; see record)."""
    c0, s0 = cos_deg(miv(v0)), sin_deg(miv(v0))
    return (cos_deg(v) * sin_deg((miv(90) - v) * Fr(1, a)) * c0
            - sin_deg(v) * cos_deg((miv(90) - v) * Fr(1, a)) * c0 * Fr(1, a)
            - s0 * cos_deg(v * Fr(b, a)) * cos_deg(v) * Fr(b, a)
            + s0 * sin_deg(v * Fr(b, a)) * sin_deg(v))


def cmd_lemmac(args):
    a, b = args.a, args.b
    v0 = Fr(90, b + 1)
    report = {"command": "lemmac", "a": a, "b": b, "v0": str(v0)}
    # H2(v0) = 0 needs (90 - v0)/a == b*v0/a: exact rational check
    assert (90 - v0) == b * v0, "endpoint identity fails"
    report["endpoint_identity_90_minus_v0_eq_b_v0"] = True

    # leg (iii): H2' < 0 on [v0 - d2, v0]
    d2 = v0 / 32
    while d2 > v0 / 2 ** 24:
        e = h2p_iv(MyIv(v0 - d2, v0), a, b, v0)
        if e.hi < 0:
            break
        d2 /= 2
    else:
        report["ok"] = False
        report["fail"] = "H2' sign leg"
        print(json.dumps(report, indent=1))
        sys.exit(1)
    report["d2"] = str(d2)
    report["h2p_hi_on_zone"] = float(h2p_iv(MyIv(v0 - d2, v0), a, b, v0).hi)

    # leg (i): for v in (0, d1]:  sin(b v /a deg) <= (b/a) tan(v deg)
    # (sin x <= x <= tan x in radians), so
    # H2 >= sin(v)[ sin((90-v)/a) cos(v0) - (b/a) sin(v0)/cos(v) ]
    #    >= sin(v)[ sin((90-d1)/a) cos(v0) - (b/a) sin(v0)/cos(d1) ]
    d1 = v0 / 8
    while d1 > v0 / 2 ** 24:
        br = (sin_deg((miv(90) - miv(d1)) * Fr(1, a)) * cos_deg(miv(v0))
              - sin_deg(miv(v0)) * miv(Fr(b, a)) / cos_deg(miv(d1)))
        if br.lo > 0:
            break
        d1 /= 2
    else:
        report["ok"] = False
        report["fail"] = "near-0 leg"
        print(json.dumps(report, indent=1))
        sys.exit(1)
    report["d1"] = str(d1)

    # leg (ii): adaptive bisection on [d1, v0 - d2]
    stack = [(d1, v0 - d2)]
    leaves = 0
    while stack:
        lo, hi = stack.pop()
        e = h2_iv(MyIv(lo, hi), a, b, v0)
        if e.lo > 0:
            leaves += 1
            continue
        if hi - lo < v0 / 2 ** 40:
            report["ok"] = False
            report["fail"] = f"bisection stuck at [{lo},{hi}]"
            print(json.dumps(report, indent=1))
            sys.exit(1)
        mid = (lo + hi) / 2
        stack += [(lo, mid), (mid, hi)]
    report["bisection_leaves"] = leaves
    report["ok"] = True
    if args.out:
        json.dump(report, open(args.out, "w"), indent=1)
    print(json.dumps(report, indent=1))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)

    c = sub.add_parser("identities")
    c.add_argument("--spots", type=int, default=20)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_identities)

    c = sub.add_parser("gridproof")
    c.add_argument("--a", type=int, default=2)
    c.add_argument("--b", type=int, default=1)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_gridproof)

    c = sub.add_parser("reduction")
    c.add_argument("--samples", type=int, default=300)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_reduction)

    c = sub.add_parser("casetree")
    c.add_argument("--samples", type=int, default=200000)
    c.add_argument("--margin", type=float, default=1e-9)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_casetree)

    c = sub.add_parser("alivecheck")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_alivecheck)

    c = sub.add_parser("death")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--uniform", type=int, default=1200)
    c.add_argument("--iters", type=int, default=46)
    c.add_argument("--full", action="store_true")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_death)

    c = sub.add_parser("deadhunt")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--uniform", type=int, default=3000)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_deadhunt)

    c = sub.add_parser("lemmac")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_lemmac)

    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
