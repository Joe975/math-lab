#!/usr/bin/env python3
"""Attempt 008: skeptic review of 005 (plaw death-law completion).

Everything here is written from scratch for this review except where a 004
skeptic asset is explicitly reused (deathlaw_skeptic's interval trig stack
and its adaptive death sampler -- 004's tools are reserved for skeptics).
Nothing imports plaw_general / plaw_suffice, and nothing trusts
deathlaw_symbolic's ring in a verdict path (it appears only as the OBJECT
under test in `segment --vs-ds` cross-prints, never as evidence).

Subcommands (all stdlib-only, deterministic, seed 20260738):

  bridge     EXACT off-torus pair-algebra test of 005's closed-form bridge:
             the half word u = R_0 (R_1 R_2)^a (R_0 R_2)^b is composed by
             DIRECT reflection application (my own maps) and compared, as an
             exact Q(i) identity at random rational off-torus points, with
             the claimed closed forms mu, w, delta, tau, m, and the
             gate-(2a+2) C-endpoint v2 = R_0(Rot_A(2a alpha) C).  Also
             checks mu = delta^2, Im(conj(delta) tau) = 0 and the glide
             action p(u(z)) = 2m - p(z) at the same exact points.
  ring       my own formal 4-variable Laurent ring Q(i)[x,y,P,Q]:
             u is RE-DERIVED symbolically by affine-map composition (the
             only lattice steps being (x^2)^a -> P^2, (y^-2)^b -> Q^-2,
             each justified by a fixed-point identity checked in the ring),
             then I1, I2, I3, I4, D1, D2 and the glide facts are re-proven
             by exact polynomial subtraction in MY ring.
  lemma      adversarial numeric stress of Lemma L1 / general Lemma C / D:
             F(v) < 0 and H2 > 0 at huge a, b = 1, b = a, v near both
             endpoints; the factorization H2 = cos v cos v0 sin(bv/a)
             (Phi - Phi(v0)) checked numerically (my own radian code).
  threshold  exact Fraction check of the generalized threshold identity
             a*w - b*v = a(b+1)(theta - theta_d) - 360*k*b.
  casetree   my own adversarial sign-pattern scan (different design and
             seed from 005's) on members with a > 2b+3 and k >= 1 reachable:
             no (N1,N2,N3) Case I/II pattern may occur at theta <= theta_d;
             positive control above theta_d must fire.
  segment    independent alive certification at sampled t on 005's
             universal segment (alpha, beta) = (90/a - t, 90(a-1)/(a(b+1))
             + 2t): the FULL 2n-gate corridor is built by my own interval
             unfolding (004's certified trig; direct projection onto the
             normal of the full-word translation tau; no glide reduction
             used), and positive width is certified at t down to 2^-30.
  trigaudit  adversarial audit of plaw_suffice's certified-trig layer:
             its sin_deg/cos_deg enclosures must contain my independently
             computed values (004 stack) at random rational degrees and
             must be valid over intervals (random interior points inside).
"""

import argparse
import json
import math
import os
import random
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import deathlaw_skeptic as dsk  # 004's skeptic stack (mine to reuse)

SEED = 20260738


def family_half(a, b):
    return (0,) + (1, 2) * a + (0, 2) * b


def family_word(a, b):
    return family_half(a, b) * 2


# =====================================================================
# Gaussian rationals (re, im) as Fraction pairs
# =====================================================================

def gadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def gsub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def gmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def gconj(u):
    return (u[0], -u[1])


def ginv(u):
    n = u[0] * u[0] + u[1] * u[1]
    return (u[0] / n, -u[1] / n)


GZERO = (F(0), F(0))
GONE = (F(1), F(0))
GI = (F(0), F(1))


# =====================================================================
# off-torus pair algebra: pair (val, star-val); star = the involution
# x -> 1/x, y -> 1/y, i -> -i, which equals complex conjugation on the
# torus.  Exact evaluation of any expression built from +, -, *, /, conj.
# =====================================================================

def pr_const(c):
    return (c, gconj(c))


def pr_add(u, v):
    return (gadd(u[0], v[0]), gadd(u[1], v[1]))


def pr_sub(u, v):
    return (gsub(u[0], v[0]), gsub(u[1], v[1]))


def pr_mul(u, v):
    return (gmul(u[0], v[0]), gmul(u[1], v[1]))


def pr_div(u, v):
    return (gmul(u[0], ginv(v[0])), gmul(u[1], ginv(v[1])))


def pr_conj(u):
    return (u[1], u[0])


def pr_im(u):
    # (u - conj u) / (2i)
    d = pr_sub(u, pr_conj(u))
    inv2i = pr_const((F(0), F(-1, 2)))
    return pr_mul(d, inv2i)


def pr_var(x0):
    """A unit-modulus variable evaluated off-torus at rational x0."""
    return ((F(x0), F(0)), (1 / F(x0), F(0)))


def pr_pow(u, n):
    if n < 0:
        u = (ginv(u[0]), ginv(u[1]))
        n = -n
    r = (GONE, GONE)
    for _ in range(n):
        r = pr_mul(r, u)
    return r


def bridge_member(a, b, rng):
    """Exact pair-algebra comparison for one (a, b) at one random point."""
    # random off-torus rational points for x = e^{i alpha}, y = e^{i beta}
    def rq():
        n = rng.randint(2, 60)
        d = rng.randint(2, 60)
        while n == d:
            d = rng.randint(2, 60)
        return F(n, d) * rng.choice([1, -1])

    X, Y = pr_var(rq()), pr_var(rq())
    # geometry: A = 0, B = sin(alpha+beta), C = sin(beta) e^{i alpha}
    inv2i = pr_const((F(0), F(-1, 2)))
    XY = pr_mul(X, Y)
    B = pr_mul(pr_sub(XY, pr_div((GONE, GONE), XY)), inv2i)
    C = pr_mul(pr_mul(pr_sub(Y, pr_div((GONE, GONE), Y)), inv2i), X)
    Y2i = pr_pow(Y, -2)   # e^{-2 i beta}
    X2 = pr_pow(X, 2)

    # base reflections, re-derived by hand in the record:
    def R2(z):
        return pr_conj(z)

    def R1(z):
        return pr_mul(X2, pr_conj(z))

    def R0(z):
        return pr_add(B, pr_mul(Y2i, pr_conj(pr_sub(z, B))))

    R = {0: R0, 1: R1, 2: R2}
    half = family_half(a, b)

    def apply_chain(letters, z):
        # M = R_{s1} o ... o R_{sk};  M(z) applies R_{sk} first
        for s in reversed(letters):
            z = R[s](z)
        return z

    z1 = pr_const((F(3, 7), F(2, 5)))     # generic constant point
    w_dir = apply_chain(half, pr_const(GZERO))         # u(0) = w
    uz1 = apply_chain(half, z1)
    mu_dir = pr_div(pr_sub(uz1, w_dir), pr_conj(z1))

    # closed forms with P -> X^a, Q -> Y^b (the bridge under test)
    Pm2 = pr_pow(X, -2 * a)      # P^-2
    Q2 = pr_pow(Y, 2 * b)        # Q^2
    mu_cf = pr_mul(pr_mul(Pm2, Q2), Y2i)
    w_cf = pr_mul(B, pr_add(pr_sub((GONE, GONE), Y2i),
                            pr_sub(pr_mul(Pm2, Y2i), mu_cf)))
    ok_mu = mu_dir == mu_cf
    ok_w = w_dir == w_cf

    # v2 = R_0(Rot_A(2 a alpha) C) must equal the gate-(2a+2) C-endpoint
    # of the direct chain M_{2a+1} = R_0 (R_1 R_2)^a applied to C
    v2_dir = apply_chain(half[:2 * a + 1], C)
    v2_cf = pr_add(B, pr_mul(Y2i, pr_sub(pr_mul(Pm2, pr_conj(C)), B)))
    ok_v2 = v2_dir == v2_cf

    # glide facts at the exact point
    delta = pr_mul(pr_pow(X, -a), pr_pow(Y, b - 1))
    ok_delta = pr_mul(delta, delta) == mu_cf
    tau = pr_add(w_dir, pr_mul(mu_dir, pr_conj(w_dir)))
    dconj = pr_conj(delta)
    ok_tau_axis = pr_im(pr_mul(dconj, tau)) == (GZERO, GZERO)

    def proj(z):
        return pr_im(pr_mul(dconj, z))

    half_tau = pr_mul(tau, pr_const((F(1, 2), F(0))))
    m = pr_mul(pr_im(pr_mul(dconj, pr_sub(w_dir, half_tau))),
               pr_const((F(1, 2), F(0))))
    # glide action p(u(z)) + p(z) - 2m = 0 at the generic point z1
    lhs = pr_sub(pr_add(proj(uz1), proj(z1)),
                 pr_mul(m, pr_const((F(2), F(0)))))
    ok_glide = lhs == (GZERO, GZERO)
    return {"mu": ok_mu, "w": ok_w, "v2": ok_v2, "delta_sq": ok_delta,
            "tau_axis": ok_tau_axis, "glide_action": ok_glide,
            "ok": all([ok_mu, ok_w, ok_v2, ok_delta, ok_tau_axis,
                       ok_glide])}


def cmd_bridge(args):
    rng = random.Random(SEED)
    members = [(2, 1), (2, 2), (3, 2), (3, 3), (4, 2), (4, 3), (4, 4),
               (5, 3), (5, 4), (5, 5), (6, 6), (7, 7), (8, 8),
               (10, 9), (12, 12), (6, 3), (6, 1), (8, 2), (8, 1), (7, 1),
               (11, 2), (14, 4), (17, 9), (1, 1), (40, 1), (40, 40),
               (37, 13), (25, 24)]
    while len(members) < args.members:
        a = rng.randint(1, 40)
        b = rng.randint(1, a)
        members.append((a, b))
    results = {}
    ok = True
    for (a, b) in members:
        rr = [bridge_member(a, b, rng) for _ in range(args.points)]
        good = all(r["ok"] for r in rr)
        results[f"W({a},{b})"] = {"points": len(rr), "ok": good,
                                  "detail": rr[0]}
        ok &= good
        print(f"W({a},{b}): {'MATCH (exact, %d pts)' % len(rr) if good else 'MISMATCH ' + str(rr)}")
    out = {"seed": SEED, "n_members": len(members), "ok": ok,
           "meaning": ("direct reflection-composition of the half word "
                       "agrees EXACTLY (off-torus pair algebra, my own "
                       "maps) with 005's closed forms mu, w, v2 and glide "
                       "facts at random rational points")}
    print("BRIDGE", "OK" if ok else "FAILED")
    if args.out:
        json.dump({"results": results, **out}, open(args.out, "w"), indent=1)
    sys.exit(0 if ok else 1)


# =====================================================================
# my own formal 4-variable Laurent ring Q(i)[x^, y^, P^, Q^]
# =====================================================================

class R4:
    """Laurent polys as dict {(ex,ey,eP,eQ): Gaussian}."""

    @staticmethod
    def clean(p):
        return {k: c for k, c in p.items() if c != GZERO}

    @staticmethod
    def add(p, q):
        r = dict(p)
        for k, c in q.items():
            r[k] = gadd(r.get(k, GZERO), c)
        return R4.clean(r)

    @staticmethod
    def sub(p, q):
        return R4.add(p, {k: (-c[0], -c[1]) for k, c in q.items()})

    @staticmethod
    def mul(p, q):
        r = {}
        for k1, c1 in p.items():
            for k2, c2 in q.items():
                k = (k1[0] + k2[0], k1[1] + k2[1],
                     k1[2] + k2[2], k1[3] + k2[3])
                r[k] = gadd(r.get(k, GZERO), gmul(c1, c2))
        return R4.clean(r)

    @staticmethod
    def scale(p, c):
        return R4.clean({k: gmul(v, c) for k, v in p.items()})

    @staticmethod
    def conj(p):
        return {(-k[0], -k[1], -k[2], -k[3]): gconj(c)
                for k, c in p.items()}

    @staticmethod
    def im(p):
        return R4.scale(R4.sub(p, R4.conj(p)), (F(0), F(-1, 2)))

    @staticmethod
    def is_zero(p):
        return R4.clean(p) == {}


def r4mono(ex, ey, eP, eQ, c=GONE):
    return {(ex, ey, eP, eQ): c}


def r4sin(ex, ey, eP, eQ):
    return R4.scale(R4.sub(r4mono(ex, ey, eP, eQ),
                           r4mono(-ex, -ey, -eP, -eQ)), (F(0), F(-1, 2)))


def r4cos(ex, ey, eP, eQ):
    return R4.scale(R4.add(r4mono(ex, ey, eP, eQ),
                           r4mono(-ex, -ey, -eP, -eQ)), (F(1, 2), F(0)))


# antilinear/linear affine maps over the ring: ('L'|'A', lam, off)
# 'L': z -> lam z + off ; 'A': z -> lam conj(z) + off

def map_compose(f, g):
    """f o g."""
    tf, lf, of = f
    tg, lg, og = g
    if tf == 'L':
        return (tg, R4.mul(lf, lg), R4.add(R4.mul(lf, og), of))
    # f antilinear: f(g(z)) = lf conj(lg z + og) + of
    if tg == 'L':
        return ('A', R4.mul(lf, R4.conj(lg)),
                R4.add(R4.mul(lf, R4.conj(og)), of))
    return ('L', R4.mul(lf, R4.conj(lg)),
            R4.add(R4.mul(lf, R4.conj(og)), of))


def map_apply(f, z):
    t, l, o = f
    return R4.add(R4.mul(l, R4.conj(z) if t == 'A' else z), o)


def cmd_ring(args):
    checks = {}
    B = r4sin(1, 1, 0, 0)
    C = R4.mul(r4sin(0, 1, 0, 0), r4mono(1, 0, 0, 0))
    y2i = r4mono(0, -2, 0, 0)
    # base maps
    R2m = ('A', r4mono(0, 0, 0, 0), {})
    R1m = ('A', r4mono(2, 0, 0, 0), {})
    R0m = ('A', y2i, R4.sub(B, R4.mul(y2i, B)))   # B + y^-2 conj(z - B)
    # single blocks
    rot_A1 = map_compose(R1m, R2m)      # must be ('L', x^2, 0)
    rot_B1 = map_compose(R0m, R2m)      # must be ('L', y^-2, B(1 - y^-2))
    checks["rotA_block"] = (rot_A1[0] == 'L'
                            and R4.is_zero(R4.sub(rot_A1[1],
                                                  r4mono(2, 0, 0, 0)))
                            and R4.is_zero(rot_A1[2]))
    wB = R4.sub(B, R4.mul(y2i, B))
    checks["rotB_block"] = (rot_B1[0] == 'L'
                            and R4.is_zero(R4.sub(rot_B1[1], y2i))
                            and R4.is_zero(R4.sub(rot_B1[2], wB)))
    # fixed-point identities that justify the a-th / b-th power lattice
    # substitution: rot_A fixes A = 0 (trivial: off = 0, checked above);
    # rot_B fixes B: lam*B + off = B in the ring
    checks["rotB_fixes_B"] = R4.is_zero(
        R4.sub(R4.add(R4.mul(y2i, B), wB), B))
    # a-th power of z -> lam(z - fix) + fix is z -> lam^a (z - fix) + fix;
    # substitute lam^a via the exponent lattice: (x^2)^a = P^2,
    # (y^-2)^b = Q^-2.  These are the ONLY lattice steps.
    P2 = r4mono(0, 0, 2, 0)
    Qm2 = r4mono(0, 0, 0, -2)
    rot_A = ('L', P2, {})                                   # fixes 0
    rot_B = ('L', Qm2, R4.sub(B, R4.mul(Qm2, B)))           # fixes B
    u = map_compose(R0m, map_compose(rot_A, rot_B))
    checks["u_antilinear"] = u[0] == 'A'
    mu, w = u[1], u[2]
    # compare with 005's claimed closed forms
    mu_claim = r4mono(0, -2, -2, 2)
    w_claim = R4.mul(B, R4.add(R4.sub(r4mono(0, 0, 0, 0),
                                      r4mono(0, -2, 0, 0)),
                               R4.sub(r4mono(0, -2, -2, 0), mu_claim)))
    checks["mu_matches_005"] = R4.is_zero(R4.sub(mu, mu_claim))
    checks["w_matches_005"] = R4.is_zero(R4.sub(w, w_claim))
    # glide facts
    delta = r4mono(0, -1, -1, 1)
    checks["mu_eq_delta_sq"] = R4.is_zero(R4.sub(mu, R4.mul(delta, delta)))
    tau = R4.add(w, R4.mul(mu, R4.conj(w)))
    dconj = R4.conj(delta)
    checks["tau_parallel_axis"] = R4.is_zero(R4.im(R4.mul(dconj, tau)))

    def proj(z):
        return R4.im(R4.mul(dconj, z))

    m = R4.scale(R4.im(R4.mul(dconj,
                              R4.sub(w, R4.scale(tau, (F(1, 2), F(0)))))),
                 (F(1, 2), F(0)))
    A1 = map_apply(R0m, {})       # R_0(A), A = 0
    v2 = map_apply(R0m, map_apply(rot_A, C))
    pB, pC, pA1, pv2 = proj(B), proj(C), proj(A1), proj(v2)
    # the identities, RHS built from the record's prose (my own encoding)
    S = R4.mul(R4.mul(r4cos(0, 0, 1, 0), r4sin(1, 0, 0, 0)),
               r4sin(0, 0, 0, 1))          # cos(a al) sin(al) sin(b be)
    T = R4.mul(R4.mul(r4sin(0, 0, 1, 0), r4sin(0, 1, 0, 0)),
               r4cos(1, 1, 0, 1))    # sin(a al) sin(be) cos(al+(b+1)be)
    checks["I1"] = R4.is_zero(R4.sub(R4.sub(m, pC),
                                     R4.scale(R4.add(S, T),
                                              ((F(-1)), F(0)))))
    I2 = R4.mul(R4.mul(r4cos(0, 0, 1, 0), r4sin(0, 1, 0, 1)),
                r4sin(1, 1, 0, 0))
    checks["I2"] = R4.is_zero(R4.sub(R4.sub(pA1, m), I2))
    I3 = R4.mul(R4.mul(r4cos(0, 1, 0, 1), r4sin(0, 0, 1, 0)),
                r4sin(1, 1, 0, 0))
    checks["I3"] = R4.is_zero(R4.sub(R4.sub(pB, m), I3))
    checks["I4"] = R4.is_zero(R4.sub(R4.sub(pv2, m), R4.sub(S, T)))
    D1 = R4.mul(r4sin(1, 0, 0, 0), r4sin(0, 0, 1, -1))
    checks["D1"] = R4.is_zero(R4.sub(R4.sub(pB, pC), D1))
    D2 = R4.scale(R4.mul(r4sin(0, 1, 0, 0), r4sin(-1, -1, 1, -1)),
                  (F(-1), F(0)))
    checks["D2"] = R4.is_zero(R4.sub(R4.sub(pA1, pC), D2))
    ok = all(checks.values())
    print(json.dumps(checks, indent=1))
    print("RING", "OK (all identities re-proven in my own formal ring)"
          if ok else "FAILED")
    if args.out:
        json.dump({"checks": checks, "ok": ok, "seed": SEED},
                  open(args.out, "w"), indent=1)
    sys.exit(0 if ok else 1)


# =====================================================================
# Lemma L1 / C / D adversarial numerics (radians, my own code)
# =====================================================================

def cmd_lemma(args):
    rng = random.Random(SEED)
    bad = []
    n = 0
    # L1: d/dc [c cot(cs)] < 0 <=> sin(2cs) < 2cs; and monotonicity spot
    for _ in range(4000):
        s = rng.uniform(1e-6, math.pi / 2)
        c1 = rng.uniform(1e-9, 1.0)
        c2 = rng.uniform(1e-9, 1.0)
        if c1 == c2:
            continue
        c1, c2 = min(c1, c2), max(c1, c2)
        f1 = c1 / math.tan(c1 * s)
        f2 = c2 / math.tan(c2 * s)
        n += 1
        if not f1 > f2:
            bad.append(("L1", s, c1, c2, f1 - f2))
    # F(v) = 2/sin 2v - (1/a)cot((pi/2-v)/a) - (b/a)cot(bv/a) < 0
    # on (0, pi/2), a >= 2, 1 <= b <= a   (radians)
    def Fv(v, a, b):
        h = math.pi / 2
        return (2.0 / math.sin(2 * v) - (1.0 / a) / math.tan((h - v) / a)
                - (b / a) / math.tan(b * v / a))

    a_list = [2, 3, 4, 7, 20, 100, 1000, 10 ** 6]
    for a in a_list:
        for b in sorted({1, 2, max(1, a // 2), a - 1 if a > 1 else 1, a}):
            if b < 1 or b > a:
                continue
            v0 = math.pi / (2 * (b + 1))
            vs = [1e-9, 1e-6, 1e-3, 0.01, v0 * 0.5, v0 * (1 - 1e-9),
                  v0, min(math.pi / 2 - 1e-9, v0 * 1.5),
                  math.pi / 4, math.pi / 2 * (1 - 1e-6),
                  math.pi / 2 * (1 - 1e-12)]
            vs += [rng.uniform(1e-8, math.pi / 2 - 1e-12)
                   for _ in range(200)]
            for v in vs:
                if not (0 < v < math.pi / 2):
                    continue
                n += 1
                fv = Fv(v, a, b)
                if not fv < 0:
                    bad.append(("F", a, b, v, fv))
            # H2 > 0 on (0, v0), via the FACTORED form (cancellation-free)
            for v in vs:
                if not (0 < v < v0):
                    continue
                n += 1
                phi = (math.tan(v) * math.sin((math.pi / 2 - v) / a)
                       / math.sin(b * v / a))
                phi0 = math.tan(v0)
                h2f = (math.cos(v) * math.cos(v0) * math.sin(b * v / a)
                       * (phi - phi0))
                if not h2f > 0:
                    bad.append(("H2fact", a, b, v, h2f))
                # direct H2 must agree with the factored form (identity)
                h2d = (math.sin(v) * math.sin((math.pi / 2 - v) / a)
                       * math.cos(v0)
                       - math.sin(v0) * math.sin(b * v / a) * math.cos(v))
                if abs(h2d - h2f) > 1e-9 * max(1.0, abs(h2d), abs(h2f)):
                    bad.append(("H2ident", a, b, v, h2d - h2f))
    # Lemma D: b cot(bx) - cot x < 0 on (0, pi/(2b)], b >= 2.  Tested in
    # the equivalent positive-quantity form tan(bx) > b tan(x) (valid
    # since both tangents are positive there), sampling x relative to the
    # domain end pi/(2b); x below ~1e-3 of the domain is float-untestable
    # (difference ~ (bx)^3/3 underflows relative precision) and is covered
    # by the hand proof via L1.
    for b in [2, 3, 5, 17, 400, 10 ** 5]:
        dom = math.pi / (2 * b)
        for lam in ([1e-3, 1e-2, 0.1, 0.5, 0.9, 0.999, 1.0]
                    + [rng.uniform(1e-3, 1.0) for _ in range(200)]):
            x = lam * dom
            n += 1
            g = math.tan(min(b * x, math.pi / 2 * (1 - 1e-12))) \
                - b * math.tan(x)
            if not g > 0:
                bad.append(("D", b, x, g))
    out = {"seed": SEED, "n_checks": n, "violations": bad[:20],
           "n_violations": len(bad), "ok": not bad}
    print(json.dumps(out))
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    sys.exit(0 if not bad else 1)


# =====================================================================
# threshold identity, exact
# =====================================================================

def cmd_threshold(args):
    rng = random.Random(SEED)
    ok = True
    for _ in range(2000):
        a = rng.randint(1, 60)
        b = rng.randint(1, a)
        k = rng.randint(0, 5)
        al = F(rng.randint(1, 10 ** 6), rng.randint(1, 997))
        be = F(rng.randint(1, 10 ** 6), rng.randint(1, 997))
        theta = al + be
        theta_d = F(90 * (a + b), a * (b + 1))
        v = 360 * k + 90 - a * al
        w = al + (b + 1) * be - 90
        ok &= (a * w - b * v == a * (b + 1) * (theta - theta_d)
               - 360 * k * b)
    print("THRESHOLD identity a*w - b*v == a(b+1)(theta-theta_d) - 360kb:",
          "OK (2000 exact rational checks)" if ok else "FAILED")
    if args.out:
        json.dump({"ok": ok, "n": 2000, "seed": SEED},
                  open(args.out, "w"), indent=1)
    sys.exit(0 if ok else 1)


# =====================================================================
# my own case-tree adversarial scan
# =====================================================================

def signs(al, be, a, b):
    N1 = -(math.cos(a * al) * math.sin(al) * math.sin(b * be)
           + math.sin(a * al) * math.sin(be)
           * math.cos(al + (b + 1) * be))
    N2 = math.cos(a * al) * math.sin((b + 1) * be) * math.sin(al + be)
    N3 = math.cos((b + 1) * be) * math.sin(a * al) * math.sin(al + be)
    return N1, N2, N3


def cmd_casetree(args):
    a, b = args.a, args.b
    rng = random.Random(SEED + 1000 * a + b)
    theta_d = 90.0 * (a + b) / (a * (b + 1))
    d2r = math.pi / 180.0
    hits, nearhits = [], []
    n = 0
    thetas = [theta_d * (1 - g)
              for g in (0.0, 1e-12, 1e-9, 1e-6, 1e-4, 1e-3,
                        0.01, 0.05, 0.2, 0.5, 0.8)]
    for th in thetas:
        cands = []
        # a*alpha near every multiple of 90 (all k, incl. aal > 360)
        kmax = int(a * th / 90.0) + 2
        for k in range(1, kmax + 1):
            base = 90.0 * k / a
            for eps in (1e-10, 1e-7, 1e-4, 1e-2, 0.5):
                cands += [base - eps, base + eps]
        # (b+1)*beta near every multiple of 90
        jmax = int((b + 1) * th / 90.0) + 2
        for j in range(1, jmax + 1):
            be = 90.0 * j / (b + 1)
            for eps in (1e-10, 1e-7, 1e-4, 1e-2, 0.5):
                cands += [th - (be - eps), th - (be + eps)]
        # death corner accumulation
        for eps in (1e-12, 1e-9, 1e-6, 1e-3, 0.1):
            cands += [90.0 / a - eps, 90.0 / a + eps]
        # random fill
        cands += [rng.uniform(1e-8, th - 1e-8) for _ in range(args.random)]
        for al in cands:
            be = th - al
            if not (al > 1e-12 and be > 1e-12):
                continue
            n += 1
            N1, N2, N3 = signs(al * d2r, be * d2r, a, b)
            tol = args.margin
            if (N1 > tol and N2 > tol and N3 > tol) or \
               (N1 < -tol and N2 < -tol and N3 < -tol):
                hits.append((th, al, N1, N2, N3))
            elif (N1 > 0 and N2 > 0 and N3 > 0) or \
                 (N1 < 0 and N2 < 0 and N3 < 0):
                nearhits.append((th, al, N1, N2, N3))
    # positive control above theta_d
    alc = 90.0 / a
    thc = theta_d + 1e-4
    c = signs((alc - 1e-5) * d2r, (thc - alc + 1e-5) * d2r, a, b)
    ctrl = c[0] > 0 and c[1] > 0 and c[2] > 0
    out = {"a": a, "b": b, "theta_d": theta_d, "seed": SEED,
           "n_points": n, "n_hits": len(hits), "hits": hits[:10],
           "n_nearhits_below_margin": len(nearhits),
           "nearhits": nearhits[:10],
           "positive_control_fires": ctrl,
           "ok": (not hits) and ctrl}
    print(json.dumps({k: v for k, v in out.items() if k != "nearhits"}))
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    sys.exit(0 if out["ok"] else 1)


# =====================================================================
# independent segment alive certification (interval unfolding, full word)
# =====================================================================

CHOP = 2 ** 320       # outward dyadic rounding grid (sound: only widens)


def chop_iv(x):
    lo = F(math.floor(x.lo * CHOP), CHOP)
    hi = F(math.ceil(x.hi * CHOP), CHOP)
    return dsk.MyIv(lo, hi)


def iv_c(re, im):
    return (re, im)     # complex interval = (MyIv, MyIv)


def ivc_chop(u):
    return (chop_iv(u[0]), chop_iv(u[1]))


def ivc_add(u, v):
    return (u[0] + v[0], u[1] + v[1])


def ivc_sub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def ivc_mul(u, v):
    return ivc_chop((u[0] * v[0] - u[1] * v[1],
                     u[0] * v[1] + u[1] * v[0]))


def ivc_conj(u):
    return (u[0], -u[1])


def certify_alive_at(a, b, t):
    """Certify positive full-corridor width of W(a,b) at
    (alpha, beta) = (90/a - t, 90(a-1)/(a(b+1)) + 2t), rational degrees,
    with my own interval unfolding (004 trig stack).  Returns
    (certified_bool, width_lo_float)."""
    MyIv, miv = dsk.MyIv, dsk.miv
    ac, bc = F(90, a), F(90 * (a - 1), a * (b + 1))
    al, be = ac - t, bc + 2 * t
    th = al + be
    sd = lambda d: chop_iv(dsk.sin_deg(d))
    cd = lambda d: chop_iv(dsk.cos_deg(d))
    # vertices (circumdiameter 1): A = 0, B = sin th, C = sin be e^{i al}
    zero = (miv(0), miv(0))
    A = zero
    B = (sd(th), miv(0))
    C = (sd(be) * cd(al), sd(be) * sd(al))
    V = [A, B, C]
    # direction-squared of sides: E0 = e^{-2i be}, E1 = e^{2i al}, E2 = 1
    Es = {0: (cd(2 * be), -sd(2 * be)),
          1: (cd(2 * al), sd(2 * al)),
          2: (miv(1), miv(0))}
    SIDE_ENDS = {0: (1, 2), 1: (2, 0), 2: (0, 1)}
    cs = {s: ivc_sub(V[SIDE_ENDS[s][0]],
                     ivc_mul(Es[s], ivc_conj(V[SIDE_ENDS[s][0]])))
          for s in range(3)}
    word = family_word(a, b)
    sigma, L, toff = 1, (miv(1), miv(0)), zero
    gates = []
    for s in word:
        i, j = SIDE_ENDS[s]
        if sigma == 1:
            gates.append((ivc_add(ivc_mul(L, V[i]), toff),
                          ivc_add(ivc_mul(L, V[j]), toff)))
            L, toff = ivc_mul(L, Es[s]), ivc_add(ivc_mul(L, cs[s]), toff)
            sigma = -1
        else:
            gates.append((ivc_add(ivc_mul(L, ivc_conj(V[i])), toff),
                          ivc_add(ivc_mul(L, ivc_conj(V[j])), toff)))
            newL = ivc_mul(L, ivc_conj(Es[s]))
            toff = ivc_add(ivc_mul(L, ivc_conj(cs[s])), toff)
            L = newL
            sigma = 1
    assert sigma == 1
    # translation word: L must enclose 1
    assert L[0].lo <= 1 <= L[0].hi and L[1].lo <= 0 <= L[1].hi
    tau = toff
    # tau certified nonzero?
    tau_nonzero = (tau[0].lo > 0 or tau[0].hi < 0
                   or tau[1].lo > 0 or tau[1].hi < 0)
    # projection p(z) = Im(conj(tau) z)
    lo_hi, hi_lo = None, None      # sup of interval-lo of gate max, etc.
    LO_hi = None                   # upper bound on corridor lo
    HI_lo = None                   # lower bound on corridor hi
    LO_lo = None
    HI_hi = None
    for (zu, zv) in gates:
        pu = ivc_mul(ivc_conj(tau), zu)[1]
        pv = ivc_mul(ivc_conj(tau), zv)[1]
        gmin_hi = min(pu.hi, pv.hi)
        gmin_lo = min(pu.lo, pv.lo)
        gmax_lo = max(pu.lo, pv.lo)
        gmax_hi = max(pu.hi, pv.hi)
        LO_hi = gmin_hi if LO_hi is None else max(LO_hi, gmin_hi)
        LO_lo = gmin_lo if LO_lo is None else max(LO_lo, gmin_lo)
        HI_lo = gmax_lo if HI_lo is None else min(HI_lo, gmax_lo)
        HI_hi = gmax_hi if HI_hi is None else min(HI_hi, gmax_hi)
    width_lo = HI_lo - LO_hi         # certified lower bound on width
    width_hi = HI_hi - LO_lo
    return (tau_nonzero and width_lo > 0), float(width_lo), float(width_hi)


def cmd_segment(args):
    members = [(2, 1), (6, 3), (6, 1), (8, 2), (8, 1), (7, 1), (12, 12)]
    ts = [F(1, 4), F(1, 8), F(1, 97), F(1, 1024), F(1, 2 ** 20),
          F(1, 2 ** 30)]
    rows = []
    ok = True
    from skeptic_family import width as sk_width
    for (a, b) in members:
        ac, bc = F(90, a), F(90 * (a - 1), a * (b + 1))
        for t in ts:
            good, wlo, whi = certify_alive_at(a, b, t)
            # float cross-check at the apex (002's corridor)
            al = math.radians(float(ac - t))
            be = math.radians(float(bc + 2 * t))
            ta, tb = math.tan(al), math.tan(be)
            x, y = tb / (ta + tb), ta * tb / (ta + tb)
            wfl = sk_width(family_word(a, b), x, y)
            rows.append({"a": a, "b": b, "t": str(t),
                         "gamma_minus_gamma_d": float(-t),
                         "certified_alive": good,
                         "width_iv": [wlo, whi], "float_width": wfl})
            ok &= good
            print(f"W({a},{b}) t={t}: certified_alive={good} "
                  f"width>[{wlo:.3e}] float={wfl:.3e}")
    out = {"rows": rows, "ok": ok, "seed": SEED,
           "meaning": ("full 2n-gate corridor built by my own interval "
                       "unfolding (no glide reduction), certified "
                       "positive width at every sampled t on 005's "
                       "universal segment incl. t = 2^-30 "
                       "(gamma within 1e-9 of gamma_d)")}
    print("SEGMENT", "OK" if ok else "FAILED")
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    sys.exit(0 if ok else 1)


# =====================================================================
# audit of plaw_suffice's certified trig layer
# =====================================================================

def reduce360(d):
    """Exact reduction of rational degrees into [-180, 180)."""
    d = F(d)
    return d - 360 * ((d + 180) // 360)


def cmd_trigaudit(args):
    """plaw_suffice's sin_deg/cos_deg enclosures must CONTAIN the true
    value.  Reference: 004's independent stack (dsk.sin_deg/cos_deg,
    Machin pi, own series) -- valid only for small arguments, so the
    rational degree argument is exactly reduced mod 360 first (sin/cos
    are 360-periodic, so the true value is unchanged; plaw_suffice's own
    internal reduction is part of what is under test, so IT gets the raw
    argument).  My reduced enclosure has width ~1e-84, so requiring it to
    sit inside plaw's enclosure certifies true-value containment."""
    import plaw_suffice as ps      # the code UNDER TEST
    rng = random.Random(SEED)
    ok = True
    worst = 0.0
    n = 0
    for _ in range(args.trials):
        num = rng.randint(-10 ** 7, 10 ** 7)
        den = rng.randint(1, 9973)
        d = F(num, den)
        dr = reduce360(d)
        for (name, theirs, mine) in (
                ("sin", ps.sin_deg(d), dsk.sin_deg(dsk.MyIv(dr))),
                ("cos", ps.cos_deg(d), dsk.cos_deg(dsk.MyIv(dr)))):
            n += 1
            if not (theirs.lo <= mine.lo and mine.hi <= theirs.hi):
                ok = False
                print(f"{name} POINT violation at d={d}")
            worst = max(worst, float(theirs.hi - theirs.lo))
    # interval validity: for random degree intervals, every sampled
    # interior rational point's true value must lie inside the enclosure
    for _ in range(args.trials // 4):
        lo = F(rng.randint(-10 ** 5, 10 ** 5), rng.randint(1, 997))
        hi = lo + F(rng.randint(1, 1000), rng.randint(1, 997))
        e_s = ps.sin_deg(lo, hi)
        e_c = ps.cos_deg(lo, hi)
        for _ in range(4):
            lam = F(rng.randint(0, 1000), 1000)
            d = lo + (hi - lo) * lam
            dr = reduce360(d)
            mys = dsk.sin_deg(dsk.MyIv(dr))
            myc = dsk.cos_deg(dsk.MyIv(dr))
            n += 2
            if not (e_s.lo <= mys.lo and mys.hi <= e_s.hi):
                ok = False
                print("SIN interval violation", lo, hi, d)
            if not (e_c.lo <= myc.lo and myc.hi <= e_c.hi):
                ok = False
                print("COS interval violation", lo, hi, d)
    out = {"n_checks": n, "worst_point_width": worst, "ok": ok,
           "seed": SEED,
           "note": ("earlier FAILED run used an UNREDUCED reference "
                    "argument (my stack has no mod-360 reduction and its "
                    "series is invalid past ~660 deg) -- an artifact of "
                    "my audit, not of plaw_suffice; see record")}
    print(json.dumps(out))
    print("TRIGAUDIT", "OK" if ok else "FAILED")
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    sys.exit(0 if ok else 1)


def miv_frac(d):
    return dsk.MyIv(F(d))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("bridge")
    c.add_argument("--members", type=int, default=40)
    c.add_argument("--points", type=int, default=2)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_bridge)
    c = sub.add_parser("ring")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_ring)
    c = sub.add_parser("lemma")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_lemma)
    c = sub.add_parser("threshold")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_threshold)
    c = sub.add_parser("casetree")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--random", type=int, default=4000)
    c.add_argument("--margin", type=float, default=1e-9)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_casetree)
    c = sub.add_parser("segment")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_segment)
    c = sub.add_parser("trigaudit")
    c.add_argument("--trials", type=int, default=200)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_trigaudit)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
