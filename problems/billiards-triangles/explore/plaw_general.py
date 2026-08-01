#!/usr/bin/env python3
"""Attempt 005: PARAMETRIC proofs of 003's per-member inputs.

003 proved the death-law necessity theorem for W(a,b) = (0(12)^a(02)^b)^2
modulo two inputs that were machine-certified per member only:

  * the binding identities I1-I3 (+ glide facts) -- exact ring identities,
    verified for 15 members;
  * Lemma C -- a 1-D inequality, certified per member.

This tool closes the first gap FOR ALL INTEGERS a, b >= 1 by a change of
viewpoint: the half word composes to

    u = R_0 . (R_1 R_2)^a . (R_0 R_2)^b
      = R_0 . Rot_A(2a alpha) . Rot_B(-2b beta)

(composed-map unfolding, reflections across BASE sides; R_1 R_2 and R_0 R_2
collapse to rotations because the side pairs share a vertex).  Hence every
quantity in I1-I3 is a Laurent polynomial in the FOUR variables

    x = e^{i alpha},  y = e^{i beta},  P = e^{i a alpha},  Q = e^{i b beta},

with integer coefficients in Q(i), where a and b appear ONLY through P, Q.
An identity that holds formally in Q(i)[x^, y^, P^, Q^] therefore holds
for every integer a, b >= 1 and all real angles simultaneously: specializing
P -> x^a, Q -> y^b and then x -> e^{i alpha}, y -> e^{i beta} is a ring
homomorphism.  The finite formal check below IS the general proof of I1-I3,
modulo the (short, hand-written in the record) derivation of the closed
forms; that derivation is itself cross-checked exactly, member by member,
against deathlaw_symbolic's independent gate-by-gate unfolding.

Subcommands (all stdlib-only, deterministic):
  formal                formal 4-variable proof of I1, I2, I3, D1, D2 and
                        the glide facts (mu = delta^2, Im(conj(delta) tau)=0)
  specialize            exact per-member cross-check of the closed forms
                        (mu, w, delta, tau, m, p(B), p(C), p(A1)) against
                        deathlaw_symbolic.unfold_sym / glide_data
  floatcheck            float spot-check of the closed-form u against
                        composed reflections for random (a,b) up to 40
  lemmac                sanity grid for the GENERAL Lemma C proof (the proof
                        itself is elementary and lives in the record):
                        checks F(v) < 0 and H2 > 0 on dense grids
  casetree              adversarial sign-pattern search for the EXTENDED
                        case tree (scope a <= 3b+4): no Case I/II pattern
                        at theta <= theta_d for the given member
"""

import argparse
import cmath
import json
import math
import os
import random
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import deathlaw_symbolic as ds  # noqa: E402

# ---------------------------------------------------------- Q(i) scalars

ZERO = (F(0), F(0))
ONE = (F(1), F(0))
I = (F(0), F(1))
HALF = (F(1, 2), F(0))
HALF_OVER_I = (F(0), F(-1, 2))


def cadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def cmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def cneg(u):
    return (-u[0], -u[1])


def cconj(u):
    return (u[0], -u[1])


# ------------------------------- formal ring Q(i)[x^, y^, P^, Q^]
# poly: dict {(ex, ey, eP, eQ): (re, im)}; the monomial means
# x^ex y^ey P^eP Q^eQ with x=e^{i alpha}, y=e^{i beta}, P=e^{i a alpha},
# Q=e^{i b beta}.

def clean(p):
    return {k: c for k, c in p.items() if c != ZERO}


def add(p, q):
    r = dict(p)
    for k, c in q.items():
        r[k] = cadd(r.get(k, ZERO), c)
    return clean(r)


def neg(p):
    return {k: cneg(c) for k, c in p.items()}


def sub(p, q):
    return add(p, neg(q))


def mul(p, q):
    r = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = tuple(a + b for a, b in zip(k1, k2))
            r[k] = cadd(r.get(k, ZERO), cmul(c1, c2))
    return clean(r)


def scale(p, c):
    return clean({k: cmul(v, c) for k, v in p.items()})


def conj(p):
    """Complex conjugation: all four variables are unit complex numbers,
    so conj is exponent negation plus coefficient conjugation."""
    return {tuple(-e for e in k): cconj(c) for k, c in p.items()}


def im(p):
    return scale(sub(p, conj(p)), HALF_OVER_I)


def mono(ex, ey, eP, eQ, c=ONE):
    return {(ex, ey, eP, eQ): c}


def fsin(ex, ey, eP, eQ):
    """sin(ex*alpha + ey*beta + eP*a*alpha + eQ*b*beta) as a ring element."""
    return scale(sub(mono(ex, ey, eP, eQ), mono(-ex, -ey, -eP, -eQ)),
                 HALF_OVER_I)


def fcos(ex, ey, eP, eQ):
    return scale(add(mono(ex, ey, eP, eQ), mono(-ex, -ey, -eP, -eQ)), HALF)


def feval(p, alpha, beta, a, b):
    """Numeric value at radian angles for concrete integers a, b."""
    z = 0j
    for (ex, ey, eP, eQ), c in p.items():
        ang = ex * alpha + ey * beta + eP * a * alpha + eQ * b * beta
        z += complex(c[0], c[1]) * cmath.exp(1j * ang)
    return z


# --------------------------------------------- closed forms (the theorem)
# Derivation in the record (attempt 005): with A = 0, B = sin(alpha+beta),
# C = sin(beta) e^{i alpha} (circumdiameter 1) and base-side reflections
#   R_2: z -> conj z            (side AB = real axis)
#   R_1: z -> e^{2i alpha} conj z
#   R_0: z -> B + e^{-2i beta} conj(z - B)
# the composed-map unfolding gives u = R_0 (R_1 R_2)^a (R_0 R_2)^b
# = R_0 . Rot_A(2a alpha) . Rot_B(-2b beta), whence u(z) = mu conj(z) + w:

def closed_forms():
    B = fsin(1, 1, 0, 0)                        # sin(alpha+beta), real
    C = mul(fsin(0, 1, 0, 0), mono(1, 0, 0, 0))  # sin(beta) e^{i alpha}
    mu = mono(0, -2, -2, 2)                     # e^{-2i a alpha + 2i(b-1)beta}
    delta = mono(0, -1, -1, 1)                  # axis direction, delta^2 = mu
    #   w = B [1 + e^{-2i(a alpha + beta)} - mu - e^{-2i beta}]
    w = mul(B, add(sub(mono(0, 0, 0, 0), mono(0, -2, 0, 0)),
                   sub(mono(0, -2, -2, 0), mu)))
    tau = add(w, mul(mu, conj(w)))
    A1 = mul(B, sub(mono(0, 0, 0, 0), mono(0, -2, 0, 0)))  # R_0(A)=B(1-y^-2)
    dconj = conj(delta)
    proj = lambda z: im(mul(dconj, z))
    m = scale(im(mul(dconj, sub(w, scale(tau, HALF)))), HALF)
    return {"B": B, "C": C, "A1": A1, "mu": mu, "delta": delta, "w": w,
            "tau": tau, "m": m,
            "pB": proj(B), "pC": proj(C), "pA1": proj(A1)}


def cmd_formal(args):
    cf = closed_forms()
    checks = {}
    # glide facts
    checks["mu_eq_delta_sq"] = clean(
        sub(cf["mu"], mul(cf["delta"], cf["delta"]))) == {}
    checks["tau_parallel_axis"] = clean(
        im(mul(conj(cf["delta"]), cf["tau"]))) == {}
    # I1: m - p(C) = -[cos(a al) sin(al) sin(b be)
    #                  + sin(a al) sin(be) cos(al + (b+1) be)]
    I1 = neg(add(mul(mul(fcos(0, 0, 1, 0), fsin(1, 0, 0, 0)),
                     fsin(0, 0, 0, 1)),
                 mul(mul(fsin(0, 0, 1, 0), fsin(0, 1, 0, 0)),
                     fcos(1, 1, 0, 1))))
    checks["I1"] = clean(sub(sub(cf["m"], cf["pC"]), I1)) == {}
    # I2: p(A1) - m = cos(a al) sin((b+1) be) sin(al+be)
    I2 = mul(mul(fcos(0, 0, 1, 0), fsin(0, 1, 0, 1)), fsin(1, 1, 0, 0))
    checks["I2"] = clean(sub(sub(cf["pA1"], cf["m"]), I2)) == {}
    # I3: p(B) - m = cos((b+1) be) sin(a al) sin(al+be)
    I3 = mul(mul(fcos(0, 1, 0, 1), fsin(0, 0, 1, 0)), fsin(1, 1, 0, 0))
    checks["I3"] = clean(sub(sub(cf["pB"], cf["m"]), I3)) == {}
    # context identities (not load-bearing):
    # D2: p(A1) - p(C) = -sin(be) sin((a-1)al - (b+1)be)
    D2 = neg(mul(fsin(0, 1, 0, 0), fsin(-1, -1, 1, -1)))
    checks["D2"] = clean(sub(sub(cf["pA1"], cf["pC"]), D2)) == {}
    # D1: p(B) - p(C) = sin(al) sin(a al - b be)
    D1 = mul(fsin(1, 0, 0, 0), fsin(0, 0, 1, -1))
    checks["D1"] = clean(sub(sub(cf["pB"], cf["pC"]), D1)) == {}
    # I4 (needed for sufficiency): the C-image endpoint of gate 2a+2,
    # v2 = R_0(Rot_A(2a al) C) = B + y^-2 (P^-2 conj(C) - B), satisfies
    #   p(v2) - m = cos(a al) sin(al) sin(b be)
    #               - sin(a al) sin(be) cos(al + (b+1) be)
    v2 = add(cf["B"], mul(mono(0, -2, 0, 0),
                          sub(mul(mono(0, 0, -2, 0), conj(cf["C"])),
                              cf["B"])))
    pv2 = im(mul(conj(cf["delta"]), v2))
    I4 = sub(mul(mul(fcos(0, 0, 1, 0), fsin(1, 0, 0, 0)),
                 fsin(0, 0, 0, 1)),
             mul(mul(fsin(0, 0, 1, 0), fsin(0, 1, 0, 0)),
                 fcos(1, 1, 0, 1)))
    checks["I4"] = clean(sub(sub(pv2, cf["m"]), I4)) == {}
    sizes = {k: len(v) for k, v in cf.items()}
    ok = all(checks.values())
    out = {"checks": checks, "term_counts": sizes, "ok": ok,
           "meaning": ("each True is an exact identity in "
                       "Q(i)[x^,y^,P^,Q^]; holds for ALL integer a,b and "
                       "all real angles by specialization")}
    print(json.dumps(out, indent=1))
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    sys.exit(0 if ok else 1)


# ------------------------------------------- per-member exact cross-check

def specialize(p, a, b):
    """Ring hom to deathlaw_symbolic's half-angle ring: P -> x^a, Q -> y^b,
    then keys (m, n) meaning e^{i(m alpha + n beta)/2}."""
    r = {}
    for (ex, ey, eP, eQ), c in p.items():
        k = (2 * (ex + a * eP), 2 * (ey + b * eQ))
        r[k] = cadd(r.get(k, ZERO), c)
    return ds.pclean(r)


def check_member(a, b):
    from skeptic_family import family_word
    cf = closed_forms()
    word = family_word(a, b)
    half = word[:len(word) // 2]
    U = ds.unfold_sym(half)
    orient, lin, off = U["maps"][-1]
    gd = ds.glide_data(a, b)
    dconj = ds.pconj(gd["delta"])
    proj = lambda z: ds.pim(ds.pmul(dconj, z))
    g1u, g1v = gd["gates"][0]
    g2u, g2v = gd["gates"][1]
    res = {
        "orient_reversing": orient == -1,
        "mu": ds.pclean(ds.psub(specialize(cf["mu"], a, b), lin)) == {},
        "w": ds.pclean(ds.psub(specialize(cf["w"], a, b), off)) == {},
        "delta": ds.pclean(
            ds.psub(specialize(cf["delta"], a, b), gd["delta"])) == {},
        "tau": ds.pclean(ds.psub(specialize(cf["tau"], a, b),
                                 gd["tau"])) == {},
        "m": ds.pclean(ds.psub(specialize(cf["m"], a, b), gd["m"])) == {},
        "pB": ds.pclean(ds.psub(specialize(cf["pB"], a, b),
                                proj(g1u))) == {},
        "pC": ds.pclean(ds.psub(specialize(cf["pC"], a, b),
                                proj(g1v))) == {},
        "pA1": ds.pclean(ds.psub(specialize(cf["pA1"], a, b),
                                 proj(g2v))) == {},
    }
    # gate 2a+2's C-image endpoint is v2 = R_0(Rot_A(2a alpha) C)
    v2 = add(cf["B"], mul(mono(0, -2, 0, 0),
                          sub(mul(mono(0, 0, -2, 0), conj(cf["C"])),
                              cf["B"])))
    pv2 = im(mul(conj(cf["delta"]), v2))
    res["p_gate2a2_v"] = ds.pclean(
        ds.psub(specialize(pv2, a, b),
                proj(gd["gates"][2 * a + 1][1]))) == {}
    res["ok"] = all(res.values())
    return res


MEMBERS_003 = [(2, 1), (2, 2), (3, 2), (3, 3), (4, 2), (4, 3), (4, 4),
               (5, 3), (5, 4), (5, 5), (6, 6), (7, 7), (8, 8),
               (10, 9), (12, 12)]
MEMBERS_NEW = [(6, 3), (6, 1), (8, 2), (8, 1), (7, 1), (11, 2), (14, 4),
               (17, 9)]


def cmd_specialize(args):
    results = {}
    ok = True
    for (a, b) in MEMBERS_003 + MEMBERS_NEW:
        r = check_member(a, b)
        results[f"W({a},{b})"] = r
        ok &= r["ok"]
        print(f"W({a},{b}): {'MATCH' if r['ok'] else 'MISMATCH ' + str(r)}")
    out = {"results": results, "ok": ok,
           "meaning": ("closed forms of u = R_0 Rot_A(2a al) Rot_B(-2b be) "
                       "agree EXACTLY with the gate-by-gate symbolic "
                       "unfolding for every listed member")}
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    print("SPECIALIZE", "OK" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def cmd_floatcheck(args):
    """Closed-form u(z) vs composed float reflections, random (a,b)."""
    rng = random.Random(20260731)
    cf = closed_forms()
    worst = 0.0
    for trial in range(args.trials):
        a = rng.randint(1, 40)
        b = rng.randint(1, a)
        al = rng.uniform(0.02, 1.2)
        be = rng.uniform(0.02, min(1.2, math.pi - 0.05 - al))
        B = math.sin(al + be)
        z = complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
        # composed reflections R_0 (R_1 R_2)^a (R_0 R_2)^b applied to z
        zz = z
        for _ in range(b):                      # (R_0 R_2): rot about B
            zz = B + cmath.exp(-2j * be) * (zz - B)
        zz = cmath.exp(2j * a * al) * zz        # (R_1 R_2)^a: rot about A
        zz = B + cmath.exp(-2j * be) * (zz - B).conjugate()   # R_0
        mu = feval(cf["mu"], al, be, a, b)
        w = feval(cf["w"], al, be, a, b)
        err = abs(mu * z.conjugate() + w - zz)
        worst = max(worst, err)
    print(json.dumps({"trials": args.trials, "worst_error": worst,
                      "ok": worst < 1e-9}))
    sys.exit(0 if worst < 1e-9 else 1)


# ------------------------------------------------- Lemma C, general form

def F_logderiv(v_deg, a, b):
    """F(v) = 2/sin 2v - (1/a) cot((90-v)/a) - (b/a) cot(bv/a), degrees.
    The general proof (in the record) shows F < 0 for all a >= 2, b >= 1,
    v in (0, 90); this is a float sanity net only."""
    v = math.radians(v_deg)
    h = math.pi / 2
    return (2.0 / math.sin(2 * v)
            - (1.0 / a) / math.tan((h - v) / a)
            - (b / a) / math.tan(b * v / a))


def H2(v_deg, a, b):
    v0 = 90.0 / (b + 1)
    v, vv0 = math.radians(v_deg), math.radians(v0)
    return (math.sin(v) * math.sin(math.radians((90 - v_deg) / a))
            * math.cos(vv0)
            - math.sin(vv0) * math.sin(math.radians(b * v_deg / a))
            * math.cos(v))


def cmd_lemmac(args):
    bad = []
    n = 0
    for a in range(2, args.amax + 1):
        for b in range(1, a + 1):
            v0 = 90.0 / (b + 1)
            for k in range(args.grid):
                v = 90.0 * (k + 0.5) / args.grid      # F < 0 on (0, 90)
                n += 1
                if F_logderiv(v, a, b) >= 0:
                    bad.append(("F", a, b, v, F_logderiv(v, a, b)))
                if v < v0 and H2(v, a, b) <= 0:       # H2 > 0 on (0, v0)
                    bad.append(("H2", a, b, v, H2(v, a, b)))
    out = {"amax": args.amax, "grid": args.grid, "n_samples": n,
           "violations": bad[:20], "n_violations": len(bad),
           "ok": not bad}
    print(json.dumps(out))
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    sys.exit(0 if not bad else 1)


# --------------------------------- extended case tree adversarial search

def signs_NNN(al_deg, be_deg, a, b):
    """(N1, N2, N3) from the closed forms I1-I3 (floats)."""
    al, be = math.radians(al_deg), math.radians(be_deg)
    th = al + be
    N1 = -(math.cos(a * al) * math.sin(al) * math.sin(b * be)
           + math.sin(a * al) * math.sin(be) * math.cos(al + (b + 1) * be))
    N2 = math.cos(a * al) * math.sin((b + 1) * be) * math.sin(th)
    N3 = math.cos((b + 1) * be) * math.sin(a * al) * math.sin(th)
    return N1, N2, N3


def cmd_casetree(args):
    """No Case I (+++) or Case II (---) pattern may occur at
    theta <= theta_d.  Boundary-targeted + random, seeded."""
    a, b = args.a, args.b
    rng = random.Random(20260731 + 100 * a + b)
    theta_d = 90.0 * (a + b) / (a * (b + 1))
    tol = 1e-9
    hits = []
    n = 0
    # targeted alpha values: multiples of 90/a (branch bounds), corner,
    # plus (b+1)beta boundaries mapped to alpha = theta - beta
    targets = [90.0 * k / a for k in range(1, 4 * a) if 90.0 * k / a < 180]
    tsets = []
    for frac in (1.0, 1 - 1e-6, 1 - 1e-3, 0.9, 0.5, 0.1):
        th = theta_d * frac
        for t in targets:
            for eps in (-1e-7, -1e-4, 1e-4, 1e-7):
                al = t + eps
                if 0 < al < th:
                    tsets.append((al, th - al))
        for k in range(1, 4 * (b + 1)):
            be = 90.0 * k / (b + 1)
            for eps in (-1e-7, -1e-4, 1e-4, 1e-7):
                if 0 < be + eps < th:
                    tsets.append((th - be - eps, be + eps))
        for _ in range(args.random):
            al = rng.uniform(1e-6, th - 1e-6)
            tsets.append((al, th - al))
    for (al, be) in tsets:
        if not (al > 0 and be > 0):
            continue
        n += 1
        N1, N2, N3 = signs_NNN(al, be, a, b)
        if (N1 > tol and N2 > tol and N3 > tol) or \
           (N1 < -tol and N2 < -tol and N3 < -tol):
            hits.append((al, be, N1, N2, N3))
    # positive control: just above theta_d at the death corner
    alc = 90.0 / a
    thc = theta_d + 1e-4
    ctrl = signs_NNN(alc - 1e-5, thc - alc + 1e-5, a, b)
    ctrl_fires = ctrl[0] > 0 and ctrl[1] > 0 and ctrl[2] > 0
    out = {"a": a, "b": b, "theta_d": theta_d, "n_points": n,
           "case_hits_at_or_below_theta_d": hits[:10],
           "n_hits": len(hits),
           "positive_control_above_theta_d_fires": ctrl_fires,
           "ok": (not hits) and ctrl_fires}
    print(json.dumps(out))
    if args.out:
        json.dump(out, open(args.out, "w"), indent=1)
    sys.exit(0 if out["ok"] else 1)


def main():
    ap = argparse.ArgumentParser()
    sp = ap.add_subparsers(dest="cmd", required=True)
    c = sp.add_parser("formal")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_formal)
    c = sp.add_parser("specialize")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_specialize)
    c = sp.add_parser("floatcheck")
    c.add_argument("--trials", type=int, default=400)
    c.set_defaults(fn=cmd_floatcheck)
    c = sp.add_parser("lemmac")
    c.add_argument("--amax", type=int, default=60)
    c.add_argument("--grid", type=int, default=500)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_lemmac)
    c = sp.add_parser("casetree")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--random", type=int, default=20000)
    c.add_argument("--out")
    c.set_defaults(fn=cmd_casetree)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
