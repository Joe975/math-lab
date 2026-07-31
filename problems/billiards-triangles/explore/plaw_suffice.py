#!/usr/bin/env python3
"""Attempt 005: certified SUFFICIENCY for the W(a,b) death law, per member.

003 proved (necessity): corridor alive => gamma < gamma_d(a,b).  This tool
proves the matching lower half per member: an explicit rational SEGMENT of
triangles, ending at the death corner, along which the corridor is certified
alive at every point, with gamma(t) -> gamma_d.  Together: death = gamma_d
EXACTLY (sup of alive gamma; the sup is not attained).

Construction.  Corner P0 = (alpha_c, beta_c) = (90/a, 90(a-1)/(a(b+1)))
(rational degrees; this is where 003 measured the last-alive apexes).  Take

    alpha(t) = alpha_c - t,   beta(t) = beta_c + rho * t,   t in (0, t0],

rho an integer >= 2, t0 rational.  Then theta(t) = theta_d + (rho-1) t and
gamma(t) = gamma_d - (rho-1) t, so gamma(t) increases to gamma_d as t -> 0+.

Alive criterion used (glide reduction, 003 sec. T3 + 004 sec. 1, converse
direction proved in the 005 record): the corridor has positive width iff the
axis offset m lies STRICTLY inside the projection interval of every
half-word gate (and tau != 0).  So we certify, for every distinct
gate-endpoint projection difference D = p(endpoint) - m, a constant sign on
the whole segment, plus per-gate opposite signs, plus Re(conj(delta) tau)
nonzero.  Two differences vanish AT the corner (t = 0) and need care:

  * p(A1) - m = N2 = cos(a alpha) sin((b+1) beta) sin(theta)   [identity I2]
    On the segment a*alpha = 90 - a t, so N2(t) = sin(a t) *
    sin((b+1)beta(t)) * sin(theta(t)): positive for t in (0, t0] by pure
    RATIONAL range checks (all three arguments stay inside (0, 180)).
  * m - p(C) = N1 = cos(at) sin(beta) sin(wh t) - sin(at) sin(alpha)
    sin(b beta)   [identity I1 with v = 90 - a alpha = at,
    w = alpha + (b+1)beta - 90 = wh t, wh = (b+1)rho - 1]
    N1(0) = 0 exactly; we certify N1' > 0 on [0, t1] (so N1 > 0 on (0, t1])
    and N1 > 0 on [t1, t0] by adaptive interval bisection.  rho is chosen so
    that N1'(0) = wh sin(beta_c) - a sin(alpha_c) sin(b beta_c) > 0
    (certified).
  * p(v2) - m = N4 = cos(a alpha) sin(alpha) sin(b beta) - sin(a alpha)
    sin(beta) cos(alpha + (b+1) beta)  [identity I4, attempt 005], where
    v2 = R_0(Rot_A(2a alpha) C) is the C-image endpoint of gate 2a+2.  On
    the segment N4(t) = sin(at) sin(alpha) sin(b beta) + cos(at) sin(beta)
    sin(wh t): BOTH products positive, again by rational range checks.

The three identities I1, I2, I4 are re-verified exactly (ring subtraction
in Q(i)[e^{i alpha/2}, e^{i beta/2}]) for the member before use; they are
also proven for ALL (a,b) in plaw_general.py.

All other differences are nonzero at t = 0 and are certified by adaptive
bisection over [0, t0].  Every enclosure is a single (non-iterated)
interval evaluation of an explicit trig polynomial at rational-degree
arguments (the affine-form lesson does not bite).  Certified trig uses
exact mod-360 degree reduction, dyadic pre-rounding, alternating-series
remainders, and the harness's rational pi enclosure; it is cross-checked
against harness unfold.sin_iv/cos_iv in --selftest.

Failure is loud: if any leg cannot be certified the member FAILS (t0 is
adaptively halved a bounded number of times first).

Usage:
  plaw_suffice.py selftest
  plaw_suffice.py member --a A --b B [--out FILE]
  plaw_suffice.py all [--out FILE]      (the 15 members of 003 + 5 new)
"""

import argparse
import json
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..",
                                "harness", "billiards-triangles"))
import unfold  # noqa: E402
import deathlaw_symbolic as ds  # noqa: E402
from skeptic_family import width as sk_width, family_word  # noqa: E402

Iv = unfold.Iv
PI = Iv(unfold.PI_LO, unfold.PI_HI)
DY = 2 ** 72                      # dyadic pre-rounding grid


# ------------------------------------------------- certified trig, degrees

_trig_cache = {}


def _sin_point_iv(x, terms=14):
    """Alternating-series enclosure of sin at rational x, |x| <= 4."""
    total, term = F(0), F(x)
    for k in range(terms):
        total += term
        term = -term * x * x / ((2 * k + 2) * (2 * k + 3))
    r = abs(term)
    return total - r, total + r


def _cos_point_iv(x, terms=14):
    total, term = F(0), F(1)
    for k in range(terms):
        total += term
        term = -term * x * x / ((2 * k + 1) * (2 * k + 2))
    r = abs(term)
    return total - r, total + r


def _trig_deg(fn_point, dlo, dhi):
    """Enclosure of sin/cos over the degree interval [dlo, dhi]."""
    key = (fn_point is _sin_point_iv, dlo, dhi)
    if key in _trig_cache:
        return _trig_cache[key]
    dm = (dlo + dhi) / 2
    # exact reduction of the midpoint mod 360 into [-180, 180]
    dm_red = dm - 360 * ((dm + 180) // 360)
    h = (dhi - dlo) / 2
    # radian midpoint interval (pi enclosure), then dyadic pre-round
    xlo = dm_red * (unfold.PI_LO if dm_red >= 0 else unfold.PI_HI) / 180
    xhi = dm_red * (unfold.PI_HI if dm_red >= 0 else unfold.PI_LO) / 180
    xm = F((xlo + xhi).numerator * DY // (xlo + xhi).denominator // 2, DY)
    slack = (xhi - xlo) / 2 + abs((xlo + xhi) / 2 - xm) + F(1, DY)
    plo, phi = fn_point(xm)
    # 1-Lipschitz in radians: widen by half-width h (degrees) in radians
    wid = slack + h * unfold.PI_HI / 180
    res = Iv(plo - wid, phi + wid)
    _trig_cache[key] = res
    return res


def sin_deg(dlo, dhi=None):
    return _trig_deg(_sin_point_iv, F(dlo), F(dlo if dhi is None else dhi))


def cos_deg(dlo, dhi=None):
    return _trig_deg(_cos_point_iv, F(dlo), F(dlo if dhi is None else dhi))


# ------------------------------------------------- ring element -> segment

def poly_terms(p):
    """Real (conj-invariant) ring element -> [(const?, m, n, re2, im2)].
    Value = c00 + sum over (m,n)>0 of [re2 cos(ang) + im2 sin(ang)],
    ang = (m alpha + n beta)/2 degrees, re2 = 2 Re c, im2 = 2 Im(-c)...
    We use: c e^{i ang} + conj(c) e^{-i ang} = 2 Re(c) cos - 2 Im(c) sin."""
    assert ds.pclean(ds.psub(p, ds.pconj(p))) == {}, "not a real element"
    const = F(0)
    terms = []
    for (m, n), c in p.items():
        if (m, n) == (0, 0):
            const += c[0]
        elif (m, n) > (0, 0):
            terms.append((m, n, 2 * c[0], -2 * c[1]))
    return const, terms


class SegPoly:
    """A real ring element restricted to the segment
    alpha = ac - t, beta = bc + rho t; evaluates over t-intervals."""

    def __init__(self, p, ac, bc, rho):
        self.const, raw = poly_terms(p)
        self.terms = []
        for (m, n, re2, im2) in raw:
            c0 = (m * ac + n * bc) / 2          # degrees at t = 0
            c1 = F(-m + n * rho, 2)             # slope in t
            self.terms.append((c0, c1, re2, im2))

    def eval_iv(self, tlo, thi):
        tot = Iv(self.const)
        for (c0, c1, re2, im2) in self.terms:
            a1, a2 = c0 + c1 * tlo, c0 + c1 * thi
            if a1 > a2:
                a1, a2 = a2, a1
            if re2:
                tot = tot + cos_deg(a1, a2) * re2
            if im2:
                tot = tot + sin_deg(a1, a2) * im2
        return tot


def certify_sign(evalf, sign, tlo, thi, max_depth=42):
    """Certify sign(evalf) == sign on [tlo, thi] by adaptive bisection.
    Returns leaf count or None on failure."""
    stack = [(F(tlo), F(thi))]
    leaves = 0
    span = F(thi) - F(tlo)
    while stack:
        lo, hi = stack.pop()
        e = evalf(lo, hi)
        if (sign > 0 and e.lo > 0) or (sign < 0 and e.hi < 0):
            leaves += 1
            continue
        if hi - lo < span / 2 ** max_depth:
            return None
        mid = (lo + hi) / 2
        stack.append((lo, mid))
        stack.append((mid, hi))
    return leaves


# ------------------------------------------------------------ the member

def n1_iv(a, b, ac, bc, rho, wh, tlo, thi):
    """N1(t) = cos(at) sin(bc+rho t) sin(wh t)
             - sin(at) sin(ac-t) sin(b bc + b rho t), certified."""
    t = (tlo, thi)

    def lin(c0, c1):
        a1, a2 = c0 + c1 * tlo, c0 + c1 * thi
        return (a1, a2) if a1 <= a2 else (a2, a1)

    return (cos_deg(*lin(F(0), F(a))) * sin_deg(*lin(bc, F(rho)))
            * sin_deg(*lin(F(0), wh))
            - sin_deg(*lin(F(0), F(a))) * sin_deg(*lin(ac, F(-1)))
            * sin_deg(*lin(b * bc, F(b * rho))))


def n1p_iv(a, b, ac, bc, rho, wh, tlo, thi):
    """dN1/dt (with the global pi/180 factor dropped: every term of the
    derivative carries exactly one such factor, so the sign is unchanged)."""
    def lin(c0, c1):
        a1, a2 = c0 + c1 * tlo, c0 + c1 * thi
        return (a1, a2) if a1 <= a2 else (a2, a1)

    at = lin(F(0), F(a))
    be = lin(bc, F(rho))
    wt = lin(F(0), wh)
    al = lin(ac, F(-1))
    bb = lin(b * bc, F(b * rho))
    t1p = (-a * sin_deg(*at) * sin_deg(*be) * sin_deg(*wt)
           + rho * cos_deg(*at) * cos_deg(*be) * sin_deg(*wt)
           + wh * cos_deg(*at) * sin_deg(*be) * cos_deg(*wt))
    t2p = (a * cos_deg(*at) * sin_deg(*al) * sin_deg(*bb)
           - sin_deg(*at) * cos_deg(*al) * sin_deg(*bb)
           + b * rho * sin_deg(*at) * sin_deg(*al) * cos_deg(*bb))
    return t1p - t2p


def run_member(a, b, verbose=True):
    theta_d = F(90 * (a + b), a * (b + 1))
    gamma_d = 180 - theta_d
    ac, bc = F(90, a), F(90 * (a - 1), a * (b + 1))
    rep = {"a": a, "b": b, "len": 4 * a + 4 * b + 2,
           "theta_d": str(theta_d), "gamma_d": str(gamma_d),
           "corner": [str(ac), str(bc)]}

    gd = ds.glide_data(a, b)          # asserts glide facts (exact, ring)
    dconj = ds.pconj(gd["delta"])
    proj = lambda z: ds.pim(ds.pmul(dconj, z))
    m = gd["m"]

    # distinct difference polynomials D = p(endpoint) - m, gate structure
    pC_m = ds.pclean(ds.psub(proj(gd["gates"][0][1]), m))   # p(C) - m = -N1
    pA1_m = ds.pclean(ds.psub(proj(gd["gates"][1][1]), m))  # p(A1) - m = N2
    pv2_m = ds.pclean(ds.psub(proj(gd["gates"][2 * a + 1][1]), m))  # = N4
    # exact per-member re-verification of I1, I2, I4 before use
    cos_of = lambda mm, nn: ds.pscale(
        ds.padd(ds.mono(mm, nn), ds.mono(-mm, -nn)), (F(1, 2), F(0)))
    S = ds.pmul(ds.pmul(cos_of(2 * a, 0), ds.sin_of(2, 0)),
                ds.sin_of(0, 2 * b))
    T = ds.pmul(ds.pmul(ds.sin_of(2 * a, 0), ds.sin_of(0, 2)),
                cos_of(2, 2 * (b + 1)))
    I1ok = ds.pclean(ds.psub(pC_m, ds.padd(S, T))) == {}   # pC-m = S+T=-N1
    I2m = ds.pmul(ds.pmul(cos_of(2 * a, 0), ds.sin_of(0, 2 * (b + 1))),
                  ds.sin_of(2, 2))
    I2ok = ds.pclean(ds.psub(pA1_m, I2m)) == {}
    I4ok = ds.pclean(ds.psub(pv2_m, ds.psub(S, T))) == {}  # pv2-m = S-T
    if not (I1ok and I2ok and I4ok):
        rep["ok"] = False
        rep["fail"] = f"identity re-check: I1={I1ok} I2={I2ok} I4={I4ok}"
        return rep
    diffs = {}                        # key -> poly
    gates_keys = []
    for k, (u, v) in enumerate(gd["gates"]):
        pair = []
        for z in (u, v):
            d = ds.pclean(ds.psub(proj(z), m))
            key = json.dumps(sorted((mn, (str(c[0]), str(c[1])))
                                    for mn, c in d.items()))
            diffs.setdefault(key, d)
            pair.append(key)
        gates_keys.append(pair)
    keyof = lambda p: json.dumps(sorted((mn, (str(c[0]), str(c[1])))
                                        for mn, c in p.items()))
    key_N1, key_N2, key_N4 = keyof(pC_m), keyof(pA1_m), keyof(pv2_m)
    assert key_N1 in diffs and key_N2 in diffs and key_N4 in diffs
    rep["n_distinct_diffs"] = len(diffs)

    # tau nondegeneracy: r(t) = Re(conj(delta) tau) as an extra function
    tau_r = ds.pre(ds.pmul(dconj, gd["tau"]))

    # choose rho: certified N1'(0) = wh sin(bc) - a sin(ac) sin(b bc) > 0
    rho = None
    for r in range(2, 65):
        wh = F((b + 1) * r - 1)
        e = (wh * sin_deg(bc) - a * sin_deg(ac) * sin_deg(b * bc))
        if e.lo > 0:
            rho = r
            break
    if rho is None:
        rep["ok"] = False
        rep["fail"] = "no rho with certified N1'(0) > 0"
        return rep
    wh = F((b + 1) * rho - 1)
    rep["rho"] = rho

    # initial t0 from the rational range constraints (N2 legs + obtuse)
    t0 = min(F(45, a), F(90, (b + 1) * rho),
             (90 - theta_d) / (2 * (rho - 1)), F(1, 4))

    for attempt in range(14):
        _trig_cache.clear()
        ok, detail = certify_segment(a, b, ac, bc, rho, wh, t0, diffs,
                                     gates_keys, key_N1, key_N2, key_N4,
                                     tau_r)
        if ok:
            rep.update(detail)
            rep["t0"] = str(t0)
            rep["ok"] = True
            # float cross-check via 002's independent corridor
            wd = []
            for frac in (F(1, 2), F(1, 97)):
                t = t0 * frac
                al, be = float(ac - t), float(bc + rho * t)
                ta = math.tan(math.radians(al))
                tb = math.tan(math.radians(be))
                x, y = tb / (ta + tb), ta * tb / (ta + tb)
                wd.append(sk_width(family_word(a, b), x, y))
            rep["float_widths_at_t0/2_t0/97"] = wd
            rep["float_alive_crosscheck"] = all(w > 0 for w in wd)
            rep["gamma_range_alive"] = [str(gamma_d - (rho - 1) * t0),
                                        f"-> {gamma_d} (open)"]
            rep["conclusion"] = (
                f"alive for all t in (0, {t0}]: gamma in "
                f"[{float(gamma_d - (rho - 1) * t0):.6f}, {float(gamma_d)})"
                f" all realized; with the parametric necessity theorem "
                f"(003 + 005 extended case tree): "
                f"death(W({a},{b})) = {gamma_d} EXACTLY")
            return rep
        if verbose:
            print(f"  W({a},{b}): t0={t0} failed ({detail}); halving")
        t0 /= 2
    rep["ok"] = False
    rep["fail"] = f"could not certify segment down to t0={t0}"
    return rep


def certify_segment(a, b, ac, bc, rho, wh, t0, diffs, gates_keys,
                    key_N1, key_N2, key_N4, tau_r):
    leaves = {}
    signs = {}
    # --- N2 leg: rational range checks only ---
    # at in (0, a t0] subset (0,90); (b+1)beta(t) in (90-90/a, +] subset
    # (0,180); theta(t) in (theta_d, theta_d+(rho-1)t0] subset (0,90)
    if not (a * t0 <= 45 and 90 - F(90, a) > 0
            and 90 - F(90, a) + (b + 1) * rho * t0 < 180
            and F(90 * (a + b), a * (b + 1)) + (rho - 1) * t0 < 90):
        return False, "N2 rational range check"
    signs[key_N2] = +1
    leaves["N2"] = 0

    # --- N4 leg: rational range checks only ---
    # N4(t) = sin(at) sin(alpha) sin(b beta) + cos(at) sin(beta) sin(wh t)
    # sin(at), cos(at) > 0 since at in (0,45]; alpha in (0, 90/a];
    # b beta in (b bc, b bc + b rho t0] subset (0,180); beta subset (0,180);
    # wh t in (0, wh t0] subset (0,180)
    if not (b * bc + b * rho * t0 < 180 and bc + rho * t0 < 180
            and wh * t0 < 180 and t0 < ac):
        return False, "N4 rational range check"
    signs[key_N4] = +1
    leaves["N4"] = 0

    # --- N1 legs ---
    t1 = t0 / 8
    for _ in range(20):
        n = certify_sign(lambda lo, hi: n1p_iv(a, b, ac, bc, rho, wh,
                                               lo, hi), +1, 0, t1)
        if n is not None:
            break
        t1 /= 2
    else:
        return False, "N1' > 0 near 0"
    leaves["N1_deriv"] = n
    n = certify_sign(lambda lo, hi: n1_iv(a, b, ac, bc, rho, wh, lo, hi),
                     +1, t1, t0)
    if n is None:
        return False, f"N1 > 0 on [{t1}, {t0}]"
    leaves["N1_bisect"] = n
    signs[key_N1] = -1                # p(C) - m = -N1 < 0

    # --- generic differences ---
    for key, p in diffs.items():
        if key in (key_N1, key_N2, key_N4):
            continue
        sp = SegPoly(p, ac, bc, rho)
        guess = sp.eval_iv(t0 / 2, t0 / 2)
        s = 1 if guess.lo > 0 else (-1 if guess.hi < 0 else 0)
        if s == 0:
            return False, "sign guess ambiguous"
        n = certify_sign(sp.eval_iv, s, 0, t0)
        if n is None:
            return False, f"generic diff sign on [0,{t0}]"
        signs[key] = s
        leaves.setdefault("generic", 0)
        leaves["generic"] += n

    # --- per-gate: strictly one endpoint above m, one below ---
    for pair in gates_keys:
        if signs[pair[0]] * signs[pair[1]] != -1:
            return False, "gate endpoints not straddling m"

    # --- tau nondegenerate along the segment ---
    sp = SegPoly(tau_r, ac, bc, rho)
    guess = sp.eval_iv(t0 / 2, t0 / 2)
    s = 1 if guess.lo > 0 else (-1 if guess.hi < 0 else 0)
    n = certify_sign(sp.eval_iv, s, 0, t0) if s else None
    if n is None:
        return False, "tau_r sign"
    leaves["tau_r"] = n
    return True, {"leaves": leaves}


# ------------------------------------------------------------------ CLI

MEMBERS = [(2, 1), (2, 2), (3, 2), (3, 3), (4, 2), (4, 3), (4, 4),
           (5, 3), (5, 4), (5, 5), (6, 6), (7, 7), (8, 8),
           (10, 9), (12, 12),
           (6, 3), (6, 1), (8, 2), (8, 1), (7, 1)]


def cmd_selftest(_args):
    ok = True
    import random
    rng = random.Random(20260731)
    for _ in range(300):
        d = F(rng.randint(-500000, 500000), rng.randint(1, 997))
        for mine, ref, fl in ((sin_deg(d), unfold.sin_iv(
                unfold.Iv(d) * PI / 180), math.sin(math.radians(float(d)))),
                (cos_deg(d), unfold.cos_iv(unfold.Iv(d) * PI / 180),
                 math.cos(math.radians(float(d))))):
            ok &= mine.lo <= ref.hi and ref.lo <= mine.hi  # overlap
            ok &= float(mine.lo) - 1e-12 <= fl <= float(mine.hi) + 1e-12
            ok &= float(mine.hi - mine.lo) < 1e-15
    print("[trig vs unfold + float x300]", "OK" if ok else "FAIL")
    # SegPoly vs direct float evaluation on a member
    gd = ds.glide_data(3, 2)
    ac, bc, rho = F(90, 3), F(90 * 2, 3 * 3), 3
    dconj = ds.pconj(gd["delta"])
    p = ds.pim(ds.pmul(dconj, gd["gates"][4][0]))
    pm = ds.pclean(ds.psub(p, gd["m"]))
    sp = SegPoly(pm, ac, bc, rho)
    ok2 = True
    for tt in (F(1, 100), F(1, 7), F(1, 3)):
        e = sp.eval_iv(tt, tt)
        al = math.radians(float(ac - tt))
        be = math.radians(float(bc + rho * tt))
        v = (ds.peval(p, al, be) - ds.peval(gd["m"], al, be)).real
        ok2 &= float(e.lo) - 1e-11 <= v <= float(e.hi) + 1e-11
    print("[SegPoly vs peval]", "OK" if ok2 else "FAIL")
    ok &= ok2
    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)
    c = sub.add_parser("member")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--out")
    c.set_defaults(fn=None)
    c = sub.add_parser("all")
    c.add_argument("--out")
    c.set_defaults(fn=None)
    args = ap.parse_args()
    if args.cmd == "selftest":
        args.fn(args)
        return
    if args.cmd == "member":
        results = [run_member(args.a, args.b)]
        print(json.dumps(results[0], indent=1))
    else:
        results = []
        for (a, b) in MEMBERS:
            r = run_member(a, b)
            results.append(r)
            print(f"W({a},{b}): rho={r.get('rho')} t0={r.get('t0')} "
                  f"diffs={r.get('n_distinct_diffs')} "
                  f"float_ok={r.get('float_alive_crosscheck')} "
                  f"=> {'DEATH = ' + r['gamma_d'] + ' EXACT' if r['ok'] else 'FAILED: ' + str(r.get('fail'))}")
            if args.out:                      # checkpoint after each member
                json.dump(results, open(args.out, "w"), indent=1)
    if args.out:
        json.dump(results, open(args.out, "w"), indent=1)
    ok = all(r["ok"] for r in results)
    print("ALL MEMBERS CERTIFIED" if ok else "SOME MEMBERS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
