#!/usr/bin/env python3
"""Closed forms for ALL W(a,b) gate margins (queue 11, from 005 leads 1+3).

005 proved (Step 3 + I1-I4) that the half-word prefix maps are
R_0 Rot_A(2s alpha) and R_0 Rot_A(2a alpha) Rot_B(-2j beta), and certified
the corridor per member by interval bisection on the ~2a+2b+1 gate-endpoint
projection differences.  This module closes the structural gap: every
distinct difference D = p(endpoint) - m has an explicit closed form,

    A-fan pivot   N2        = p(A1) - m = cos(a al) sin((b+1) be) sin(th)
    gate-1 far    PB        = p(B)  - m = cos((b+1) be) sin(a al) sin(th)
    A-fan C-image (s=0..a-1) N2 + sin(be) sin((a-2s-1) al - (b+1) be)
    A-fan B-image (s=1..a)   N2 + sin(th) sin((a-2s)   al - (b+1) be)
    B-fan C-image (j=0..b-1) -PB + sin(al) sin(a al + (b-2j)   be)
    B-fan A-image (j=1..b)   -PB + sin(th) sin(a al + (b+1-2j) be)

with th = al + be, and the B-fan pivot p(B1) - m = -PB (an exact
antisymmetry: p(B) + p(B1) = 2m).  The B-image s=a IS the B-fan pivot
(same unfolded point), which is the sin-addition identity
-I3 = I2 - sin(th) sin(a al + (b+1) be).

Derivation (checked exactly in the ring by `selftest`): the far endpoint of
an A-fan gate is z_s = R_0 Rot_A(2s al) X, X in {C, B}, so
z_s - A1 = e0 e^{-2is al} conj(X) (A = 0, R_0 w = e0 conj(w) + c0,
e0 = e^{-2i be}); projecting with delta = e^{i(-a al + (b-1) be)} gives the
sinusoid.  B-fan far endpoints factor as z'_j - B1 =
-sin(al) e^{i(-2a al + (2j-1) be)} (C-images) and
-sin(th) e^{i(-2a al + 2(j-1) be)} (A-images).

Alive criterion (glide reduction, 003/004/005: positive corridor width
<=> m strictly inside every half-word gate projection, tau != 0):
    PB > 0, N2 > 0, all A-fan far margins < 0, all B-fan far margins > 0.
(N1 > 0 and N4 > 0 are the s=0 C-image < 0 and j=0 C-image > 0 rows.)

Commands (stdlib only, deterministic):
    selftest              exact ring verification of every closed form for
                          a sweep of members + float corridor cross-check
    segment --amax A      margin sign stress over the universal segment
    alive --a A --b B --alpha X --beta Y     evaluate margins at a point
"""

import argparse
import math
import sys

D2R = math.pi / 180.0


def margins(a, b, al, be):
    """All distinct gate margins at (al, be) in DEGREES.

    Returns list of (label, value, need_sign) with need_sign in {+1, -1}:
    alive at (al, be)  <=>  every value has its need_sign strictly.
    """
    th = al + be
    s_al, s_be, s_th = (math.sin(x * D2R) for x in (al, be, th))
    N2 = (math.cos(a * al * D2R) * math.sin((b + 1) * be * D2R) * s_th)
    PB = (math.cos((b + 1) * be * D2R) * math.sin(a * al * D2R) * s_th)
    out = [("PB", PB, +1), ("N2", N2, +1)]
    for s in range(0, a):
        arg = ((a - 2 * s - 1) * al - (b + 1) * be) * D2R
        out.append((f"AC{s}", N2 + s_be * math.sin(arg), -1))
    for s in range(1, a + 1):
        arg = ((a - 2 * s) * al - (b + 1) * be) * D2R
        out.append((f"AB{s}", N2 + s_th * math.sin(arg), -1))
    for j in range(0, b):
        arg = (a * al + (b - 2 * j) * be) * D2R
        out.append((f"BC{j}", -PB + s_al * math.sin(arg), +1))
    for j in range(1, b + 1):
        arg = (a * al + (b + 1 - 2 * j) * be) * D2R
        out.append((f"BA{j}", -PB + s_th * math.sin(arg), +1))
    return out


def alive(a, b, al, be):
    """Strictly-inside test via the closed forms (float)."""
    return all(v * s > 0 for _, v, s in margins(a, b, al, be))


def min_signed_margin(a, b, al, be):
    """min over margins of (need_sign * value): positive iff alive."""
    return min(v * s for _, v, s in margins(a, b, al, be))


# --------------------------------------------------------------- selftest

def _ring_check(a, b):
    """Exact verification of every closed form against deathlaw_symbolic."""
    import deathlaw_symbolic as ds
    gd = ds.glide_data(a, b)
    dconj = ds.pconj(gd["delta"])
    proj = lambda z: ds.pim(ds.pmul(dconj, z))
    m = gd["m"]
    diffs = []
    seen = set()
    import json as _json
    for (u, v) in gd["gates"]:
        for z in (u, v):
            d = ds.pclean(ds.psub(proj(z), m))
            key = _json.dumps(sorted((mn, (str(c[0]), str(c[1])))
                                     for mn, c in d.items()))
            if key not in seen:
                seen.add(key)
                diffs.append(d)

    cos_of = lambda mm, nn: ds.pscale(
        ds.padd(ds.mono(mm, nn), ds.mono(-mm, -nn)), (0.5).__class__ and
        __import__("fractions").Fraction and (
            __import__("fractions").Fraction(1, 2),
            __import__("fractions").Fraction(0)))
    from fractions import Fraction as F
    cos_of = lambda mm, nn: ds.pscale(
        ds.padd(ds.mono(mm, nn), ds.mono(-mm, -nn)), (F(1, 2), F(0)))
    sin_of = ds.sin_of

    N2 = ds.pmul(ds.pmul(cos_of(2 * a, 0), sin_of(0, 2 * (b + 1))),
                 sin_of(2, 2))
    PB = ds.pmul(ds.pmul(cos_of(0, 2 * (b + 1)), sin_of(2 * a, 0)),
                 sin_of(2, 2))
    forms = [N2, PB]
    for s in range(0, a):
        forms.append(ds.padd(N2, ds.pmul(
            sin_of(0, 2), sin_of(2 * (a - 2 * s - 1), -2 * (b + 1)))))
    for s in range(1, a + 1):
        forms.append(ds.padd(N2, ds.pmul(
            sin_of(2, 2), sin_of(2 * (a - 2 * s), -2 * (b + 1)))))
    for j in range(0, b):
        forms.append(ds.psub(ds.pmul(
            sin_of(2, 0), sin_of(2 * a, 2 * (b - 2 * j))), PB))
    for j in range(1, b + 1):
        forms.append(ds.psub(ds.pmul(
            sin_of(2, 2), sin_of(2 * a, 2 * (b + 1 - 2 * j))), PB))

    # every closed form must be one of the engine's distinct diffs, the
    # match must be a BIJECTION, and the count must be 2a+2b+2
    unmatched = 0
    hit = set()
    for f in forms:
        found = [i for i, d in enumerate(diffs)
                 if ds.pclean(ds.psub(f, d)) == {}]
        if len(found) != 1:
            unmatched += 1
        else:
            hit.add(found[0])
    count_ok = (len(diffs) == 2 * a + 2 * b + 2
                and len(forms) == len(diffs) and len(hit) == len(diffs))
    # pivot antisymmetry p(B) + p(B1) = 2m  <=>  engine diff -PB present
    anti_ok = any(ds.pclean(ds.padd(f, PB)) == {} for f in diffs)
    return unmatched == 0 and count_ok and anti_ok


def cmd_selftest(_args):
    sys.path.insert(0, __import__("os").path.dirname(
        __import__("os").path.abspath(__file__)))
    ok = True
    members = [(2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2),
               (4, 3), (4, 4), (5, 2), (5, 5), (6, 1), (6, 3), (6, 4),
               (7, 1), (7, 7), (8, 2), (8, 8), (9, 4), (10, 9), (11, 3),
               (12, 5), (12, 12), (13, 2), (15, 7)]
    for (a, b) in members:
        r = _ring_check(a, b)
        ok &= r
        print(f"[ring] W({a},{b}): closed forms exact = {r}")

    # float cross-check: alive() vs the independent corridor width of 002
    from skeptic_family import width as sk_width, family_word
    import random
    rng = random.Random(20260816)
    agree = tried = 0
    for _ in range(400):
        a = rng.randint(2, 9)
        b = rng.randint(1, a)
        ac, bc = 90.0 / a, 90.0 * (a - 1) / (a * (b + 1))
        t = rng.uniform(-0.6, 0.6)
        al, be = ac - t, bc + 2 * t + rng.uniform(-0.3, 0.3)
        if not (0 < al and 0 < be and al + be < 90):
            continue
        ta = math.tan(al * D2R)
        tb = math.tan(be * D2R)
        x, y = tb / (ta + tb), ta * tb / (ta + tb)
        w = sk_width(family_word(a, b), x, y)
        mine = min_signed_margin(a, b, al, be)
        if abs(w) < 1e-9 or abs(mine) < 1e-9:
            continue  # too close to a boundary for float agreement
        tried += 1
        agree += ((w > 0) == (mine > 0))
    print(f"[corridor] alive() vs skeptic_family.width: "
          f"{agree}/{tried} agree")
    ok &= agree == tried
    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def cmd_segment(args):
    """Margin sign stress along the universal segment for all a <= amax."""
    worst = []
    for a in range(2, args.amax + 1):
        for b in range(1, a + 1):
            if (a, b) == (1, 1):
                continue
            ac, bc = 90.0 / a, 90.0 * (a - 1) / (a * (b + 1))
            t0 = min(45.0 / a, 90.0 / (2 * (b + 1)), 0.25)
            bad = None
            mn = math.inf
            for i in range(1, args.steps + 1):
                t = t0 * i / args.steps
                mm = min_signed_margin(a, b, ac - t, bc + 2 * t)
                mn = min(mn, mm)
                if mm <= 0 and bad is None:
                    bad = t
            worst.append((mn, a, b, bad))
    worst.sort()
    for mn, a, b, bad in worst[:12]:
        print(f"W({a},{b}): min signed margin {mn:+.3e}"
              + (f"  FIRST NONALIVE t={bad:.4f}" if bad else ""))
    n_bad = sum(1 for mn, _, _, bad in worst if bad is not None)
    print(f"{len(worst)} members, {n_bad} with a non-alive segment point")


def cmd_alive(args):
    for lbl, v, s in margins(args.a, args.b, args.alpha, args.beta):
        flag = "ok" if v * s > 0 else "VIOLATED"
        print(f"{lbl:6s} {v:+.9f} need {'+' if s > 0 else '-'}  {flag}")
    print("alive:", alive(args.a, args.b, args.alpha, args.beta))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)
    c = sub.add_parser("segment")
    c.add_argument("--amax", type=int, default=40)
    c.add_argument("--steps", type=int, default=200)
    c.set_defaults(fn=cmd_segment)
    c = sub.add_parser("alive")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--alpha", type=float, required=True)
    c.add_argument("--beta", type=float, required=True)
    c.set_defaults(fn=cmd_alive)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
