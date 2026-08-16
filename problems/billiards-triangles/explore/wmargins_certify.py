#!/usr/bin/env python3
"""Certified (exact interval) alive-segments for the W(a,b) window theorem.

For each member this certifies, in exact rational interval arithmetic (the
certified-trig layer of plaw_suffice, attempt 005), that BOTH canonical
segments are alive:

  DEATH segment  (alpha, beta) = (90/a - t, 90(a-1)/(a(b+1)) + 2t),
                 t in (0, 22.5/(a(b+1))]   (= half the distance to the
                 PB = 0 wall at t_PB = 45/(a(b+1)); theta = theta_d + t)
  BIRTH segment  (alpha, beta) = (90/a - t, 90/(b+1) - t),
                 t in (0, 22.5/(a(b+1))]   (theta = theta_birth - 2t)

Together the two segments realize EVERY theta in (theta_d, theta_mid] u
[theta_mid, theta_birth) = (theta_d, theta_birth), theta_mid = theta_d +
45/(a(b+1)): the window of W(a,b) is fully swept (the sufficiency half of
the window law; necessity is 005's theorem on the death side and open on
the birth side, see the 015 record).

Margins are the closed forms of wmargins.py (ring-verified there):
  N2, PB           positive by rational range checks on both segments
  A-fan / B-fan    sinusoid margins, certified by adaptive interval
                   bisection (they are nonzero at t = 0 on each segment,
                   EXCEPT the three corner-degenerate ones below)
  N1 (death)       vanishes at t = 0: N1' > 0 on [0, t1] + bisection on
                   [t1, T]  (exactly 005's treatment, same code path)
  AB_a (birth)     equals cos(at) sin((b+1)t) sin(theta) > 0 by the sine
                   addition formula: rational range check only
  tau != 0         Re(conj(delta) tau) certified sign via SegPoly

Usage:
  wmargins_certify.py member --a A --b B
  wmargins_certify.py sweep --amax 12 [--out FILE]
"""

import argparse
import json
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "..", "..",
                                "harness", "billiards-triangles"))
import deathlaw_symbolic as ds  # noqa: E402
from plaw_suffice import (sin_deg, cos_deg, certify_sign, SegPoly,  # noqa: E402
                          n1_iv, n1p_iv, _trig_cache)


def lin(c0, c1, tlo, thi):
    a1, a2 = c0 + c1 * tlo, c0 + c1 * thi
    return (a1, a2) if a1 <= a2 else (a2, a1)


def death_margins(a, b):
    """[(label, evalf(tlo,thi) -> Iv, need_sign)] for the death segment,
    excluding N1/N2/PB which are handled separately."""
    A, bc = F(90, a), F(90 * (a - 1), a * (b + 1))
    out = []

    def n2(tlo, thi):
        return (sin_deg(*lin(F(0), F(a), tlo, thi))
                * cos_deg(*lin(A, F(-2 * (b + 1)), tlo, thi))
                * sin_deg(*lin(A + bc, F(1), tlo, thi)))

    for s in range(1, a):
        c0, c1 = 2 * s * A, F(a + 2 * b + 1 - 2 * s)
        out.append((f"AC{s}", (lambda s=s, c0=c0, c1=c1: lambda tlo, thi:
                    n2(tlo, thi) - sin_deg(*lin(bc, F(2), tlo, thi))
                    * sin_deg(*lin(c0, c1, tlo, thi)))(), -1))
    for s in range(1, a):
        c0, c1 = (2 * s - 1) * A, F(a + 2 * b + 2 - 2 * s)
        out.append((f"AB{s}", (lambda s=s, c0=c0, c1=c1: lambda tlo, thi:
                    n2(tlo, thi) - sin_deg(*lin(A + bc, F(1), tlo, thi))
                    * sin_deg(*lin(c0, c1, tlo, thi)))(), -1))

    def pb(tlo, thi):
        return (sin_deg(*lin(A, F(-2 * (b + 1)), tlo, thi))
                * cos_deg(*lin(F(0), F(a), tlo, thi))
                * sin_deg(*lin(A + bc, F(1), tlo, thi)))

    for j in range(1, b):
        c0 = 90 + (b - 2 * j) * bc
        c1 = F(-a + 2 * (b - 2 * j))
        out.append((f"BC{j}", (lambda j=j, c0=c0, c1=c1: lambda tlo, thi:
                    sin_deg(*lin(A, F(-1), tlo, thi))
                    * sin_deg(*lin(c0, c1, tlo, thi)) - pb(tlo, thi))(), +1))
    for j in range(1, b + 1):
        c0 = 90 + (b + 1 - 2 * j) * bc
        c1 = F(-a + 2 * (b + 1 - 2 * j))
        out.append((f"BA{j}", (lambda j=j, c0=c0, c1=c1: lambda tlo, thi:
                    sin_deg(*lin(A + bc, F(1), tlo, thi))
                    * sin_deg(*lin(c0, c1, tlo, thi)) - pb(tlo, thi))(), +1))
    return out


def birth_margins(a, b):
    """Same for the birth segment (alpha, beta) = (A - t, Bh - t);
    excludes N2, PB, AB_a (rational range checks)."""
    A, Bh = F(90, a), F(90, b + 1)
    th0 = A + Bh
    out = []

    def n2(tlo, thi):
        return (sin_deg(*lin(F(0), F(a), tlo, thi))
                * cos_deg(*lin(F(0), F(b + 1), tlo, thi))
                * sin_deg(*lin(th0, F(-2), tlo, thi)))

    for s in range(0, a):
        c0, c1 = (2 * s + 1) * A, F(a - b - 2 - 2 * s)
        out.append((f"AC{s}", (lambda s=s, c0=c0, c1=c1: lambda tlo, thi:
                    n2(tlo, thi) - sin_deg(*lin(Bh, F(-1), tlo, thi))
                    * sin_deg(*lin(c0, c1, tlo, thi)))(), -1))
    for s in range(1, a):
        c0, c1 = 2 * s * A, F(a - b - 1 - 2 * s)
        out.append((f"AB{s}", (lambda s=s, c0=c0, c1=c1: lambda tlo, thi:
                    n2(tlo, thi) - sin_deg(*lin(th0, F(-2), tlo, thi))
                    * sin_deg(*lin(c0, c1, tlo, thi)))(), -1))

    def pb(tlo, thi):
        return (sin_deg(*lin(F(0), F(b + 1), tlo, thi))
                * cos_deg(*lin(F(0), F(a), tlo, thi))
                * sin_deg(*lin(th0, F(-2), tlo, thi)))

    for j in range(0, b):
        c0, c1 = 90 + (b - 2 * j) * Bh, F(-a - (b - 2 * j))
        out.append((f"BC{j}", (lambda j=j, c0=c0, c1=c1: lambda tlo, thi:
                    sin_deg(*lin(A, F(-1), tlo, thi))
                    * sin_deg(*lin(c0, c1, tlo, thi)) - pb(tlo, thi))(), +1))
    for j in range(1, b + 1):
        c0, c1 = 90 + (b + 1 - 2 * j) * Bh, F(-a - (b + 1 - 2 * j))
        out.append((f"BA{j}", (lambda j=j, c0=c0, c1=c1: lambda tlo, thi:
                    sin_deg(*lin(th0, F(-2), tlo, thi))
                    * sin_deg(*lin(c0, c1, tlo, thi)) - pb(tlo, thi))(), +1))
    return out


def rational_checks(a, b, T):
    """The margins provable by rational range comparisons alone.
    Death: N2 > 0, PB > 0.  Birth: N2 > 0, PB > 0, AB_a > 0.
    Returns list of failed checks (empty = all pass)."""
    fails = []
    A, bc, Bh = F(90, a), F(90 * (a - 1), a * (b + 1)), F(90, b + 1)
    th_d = A + bc
    # death: at in (0, aT] within (0,90); (b+1)beta in (90-A, 90) strictly
    # below 90 iff T < 45/(a(b+1)); theta in (th_d, th_d+T] within (0,90)
    if not (a * T < 90):
        fails.append("death: aT < 90")
    if not (2 * (b + 1) * T < A):
        fails.append("death: (b+1)beta < 90 up to T")
    if not (th_d + T < 90):
        fails.append("death: theta < 90")
    # birth: at in (0, aT]; (b+1)t in (0, (b+1)T]; theta = A+Bh-2t in (0,90]
    # with theta < 180 always and sin > 0; A+Bh <= 90 (equality (2,1) only,
    # then theta < 90 for t > 0)
    if not (a * T < 90 and (b + 1) * T < 90):
        fails.append("birth: at,(b+1)t < 90")
    if not (A + Bh <= 90):
        fails.append("birth: theta_birth <= 90")
    return fails


def tau_polys(a, b):
    gd = ds.glide_data(a, b)
    return ds.pre(ds.pmul(ds.pconj(gd["delta"]), gd["tau"]))


def certify_member(a, b, max_depth=40):
    _trig_cache.clear()
    T = F(45, 2 * a * (b + 1))            # 22.5/(a(b+1))
    rep = {"a": a, "b": b, "T": str(T)}
    fails = rational_checks(a, b, T)
    if fails:
        rep["ok"] = False
        rep["fail"] = fails
        return rep

    leaves = {}
    # death N1: derivative near 0, bisection beyond (005's treatment, rho=2)
    ac, bc = F(90, a), F(90 * (a - 1), a * (b + 1))
    wh = F(2 * b + 1)
    t1 = T / 8
    for _ in range(24):
        n = certify_sign(lambda lo, hi: n1p_iv(a, b, ac, bc, 2, wh, lo, hi),
                         +1, 0, t1, max_depth)
        if n is not None:
            break
        t1 /= 2
    else:
        rep["ok"] = False
        rep["fail"] = "death N1' near 0"
        return rep
    leaves["death_N1_deriv"] = n
    n = certify_sign(lambda lo, hi: n1_iv(a, b, ac, bc, 2, wh, lo, hi),
                     +1, t1, T, max_depth)
    if n is None:
        rep["ok"] = False
        rep["fail"] = "death N1 bisection"
        return rep
    leaves["death_N1"] = n

    for seg, margin_list in (("death", death_margins(a, b)),
                             ("birth", birth_margins(a, b))):
        tot = 0
        for lbl, evalf, sign in margin_list:
            n = certify_sign(evalf, sign, 0, T, max_depth)
            if n is None:
                rep["ok"] = False
                rep["fail"] = f"{seg} {lbl}"
                return rep
            tot += n
        leaves[f"{seg}_fan_leaves"] = tot

    # tau nondegenerate on both segments
    taur = tau_polys(a, b)
    for seg, (c_a, c_b, rho) in (("death", (ac, bc, 2)),
                                 ("birth", (ac, F(90, b + 1), -1))):
        sp = SegPoly(taur, c_a, c_b, rho)
        mid = sp.eval_iv(T / 2, T / 2)
        s = 1 if mid.lo > 0 else (-1 if mid.hi < 0 else 0)
        n = certify_sign(sp.eval_iv, s, 0, T, max_depth) if s else None
        if n is None:
            rep["ok"] = False
            rep["fail"] = f"{seg} tau"
            return rep
        leaves[f"{seg}_tau"] = n

    rep["ok"] = True
    rep["leaves"] = leaves
    rep["theta_covered"] = (
        f"({90*(a+b)}/{a*(b+1)}, {90*(a+b)+45}/{a*(b+1)}] u "
        f"[{90*(a+b)+45}/{a*(b+1)}, {90*(a+b+1)}/{a*(b+1)})")
    return rep


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("member")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c = sub.add_parser("sweep")
    c.add_argument("--amax", type=int, default=12)
    c.add_argument("--out")
    args = ap.parse_args()
    if args.cmd == "member":
        print(json.dumps(certify_member(args.a, args.b), indent=1))
        return
    results = []
    for a in range(2, args.amax + 1):
        for b in range(1, a + 1):
            r = certify_member(a, b)
            results.append(r)
            print(f"W({a},{b}): {'CERTIFIED both segments' if r['ok'] else 'FAIL ' + str(r.get('fail'))}")
            if args.out:
                json.dump(results, open(args.out, "w"), indent=1)
    ok = all(r["ok"] for r in results)
    print(f"{len(results)} members:",
          "ALL CERTIFIED" if ok else "SOME FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
