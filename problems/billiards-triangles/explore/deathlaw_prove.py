#!/usr/bin/env python3
"""Attempt 003: per-member certified proof obligations for the W(a,b)
death-angle theorem.

The paper argument (recorded in attempts/003-death-angle-laws.md) proves:

  THEOREM (necessity).  Let a >= b >= 1, a >= 2, a <= 2b+3.  If the
  unfolding corridor of W(a,b) = (0 (12)^a (02)^b)^2 has positive width at a
  triangle with angles alpha (at A), beta (at B), then
      alpha + beta > theta_d(a,b) = 90 (a+b) / (a (b+1))  degrees,
  i.e. the obtuse angle satisfies gamma < gamma_d(a,b) = 180 - theta_d.

modulo exactly three machine-checkable inputs per member (a,b), which this
tool certifies:

  [ID]  Exact ring identities (in Q(i)[e^{i alpha/2}, e^{i beta/2}]) for the
        glide-axis projections of the first two gate endpoints and the glide
        midline m of the half word u = 0(12)^a(02)^b:
          I1:  m - p(C)   = -[cos(a alpha) sin(alpha) sin(b beta)
                              + sin(a alpha) sin(beta) cos(alpha+(b+1)beta)]
          I2:  p(A1) - m  = cos(a alpha) sin((b+1) beta) sin(alpha+beta)
          I3:  p(B) - m   = cos((b+1) beta) sin(a alpha) sin(alpha+beta)
        plus the glide facts: the half-word isometry is orientation-
        reversing with monomial linear part mu = delta^2 (delta the axis
        direction) and translation tau parallel to the axis
        (Im(conj(delta) tau) = 0 in the ring).  These four facts make
        "corridor has positive width  =>  m lies strictly inside the
        projection intervals of gate 1 = [B, C] and gate 2 = [C, A1]"
        a theorem (gate n'+k is the glide image of gate k, and a glide
        reflection acts on axis-normal projections by p -> 2m - p).

  [LC]  Lemma C: with v0 = 90/(b+1) deg, for all v in (0, v0):
          H2(v) := sin(v) sin((90-v)/a) cos(v0)
                   - sin(v0) sin(b v/a) cos(v)  > 0,
        and H2(v0) = 0 identically (algebra: (90-v0)/a = b v0 / a).
        Certified here by: (i) an explicit positive lower bound on (0, d1]
        via sin x <= x <= sin x / cos x; (ii) adaptive interval bisection
        on [d1, v0 - d2]; (iii) a certified sign H2' < 0 on [v0 - d2, v0],
        which forces H2 > 0 there since H2(v0) = 0.  All enclosures are
        unfold.py's rational-interval sin/cos with the rational pi
        enclosure; single (non-iterated) evaluations, so plain intervals
        are the right tool (the affine-form lesson applies to iterated
        maps).

  [SANE]  a numeric spot-check that the projections evaluate to real
        numbers at a generic angle pair (guards the ring plumbing).

  Also checked exactly, context rather than load-bearing:
          D2:  p(A1) - p(C) = -sin(beta) sin((a-1)alpha - (b+1)beta)
          D1:  p(B)  - p(C) =  sin(alpha) sin(a alpha - b beta)
        (the collapsing-gate loci the float layer sees as the pinch).

Usage:
  deathlaw_prove.py member --a A --b B [--out FILE]
  deathlaw_prove.py all    [--out FILE]     (the members recorded in 003)
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

Iv = unfold.Iv


# ------------------------------------------------------------------ [ID]

def sin_of(m, n):
    return ds.sin_of(m, n)


def cos_of(m, n):
    return ds.pscale(ds.padd(ds.mono(m, n), ds.mono(-m, -n)),
                     (F(1, 2), F(0)))


def check_identities(a, b):
    """Exact verification of I1, I2, I3 and the glide facts.  Returns dict
    of booleans; glide facts are enforced by assertions in glide_data."""
    gd = ds.glide_data(a, b)          # asserts glide facts internally
    dconj = ds.pconj(gd["delta"])
    proj = lambda z: ds.pim(ds.pmul(dconj, z))
    g1u, g1v = gd["gates"][0]         # gate 1 = (B, C)
    g2u, g2v = gd["gates"][1]         # gate 2 = (C, A1)
    pB, pC, pP = proj(g1u), proj(g1v), proj(g2v)
    m = gd["m"]

    I1 = ds.pneg(ds.padd(
        ds.pmul(ds.pmul(cos_of(2 * a, 0), sin_of(2, 0)), sin_of(0, 2 * b)),
        ds.pmul(ds.pmul(sin_of(2 * a, 0), sin_of(0, 2)),
                cos_of(2, 2 * (b + 1)))))
    I2 = ds.pmul(ds.pmul(cos_of(2 * a, 0), sin_of(0, 2 * (b + 1))),
                 sin_of(2, 2))
    I3 = ds.pmul(ds.pmul(cos_of(0, 2 * (b + 1)), sin_of(2 * a, 0)),
                 sin_of(2, 2))

    ok1 = ds.pclean(ds.psub(ds.psub(m, pC), I1)) == {}
    ok2 = ds.pclean(ds.psub(ds.psub(pP, m), I2)) == {}
    ok3 = ds.pclean(ds.psub(ds.psub(pB, m), I3)) == {}
    # collapsing-gate identities (context, not load-bearing for the theorem):
    #   D2: p(A1) - p(C) = -sin(beta) sin((a-1)alpha - (b+1)beta)
    #   D1: p(B)  - p(C) =  sin(alpha) sin(a alpha - b beta)
    D2 = ds.pneg(ds.pmul(sin_of(0, 2), sin_of(2 * (a - 1), -2 * (b + 1))))
    D1 = ds.pmul(sin_of(2, 0), sin_of(2 * a, -2 * b))
    okd2 = ds.pclean(ds.psub(ds.psub(pP, pC), D2)) == {}
    okd1 = ds.pclean(ds.psub(ds.psub(pB, pC), D1)) == {}
    # numeric sanity: identities really describe the same functions floats see
    al, be = 0.37, 0.24
    sane = all(abs(ds.peval(p, al, be).imag) < 1e-9 for p in (pB, pC, pP, m))
    return {"I1": ok1, "I2": ok2, "I3": ok3, "D1": okd1, "D2": okd2,
            "glide_facts": True, "numeric_sane": sane}


# ------------------------------------------------------------------ [LC]

def sin_deg(x):
    """Certified enclosure of sin(x degrees), x an Iv of rationals."""
    rad = x * Iv(unfold.PI_LO, unfold.PI_HI) / 180
    return unfold.sin_iv(rad)


def cos_deg(x):
    rad = x * Iv(unfold.PI_LO, unfold.PI_HI) / 180
    return unfold.cos_iv(rad)


def h2_iv(v, a, b, v0):
    """Enclosure of H2 over the rational-degree interval v."""
    return (sin_deg(v) * sin_deg((Iv(90) - v) * F(1, a)) * cos_deg(Iv(v0))
            - sin_deg(Iv(v0)) * sin_deg(v * F(b, a)) * cos_deg(v))


def h2p_iv(v, a, b, v0):
    """Enclosure of dH2/dv (per degree; the pi/180 factor is positive and
    dropped, which does not change the sign)."""
    c0, s0 = cos_deg(Iv(v0)), sin_deg(Iv(v0))
    t1 = cos_deg(v) * sin_deg((Iv(90) - v) * F(1, a)) * c0
    t2 = sin_deg(v) * cos_deg((Iv(90) - v) * F(1, a)) * c0 * F(1, a)
    t3 = s0 * cos_deg(v * F(b, a)) * cos_deg(v) * F(b, a)
    t4 = s0 * sin_deg(v * F(b, a)) * sin_deg(v)
    return t1 - t2 - t3 + t4


def certify_lemma_c(a, b, max_depth=40):
    """Certify H2 > 0 on (0, v0) with v0 = 90/(b+1).  Returns dict."""
    v0 = F(90, b + 1)
    report = {"v0": str(v0)}

    # (iii) near v0: certify H2' < 0 on [v0 - d2, v0]; then for v in that
    # interval H2(v) = -int_v^{v0} H2' > 0 because H2(v0) = 0 exactly.
    d2 = v0 / 64
    while d2 > v0 / 2 ** 20:
        e = h2p_iv(Iv(v0 - d2, v0), a, b, v0)
        if e.hi < 0:
            break
        d2 /= 2
    else:
        return {"ok": False, "fail": "H2' sign near v0"}
    report["d2"] = str(d2)

    # (i) near 0: for v in (0, d1]:
    #   sin(b v/a) <= sin(v) * (b/a) / cos(v)   (0 < bv/a <= v < 90)
    # so H2 >= sin(v) [ sin((90-d1)/a) cos(v0) - sin(v0) (b/a) / cos(d1) ]
    # and the bracket is certified positive at a concrete d1.
    d1 = v0 / 8
    while d1 > v0 / 2 ** 20:
        br = (sin_deg((Iv(90) - Iv(d1)) * F(1, a)) * cos_deg(Iv(v0))
              - sin_deg(Iv(v0)) * Iv(F(b, a)) / cos_deg(Iv(d1)))
        if br.lo > 0:
            break
        d1 /= 2
    else:
        return {"ok": False, "fail": "near-0 bracket"}
    report["d1"] = str(d1)

    # (ii) adaptive bisection on [d1, v0 - d2]
    stack = [(d1, v0 - d2)]
    n_leaves = 0
    while stack:
        lo, hi = stack.pop()
        e = h2_iv(Iv(lo, hi), a, b, v0)
        if e.lo > 0:
            n_leaves += 1
            continue
        if (hi - lo) < v0 / 2 ** max_depth:
            return {"ok": False, "fail": f"bisection stuck at [{lo},{hi}]"}
        mid = (lo + hi) / 2
        stack.append((lo, mid))
        stack.append((mid, hi))
    report["ok"] = True
    report["bisection_leaves"] = n_leaves
    return report


# ------------------------------------------------------------------ main

def gamma_d(a, b):
    return 180 - F(90 * (a + b), a * (b + 1))


def run_member(a, b):
    res = {"a": a, "b": b, "len": 4 * a + 4 * b + 2,
           "theta_d": str(F(90 * (a + b), a * (b + 1))),
           "gamma_d": str(gamma_d(a, b)),
           "scope_ok": (a >= b >= 1 and a >= 2 and a <= 2 * b + 3)}
    ids = check_identities(a, b)
    res["identities"] = ids
    lc = certify_lemma_c(a, b)
    res["lemma_c"] = lc
    res["ok"] = (all(ids.values()) and lc.get("ok", False)
                 and res["scope_ok"])
    return res


MEMBERS = [(2, 1), (2, 2), (3, 2), (3, 3), (4, 2), (4, 3), (4, 4),
           (5, 3), (5, 4), (5, 5), (6, 6), (7, 7), (8, 8),
           (10, 9), (12, 12)]


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("member")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.add_argument("--out")
    c = sub.add_parser("all")
    c.add_argument("--out")
    args = ap.parse_args()

    if args.cmd == "member":
        results = [run_member(args.a, args.b)]
    else:
        results = []
        for (a, b) in MEMBERS:
            r = run_member(a, b)
            results.append(r)
            print(f"W({a},{b}): gamma_d={r['gamma_d']} "
                  f"identities={all(r['identities'].values())} "
                  f"lemmaC={r['lemma_c'].get('ok')} "
                  f"leaves={r['lemma_c'].get('bisection_leaves')} "
                  f"=> {'PROVEN' if r['ok'] else 'FAILED'}")
    if args.out:
        json.dump(results, open(args.out, "w"), indent=1)
    if args.cmd == "member":
        print(json.dumps(results[0], indent=1))
    ok = all(r["ok"] for r in results)
    print("ALL OBLIGATIONS CERTIFIED" if ok else "SOME OBLIGATIONS FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
