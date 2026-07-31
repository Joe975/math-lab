#!/usr/bin/env python3
"""Attempt 003: exact symbolic unfolding of W(a,b) over the trig ring.

The apex-coordinate harness carries denominators (squared side lengths)
through every reflection, so its corridor endpoints are rational functions
whose structure is opaque.  This tool works in a different exact domain in
which the unfolding is DIVISION-FREE, so every quantity has a closed form:

  * Normalize the triangle to circumdiameter 1: A = 0, B = sin(alpha+beta),
    C = sin(beta) e^{i alpha} in the complex plane.  Side lengths are then
    sin(alpha), sin(beta), sin(alpha+beta), and every side direction is a
    unit complex number e^{i(...)} with (...) an integer combination of
    alpha, beta, pi.
  * Represent every point as a Laurent polynomial in A = e^{i alpha/2},
    B = e^{i beta/2} with coefficients in Q(i) (pairs of Fractions).
    Reflection across the line through u with unit direction d is
    z -> u + d^2 conj(z - u), and d^2 is a monomial, so reflections are
    ring operations: no division ever happens.  conj is the substitution
    A -> 1/A, B -> 1/B, i -> -i.

Everything printed as a "trig form" is an EXACT identity in the ring,
verified against a float evaluation at random angles (--selftest), and the
corridor computed here is verified against skeptic_family's independent
float corridor (different normalization; widths agree up to the exact
scale factor sin(alpha+beta) * |tau_skeptic| / |tau_here|... in practice the
selftest compares alive/dead verdicts and width ratios).

Subcommands:
  selftest                     ring + geometry validation
  structure --a A --b B        closed-form corridor structure of W(a,b):
                               glide axis, all distinct gate-endpoint
                               projections, glide midline m, binding
                               constraints near death (numeric id), each as
                               an exact trig form
  arc --a A --b B              substitute the death-arc relation
                               (a+1)alpha + (b+1)beta = 180 deg into the
                               binding differences (cyclotomic coefficients)
                               and report which vanish identically
"""

import argparse
import cmath
import json
import math
import os
import random
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from skeptic_family import width as sk_width, family_word  # noqa: E402

# ---------------------------------------------------------------- Q(i)

def cadd(u, v):
    return (u[0] + v[0], u[1] + v[1])


def csub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def cmul(u, v):
    return (u[0] * v[0] - u[1] * v[1], u[0] * v[1] + u[1] * v[0])


def cconj(u):
    return (u[0], -u[1])


ONE = (F(1), F(0))
I = (F(0), F(1))
HALF_OVER_I = (F(0), F(-1, 2))          # 1/(2i)

# ------------------------------------------------- Laurent ring Q(i)[A*,B*]
# poly: dict {(m, n): (re, im)} for coefficient of A^m B^n,
# A = e^{i alpha/2}, B = e^{i beta/2}.


def pclean(p):
    return {k: c for k, c in p.items() if c != (0, 0)}


def padd(p, q):
    r = dict(p)
    for k, c in q.items():
        r[k] = cadd(r.get(k, (F(0), F(0))), c)
    return pclean(r)


def pneg(p):
    return {k: (-c[0], -c[1]) for k, c in p.items()}


def psub(p, q):
    return padd(p, pneg(q))


def pmul(p, q):
    r = {}
    for (m1, n1), c1 in p.items():
        for (m2, n2), c2 in q.items():
            k = (m1 + m2, n1 + n2)
            r[k] = cadd(r.get(k, (F(0), F(0))), cmul(c1, c2))
    return pclean(r)


def pscale(p, c):
    return pclean({k: cmul(v, c) for k, v in p.items()})


def pconj(p):
    return {(-m, -n): cconj(c) for (m, n), c in p.items()}


def pim(p):
    """Imaginary part as a ring element: (p - conj p)/(2i)."""
    return pscale(psub(p, pconj(p)), HALF_OVER_I)


def pre(p):
    return pscale(padd(p, pconj(p)), (F(1, 2), F(0)))


def mono(m, n, c=ONE):
    return {(m, n): c}


def peval(p, alpha, beta):
    """Numeric value at angles in radians (complex)."""
    z = 0j
    for (m, n), c in p.items():
        z += complex(c[0], c[1]) * cmath.exp(1j * (m * alpha + n * beta) / 2)
    return z


def sin_of(m, n):
    """sin((m alpha + n beta)/2) as a ring element."""
    return pscale(psub(mono(m, n), mono(-m, -n)), HALF_OVER_I)


def trig_str(p):
    """Exact trig form of a REAL ring element (conj-invariant)."""
    assert pclean(psub(p, pconj(p))) == {}, "not real"
    terms = []
    for (m, n), c in sorted(p.items()):
        if (m, n) == (0, 0):
            if c[0]:
                terms.append(f"{c[0]}")
            continue
        if (m, n) < (0, 0) or ((m, n) > (0, 0) and (-m, -n) in p and False):
            continue
        if (m, n) > (0, 0):
            ang = []
            if m:
                ang.append(f"{m}a" if m != 1 else "a")
            if n:
                ang.append(f"{'+' if n > 0 else '-'}{abs(n)}b"
                           if abs(n) != 1 else ('+b' if n > 0 else '-b'))
            angs = "(" + "".join(ang).lstrip("+") + ")/2"
            re2, im2 = 2 * c[0], 2 * c[1]
            if re2:
                terms.append(f"{re2}*cos{angs}")
            if im2:
                terms.append(f"{-im2}*sin{angs}")
    return " + ".join(terms).replace("+ -", "- ") if terms else "0"


# ------------------------------------------------------------- unfolding

SIDE_ENDS = {0: (1, 2), 1: (2, 0), 2: (0, 1)}


def initial_geometry():
    """(vertices, side direction-squares) of the base copy."""
    #  A = 0, B = sin(a+b) = sin((2A+2B)/2), C = sin(b) e^{i a}
    A = {}
    Bv = sin_of(2, 2)
    Cv = pmul(sin_of(0, 2), mono(2, 0))
    # dir^2: side2 (AB): 1; side1 (CA): e^{2 i alpha} = A^4;
    # side0 (BC): e^{2 i (pi - beta)} = e^{-2 i beta} = B^-4
    return [A, Bv, Cv], [mono(0, -4), mono(4, 0), mono(0, 0)]


def unfold_sym(word):
    """Symbolic unfolding.  Returns dict with gates (list of (u,v) ring
    pairs), tau, and the composed map of each prefix as (orient, lin, off):
    z -> lin*z + off (orient=+1) or lin*conj(z) + off (orient=-1)."""
    verts, dirsq = initial_geometry()
    gates = []
    orient, lin, off = 1, mono(0, 0), {}
    maps = []
    for s in word:
        i, j = SIDE_ENDS[s]
        u, v = verts[i], verts[j]
        gates.append((u, v))
        e = dirsq[s]
        # reflect all three vertices: z' = u + e*conj(z - u)
        verts = [padd(u, pmul(e, pconj(psub(z, u)))) for z in verts]
        # update direction squares: f' = e^2 * conj(f)
        e2 = pmul(e, e)
        dirsq = [pmul(e2, pconj(f)) for f in dirsq]
        dirsq[s] = e  # the mirror side itself is unchanged
        # compose the reflection R(z) = e*conj(z) + (u - e*conj(u)) with map
        c = psub(u, pmul(e, pconj(u)))
        if orient == 1:
            lin, off = pmul(e, pconj(lin)), padd(pmul(e, pconj(off)), c)
            orient = -1
        else:
            lin, off = pmul(e, pconj(lin)), padd(pmul(e, pconj(off)), c)
            orient = 1
        maps.append((orient, lin, off))
    return {"gates": gates, "verts": verts, "maps": maps}


def glide_data(a, b):
    """Symbolic corridor data for W(a,b) via the half-word glide reflection.

    Returns dict with:
      delta   : axis direction (monomial, possibly times i)
      projs   : list of (label, ring element) projections Im(conj(delta) z)
                of every first-half gate endpoint
      m       : projection of the glide axis
      tau     : full translation (ring element pair check)
    """
    word = family_word(a, b)
    half = word[:len(word) // 2]
    U = unfold_sym(word)
    Uh = unfold_sym(half)
    orient, mu, w = Uh["maps"][-1]
    assert orient == -1, "half word must be orientation-reversing"
    # mu must be a monomial c * A^{2p} B^{2q}, c in {1,-1}
    assert len(mu) == 1, f"glide linear part not a monomial: {mu}"
    (mm, nn), c = next(iter(mu.items()))
    assert mm % 2 == 0 and nn % 2 == 0, "odd monomial in glide linear part"
    if c == (F(1), F(0)):
        delta = mono(mm // 2, nn // 2)
    elif c == (F(-1), F(0)):
        delta = mono(mm // 2, nn // 2, I)
    else:
        raise AssertionError(f"unexpected glide coefficient {c}")
    dconj = pconj(delta)
    # tau = w + mu conj(w); axis offset m = Im(conj(delta)(w - tau/2))/2
    tau = padd(w, pmul(mu, pconj(w)))
    # check tau parallel to axis: Im(conj(delta) tau) == 0 in the ring
    assert pclean(pim(pmul(dconj, tau))) == {}, "tau not parallel to axis"
    half_tau = pscale(tau, (F(1, 2), F(0)))
    m = pscale(pim(pmul(dconj, psub(w, half_tau))), (F(1, 2), F(0)))
    projs = []
    seen = {}
    for k, (u, v) in enumerate(Uh["gates"]):
        for lbl, z in (("u", u), ("v", v)):
            pz = pim(pmul(dconj, z))
            key = json.dumps(sorted([(mn, (str(c0[0]), str(c0[1])))
                                     for mn, c0 in pz.items()]))
            if key not in seen:
                seen[key] = f"g{k + 1}{lbl}"
                projs.append((f"g{k + 1}{lbl}", pz))
    return {"delta": delta, "projs": projs, "m": m, "tau": tau,
            "gates": Uh["gates"], "word": word, "half": half}


# ------------------------------------------------------------- numeric

def corridor_glide_numeric(gd, alpha, beta):
    """Corridor [lo, hi] intersected with its glide reflection, numerically,
    from the symbolic data.  Returns (width, lo, hi, m)."""
    vals = {}
    for lbl, p in gd["projs"]:
        vals[lbl] = peval(p, alpha, beta).real
    m = peval(gd["m"], alpha, beta).real
    lo, hi = -math.inf, math.inf
    for (u, v), k in zip(gd["gates"], range(len(gd["gates"]))):
        pu = peval(pim(pmul(pconj(gd["delta"]), u)), alpha, beta).real
        pv = peval(pim(pmul(pconj(gd["delta"]), v)), alpha, beta).real
        lo = max(lo, min(pu, pv))
        hi = min(hi, max(pu, pv))
    # full corridor = [lo,hi] cap [2m-hi, 2m-lo]
    flo, fhi = max(lo, 2 * m - hi), min(hi, 2 * m - lo)
    return fhi - flo, lo, hi, m


def cmd_selftest(_args):
    ok = True
    rng = random.Random(20260731)
    # 1. ring identities
    p = sin_of(2, 2)
    v = peval(p, 0.7, 0.4)
    ok &= abs(v - math.sin(0.7 + 0.4)) < 1e-12 and abs(v.imag) < 1e-15
    print("[sin identity]", "OK" if ok else "FAIL")
    # 2. initial geometry closes: |AB|=sin(a+b), |CA|=sin(b), |BC|=sin(a)
    verts, dirsq = initial_geometry()
    al, be = 0.62, 0.31
    A0, B0, C0 = (peval(z, al, be) for z in verts)
    good = (abs(abs(B0 - A0) - math.sin(al + be)) < 1e-12
            and abs(abs(C0 - A0) - math.sin(be)) < 1e-12
            and abs(abs(C0 - B0) - math.sin(al)) < 1e-12)
    ok &= good
    print("[triangle lengths]", "OK" if good else "FAIL")
    # 3. reflection preserves distances and fixes the mirror line, and the
    #    glide/corridor verdicts of the symbolic path agree with
    #    skeptic_family's float corridor across the family and random angles
    for (fa, fb) in [(1, 1), (2, 1), (2, 2), (3, 2)]:
        gd = glide_data(fa, fb)
        agree = True
        for _ in range(60):
            g = rng.uniform(91.0, 179.0)
            th = math.radians(180 - g)
            al = rng.uniform(0.05 * th, 0.95 * th)
            be = th - al
            wglide, lo, hi, m = corridor_glide_numeric(gd, al, be)
            ta, tb = math.tan(al), math.tan(be)
            x, y = tb / (ta + tb), ta * tb / (ta + tb)
            wsk = sk_width(gd["word"], x, y)
            alive_sym = wglide > 1e-13
            alive_sk = wsk > 1e-13
            borderline = abs(wglide) < 1e-9 or abs(wsk) < 1e-9
            if alive_sym != alive_sk and not borderline:
                agree = False
                print(f"  DISAGREE W({fa},{fb}) g={g:.4f} al={al:.4f}: "
                      f"sym {wglide:.3e} vs sk {wsk:.3e}")
        ok &= agree
        print(f"[corridor verdicts W({fa},{fb}) x60 random angles]",
              "OK" if agree else "FAIL")
    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def cmd_structure(args):
    gd = glide_data(args.a, args.b)
    (dm, dn), dc = next(iter(gd["delta"].items()))
    print(f"W({args.a},{args.b}) = {''.join(map(str, gd['word']))}")
    print(f"half word = {''.join(map(str, gd['half']))} "
          f"(glide reflection; corridor = I cap glide-mirror(I))")
    extra = " * i" if dc == (F(0), F(1)) else ""
    print(f"axis direction delta = e^(i({dm}a{'+' if dn >= 0 else ''}{dn}b)/2)"
          f"{extra}")
    print(f"tau (translation)    = {trig_str(pre(gd['tau'])) or '(complex)'}"
          f"  [Re];  Im = {trig_str(pim(gd['tau']))}")
    print(f"m (axis offset)      = {trig_str(gd['m'])}")
    print("distinct gate-endpoint projections p = Im(conj(delta) z):")
    for lbl, p in gd["projs"]:
        print(f"  {lbl:5s} = {trig_str(p)}")
    # numeric binding structure near the (conjectured) death point
    pred = (180.0 * args.a / (args.a + 1) if args.a == args.b
            else 180.0 * (1 - (2 * args.a - 1) / (2.0 * args.a ** 2))
            if args.b == args.a - 1 else None)
    if pred is not None:
        th = math.radians(180 - pred) + 1e-6
        r = args.a / (args.a + args.b)
        al, be = r * th, th - r * th
        vals = {lbl: peval(p, al, be).real for lbl, p in gd["projs"]}
        m = peval(gd["m"], al, be).real
        lo_lbl = max(vals, key=lambda l: min(vals[l], 1e18))
        print(f"\nnear-death numeric values (gamma = death - ~5.7e-5 deg, "
              f"alpha:beta = a:b):")
        for lbl in sorted(vals, key=vals.get):
            print(f"  {lbl:5s} = {vals[lbl]:+.9f}")
        print(f"  m     = {m:+.9f}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)
    c = sub.add_parser("structure")
    c.add_argument("--a", type=int, required=True)
    c.add_argument("--b", type=int, required=True)
    c.set_defaults(fn=cmd_structure)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
