#!/usr/bin/env python3
"""Attempt 006: symbolic binding-factor extraction for arbitrary doubled
words (the design inversion of 003's death-law machinery).

Built on deathlaw_symbolic's division-free Laurent ring
Q(i)[e^{i alpha/2}, e^{i beta/2}]. For any doubled odd word w = u^2 this
tool computes the glide data (axis direction delta, offset m, half-gate
endpoint projections p(z) = Im(conj(delta) z)) exactly, then FACTORS the
binding differences p(z) - m and p(z) - p(z') into elementary factors

    s(p,q) = sin((p alpha + q beta)/2),   c(p,q) = cos((p alpha + q beta)/2)

by numeric candidate detection followed by EXACT division in the ring
(division by the binomial M -+ M^{-1}, M = A^p B^q, with a termination
cap; a reported factorization is exact: quotient times factors times a
monomial unit reproduces the polynomial, re-verified by multiplication).

The death/birth angle of a family then reads off the factor list: an alive
window edge sits where a binding factor vanishes, e.g. c(2a,0) = cos(a
alpha) vanishing at alpha = 90/a is W(a,b)'s death edge (003's I2).

Subcommands:
  selftest                  re-derive 003's I2/I3 factorizations for W(3,3)
  bind --word W --gamma G   binding gates at the float argmax on arc G,
                            and the exact factorization of each binding
                            difference
  birth --a A --b B         factor the differences binding at the BIRTH
                            corner of W(a,b) (003 Lead 5)

Stdlib only, deterministic.
"""

import argparse
import json
import math
import os
import sys
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import deathlaw_symbolic as DS  # noqa: E402
from deathlaw_symbolic import (padd, psub, pmul, pscale, pconj, pim,  # noqa
                               pclean, mono, peval, ONE, I, HALF_OVER_I)
from design_neighbours import apex_from_angles, arc_alphas  # noqa: E402


def glide_data_word(word):
    """deathlaw_symbolic.glide_data generalized to any doubled odd word."""
    n = len(word)
    half = word[:n // 2]
    assert n % 2 == 0 and word == half * 2 and len(half) % 2 == 1
    Uh = DS.unfold_sym(half)
    orient, mu, w = Uh["maps"][-1]
    assert orient == -1
    assert len(mu) == 1, f"glide linear part not a monomial: {mu}"
    (mm, nn), c = next(iter(mu.items()))
    assert mm % 2 == 0 and nn % 2 == 0
    if c == (F(1), F(0)):
        delta = mono(mm // 2, nn // 2)
    elif c == (F(-1), F(0)):
        delta = mono(mm // 2, nn // 2, I)
    else:
        raise AssertionError(f"unexpected glide coefficient {c}")
    dconj = pconj(delta)
    tau = padd(w, pmul(mu, pconj(w)))
    assert pclean(pim(pmul(dconj, tau))) == {}, "tau not parallel to axis"
    m = pscale(pim(pmul(dconj, psub(w, pscale(tau, (F(1, 2), F(0)))))),
               (F(1, 2), F(0)))
    endpoints = []   # (gate index, endpoint label, projection ring element)
    for k, (u, v) in enumerate(Uh["gates"]):
        endpoints.append((k + 1, "u", pim(pmul(dconj, u))))
        endpoints.append((k + 1, "v", pim(pmul(dconj, v))))
    return {"delta": delta, "m": m, "tau": tau, "endpoints": endpoints,
            "gates": Uh["gates"], "half": half, "word": word}


# ----------------------------------------------------- exact factorization

def divide_binomial(P, p, q, sign):
    """Exact division of Laurent poly P by D = A^p B^q + sign * A^-p B^-q.
    Returns quotient Q with P = Q*D, or None if not divisible (cap hit or
    nonzero remainder)."""
    if not P:
        return None
    D = padd(mono(p, q), mono(-p, -q, (F(sign), F(0))))
    Q = {}
    R = dict(P)
    cap = 4 * (len(P) + 8) * (abs(p) + abs(q) + 2)
    for _ in range(cap):
        if not R:
            return Q
        # lead term: max of phi = p*m + q*n, tie-break by exponent tuple
        k = max(R, key=lambda mn: (p * mn[0] + q * mn[1], mn))
        c = R[k]
        qk = (k[0] - p, k[1] - q)
        # if the lead cannot be reduced further the division fails
        if p * k[0] + q * k[1] <= p * qk[0] + q * qk[1]:
            return None
        Q[qk] = c
        R = psub(R, pmul({qk: c}, D))
    return None


def elementary_candidates(P, pq_max=40, tol=1e-9, npts=4):
    """Numeric candidates (p, q, sign) with D = M + sign*M^{-1} dividing P:
    P must vanish on the locus  (p alpha + q beta)/2 = 0 (mod pi)   [sign -]
    or = pi/2 (mod pi) [sign +]. Sampled at npts real points per locus."""
    # scale reference
    import random
    rng = random.Random(20260731)
    ref = max(abs(peval(P, rng.uniform(0.3, 2.5), rng.uniform(0.3, 2.5)))
              for _ in range(6))
    out = []
    seen = set()
    for p in range(0, pq_max + 1):
        for q in range(-pq_max, pq_max + 1):
            if p == 0 and q <= 0:
                continue
            if math.gcd(p, abs(q)) != 1 and not (p == 0 or q == 0):
                pass  # non-primitive directions allowed: frequencies matter
            for sign, off in ((-1, 0.0), (1, math.pi)):
                key = (p, q, sign)
                if key in seen:
                    continue
                good = True
                for t in range(npts):
                    al = 0.35 + 0.53 * t
                    # (p al + q be)/2 = off/2 (mod pi): pick k spread out
                    k = t - npts // 2
                    if q != 0:
                        be = (off + 2 * math.pi * k - p * al) / q
                    else:
                        if p == 0:
                            good = False
                            break
                        al = (off + 2 * math.pi * k) / p
                        be = 0.7 + 0.31 * t
                    v = abs(peval(P, al, be))
                    if v > tol * max(ref, 1e-30):
                        good = False
                        break
                if good:
                    seen.add(key)
                    out.append((p, q, sign))
    return out


def factor_elementary(P, pq_max=40):
    """Peel exact elementary factors off P greedily. Returns
    (factors, cofactor) with P == unit-monomial-scaled product; each factor
    reported as ('sin'|'cos', p, q) meaning sin/cos((p a + q b)/2).
    Division is EXACT; the numeric step only proposes candidates.
    The reconstruction P == c * prod * cofactor is re-verified exactly."""
    factors = []
    cur = dict(P)
    changed = True
    while changed and cur:
        changed = False
        cands = elementary_candidates(cur, pq_max)
        # peel the largest factors first, so composite loci (e.g. sin(4b)
        # = product of four half-angle factors) come out in maximal form
        cands.sort(key=lambda t: -(t[0] * t[0] + t[1] * t[1]))
        for (p, q, sign) in cands:
            while True:
                Q = divide_binomial(cur, p, q, sign)
                if Q is None:
                    break
                factors.append(("cos" if sign == 1 else "sin", p, q))
                cur = Q
                changed = True
    # exact reconstruction check
    rec = {k: c for k, c in cur.items()}
    for (kind, p, q) in factors:
        sign = 1 if kind == "cos" else -1
        rec = pmul(rec, padd(mono(p, q), mono(-p, -q, (F(sign), F(0)))))
    # rec should equal P (both are the raw binomial products; the sin/cos
    # normalization constants (2i resp. 2) are carried by the caller)
    assert pclean(psub(rec, P)) == {}, "reconstruction failed"
    return factors, cur


def fmt_factors(factors, cofactor):
    def ang(p, q):
        parts = []
        if p:
            parts.append(f"{p}a" if p != 1 else "a")
        if q:
            parts.append(f"{'+' if q > 0 else '-'}{abs(q)}b"
                         if abs(q) != 1 else ("+b" if q > 0 else "-b"))
        return "(" + "".join(parts).lstrip("+") + ")/2"
    s = " * ".join(f"{k}{ang(p, q)}" for (k, p, q) in factors)
    if cofactor and pclean(cofactor) not in ({}, ):
        if len(cofactor) == 1 and list(cofactor)[0] == (0, 0):
            c = cofactor[(0, 0)]
            s += f" * ({c[0]}{'+' if c[1] >= 0 else ''}{c[1]}i)"
        else:
            s += f" * [{DS.trig_str(pre(cofactor)) if is_real(cofactor) else 'cofactor(' + str(len(cofactor)) + ' terms)'}]"
    return s or "1"


def is_real(p):
    return pclean(psub(p, pconj(p))) == {}


def pre(p):
    return DS.pre(p)


# ------------------------------------------------------------- subcommands

def binding_report(word, gamma, corner=None, out=None):
    gd = glide_data_word(tuple(word))
    th = math.radians(180.0 - gamma)
    # float argmax on the arc
    best, arg = -math.inf, None
    als = [math.radians(a) for a in arc_alphas(gamma, uniform=200, depth=14)]
    if corner is not None:
        for d in range(1, 46):
            eps = 0.5 * 2.0 ** -d
            for a in (corner - eps, corner + eps):
                aa = math.radians(a)
                if 0 < aa < th:
                    als.append(aa)
    mval_p = gd["m"]
    dconj = pconj(gd["delta"])
    for al in als:
        be = th - al
        m = peval(mval_p, al, be).real
        lo, hi = -math.inf, math.inf
        for (u, v) in gd["gates"]:
            pu = peval(pim(pmul(dconj, u)), al, be).real
            pv = peval(pim(pmul(dconj, v)), al, be).real
            if pu > pv:
                pu, pv = pv, pu
            lo, hi = max(lo, pu), min(hi, pv)
        wdt = min(hi, 2 * m - lo) - max(lo, 2 * m - hi)
        if wdt > best:
            best, arg = wdt, (al, be)
    al, be = arg
    m = peval(mval_p, al, be).real
    rows = []
    for (k, lbl, pz) in gd["endpoints"]:
        rows.append((k, lbl, pz, peval(pz, al, be).real - m))
    rows.sort(key=lambda r: abs(r[3]))
    print(f"word {''.join(map(str, word))}  gamma={gamma}  "
          f"argmax alpha={math.degrees(al):.6f} width={best:.3e}")
    res = {"word": "".join(map(str, word)), "gamma": gamma,
           "argmax_alpha_deg": math.degrees(al), "width": best,
           "bindings": []}
    seen_forms = set()
    for (k, lbl, pz, dv) in rows[:6]:
        diff = pclean(psub(pz, mval_p))
        key = json.dumps(sorted([(mn, (str(c[0]), str(c[1])))
                                 for mn, c in diff.items()]))
        if key in seen_forms:
            continue
        seen_forms.add(key)
        if not diff:
            desc = "0 (endpoint ON the axis offset identically)"
        else:
            fac, cof = factor_elementary(diff)
            desc = fmt_factors(fac, cof)
        print(f"  gate {k}{lbl}: p - m = {dv:+.3e}   factors: {desc}")
        res["bindings"].append({"gate": k, "end": lbl, "value": dv,
                                "factorization": desc})
    if out:
        json.dump(res, open(out, "w"), indent=1)
    return res


def cmd_bind(args):
    w = tuple(int(c) for c in args.word)
    binding_report(w, args.gamma, corner=args.corner, out=args.out)


def cmd_selftest(_args):
    ok = True
    # I2/I3 for W(3,3): p(A1) - m = cos(a alpha) sin((b+1) beta) sin(a+b)
    # in half-angle vocabulary: cos((2a alpha)/2) sin((2(b+1) beta)/2)
    # sin((2 alpha + 2 beta)/2), each carrying its normalization constant.
    a, b = 3, 3
    word = ((0,) + (1, 2) * a + (0, 2) * b) * 2
    gd = glide_data_word(word)
    # gate 2 endpoint v is A1 (R0(A)); find its projection minus m
    k, lbl, pz = [e for e in gd["endpoints"] if e[0] == 2 and e[1] == "v"][0]
    diff = pclean(psub(pz, gd["m"]))
    fac, cof = factor_elementary(diff)
    want = sorted([("cos", 2 * a, 0), ("sin", 0, 2 * (b + 1)),
                   ("sin", 2, 2)])
    got = sorted(fac)
    # cofactor must be the constant scale: sin*sin -> (1/2i)^2, cos -> 1/2
    const_ok = (len(cof) == 1 and (0, 0) in cof)
    good = got == want and const_ok
    ok &= good
    print(f"[I2 factorization W(3,3): {fmt_factors(fac, cof)}]",
          "OK" if good else f"FAIL (got {got}, want {want})")
    # numeric agreement of the factored form at a random angle
    al, be = 0.41, 0.35
    lhs = peval(diff, al, be)
    rhs = (math.cos(a * al) * math.sin((b + 1) * be) * math.sin(al + be))
    c = cof[(0, 0)]
    cval = complex(c[0], c[1])
    # each sin binomial carries 2i, cos carries 2:
    scale = cval * (2j) * (2j) * 2 / (2j) ** 0  # sin,sin,cos -> (2i)(2i)(2)
    good = abs(lhs - rhs * (scale / abs(scale)) * abs(scale)) < 1e-9 or \
        abs(abs(lhs) - abs(rhs)) < 1e-9
    ok &= good
    print(f"[I2 numeric magnitude {abs(lhs):.6f} vs {abs(rhs):.6f}]",
          "OK" if good else "FAIL")
    print("SELFTEST", "PASSED" if ok else "FAILED")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("selftest").set_defaults(fn=cmd_selftest)
    c = sub.add_parser("bind")
    c.add_argument("--word", required=True)
    c.add_argument("--gamma", type=float, required=True)
    c.add_argument("--corner", type=float, default=None,
                   help="corner alpha in deg for deep accumulation")
    c.add_argument("--out")
    c.set_defaults(fn=cmd_bind)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
