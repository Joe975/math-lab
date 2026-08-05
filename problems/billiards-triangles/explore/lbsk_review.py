#!/usr/bin/env python3
"""Attempt 012: skeptic pass on attempt 011 (Lean LaurentBlock.lean) --
independent re-derivation, written from the LEAN FILE's definitions (not
from plaw_general.py, whose ring this script does not import).

Layer 1 (SYMBOLIC, exact): a from-scratch Laurent-polynomial ring over Q in
SEVEN variables (x, y, P, Q, i, z1, z2) -- i is a formal unit exactly as in
the Lean statements (no i^2 = -1).  Every division in the Lean definitions
is by a unit monomial (2, 2i, generator inverses), so the whole block lives
inside the Laurent ring and each identity is checked by exact polynomial
subtraction.  This re-proves the Lean statements independently of both the
Lean kernel and 005's/plaw_general's code.

Layer 2 (FLOAT, mapping): instantiate x=e^{ia_}, y=e^{ib_}, P=e^{i a a_},
Q=e^{i b b_}, i=1j and check the Lean-side expressions against 005's OWN
prose formulas (trig right-hand sides of I1-I4/D1/D2, Step-2 reflection
formulas, Step-3 closed forms mu/delta/w/tau/m, p(v)=Im(conj delta v)).
This nails the statement mapping: the Lean theorems at the intended
instance say exactly what 005 says.

stdlib only; deterministic (fixed seed).  Usage:
    python3 lbsk_review.py [--out data/lbsk_review.json]
"""
import argparse
import cmath
import json
import math
import random
import sys
from fractions import Fraction as Fr

# ---------------------------------------------------------------- Layer 1
# poly: dict {(ex,ey,eP,eQ,ei,ez1,ez2): Fraction}, Laurent exponents.

NV = 7
X = 0; Y = 1; PP = 2; QQ = 3; II = 4; Z1 = 5; Z2 = 6

def mono(coef=Fr(1), **exps):
    e = [0]*NV
    for k, v in exps.items():
        e[{'x':X,'y':Y,'P':PP,'Q':QQ,'i':II,'z1':Z1,'z2':Z2}[k]] = v
    return {tuple(e): Fr(coef)}

ONE = mono()
def const(c): return mono(Fr(c))

def padd(p, q):
    r = dict(p)
    for k, c in q.items():
        r[k] = r.get(k, Fr(0)) + c
        if r[k] == 0: del r[k]
    return r

def pneg(p): return {k: -c for k, c in p.items()}
def psub(p, q): return padd(p, pneg(q))

def pmul(p, q):
    r = {}
    for k1, c1 in p.items():
        for k2, c2 in q.items():
            k = tuple(a+b for a, b in zip(k1, k2))
            r[k] = r.get(k, Fr(0)) + c1*c2
            if r[k] == 0: del r[k]
    return r

def is_mono(p): return len(p) == 1

def pinv(p):
    """Inverse of a unit monomial."""
    assert is_mono(p), "inverting a non-monomial: encoding broken"
    (k, c), = p.items()
    return {tuple(-e for e in k): Fr(1)/c}

def pdiv_unit(p, u):  # divide by unit monomial u
    return pmul(p, pinv(u))

def ppow(p, n):
    r = ONE
    for _ in range(n): r = pmul(r, p)
    return r

def zero(p): return p == {}

# generators
gx, gy, gP, gQ, gi = (mono(x=1), mono(y=1), mono(P=1), mono(Q=1), mono(i=1))
gz1, gz2 = mono(z1=1), mono(z2=1)
TWO = const(2)

# --- Lean definitions, transcribed line by line from LaurentBlock.lean ---
def sinF(z, i):  # (z - z^-1)/(2*i)
    return pdiv_unit(psub(z, pinv(z)), pmul(TWO, i))

def cosF(z):     # (z + z^-1)/2
    return pdiv_unit(padd(z, pinv(z)), TWO)

def Bd(x, y, i): return sinF(pmul(x, y), i)
def Cd(x, y, i): return pmul(sinF(y, i), x)
def A1d(x, y, i): return pmul(Bd(x, y, i), psub(ONE, ppow(pinv(y), 2)))
def mud(y, P, Q):
    return pmul(pmul(ppow(pinv(y), 2), ppow(pinv(P), 2)), ppow(Q, 2))
def deltad(y, P, Q): return pmul(pmul(pinv(y), pinv(P)), Q)
def wd(x, y, P, Q, i):
    inner = psub(padd(psub(ONE, ppow(pinv(y), 2)),
                      pmul(ppow(pinv(y), 2), ppow(pinv(P), 2))),
                 mud(y, P, Q))
    return pmul(Bd(x, y, i), inner)
def taud(x, y, P, Q, i):
    return padd(wd(x, y, P, Q, i),
                pmul(mud(y, P, Q),
                     wd(pinv(x), pinv(y), pinv(P), pinv(Q), pneg(i))))
def projd(y, P, Q, i, v, vc):
    num = psub(pmul(deltad(pinv(y), pinv(P), pinv(Q)), v),
               pmul(deltad(y, P, Q), vc))
    return pdiv_unit(num, pmul(TWO, i))
def md(x, y, P, Q, i):
    wc = wd(pinv(x), pinv(y), pinv(P), pinv(Q), pneg(i))
    tc = taud(pinv(x), pinv(y), pinv(P), pinv(Q), pneg(i))
    t = taud(x, y, P, Q, i)
    v = psub(wd(x, y, P, Q, i), pdiv_unit(t, TWO))
    vc = psub(wc, pdiv_unit(tc, TWO))
    return pdiv_unit(projd(y, P, Q, i, v, vc), TWO)
def v2d(x, y, P, i):
    B_ = Bd(x, y, i)
    return padd(B_, pmul(ppow(pinv(y), 2),
                         psub(pmul(ppow(pinv(P), 2),
                                   Cd(pinv(x), pinv(y), pneg(i))), B_)))

# maps on pairs
def reflAB(p): return (p[1], p[0])
def reflCA(x, p): return (pmul(ppow(x, 2), p[1]), pmul(ppow(pinv(x), 2), p[0]))
def reflBC(x, y, i, p):
    B_ = Bd(x, y, i); Bc = Bd(pinv(x), pinv(y), pneg(i))
    return (padd(B_, pmul(ppow(pinv(y), 2), psub(p[1], B_))),
            padd(Bc, pmul(ppow(pinv(pinv(y)), 2), psub(p[0], Bc))))
def rotA(P, p): return (pmul(ppow(P, 2), p[0]), pmul(ppow(pinv(P), 2), p[1]))
def rotB(x, y, Q, i, p):
    B_ = Bd(x, y, i); Bc = Bd(pinv(x), pinv(y), pneg(i))
    return (padd(B_, pmul(ppow(pinv(Q), 2), psub(p[0], B_))),
            padd(Bc, pmul(ppow(pinv(pinv(Q)), 2), psub(p[1], Bc))))

def run_symbolic():
    x, y, P, Q, i = gx, gy, gP, gQ, gi
    xc, yc, Pc, Qc, ic = pinv(x), pinv(y), pinv(P), pinv(Q), pneg(i)
    checks = {}
    m_ = md(x, y, P, Q, i)
    pr = lambda v, vc: projd(y, P, Q, i, v, vc)
    # I1
    lhs = psub(m_, pr(Cd(x, y, i), Cd(xc, yc, ic)))
    rhs = pneg(padd(pmul(pmul(cosF(P), sinF(x, i)), sinF(Q, i)),
                    pmul(pmul(sinF(P, i), sinF(y, i)),
                         cosF(pmul(pmul(x, y), Q)))))
    checks['I1'] = zero(psub(lhs, rhs))
    # I2
    lhs = psub(pr(A1d(x, y, i), A1d(xc, yc, ic)), m_)
    rhs = pmul(pmul(cosF(P), sinF(pmul(y, Q), i)), sinF(pmul(x, y), i))
    checks['I2'] = zero(psub(lhs, rhs))
    # I3
    lhs = psub(pr(Bd(x, y, i), Bd(xc, yc, ic)), m_)
    rhs = pmul(pmul(cosF(pmul(y, Q)), sinF(P, i)), sinF(pmul(x, y), i))
    checks['I3'] = zero(psub(lhs, rhs))
    # I4
    lhs = psub(pr(v2d(x, y, P, i), v2d(xc, yc, Pc, ic)), m_)
    rhs = psub(pmul(pmul(cosF(P), sinF(x, i)), sinF(Q, i)),
               pmul(pmul(sinF(P, i), sinF(y, i)),
                    cosF(pmul(pmul(x, y), Q))))
    checks['I4'] = zero(psub(lhs, rhs))
    # D1
    lhs = psub(pr(Bd(x, y, i), Bd(xc, yc, ic)),
               pr(Cd(x, y, i), Cd(xc, yc, ic)))
    rhs = pmul(sinF(x, i), sinF(pmul(P, pinv(Q)), i))
    checks['D1'] = zero(psub(lhs, rhs))
    # D2
    lhs = psub(pr(A1d(x, y, i), A1d(xc, yc, ic)),
               pr(Cd(x, y, i), Cd(xc, yc, ic)))
    rhs = pneg(pmul(sinF(y, i),
                    sinF(pmul(pmul(pmul(pinv(x), pinv(y)), P), pinv(Q)), i)))
    checks['D2'] = zero(psub(lhs, rhs))
    # glide facts
    checks['mu_eq_delta_sq'] = zero(psub(mud(y, P, Q),
                                         ppow(deltad(y, P, Q), 2)))
    checks['tau_parallel_axis'] = zero(
        pr(taud(x, y, P, Q, i), taud(xc, yc, Pc, Qc, ic)))
    checks['mu_mul_conj_mu'] = zero(psub(
        pmul(mud(y, P, Q), mud(yc, Pc, Qc)), ONE))
    # half_word_sq_is_translation (z1 formal)
    wc = wd(xc, yc, Pc, Qc, ic)
    lhs = padd(pmul(mud(y, P, Q), padd(pmul(mud(yc, Pc, Qc), gz1), wc)),
               wd(x, y, P, Q, i))
    checks['half_word_sq'] = zero(psub(lhs, padd(gz1, taud(x, y, P, Q, i))))
    # glide_action (z1, z2 formal)
    lhs = pr(padd(pmul(mud(y, P, Q), gz2), wd(x, y, P, Q, i)),
             padd(pmul(mud(yc, Pc, Qc), gz1), wc))
    rhs = psub(pmul(TWO, m_), pr(gz1, gz2))
    checks['glide_action'] = zero(psub(lhs, rhs))
    # half_word_closed_form (pair (z1,z2) formal)
    got = reflBC(x, y, i, rotA(P, rotB(x, y, Q, i, (gz1, gz2))))
    want = (padd(pmul(mud(y, P, Q), gz2), wd(x, y, P, Q, i)),
            padd(pmul(mud(yc, Pc, Qc), gz1), wc))
    checks['half_word_closed_form.1'] = zero(psub(got[0], want[0]))
    checks['half_word_closed_form.2'] = zero(psub(got[1], want[1]))
    # pair collapses
    got = reflCA(x, reflAB((gz1, gz2)))
    want = rotA(x, (gz1, gz2))
    checks['pair_collapse_A'] = (zero(psub(got[0], want[0]))
                                 and zero(psub(got[1], want[1])))
    got = reflBC(x, y, i, reflAB((gz1, gz2)))
    want = rotB(x, y, y, i, (gz1, gz2))
    checks['pair_collapse_B'] = (zero(psub(got[0], want[0]))
                                 and zero(psub(got[1], want[1])))
    # half_word_letters for a,b in 0..3 (iterates computed symbolically)
    for a in range(4):
        for b in range(4):
            p = (gz1, gz2)
            for _ in range(b):
                p = reflBC(x, y, i, reflAB(p))
            for _ in range(a):
                p = reflCA(x, reflAB(p))
            p = reflBC(x, y, i, p)
            Pa, Qb = ppow(x, a), ppow(y, b)
            Pac, Qbc = pinv(Pa), pinv(Qb)
            want = (padd(pmul(mud(y, Pa, Qb), gz2), wd(x, y, Pa, Qb, i)),
                    padd(pmul(mud(yc, Pac, Qbc), gz1),
                         wd(xc, yc, Pac, Qbc, ic)))
            checks[f'half_word_letters[{a},{b}]'] = (
                zero(psub(p[0], want[0])) and zero(psub(p[1], want[1])))
    return checks

# ---------------------------------------------------------------- Layer 2
def crun(f, *args): return f(*args)

def run_float(trials=200, seed=20260804):
    rng = random.Random(seed)
    worst = {}
    def upd(tag, err):
        worst[tag] = max(worst.get(tag, 0.0), err)

    # numeric versions of the Lean defs (complex arithmetic, i passed in)
    def nsinF(z, i): return (z - 1/z) / (2*i)
    def ncosF(z): return (z + 1/z) / 2
    def nB(x, y, i): return nsinF(x*y, i)
    def nC(x, y, i): return nsinF(y, i) * x
    def nA1(x, y, i): return nB(x, y, i) * (1 - (1/y)**2)
    def nmu(y, P, Q): return (1/y)**2 * (1/P)**2 * Q**2
    def ndelta(y, P, Q): return (1/y) * (1/P) * Q
    def nw(x, y, P, Q, i):
        return nB(x, y, i) * (1 - y**-2 + y**-2 * P**-2 - nmu(y, P, Q))
    def ntau(x, y, P, Q, i):
        return (nw(x, y, P, Q, i)
                + nmu(y, P, Q) * nw(1/x, 1/y, 1/P, 1/Q, -i))
    def nproj(y, P, Q, i, v, vc):
        return (ndelta(1/y, 1/P, 1/Q)*v - ndelta(y, P, Q)*vc) / (2*i)
    def nm(x, y, P, Q, i):
        wc = nw(1/x, 1/y, 1/P, 1/Q, -i)
        tc = ntau(1/x, 1/y, 1/P, 1/Q, -i)
        return nproj(y, P, Q, i,
                     nw(x, y, P, Q, i) - ntau(x, y, P, Q, i)/2,
                     wc - tc/2) / 2
    def nv2(x, y, P, i):
        return (nB(x, y, i)
                + y**-2 * (P**-2 * nC(1/x, 1/y, -i) - nB(x, y, i)))

    for _ in range(trials):
        a = rng.randint(1, 12)
        b = rng.randint(1, a)
        al = rng.uniform(0.05, 1.3)
        be = rng.uniform(0.05, min(1.3, math.pi - 0.1 - al))
        x = cmath.exp(1j*al); y = cmath.exp(1j*be)
        P = cmath.exp(1j*a*al); Q = cmath.exp(1j*b*be)
        i = 1j
        Bg = math.sin(al+be)          # 005: B = sin(alpha+beta), real
        # -- the trig dictionary of the Lean statements vs 005's prose
        upd('cosF_P=cos(a.al)', abs(ncosF(P) - math.cos(a*al)))
        upd('sinF_P=sin(a.al)', abs(nsinF(P, i) - math.sin(a*al)))
        upd('sinF_Q=sin(b.be)', abs(nsinF(Q, i) - math.sin(b*be)))
        upd('sinF_x=sin(al)', abs(nsinF(x, i) - math.sin(al)))
        upd('sinF_y=sin(be)', abs(nsinF(y, i) - math.sin(be)))
        upd('sinF_xy=sin(al+be)', abs(nsinF(x*y, i) - math.sin(al+be)))
        upd('cosF_xyQ=cos(al+(b+1)be)',
            abs(ncosF(x*y*Q) - math.cos(al+(b+1)*be)))
        upd('sinF_yQ=sin((b+1)be)', abs(nsinF(y*Q, i) - math.sin((b+1)*be)))
        upd('cosF_yQ=cos((b+1)be)', abs(ncosF(y*Q) - math.cos((b+1)*be)))
        upd('sinF_PQinv=sin(a.al-b.be)',
            abs(nsinF(P/Q, i) - math.sin(a*al-b*be)))
        upd('sinF_D2mono=sin((a-1)al-(b+1)be)',
            abs(nsinF(P/(x*y*Q), i) - math.sin((a-1)*al-(b+1)*be)))
        # -- 005 Step 2/3 objects vs Lean defs at the instance
        upd('B=sin(al+be)', abs(nB(x, y, i) - Bg))
        upd('C=sin(be)e^{i.al}', abs(nC(x, y, i)
                                     - math.sin(be)*cmath.exp(1j*al)))
        mu5 = cmath.exp(-2j*a*al + 2j*(b-1)*be)
        upd('mu=e^{-2ia.al+2i(b-1)be}', abs(nmu(y, P, Q) - mu5))
        de5 = cmath.exp(1j*(-a*al + (b-1)*be))
        upd('delta', abs(ndelta(y, P, Q) - de5))
        w5 = Bg * (1 + cmath.exp(-2j*(a*al+be)) - cmath.exp(-2j*be) - mu5)
        upd('w(005 Step3)', abs(nw(x, y, P, Q, i) - w5))
        tau5 = w5 + mu5 * w5.conjugate()
        upd('tau', abs(ntau(x, y, P, Q, i) - tau5))
        m5 = ((de5.conjugate() * (w5 - tau5/2)).imag) / 2
        upd('m=Im(conj.delta(w-tau/2))/2', abs(nm(x, y, P, Q, i) - m5))
        # p(v) = Im(conj delta . v) at a random point, pair = (v, conj v)
        v = complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
        upd('proj=Im(conj.delta.v)',
            abs(nproj(y, P, Q, i, v, v.conjugate())
                - (de5.conjugate()*v).imag))
        # v2 = R0(RotA(2a.al) C) computed from 005's geometric maps
        R0 = lambda z: Bg + cmath.exp(-2j*be) * (z - Bg).conjugate()
        v2g = R0(cmath.exp(2j*a*al) * math.sin(be)*cmath.exp(1j*al))
        upd('v2=R0(RotA C)', abs(nv2(x, y, P, i) - v2g))
        # closed form u(z) vs direct composition of 005's maps
        z = complex(rng.uniform(-1, 1), rng.uniform(-1, 1))
        zz = Bg + cmath.exp(-2j*b*be) * (z - Bg)         # Rot_B(-2b.be)
        zz = cmath.exp(2j*a*al) * zz                     # Rot_A(2a.al)
        zz = R0(zz)                                      # R_0
        upd('u=mu.conj z+w', abs(nmu(y, P, Q)*z.conjugate()
                                 + nw(x, y, P, Q, i) - zz))
        # pair second components genuinely conjugate at the instance
        upd('pair C', abs(nC(1/x, 1/y, -i) - nC(x, y, i).conjugate()))
        upd('pair w', abs(nw(1/x, 1/y, 1/P, 1/Q, -i)
                          - nw(x, y, P, Q, i).conjugate()))
        upd('pair v2', abs(nv2(1/x, 1/y, 1/P, -i)
                           - nv2(x, y, P, i).conjugate()))
        # -- 005's I1-I4/D1/D2 trig RHS vs the Lean-statement LHS
        m_ = nm(x, y, P, Q, i)
        pC = nproj(y, P, Q, i, nC(x, y, i), nC(1/x, 1/y, -i))
        pB = nproj(y, P, Q, i, nB(x, y, i), nB(1/x, 1/y, -i))
        pA1 = nproj(y, P, Q, i, nA1(x, y, i), nA1(1/x, 1/y, -i))
        pv2 = nproj(y, P, Q, i, nv2(x, y, P, i), nv2(1/x, 1/y, 1/P, -i))
        S = math.cos(a*al)*math.sin(al)*math.sin(b*be)
        T = math.sin(a*al)*math.sin(be)*math.cos(al+(b+1)*be)
        upd('I1 trig', abs((m_ - pC) - (-(S + T))))
        upd('I2 trig', abs((pA1 - m_)
                           - math.cos(a*al)*math.sin((b+1)*be)
                           * math.sin(al+be)))
        upd('I3 trig', abs((pB - m_)
                           - math.cos((b+1)*be)*math.sin(a*al)
                           * math.sin(al+be)))
        upd('I4 trig', abs((pv2 - m_) - (S - T)))
        upd('D1 trig', abs((pB - pC)
                           - math.sin(al)*math.sin(a*al-b*be)))
        upd('D2 trig', abs((pA1 - pC)
                           - (-math.sin(be)
                              * math.sin((a-1)*al-(b+1)*be))))
    return worst

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out")
    args = ap.parse_args()
    print("=== Layer 1: symbolic (exact, 7-var Laurent ring over Q) ===")
    checks = run_symbolic()
    ok = True
    for k, v in checks.items():
        if not v:
            print(f"  FAIL {k}")
            ok = False
    n = len(checks)
    print(f"  {sum(checks.values())}/{n} identities are exact zero "
          f"polynomials (i FORMAL, no i^2=-1)")
    print("=== Layer 2: float mapping to 005's prose (200 random "
          "(a,b,al,be), a<=12) ===")
    worst = run_float()
    w = max(worst.values())
    for k in sorted(worst, key=worst.get, reverse=True)[:5]:
        print(f"  worst[{k}] = {worst[k]:.3e}")
    print(f"  overall worst error = {w:.3e} over {len(worst)} mapped "
          f"quantities")
    ok = ok and w < 1e-9
    print("SKEPTIC-011 SYMBOLIC+MAPPING:", "ALL OK" if ok else "FAILURES")
    if args.out:
        json.dump({"symbolic_checks": checks,
                   "n_symbolic": n,
                   "n_symbolic_pass": sum(checks.values()),
                   "float_worst_errors": worst,
                   "float_overall_worst": w,
                   "ok": ok,
                   "meaning": ("Layer 1: every Lean theorem statement of "
                               "011 re-proven as an exact zero polynomial "
                               "in an independent 7-variable Laurent ring "
                               "over Q with i FORMAL; Layer 2: the Lean "
                               "statements at the torus instance match "
                               "005's prose quantities to float "
                               "precision")},
                  open(args.out, "w"), indent=1)
    sys.exit(0 if ok else 1)

if __name__ == "__main__":
    main()
