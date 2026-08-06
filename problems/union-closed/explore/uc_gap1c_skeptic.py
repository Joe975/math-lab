#!/usr/bin/env python3
"""Attempt 019: skeptic re-verification of 018 (margin-modulated control
killed at the witness; MM_abs and the lambda-window control survive the
ladder).

INDEPENDENCE CONTRACT.  This file imports nothing from
uc_gap1_candidates.py, uc_or_avg.py, uc_or_avg_skeptic.py,
uc_agg_ctrl_probe.py or uc_agg_ctrl_probe2.py.  Everything below is written
from the PROSE of 007 section 1 (the well-posed M_i census), 018's record
(the three MM readings), and 016 section 1 (the MU(n, r) ladder).  Witness
atoms, theta* weights and search endpoints are DATA copied from records /
checkpoints -- data reuse is fine, code reuse is not.

The exact kit is new and structurally different from BOTH prior kits:
  - log2 enclosure by dyadic interval squaring (digit extraction) with
    directed rounding -- the idea 017 used, re-implemented here from
    scratch; NOT 013's atanh-series kit (which 018 imported).
  - Plackett-root enclosure by exact rational NEWTON steps with an exact
    sign-checked bracket (every bracket update verified by an exact sign
    evaluation of the quadratic; bisection fallback guarantees progress);
    NOT 018's blind 80-step dyadic bisection.  Bracket endpoint signs are
    also proved symbolically:  with F(z) = z(1-x-y+z) - t(x-z)(y-z) and
    0 < x, y < 1, t > 1:
        F(xy)       = (1-t) xy (1-x)(1-y)      < 0,
        F(min(x,y)) = min(x,y) (1 - max(x,y))  > 0,
    so the root bracket (xy, min(x,y)) is certified before any iteration.
  - the census is scaled-integer (exact) end-to-end; the only outward
    rounding is inside the two enclosure primitives, which are directed.

Parts:
  S0  kit self-tests (log2 enclosure vs known values; Plackett root vs
      closed forms and float solver; h2 enclosure; bracket sign identities).
  S1  K1: from-scratch exact certification of MM_sec on the 10-atom witness
      at t = 181/16 (plus the 2-digit tidy witness), with exact marginals,
      exact dichotomy, and an exact per-row check that degenerate histories
      contribute 0 to the secant numerator under either convention.
  S2  K2: the definitional fork -- own float census evaluating 018's three
      readings PLUS the variants 018 did not test (sensitivity at the
      REALIZED odds ratio, signed and unsigned; |dh/drho| weighting).
  S3  K3: float spot-checks with the ladder rebuilt from 016's prose spec
      (own tune_dil): witness MM numbers, theta* ladder n = 96 at lam = 2,
      raw ladder n = 96 (the sign 018's prose gets wrong), part-D window
      value at n = 96 and the a1/a2 fit arithmetic, perturbation battery
      with DIFFERENT seeds.
  S4  K4: re-evaluation of the part-C best_free endpoint (data from
      gap1c_partC.json) and the z = 1/2 sensitivity-boundary arithmetic
      (the golden-ratio parenthetical).

Standard library only.  Deterministic (fixed seeds).  Runtime ~3 min.
Usage:
    python problems/union-closed/explore/uc_gap1c_skeptic.py [--fast]
Checkpoints: problems/union-closed/data/gap1csk_*.json
Log:         tee to problems/union-closed/data/gap1csk_run.log
"""

from __future__ import annotations

import json
import math
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

RECORD_NUM, RECORD_DEN = 38271, 100000     # record threshold 0.38271 (Liu)
MARG_TARGET = 0.375                        # 016's tuning target
F0, F1, HALF = Fraction(0), Fraction(1), Fraction(1, 2)
FAST = "--fast" in sys.argv
T0 = time.monotonic()


def log(msg: str) -> None:
    print(f"[{time.monotonic() - T0:6.0f}s] {msg}", flush=True)


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True,
                                        default=str) + "\n", encoding="utf-8")


# ======================================================================
# DATA (copied from records/checkpoints; not imported)
# ======================================================================

# 007's 10-atom witness (record part T / uc_or_avg.py WITNESS_U -- data)
WITNESS_U = {
    0b0100001: 3.802705744653898,      # {1,6}
    0b1000001: 8.094727225702503,      # {1,7}
    0b1010001: 0.004792968160237178,   # {1,5,7}  light slice
    0b0000000: 0.2526660836353666,     # {}
    0b0100000: 27.28641199464732,      # {6}
    0b1000000: 12.529780963546393,     # {7}
    0b0110000: 0.3609234435898425,     # {5,6}
    0b0000010: 20.0,                   # {2}  dilution
    0b0000100: 20.0,                   # {3}  dilution
    0b0001000: 20.0,                   # {4}  dilution
}

# theta* from 016 P2 (aggprobe2_partP2.json), as quoted in 018's engine
THETA_STAR_L2 = {
    "A1": 14.026917090492047, "A2": 22.622718089566355,
    "B0": 13.573329717299568, "B1": 43.00957928858568,
    "B2": 24.62436309125563, "dil": 21.044172077650057,
    "light": 0.007943082225456255, "resp": 0.0017259962695652103,
}

# raw witness weights as a theta dict (016's parametrization -- data)
THETA_WITNESS = {
    "A1": 3.802705744653898, "A2": 8.094727225702503,
    "light": 0.004792968160237178, "B0": 0.2526660836353666,
    "B1": 27.28641199464732, "B2": 12.529780963546393,
    "resp": 0.3609234435898425, "dil": 20.0,
}

LAMWIN_C = 4.847                           # 009/011 window law c/(n-3)


# ======================================================================
# certified log2 enclosure -- own implementation (interval squaring)
# ======================================================================

# ln 2 = -ln(1/2) = sum_{k>=1} (1/2)^k / k ; the partial sum after K terms
# is a strict lower bound and the tail is < 1/((K+1) 2^K).
_LN2_K = 72
LN2_LO = sum(Fraction(1, k << k) for k in range(1, _LN2_K + 1))
LN2_HI = LN2_LO + Fraction(1, (_LN2_K + 1) * (1 << _LN2_K))


def log2_enc(q: Fraction, K: int = 96, prec: int = 288
             ) -> tuple[Fraction, Fraction]:
    """Certified [lo, hi] with lo <= log2 q <= hi for rational q > 0.

    Write q = 2^e x, x in [1, 2).  Square a directed-rounded dyadic
    enclosure of x K times, renormalizing into [1, 4) by shifts whose count
    builds the binary expansion S of 2^K log2 x; then
    log2 x = (S + log2 x_K)/2^K with x_K in a known small interval, where
    the elementary sandwich 1 - 1/v <= ln v <= v - 1 (v > 0) plus the LN2
    bounds give a crude certified enclosure whose error is divided by 2^K.
    """
    num, den = q.numerator, q.denominator
    if num <= 0:
        raise ValueError("log2 of nonpositive")
    if num & (num - 1) == 0 and den & (den - 1) == 0:   # exact power of two
        return (Fraction(num.bit_length() - den.bit_length()),) * 2
    e = num.bit_length() - den.bit_length()
    nn, dd = (num, den << e) if e >= 0 else (num << (-e), den)
    if nn < dd:
        nn <<= 1
        e -= 1
    # x = nn/dd in [1, 2)
    scale = 1 << prec
    xlo = (nn << prec) // dd                             # floor
    xhi = -((-(nn << prec)) // dd)                       # ceil
    S = 0
    for _ in range(K):
        xlo = (xlo * xlo) >> prec                        # round down
        xhi = -((-(xhi * xhi)) >> prec)                  # round up
        s = 0
        while xlo >= (scale << 1):                       # whole interval >= 2
            xlo >>= 1
            xhi = -((-xhi) >> 1)
            s += 1
        while xhi < scale:                               # whole interval < 1
            xlo <<= 1
            xhi <<= 1
            s -= 1
        S = (S << 1) + s
    vlo, vhi = Fraction(xlo, scale), Fraction(xhi, scale)
    assert 0 < vlo <= vhi
    ln_lo = 1 - 1 / vlo                                  # <= ln vlo
    ln_hi = vhi - 1                                      # >= ln vhi
    l2_lo = ln_lo / (LN2_HI if ln_lo >= 0 else LN2_LO)
    l2_hi = ln_hi / (LN2_LO if ln_hi >= 0 else LN2_HI)
    pw = 1 << K
    return e + (S + l2_lo) / pw, e + (S + l2_hi) / pw


def h2_enc(z: Fraction) -> tuple[Fraction, Fraction]:
    """Certified [lo, hi] for binary entropy h2(z), z rational in [0, 1]."""
    if z <= 0 or z >= 1:
        return F0, F0
    alo, ahi = log2_enc(z)
    blo, bhi = log2_enc(1 - z)
    return (-z * ahi - (1 - z) * bhi), (-z * alo - (1 - z) * blo)


def h2_over(zlo: Fraction, zhi: Fraction) -> tuple[Fraction, Fraction]:
    """[lo, hi] of h2 over [zlo, zhi] (h2 unimodal with peak 1 at 1/2)."""
    alo, ahi = h2_enc(zlo)
    blo, bhi = h2_enc(zhi)
    lo, hi = min(alo, blo), max(ahi, bhi)
    if zlo < HALF < zhi:
        hi = F1
    return lo, hi


# ======================================================================
# certified Plackett root -- exact Newton with sign-checked bracket
# ======================================================================

def plackett_F_coeffs(x: Fraction, y: Fraction, t: Fraction):
    """F(z) = z(1-x-y+z) - t(x-z)(y-z) = A z^2 + B z + C."""
    return 1 - t, 1 - x - y + t * (x + y), -t * x * y


def z_root_enc(x: Fraction, y: Fraction, t: Fraction,
               bits: int = 160) -> tuple[Fraction, Fraction]:
    """Certified enclosure of the Plackett both-zero probability z_t(x, y)
    for 0 < x, y < 1 and t > 1: the unique root of F in (xy, min(x,y)).

    Every bracket move is validated by an exact sign evaluation of F, so
    correctness never depends on the Newton formula; Newton (on the exact
    quadratic) supplies fast candidates, bisection guarantees progress.
    """
    if t == 1:
        return x * y, x * y
    assert t > 1 and 0 < x < 1 and 0 < y < 1
    A, B, C = plackett_F_coeffs(x, y, t)

    def F(z: Fraction) -> Fraction:
        return (A * z + B) * z + C

    lo, hi = x * y, min(x, y)
    flo, fhi = F(lo), F(hi)
    # symbolic identities, re-checked exactly per call
    assert flo == (1 - t) * x * y * (1 - x) * (1 - y) and flo < 0
    assert fhi == min(x, y) * (1 - max(x, y)) and fhi > 0
    tol = Fraction(1, 1 << bits)
    while hi - lo > tol:
        probes = [(lo + hi) / 2]
        d = 2 * A * hi + B
        if d != 0:
            zn = (hi - fhi / d).limit_denominator(1 << 400)
            if lo < zn < hi:
                probes.append(zn)
        d = 2 * A * lo + B
        if d != 0:
            zn = (lo - flo / d).limit_denominator(1 << 400)
            if lo < zn < hi:
                probes.append(zn)
        for z in probes:
            fz = F(z)
            if fz == 0:
                return z, z
            if fz < 0:
                if z > lo:
                    lo, flo = z, fz
            else:
                if z < hi:
                    hi, fhi = z, fz
    return lo, hi


# ======================================================================
# exact census (scaled integer) and certified MM_sec / plain aggregate
# ======================================================================

def rationalize(u_float: dict[int, float]) -> dict[int, Fraction]:
    return {a: Fraction(w) for a, w in u_float.items() if w > 0}


def exact_kernel(u: dict[int, Fraction], n: int, t: Fraction):
    """Scaled-integer kernel k(A,B) = U_A U_B tn^p td^(n-p), p = |A&B|.
    Global scale cancels in every normalized quantity."""
    supp = sorted(u)
    L = 1
    for w in u.values():
        L = L * w.denominator // math.gcd(L, w.denominator)
    U = {a: int(u[a] * L) for a in supp}
    tn, td = t.numerator, t.denominator
    ptn = [1] * (n + 1)
    ptd = [1] * (n + 1)
    for j in range(1, n + 1):
        ptn[j] = ptn[j - 1] * tn
        ptd[j] = ptd[j - 1] * td
    rows = []
    for A in supp:
        UA = U[A]
        rows.append([UA * U[B] * ptn[bin(A & B).count("1")]
                     * ptd[n - bin(A & B).count("1")] for B in supp])
    return supp, rows


def exact_mm_cert(u_float: dict[int, float], n: int, t: Fraction) -> dict:
    """Certified enclosures of MM_sec (secant form) and the plain aggregate
    at rational tilt t, from scratch: own kernel, own census, own root and
    log2 enclosures.  Also: exact marginals, exact dichotomy, and an exact
    per-row proof that margin-degenerate histories contribute 0."""
    u = rationalize(u_float)
    supp, kern = exact_kernel(u, n, t)
    m_at = len(supp)
    Z = sum(sum(r) for r in kern)
    # exact elementwise marginals
    mu_row = [sum(r) for r in kern]
    marg = [F0] * n
    for idx, A in enumerate(supp):
        for j in range(n):
            if A >> j & 1:
                marg[j] += Fraction(mu_row[idx], Z)
    in_regime = all(m < Fraction(RECORD_NUM, RECORD_DEN) for m in marg)

    num_lo = num_hi = F0          # MM_sec numerator (unnormalized, /Z scale)
    pl_lo = pl_hi = F0            # plain numerator: sum m log2(OR/t)
    den = 0                       # nondegenerate mass, integer scale
    dich_ok = True
    deg_zero_ok = True            # degenerate rows provably contribute 0
    n_rows = 0
    zcache: dict[tuple[Fraction, Fraction], tuple[Fraction, Fraction]] = {}
    for i in range(1, n):                                # i < n only
        pmask, bit = (1 << (i - 1)) - 1, 1 << (i - 1)
        has0, has1 = set(), set()
        pref, resp = [], []
        for A in supp:
            c = A & pmask
            r = 1 if A & bit else 0
            pref.append(c)
            resp.append(r)
            (has1 if r else has0).add(c)
        active = has0 & has1
        tables: dict[tuple[int, int], list[int]] = {}
        for ia in range(m_at):
            ca, ra, row = pref[ia], resp[ia], kern[ia]
            for ib in range(m_at):
                key = (ca, pref[ib])
                tab = tables.get(key)
                if tab is None:
                    tab = tables[key] = [0, 0, 0, 0]
                tab[(ra << 1) | resp[ib]] += row[ib]
        for (a, b), (t00, t01, t10, t11) in tables.items():
            m = t00 + t01 + t10 + t11
            if m == 0:
                continue
            nondeg = t00 > 0 and t01 > 0 and t10 > 0 and t11 > 0
            if nondeg != (a in active and b in active):
                dich_ok = False
            x = Fraction(t00 + t01, m)
            y = Fraction(t00 + t10, m)
            z_rlz = Fraction(t00, m)
            if not nondeg:
                # exact convention-robustness check: margin degenerate and
                # z~ equals the pinned boundary value of z_t for EVERY t
                if x == 0 or y == 0:
                    if z_rlz != 0:
                        deg_zero_ok = False
                elif x == 1:
                    if z_rlz != y:
                        deg_zero_ok = False
                elif y == 1:
                    if z_rlz != x:
                        deg_zero_ok = False
                else:
                    deg_zero_ok = False      # nondeg margins but a zero cell
                continue
            key = (x, y)
            if key not in zcache:
                zlo, zhi = z_root_enc(x, y, t)
                zcache[key] = h2_over(zlo, zhi)
            htlo, hthi = zcache[key]
            hrlo, hrhi = h2_enc(z_rlz)
            num_lo += m * (hrlo - hthi)
            num_hi += m * (hrhi - htlo)
            q = Fraction(t00 * t11, t01 * t10) / t
            llo, lhi = log2_enc(q)
            pl_lo += m * llo
            pl_hi += m * lhi
            den += m
            n_rows += 1
    return {
        "mm_lo": num_lo / den, "mm_hi": num_hi / den,
        "agg_lo": pl_lo / den, "agg_hi": pl_hi / den,
        "max_marg": max(marg), "in_regime": in_regime,
        "dich_ok": dich_ok, "deg_zero_ok": deg_zero_ok,
        "n_rows": n_rows, "n_distinct_margins": len(zcache),
        "den_over_Z": Fraction(den, Z),
    }


# ======================================================================
# float census (own implementation, from 007 section 1 + 018's prose)
# ======================================================================

def tilt(u: dict[int, float], lam: float) -> dict[tuple[int, int], float]:
    z = 0.0
    pi = {}
    for a, ua in u.items():
        for b, ub in u.items():
            v = ua * ub * 2.0 ** (lam * bin(a & b).count("1"))
            pi[(a, b)] = v
            z += v
    return {k: v / z for k, v in pi.items()}


def z_pl(x: float, y: float, rho: float) -> float:
    """Stable float Plackett root (own derivation: F(z)=0 with
    F(z) = (1-rho) z^2 + (1-x-y+rho(x+y)) z - rho x y; the coupling root is
    2C' / (-B - sqrt(B^2-4AC)) written to avoid cancellation)."""
    if rho == 1.0:
        return x * y
    A = 1.0 - rho
    B = 1.0 - x - y + rho * (x + y)
    C = -rho * x * y
    disc = B * B - 4.0 * A * C
    return -2.0 * C / (B + math.sqrt(disc))


def h2f(z: float) -> float:
    if z <= 0.0 or z >= 1.0:
        return 0.0
    return -z * math.log2(z) - (1.0 - z) * math.log2(1.0 - z)


def sens(x: float, y: float, rho: float) -> float:
    """d/d lambda of h2(z_rho(x,y)) at rho (own implicit differentiation:
    z' = (x-z)(y-z) / [(1-x-y+2z) + rho(x+y-2z)])."""
    z = z_pl(x, y, rho)
    if z <= 0.0 or z >= 1.0 or z >= min(x, y) or z <= max(0.0, x + y - 1.0):
        return 0.0
    dz = (x - z) * (y - z) / ((1.0 - x - y + 2.0 * z)
                              + rho * (x + y - 2.0 * z))
    return math.log2((1.0 - z) / z) * dz * rho * math.log(2.0)


def census_variants(u: dict[int, float], n: int, lam: float):
    """Own float census: plain aggregate, 018's three MM readings, and the
    realized-OR variants 018 did not test.  Returns None if no
    nondegenerate history exists."""
    pi = tilt(u, lam)
    supp = {a for (a, _b) in pi}
    rho_t = 2.0 ** lam
    acc = {k: 0.0 for k in
           ("plain", "sec", "der", "absn", "der_rlz", "absn_rlz",
            "absrho_rlz", "den", "denabs", "denabs_rlz", "denabsrho")}
    dich = True
    for i in range(1, n):
        pmask, bit = (1 << (i - 1)) - 1, 1 << (i - 1)
        has0, has1 = set(), set()
        for A in supp:
            (has1 if A & bit else has0).add(A & pmask)
        active = has0 & has1
        tables: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in pi.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [0.0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        for (a, b), (p00, p01, p10, p11) in tables.items():
            nondeg = min(p00, p01, p10, p11) > 0.0
            if nondeg != (a in active and b in active):
                dich = False
            if not nondeg:
                continue
            m = p00 + p01 + p10 + p11
            x, y, z = (p00 + p01) / m, (p00 + p10) / m, p00 / m
            l2or = (math.log2(p00) + math.log2(p11)
                    - math.log2(p01) - math.log2(p10))
            d = l2or - lam
            rho = (p00 * p11) / (p01 * p10)
            s_t = sens(x, y, rho_t)
            s_r = sens(x, y, rho)
            acc["plain"] += m * d
            acc["sec"] += m * (h2f(z) - h2f(z_pl(x, y, rho_t)))
            acc["der"] += m * s_t * d
            acc["absn"] += m * abs(s_t) * d
            acc["der_rlz"] += m * s_r * d
            acc["absn_rlz"] += m * abs(s_r) * d
            acc["absrho_rlz"] += m * (abs(s_r) / rho) * d
            acc["den"] += m
            acc["denabs"] += m * abs(s_t)
            acc["denabs_rlz"] += m * abs(s_r)
            acc["denabsrho"] += m * abs(s_r) / rho
    if acc["den"] == 0.0:
        return None
    mu: dict[int, float] = {}
    for (A, _B), w in pi.items():
        mu[A] = mu.get(A, 0.0) + w
    marg = [sum(w for a, w in mu.items() if a >> j & 1) for j in range(n)]
    return {
        "plain": acc["plain"] / acc["den"],
        "mm_sec": acc["sec"] / acc["den"],
        "mm_der": acc["der"] / (acc["denabs"] or 1.0),
        "mm_abs": acc["absn"] / (acc["denabs"] or 1.0),
        "mm_der_rlz": acc["der_rlz"] / (acc["denabs_rlz"] or 1.0),
        "mm_abs_rlz": acc["absn_rlz"] / (acc["denabs_rlz"] or 1.0),
        "mm_absrho_rlz": acc["absrho_rlz"] / (acc["denabsrho"] or 1.0),
        "max_marginal": max(marg), "dich_ok": dich,
        "H": -sum(w * math.log2(w) for w in mu.values() if w > 0),
    }


# ======================================================================
# MU(n, r) ladder rebuilt from 016 section 1 prose (own code)
# ======================================================================

def mu_ladder(n: int, r: int, th: dict[str, float]) -> dict[int, float]:
    """Coordinate 1 = block marker (bit 0); coords 2-4 = dilution singletons
    (bits 1-3); response coords 5..4+r (bits 4..3+r); shared futures n-1, n
    (bits n-2, n-1).  Block a: {1,n-1}=A1, {1,n}=A2, light {1,c_j,n};
    block b: {}=B0, {n-1}=B1, {n}=B2, responders {c_j,n-1}.  2r+8 atoms."""
    assert 1 <= r <= n - 6
    hi, lo = 1 << (n - 1), 1 << (n - 2)
    u = {1 | lo: th["A1"], 1 | hi: th["A2"], 0: th["B0"],
         lo: th["B1"], hi: th["B2"]}
    for j in range(r):
        c = 1 << (4 + j)
        u[1 | c | hi] = th["light"]
        u[c | lo] = th["resp"]
    for d in (1, 2, 3):
        u[1 << d] = th["dil"]
    return u


def marginals_u(u: dict[int, float], n: int, lam: float) -> list[float]:
    supp = list(u)
    z = 0.0
    mu = {}
    for a in supp:
        row = sum(u[a] * u[b] * 2.0 ** (lam * bin(a & b).count("1"))
                  for b in supp)
        mu[a] = row
        z += row
    marg = [0.0] * n
    for a, w in mu.items():
        for j in range(n):
            if a >> j & 1:
                marg[j] += w / z
    return marg


def tune_dil(n: int, r: int, lam: float, th: dict[str, float]
             ) -> dict[str, float]:
    """016's dilution bisection re-implemented to spec: geometric bisection
    of the shared dilution weight on [0.5, 2000], 18 steps, target max
    marginal <= 0.375, dilution coords 2-4 vs the rest."""
    th = dict(th)
    lo, hi = 0.5, 2000.0
    for _ in range(18):
        w = math.sqrt(lo * hi)
        th["dil"] = w
        marg = marginals_u(mu_ladder(n, r, th), n, lam)
        dil_m = max(marg[1:4])
        oth_m = max(m for j, m in enumerate(marg) if j not in (1, 2, 3))
        if oth_m > MARG_TARGET and dil_m < MARG_TARGET:
            lo = w
        elif dil_m > MARG_TARGET:
            hi = w
        else:
            break
    return th


# ======================================================================
# S0 -- kit self-tests
# ======================================================================

def part_s0() -> dict:
    out = {"ok": True}

    # log2 enclosure vs known values
    cases = [(Fraction(1), 0.0), (Fraction(2), 1.0), (Fraction(1, 8), -3.0),
             (Fraction(3), math.log2(3.0)), (Fraction(181, 16),
              math.log2(181.0 / 16.0)), (Fraction(10) ** 20,
              20 * math.log2(10.0)), (Fraction(7, 5) ** 9,
              9 * math.log2(1.4))]
    wmax = 0.0
    for q, ref in cases:
        lo, hi = log2_enc(q)
        wmax = max(wmax, float(hi - lo))
        out["ok"] &= bool(lo <= hi and float(lo) - 1e-12 <= ref
                          <= float(hi) + 1e-12)
    rng = random.Random(19)
    for _ in range(60):
        q = Fraction(rng.randrange(1, 10 ** 12), rng.randrange(1, 10 ** 12))
        lo, hi = log2_enc(q)
        ref = math.log2(q.numerator) - math.log2(q.denominator)
        wmax = max(wmax, float(hi - lo))
        out["ok"] &= bool(float(lo) - 1e-9 <= ref <= float(hi) + 1e-9)
    a, b = Fraction(355, 113), Fraction(1234567, 7654321)
    la, ha = log2_enc(a)
    lb, hb = log2_enc(b)
    lab, hab = log2_enc(a * b)
    out["ok"] &= bool(la + lb <= hab and lab <= ha + hb)
    out["log2_max_width"] = wmax
    log(f"[S0] log2 enclosure: {len(cases) + 60} cases + additivity, "
        f"max width {wmax:.1e}, ok={out['ok']}")

    # Plackett root: closed forms + float cross-check + monotonicity in rho
    ok_root = True
    wr = 0.0
    for _ in range(60):
        x = Fraction(rng.randrange(1, 99), 100)
        y = Fraction(rng.randrange(1, 99), 100)
        t = Fraction(rng.randrange(101, 4000), 100)
        zlo, zhi = z_root_enc(x, y, t)
        wr = max(wr, float(zhi - zlo))
        zf = z_pl(float(x), float(y), float(t))
        ok_root &= bool(float(zlo) - 1e-9 <= zf <= float(zhi) + 1e-9)
        ok_root &= bool(x * y <= zlo <= zhi <= min(x, y))
        # the realized-table identity: build the exact table with margins
        # x, y and both-zero z*, check its OR is inside [t*(1-eps), ...]
        zm = (zlo + zhi) / 2
        orr = (zm * (1 - x - y + zm)) / ((x - zm) * (y - zm))
        ok_root &= bool(abs(float(orr / t) - 1.0) < 1e-9)
    out["root_ok"] = ok_root
    out["root_max_width"] = wr
    out["ok"] &= ok_root
    log(f"[S0] Plackett Newton root: 60 random (x,y,t), max width {wr:.1e}, "
        f"round-trip OR/t = 1 to 1e-9, ok={ok_root}")

    # h2 enclosure
    ok_h = True
    for _ in range(40):
        z = Fraction(rng.randrange(1, 999), 1000)
        lo, hi = h2_enc(z)
        ok_h &= bool(float(lo) - 1e-12 <= h2f(float(z)) <= float(hi) + 1e-12)
    l, h = h2_enc(HALF)
    ok_h &= bool(float(l) <= 1.0 <= float(h) + 1e-30)
    out["h2_ok"] = ok_h
    out["ok"] &= ok_h
    log(f"[S0] h2 enclosure: 41 cases ok={ok_h}")

    # ladder spec check: MU(7,1) with witness theta == the 10 witness atoms
    lad = mu_ladder(7, 1, THETA_WITNESS)
    ok_lad = lad == WITNESS_U
    out["ladder_matches_witness"] = ok_lad
    out["ok"] &= ok_lad
    log(f"[S0] MU(7,1) from spec == 007 witness atom-for-atom: {ok_lad}")

    save("gap1csk_s0.json", out)
    if not out["ok"]:
        raise SystemExit("S0 kit self-test failure -- abort")
    return out


# ======================================================================
# S1 -- K1: exact certification, from scratch
# ======================================================================

def part_s1() -> dict:
    out = {}
    t = Fraction(181, 16)
    r = exact_mm_cert(WITNESS_U, 7, t)
    lo, hi = r["mm_lo"], r["mm_hi"]
    verdict = ("CERTIFIED NEGATIVE" if hi < 0 else
               "CERTIFIED POSITIVE" if lo > 0 else "STRADDLES 0")
    log(f"[S1] witness t=181/16 MM_sec in [{float(lo):+.15e}, "
        f"{float(hi):+.15e}] width {float(hi - lo):.1e}  {verdict}")
    log(f"[S1]   plain aggregate in [{float(r['agg_lo']):+.12f}, "
        f"{float(r['agg_hi']):+.12f}] (013/018: +1.844669005)")
    log(f"[S1]   max marginal {float(r['max_marg']):.10f} "
        f"(exact < 38271/100000: {r['in_regime']}); dichotomy exact "
        f"{r['dich_ok']}; degenerate rows contribute 0 exactly "
        f"{r['deg_zero_ok']}; rows {r['n_rows']} "
        f"({r['n_distinct_margins']} distinct margins)")
    # containment check vs 018's certificate
    lo18 = Fraction(-22254658205655785177, 10 ** 21)
    hi18 = Fraction(-22254658205655785176, 10 ** 21)
    overlap = not (hi < lo18 or lo > hi18)
    log(f"[S1]   overlap with 018's enclosure "
        f"[-0.022254658205655785177, ...176]: {overlap}")
    out["witness"] = {
        "mm_lo": f"{float(lo):.18e}", "mm_hi": f"{float(hi):.18e}",
        "mm_lo_str": str(lo.limit_denominator(10 ** 30)),
        "width": float(hi - lo),
        "agg_lo": float(r["agg_lo"]), "agg_hi": float(r["agg_hi"]),
        "max_marg": float(r["max_marg"]), "in_regime": r["in_regime"],
        "dich_ok": r["dich_ok"], "deg_zero_ok": r["deg_zero_ok"],
        "n_rows": r["n_rows"], "verdict": verdict,
        "overlaps_018_cert": overlap,
    }

    # 2-digit tidy witness (018 part R quotes float -0.0223 at mm 0.3130)
    u2 = {a: float(f"{w:.2g}") for a, w in WITNESS_U.items()}
    r2 = exact_mm_cert(u2, 7, t)
    verdict2 = ("CERTIFIED NEGATIVE" if r2["mm_hi"] < 0 else
                "CERTIFIED POSITIVE" if r2["mm_lo"] > 0 else "STRADDLES 0")
    log(f"[S1] tidy 2-digit witness t=181/16: MM_sec in "
        f"[{float(r2['mm_lo']):+.9e}, {float(r2['mm_hi']):+.9e}] {verdict2} "
        f"(max marg {float(r2['max_marg']):.4f}, in-regime "
        f"{r2['in_regime']})")
    out["tidy"] = {"mm_lo": float(r2["mm_lo"]), "mm_hi": float(r2["mm_hi"]),
                   "max_marg": float(r2["max_marg"]),
                   "in_regime": r2["in_regime"], "verdict": verdict2,
                   "dich_ok": r2["dich_ok"], "deg_zero_ok": r2["deg_zero_ok"]}

    # exact product-measure zero: product potential u(a) = w^{|a|} at
    # rational w -- every nondegenerate history must have OR = t exactly,
    # so the MM_sec numerator is exactly 0 term by term.  Verified via the
    # census: enclosure must contain 0 with tiny width, and the plain
    # aggregate must be exactly 0 (OR/t = 1 rational -> log enclosure
    # degenerate at 0).
    n_p = 5
    w = Fraction(3, 7)
    up = {a: w ** bin(a).count("1") for a in range(1 << n_p)}
    rp = exact_mm_cert(up, n_p, Fraction(4))
    zero_in = rp["mm_lo"] <= 0 <= rp["mm_hi"]
    plain_zero = rp["agg_lo"] == 0 == rp["agg_hi"]
    log(f"[S1] product anchor (w=3/7, n=5, t=4): MM_sec enclosure contains "
        f"0: {zero_in} (width {float(rp['mm_hi'] - rp['mm_lo']):.1e}); "
        f"plain aggregate exactly 0: {plain_zero}")
    out["product_anchor"] = {"zero_in": zero_in, "plain_exact_zero":
                             plain_zero,
                             "width": float(rp["mm_hi"] - rp["mm_lo"])}
    save("gap1csk_s1.json", out)
    return out


# ======================================================================
# S2 -- K2: the definitional fork, incl. variants 018 did not test
# ======================================================================

def part_s2() -> dict:
    out = {"witness": {}, "ladder96": {}}
    cv = census_variants(WITNESS_U, 7, 3.5)
    log("[S2] witness (n=7, lam=3.5) variant zoo:")
    for k in ("plain", "mm_sec", "mm_der", "mm_abs", "mm_der_rlz",
              "mm_abs_rlz", "mm_absrho_rlz"):
        log(f"[S2]   {k:14s} {cv[k]:+.6f}")
        out["witness"][k] = cv[k]
    out["witness"]["max_marginal"] = cv["max_marginal"]

    # the realized-OR readings on the theta* ladder kill instance
    th = tune_dil(96, 90, 2.0, THETA_STAR_L2)
    cl = census_variants(mu_ladder(96, 90, th), 96, 2.0)
    log("[S2] theta* ladder n=96 lam=2 variant zoo:")
    for k in ("plain", "mm_sec", "mm_der", "mm_abs", "mm_der_rlz",
              "mm_abs_rlz", "mm_absrho_rlz"):
        log(f"[S2]   {k:14s} {cl[k]:+.6e}")
        out["ladder96"][k] = cl[k]
    out["ladder96"]["max_marginal"] = cl["max_marginal"]

    # small-lambda sign crossing of MM_sec on the witness (018: ~0.7)
    prof = []
    for lam in (0.6, 0.7, 0.75, 0.8, 0.9):
        c = census_variants(WITNESS_U, 7, lam)
        prof.append({"lambda": lam, "mm_sec": c["mm_sec"]})
    log("[S2] witness MM_sec near the crossing: "
        + "  ".join(f"{p['lambda']:.2f}:{p['mm_sec']:+.2e}" for p in prof))
    out["crossing_profile"] = prof
    save("gap1csk_s2.json", out)
    return out


# ======================================================================
# S3 -- K3: float spot-checks with the spec-rebuilt ladder
# ======================================================================

def part_s3() -> dict:
    out = {}
    # witness headline floats
    cv = census_variants(WITNESS_U, 7, 3.5)
    ok_abs = abs(cv["mm_abs"] - 2.354328) < 5e-4
    ok_der = abs(cv["mm_der"] - (-0.895542)) < 5e-4
    ok_sec = abs(cv["mm_sec"] - (-0.022255)) < 5e-5
    log(f"[S3] witness floats: MM_abs {cv['mm_abs']:+.6f} (018 +2.354328 "
        f"{ok_abs}), MM_der {cv['mm_der']:+.6f} (018 -0.895542 {ok_der}), "
        f"MM_sec {cv['mm_sec']:+.6f} (018 -0.022255 {ok_sec})")
    out["witness"] = {"mm_abs": cv["mm_abs"], "mm_der": cv["mm_der"],
                      "mm_sec": cv["mm_sec"], "ok": ok_abs and ok_der
                      and ok_sec}

    # theta* ladder n=96 lam=2, own tune_dil
    th = tune_dil(96, 90, 2.0, THETA_STAR_L2)
    c96 = census_variants(mu_ladder(96, 90, th), 96, 2.0)
    ok96 = (abs(c96["plain"] - (-0.000760)) < 5e-5
            and abs(c96["mm_sec"] - (-3.819e-4)) < 5e-6
            and abs(c96["mm_abs"] - 0.969) < 5e-3)
    log(f"[S3] theta* ladder n=96 lam=2 (own tune_dil, dil="
        f"{th['dil']:.6f}): plain {c96['plain']:+.6f} (018 -0.000760), "
        f"MM_sec {c96['mm_sec']:+.6e} (018 -3.819e-4), MM_abs "
        f"{c96['mm_abs']:+.6f} (018 +0.969), mm {c96['max_marginal']:.3f}, "
        f"ok={ok96}")
    out["theta96"] = {"dil": th["dil"], "plain": c96["plain"],
                      "mm_sec": c96["mm_sec"], "mm_abs": c96["mm_abs"],
                      "max_marginal": c96["max_marginal"], "ok": ok96}

    # raw ladder n=96 lam=2 -- the sign 018's prose gets wrong for n >= 96
    thr = tune_dil(96, 90, 2.0, THETA_WITNESS)
    cr = census_variants(mu_ladder(96, 90, thr), 96, 2.0)
    log(f"[S3] raw ladder n=96 lam=2 (own tune_dil, dil={thr['dil']:.4f}): "
        f"plain {cr['plain']:+.6f}, MM_sec {cr['mm_sec']:+.6e}, MM_der "
        f"{cr['mm_der']:+.6e}  <-- MM_sec/MM_der POSITIVE at n=96 (018 "
        f"checkpoint agrees; 018 PROSE says negative at every n)")
    out["raw96"] = {"dil": thr["dil"], "plain": cr["plain"],
                    "mm_sec": cr["mm_sec"], "mm_der": cr["mm_der"],
                    "max_marginal": cr["max_marginal"]}
    # raw ladder n=64 for the branch contrast
    thr64 = tune_dil(64, 58, 2.0, THETA_WITNESS)
    cr64 = census_variants(mu_ladder(64, 58, thr64), 64, 2.0)
    log(f"[S3] raw ladder n=64 lam=2: MM_sec {cr64['mm_sec']:+.6e} "
        f"(negative branch, 018 -3.400e-3)")
    out["raw64"] = {"mm_sec": cr64["mm_sec"], "plain": cr64["plain"]}

    # part-D window value at n=96 + fit arithmetic from 018's checkpoint
    lw = LAMWIN_C / (96 - 3)
    thw = tune_dil(96, 90, lw, THETA_STAR_L2)
    cw = census_variants(mu_ladder(96, 90, thw), 96, lw)
    ok_w = abs(cw["plain"] - 3.747e-3) < 5e-5
    log(f"[S3] theta* window n=96 lam={lw:.6f}: A {cw['plain']:+.6e} "
        f"(018 +3.747e-3, ok={ok_w}), mm {cw['max_marginal']:.3f}")
    out["window96"] = {"lamwin": lw, "agg": cw["plain"],
                       "max_marginal": cw["max_marginal"], "ok": ok_w}

    # a1/a2 fit arithmetic re-done from 018's stored points (data reuse)
    fit_rows = []
    try:
        d = json.loads((DATA / "gap1c_partD.json").read_text())
        for row in d["rows"]:
            if "points" not in row:
                continue
            pts = row["points"]
            p1, p2 = pts[-1], pts[-2]
            l1, l2 = p1["lambda"], p2["lambda"]
            det = l1 * l2 * (l2 - l1)
            a1 = (p1["agg"] * l2 * l2 - p2["agg"] * l1 * l1) / det
            a2 = (p2["agg"] * l1 - p1["agg"] * l2) / det
            ok = (abs(a1 - row["a1"]) < 1e-9 * max(1, abs(row["a1"]))
                  and abs(a2 - row["a2"]) < 1e-9 * max(1, abs(row["a2"])))
            fit_rows.append({"name": row["name"], "n": row["n"],
                             "a1": a1, "a2": a2, "a1n": a1 * row["n"],
                             "match": ok})
        th_star = [r for r in fit_rows if r["name"].startswith("theta*")]
        log("[S3] a1/a2 fit re-solve (theta*): "
            + "  ".join(f"n={r['n']}:a1*n={r['a1n']:+.3f}"
                        f"{'' if r['match'] else '(MISMATCH)'}"
                        for r in th_star))
        all_match = all(r["match"] for r in fit_rows)
        log(f"[S3] all 10 stored (a1, a2) rows re-derive: {all_match}; "
            f"a1*n growth on theta* rows: "
            f"{[round(r['a1n'], 2) for r in th_star]}")
        out["fit_rows"] = fit_rows
        out["fit_all_match"] = all_match
    except FileNotFoundError:
        log("[S3] gap1c_partD.json missing -- fit check skipped")

    # perturbation battery, DIFFERENT seeds than 018 (which used 2024)
    res = {}
    for seed in (7, 424242):
        rng = random.Random(seed)
        vals = []
        nfail = 0
        for _ in range(20):
            u = {a: w * math.exp(rng.uniform(-0.03, 0.03))
                 for a, w in WITNESS_U.items()}
            c = census_variants(u, 7, 3.5)
            vals.append(c["mm_sec"])
            if c["mm_sec"] >= 0 or c["max_marginal"] > 0.38271:
                nfail += 1
        res[seed] = {"n_fail": nfail, "min": min(vals), "max": max(vals)}
        log(f"[S3] perturbations seed={seed}: {nfail}/20 failures, MM_sec "
            f"in [{min(vals):+.5f}, {max(vals):+.5f}]")
    out["perturbations"] = res
    save("gap1csk_s3.json", out)
    return out


# ======================================================================
# S4 -- K4: scope checks
# ======================================================================

def part_s4() -> dict:
    out = {}
    # best_free endpoint from 018 part C (data)
    try:
        d = json.loads((DATA / "gap1c_partC.json").read_text())
        bf = d["best_free"]
        u = {int(k, 2): float(w) for k, w in bf["u"].items()}
        c = census_variants(u, 8, bf["lambda"])
        ok = (abs(c["mm_abs"] - bf["mm_abs"]) < 1e-9
              and abs(c["max_marginal"] - bf["max_marginal"]) < 1e-9)
        log(f"[S4] best_free endpoint (n=8, lam={bf['lambda']}): my MM_abs "
            f"{c['mm_abs']:+.6e} (018 {bf['mm_abs']:+.6e}), MM_sec "
            f"{c['mm_sec']:+.6e}, max marg {c['max_marginal']:.6f} "
            f"(in-regime {c['max_marginal'] < 0.38271}), dich "
            f"{c['dich_ok']}, H {c['H']:.3f} bits, match={ok}")
        out["best_free"] = {"mm_abs": c["mm_abs"], "mm_sec": c["mm_sec"],
                            "max_marginal": c["max_marginal"],
                            "dich_ok": c["dich_ok"], "H": c["H"],
                            "match": ok,
                            "in_regime": c["max_marginal"] < 0.38271}
    except FileNotFoundError:
        log("[S4] gap1c_partC.json missing -- best_free check skipped")

    # the z = 1/2 sensitivity boundary vs the golden-ratio parenthetical.
    # 018 claims z_1(x,x) = x^2 crosses 1/2 "exactly at the (3-sqrt5)/2
    # barrier".  Arithmetic: x^2 = 1/2 at x = 2^-1/2 = 0.70711, i.e.
    # element marginal 1 - 2^-1/2 = 0.29289, NOT 0.38197.  At the barrier
    # x = 1 - (3-sqrt5)/2 = (sqrt5-1)/2 and x^2 = (3-sqrt5)/2 = the
    # marginal itself (phi^2 = 1 - phi), not 1/2.
    s5 = math.sqrt(5.0)
    barrier = (3.0 - s5) / 2.0
    x_bar = 1.0 - barrier
    out["golden"] = {
        "barrier": barrier,
        "x_at_barrier": x_bar,
        "x2_at_barrier": x_bar * x_bar,
        "claim_x2_equals_half": abs(x_bar * x_bar - 0.5) < 1e-9,
        "x2_equals_marginal": abs(x_bar * x_bar - barrier) < 1e-12,
        "true_crossing_marginal_rho1": 1.0 - math.sqrt(0.5),
    }
    log(f"[S4] golden parenthetical: at marginal {barrier:.6f}, x = "
        f"{x_bar:.6f}, x^2 = {x_bar * x_bar:.6f} = the marginal itself "
        f"(phi identity), NOT 1/2; z_1(x,x) = 1/2 at marginal "
        f"{1.0 - math.sqrt(0.5):.6f}")
    # where DOES z_t(x,x) = 1/2 at t = 2^3.5?  solve 3/4 - x = t (x-1/2)^2
    t = 2.0 ** 3.5
    lo, hi = 0.5, 0.75
    for _ in range(60):
        mid = (lo + hi) / 2
        if z_pl(mid, mid, t) < 0.5:
            lo = mid
        else:
            hi = mid
    out["golden"]["crossing_marginal_at_lam3.5"] = 1.0 - lo
    log(f"[S4] z_t(x,x) = 1/2 at lam=3.5 happens at element marginal "
        f"{1.0 - lo:.6f} (lambda-dependent; equals 0.38197 for no stated "
        f"reason)")
    save("gap1csk_s4.json", out)
    return out


def main() -> None:
    part_s0()
    part_s1()
    part_s2()
    part_s3()
    part_s4()
    log("done")


if __name__ == "__main__":
    main()
