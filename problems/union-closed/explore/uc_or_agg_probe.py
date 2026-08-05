#!/usr/bin/env python3
"""Attempt 014: large-n certified probe of the i-AGGREGATED odds-ratio control.

Object under test (007 SS5 "what survives" item 1, exactly as 013 part C
certified it at n = 7): for the tilt coupling pi_lambda(A,B) prop
u(A) u(B) 2^{lambda |A&B|} with both marginals mu, and M_i the well-posed
per-coordinate average of 007 SS1 (mass-weighted mean of log2 OR_i over
N_i x N_i histories; N_i = reachable (i-1)-bit prefixes with both
continuations reachable in supp(mu); i = n excluded since M_n = lambda
identically),

    A(mu, lambda) = sum_{i<n, N_i != empty} w_i (M_i - lambda) / sum w_i  >= 0
    with w_i = pi-mass of N_i x N_i (the "defined mass"),

for all mu in the record regime (every elementwise marginal < 0.38271) and
all lambda >= 0?  Before this attempt the aggregate was certified positive
only at n <= 7 (013 part C: +1.844669 on the 10-atom witness).  Per the
standing warning from 012 (small-n flattening reversed past n ~ 22) this
probe extends the census to n in {10, 14, 20, 24, 28, 32} BEFORE any proof
effort.

Sign convention note: under the conjecture-friendliest alternative
bookkeeping (score every margin-degenerate history at exactly lambda, 013's
"surrogate-lambda" convention) the aggregate's NUMERATOR is unchanged and
its denominator only grows, so the SIGN of A is the same under either
convention.  A sign probe of the conditioning-out aggregate is therefore
convention-robust for free.

Engines:
  * sparse engine -- direct atom-pair census from a potential (005's trick:
    any positive u is its own marginal's Sinkhorn potential).  Cost
    O(m^2 n) for m atoms, so it reaches ANY n on sparse supports.  Fresh
    code, conventions copied from 007 SS1 / 013 (cross-checked in part X).
  * orbit engine -- for measures invariant under permutations of coords
    {3..n} (all the structured families: Bernoulli/Sawin mixtures, slices,
    delta_0+Bern, crash mixtures).  Histories are collapsed to orbit
    classes ((a1,a2) special-coordinate pattern, tail prefix weights,
    overlap), the same S_{n-2} symmetry 012's budget engine used to reach
    n = 32.  Cells per class are exact orbit sums via tail pair-kernels
    PK[s][l] = C(f,s) * sum_j C(s,j) C(f-s,l-j) t^j.
  * exact rational certifier -- both engines re-run end-to-end in
    fractions.Fraction at rational tilt t, with every log2 replaced by a
    certified rational enclosure (atanh series, directed rounding, explicit
    tail bound; own implementation following the 2026-07-31 lab standard:
    float agreement is NOT certification).  Potentials are dyadic
    rationals (Fraction(float)), so all table entries stay dyadic and no
    denominator blow-up occurs.

Parts (run `python uc_or_agg_probe.py X O S H E` or any subset; default all):
  X  cross-checks: (X1) sparse engine reproduces 007's witness numbers
     (M_5 - lambda = -0.122033, aggregate +1.8447) and 013's exact
     aggregate enclosure at t = 181/16; (X2) orbit engine vs direct
     sparse+Sinkhorn on full supports at n = 6, 7 (slice, Sawin mixture,
     eps-crash, d0+Bern); (X3) exact orbit rows vs float orbit rows.
  O  orbit family census over n in {10,14,20,24,28,32} x lambda grid
     (incl. the 011 window law 4.847/(n-3)): Bernoulli cells, Sawin
     geometric mixtures, delta0+Bern, Chase-Lovett slice + smoothed slice,
     eps-mixed crash, 012 crash-mixture cells (p,d).
  S  sparse structured census at the same n: pure crash / mirror, the
     10-atom witness genre re-laid-out at each n with its dilution-weight
     sweep (tightest in-regime point per n), multi-gadget stacks (K
     light-slice gadgets, 007 lead 1 generalized).
  H  seeded adversarial hill-climbs on the AGGREGATE objective (marginal-
     penalized, dynamic-range-clamped): free sparse supports, witness-
     seeded, multi-gadget-parametrized, at n = 10, 14, 20 (+ spot 24).
  E  exact rational certification of the tightest points found (per-n
     minimum and the global minimum), plus exact in-regime marginals.

Standard library only.  Deterministic (fixed seeds).  Runtime: X+S+H ~ 15
min, O ~ 40-80 min (n = 28, 32 dominate), E ~ 10-30 min.
Checkpoints: problems/union-closed/data/or_agg_probe_part[XOSHE]*.json
Log:         problems/union-closed/data/or_agg_probe_run.log  (tee by hand)

Usage:
    python problems/union-closed/explore/uc_or_agg_probe.py [parts] 2>&1 |
        tee -a problems/union-closed/data/or_agg_probe_run.log
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

RECORD = 0.38271
F0 = Fraction(0)
NS_GRID = (10, 14, 20, 24, 28, 32)


def popc(x: int) -> int:
    return bin(x).count("1")


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(
        json.dumps(obj, indent=1, sort_keys=True, default=str) + "\n",
        encoding="utf-8")


def lam_grid(n: int) -> list[float]:
    """lambda grid per n: tilt-recipe cells log2(1.06), log2(1.20) from
    008/012's calibration, O(1) values, and the 011 window law
    lam_max ~ 4.847/(n-3) (half / at / twice)."""
    base = [math.log2(1.06), math.log2(1.20), 0.5, 1.0, 2.0, 3.5, 5.0]
    win = [2.4 / (n - 3), 4.847 / (n - 3), 9.7 / (n - 3)]
    return sorted(set(round(v, 6) for v in base + win))


# ======================================================================
# certified log2 enclosures (rational in, rational bounds out).
# Own implementation of the 013 standard: ln s = 2 atanh((s-1)/(s+1)),
# partial sum is a lower bound (all terms positive), tail after J terms
# <= 2 z^{2J+1} / ((2J+1)(1-z^2)).  Inputs are directed-rounded to
# bounded precision first so Fraction sizes stay tame.
# ======================================================================

def _ln_1_to_2(s: Fraction, terms: int) -> tuple[Fraction, Fraction]:
    if s == 1:
        return F0, F0
    z = (s - 1) / (s + 1)
    z2 = z * z
    term = z
    acc = F0
    for j in range(terms):
        acc += term / (2 * j + 1)
        term *= z2
    lo = 2 * acc
    tail = 2 * term / ((2 * terms + 1) * (1 - z2))
    return lo, lo + tail


LN2_LO, LN2_HI = _ln_1_to_2(Fraction(2), 48)   # width ~ 1e-47


def _dround(x: Fraction, up: bool, bits: int = 120) -> Fraction:
    """Directed rounding of x to a dyadic rational with a <= bits-bit
    numerator: result <= x if not up, result >= x if up.  Keeping every
    intermediate dyadic avoids the gcd churn that makes long certified
    sums explode."""
    if x == 0:
        return F0
    n, d = x.numerator, x.denominator
    shift = bits - (abs(n).bit_length() - d.bit_length())
    if shift >= 0:
        num, den = n << shift, d
    else:
        num, den = n, d << (-shift)
    q = num // den                       # floor (rounds down, also for <0)
    if up and q * den != num:
        q += 1
    if shift >= 0:
        return Fraction(q, 1 << shift)
    return Fraction(q << (-shift))


def _ln_1_to_2_fp(s: Fraction, up: bool, terms: int,
                  bits: int = 120) -> Fraction:
    """Directed bound on ln s for 1 <= s < 2, dyadic fixed-point: every
    operation is rounded in the safe direction, and the atanh-series tail
    2 z^{2J+1}/((2J+1)(1-z^2)) is added on the upper side."""
    if s == 1:
        return F0
    z = _dround((s - 1) / (s + 1), up, bits)
    z2 = _dround(z * z, up, bits)
    term = z
    acc = F0
    for j in range(terms):
        acc = _dround(acc + term / (2 * j + 1), up, bits)
        term = _dround(term * z2, up, bits)
    val = 2 * acc
    if up:
        tail = _dround(2 * term / ((2 * terms + 1) * (1 - z2)), True, bits)
        val = _dround(val + tail, True, bits)
    return val


LN2D_LO = _dround(LN2_LO, False)
LN2D_HI = _dround(LN2_HI, True)


def log2_enclosure(r: Fraction, terms: int = 12,
                   bits: int = 120) -> tuple[Fraction, Fraction]:
    """[lo, hi] with lo <= log2 r <= hi, dyadic rationals; r > 0.
    terms=12 gives width ~ 1e-13 (z <= 1/3 after range reduction)."""
    assert r > 0
    out = []
    for up in (False, True):
        rr = _dround(r, up, bits)
        if rr <= 0:
            rr = r
        k = 0
        while rr >= 2:
            rr /= 2
            k += 1
        while rr < 1:
            rr *= 2
            k -= 1
        ln_dir = _ln_1_to_2_fp(rr, up, terms, bits)
        if up:
            out.append(_dround(k + ln_dir / LN2D_LO, True, bits)
                       if ln_dir >= 0 else
                       _dround(k + ln_dir / LN2D_HI, True, bits))
        else:
            out.append(_dround(k + ln_dir / LN2D_HI, False, bits)
                       if ln_dir >= 0 else
                       _dround(k + ln_dir / LN2D_LO, False, bits))
    return out[0], out[1]


# ======================================================================
# sparse engine (float): census + aggregate from a potential
# ======================================================================

def sparse_rows(u: dict[int, float], lam: float, n: int):
    """Per i = 1..n: dict(rows=[(mass, log2 OR)], def_mass, n_active,
    cell_floor).  Masses unnormalized (one common Z for all i, which
    cancels in every ratio reported).  N_i membership is structural
    (from supp(u)); a structurally nondegenerate history with a zero
    float cell is an underflow and is flagged."""
    t = 2.0 ** lam
    supp = sorted(u)
    kern = {}
    for a in supp:
        ua = u[a]
        for b in supp:
            kern[(a, b)] = ua * u[b] * t ** popc(a & b)
    out = []
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        prefixes = {A & pmask for A in supp}
        active = {c for c in prefixes
                  if any(A & pmask == c and A & bit for A in supp)
                  and any(A & pmask == c and not A & bit for A in supp)}
        tables: dict[tuple[int, int], list[float]] = {}
        for (A, B), w in kern.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [0.0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        rows = []
        dm = 0.0
        floor = None
        underflow = False
        for (a, b), (p00, p01, p10, p11) in tables.items():
            if a in active and b in active:
                if min(p00, p01, p10, p11) <= 0.0:
                    underflow = True          # structural nondeg, float zero
                    continue
                mass = p00 + p01 + p10 + p11
                v = (math.log2(p11) + math.log2(p00)
                     - math.log2(p10) - math.log2(p01))
                rows.append((a, b, mass, v))
                dm += mass
                fl = min(p00, p01, p10, p11) / mass
                floor = fl if floor is None else min(floor, fl)
        out.append({"i": i, "rows": rows, "def_mass": dm,
                    "n_active": len(active), "cell_floor": floor,
                    "underflow": underflow})
    return out


def agg_of_rows(per_i, lam: float, n: int):
    """(aggregate, min per-i excess over i<n, argmin i, def-mass share).
    aggregate = sum_{i<n} w_i (M_i - lam) / sum_{i<n} w_i."""
    num = den = 0.0
    worst = None
    worst_i = None
    for c in per_i:
        if c["i"] >= n or not c["rows"]:
            continue
        w = c["def_mass"]
        mi = sum(m * v for (_a, _b, m, v) in c["rows"]) / w
        num += w * (mi - lam)
        den += w
        if worst is None or mi - lam < worst:
            worst, worst_i = mi - lam, c["i"]
    if den == 0.0:
        return None, None, None
    return num / den, worst, worst_i


def sparse_marginals(u: dict[int, float], lam: float, n: int):
    t = 2.0 ** lam
    supp = sorted(u)
    z = 0.0
    mu = {}
    for a in supp:
        row = sum(u[a] * u[b] * t ** popc(a & b) for b in supp)
        mu[a] = row
        z += row
    mu = {a: v / z for a, v in mu.items()}
    return mu, [sum(v for a, v in mu.items() if a >> j & 1) for j in range(n)]


def sparse_eval(u: dict[int, float], lam: float, n: int):
    """(agg, worst per-i excess, worst i, max marginal, underflow flag)."""
    per_i = sparse_rows(u, lam, n)
    agg, worst, wi = agg_of_rows(per_i, lam, n)
    _mu, marg = sparse_marginals(u, lam, n)
    return agg, worst, wi, max(marg), any(c["underflow"] for c in per_i)


def sinkhorn_atoms(mu: dict[int, float], lam: float, iters: int = 30000,
                   tol: float = 1e-13):
    """Symmetric geometric-damped scaling on atom space (013's iteration
    style); for mu-specified sparse instances (crash) and small-n direct
    cross-checks of the orbit engine."""
    supp = sorted(mu)
    N = len(supp)
    K = [[2.0 ** (lam * popc(supp[a] & supp[b])) for b in range(N)]
         for a in range(N)]
    tgt = [mu[s] for s in supp]
    u = [math.sqrt(max(v, 1e-300)) for v in tgt]
    resid = float("inf")
    for _ in range(iters):
        cur = [u[a] * sum(u[b] * K[a][b] for b in range(N)) for a in range(N)]
        z = sum(cur)
        resid = max(abs(c / z - g) for c, g in zip(cur, tgt))
        if resid < tol:
            break
        u = [ua * math.sqrt(g * z / c) for ua, c, g in zip(u, cur, tgt)]
    return {supp[a]: u[a] for a in range(N)}, resid


# ======================================================================
# sparse engine (exact): same census with Fractions; certified aggregate
# ======================================================================

def sparse_rows_exact(u: dict[int, Fraction], t: Fraction, n: int):
    supp = sorted(u)
    kern = {}
    for a in supp:
        ua = u[a]
        for b in supp:
            kern[(a, b)] = ua * u[b] * t ** popc(a & b)
    out = []
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        prefixes = {A & pmask for A in supp}
        active = {c for c in prefixes
                  if any(A & pmask == c and A & bit for A in supp)
                  and any(A & pmask == c and not A & bit for A in supp)}
        tables: dict[tuple[int, int], list[Fraction]] = {}
        for (A, B), w in kern.items():
            key = (A & pmask, B & pmask)
            tab = tables.setdefault(key, [F0] * 4)
            tab[(2 if A & bit else 0) + (1 if B & bit else 0)] += w
        rows = []
        for (a, b), (p00, p01, p10, p11) in tables.items():
            nondeg = min(p00, p01, p10, p11) > 0
            assert nondeg == (a in active and b in active), \
                f"dichotomy failure at i={i}"      # exact: must never fire
            if nondeg:
                rows.append((p00 + p01 + p10 + p11,
                             (p11 * p00) / (p10 * p01)))
        out.append(rows)
    return out


def sparse_marginals_exact(u: dict[int, Fraction], t: Fraction, n: int):
    supp = sorted(u)
    z = F0
    mu = {}
    for a in supp:
        row = sum(u[a] * u[b] * t ** popc(a & b) for b in supp)
        mu[a] = row
        z += row
    mu = {a: v / z for a, v in mu.items()}
    return [sum(v for a, v in mu.items() if a >> j & 1) for j in range(n)]


def agg_enclosure(rows_by_i, t: Fraction, n: int, terms: int = 12):
    """Certified [lo, hi] for A = sum_{i<n} sum_rows m log2(OR/t) / sum m.
    Also returns per-i M_i - log2 t enclosures.  Running sums are
    compacted to 120-bit dyadics with directed rounding, so the whole
    accumulation stays certified AND cheap."""
    num_lo = num_hi = den_lo = den_hi = F0
    per_i = []
    memo: dict[Fraction, tuple[Fraction, Fraction]] = {}
    for i, rows in enumerate(rows_by_i, start=1):
        if i >= n or not rows:
            per_i.append(None)
            continue
        slo = shi = w_lo = w_hi = F0
        for (m, orr) in rows:
            q = orr / t
            enc = memo.get(q)
            if enc is None:
                enc = log2_enclosure(q, terms=terms)
                memo[q] = enc
            mlo = _dround(m, False)
            mhi = _dround(m, True)
            # m > 0: directed products depend on the sign of the bound
            slo = _dround(slo + (mlo if enc[0] >= 0 else mhi) * enc[0],
                          False)
            shi = _dround(shi + (mhi if enc[1] >= 0 else mlo) * enc[1],
                          True)
            w_lo = _dround(w_lo + mlo, False)
            w_hi = _dround(w_hi + mhi, True)
        per_i.append((slo / w_hi if slo >= 0 else slo / w_lo,
                      shi / w_lo if shi >= 0 else shi / w_hi,
                      (w_lo + w_hi) / 2))
        num_lo = _dround(num_lo + slo, False)
        num_hi = _dround(num_hi + shi, True)
        den_lo = _dround(den_lo + w_lo, False)
        den_hi = _dround(den_hi + w_hi, True)
    a_lo = num_lo / den_hi if num_lo >= 0 else num_lo / den_lo
    a_hi = num_hi / den_lo if num_hi >= 0 else num_hi / den_hi
    return a_lo, a_hi, per_i


# ======================================================================
# orbit engine: measures invariant under permutations of coords {3..n}.
# States s = (a1, a2, k): membership of coords 1, 2 and tail count
# k = |A & {3..n}|.  Same symmetry class as 012's engine; fresh code.
# ======================================================================

def orbit_states(n: int):
    r = n - 2
    return [(a1, a2, k) for a1 in (0, 1) for a2 in (0, 1)
            for k in range(r + 1)]


def wtail_gen(r: int, t):
    """Wt[k][l] = sum over tail sets T_B of size l of t^{|T_A & T_B|},
    T_A fixed of size k.  Arithmetic-generic (float or Fraction t)."""
    tp = [t ** j for j in range(r + 1)]
    Wt = [[None] * (r + 1) for _ in range(r + 1)]
    for k in range(r + 1):
        for l in range(r + 1):
            s = 0
            for j in range(max(0, k + l - r), min(k, l) + 1):
                s += math.comb(k, j) * math.comb(r - k, l - j) * tp[j]
            Wt[k][l] = s
    return Wt


def pair_kernel(f: int, t):
    """PK[s][l] = C(f,s) * Wt_f[s][l] = sum over PAIRS of future sets
    (|X_A|=s, |X_B|=l) of t^{|X_A & X_B|}."""
    Wt = wtail_gen(f, t)
    return [[math.comb(f, s) * Wt[s][l] for l in range(f + 1)]
            for s in range(f + 1)]


def orbit_kernel_matrix(n: int, t, states):
    """M[s][t2] = t^{s1 t1 + s2 t2} * Wt[k_s][k_t]: for a FIXED atom in
    state s, the total kernel weight against the whole orbit t2."""
    r = n - 2
    Wt = wtail_gen(r, t)
    return [[(t ** (s[0] * t2[0] + s[1] * t2[1])) * Wt[s[2]][t2[2]]
             for t2 in states] for s in states]


def orbit_sinkhorn(n: int, tgt: dict, lam: float, iters: int = 12000,
                   tol: float = 1e-11):
    """Symmetric scaling g <- g * sqrt(tgt/cur) on orbit space; tgt is the
    PER-ATOM target mass by state (zero-mass states excluded from the
    support).  Returns (g, relative residual)."""
    states = [s for s in orbit_states(n) if tgt.get(s, 0.0) > 0.0]
    cnt = [math.comb(n - 2, s[2]) for s in states]
    M = orbit_kernel_matrix(n, 2.0 ** lam, states)
    NS = len(states)
    tv = [tgt[s] for s in states]
    g = [math.sqrt(tv[a]) for a in range(NS)]
    resid = float("inf")
    for it in range(iters):
        cur = [g[a] * sum(g[b] * M[a][b] for b in range(NS))
               for a in range(NS)]
        z = sum(cnt[a] * cur[a] for a in range(NS))
        resid = max(abs(cur[a] / (z * tv[a]) - 1.0) for a in range(NS))
        if resid < tol:
            break
        g = [g[a] * math.sqrt(tv[a] * z / cur[a]) for a in range(NS)]
    return {s: g[a] for a, s in enumerate(states)}, resid


def orbit_rows(n: int, t, g: dict):
    """Per i = 1..n: list of (mass, OR) over nondegenerate history orbit
    classes, masses unnormalized.  Arithmetic-generic: t and the values of
    g may be floats or Fractions (math.comb counts are ints).

    i = 1: response coord 1; prefix empty.
    i = 2: response coord 2; prefix = coord-1 value.
    i >= 3: response = a tail coord; prefix = (a1, a2) + q = i-3 tail
    coords; future f = n-i tail coords.  History classes are
    (oa, ob, ja, jb) with the prefix-overlap j summed into the mass factor
    JW[ja][jb] = sum_j Ncl(ja,jb,j) t^j; the OR does not depend on j."""
    r = n - 2
    zero = 0 * t

    def gv(o, c):
        return g.get((o[0], o[1], c), zero)

    orb4 = [(0, 0), (0, 1), (1, 0), (1, 1)]
    out = []

    # ---- i = 1
    PK = pair_kernel(r, t)

    def quad(o_a, ca, o_b, cb, PKm, f):
        acc = zero
        for s in range(f + 1):
            ga = gv(o_a, ca + s)
            if ga == 0:
                continue
            row = PKm[s]
            sub = zero
            for l in range(f + 1):
                gb = gv(o_b, cb + l)
                if gb != 0:
                    sub += row[l] * gb
            acc += ga * sub
        return acc

    # i=1: cells over (alpha, beta) = coord-1 values; sum over coord-2
    reach1 = {al: any(g.get((al, a2, k), zero) != 0
                      for a2 in (0, 1) for k in range(r + 1))
              for al in (0, 1)}
    rows1 = []
    if reach1[0] and reach1[1]:
        cell = {}
        for al in (0, 1):
            for be in (0, 1):
                acc = zero
                for a2 in (0, 1):
                    for b2 in (0, 1):
                        acc += (t ** (a2 * b2)) * quad((al, a2), 0,
                                                       (be, b2), 0, PK, r)
                cell[(al, be)] = (t ** (al * be)) * acc
        mass = cell[(0, 0)] + cell[(0, 1)] + cell[(1, 0)] + cell[(1, 1)]
        orr = (cell[(1, 1)] * cell[(0, 0)]) / (cell[(1, 0)] * cell[(0, 1)])
        rows1.append((mass, orr))
    out.append(rows1)

    # i=2: classes (a, b) = coord-1 prefix values
    reach2 = {(a, al): any(g.get((a, al, k), zero) != 0
                           for k in range(r + 1))
              for a in (0, 1) for al in (0, 1)}
    act2 = {a for a in (0, 1) if reach2[(a, 0)] and reach2[(a, 1)]}
    rows2 = []
    for a in act2:
        for b in act2:
            cell = {}
            for al in (0, 1):
                for be in (0, 1):
                    cell[(al, be)] = (t ** (a * b + al * be)) * \
                        quad((a, al), 0, (b, be), 0, PK, r)
            mass = (cell[(0, 0)] + cell[(0, 1)] + cell[(1, 0)]
                    + cell[(1, 1)])
            orr = (cell[(1, 1)] * cell[(0, 0)]) / (cell[(1, 0)]
                                                   * cell[(0, 1)])
            rows2.append((mass, orr))
    out.append(rows2)

    # i >= 3
    for i in range(3, n + 1):
        q = i - 3
        f = n - i
        PKf = pair_kernel(f, t)
        # S matrices per (oa, ob): S[x][y], x,y in 0..q+1
        # via Gmat_o[x][s] = g(o, x+s)
        Gm = {}
        for o in orb4:
            Gm[o] = [[gv(o, x + s) if x + s <= r else zero
                      for s in range(f + 1)] for x in range(q + 2)]
        # Half_o[l][y] = sum_s PKf[s][l] ... compute H_o = PKf^T Gm_o^T:
        # H_o[l][y] = sum_s PKf[s][l] * Gm_o[y][s]  (PK symmetric would
        # allow PKf[l][s]; keep general)
        Ho = {}
        for o in orb4:
            Go = Gm[o]
            H = [[zero] * (q + 2) for _ in range(f + 1)]
            for l in range(f + 1):
                Hl = H[l]
                for s in range(f + 1):
                    ps = PKf[s][l]
                    if ps == 0:
                        continue
                    for y in range(q + 2):
                        gy = Go[y][s]
                        if gy != 0:
                            Hl[y] += ps * gy
            Ho[o] = H
        # JW[ja][jb] = sum_j Ncl * t^j
        tp = [t ** j for j in range(q + 1)]
        JW = [[zero] * (q + 1) for _ in range(q + 1)]
        for ja in range(q + 1):
            for jb in range(q + 1):
                acc = zero
                for j in range(max(0, ja + jb - q), min(ja, jb) + 1):
                    acc += (math.comb(q, j) * math.comb(q - j, ja - j)
                            * math.comb(q - ja, jb - j)) * tp[j]
                JW[ja][jb] = acc
        # reachability of continuations c = 0..q+1 per orbit pattern
        reach = {o: [any(gv(o, c + s) != 0 for s in range(f + 1))
                     for c in range(q + 3)] for o in orb4}
        rows = []
        for oa in orb4:
            Ra = reach[oa]
            Ga = Gm[oa]
            for ob in orb4:
                Rb = reach[ob]
                Hb = Ho[ob]
                base = t ** (oa[0] * ob[0] + oa[1] * ob[1])
                # S[x][y] = sum_s Ga[x][s] * sum_l PKf[s][l] Gb[y][l]
                #         = sum_l ... use Hb: S[x][y] = sum_s Ga[x][s]*..
                # Hb[l][y] = sum_s PKf[s][l] Gb[y][s] -- careful: Hb was
                # built from Gm[ob], i.e. Hb[l][y] = sum_s PK[s][l] Gb[y][s]
                # We need S[x][y] = sum_{s,l} Ga[x][s] PK[s][l] Gb[y][l]
                #                = sum_l (sum_s Ga[x][s] PK[s][l]) Gb[y][l]
                # PK symmetric (pair count), so sum_s Ga[x][s] PK[s][l]
                # = (H built from Ga)[l][x].  Use Ho[oa].
                Ha = Ho[oa]
                S = [[zero] * (q + 2) for _ in range(q + 2)]
                Gb = Gm[ob]
                for x in range(q + 2):
                    Sx = S[x]
                    for l in range(f + 1):
                        hax = Ha[l][x]
                        if hax == 0:
                            continue
                        for y in range(q + 2):
                            gby = Gb[y][l]
                            if gby != 0:
                                Sx[y] += hax * gby
                for ja in range(q + 1):
                    if not (Ra[ja] and Ra[ja + 1]):
                        continue
                    for jb in range(q + 1):
                        if not (Rb[jb] and Rb[jb + 1]):
                            continue
                        s00 = S[ja][jb]
                        s01 = S[ja][jb + 1]
                        s10 = S[ja + 1][jb]
                        s11 = S[ja + 1][jb + 1]
                        if 0 in (s00, s01, s10, s11):
                            continue      # unreachable class (zero count)
                        jw = JW[ja][jb]
                        if jw == 0:
                            continue
                        orr = t * (s00 * s11) / (s01 * s10)
                        mass = base * jw * (s00 + s01 + s10 + t * s11)
                        rows.append((mass, orr))
        out.append(rows)
    return out


def orbit_agg(rows_by_i, lam: float, n: int):
    """Float aggregate + per-i excess from orbit rows."""
    num = den = 0.0
    worst = worst_i = None
    per_i = []
    for i, rows in enumerate(rows_by_i, start=1):
        if not rows:
            per_i.append(None)
            continue
        W = sum(m for (m, _o) in rows)
        mi = sum(m * math.log2(o) for (m, o) in rows) / W
        per_i.append({"i": i, "M_i": mi, "def_mass": W})
        if i < n:
            num += W * (mi - lam)
            den += W
            if worst is None or mi - lam < worst:
                worst, worst_i = mi - lam, i
    if den == 0.0:
        return None, None, None, per_i
    return num / den, worst, worst_i, per_i


def orbit_marginals(n: int, t, g: dict):
    """Elementwise marginals [coord1, coord2, tail coord] for the tilt
    coupling of orbit potential g.  Arithmetic-generic."""
    r = n - 2
    zero = 0 * t
    states = [s for s in orbit_states(n) if g.get(s, zero) != zero]
    M = orbit_kernel_matrix(n, t, states)
    z = zero
    m1 = m2 = mt = zero
    for a, s in enumerate(states):
        row = zero
        for b, s2 in enumerate(states):
            row += g[s2] * M[a][b]
        w = math.comb(r, s[2]) * g[s] * row
        z += w
        if s[0]:
            m1 += w
        if s[1]:
            m2 += w
        mt += w * Fraction(s[2], r) if isinstance(t, Fraction) \
            else w * (s[2] / r)
    return m1 / z, m2 / z, mt / z


# ----------------------------------------------------------------------
# orbit family targets: PER-ATOM mass by state (a1, a2, k)
# ----------------------------------------------------------------------

def st_bern(n: int, p: float) -> dict:
    return {(a1, a2, k): p ** (a1 + a2 + k) * (1 - p) ** (n - a1 - a2 - k)
            for (a1, a2, k) in orbit_states(n)}


def st_mixture(n: int, comps) -> dict:
    out = {s: 0.0 for s in orbit_states(n)}
    for w, p in comps:
        for s, v in st_bern(n, p).items():
            out[s] += w * v
    return {s: v for s, v in out.items() if v > 0}


def st_slice(n: int, w0: int) -> dict:
    tot = math.comb(n, w0)
    return {(a1, a2, k): 1.0 / tot
            for (a1, a2, k) in orbit_states(n) if a1 + a2 + k == w0}


def st_smslice(n: int, w0: int, tau: float) -> dict:
    out = {s: (1 - tau) * v for s, v in st_slice(n, w0).items()}
    for s, v in st_bern(n, w0 / n).items():
        out[s] = out.get(s, 0.0) + tau * v
    return out


def st_d0bern(n: int, pl: float, ph: float) -> dict:
    out = {s: (1 - pl) * v for s, v in st_bern(n, ph).items()}
    out[(0, 0, 0)] = out.get((0, 0, 0), 0.0) + pl
    return out


def st_crashmix(n: int, p: float, d: float) -> dict:
    """(1-d) Bern(p)^n + d * crash (005's 4-atom family), 012's cells."""
    r = n - 2
    out = {s: (1 - d) * v for s, v in st_bern(n, p).items()}
    out[(0, 1, 0)] = out.get((0, 1, 0), 0.0) + d * 0.33
    out[(0, 0, r)] = out.get((0, 0, r), 0.0) + d * 0.33 / math.comb(r, r)
    out[(1, 1, r)] = out.get((1, 1, r), 0.0) + d * 0.04
    out[(1, 0, 0)] = out.get((1, 0, 0), 0.0) + d * 0.30
    return out


def st_epscrash(n: int, eps: float) -> dict:
    """(1-eps) crash + eps * uniform (007 part K2 genre)."""
    out = {s: eps * 2.0 ** (-n) for s in orbit_states(n)}
    r = n - 2
    out[(0, 1, 0)] += (1 - eps) * 0.33
    out[(0, 0, r)] += (1 - eps) * 0.33
    out[(1, 1, r)] += (1 - eps) * 0.04
    out[(1, 0, 0)] += (1 - eps) * 0.30
    return out


def st_mix_dir(n: int, p: float, nu: str, d: float) -> dict:
    """(1-d) Bern(p)^n + d * direction, for the perturbative sweep around
    the product EQUALITY BOUNDARY (aggregate == 0 identically at d = 0).
    Directions: crash (005), empty (delta_0), full (delta_[n]),
    slice (uniform on the in-regime layer)."""
    out = {s: (1 - d) * v for s, v in st_bern(n, p).items()}
    r = n - 2
    if nu == "crash":
        out[(0, 1, 0)] += d * 0.33
        out[(0, 0, r)] += d * 0.33
        out[(1, 1, r)] += d * 0.04
        out[(1, 0, 0)] += d * 0.30
    elif nu == "empty":
        out[(0, 0, 0)] += d
    elif nu == "full":
        out[(1, 1, r)] += d
    elif nu == "slice":
        w0 = int(RECORD * n)
        for s, v in st_slice(n, w0).items():
            out[s] = out.get(s, 0.0) + d * v
    elif nu == "slice_matched":
        w0 = round(p * n)                 # the layer under the Bern bulk
        for s, v in st_slice(n, w0).items():
            out[s] = out.get(s, 0.0) + d * v
    else:
        raise ValueError(nu)
    return out


def orbit_families(n: int) -> dict:
    """The structured adversary genres, per-atom orbit targets."""
    w0 = int(RECORD * n)          # slice weight, in-regime by construction
    fams = {
        "bern_p3823": st_bern(n, 0.3823),
        "bern_p30": st_bern(n, 0.30),
        "sawin_geo_a": st_mixture(n, [(0.80, 0.28), (0.12, 0.55),
                                      (0.06, 0.75), (0.02, 0.92)]),
        "sawin_geo_b": st_mixture(n, [(0.70, 0.32), (0.20, 0.50),
                                      (0.08, 0.70), (0.02, 0.90)]),
        "sawin_geo_c": st_mixture(n, [(0.60, 0.30), (0.25, 0.45),
                                      (0.10, 0.60), (0.05, 0.80)]),
        "partE_mix": st_mixture(n, [(0.7, 0.30), (0.3, 0.05)]),
        "d0bern_30_51": st_d0bern(n, 0.30, 0.51),
        "d0bern_30_55": st_d0bern(n, 0.30, 0.55),
        f"slice_w{w0}": st_slice(n, w0),
        f"smslice_w{w0}": st_smslice(n, w0, 2.0 ** -8),
        "crash_eps2": st_epscrash(n, 1e-2),
        "crash_eps3": st_epscrash(n, 1e-3),
        "crashmix_383_05": st_crashmix(n, 0.383, 0.05),
        "crashmix_30_05": st_crashmix(n, 0.30, 0.05),
    }
    return fams


def orbit_instance(name: str, n: int, tgt: dict, lam: float):
    """Sinkhorn-fit + census + marginals; returns result dict + g."""
    t0 = time.time()
    g, resid = orbit_sinkhorn(n, tgt, lam)
    rows = orbit_rows(n, 2.0 ** lam, g)
    agg, worst, wi, per_i = orbit_agg(rows, lam, n)
    m1, m2, mt = orbit_marginals(n, 2.0 ** lam, g)
    mm = max(m1, m2, mt)
    # sanity: M_n = lambda identically (Prop 1)
    mn_err = (abs(per_i[n - 1]["M_i"] - lam)
              if per_i[n - 1] is not None else None)
    return {
        "name": name, "n": n, "lambda": lam, "aggregate": agg,
        "min_excess_lt_n": worst, "argmin_i": wi, "max_marginal": mm,
        "in_regime": mm < RECORD, "sinkhorn_resid": resid,
        "Mn_minus_lam_abs": mn_err, "seconds": round(time.time() - t0, 2),
    }, g


# ======================================================================
# sparse structured builders: witness genre at any n, multi-gadget stacks
# ======================================================================

# 007's 10-atom witness, original layout (n=7: marker 1, dilution 2-4,
# response 5, futures 6-7), bit j-1 = coord j.
WITNESS_U7 = {
    0b0100001: 3.802705744653898,      # {1,6}
    0b1000001: 8.094727225702503,      # {1,7}
    0b1010001: 0.004792968160237178,   # {1,5,7}  (light slice)
    0b0000000: 0.2526660836353666,     # {}
    0b0100000: 27.28641199464732,      # {6}
    0b1000000: 12.529780963546393,     # {7}
    0b0110000: 0.3609234435898425,     # {5,6}
    0b0000010: 20.0, 0b0000100: 20.0, 0b0001000: 20.0,
}


def witness_at_n(n: int, wdil: float = 20.0) -> tuple[dict, int]:
    """The witness genre re-laid-out at size n: coord 1 = block marker,
    coords 2..n-3 = fresh-prefix dilution singletons (weight wdil),
    coord n-2 = response, coords n-1, n = futures.  By 007's dilution
    invariance M_response = -0.122 independent of the number of dilution
    atoms; what changes with n is the AGGREGATE bookkeeping (new
    coordinates bring their own w_i and M_i).  Returns (u, response_i)."""
    i_resp = n - 2
    b_resp = 1 << (i_resp - 1)
    b_f1 = 1 << (n - 2)
    b_f2 = 1 << (n - 1)
    b_mark = 1
    u = {
        b_mark | b_f1: 3.802705744653898,             # {1, f1}
        b_mark | b_f2: 8.094727225702503,             # {1, f2}
        b_mark | b_resp | b_f2: 0.004792968160237178,  # light slice
        0: 0.2526660836353666,                        # {}
        b_f1: 27.28641199464732,
        b_f2: 12.529780963546393,
        b_resp | b_f1: 0.3609234435898425,
    }
    for j in range(2, n - 2):                          # dilution singletons
        u[1 << (j - 1)] = wdil
    return u, i_resp


def multi_gadget(n: int, K: int, m: int, slices, wdil: float,
                 sub_w: float = 1.0):
    """K light-slice gadgets stacked: gadget k occupies a marker coord, a
    sub-marker coord, a response coord and m future coords; remaining
    coords are fresh-prefix dilution singletons.  Every atom of gadget k
    carries marker_k, so at gadget k's response coordinate all other
    gadgets' atoms sit on fresh prefixes (deterministic response) and do
    not pollute the active blocks -- the multi-coordinate generalization
    007 lead 1 asked for.  slices = list of K tuples (F0, F1, G0, G1),
    each a vector over 2^m future patterns:
      a-block atoms: marker + sub-marker + response alpha + future x,
      b-block atoms: marker + response alpha + future x.
    Returns (u, [response coords]) or None if it does not fit in n."""
    span = K * (3 + m)
    ndil = n - span
    if ndil < 0:
        return None
    u: dict[int, float] = {}
    resp_coords = []
    pos = 0
    for k in range(K):
        b_mark = 1 << pos
        b_sub = 1 << (pos + 1)
        i_resp = pos + 3
        b_resp = 1 << (pos + 2)
        fut0 = pos + 3
        F0, F1, G0, G1 = slices[k]
        for al, (FA, GA) in enumerate(((F0, G0), (F1, G1))):
            for x in range(1 << m):
                fbits = 0
                for jj in range(m):
                    if x >> jj & 1:
                        fbits |= 1 << (fut0 + jj)
                if FA[x] > 0:
                    key = b_mark | b_sub | (b_resp if al else 0) | fbits
                    u[key] = u.get(key, 0.0) + FA[x] * sub_w
                if GA[x] > 0:
                    key = b_mark | (b_resp if al else 0) | fbits
                    u[key] = u.get(key, 0.0) + GA[x]
        resp_coords.append(i_resp)
        pos += 3 + m
    for j in range(ndil):
        u[1 << (pos + j)] = wdil
    return u, resp_coords


# The witness's slices in (F0, F1, G0, G1) form over 2 future coords
# (x bit0 = first future, bit1 = second future):
WITNESS_SLICES = (
    [0.0, 3.802705744653898, 8.094727225702503, 0.0],    # F0: {f1}, {f2}
    [0.0, 0.0, 0.004792968160237178, 0.0],               # F1: {f2} light
    [0.2526660836353666, 27.28641199464732, 12.529780963546393, 0.0],
    [0.0, 0.3609234435898425, 0.0, 0.0],                 # G1: {f1}
)


def rand_slice(rng, dim: int):
    if rng.random() < 0.5:
        v = [0.0] * dim
        for _ in range(rng.randrange(1, 3)):
            v[rng.randrange(dim)] = math.exp(rng.uniform(-2, 2))
        if sum(v) == 0:
            v[rng.randrange(dim)] = 1.0
        return v
    return [math.exp(rng.uniform(-3, 3)) if rng.random() < 0.7 else 0.0
            for _ in range(dim)]


def clamp_pot(u: dict[int, float], rel: float = 1e-9) -> dict[int, float]:
    """Dynamic-range clamp (007's underflow guard): drop atoms below
    rel * max weight.  Honest support changes, handled by N_i."""
    mx = max(u.values())
    return {k: v for k, v in u.items() if v >= mx * rel}


# ======================================================================
# part X -- cross-checks
# ======================================================================

def part_x() -> dict:
    out: dict = {}
    print("== part X: engine cross-checks ==")

    # X1a: sparse engine reproduces 007's witness numbers at lambda = 3.5
    lam, n = 3.5, 7
    per_i = sparse_rows(WITNESS_U7, lam, n)
    agg, worst, wi, per = agg_of_rows(per_i, lam, n), None, None, None
    a, w, wi2 = agg
    _mu, marg = sparse_marginals(WITNESS_U7, lam, n)
    w5 = per_i[4]
    m5 = sum(m * v for (_a, _b, m, v) in w5["rows"]) / w5["def_mass"]
    print(f"[X1] witness n=7 lam=3.5: M_5 - lam = {m5 - lam:+.6f} "
          f"(007: -0.122033), aggregate = {a:+.6f} (007: +1.84), "
          f"max marginal = {max(marg):.4f} (007: 0.3178)")
    ok1 = abs(m5 - lam + 0.122033) < 5e-6 and abs(a - 1.8447) < 5e-4
    out["x1"] = {"M5_minus_lam": m5 - lam, "aggregate": a,
                 "max_marginal": max(marg), "pass": ok1}

    # X1b: exact aggregate at t = 181/16 vs 013 part C (+1.844669)
    tF = Fraction(181, 16)
    uF = {k: Fraction(v) for k, v in WITNESS_U7.items()}
    rows = sparse_rows_exact(uF, tF, n)
    lo, hi, _pi = agg_enclosure(rows, tF, n)
    print(f"[X1] exact aggregate at t=181/16 in "
          f"[{float(lo):+.9f}, {float(hi):+.9f}] (013: +1.844669), "
          f"width {float(hi - lo):.1e}")
    ok1b = lo > 0 and abs(float(lo) - 1.844669) < 1e-5
    out["x1_exact"] = {"lo": float(lo), "hi": float(hi), "pass": ok1b}

    # X2: orbit engine vs direct sparse+Sinkhorn on full supports, n=6,7.
    # Direct path: atom-space Sinkhorn on the full 2^n support, then the
    # sparse census; orbit path: orbit Sinkhorn + orbit rows.
    max_dv = 0.0
    for (name, n2, tgt_orbit, mu_full) in x2_instances():
        lam2 = 1.0
        gg, res_o = orbit_sinkhorn(n2, tgt_orbit, lam2)
        rows_o = orbit_rows(n2, 2.0 ** lam2, gg)
        agg_o, w_o, _wi, per_o = orbit_agg(rows_o, lam2, n2)
        u_d, res_d = sinkhorn_atoms(mu_full, lam2)
        per_d = sparse_rows(u_d, lam2, n2)
        agg_d, w_d, _x = agg_of_rows(per_d, lam2, n2)
        dv = abs(agg_o - agg_d)
        # per-i comparison
        for i in range(1, n2 + 1):
            po = per_o[i - 1]
            rows_i = per_d[i - 1]["rows"]
            if po is None or not rows_i:
                assert po is None and not rows_i, (name, i)
                continue
            mid = sum(m * v for (_a, _b, m, v) in rows_i) \
                / per_d[i - 1]["def_mass"]
            dv = max(dv, abs(po["M_i"] - mid))
        max_dv = max(max_dv, dv)
        print(f"[X2] orbit vs direct, {name} n={n2}: max |diff| = {dv:.2e} "
              f"(resids {res_o:.1e}/{res_d:.1e})")
    out["x2_max_diff"] = max_dv
    out["x2_pass"] = max_dv < 5e-7

    # X3: exact orbit rows vs float orbit rows on one instance
    n3, lam3 = 8, 1.0
    tgt = st_mixture(n3, [(0.7, 0.30), (0.3, 0.05)])
    g3, _res = orbit_sinkhorn(n3, tgt, lam3)
    rows_f = orbit_rows(n3, 2.0 ** lam3, g3)
    agg_f, _w, _i, _p = orbit_agg(rows_f, lam3, n3)
    g3x = {s: Fraction(v) for s, v in g3.items()}
    t3x = Fraction(2.0 ** lam3)
    rows_x = orbit_rows(n3, t3x, g3x)
    lo3, hi3, _pp = agg_enclosure(rows_x, t3x, n3)
    d3 = max(abs(agg_f - float(lo3)), abs(agg_f - float(hi3)))
    print(f"[X3] exact vs float orbit rows (partE_mix n=8): float "
          f"{agg_f:+.9f}, exact [{float(lo3):+.9f}, {float(hi3):+.9f}], "
          f"diff {d3:.1e}")
    out["x3_diff"] = d3
    out["x3_pass"] = d3 < 1e-7

    save("or_agg_probe_partX.json", out)
    print(f"[X] ALL {'PASS' if all(v for k, v in out.items() if str(k).endswith('pass')) else 'FAIL'}")
    return out


def x2_instances():
    """(name, n, orbit target, full-support mu) for the X2 cross-check."""
    def bern_full(n, p):
        return {m: p ** popc(m) * (1 - p) ** (n - popc(m))
                for m in range(1 << n)}

    def mix_full(n, comps):
        mu = {m: 0.0 for m in range(1 << n)}
        for w, p in comps:
            for m2, v in bern_full(n, p).items():
                mu[m2] += w * v
        return {m: v for m, v in mu.items() if v > 0}

    insts = []
    n = 6
    w0 = 2
    mu = {m: 1.0 / math.comb(n, w0) for m in range(1 << n)
          if popc(m) == w0}
    insts.append((f"slice_w{w0}", n, st_slice(n, w0), mu))
    insts.append(("sawin_geo_a", n,
                  st_mixture(n, [(0.80, 0.28), (0.12, 0.55),
                                 (0.06, 0.75), (0.02, 0.92)]),
                  mix_full(n, [(0.80, 0.28), (0.12, 0.55),
                               (0.06, 0.75), (0.02, 0.92)])))
    # eps-crash n=6
    eps = 1e-2
    mu_ec = {m: eps / (1 << n) for m in range(1 << n)}
    full = (1 << n) - 1
    for a, wgt in ((0b10, 0.33), (full ^ 0b11, 0.33), (full, 0.04),
                   (0b01, 0.30)):
        mu_ec[a] += (1 - eps) * wgt
    insts.append(("crash_eps2", n, st_epscrash(n, eps), mu_ec))
    n = 7
    insts.append(("d0bern_30_51", n, st_d0bern(n, 0.30, 0.51),
                  mix_full(n, [(0.30, 0.0), (0.70, 0.51)])))
    return insts


# ======================================================================
# part O -- orbit family census over the n grid
# ======================================================================

def part_o(ns=NS_GRID) -> dict:
    print("== part O: orbit family census ==")
    all_res = {}
    tightest = {}
    for n in ns:
        res = []
        fams = orbit_families(n)
        lams = lam_grid(n)
        t_n = time.time()
        for name, tgt in fams.items():
            for lam in lams:
                rec, _g = orbit_instance(name, n, tgt, lam)
                res.append(rec)
        # per-n summary over IN-REGIME instances
        inreg = [r for r in res if r["in_regime"] and r["aggregate"] is not None]
        n_viol = sum(1 for r in inreg if r["aggregate"] < -1e-9)
        worst = min(inreg, key=lambda r: r["aggregate"]) if inreg else None
        max_resid = max(r["sinkhorn_resid"] for r in res)
        max_mn = max(r["Mn_minus_lam_abs"] for r in res
                     if r["Mn_minus_lam_abs"] is not None)
        print(f"[O] n={n}: {len(res)} instances ({len(inreg)} in-regime), "
              f"aggregate violations: {n_viol}; min in-regime aggregate = "
              f"{worst['aggregate']:+.6f} at {worst['name']} "
              f"lam={worst['lambda']:g} (i* = {worst['argmin_i']}, "
              f"min per-i excess {worst['min_excess_lt_n']:+.4f}); "
              f"max resid {max_resid:.1e}, max |M_n - lam| {max_mn:.1e}; "
              f"{time.time() - t_n:.0f}s")
        all_res[n] = res
        tightest[n] = worst
        save(f"or_agg_probe_partO_n{n}.json",
             {"n": n, "results": res, "n_in_regime_violations": n_viol,
              "min_in_regime_aggregate": worst})
    save("or_agg_probe_partO_summary.json",
         {str(n): tightest[n] for n in all_res})
    return {"results": all_res, "tightest": tightest}


# ======================================================================
# part P -- perturbative sweep around the product equality boundary.
# Product measures have proportional slices, hence OR == 2^lambda
# pointwise and aggregate == 0 IDENTICALLY: they are the equality set of
# the aggregated control.  A violation therefore only needs the
# d-derivative (or curvature) of the aggregate to be negative along SOME
# direction at SOME (n, p, lambda) -- and 012 showed the perturbative
# downward-budget grows with n on exactly this genre.  Sweep d and fit
# the small-d scaling.
# ======================================================================

def part_p(ns=(10, 20, 32)) -> dict:
    print("== part P: perturbative d-sweep around product boundary ==")
    res = []
    dirs = ("crash", "empty", "full", "slice")
    dgrid = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
    for n in ns:
        lams = sorted({round(math.log2(1.06), 6), 0.5,
                       round(4.847 / (n - 3), 6)})
        for p in (0.3823, 0.30):
            for nu in dirs:
                for lam in lams:
                    aggs = []
                    for d in dgrid:
                        tgt = st_mix_dir(n, p, nu, d)
                        rec, _g = orbit_instance(
                            f"mixdir_{nu}_p{p:g}", n, tgt, lam)
                        rec["d"] = d
                        res.append(rec)
                        aggs.append((d, rec["aggregate"],
                                     rec["in_regime"]))
                    # small-d scaling: ratio agg(2d)/agg(d) ~ 2^order
                    a1, a2 = aggs[0][1], aggs[1][1]
                    order = (math.log2(a2 / a1)
                             if a1 and a2 and a1 > 0 and a2 > 0 else None)
                    neg = [a for a in aggs if a[1] is not None
                           and a[1] < -1e-12 and a[2]]
                    print(f"[P] n={n} p={p:g} dir={nu:6s} lam={lam:g}: "
                          f"agg(d) = "
                          + " ".join(f"{d:g}:{a:+.2e}" for d, a, _r in aggs)
                          + (f"  order~d^{order:.2f}" if order else "")
                          + ("  <-- NEGATIVE" if neg else ""))
    inreg = [r for r in res if r["in_regime"] and r["aggregate"] is not None]
    nv = sum(1 for r in inreg if r["aggregate"] < -1e-12)
    worst = min(inreg, key=lambda r: r["aggregate"])
    print(f"[P] {len(res)} instances ({len(inreg)} in-regime); negative "
          f"aggregates: {nv}; min in-regime aggregate = "
          f"{worst['aggregate']:+.3e} at n={worst['n']} {worst['name']} "
          f"d={worst['d']:g} lam={worst['lambda']:g}")
    save("or_agg_probe_partP.json",
         {"results": res, "n_negative": nv, "min_in_regime": worst})
    return {"results": res, "min_in_regime": worst}


# ======================================================================
# part T -- trend extension past n = 32 on the softest observed
# directions (part P's slice-direction curvature at p = 0.30 decays with
# n: 1.35e-5 -> 6.2e-7 over n = 10 -> 32).  The orbit engine is cheap
# enough to chase the trend to n = 64: does the margin cross zero, decay
# to 0 polynomially, or plateau?  Also tracks the tightest census family
# (crashmix_383_05) and the crash direction at the window lambda.
# ======================================================================

def part_t(ns=(32, 40, 48, 56, 64)) -> dict:
    print("== part T: margin trend past n = 32 ==")
    res = []
    for n in ns:
        for (p, nu, lam, d) in (
                (0.30, "slice", 0.084, 0.05),
                (0.30, "slice", 0.5, 0.05),
                (0.30, "slice_matched", 0.084, 0.05),
                (0.3823, "slice", 0.084, 0.05),
                (0.3823, "slice_matched", 0.084, 0.05),
                (0.3823, "crash", 0.084, 0.05),
                (0.3823, "crash", round(4.847 / (n - 3), 6), 0.05),
        ):
            tgt = st_mix_dir(n, p, nu, d)
            rec, _g = orbit_instance(f"mixdir_{nu}_p{p:g}_d{d:g}", n, tgt,
                                     lam)
            rec["d"] = d
            res.append(rec)
        # the tightest census family, at its tightest lambda
        tgt = st_crashmix(n, 0.383, 0.05)
        lamw = round(math.log2(1.06), 6)
        rec, _g = orbit_instance("crashmix_383_05", n, tgt, lamw)
        res.append(rec)
        row = [r for r in res if r["n"] == n]
        print(f"[T] n={n}: "
              + "  ".join(f"{r['name']}@{r['lambda']:g}:"
                          f"{r['aggregate']:+.3e}" for r in row))
    neg = [r for r in res if r["in_regime"] and r["aggregate"] is not None
           and r["aggregate"] < -1e-12]
    print(f"[T] negative in-regime aggregates: {len(neg)}")
    save("or_agg_probe_partT.json", {"results": res, "n_negative": len(neg)})
    return {"results": res}


# ======================================================================
# part S -- sparse structured census: crash/mirror, witness genre,
# multi-gadget stacks
# ======================================================================

def part_s(ns=NS_GRID) -> dict:
    print("== part S: sparse structured census ==")
    res = []

    def run_pot(name, u, n, lam):
        agg, worst, wi, mm, uf = sparse_eval(u, lam, n)
        rec = {"name": name, "n": n, "lambda": lam, "aggregate": agg,
               "min_excess_lt_n": worst, "argmin_i": wi,
               "max_marginal": mm, "in_regime": mm < RECORD,
               "underflow": uf}
        res.append(rec)
        return rec

    full_atoms = lambda n: (1 << n) - 1

    for n in ns:
        lams = lam_grid(n)
        # crash / mirror (4 atoms; Theorem B territory -- margin trend)
        full = full_atoms(n)
        crash = {0b10: 0.33, full ^ 0b11: 0.33, full: 0.04, 0b01: 0.30}
        mirror = {0b10 | (full ^ 0b11): 0.3, 0b00: 0.2, full: 0.3,
                  0b01: 0.2}
        for lam in lams:
            for fname, mu in (("crash", crash), ("mirror", mirror)):
                u, resid = sinkhorn_atoms(mu, lam)
                rec = run_pot(f"{fname}", u, n, lam)
                rec["sinkhorn_resid"] = resid

        # witness genre re-laid-out at n, dilution-weight sweep: the
        # tightest in-regime aggregate the genre can reach at this n
        for lam in (0.5, 1.0, 2.0, 3.5, 5.0):
            best = None
            for wdil in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 20.0, 40.0):
                u, i_resp = witness_at_n(n, wdil)
                agg, worst, wi, mm, uf = sparse_eval(u, lam, n)
                rec = {"name": f"witness_wdil{wdil:g}", "n": n,
                       "lambda": lam, "aggregate": agg,
                       "min_excess_lt_n": worst, "argmin_i": wi,
                       "max_marginal": mm, "in_regime": mm < RECORD,
                       "underflow": uf}
                res.append(rec)
                if rec["in_regime"] and (best is None
                                         or agg < best["aggregate"]):
                    best = rec
            if best is not None:
                print(f"[S] n={n} lam={lam:g}: witness-genre tightest "
                      f"in-regime aggregate {best['aggregate']:+.6f} "
                      f"({best['name']}, min per-i excess "
                      f"{best['min_excess_lt_n']:+.4f} at i="
                      f"{best['argmin_i']}, marg {best['max_marginal']:.4f})")

        # multi-gadget stacks seeded with the witness slices: K gadgets,
        # m = 2 futures each (K as large as fits), sub-marker weight and
        # dilution weight swept coarsely
        for lam in (2.0, 3.5):
            K = (n - 2) // 5
            if K < 1:
                continue
            for wdil in (1.0, 4.0, 16.0):
                built = multi_gadget(n, K, 2, [WITNESS_SLICES] * K, wdil)
                if built is None:
                    continue
                u, resp = built
                u = clamp_pot(u)
                rec = run_pot(f"stack_K{K}_wdil{wdil:g}", u, n, lam)
        stack = [r for r in res if r["n"] == n and r["name"].startswith("stack")]
        if stack:
            b = min((r for r in stack if r["in_regime"]),
                    key=lambda r: r["aggregate"], default=None)
            if b:
                print(f"[S] n={n}: tightest in-regime witness-slice stack "
                      f"aggregate {b['aggregate']:+.6f} ({b['name']} "
                      f"lam={b['lambda']:g})")

    inreg = [r for r in res if r["in_regime"] and r["aggregate"] is not None]
    n_viol = sum(1 for r in inreg if r["aggregate"] < -1e-9)
    worst = min(inreg, key=lambda r: r["aggregate"])
    print(f"[S] total {len(res)} instances ({len(inreg)} in-regime); "
          f"aggregate violations: {n_viol}; min in-regime aggregate = "
          f"{worst['aggregate']:+.6f} at n={worst['n']} {worst['name']} "
          f"lam={worst['lambda']:g}")
    save("or_agg_probe_partS.json",
         {"results": res, "n_in_regime_violations": n_viol,
          "min_in_regime": worst})
    return {"results": res, "min_in_regime": worst}


# ======================================================================
# part H -- adversarial hill-climbs on the aggregate objective
# ======================================================================

MARG_TARGET = 0.375


def climb(u0: dict[int, float], n: int, lam: float, rng, steps: int = 400):
    """Multiplicative hill-climb minimizing aggregate + marginal penalty,
    dynamic-range-clamped.  Returns (best value, best u, max marginal)."""
    def obj(u):
        u = clamp_pot(u)
        agg, _w, _i, mm, uf = sparse_eval(u, lam, n)
        if agg is None or uf:
            return float("inf"), None, None
        return agg + 25.0 * max(0.0, mm - MARG_TARGET), agg, mm

    cur, cagg, cmm = obj(u0)
    u = dict(u0)
    if not math.isfinite(cur):
        return None, None, None
    for step in range(steps):
        sigma = 1.0 * (1 - step / steps) + 0.05
        cand = dict(u)
        r = rng.random()
        if r < 0.12 and len(cand) > 5:            # drop an atom
            cand.pop(rng.choice(sorted(cand)))
        elif r < 0.24:                            # add a random atom
            a = rng.randrange(1 << n)
            cand[a] = math.exp(rng.uniform(-2, 2))
        else:                                     # perturb one weight
            a = rng.choice(sorted(cand))
            cand[a] = max(cand[a], 1e-300) * math.exp(rng.gauss(0, sigma))
        v, agg, mm = obj(cand)
        if v < cur:
            u, cur, cagg, cmm = clamp_pot(cand), v, agg, mm
    return cagg, u, cmm


def part_h() -> dict:
    print("== part H: adversarial hill-climbs on the aggregate ==")
    rng = random.Random(140001)
    endpoints = []
    best_overall = None

    def record(tag, n, lam, cagg, u, cmm):
        nonlocal best_overall
        e = {"tag": tag, "n": n, "lambda": lam, "aggregate": cagg,
             "max_marginal": cmm, "n_atoms": len(u) if u else None}
        endpoints.append(e)
        if (cagg is not None and cmm is not None and cmm < RECORD
                and (best_overall is None
                     or cagg < best_overall["aggregate"])):
            best_overall = dict(e)
            best_overall["u"] = {format(k, "b"): v for k, v in u.items()}

    # H-a: free sparse supports at n = 10, 14 (8..18 atoms)
    for n in (10, 14):
        for restart in range(10):
            k = rng.randrange(8, 19)
            atoms = rng.sample(range(1 << n), k)
            u0 = {a: math.exp(rng.uniform(-2, 2)) for a in atoms}
            lam = rng.choice([0.5, 1.0, 2.0, 3.5])
            cagg, u, cmm = climb(u0, n, lam, rng, steps=400)
            record("free", n, lam, cagg, u, cmm)

    # H-b: witness-seeded at n = 10, 14, 20
    for n in (10, 14, 20):
        for lam in (2.0, 3.5):
            for restart in range(3):
                u0, _i = witness_at_n(n, wdil=rng.choice([1.0, 4.0, 20.0]))
                if restart:                        # jiggle the seed
                    u0 = {k: v * math.exp(rng.gauss(0, 0.1))
                          for k, v in u0.items()}
                cagg, u, cmm = climb(u0, n, lam, rng, steps=350)
                record("witness_seed", n, lam, cagg, u, cmm)

    # H-c: multi-gadget-parametrized at n = 14, 20, 24: witness slices +
    # random slices, K gadgets, climbing atom weights freely
    for n in (14, 20, 24):
        K = (n - 2) // 5
        for restart in range(4):
            m = 2
            slices = []
            for k in range(K):
                if rng.random() < 0.5:
                    slices.append(WITNESS_SLICES)
                else:
                    slices.append(tuple(rand_slice(rng, 1 << m)
                                        for _ in range(4)))
            built = multi_gadget(n, K, m, slices, wdil=4.0)
            if built is None:
                continue
            u0, _resp = built
            lam = rng.choice([2.0, 3.5])
            cagg, u, cmm = climb(clamp_pot(u0), n, lam, rng, steps=350)
            record("stack_seed", n, lam, cagg, u, cmm)

    # H-d: small-lambda climbs (Theorem C makes the first order in lambda
    # a perfect square per i, so a small-lambda kill needs all D-sums ~ 0
    # with negative curvature -- let the search try, near the window law)
    for n in (10, 14):
        for lam in (0.084, 0.15, 0.3, 4.847 / (n - 3)):
            for restart in range(4):
                if restart % 2 == 0:
                    k = rng.randrange(8, 16)
                    u0 = {a: math.exp(rng.uniform(-2, 2))
                          for a in rng.sample(range(1 << n), k)}
                else:
                    u0, _i = witness_at_n(n, wdil=4.0)
                    u0 = {k2: v * math.exp(rng.gauss(0, 0.1))
                          for k2, v in u0.items()}
                cagg, u, cmm = climb(u0, n, lam, rng, steps=350)
                record("small_lam", n, lam, cagg, u, cmm)

    ok = [e for e in endpoints if e["aggregate"] is not None]
    nv = sum(1 for e in ok if e["aggregate"] < -1e-9
             and e["max_marginal"] < RECORD)
    lows = sorted(e["aggregate"] for e in ok
                  if e["max_marginal"] < RECORD)[:5]
    print(f"[H] {len(ok)} climb endpoints; in-regime aggregate violations: "
          f"{nv}; five lowest in-regime endpoints: "
          + " ".join(f"{v:+.5f}" for v in lows))
    if best_overall:
        print(f"[H] lowest in-regime endpoint: {best_overall['aggregate']:+.6f} "
              f"({best_overall['tag']} n={best_overall['n']} "
              f"lam={best_overall['lambda']:g}, "
              f"{best_overall['n_atoms']} atoms)")
    save("or_agg_probe_partH.json",
         {"endpoints": endpoints, "n_in_regime_violations": nv,
          "best": best_overall})
    return {"endpoints": endpoints, "best": best_overall}


# ======================================================================
# part E -- exact rational certification of the tightest points
# ======================================================================

def certify_sparse(name: str, u: dict[int, float], lam_num: int,
                   lam_den: int, n: int, use_t: Fraction | None = None):
    """Certify the aggregate sign for the exact dyadic-rational potential
    Fraction(u) at exact rational tilt.  If use_t is given it IS the tilt
    t (lambda = log2 t enclosed); else t = 2^{lam_num/lam_den} is
    float-rationalized (dyadic).  Returns dict."""
    t = use_t if use_t is not None else Fraction(2.0 ** (lam_num / lam_den))
    uF = {k: Fraction(v) for k, v in u.items()}
    t0 = time.time()
    rows = sparse_rows_exact(uF, t, n)
    lo, hi, per_i = agg_enclosure(rows, t, n)
    marg = sparse_marginals_exact(uF, t, n)
    mm = max(marg)
    rec = {"name": name, "n": n, "t": str(t), "engine": "sparse",
           "agg_lo": float(lo), "agg_hi": float(hi),
           "enclosure_width": float(hi - lo),
           "certified_sign": ("positive" if lo > 0 else
                              "NEGATIVE" if hi < 0 else "indeterminate"),
           "max_marginal_exact_lt_record": bool(mm < Fraction(38271, 100000)),
           "max_marginal": float(mm), "seconds": round(time.time() - t0, 1)}
    print(f"[E] {name} (sparse, n={n}, t={t.numerator}/{t.denominator}"
          f"{'' if t.denominator < 10**6 else ' ~ ' + format(float(t), '.6g')}):"
          f" aggregate in [{float(lo):+.9f}, {float(hi):+.9f}] -> "
          f"{rec['certified_sign']}; exact max marginal "
          f"{float(mm):.6f} {'<' if rec['max_marginal_exact_lt_record'] else '>='} 0.38271; {rec['seconds']}s")
    return rec


def _quant(v: float, bits: int = 36) -> Fraction:
    """Nearest dyadic with a bits-bit mantissa (exact, sign-preserving)."""
    if v == 0.0:
        return F0
    m, e = math.frexp(v)
    return Fraction(round(m * (1 << bits)), 1 << bits) * Fraction(2) ** e


def certify_orbit(name: str, n: int, tgt: dict, lam: float,
                  terms: int = 10, use_t: Fraction | None = None,
                  quant_bits: int = 36):
    """Sinkhorn-fit at float lambda, rationalize the potential (dyadic,
    quantized to quant_bits-bit mantissas to keep exact arithmetic fast)
    and the tilt (use_t if given -- prefer a short rational like 271/256
    -- else the dyadic float-rationalization of 2^lambda), then certify
    the aggregate sign for that EXACT rational measure (013's
    rationalized-refit scope: the certificate is about a mu-tilde within
    ~Sinkhorn-residual + 2^-quant_bits relative of the named family,
    which suffices for the universally-quantified claim)."""
    t0 = time.time()
    if use_t is not None:
        t = use_t
        lam = math.log2(float(t.numerator) / float(t.denominator))
    else:
        t = Fraction(2.0 ** lam)
    g, resid = orbit_sinkhorn(n, tgt, lam)
    gF = {s: _quant(v, quant_bits) for s, v in g.items()}
    rows = orbit_rows(n, t, gF)
    lo, hi, _per = agg_enclosure(rows, t, n, terms=terms)
    m1, m2, mt = orbit_marginals(n, t, gF)
    mm = max(m1, m2, mt)
    rec = {"name": name, "n": n, "lambda": lam, "t": str(t),
           "engine": "orbit", "sinkhorn_resid": resid,
           "agg_lo": float(lo), "agg_hi": float(hi),
           "enclosure_width": float(hi - lo),
           "certified_sign": ("positive" if lo > 0 else
                              "NEGATIVE" if hi < 0 else "indeterminate"),
           "max_marginal_exact_lt_record": bool(mm < Fraction(38271, 100000)),
           "max_marginal": float(mm), "seconds": round(time.time() - t0, 1)}
    print(f"[E] {name} (orbit, n={n}, lam~{lam:g}): aggregate in "
          f"[{float(lo):+.9f}, {float(hi):+.9f}] -> {rec['certified_sign']};"
          f" exact max marginal {float(mm):.6f}; resid {resid:.1e}; "
          f"{rec['seconds']}s")
    return rec


def part_e(o_tightest=None) -> dict:
    """Certify: (a) the witness-genre tightest in-regime point per n
    (sparse, cheap, all n up to 32); (b) the tightest orbit point per n
    (from part O's summary if available on disk); (c) rational-lambda
    witness-genre points (t = 4 exactly, lambda = 2)."""
    print("== part E: exact rational certification ==")
    out = {"sparse": [], "orbit": []}

    # (a) witness genre at every n: the wdil giving the tightest in-regime
    # aggregate at lam = 3.5 (float pre-scan), certified at t = 181/16
    for n in NS_GRID:
        best = None
        for wdil in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 20.0):
            u, _i = witness_at_n(n, wdil)
            agg, _w, _wi, mm, _uf = sparse_eval(u, 3.5, n)
            if mm < RECORD and (best is None or agg < best[0]):
                best = (agg, wdil, u)
        if best is None:
            continue
        _agg, wdil, u = best
        rec = certify_sparse(f"witness_n{n}_wdil{wdil:g}", u, 7, 2, n,
                             use_t=Fraction(181, 16))
        out["sparse"].append(rec)

    # (c) fully rational statement at lambda = 2 exactly (t = 4): the
    # witness genre at n = 20 (wdil = 20 keeps it in-regime at t = 4)
    u20, _i = witness_at_n(20, 20.0)
    rec = certify_sparse("witness_n20_t4_lam2exact", u20, 2, 1, 20,
                         use_t=Fraction(4))
    out["sparse"].append(rec)

    # (b) tightest NON-PRODUCT orbit point per n, read from part O
    # checkpoints.  Product (bern) cells are the exact equality boundary
    # of the aggregated control (proportional slices => OR == 2^lambda
    # pointwise => aggregate == 0 identically); an enclosure there
    # straddles 0 by construction and certifies nothing, so they are
    # excluded from the "tightest point" and recorded as boundary.
    for n in NS_GRID:
        p = DATA / f"or_agg_probe_partO_n{n}.json"
        if not p.exists():
            continue
        allres = json.loads(p.read_text())["results"]
        cand = [r for r in allres if r["in_regime"]
                and r["aggregate"] is not None
                and not r["name"].startswith("bern")]
        if not cand:
            continue
        worst = min(cand, key=lambda r: r["aggregate"])
        name = worst["name"]
        fams = orbit_families(n)
        if name not in fams:
            continue
        # census min sits at the smallest grid lambda (~0.083-0.084);
        # certify at the nearby SHORT exact rational tilt t = 271/256
        # (lambda = log2(271/256) ~ 0.08215), 013's clean-rational move
        rec = certify_orbit(f"{name}_n{n}", n, fams[name], 0.0,
                            terms=10, use_t=Fraction(271, 256))
        out["orbit"].append(rec)

    # (d) the softest direction found (part P/T): slice-direction
    # perturbation of Bern(0.30) at d = 0.05, small lambda -- the
    # smallest nonzero margin in the whole probe (~6e-7 at n = 32).
    # Certify at n = 20 and 32 (enclosure width ~1e-10 << margin).
    for n_soft in (20, 32):
        tgt = st_mix_dir(n_soft, 0.30, "slice", 0.05)
        rec = certify_orbit(f"mixdir_slice_p0.3_d0.05_n{n_soft}", n_soft,
                            tgt, 0.0, terms=10, use_t=Fraction(271, 256))
        out["orbit"].append(rec)

    save("or_agg_probe_partE.json", out)
    return out


def main() -> None:
    parts = [p.upper() for p in sys.argv[1:] if p.upper() in "XOPTSHE"] \
        or list("XOPTSHE")
    t0 = time.time()
    if "X" in parts:
        part_x()
    if "O" in parts:
        part_o()
    if "P" in parts:
        part_p()
    if "T" in parts:
        part_t()
    if "S" in parts:
        part_s()
    if "H" in parts:
        part_h()
    if "E" in parts:
        part_e()
    print(f"done in {time.time() - t0:.0f}s; checkpoints in {DATA}")


if __name__ == "__main__":
    main()
