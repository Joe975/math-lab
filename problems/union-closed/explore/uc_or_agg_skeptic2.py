#!/usr/bin/env python3
"""Attempt 015: adversarial skeptic review of 014 (aggregated OR probe).

Independent re-implementation, from scratch, of everything 014's headline
numbers rest on.  NOTHING on the critical path is shared with
uc_or_agg_probe.py:

  * log2 certification: binary DIGIT-EXTRACTION (square-and-compare) with
    directed integer shifts -- a different algorithm from 014's atanh
    series + dyadic compaction.  Soundness proof in log2_bracket's
    docstring.
  * exact aggregate: integer-scaled row builders (all masses exact
    integers on one common scale; only the log2 values are intervals), so
    the denominator of A is EXACT and only the numerator is an interval.
  * own direct full-support census evaluator (float + the exact sparse
    one), own tail-exchangeable orbit evaluator re-derived by hand, own
    Sinkhorn (exponent-0.4-damped scaling, not 014's sqrt damping).

uc_or_agg_probe.py is imported ONLY (a) to reproduce the exact rational
refit measures mu-tilde its certificates are ABOUT (part C re-certifies
the same stated inputs -- the measure spec is data, not verification
logic), and (b) to compare its float orbit rows against mine (part X).

Parts:
  K  attack 014's enclosure kit: directional unit tests of _dround /
     _ln_1_to_2_fp / log2_enclosure against my bracket enclosure on
     adversarial inputs (near powers of 2, tiny/huge, tight rationals);
     brute-force validation of the orbit combinatorics (Ncl, PK) that my
     hand re-derivation predicts.
  X  formula cross-checks: my direct evaluator vs my orbit evaluator at
     n = 8 AND 9 (014 only did 6, 7), on crashmix / sawin / slice-dir /
     d0bern; my direct vs 014's sparse engine on the witness; explicit
     j-independence check of class ORs at n = 9.
  C  independent re-certification of five of 014's fifteen certificates
     (incl. the softest, +5.99e-7 at n = 32, and the fully rational
     t = 4 point), on the identical exact rational inputs, plus one
     fresh certificate at the softest point from MY OWN Sinkhorn+quant.
  T  trend attack: the n ~ 48 dip re-computed with my engine, off-grid
     lambda (below 014's grid floor 0.079 and between its big-grid
     points), and the census-minimum's lambda-scaling.
  V  new-adversary hunt in regions 014's search could not see:
     witness-DIRECTION perturbation of a dense product bulk (non-tail-
     symmetric direction off the equality boundary, n = 10),
     boundary-of-regime marginals (max marginal pushed to ~0.3827),
     tiny lambda.

Standard library only.  Deterministic.  Runtime ~10-25 min total.
Checkpoints: problems/union-closed/data/oragg_sk2_part[KXCTV].json
Usage: python problems/union-closed/explore/uc_or_agg_skeptic2.py [parts]
"""

from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"
sys.path.insert(0, str(Path(__file__).resolve().parent))

RECORD = Fraction(38271, 100000)
RECORD_F = 0.38271


def popc(x: int) -> int:
    return bin(x).count("1")


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)
    (DATA / name).write_text(
        json.dumps(obj, indent=1, sort_keys=True, default=str) + "\n",
        encoding="utf-8")


# ======================================================================
# my certified log2: binary digit extraction (square-and-compare).
#
# For y in [1, 2) kept as an integer mantissa m with y = m / 2^WB:
# repeat B times: y <- y^2; if y >= 2 emit digit 1 and halve else emit 0.
# Then log2(y_0) = sum_j d_j 2^{-j} + 2^{-B} log2(y_B) with y_B in [1,2).
#
# Soundness with directed integer shifts (proved, used below):
#   DOWN pass: m <- (m*m) >> (WB + d)  (floor).  By induction the tracked
#   value y_j = m_j/2^WB satisfies y_j <= v_j, where v_j is the EXACT
#   process that uses the algorithm's emitted digits d_j
#   (v_{j+1} = v_j^2/2^{d_j}, v_0 = true input).  Digits stay in sync by
#   construction (they are defined by the tracked m).  Unrolling
#   log2 v_{j+1} = 2 log2 v_j - d_j gives
#   log2 v_0 = sum d_j 2^{-(j+1)} + 2^{-B} log2 v_B, and v_B >= y_B >= 1
#   so log2 v_B >= 0: acc/2^B <= log2(input).
#   UP pass: ceil shifts give y_j >= v_j and y stays <= 2 (squaring a
#   value <= 2 then halving on digit 1 keeps it <= 2; ceil adds < 1 ulp),
#   so v_B <= y_B <= 2 and log2 v_B <= 1: log2(input) <= (acc+1)/2^B.
# The initial mantissa is floor/ceil of q * 2^(WB-e), which only widens
# the bracket in the safe direction (log2 is increasing).
# ======================================================================

WB = 80    # working mantissa bits
BB = 48    # extracted fraction bits


def _digits_dir(m: int, up: bool, wb: int = WB, bits: int = BB) -> int:
    """Truncated binary expansion of log2(m/2^wb), m/2^wb in [1,2];
    floor variant if not up (acc/2^bits <= log2), ceil variant if up
    (log2 <= (acc+1)/2^bits)."""
    acc = 0
    two = 1 << (2 * wb + 1)          # value 2.0 at the squared scale
    for _ in range(bits):
        m *= m
        acc <<= 1
        if m >= two:
            acc |= 1
            sh = wb + 1
        else:
            sh = wb
        if up:
            m = -((-m) >> sh)         # ceil
        else:
            m >>= sh                  # floor
    return acc


def log2_bracket(q: Fraction) -> tuple[Fraction, Fraction]:
    """[lo, hi] dyadic with lo <= log2 q <= hi, width ~ 3*2^-BB; q > 0.

    Safety margins: the tracked mantissa can drift below 1 (floor pass)
    or above 2 (ceil pass) by < 2^-(WB-2), so log2 v_B lies in
    [-2^-(WB-1), 1 + 2^-(WB-1)] instead of [0, 1]; the -1/+2 ulps on the
    extracted digits absorb that with room to spare."""
    assert q > 0
    num, den = q.numerator, q.denominator
    e = num.bit_length() - den.bit_length()
    # ensure 2^e <= q < 2^(e+1)
    ge = num >= (den << e) if e >= 0 else (num << (-e)) >= den
    if not ge:
        e -= 1
    # mantissa q/2^e in [1,2): integer bracket at WB bits
    sh = WB - e
    if sh >= 0:
        m_lo = (num << sh) // den
        m_hi = -((-(num << sh)) // den)
    else:
        m_lo = num // (den << (-sh))
        m_hi = -((-num) // (den << (-sh)))
    lo = Fraction(_digits_dir(m_lo, False) - 1, 1 << BB) + e
    hi = Fraction(_digits_dir(m_hi, True) + 2, 1 << BB) + e
    return lo, hi


# ======================================================================
# my direct full-support / sparse census evaluator (float), from the
# 007 SS1 definition, written from scratch: N_i = reachable prefixes with
# both continuations reachable in supp(u); M_i = mass-weighted mean of
# log2 OR over N_i x N_i; A = sum_{i<n} w_i (M_i - lam) / sum w_i.
# ======================================================================

def my_eval_float(u: dict[int, float], lam: float, n: int):
    """Returns (aggregate, per_i list of (M_i, w_i) or None, max marginal,
    min excess over i<n, argmin i)."""
    t = 2.0 ** lam
    atoms = sorted(u)
    wts = [u[a] for a in atoms]
    m = len(atoms)
    kern = [[wts[x] * wts[y] * t ** popc(atoms[x] & atoms[y])
             for y in range(m)] for x in range(m)]
    # marginals of mu
    rowsum = [sum(kern[x]) for x in range(m)]
    Z = sum(rowsum)
    marg = [sum(rowsum[x] for x in range(m) if atoms[x] >> j & 1) / Z
            for j in range(n)]
    per_i = []
    num = den = 0.0
    worst = worst_i = None
    for i in range(1, n + 1):
        mask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        byp: dict[int, list[int]] = {}
        for x in range(m):
            byp.setdefault(atoms[x] & mask, []).append(x)
        act = [p for p, xs in byp.items()
               if any(atoms[x] & bit for x in xs)
               and any(not atoms[x] & bit for x in xs)]
        rows_w = 0.0
        rows_s = 0.0
        any_row = False
        for pa in act:
            for pb in act:
                c = [0.0, 0.0, 0.0, 0.0]
                for x in byp[pa]:
                    ax = 2 if atoms[x] & bit else 0
                    kx = kern[x]
                    for y in byp[pb]:
                        c[ax + (1 if atoms[y] & bit else 0)] += kx[y]
                if min(c) <= 0.0:
                    raise RuntimeError(f"underflow/degeneracy i={i}")
                mass = c[0] + c[1] + c[2] + c[3]
                v = math.log2(c[3] * c[0] / (c[2] * c[1]))
                rows_w += mass
                rows_s += mass * v
                any_row = True
        if not any_row:
            per_i.append(None)
            continue
        mi = rows_s / rows_w
        per_i.append((mi, rows_w))
        if i < n:
            num += rows_w * (mi - lam)
            den += rows_w
            if worst is None or mi - lam < worst:
                worst, worst_i = mi - lam, i
    agg = num / den if den > 0 else None
    return agg, per_i, max(marg), worst, worst_i


# ======================================================================
# my exact sparse certifier: integer-scaled tables + log2_bracket.
# t = tn/td rational; u values Fractions -> common-denominator integers.
# All kernel entries are scaled by td^n and the u common denominator:
# one global scale, so masses are exact ints and A's denominator exact.
# ======================================================================

def my_certify_sparse(u: dict[int, Fraction], tn: int, td: int, n: int):
    atoms = sorted(u)
    m = len(atoms)
    cd = 1
    for a in atoms:
        cd = cd * u[a].denominator // math.gcd(cd, u[a].denominator)
    ui = [int(u[a] * cd) for a in atoms]
    tp = [tn ** k * td ** (n - k) for k in range(n + 1)]
    kern = [[ui[x] * ui[y] * tp[popc(atoms[x] & atoms[y])]
             for y in range(m)] for x in range(m)]
    # exact marginals vs RECORD
    rowsum = [sum(kern[x]) for x in range(m)]
    Z = sum(rowsum)
    margs = [Fraction(sum(rowsum[x] for x in range(m)
                          if atoms[x] >> j & 1), Z) for j in range(n)]
    mm = max(margs)
    num_lo = num_hi = Fraction(0)
    den = 0
    memo: dict[tuple[int, int], tuple[Fraction, Fraction]] = {}
    for i in range(1, n):                       # i = n excluded
        mask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        byp: dict[int, list[int]] = {}
        for x in range(m):
            byp.setdefault(atoms[x] & mask, []).append(x)
        act = [p for p, xs in byp.items()
               if any(atoms[x] & bit for x in xs)
               and any(not atoms[x] & bit for x in xs)]
        for pa in act:
            for pb in act:
                c = [0, 0, 0, 0]
                for x in byp[pa]:
                    ax = 2 if atoms[x] & bit else 0
                    kx = kern[x]
                    for y in byp[pb]:
                        c[ax + (1 if atoms[y] & bit else 0)] += kx[y]
                assert min(c) > 0, f"degenerate cell in active pair i={i}"
                mass = c[0] + c[1] + c[2] + c[3]
                # log2(OR/t) = log2( (c11 c00 td) / (c01 c10 tn) )
                qn, qd = c[3] * c[0] * td, c[2] * c[1] * tn
                g = math.gcd(qn, qd)
                key = (qn // g, qd // g)
                enc = memo.get(key)
                if enc is None:
                    enc = log2_bracket(Fraction(key[0], key[1]))
                    memo[key] = enc
                num_lo += mass * enc[0]
                num_hi += mass * enc[1]
                den += mass
    assert den > 0
    return num_lo / den, num_hi / den, mm


# ======================================================================
# my orbit evaluator for tail-exchangeable measures (own re-derivation).
#
# Atom state s = (a1, a2, k), k = |A cap {3..n}|, r = n - 2 tail coords,
# potential g by state.  For coordinate i >= 3 (response = tail coord i):
# q = i - 3 prefix tail coords, f = n - i future tail coords.  History
# class of a prefix pair: (oa, ob) in {0,1}^2 x {0,1}^2 (special coords),
# (ja, jb) prefix tail counts, j = |prefix tails of a cap of b|.  The
# 2x2 table of (a, b) is
#     cell(al, be) = t^{oa.ob + j} * t^{al.be} * S[ja+al][jb+be],
#     S[x][y] = sum_{s,l} g(oa, x+s) PK_f[s][l] g(ob, y+l),
#     PK_f[s][l] = C(f,s) sum_w C(s,w) C(f-s,l-w) t^w
#       (= sum over ordered pairs of future sets, |XA| = s, |XB| = l, of
#        t^{|XA cap XB|}),
# so OR = t * S[ja][jb] S[ja+1][jb+1] / (S[ja][jb+1] S[ja+1][jb]) is
# j-independent, and the class mass sums the j-factor into
#     JW[ja][jb] = sum_j Ncl(q, ja, jb, j) t^j,
#     Ncl = C(q,j) C(q-j, ja-j) C(q-ja, jb-j)  (ordered prefix pairs).
# i = 1, 2 are handled directly over the <= 4 special-coordinate
# histories with the full r-tail pair kernel.  Membership in N_i is
# structural: pattern (o, c) is active iff both continuations c, c+1
# reach supp(g), i.e. exists s <= f with g(o, c + al + s) > 0.
#
# Exact mode: g values as integers over a common denominator 2^G (or cd),
# t-powers as tn^w td^(f-w): every S entry, JW and mass is an exact
# INTEGER on a per-i scale cd^2 * td^(n-3+extras); the extras are aligned
# across i so the aggregate denominator is one exact integer.
# ======================================================================

def _pairkernel_int(f: int, tn: int, td: int):
    """PK_f[s][l] scaled by td^f (exact ints)."""
    tp = [tn ** w * td ** (f - w) for w in range(f + 1)]
    PK = [[0] * (f + 1) for _ in range(f + 1)]
    for s in range(f + 1):
        cfs = math.comb(f, s)
        for l in range(f + 1):
            acc = 0
            for w in range(max(0, s + l - f), min(s, l) + 1):
                acc += math.comb(s, w) * math.comb(f - s, l - w) * tp[w]
            PK[s][l] = cfs * acc
    return PK


def _jw_int(q: int, tn: int, td: int):
    """JW[ja][jb] scaled by td^q (exact ints)."""
    tp = [tn ** j * td ** (q - j) for j in range(q + 1)]
    JW = [[0] * (q + 1) for _ in range(q + 1)]
    for ja in range(q + 1):
        for jb in range(q + 1):
            acc = 0
            for j in range(max(0, ja + jb - q), min(ja, jb) + 1):
                acc += (math.comb(q, j) * math.comb(q - j, ja - j)
                        * math.comb(q - ja, jb - j)) * tp[j]
            JW[ja][jb] = acc
    return JW


ORB4 = ((0, 0), (0, 1), (1, 0), (1, 1))


def my_orbit_rows_int(n: int, gi: dict, tn: int, td: int):
    """Exact integer orbit rows.  gi: state -> int (common denominator
    folded out; it cancels in OR and is a common mass factor).  Returns
    per i in 1..n a list of (mass_int, or_num, or_den) with mass_int on
    ONE common scale across i and or_num/or_den = OR/t exactly."""
    r = n - 2

    def g(o, c):
        return gi.get((o[0], o[1], c), 0)

    out = []
    # ---- i = 1, 2: special coords.  Q[(oa, ob)][x][y] over continuation
    # counts is the full-tail pair kernel contraction.
    PKr = _pairkernel_int(r, tn, td)

    def quad(oa, ob):
        acc = 0
        for s in range(r + 1):
            ga = g(oa, s)
            if ga:
                row = PKr[s]
                acc += ga * sum(row[l] * gb for l, gb in
                                ((l, g(ob, l)) for l in range(r + 1)) if gb)
        return acc

    QQ = {(oa, ob): quad(oa, ob) for oa in ORB4 for ob in ORB4}
    # i = 1: response coord 1; history (empty, empty) if both values reach
    reach_c1 = {v: any(g((v, a2), k) for a2 in (0, 1) for k in range(r + 1))
                for v in (0, 1)}
    rows1 = []
    if reach_c1[0] and reach_c1[1]:
        # cell(al, be) = t^{al be} sum_{a2 b2} t^{a2 b2} QQ[(al,a2),(be,b2)]
        # scale: td^r (QQ) * td^2 (two explicit t powers) -> mult by td^2
        cell = {}
        for al in (0, 1):
            for be in (0, 1):
                acc = 0
                for a2 in (0, 1):
                    for b2 in (0, 1):
                        acc += (tn if a2 * b2 else td) * QQ[((al, a2),
                                                             (be, b2))]
                cell[(al, be)] = (tn if al * be else td) * acc
        mass = cell[(0, 0)] + cell[(0, 1)] + cell[(1, 0)] + cell[(1, 1)]
        rows1.append((mass, cell[(1, 1)] * cell[(0, 0)] * td,
                      cell[(1, 0)] * cell[(0, 1)] * tn))
    out.append(rows1)
    # i = 2: response coord 2; prefixes = coord-1 values
    reach_12 = {(v, w): any(g((v, w), k) for k in range(r + 1))
                for v in (0, 1) for w in (0, 1)}
    act2 = [v for v in (0, 1) if reach_12[(v, 0)] and reach_12[(v, 1)]]
    rows2 = []
    for a in act2:
        for b in act2:
            cell = {}
            for al in (0, 1):
                for be in (0, 1):
                    cell[(al, be)] = ((tn if a * b else td)
                                      * (tn if al * be else td)
                                      * QQ[((a, al), (b, be))])
            mass = (cell[(0, 0)] + cell[(0, 1)] + cell[(1, 0)]
                    + cell[(1, 1)])
            rows2.append((mass, cell[(1, 1)] * cell[(0, 0)] * td,
                          cell[(1, 0)] * cell[(0, 1)] * tn))
    out.append(rows2)
    # scale so far: i=1: td^{r+2}; i=2: td^{r+2}.  For i >= 3 the natural
    # scale is td^{q + f + 2} = td^{n - 3} (JW) * (cells: td^{f} * td^2)
    # ... one td short of the i<=2 scale (r + 2 = n).  Track exponents
    # and align at the end.
    scale12 = r + 2                    # td exponent used above
    per_i_scale = [scale12, scale12]
    for i in range(3, n + 1):
        q, f = i - 3, n - i
        PK = _pairkernel_int(f, tn, td)
        JW = _jw_int(q, tn, td)
        rows = []
        # T[o][x][s] = g(o, x+s); S per (oa, ob)
        for oa in ORB4:
            ga = [[g(oa, x + s) if x + s <= r else 0
                   for s in range(f + 1)] for x in range(q + 2)]
            if not any(any(row) for row in ga):
                continue
            # A[x][l] = sum_s ga[x][s] PK[s][l]
            Amat = [[sum(ga[x][s] * PK[s][l] for s in range(f + 1)
                         if ga[x][s]) for l in range(f + 1)]
                    for x in range(q + 2)]
            reach_a = [any(ga[c]) for c in range(q + 2)]
            for ob in ORB4:
                gb = [[g(ob, y + l) if y + l <= r else 0
                       for l in range(f + 1)] for y in range(q + 2)]
                if not any(any(row) for row in gb):
                    continue
                reach_b = [any(gb[c]) for c in range(q + 2)]
                S = [[sum(Amat[x][l] * gb[y][l] for l in range(f + 1)
                          if gb[y][l]) for y in range(q + 2)]
                     for x in range(q + 2)]
                e = oa[0] * ob[0] + oa[1] * ob[1]     # in {0, 1, 2}
                base = tn ** e * td ** (2 - e)
                for ja in range(q + 1):
                    if not (reach_a[ja] and reach_a[ja + 1]):
                        continue
                    for jb in range(q + 1):
                        if not (reach_b[jb] and reach_b[jb + 1]):
                            continue
                        s00, s01 = S[ja][jb], S[ja][jb + 1]
                        s10, s11 = S[ja + 1][jb], S[ja + 1][jb + 1]
                        if not (s00 and s01 and s10 and s11):
                            continue
                        mass = base * JW[ja][jb] * (
                            td * (s00 + s01 + s10) + tn * s11)
                        rows.append((mass, s00 * s11, s01 * s10))
        out.append(rows)
        per_i_scale.append(q + f + 3)  # td^{2-e}*tn^e (2) + JW q + PK f + 1
    # align all masses to the max scale
    mx = max(per_i_scale)
    aligned = []
    for sc, rows in zip(per_i_scale, out):
        mult = td ** (mx - sc)
        aligned.append([(m * mult, on, od) for (m, on, od) in rows]
                       if mult != 1 else rows)
    return aligned


def my_orbit_agg_float(n: int, gf: dict, lam: float):
    """Float aggregate via the integer-mode row builder run on floats is
    unsafe (overflow); instead run the same formulas in float directly."""
    t = 2.0 ** lam
    rows_by_i = my_orbit_rows_generic(n, gf, t)
    num = den = 0.0
    worst = worst_i = None
    for i, rows in enumerate(rows_by_i, start=1):
        if i >= n or not rows:
            continue
        W = sum(m for m, _o in rows)
        mi = sum(m * math.log2(o) for m, o in rows) / W
        num += W * (mi - lam)
        den += W
        if worst is None or mi - lam < worst:
            worst, worst_i = mi - lam, i
    return (num / den if den else None), worst, worst_i


def my_orbit_rows_generic(n: int, g: dict, t: float):
    """Float rows (mass, OR).  Same hand-derived formulas, float t."""
    r = n - 2

    def gv(o, c):
        return g.get((o[0], o[1], c), 0.0)

    def pk(f):
        return [[math.comb(f, s) * sum(
            math.comb(s, w) * math.comb(f - s, l - w) * t ** w
            for w in range(max(0, s + l - f), min(s, l) + 1))
            for l in range(f + 1)] for s in range(f + 1)]

    out = []
    PKr = pk(r)

    def quad(oa, ob):
        acc = 0.0
        for s in range(r + 1):
            ga = gv(oa, s)
            if ga:
                acc += ga * sum(PKr[s][l] * gv(ob, l)
                                for l in range(r + 1) if gv(ob, l))
        return acc

    QQ = {(oa, ob): quad(oa, ob) for oa in ORB4 for ob in ORB4}
    reach_c1 = {v: any(gv((v, a2), k) for a2 in (0, 1)
                       for k in range(r + 1)) for v in (0, 1)}
    rows1 = []
    if reach_c1[0] and reach_c1[1]:
        cell = {}
        for al in (0, 1):
            for be in (0, 1):
                cell[(al, be)] = t ** (al * be) * sum(
                    t ** (a2 * b2) * QQ[((al, a2), (be, b2))]
                    for a2 in (0, 1) for b2 in (0, 1))
        mass = sum(cell.values())
        rows1.append((mass, cell[(1, 1)] * cell[(0, 0)]
                      / (cell[(1, 0)] * cell[(0, 1)])))
    out.append(rows1)
    reach_12 = {(v, w): any(gv((v, w), k) for k in range(r + 1))
                for v in (0, 1) for w in (0, 1)}
    act2 = [v for v in (0, 1) if reach_12[(v, 0)] and reach_12[(v, 1)]]
    rows2 = []
    for a in act2:
        for b in act2:
            cell = {(al, be): t ** (a * b + al * be) * QQ[((a, al), (b, be))]
                    for al in (0, 1) for be in (0, 1)}
            rows2.append((sum(cell.values()),
                          cell[(1, 1)] * cell[(0, 0)]
                          / (cell[(1, 0)] * cell[(0, 1)])))
    out.append(rows2)
    for i in range(3, n + 1):
        q, f = i - 3, n - i
        PK = pk(f)
        tpj = [t ** j for j in range(q + 1)]
        JW = [[sum(math.comb(q, j) * math.comb(q - j, ja - j)
                   * math.comb(q - ja, jb - j) * tpj[j]
                   for j in range(max(0, ja + jb - q), min(ja, jb) + 1))
               for jb in range(q + 1)] for ja in range(q + 1)]
        rows = []
        for oa in ORB4:
            ga = [[gv(oa, x + s) if x + s <= r else 0.0
                   for s in range(f + 1)] for x in range(q + 2)]
            if not any(any(row) for row in ga):
                continue
            Amat = [[sum(ga[x][s] * PK[s][l] for s in range(f + 1)
                         if ga[x][s]) for l in range(f + 1)]
                    for x in range(q + 2)]
            reach_a = [any(ga[c]) for c in range(q + 2)]
            for ob in ORB4:
                gb = [[gv(ob, y + l) if y + l <= r else 0.0
                       for l in range(f + 1)] for y in range(q + 2)]
                if not any(any(row) for row in gb):
                    continue
                reach_b = [any(gb[c]) for c in range(q + 2)]
                S = [[sum(Amat[x][l] * gb[y][l] for l in range(f + 1)
                          if gb[y][l]) for y in range(q + 2)]
                     for x in range(q + 2)]
                base = t ** (oa[0] * ob[0] + oa[1] * ob[1])
                for ja in range(q + 1):
                    if not (reach_a[ja] and reach_a[ja + 1]):
                        continue
                    for jb in range(q + 1):
                        if not (reach_b[jb] and reach_b[jb + 1]):
                            continue
                        s00, s01 = S[ja][jb], S[ja][jb + 1]
                        s10, s11 = S[ja + 1][jb], S[ja + 1][jb + 1]
                        if 0.0 in (s00, s01, s10, s11):
                            continue
                        rows.append((base * JW[ja][jb]
                                     * (s00 + s01 + s10 + t * s11),
                                     t * s00 * s11 / (s01 * s10)))
        out.append(rows)
    return out


def my_orbit_marginals(n: int, g: dict, t: float):
    r = n - 2
    states = [s for s in g if g[s] > 0]

    def wt(k, l):
        return sum(math.comb(k, j) * math.comb(r - k, l - j) * t ** j
                   for j in range(max(0, k + l - r), min(k, l) + 1))

    z = m1 = m2 = mt = 0.0
    for s in states:
        row = sum(g[s2] * t ** (s[0] * s2[0] + s[1] * s2[1])
                  * wt(s[2], s2[2]) for s2 in states)
        w = math.comb(r, s[2]) * g[s] * row
        z += w
        m1 += w * s[0]
        m2 += w * s[1]
        mt += w * s[2] / r
    return m1 / z, m2 / z, mt / z


def my_orbit_sinkhorn(n: int, tgt: dict, lam: float, iters: int = 40000,
                      tol: float = 1e-11):
    """My scaling iteration: g <- g * (tgt*z/cur)^0.4 (different damping
    exponent from 014's 0.5)."""
    r = n - 2
    t = 2.0 ** lam
    states = [s for s in sorted(tgt) if tgt[s] > 0]

    def wt(k, l):
        return sum(math.comb(k, j) * math.comb(r - k, l - j) * t ** j
                   for j in range(max(0, k + l - r), min(k, l) + 1))

    M = [[t ** (a[0] * b[0] + a[1] * b[1]) * wt(a[2], b[2])
          for b in states] for a in states]
    cnt = [math.comb(r, s[2]) for s in states]
    tv = [tgt[s] for s in states]
    N = len(states)
    g = [tv[a] ** 0.5 for a in range(N)]
    resid = float("inf")
    for _ in range(iters):
        cur = [g[a] * sum(g[b] * M[a][b] for b in range(N))
               for a in range(N)]
        z = sum(cnt[a] * cur[a] for a in range(N))
        resid = max(abs(cur[a] / (z * tv[a]) - 1.0) for a in range(N))
        if resid < tol:
            break
        g = [g[a] * (tv[a] * z / cur[a]) ** 0.4 for a in range(N)]
    return {s: g[a] for a, s in enumerate(states)}, resid


def my_certify_orbit_rows(n: int, gF: dict, tn: int, td: int):
    """Exact aggregate enclosure from integer orbit rows + log2_bracket.
    gF: state -> Fraction (dyadic).  Returns (lo, hi, den>0)."""
    cd = 1
    for v in gF.values():
        cd = cd * v.denominator // math.gcd(cd, v.denominator)
    gi = {s: int(v * cd) for s, v in gF.items() if v != 0}
    rows_by_i = my_orbit_rows_int(n, gi, tn, td)
    num_lo = num_hi = Fraction(0)
    den = 0
    memo: dict[tuple[int, int], tuple[Fraction, Fraction]] = {}
    for i, rows in enumerate(rows_by_i, start=1):
        if i >= n:
            continue
        for (mass, on, od) in rows:
            # OR/t = (on*td.. already exact: on/od * ... ) -- rows carry
            # or_num/or_den = OR/t exactly (see builder)
            g2 = math.gcd(on, od)
            key = (on // g2, od // g2)
            enc = memo.get(key)
            if enc is None:
                enc = log2_bracket(Fraction(key[0], key[1]))
                memo[key] = enc
            num_lo += mass * enc[0]
            num_hi += mass * enc[1]
            den += mass
    assert den > 0
    return num_lo / den, num_hi / den


def my_orbit_marginals_exact(n: int, gF: dict, tn: int, td: int):
    r = n - 2
    states = [s for s in gF if gF[s] > 0]
    cd = 1
    for v in gF.values():
        cd = cd * v.denominator // math.gcd(cd, v.denominator)
    gi = {s: int(gF[s] * cd) for s in states}

    def wt(k, l):
        return sum(math.comb(k, j) * math.comb(r - k, l - j)
                   * tn ** j * td ** (r - j)
                   for j in range(max(0, k + l - r), min(k, l) + 1))

    z = m1 = m2 = Fraction(0)
    mt = Fraction(0)
    for a in states:
        row = sum(gi[b] * tn ** (a[0] * b[0] + a[1] * b[1])
                  * td ** (2 - a[0] * b[0] - a[1] * b[1])
                  * wt(a[2], b[2]) for b in states)
        w = math.comb(r, a[2]) * gi[a] * row
        z += w
        m1 += w * a[0]
        m2 += w * a[1]
        mt += Fraction(w * a[2], r)
    return m1 / z, m2 / z, mt / z


# ======================================================================
# family targets (own code, same definitions as the record states them)
# ======================================================================

def orbit_states(n: int):
    return [(a1, a2, k) for a1 in (0, 1) for a2 in (0, 1)
            for k in range(n - 1)]


def st_bern(n: int, p: float) -> dict:
    return {(a1, a2, k): p ** (a1 + a2 + k) * (1 - p) ** (n - a1 - a2 - k)
            for a1 in (0, 1) for a2 in (0, 1) for k in range(n - 1)}


def st_mixture(n: int, comps) -> dict:
    out = {}
    for w, p in comps:
        for s, v in st_bern(n, p).items():
            out[s] = out.get(s, 0.0) + w * v
    return {s: v for s, v in out.items() if v > 0}


def st_slice(n: int, w0: int) -> dict:
    tot = math.comb(n, w0)
    return {(a1, a2, k): 1.0 / tot for a1 in (0, 1) for a2 in (0, 1)
            for k in range(n - 1) if a1 + a2 + k == w0}


def st_crashmix(n: int, p: float, d: float) -> dict:
    r = n - 2
    out = {s: (1 - d) * v for s, v in st_bern(n, p).items()}
    for s, w in (((0, 1, 0), 0.33), ((0, 0, r), 0.33), ((1, 1, r), 0.04),
                 ((1, 0, 0), 0.30)):
        out[s] = out.get(s, 0.0) + d * w
    return out


def st_mix_dir(n: int, p: float, nu: str, d: float) -> dict:
    out = {s: (1 - d) * v for s, v in st_bern(n, p).items()}
    r = n - 2
    if nu == "crash":
        for s, w in (((0, 1, 0), 0.33), ((0, 0, r), 0.33),
                     ((1, 1, r), 0.04), ((1, 0, 0), 0.30)):
            out[s] = out.get(s, 0.0) + d * w
    elif nu == "slice":
        for s, v in st_slice(n, int(RECORD_F * n)).items():
            out[s] = out.get(s, 0.0) + d * v
    elif nu == "slice_matched":
        for s, v in st_slice(n, round(p * n)).items():
            out[s] = out.get(s, 0.0) + d * v
    else:
        raise ValueError(nu)
    return out


def st_d0bern(n: int, pl: float, ph: float) -> dict:
    out = {s: (1 - pl) * v for s, v in st_bern(n, ph).items()}
    out[(0, 0, 0)] = out.get((0, 0, 0), 0.0) + pl
    return out


def orbit_to_full(n: int, g: dict) -> dict[int, float]:
    """Expand an orbit potential to an atom potential on 2^[n]."""
    u = {}
    for A in range(1 << n):
        s = (A & 1, (A >> 1) & 1, popc(A >> 2))
        v = g.get(s, 0.0)
        if v > 0:
            u[A] = v
    return u


# 007's 10-atom witness weights (the measure definition, data):
W_LIGHT = 0.004792968160237178
W_ATOMS = (3.802705744653898, 8.094727225702503, W_LIGHT,
           0.2526660836353666, 27.28641199464732, 12.529780963546393,
           0.3609234435898425)


def my_witness_at_n(n: int, wdil: float) -> dict[int, float]:
    """Own implementation of the record's witness-at-n layout: marker
    coord 1, dilution singletons 2..n-3, response n-2, futures n-1, n."""
    bm, br = 1, 1 << (n - 3)
    f1, f2 = 1 << (n - 2), 1 << (n - 1)
    a = W_ATOMS
    u = {bm | f1: a[0], bm | f2: a[1], bm | br | f2: a[2], 0: a[3],
         f1: a[4], f2: a[5], br | f1: a[6]}
    for j in range(2, n - 2):
        u[1 << (j - 1)] = wdil
    return u


# ======================================================================
# part K -- attack the enclosure kit + combinatorics
# ======================================================================

def part_k() -> dict:
    print("== part K: enclosure-kit attack ==")
    import uc_or_agg_probe as probe
    out = {}
    import random
    rng = random.Random(20260804)

    # K1: their log2_enclosure vs my log2_bracket on adversarial rationals
    cases = []
    for _ in range(400):
        num = rng.randrange(1, 1 << rng.randrange(4, 200))
        den = rng.randrange(1, 1 << rng.randrange(4, 200))
        cases.append(Fraction(num, den))
    # near powers of two (both sides), tiny, huge, exact powers
    for k in (-100, -7, -1, 0, 1, 7, 100):
        for eps_num, eps_den in ((0, 1), (1, 10 ** 25), (-1, 10 ** 25),
                                 (1, 3), (-1, 3)):
            q = Fraction(2) ** k * (1 + Fraction(eps_num, eps_den))
            if q > 0:
                cases.append(q)
    cases.append(Fraction(181, 16))
    cases.append(Fraction(271, 256))
    bad = []
    maxw_theirs = maxw_mine = 0.0
    for q in cases:
        tlo, thi = probe.log2_enclosure(q)
        mlo, mhi = log2_bracket(q)
        # both must contain the true value: cross-containment demands
        # overlap, and each interval must contain the other's midpoint
        # up to widths; the decisive check: tlo <= mhi and mlo <= thi
        if not (tlo <= mhi and mlo <= thi):
            bad.append((str(q), float(tlo), float(thi), float(mlo),
                        float(mhi)))
        maxw_theirs = max(maxw_theirs, float(thi - tlo))
        maxw_mine = max(maxw_mine, float(mhi - mlo))
        # sanity anchor: floor consistency with float log2 where safe
        fl = math.log2(float(q)) if 1e-300 < float(q) < 1e300 else None
        if fl is not None and (fl < float(mlo) - 1e-6
                               or fl > float(mhi) + 1e-6):
            bad.append(("float-anchor", str(q), fl, float(mlo), float(mhi)))
    print(f"[K1] {len(cases)} rationals: disjoint-enclosure failures = "
          f"{len(bad)}; max widths theirs {maxw_theirs:.1e} / mine "
          f"{maxw_mine:.1e}")
    out["k1_cases"] = len(cases)
    out["k1_failures"] = bad
    out["k1_pass"] = not bad

    # K2: their _dround directional soundness, random Fractions inc. negatives
    dbad = 0
    for _ in range(2000):
        x = Fraction(rng.randrange(-10 ** 30, 10 ** 30),
                     rng.randrange(1, 10 ** 30))
        lo = probe._dround(x, False)
        hi = probe._dround(x, True)
        if not (lo <= x <= hi):
            dbad += 1
    print(f"[K2] _dround directionality: {dbad} failures / 2000")
    out["k2_failures"] = dbad
    out["k2_pass"] = dbad == 0

    # K3: their _ln_1_to_2_fp directed bounds -- sanity anchor against my
    # bracket via ln s = log2(s) * ln 2 (float ln 2; the RIGOROUS
    # end-to-end test of their composed log2_enclosure is K1)
    lbad = 0
    for _ in range(300):
        s = 1 + Fraction(rng.randrange(0, 10 ** 18), 10 ** 18)
        if s >= 2 or s == 1:
            continue
        lo = probe._ln_1_to_2_fp(s, False, 12)
        hi = probe._ln_1_to_2_fp(s, True, 12)
        l2lo, l2hi = log2_bracket(s)
        mid = (float(l2lo) + float(l2hi)) / 2 * math.log(2.0)
        if not (float(lo) - 1e-12 <= mid <= float(hi) + 1e-12):
            lbad += 1
    print(f"[K3] _ln_1_to_2_fp anchor vs my bracket*ln2 (slack 1e-12): "
          f"{lbad} failures / 300")
    out["k3_failures"] = lbad
    out["k3_pass"] = lbad == 0

    # K4: brute-force the orbit combinatorics my derivation predicts
    #  Ncl(q,ja,jb,j) counts ordered prefix pairs; PK_f pair kernel
    ncl_bad = pk_bad = 0
    for q in range(0, 7):
        for ja in range(q + 1):
            for jb in range(q + 1):
                cnt = {}
                for a in range(1 << q):
                    if popc(a) != ja:
                        continue
                    for b in range(1 << q):
                        if popc(b) != jb:
                            continue
                        cnt[popc(a & b)] = cnt.get(popc(a & b), 0) + 1
                for j in range(q + 1):
                    pred = (math.comb(q, j) * math.comb(q - j, ja - j)
                            * math.comb(q - ja, jb - j)
                            if j <= min(ja, jb) and ja + jb - j <= q else 0)
                    if cnt.get(j, 0) != pred:
                        ncl_bad += 1
    t_test = Fraction(7, 3)
    for f in range(0, 6):
        PK = _pairkernel_int(f, 7, 3)
        for s in range(f + 1):
            for l in range(f + 1):
                acc = Fraction(0)
                for XA in range(1 << f):
                    if popc(XA) != s:
                        continue
                    for XB in range(1 << f):
                        if popc(XB) != l:
                            continue
                        acc += t_test ** popc(XA & XB)
                if acc * 3 ** f != PK[s][l]:
                    pk_bad += 1
    print(f"[K4] combinatorics brute force: Ncl failures {ncl_bad}, "
          f"PK failures {pk_bad}")
    out["k4_pass"] = ncl_bad == 0 and pk_bad == 0

    save("oragg_sk2_partK.json", out)
    return out


# ======================================================================
# part X -- formula cross-checks (n = 8, 9; j-independence; witness)
# ======================================================================

def part_x() -> dict:
    print("== part X: formula cross-checks ==")
    import uc_or_agg_probe as probe
    out = {"checks": []}

    # X1: my direct evaluator on the witness at n = 7 (007 numbers)
    u7 = my_witness_at_n(7, 20.0)
    agg, per_i, mm, worst, wi = my_eval_float(u7, 3.5, 7)
    print(f"[X1] my direct on witness n=7 lam=3.5: M_5-lam = "
          f"{per_i[4][0] - 3.5:+.6f} (007: -0.122033), agg = {agg:+.6f} "
          f"(014: +1.844754), marg = {mm:.4f} (0.3178)")
    ok = (abs(per_i[4][0] - 3.5 + 0.122033) < 5e-6
          and abs(agg - 1.844754) < 5e-5 and abs(mm - 0.317769) < 5e-6)
    out["checks"].append({"name": "witness_n7", "pass": bool(ok),
                          "agg": agg, "M5_minus_lam": per_i[4][0] - 3.5})

    # X2: my orbit vs my direct vs 014's orbit at n = 8, 9
    fams = [("crashmix_383_05", lambda n: st_crashmix(n, 0.383, 0.05)),
            ("sawin_geo_a", lambda n: st_mixture(
                n, [(0.80, 0.28), (0.12, 0.55), (0.06, 0.75),
                    (0.02, 0.92)])),
            ("mixdir_slice_p.30_d.05", lambda n: st_mix_dir(
                n, 0.30, "slice", 0.05)),
            ("d0bern_30_51", lambda n: st_d0bern(n, 0.30, 0.51))]
    worst_dd = 0.0
    for n in (8, 9):
        for lam in (0.5, 1.0):
            for name, mk in fams:
                tgt = mk(n)
                g, res = my_orbit_sinkhorn(n, tgt, lam)
                agg_o, w_o, _ = my_orbit_agg_float(n, g, lam)
                # direct on the SAME potential expanded to atoms
                u_full = orbit_to_full(n, g)
                agg_d, per_d, mm_d, w_d, _ = my_eval_float(u_full, lam, n)
                # 014's orbit rows on the same g
                rows_p = probe.orbit_rows(n, 2.0 ** lam, g)
                agg_p, w_p, _wi, _pp = probe.orbit_agg(rows_p, lam, n)
                dd = max(abs(agg_o - agg_d), abs(agg_o - agg_p))
                worst_dd = max(worst_dd, dd)
                out["checks"].append(
                    {"name": f"x2_{name}_n{n}_lam{lam}", "agg_mine_orbit":
                     agg_o, "agg_mine_direct": agg_d, "agg_014_orbit":
                     agg_p, "maxdiff": dd, "resid": res})
                print(f"[X2] {name} n={n} lam={lam}: mine-orbit "
                      f"{agg_o:+.9f} mine-direct {agg_d:+.9f} 014-orbit "
                      f"{agg_p:+.9f} (maxdiff {dd:.2e})")
    out["x2_maxdiff"] = worst_dd
    out["x2_pass"] = worst_dd < 1e-8

    # X3: explicit j-independence of class ORs at n = 9, i = 6:
    # per-history ORs from my direct evaluator grouped by orbit class
    n, lam = 9, 1.0
    g, _ = my_orbit_sinkhorn(n, st_crashmix(n, 0.383, 0.05), lam)
    u = orbit_to_full(n, g)
    t = 2.0 ** lam
    atoms = sorted(u)
    i = 6
    mask, bit = (1 << (i - 1)) - 1, 1 << (i - 1)
    byp: dict[int, list[int]] = {}
    for a in atoms:
        byp.setdefault(a & mask, []).append(a)
    act = [p for p, xs in byp.items()
           if any(x & bit for x in xs) and any(not x & bit for x in xs)]
    classes: dict[tuple, list[float]] = {}
    for pa in act:
        for pb in act:
            c = [0.0] * 4
            for A in byp[pa]:
                for B in byp[pb]:
                    c[(2 if A & bit else 0) + (1 if B & bit else 0)] += \
                        u[A] * u[B] * t ** popc(A & B)
            orr = c[3] * c[0] / (c[2] * c[1])
            key = (pa & 1, (pa >> 1) & 1, pb & 1, (pb >> 1) & 1,
                   popc(pa >> 2), popc(pb >> 2))
            classes.setdefault(key, []).append(orr)
    spread = max((max(v) - min(v)) / max(v) for v in classes.values())
    print(f"[X3] j-independence at n=9 i=6: {len(classes)} classes, max "
          f"rel spread of within-class ORs = {spread:.2e}")
    out["x3_class_or_rel_spread"] = spread
    out["x3_pass"] = spread < 1e-9

    ok_all = all(c.get("pass", True) for c in out["checks"]) and \
        out["x2_pass"] and out["x3_pass"]
    print(f"[X] {'ALL PASS' if ok_all else 'FAILURES PRESENT'}")
    save("oragg_sk2_partX.json", out)
    return out


# ======================================================================
# part C -- independent re-certification
# ======================================================================

def part_c() -> dict:
    print("== part C: independent re-certification ==")
    import uc_or_agg_probe as probe
    out = {"certs": []}
    E = json.loads((DATA / "or_agg_probe_partE.json").read_text())
    theirs = {r["name"]: r for r in E["sparse"] + E["orbit"]}

    def report(name, lo, hi, mm, ref):
        rec = {"name": name, "my_lo": float(lo), "my_hi": float(hi),
               "my_width": float(hi - lo),
               "my_sign": "positive" if lo > 0 else
               ("NEGATIVE" if hi < 0 else "indeterminate"),
               "max_marginal": float(mm),
               "in_regime_exact": bool(mm < RECORD),
               "their_lo": ref["agg_lo"], "their_hi": ref["agg_hi"]}
        # decisive comparison: intervals must overlap (mine is proven
        # sound with width ~1e-14, so non-overlap convicts theirs)
        overlap = not (float(hi) < ref["agg_lo"] - 1e-15
                       or float(lo) > ref["agg_hi"] + 1e-15)
        rec["intervals_overlap"] = bool(overlap)
        rec["sign_agrees"] = rec["my_sign"] == ref["certified_sign"]
        out["certs"].append(rec)
        print(f"[C] {name}: mine [{float(lo):+.12g}, {float(hi):+.12g}] "
              f"(w {float(hi - lo):.1e}) vs theirs [{ref['agg_lo']:+.12g},"
              f" {ref['agg_hi']:+.12g}] -> overlap={overlap} "
              f"sign={rec['my_sign']} marg={float(mm):.6f}")
        return rec

    # C1: witness n=10 wdil=16 at t = 181/16 (sparse; my own builder)
    t0 = time.time()
    u = {k: Fraction(v) for k, v in my_witness_at_n(10, 16.0).items()}
    lo, hi, mm = my_certify_sparse(u, 181, 16, 10)
    report("witness_n10_wdil16", lo, hi, mm, theirs["witness_n10_wdil16"])

    # C2: the fully rational t = 4 point (n = 20, wdil = 20, lambda = 2)
    u = {k: Fraction(v) for k, v in my_witness_at_n(20, 20.0).items()}
    lo, hi, mm = my_certify_sparse(u, 4, 1, 20)
    report("witness_n20_t4_lam2exact", lo, hi, mm,
           theirs["witness_n20_t4_lam2exact"])

    # C3: witness n=32 at t = 181/16
    u = {k: Fraction(v) for k, v in my_witness_at_n(32, 20.0).items()}
    lo, hi, mm = my_certify_sparse(u, 181, 16, 32)
    report("witness_n32_wdil20", lo, hi, mm, theirs["witness_n32_wdil20"])
    print(f"    (sparse certs {time.time() - t0:.0f}s)")

    # C4: crashmix_383_05 n=32 at t = 271/256 -- SAME mu-tilde as 014:
    # their float Sinkhorn + their 36-bit quantization (input spec),
    # my exact rows + my log2 from there on.
    t0 = time.time()
    n = 32
    lamf = math.log2(271 / 256)
    g, resid = probe.orbit_sinkhorn(n, probe.st_crashmix(n, 0.383, 0.05),
                                    lamf)
    gF = {s: probe._quant(v, 36) for s, v in g.items()}
    lo, hi = my_certify_orbit_rows(n, gF, 271, 256)
    m1, m2, mt = my_orbit_marginals_exact(n, gF, 271, 256)
    report("crashmix_383_05_n32", lo, hi, max(m1, m2, mt),
           theirs["crashmix_383_05_n32"])
    print(f"    ({time.time() - t0:.0f}s, their-resid {resid:.1e})")

    # C5: the softest point, mixdir_slice p=0.30 d=0.05 n=32, same mu-tilde
    t0 = time.time()
    tgt = probe.st_mix_dir(32, 0.30, "slice", 0.05)
    g, resid = probe.orbit_sinkhorn(32, tgt, lamf)
    gF = {s: probe._quant(v, 36) for s, v in g.items()}
    lo, hi = my_certify_orbit_rows(32, gF, 271, 256)
    m1, m2, mt = my_orbit_marginals_exact(32, gF, 271, 256)
    report("mixdir_slice_p0.3_d0.05_n32", lo, hi, max(m1, m2, mt),
           theirs["mixdir_slice_p0.3_d0.05_n32"])
    print(f"    ({time.time() - t0:.0f}s)")

    # C6: fresh certificate at the softest point from MY OWN pipeline
    # (my Sinkhorn damping, my 40-bit quantization): an independent
    # mu-tilde' near the same family point -- new evidence, not a rerun.
    t0 = time.time()
    tgt = st_mix_dir(32, 0.30, "slice", 0.05)
    g, resid = my_orbit_sinkhorn(32, tgt, lamf)

    def quant40(v: float) -> Fraction:
        mant, e = math.frexp(v)
        return Fraction(round(mant * (1 << 40)), 1 << 40) * Fraction(2) ** e

    gF2 = {s: quant40(v) for s, v in g.items()}
    lo, hi = my_certify_orbit_rows(32, gF2, 271, 256)
    m1, m2, mt = my_orbit_marginals_exact(32, gF2, 271, 256)
    mm = max(m1, m2, mt)
    rec = {"name": "mixdir_slice_p0.3_d0.05_n32_MY_mu", "my_lo": float(lo),
           "my_hi": float(hi), "my_sign": "positive" if lo > 0 else
           ("NEGATIVE" if hi < 0 else "indeterminate"),
           "max_marginal": float(mm), "in_regime_exact": bool(mm < RECORD),
           "sinkhorn_resid": resid}
    out["certs"].append(rec)
    print(f"[C] MY OWN mu-tilde' at softest point: [{float(lo):+.9e}, "
          f"{float(hi):+.9e}] -> {rec['my_sign']}, marg {float(mm):.6f}, "
          f"resid {resid:.1e} ({time.time() - t0:.0f}s)")

    save("oragg_sk2_partC.json", out)
    return out


# ======================================================================
# part T -- trend attack: the dip, off-grid lambda, lambda-scaling
# ======================================================================

def run_orbit_float(n, tgt, lam):
    g, resid = my_orbit_sinkhorn(n, tgt, lam)
    agg, worst, wi = my_orbit_agg_float(n, g, lam)
    m1, m2, mt = my_orbit_marginals(n, g, 2.0 ** lam)
    return {"n": n, "lambda": lam, "aggregate": agg, "min_excess": worst,
            "argmin_i": wi, "max_marginal": max(m1, m2, mt),
            "in_regime": max(m1, m2, mt) < RECORD_F, "resid": resid}


def part_t() -> dict:
    print("== part T: trend attack ==")
    res = []

    # T1: reproduce the dip row with my engine
    for n in (32, 40, 48, 56, 64):
        r = run_orbit_float(n, st_mix_dir(n, 0.30, "slice", 0.05), 0.084)
        r["name"] = "dip_slice_p.30_d.05"
        res.append(r)
        print(f"[T1] n={n}: agg {r['aggregate']:+.3e} (014: 6.15e-7 / "
              f"2.22e-7 / 1.41e-7 / 2.33e-7 / 4.31e-7) resid {r['resid']:.1e}")

    # T2: the dip point at off-grid lambda (below 014's grid floor and
    # around it) and other d -- does anything cross zero?
    for lam in (0.005, 0.01, 0.02, 0.04, 0.06, 0.12, 0.25, 0.5):
        r = run_orbit_float(48, st_mix_dir(48, 0.30, "slice", 0.05), lam)
        r["name"] = "dip_lambda_scan"
        res.append(r)
        print(f"[T2] n=48 lam={lam:g}: agg {r['aggregate']:+.3e} "
              f"agg/lam {r['aggregate'] / lam:+.3e}")
    for d in (0.01, 0.02, 0.1, 0.2):
        r = run_orbit_float(48, st_mix_dir(48, 0.30, "slice", d), 0.084)
        r["name"] = f"dip_d{d:g}"
        res.append(r)
        print(f"[T2] n=48 d={d:g}: agg {r['aggregate']:+.3e}")
    # neighborhood in p and layer (lead-2 spot sample)
    for p in (0.26, 0.28, 0.32, 0.34):
        r = run_orbit_float(48, st_mix_dir(48, p, "slice", 0.05), 0.084)
        r["name"] = f"dip_p{p:g}"
        res.append(r)
        print(f"[T2] n=48 p={p:g} (record layer): agg {r['aggregate']:+.3e}")
    for p in (0.28, 0.30, 0.32):
        r = run_orbit_float(48, st_mix_dir(48, p, "slice_matched", 0.05),
                            0.084)
        r["name"] = f"dipm_p{p:g}"
        res.append(r)
        print(f"[T2] n=48 p={p:g} (matched layer): agg {r['aggregate']:+.3e}")

    # T3: census-minimum family at off-grid lambda: crashmix_383_05,
    # lambda scan incl. tiny and between 014's big grid points
    for n in (10, 20, 32):
        for lam in (0.005, 0.02, 0.05, 0.084, 0.15, 0.3, 0.7, 1.5,
                    2.75, 4.4):
            r = run_orbit_float(n, st_crashmix(n, 0.383, 0.05), lam)
            r["name"] = "crashmix_383_05_scan"
            res.append(r)
        row = [x for x in res if x["name"] == "crashmix_383_05_scan"
               and x["n"] == n]
        print(f"[T3] crashmix n={n}: " + " ".join(
            f"{x['lambda']:g}:{x['aggregate']:+.2e}" for x in row))

    neg = [r for r in res if r["in_regime"] and r["aggregate"] is not None
           and r["aggregate"] < -1e-12]
    print(f"[T] negative in-regime aggregates: {len(neg)}")
    for r in neg:
        print("    NEG:", r)
    save("oragg_sk2_partT.json", {"results": res, "n_negative": len(neg)})
    return {"results": res, "n_negative": len(neg)}


# ======================================================================
# part V -- new adversaries 014 could not see
# ======================================================================

def part_v() -> dict:
    print("== part V: new-adversary hunt ==")
    res = []

    # V1: witness-DIRECTION perturbation of a dense product bulk at
    # n = 10: u = x^{|A|} (product potential, an equality point) + d *
    # (witness atoms).  Non-tail-symmetric direction; 014's part P had
    # only 4 symmetric rays.  x is set so the pure-bulk marginal is ~ p0.
    n = 10

    def bulk_x(p0: float, t: float) -> float:
        # per-coordinate: p = (x + x^2 t)/(1 + 2x + x^2 t); solve for x
        lo_, hi_ = 1e-6, 10.0
        for _ in range(200):
            mid = (lo_ + hi_) / 2
            p = (mid + mid * mid * t) / (1 + 2 * mid + mid * mid * t)
            if p < p0:
                lo_ = mid
            else:
                hi_ = mid
        return (lo_ + hi_) / 2

    for lam in (0.084, 0.5, 2.0, 3.5):
        t = 2.0 ** lam
        for p0 in (0.30, 0.375):
            x = bulk_x(p0, t)
            base = {A: x ** popc(A) for A in range(1 << n)}
            uwit = my_witness_at_n(n, 4.0)
            scale = max(base.values())
            for d in (0.03, 0.3, 3.0):
                u = dict(base)
                for a, v in uwit.items():
                    u[a] = u.get(a, 0.0) + d * scale * v / max(uwit.values())
                agg, per_i, mm, worst, wi = my_eval_float(u, lam, n)
                rec = {"name": f"dense_bulk_wit_p{p0:g}_d{d:g}",
                       "n": n, "lambda": lam, "aggregate": agg,
                       "max_marginal": mm, "in_regime": mm < RECORD_F,
                       "min_excess": worst, "argmin_i": wi}
                res.append(rec)
                print(f"[V1] lam={lam:g} p0={p0:g} d={d:g}: agg "
                      f"{agg:+.6e} marg {mm:.4f} minexc {worst:+.4f}"
                      f"{' OUT' if not rec['in_regime'] else ''}")

    # V2: boundary-of-regime marginals: crashmix with p tuned so the max
    # marginal sits just under 0.38271, at small/window lambda
    for n in (20, 32):
        for lam in (0.02, 0.084, 0.5):
            p_lo, p_hi = 0.383, 0.3845
            for _ in range(20):
                p = (p_lo + p_hi) / 2
                r = run_orbit_float(n, st_crashmix(n, p, 0.05), lam)
                if r["max_marginal"] < 0.382695:
                    p_lo = p
                else:
                    p_hi = p
            r = run_orbit_float(n, st_crashmix(n, p_lo, 0.05), lam)
            r["name"] = f"crashmix_boundary_p{p_lo:.6f}"
            res.append(r)
            print(f"[V2] n={n} lam={lam:g} p={p_lo:.6f}: agg "
                  f"{r['aggregate']:+.3e} marg {r['max_marginal']:.6f}")
        # d0bern pushed to the boundary: (1-pl) ph ~ 0.38268
        for (pl, ph) in ((0.304145, 0.55), (0.48976, 0.75)):
            for lam in (0.084, 0.5):
                r = run_orbit_float(n, st_d0bern(n, pl, ph), lam)
                r["name"] = f"d0bern_boundary_{pl:g}_{ph:g}"
                res.append(r)
                print(f"[V2] n={n} d0bern({pl:g},{ph:g}) lam={lam:g}: agg "
                      f"{r['aggregate']:+.3e} marg {r['max_marginal']:.6f}")
        # sawin mixture scaled to the boundary
        comps = [(0.80, 0.30317), (0.12, 0.59551), (0.06, 0.81206),
                 (0.02, 0.99612)]
        for lam in (0.084, 0.5):
            r = run_orbit_float(n, st_mixture(n, comps), lam)
            r["name"] = "sawin_boundary"
            res.append(r)
            print(f"[V2] n={n} sawin-boundary lam={lam:g}: agg "
                  f"{r['aggregate']:+.3e} marg {r['max_marginal']:.6f}")

    # V3: tiny lambda at the tight census points (below 014's grid floor)
    for n in (20, 32):
        for lam in (0.002, 0.01, 0.04):
            r = run_orbit_float(n, st_crashmix(n, 0.383, 0.05), lam)
            r["name"] = "crashmix_tiny_lambda"
            res.append(r)
            r2 = run_orbit_float(n, st_mix_dir(n, 0.30, "slice", 0.05), lam)
            r2["name"] = "slice_dir_tiny_lambda"
            res.append(r2)
            print(f"[V3] n={n} lam={lam:g}: crashmix {r['aggregate']:+.3e} "
                  f"slice-dir {r2['aggregate']:+.3e}")

    neg = [r for r in res if r["in_regime"] and r["aggregate"] is not None
           and r["aggregate"] < -1e-12]
    print(f"[V] in-regime negatives: {len(neg)}")
    for r in neg:
        print("    NEG:", r)
    save("oragg_sk2_partV.json", {"results": res, "n_negative": len(neg)})
    return {"results": res, "n_negative": len(neg)}


def main() -> None:
    parts = [p.upper() for p in sys.argv[1:] if p.upper() in "KXCTV"] \
        or list("KXCTV")
    t0 = time.time()
    if "K" in parts:
        part_k()
    if "X" in parts:
        part_x()
    if "C" in parts:
        part_c()
    if "T" in parts:
        part_t()
    if "V" in parts:
        part_v()
    print(f"done in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
