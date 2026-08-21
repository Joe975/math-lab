#!/usr/bin/env python3
"""Attempt 049 (reviewer pass on 044-048): independent re-implementation.

Nothing here imports from uc_hu_ownconst.py, uc_hu_bestorder.py,
uc_hu_n2.py, uc_hu_n2_dichotomy.py, uc_hu_L1.py, their skeptics, or the
uc_hu_order2 / uc_hu_certify / uc_reviewer036 stack.  Everything is
rebuilt from the records' prose:

  * the half-union coupling (031 "The coupling": per history cell,
    x = P(A_i = 0 | a), y = P(B_i = 0 | b),
    z = min(max(1/2, x + y - 1), x, y), cell split
    (z, x-z, y-z, 1-x-y+z)) and the chain-rule value
    CR = sum_i E[h(z_i)] - H(mu) (009), here as an EXACT rational
    history-pair evaluator with h() evaluated in 80-digit Decimal
    (Decimal.ln(); no series/digit-extraction kit is reused);
  * the canonical order (035: greedy min H(A_i | A_S), ties -> lowest
    index) and the rollout order (037: maximise full CR with the
    remainder completed canonically, ties -> lowest index);
  * the constant c*(p) = (h(max(1/2, 1-2p)) - h(p))/h(p) (034 as
    corrected by 046 F);
  * the n = 2 closed form, the decomposition, and branch L1 as stated
    in 046/047/048.

Checks (R* fail -> exit 1).  A check whose name starts with "CORRECTION"
asserts a statement this review found to be WRONG in the reviewed
records; it passes when the error is present as described.

  R1  own evaluator vs the n = 2 closed form of 046 A, and the cell
      weights' two marginals
  R2  the hand lemmas: N2-ONE-BAD (046), N2-ONE-ABOVE (047), N2-CONC
      (046) incl. its averaging step, and the exact decomposition
      identity of 046 E in BOTH cases -- plus the CORRECTION that
      N2-CONC's "tight iff p_a = p_b" clause is false
  R3  the c* correction of 046 F, both directions: the stated form is
      wrong above 1/4 but NEGATIVE only above 1/3 (CORRECTION), and
      every engine in explore/ computes the max(1/2, 1-2p) form
  R4  the 044 certificate: witness provenance in 037's checkpoint,
      50-digit recomputation of both certified values, containment in
      both kits' enclosures, own rollout-order derivation
  R5  045: every hu_bestorder.json endpoint re-scored by full order
      enumeration through this file's evaluator; floors, kills, sharp
      violations, the three structural observations
  R6  046: order-quantifier vacuity on an own grid, the B&B volume
      bookkeeping, part E's (**) failure rate on an own sampler
  R7  047: the dichotomy census re-run on an own sampler/evaluator,
      plus the SCOPE finding that the census is Case-A only and that
      "only B can be negative" fails outside Case A
  R8  048: the boundary identity in exact arithmetic, the ratio floor,
      the two-regime minimiser table, and the CORRECTION that P3's
      quoted tight counts are a --fast run, not the checkpoint

Usage: python uc_reviewer049_reimpl.py       (exit 0 iff all checks pass)
"""
from __future__ import annotations

import itertools
import json
import math
import random
import re
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FAILS = []
LOG2 = math.log(2.0)
getcontext().prec = 80


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  {detail}" if detail
                                                    else ""))
    if not ok:
        FAILS.append(name)


# --------------------------------------------------------------- floats
def h(p):
    """Binary entropy in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1 - p) * math.log1p(-p)) / LOG2


def zstar(p):
    return max(0.5, 1.0 - 2.0 * p)


def cstar(p):
    return (h(zstar(p)) - h(p)) / h(p)


def clip(x, y):
    """031's per-cell both-zero probability."""
    return min(max(0.5, x + y - 1.0), x, y)


# ------------------------------------------------------- exact side (Dec)
_DHALF = Decimal(1) / Decimal(2)


def h_dec(q):
    """Binary entropy in bits at a rational q, to the ambient precision.
    The reflection q -> 1-q is done exactly in Fraction arithmetic so
    that h is accurate at both ends."""
    if q <= 0 or q >= 1:
        return Decimal(0)
    a = q if q <= Fraction(1, 2) else 1 - q
    b = 1 - a
    da = Decimal(a.numerator) / Decimal(a.denominator)
    db = Decimal(b.numerator) / Decimal(b.denominator)
    return -(da * da.ln() + db * db.ln()) / Decimal(2).ln()


# ----------------------------------------------- own HU evaluator (exact)
def _cond_zero(n, mu, hist, i):
    """P(A_i = 0 | A_S = hist) with hist a dict coord -> bit, and the
    weight of that history."""
    tot = Fraction(0)
    zero = Fraction(0)
    for a, w in mu.items():
        if all(((a >> c) & 1) == b for c, b in hist.items()):
            tot += w
            if not ((a >> i) & 1):
                zero += w
    return (zero / tot if tot else Fraction(0)), tot


def hu_cr(n, mu, order, dec=True):
    """CR = sum_i E[h(z_i)] - H(mu) for the half-union coupling under
    the given revelation order (031/009).  dec=True: exact rational
    cells with h() in 80-digit Decimal.  dec=False: float throughout.
    Returns (CR, H)."""
    if dec:
        ent, half, one = h_dec, Fraction(1, 2), Fraction(1)
        m = {a: (w if isinstance(w, Fraction) else Fraction(w))
             for a, w in mu.items()}
        S = Decimal(0)
        wcast = (lambda w: Decimal(w.numerator) / Decimal(w.denominator))
    else:
        ent, half, one = (lambda q: h(q)), 0.5, 1.0
        m = {a: float(w) for a, w in mu.items()}
        S = 0.0
        wcast = (lambda w: w)
    memo = {}

    def cond_zero(hist, i):
        k = (hist, i)
        if k in memo:
            return memo[k]
        tot = zero = (Fraction(0) if dec else 0.0)
        for a, w in m.items():
            if all(((a >> c) & 1) == b for c, b in hist):
                tot += w
                if not ((a >> i) & 1):
                    zero += w
        v = (zero / tot if tot else (Fraction(0) if dec else 0.0))
        memo[k] = v
        return v

    state = {((), ()): (Fraction(1) if dec else 1.0)}
    for step, i in enumerate(order):
        nxt = {}
        pref = order[:step]
        for (ha, hb), w in state.items():
            x = cond_zero(tuple(zip(pref, ha)), i)
            y = cond_zero(tuple(zip(pref, hb)), i)
            z = min(max(half, x + y - one), x, y)
            S += ent(z) * wcast(w)
            for (ba, bb), cw in (((0, 0), z), ((0, 1), x - z),
                                 ((1, 0), y - z),
                                 ((1, 1), one - x - y + z)):
                if cw <= 0:
                    continue
                k = (ha + (ba,), hb + (bb,))
                nxt[k] = nxt.get(k, 0) + w * cw
        state = nxt
    H = Decimal(0) if dec else 0.0
    for w in m.values():
        if w > 0:
            if dec:
                dw = wcast(w)
                H += -dw * dw.ln() / Decimal(2).ln()
            else:
                H += -w * math.log(w) / LOG2
    return S - H, H


def order_canon(n, mu):
    """035: greedily pick the unrevealed coordinate minimising
    H(A_i | A_S); ties -> lowest index."""
    chosen = []
    while len(chosen) < n:
        best = None
        for i in range(n):
            if i in chosen:
                continue
            ce = 0.0
            for bits in itertools.product((0, 1), repeat=len(chosen)):
                hist = dict(zip(chosen, bits))
                x, w = _cond_zero(n, mu, hist, i)
                if w > 0:
                    ce += float(w) * h(float(x))
            if best is None or ce < best[0] - 1e-12:
                best = (ce, i)
        chosen.append(best[1])
    return tuple(chosen)


def canon_complete(n, mu, prefix):
    """Canonical completion of a prefix (same greedy rule, restricted)."""
    chosen = list(prefix)
    while len(chosen) < n:
        best = None
        for i in range(n):
            if i in chosen:
                continue
            ce = 0.0
            for bits in itertools.product((0, 1), repeat=len(chosen)):
                hist = dict(zip(chosen, bits))
                x, w = _cond_zero(n, mu, hist, i)
                if w > 0:
                    ce += float(w) * h(float(x))
            if best is None or ce < best[0] - 1e-12:
                best = (ce, i)
        chosen.append(best[1])
    return tuple(chosen)


def order_roll(n, mu, dec=False):
    """037: at each step pick the coordinate maximising the FULL CR of
    the canonical completion; ties (1e-12) -> lowest index."""
    chosen = []
    while len(chosen) < n:
        scores = {}
        for i in range(n):
            if i in chosen:
                continue
            full = canon_complete(n, mu, chosen + [i])
            cr, _ = hu_cr(n, mu, full, dec=dec)
            scores[i] = float(cr)
        hi = max(scores.values())
        chosen.append(min(i for i, v in scores.items() if hi - v <= 1e-12))
    return tuple(chosen)


def marginals(n, mu):
    return [float(sum(w for a, w in mu.items() if (a >> i) & 1))
            for i in range(n)]


def as_frac(mu):
    return {a: (w if isinstance(w, Fraction) else Fraction(w))
            for a, w in mu.items()}


# ============================================================== R1
def r1():
    print("\nR1. Own evaluator vs 046 A's n = 2 closed form.")
    rng = random.Random(4901)
    worst = 0.0
    worstH = 0.0
    seen = 0
    marg_bad = 0
    while seen < 400:
        x = rng.uniform(0.5, 1.0)
        u0, u1 = rng.uniform(0, 1), rng.uniform(0, 1)
        f0, f1 = 1 - x, x * (1 - u0) + (1 - x) * (1 - u1)
        if not (0 < max(f0, f1) < 0.5):
            continue
        seen += 1
        xq = Fraction(x).limit_denominator(10 ** 6)
        u0q = Fraction(u0).limit_denominator(10 ** 6)
        u1q = Fraction(u1).limit_denominator(10 ** 6)
        mu = {0b00: xq * u0q, 0b10: xq * (1 - u0q),
              0b01: (1 - xq) * u1q, 0b11: (1 - xq) * (1 - u1q)}
        mu = {a: w for a, w in mu.items() if w > 0}
        cr, H = hu_cr(2, mu, (0, 1), dec=False)
        # the record's closed form
        xf, u0f, u1f = float(xq), float(u0q), float(u1q)
        z0 = clip(xf, xf)
        w = {(0, 0): z0, (0, 1): xf - z0, (1, 0): xf - z0,
             (1, 1): 1 - 2 * xf + z0}
        uu = (u0f, u1f)
        Sc = h(z0) + sum(wt * h(clip(uu[a], uu[b]))
                         for (a, b), wt in w.items() if wt > 0)
        Hc = h(xf) + xf * h(u0f) + (1 - xf) * h(u1f)
        worst = max(worst, abs((Sc - Hc) - cr))
        worstH = max(worstH, abs(Hc - H))
        # the cell weights have BOTH marginals (x, 1-x)  [046 E averaging]
        r0 = w[(0, 0)] + w[(0, 1)]
        c0 = w[(0, 0)] + w[(1, 0)]
        if abs(r0 - xf) > 1e-12 or abs(c0 - xf) > 1e-12:
            marg_bad += 1
    check("046 A: closed form reproduces an independent HU evaluator "
          f"on {seen} in-regime samples", worst < 1e-12 and worstH < 1e-12,
          f"max|dCR| {worst:.2e}, max|dH| {worstH:.2e}")
    check("046 E: the n = 2 cell weights have both marginals (x, 1-x)",
          marg_bad == 0)
    # t = the off-diagonal cell weight, and w00 = x - t, w11 = (1-x) - t
    bad = 0
    for x in [0.5 + 0.5 * k / 50 for k in range(1, 50)]:
        z0 = clip(x, x)
        t = min(x - 0.5, 1 - x)
        if (abs((x - z0) - t) > 1e-14 or abs(z0 - (x - t)) > 1e-14
                or abs((1 - 2 * x + z0) - (1 - x - t)) > 1e-14):
            bad += 1
    check("046 E: t = min(x-1/2, 1-x) IS the off-diagonal cell weight, "
          "w00 = x - t, w11 = (1-x) - t", bad == 0)


# ============================================================== R2
def psi(s):
    return h(min(0.5, s))


def G(p):
    """sigma(p) = h(z*(p)) - h(p) = c*(p) h(p), the diagonal surplus."""
    return h(zstar(p)) - h(p)


def sigma_046(p):
    """046 E's sigma: G(p) for p <= 1/2, 0 above."""
    return G(p) if p <= 0.5 else 0.0


def s_cell(pa, pb):
    """The 037 cell ledger s(a,b) = h(clip(a,b)) - (h(a)+h(b))/2, in
    p-coordinates (p = 1 - u)."""
    return h(clip(1 - pa, 1 - pb)) - (h(pa) + h(pb)) / 2


def r2():
    print("\nR2. The hand lemmas of 046/047.")
    rng = random.Random(4902)
    both_bad = 0
    both_above = 0
    conc_bad = 0
    conc_badB = 0
    avg_bad = 0
    dec_worst = 0.0
    decA_worst = 0.0
    dec_worst_closed = 0.0
    caseB = 0
    caseB_C_neg = 0
    n_seen = 0
    while n_seen < 20000:
        x = rng.uniform(0.5, 1.0)
        u0, u1 = rng.uniform(0, 1), rng.uniform(0, 1)
        p0, p1 = 1 - u0, 1 - u1
        f0, f1 = 1 - x, x * p0 + (1 - x) * p1
        q = max(f0, f1)
        if not (0 < q < 0.5):
            continue
        n_seen += 1
        if u0 < 0.5 and u1 < 0.5:
            both_bad += 1
        if p0 > q + 1e-12 and p1 > q + 1e-12:
            both_above += 1
        # N2-CONC, cell form -- stated for the both-good case
        if p0 <= 0.5 and p1 <= 0.5:
            if s_cell(p0, p1) < (G(p0) + G(p1)) / 2 - 1e-12:
                conc_bad += 1
        elif s_cell(p0, p1) < (G(p0) + G(p1)) / 2 - 1e-12:
            conc_badB += 1
        # N2-CONC, averaged over the cells
        z0 = clip(x, x)
        w = {(0, 0): z0, (0, 1): x - z0, (1, 0): x - z0,
             (1, 1): 1 - 2 * x + z0}
        pp = (p0, p1)
        lhs = sum(wt * s_cell(pp[a], pp[b]) for (a, b), wt in w.items()
                  if wt > 0)
        if p0 <= 0.5 and p1 <= 0.5:
            if lhs < x * G(p0) + (1 - x) * G(p1) - 1e-12:
                avg_bad += 1
        # the exact decomposition of 046 E
        F = closed_margin(x, p0, p1)
        t = min(x - 0.5, 1 - x)
        Delta = 2 * s_cell(p0, p1) - s_cell(p0, p0) - s_cell(p1, p1)
        rhs = ((cstar(f0) - cstar(q)) * h(f0)
               + x * sigma_046(p0) + (1 - x) * sigma_046(p1)
               + t * Delta
               - cstar(q) * (x * h(p0) + (1 - x) * h(p1)))
        dec_worst = max(dec_worst, abs(F - rhs))
        if p0 <= 0.5 and p1 <= 0.5:
            decA_worst = max(decA_worst, abs(F - rhs))
        else:
            caseB += 1
            if t * Delta < -1e-12:
                caseB_C_neg += 1
            # 047's closed-form Delta in place of the ledger Delta
            Dc = 2 * psi(p0 + p1) - psi(2 * p0) - psi(2 * p1)
            rhs2 = ((cstar(f0) - cstar(q)) * h(f0)
                    + x * G(p0) + (1 - x) * G(p1) + t * Dc
                    - cstar(q) * (x * h(p0) + (1 - x) * h(p1)))
            dec_worst_closed = max(dec_worst_closed, abs(F - rhs2))
    check("046 N2-ONE-BAD: no in-regime sample has both conditionals "
          f"out of regime ({n_seen} samples)", both_bad == 0)
    check("047 N2-ONE-ABOVE: no in-regime sample has both conditional "
          "marginals above q", both_above == 0)
    check("046 N2-CONC (cell form): s(p_a,p_b) >= (G(p_a)+G(p_b))/2 in "
          "the both-good case", conc_bad == 0)
    check("SCOPE 046: the both-good hypothesis of N2-CONC is essential - "
          "the cell bound fails outside it", conc_badB > 0,
          f"{conc_badB} Case-B samples violate the cell bound")
    check("046 N2-CONC (averaging step): sum_ab w_ab s >= x G(p0) + "
          "(1-x) G(p1)", avg_bad == 0)
    check("046 E: the decomposition identity is EXACT over the whole "
          "regime (both cases), with sigma = 0 above 1/2 and Delta the "
          "ledger form", dec_worst < 1e-12,
          f"max |F - decomposition| {dec_worst:.2e} "
          f"(Case A alone {decA_worst:.2e}; {caseB} Case-B samples)")
    check("CORRECTION 047: outside Case A the closed-form Delta is NOT "
          "the ledger Delta, so 047's A+B+C is not the margin there",
          dec_worst_closed > 1e-6,
          f"max |F - (A+B+C_closed)| on Case B = {dec_worst_closed:.2e}")
    check("CORRECTION 047: C = t*Delta CAN be negative outside Case A, "
          "so \"only B can be negative\" holds in Case A only",
          caseB_C_neg > 0,
          f"{caseB_C_neg} of {caseB} Case-B samples have C < 0")
    # N2-CONC's tightness clause
    eq_unequal = []
    for pa, pb in ((0.30, 0.40), (0.26, 0.49), (0.25, 0.50)):
        if abs(s_cell(pa, pb) - (G(pa) + G(pb)) / 2) < 1e-15 and pa != pb:
            eq_unequal.append((pa, pb))
    check("CORRECTION 046: N2-CONC's \"equality iff p_a = p_b\" is FALSE "
          "- equality also holds whenever min(p_a,p_b) >= 1/4 (psi is "
          "constant there)", len(eq_unequal) == 3, f"{eq_unequal}")
    # and the corrected clause holds
    rng2 = random.Random(77)
    bad = 0
    for _ in range(200000):
        pa, pb = rng2.uniform(0, 0.5), rng2.uniform(0, 0.5)
        tight = abs(s_cell(pa, pb) - (G(pa) + G(pb)) / 2) < 1e-14
        expect = abs(pa - pb) < 1e-14 or min(pa, pb) >= 0.25
        if tight != expect:
            bad += 1
    check("corrected clause: equality iff p_a = p_b OR min(p_a,p_b) >= 1/4",
          bad == 0)


def closed_margin(x, p0, p1):
    """F = CR - c*(q) H at n = 2, from the closed form (046 A)."""
    u0, u1 = 1 - p0, 1 - p1
    f0, f1 = 1 - x, x * p0 + (1 - x) * p1
    q = max(f0, f1)
    z0 = clip(x, x)
    w = {(0, 0): z0, (0, 1): x - z0, (1, 0): x - z0,
         (1, 1): 1 - 2 * x + z0}
    uu = (u0, u1)
    S = h(z0) + sum(wt * h(clip(uu[a], uu[b])) for (a, b), wt in w.items()
                    if wt > 0)
    H = h(x) + x * h(u0) + (1 - x) * h(u1)
    return (S - H) - cstar(q) * H


# ============================================================== R3
def r3():
    print("\nR3. 046 F: the c* correction, both directions.")

    def stated(p):
        return (h(min(2 * p, 1.0)) - h(p)) / h(p)

    # (a) the stated form agrees below 1/4 and is wrong above
    agree = all(abs(stated(p) - cstar(p)) < 1e-14
                for p in [0.001 + 0.249 * k / 200 for k in range(201)])
    check("046 F: the stated form equals the correct one for p <= 1/4",
          agree)
    diff = all(stated(p) < cstar(p) - 1e-9
               for p in [0.25 + 0.25 * k / 200 for k in range(1, 200)])
    check("046 F: the stated form is strictly SMALLER than the correct "
          "one for every p in (1/4, 1/2)", diff)
    # the record's quoted numbers
    quoted = {0.38271: (-0.181287, 0.041739), 0.45: (-0.527591, 0.007278),
              0.49: (-0.858519, 0.000289), 0.497: (None, 0.0000260)}
    bad = []
    for p, (st, co) in quoted.items():
        if st is not None and abs(stated(p) - st) > 5e-6:
            bad.append(("stated", p, stated(p), st))
        if abs(cstar(p) - co) > 5e-6 * max(1.0, abs(co) * 100):
            bad.append(("correct", p, cstar(p), co))
    check("046 F: every quoted value of both forms reproduces", not bad,
          f"{bad}")
    # (b) BUT: negative only above 1/3, not above 1/4
    pos = [p for p in [0.25 + (1 / 3 - 0.25) * k / 100 for k in range(1, 100)]
           if stated(p) > 0]
    neg = [p for p in [1 / 3 + (0.5 - 1 / 3) * k / 100
                       for k in range(1, 100)] if stated(p) < 0]
    root = abs(stated(1 / 3)) < 1e-12
    check("CORRECTION 046 F: the stated form is NEGATIVE only above "
          "p = 1/3 (h(2p) = h(p) iff 3p = 1), not above p = 1/4; it is "
          "positive-but-wrong on (1/4, 1/3)",
          len(pos) == 99 and len(neg) == 99 and root,
          f"stated(0.30) = {stated(0.30):+.6f} > 0, "
          f"stated(1/3) = {stated(1/3):+.1e}, "
          f"stated(0.40) = {stated(0.40):+.6f}")
    # (c) every engine computes the max(1/2, 1-2p) form
    good = re.compile(r"max\(\s*(?:0\.5|HALF|Fraction\(1,\s*2\))\s*,\s*"
                      r"1(?:\.0)?\s*-\s*2(?:\.0)?\s*\*\s*p\s*\)")
    badform = re.compile(r"min\(\s*2\s*\*?\s*\*?\s*p")
    hits, offenders = 0, []
    for f in sorted(HERE.glob("*.py")):
        src = f.read_text()
        for m in re.finditer(r"def (cstar|cst|zstar|zst|c_star|sig|sigma)"
                             r"\b[^\n]*\n((?:\s+[^\n]*\n)+?)(?=\S)",
                             src):
            body = m.group(2)
            hits += 1
            if not good.search(body) and "z*(p)" not in body \
                    and "zst(" not in body and "zstar(" not in body:
                offenders.append((f.name, m.group(1)))
    check("046 F: every c*/z*/sigma definition in explore/ uses the "
          f"max(1/2, 1-2p) form ({hits} definitions scanned)",
          not offenders, f"offenders {offenders}")
    # and: where does the WRONG form appear in live code at all?
    root = HERE.parent.parent.parent
    live = []
    for f in sorted(root.rglob("*.py")):
        if f.name == "uc_reviewer049_reimpl.py":
            continue
        for ln, line in enumerate(f.read_text().splitlines(), 1):
            if re.search(r"min\(\s*2(?:\.0)?\s*\*\s*p(?:bar)?\s*,\s*1",
                         line):
                benign = ("say(" in line or "print(" in line
                          or line.lstrip().startswith("#")
                          or line.lstrip().startswith("f\""))
                live.append((str(f.relative_to(root)), ln, benign))
    check("046 F: no live computation anywhere in the repo evaluates the "
          "mis-stated h(min(2p,1)) form (the only occurrences are "
          "diagnostic prints)", all(b for _, _, b in live),
          f"{[(a, b) for a, b, _ in live]}")


# ============================================================== R4
def r4():
    print("\nR4. The 044 certificate, recomputed at 50+ digits.")
    o2 = json.loads((DATA / "hu_order2.json").read_text())
    key = "D_bestorder_cap_0.497"
    rows = [r for r in o2[key]["rows"] if r["start"] == "floor:windowkill"]
    check("044 A: the kill witness is ALREADY PRESENT in 037's committed "
          "checkpoint (hu_order2.json / D_bestorder_cap_0.497 / start "
          "floor:windowkill)", len(rows) == 1,
          f"{len(rows)} matching row(s), n = {rows[0]['n'] if rows else '-'}")
    row = rows[0]
    n = row["n"]
    muF = {int(s, 2): w for s, w in row["mu"].items()}
    muQ = {a: Fraction(w).limit_denominator(10 ** 7) for a, w in muF.items()}
    tot = sum(muQ.values())
    muQ = {a: w / tot for a, w in muQ.items()}
    check("044 B: the witness is a 9-atom n = 4 measure",
          n == 4 and len(muQ) == 9, f"n = {n}, atoms = {len(muQ)}")
    marg = max(sum(w for a, w in muQ.items() if (a >> i) & 1)
               for i in range(n))
    check("044 B: max marginal 0.49500 < 0.497 < 1/2",
          Fraction(494, 1000) < marg < Fraction(497, 1000),
          f"{float(marg):.8f}")
    # own rollout order
    seq = order_roll(n, muQ)
    check("044 B: own rollout-order derivation gives (1,0,2,3)",
          seq == (1, 0, 2, 3), f"{seq}")
    cert = json.loads((DATA / "hu_ownconst_certify.json").read_text())
    check("044 B: the certified rollout/best orders match this file's "
          "derivation / enumeration",
          tuple(cert["roll_seq"]) == seq, f"checkpoint {cert['roll_seq']}")
    cr_roll, H = hu_cr(n, muQ, seq, dec=True)
    print(f"      own CR_roll  = {cr_roll:.52}")
    print(f"      own H        = {H:.52}")
    best = max(itertools.permutations(range(n)),
               key=lambda s: float(hu_cr(n, muQ, s, dec=False)[0]))
    cr_best, _ = hu_cr(n, muQ, best, dec=True)
    print(f"      own best ord = {best}, CR_best = {cr_best:.52}")
    check("044 B: the best order is (0,3,2,1) and matches the checkpoint",
          best == (0, 3, 2, 1) == tuple(cert["best_seq"]), f"{best}")
    check("044 B: CR under rollout = -9.596936657e-4 (own 50-digit value)",
          abs(cr_roll - Decimal("-9.596936657e-4")) < Decimal("5e-14"),
          f"{cr_roll:.20}")
    check("044 B: CR under the best order = +1.610189393e-2",
          abs(cr_best - Decimal("1.610189393e-2")) < Decimal("5e-12"),
          f"{cr_best:.20}")
    check("044 B: H = 3.137 as quoted",
          abs(H - Decimal("3.1373367156258767")) < Decimal("1e-15"),
          f"{H:.20}")
    # containment in BOTH kits' enclosures
    encl = []
    for tag, val in (("K1_roll", cr_roll), ("K4_best", cr_best)):
        for kit in ("A_digit_extraction", "B_atanh_series"):
            k = cert["kits"][f"{tag}:{kit}"]
            lo, hi = Decimal(repr(k["CR_lo"])), Decimal(repr(k["CR_hi"]))
            slack = Decimal("1e-16")
            encl.append((f"{tag}:{kit}", lo - slack <= val <= hi + slack))
    check("044 B: the own 50-digit values lie inside BOTH kits' "
          "enclosures for both orders (float-rounded endpoints, 1e-16 "
          "slack)", all(ok for _, ok in encl), f"{encl}")
    check("044 B: the kill is certified negative and the control "
          "positive under each kit alone",
          all(v["certifies"] for v in cert["kits"].values()))
    # own-constant margins
    cs = cstar(float(marg))
    roll_ratio = float(cr_roll) / float(H)
    best_ratio = float(cr_best) / float(H)
    check("044 B: rollout violates its OWN constant; best order clears it",
          roll_ratio < 0 < cs < best_ratio,
          f"CR_roll/H {roll_ratio:+.3e} < 0, CR_best/H {best_ratio:+.5f}, "
          f"c*(fmax) {cs:.3e}")
    check("044 B: the record's best-order ratio +0.00513 and "
          "c*(0.495) = 7.2e-5", abs(best_ratio - 0.00513) < 5e-6
          and abs(cs - 7.21e-5) < 5e-7, f"{best_ratio:.5f}, {cs:.3e}")
    # 18 of 24 orders negative, rollout ranks 13/24
    vals = sorted(((float(hu_cr(n, muQ, s, dec=False)[0]), s)
                   for s in itertools.permutations(range(n))), reverse=True)
    negs = sum(1 for v, _ in vals if v < 0)
    rank = 1 + [s for _, s in vals].index(seq)
    rank = min(i for i, (v, _) in enumerate(vals, 1)
               if abs(v - float(cr_roll)) < 1e-15)
    tie = sum(1 for v, _ in vals if abs(v - float(cr_roll)) < 1e-15)
    check("044 B: 18 of the 24 orders are negative", negs == 18)
    check("CORRECTION 044: \"rollout ranks 13/24\" is not reproducible - "
          "rollout sits in a 6-way exact tie at ranks 7-12 by descending "
          "CR (13 is the tie block's position counting from the WORST "
          "end, under which rollout is ABOVE the median, not \"worse "
          "than the middle\")",
          rank == 7 and tie == 6 and 24 - 12 + 1 == 13,
          f"descending rank of the tie block 7-12 ({tie} orders tied at "
          f"{float(cr_roll):+.6e}); 6 orders strictly better, 12 strictly "
          "worse")
    return muQ, n


# ============================================================== R5
def r5(muQ, n4):
    print("\nR5. 045: every hu_bestorder.json endpoint re-scored.")
    bo = json.loads((DATA / "hu_bestorder.json").read_text())
    quoted = {"cap_0.495": 2.190730e-4, "cap_0.497": 2.190730e-4,
              "cap_0.499": 1.295407e-4}
    worst_dev = worst_own = 0.0
    kills = sharp = rows_seen = 0
    floors = {}
    stall_min = None
    stall_rows = 0
    fam0 = []
    kill_rows = {}
    for cap in sorted(k for k in bo if k.startswith("cap_")):
        fl = None
        for row in bo[cap]["rows"]:
            n = row["n"]
            mu = {int(s, 2): w for s, w in row["mu"].items() if w > 0}
            tot = sum(mu.values())
            mu = {a: w / tot for a, w in mu.items()}
            rows_seen += 1
            vals = [float(hu_cr(n, mu, s, dec=False)[0])
                    for s in itertools.permutations(range(n))]
            H = float(hu_cr(n, mu, tuple(range(n)), dec=False)[1])
            ratio = max(vals) / H
            fmax = max(marginals(n, mu))
            own = ratio - cstar(fmax)
            worst_dev = max(worst_dev, abs(ratio - row["floor"]))
            worst_own = max(worst_own, abs(own - row["own_margin"]))
            if ratio < 0:
                kills += 1
            if own < -1e-9:
                sharp += 1
            if 1e-12 < abs(own) < 1e-8:
                stall_rows += 1
                if stall_min is None or abs(own) < abs(stall_min):
                    stall_min = own
            if row["start"].startswith("family(('d', 4)"):
                fam0.append(own)
            if row["start"] == "044kill":
                kill_rows[cap] = ratio
            fl = ratio if fl is None else min(fl, ratio)
        floors[cap] = fl
    check(f"045: all {rows_seen} endpoint rows re-scored by full "
          "24/120-order enumeration through this file's evaluator",
          worst_dev < 1e-9 and worst_own < 1e-9,
          f"max |ratio - checkpoint floor| {worst_dev:.2e}, "
          f"max |own margin - checkpoint| {worst_own:.2e}")
    check("045: zero kills (ratio < 0) at every cap", kills == 0)
    check("045: zero sharp violations (ratio < c*(own fmax))", sharp == 0,
          "and the checkpoint's kill / sharp lists are empty: "
          + str(all(not bo[c]["kills"] and not bo[c]["sharp_violations"]
                    for c in bo)))
    check("045: the three quoted global floors reproduce",
          all(abs(floors[c] - quoted[c]) < 5e-10 for c in quoted),
          ", ".join(f"{c} {floors[c]:+.6e}" for c in sorted(floors)))
    check("045 obs 1: several endpoints saturate their own constant to "
          "+2.7e-9", stall_rows >= 5 and abs(stall_min) < 3e-9,
          f"{stall_rows} rows with 1e-12 < |own margin| < 1e-8, smallest "
          f"{stall_min:+.2e}")
    check("045 obs 2: the (d,4)-family perturbation returns to own "
          "margin EXACTLY 0 at every cap",
          bool(fam0) and all(abs(v) < 1e-15 for v in fam0),
          f"{fam0}")
    check("045 obs 3: the 044 kill witness is comfortably positive under "
          "best order; the quoted +3.07e-3 is the cap-0.499 row (it is "
          "+5.13e-3 at caps 0.495 and 0.497)",
          abs(kill_rows["cap_0.499"] - 3.071221e-3) < 1e-9
          and abs(kill_rows["cap_0.495"] - 5.132345e-3) < 1e-9,
          ", ".join(f"{c} {v:+.6e}" for c, v in sorted(kill_rows.items())))


# ============================================================== R6
def r6():
    print("\nR6. 046's numbers on an own grid/sampler.")
    d = json.loads((DATA / "hu_n2.json").read_text())
    # order quantifier vacuous: own grid over the whole box
    M = 90
    worst = None
    grid = [k / M for k in range(M + 1)]
    for xi in grid:
        x = 0.5 + 0.5 * xi
        for u0 in grid:
            for u1 in grid:
                f0, f1 = 1 - x, x * (1 - u0) + (1 - x) * (1 - u1)
                q = max(f0, f1)
                if not (0 < q < 0.5):
                    continue
                m = closed_margin(x, 1 - u0, 1 - u1)
                if worst is None or m < worst[0]:
                    worst = (m, x, u0, u1)
    check("046 A: identity-order margin >= 0 on an own "
          f"{M + 1}^3 grid of the regime box", worst[0] > -1e-12,
          f"worst {worst[0]:+.2e} at x={worst[1]:.4f}, u0={worst[2]:.4f}, "
          f"u1={worst[3]:.4f}")
    # order quantifier: identity vs best (relabelling) on samples
    rng = random.Random(4906)
    below = 0
    seen = 0
    worst_id = None
    while seen < 3000:
        x = rng.uniform(0.5, 1.0)
        u0, u1 = rng.uniform(0, 1), rng.uniform(0, 1)
        f0, f1 = 1 - x, x * (1 - u0) + (1 - x) * (1 - u1)
        q = max(f0, f1)
        if not (0 < q < 0.5):
            continue
        seen += 1
        xq = Fraction(x).limit_denominator(10 ** 6)
        u0q = Fraction(u0).limit_denominator(10 ** 6)
        u1q = Fraction(u1).limit_denominator(10 ** 6)
        mu = {0b00: xq * u0q, 0b10: xq * (1 - u0q),
              0b01: (1 - xq) * u1q, 0b11: (1 - xq) * (1 - u1q)}
        mu = {a: w for a, w in mu.items() if w > 0}
        c01, H = hu_cr(2, mu, (0, 1), dec=False)
        c10, _ = hu_cr(2, mu, (1, 0), dec=False)
        if float(c01) < float(c10) - 1e-12:
            below += 1
        m = float(c01) - cstar(q) * float(H)
        if worst_id is None or m < worst_id:
            worst_id = m
    check("046 A: the identity order is often strictly worse than the "
          "best order yet its margin never goes below 0",
          below > seen // 10 and worst_id > -1e-12,
          f"{below}/{seen} strictly worse, worst identity margin "
          f"{worst_id:+.2e}")
    # B&B volume bookkeeping
    c = d["C_bnb"]
    inreg = c["vol_root"] - c["vol_out_of_regime"]
    frac = c["vol_certified"] / inreg
    unproc = c["vol_unprocessed"] / inreg
    check("046 C: the committed B&B is the largest-box-first version and "
          "its volumes add up (certified + out-of-regime + residue + "
          "unprocessed = root)",
          abs(c["vol_certified"] + c["vol_out_of_regime"]
              + c["vol_residue"] + c["vol_unprocessed"] - c["vol_root"])
          < 1e-9 and c["vol_residue"] == 0.0,
          f"certified {100 * frac:.1f}% of in-regime, unprocessed "
          f"{100 * unproc:.1f}%, residue {c['vol_residue']}")
    check("046 C: the record's quoted 88.4% / 11.6% / zero residue match "
          "the checkpoint",
          abs(frac - 0.884) < 5e-4 and abs(unproc - 0.116) < 5e-4
          and c["residue_boxes"] == 0)
    check("046 C: min edge 1/512 and budget 2,000,000 as quoted",
          abs(c["min_edge"] - 1 / 512) < 1e-15 and c["processed"] == 2000000)
    # part E: (**) on an own sampler
    rng = random.Random(4907)
    tested = bad = 0
    worst_pt = None
    while tested < 200000:
        x = rng.uniform(0.5, 1.0)
        p0, p1 = rng.uniform(0, 0.5), rng.uniform(0, 0.5)
        f0, f1 = 1 - x, x * p0 + (1 - x) * p1
        q = max(f0, f1)
        if not (0 < q < 0.5) or min(p0, p1) <= 0:
            continue
        tested += 1
        lhs = x * G(p0) + (1 - x) * G(p1)
        rhs = cstar(q) * (x * h(p0) + (1 - x) * h(p1))
        if lhs < rhs - 1e-12:
            bad += 1
            if worst_pt is None or lhs - rhs < worst_pt[0]:
                worst_pt = (lhs - rhs, x, p0, p1, q)
    frac = bad / tested
    check("046 E: (**) fails on ~21.3% of Case-A points (own sampler)",
          abs(frac - 0.21345) < 0.01, f"{100 * frac:.1f}% of {tested}")
    m = closed_margin(worst_pt[1], worst_pt[2], worst_pt[3])
    check("046 E: at the worst (**) deficit the ORIGINAL margin is "
          "positive", m > 0,
          f"deficit {worst_pt[0]:+.3e} at x={worst_pt[1]:.4f}, "
          f"p0={worst_pt[2]:.4f}, p1={worst_pt[3]:.4f}; margin {m:+.3e}")
    check("046 E: the record's quoted worst deficit -9.72e-2 is of the "
          "same order as an own sampler's",
          abs(worst_pt[0]) > 0.05, f"own worst {worst_pt[0]:+.3e}")


# ============================================================== R7
def r7():
    print("\nR7. 047's dichotomy census, own implementation.")
    d = json.loads((DATA / "hu_n2_dichotomy.json").read_text())

    def terms(x, p0, p1):
        f0, f1 = 1 - x, x * p0 + (1 - x) * p1
        q = max(f0, f1)
        if not (0 < q < 0.5):
            return None
        A = (cstar(f0) - cstar(q)) * h(f0)
        B = (x * (G(p0) - cstar(q) * h(p0))
             + (1 - x) * (G(p1) - cstar(q) * h(p1)))
        t = min(x - 0.5, 1 - x)
        C = t * (2 * psi(p0 + p1) - psi(2 * p0) - psi(2 * p1))
        return A, B, C, q, t

    rng = random.Random(4908)
    seen = deficit = byA = byC = both_needed = fails = 0
    split = {"q=f0": [0, 0, 0], "q=f1": [0, 0, 0]}
    worst_dec = 0.0
    minL1 = None
    minmax = None
    while seen < 200000:
        x, p0, p1 = rng.uniform(0.5, 1), rng.uniform(0, 0.5), \
            rng.uniform(0, 0.5)
        r = terms(x, p0, p1)
        if r is None:
            continue
        A, B, C, q, t = r
        seen += 1
        # the decomposition is the margin (Case A)
        worst_dec = max(worst_dec, abs(A + B + C - closed_margin(x, p0, p1)))
        if B >= 0:
            continue
        deficit += 1
        okA, okC = A + B >= 0, B + C >= 0
        byA += okA
        byC += okC
        if not okA and not okC:
            if A + B + C >= 0:
                both_needed += 1
            else:
                fails += 1
        k = "q=f0" if 1 - x >= q - 1e-15 else "q=f1"
        split[k][0] += 1
        split[k][1] += (not okC)
        split[k][2] += (not okA)
        if k == "q=f0" and (minL1 is None or B + C < minL1):
            minL1 = B + C
        v = max(A, C) + B
        if minmax is None or v < minmax:
            minmax = v
    check("047: the A+B+C decomposition equals the margin (own closed "
          "form) on 200,000 Case-A samples", worst_dec < 1e-12,
          f"max |A+B+C - F| {worst_dec:.2e}")
    check("047 B: needing BOTH terms = 0 and failing outright = 0 "
          "reproduce on an own sampler",
          both_needed == 0 and fails == 0,
          f"{deficit} deficit cases, covered by A {byA}, by C {byC}")
    check("047 B: the deficit rate reproduces (checkpoint 42,575 of "
          "200,000)", abs(deficit / seen - 42575 / 200000) < 0.01,
          f"own {deficit}/{seen}")
    check("047 C: on the q = f0 branch the interaction alone NEVER fails",
          split["q=f0"][1] == 0,
          f"q=f0 {split['q=f0'][0]} cases, C fails {split['q=f0'][1]}, "
          f"A fails {split['q=f0'][2]}")
    check("047 C: on the q = f1 branch both single-term covers fail "
          "sometimes, but never together",
          split["q=f1"][1] > 0 and split["q=f1"][2] > 0,
          f"q=f1 {split['q=f1'][0]} cases, C fails {split['q=f1'][1]}, "
          f"A fails {split['q=f1'][2]}")
    check("047 D: the two tightness minima are positive and of the "
          "quoted order (+1.17e-4 / +1.66e-5)",
          minL1 > 0 and minmax > 0,
          f"own min B+C on L1 {minL1:+.2e}, own min max(A,C)+B "
          f"{minmax:+.2e}")
    # the record's own checkpoint: part B and part C use different seeds
    b, cs = d["B_dichotomy"], d["C_split"]["stats"]
    check("REPORTING 047: parts B and C are different samples - their "
          "deficit totals do not agree (42,575 vs 7,435 + 34,650 = "
          "42,085), which the record's prose does not say",
          b["deficit_cases"] != cs["q=f0"]["cases"] + cs["q=f1"]["cases"],
          f"{b['deficit_cases']} vs "
          f"{cs['q=f0']['cases'] + cs['q=f1']['cases']}")
    # SCOPE: Case B is a real part of the regime and the census skips it
    rng = random.Random(4909)
    tot = caseB = 0
    while tot < 200000:
        x = rng.uniform(0.5, 1.0)
        u0, u1 = rng.uniform(0, 1), rng.uniform(0, 1)
        f0, f1 = 1 - x, x * (1 - u0) + (1 - x) * (1 - u1)
        if not (0 < max(f0, f1) < 0.5):
            continue
        tot += 1
        if u0 < 0.5 or u1 < 0.5:
            caseB += 1
    check("SCOPE 047: Case B (one conditional out of regime) is a "
          "substantial part of the n = 2 regime that the dichotomy "
          "census never samples", caseB > 0.05 * tot,
          f"{100 * caseB / tot:.1f}% of in-regime measures are Case B")
    # but the margin itself is still nonnegative there (own check)
    rng = random.Random(4910)
    worstB = None
    seen = 0
    while seen < 200000:
        x = rng.uniform(0.5, 1.0)
        u0, u1 = rng.uniform(0, 1), rng.uniform(0, 1)
        f0, f1 = 1 - x, x * (1 - u0) + (1 - x) * (1 - u1)
        if not (0 < max(f0, f1) < 0.5) or (u0 >= 0.5 and u1 >= 0.5):
            continue
        seen += 1
        m = closed_margin(x, 1 - u0, 1 - u1)
        if worstB is None or m < worstB:
            worstB = m
    check("047: (HU-TAX) at n = 2 is not endangered in Case B - the "
          f"margin stays >= 0 over {seen} Case-B samples", worstB > -1e-12,
          f"worst Case-B margin {worstB:+.2e}")


# ============================================================== R8
def r8():
    print("\nR8. 048: branch L1 in ratio form.")
    d = json.loads((DATA / "hu_L1.json").read_text())

    def L1(q, p0, p1):
        x = 1 - q
        if not (0 < q < 0.5) or x * p0 + (1 - x) * p1 > q + 1e-15:
            return None
        negB = (x * (h(p0) * cstar(q) - G(p0))
                + (1 - x) * (h(p1) * cstar(q) - G(p1)))
        t = min(x - 0.5, 1 - x)
        return negB, t * (2 * psi(p0 + p1) - psi(2 * p0) - psi(2 * p1)), t

    # P2: the boundary identity, in EXACT rational/Decimal arithmetic.
    # Own derivation at (p0, p1) = (0, 1/2): Delta = 2h(1/2) - h(0) - h(1/2)
    # = 1, so C = t; and -B = x h(0)[c*(q)-c*(0)] + (1-x) h(1/2)[c*(q) -
    # c*(1/2)] = q c*(q) since h(0) = 0, h(1/2) = 1, c*(1/2) = 0.  Hence
    # ratio = t/(q c*(q)), which is 1/c*(q) ONLY when t = q, i.e. q <= 1/4.
    small, large = [], []
    for qs in ("1/10", "1/100", "1/1000", "1/10000", "1/100000",
               "1/1000000", "1/4", "1/5"):
        q = Fraction(qs)
        t = min(1 - q - Fraction(1, 2), q)
        cq = (h_dec(max(Fraction(1, 2), 1 - 2 * q)) - h_dec(q)) / h_dec(q)
        ratio = (Decimal(t.numerator) / Decimal(t.denominator)) / (
            (Decimal(q.numerator) / Decimal(q.denominator)) * cq)
        small.append((qs, t == q, abs(ratio - 1 / cq) < Decimal("1e-60")))
    for qs in ("7/23", "3/10", "2/5", "9/20"):
        q = Fraction(qs)
        t = min(1 - q - Fraction(1, 2), q)
        cq = (h_dec(max(Fraction(1, 2), 1 - 2 * q)) - h_dec(q)) / h_dec(q)
        ratio = (Decimal(t.numerator) / Decimal(t.denominator)) / (
            (Decimal(q.numerator) / Decimal(q.denominator)) * cq)
        pred = (Decimal((Fraction(1, 2) - q).numerator)
                / Decimal((Fraction(1, 2) - q).denominator)) / (
            (Decimal(q.numerator) / Decimal(q.denominator)) * cq)
        large.append((qs, float(ratio), float(1 / cq),
                      abs(ratio - pred) < Decimal("1e-60"), float(ratio) > 1))
    check("048 P2: ratio(q, p0=0, p1=1/2) = 1/c*(q) EXACTLY for q <= 1/4 "
          "(own derivation: C = t = q, -B = q c*(q)), to 60 digits",
          all(a and b for _, a, b in small), f"{small}")
    check("CORRECTION 048 P2: the identity is stated unqualified but "
          "holds only for q <= 1/4 - above that t = 1/2 - q < q and the "
          "ratio at (0,1/2) is (1/2-q)/(q c*(q)), not 1/c*(q); the "
          "record's table only goes up to q = 0.1",
          all(exact and abs(r - inv) > 1e-6 for _, r, inv, exact, _ in large),
          "; ".join(f"q={a}: ratio {b:.5f} vs 1/c* {c:.5f}"
                    for a, b, c, _, _ in large))
    check("048 P2: the corner configuration still satisfies L1 strictly "
          "at those q (the record's conclusion survives its proof)",
          all(g for *_, g in large))
    check("048 P2: the six quoted table rows reproduce",
          all(abs(r["ratio"] - r["inv_cstar"]) < 1e-12
              and abs(r["ratio"] - 1 / cstar(r["q"])) < 1e-9
              for r in d["P2_extremal"]))
    # c*(q) < 1 on (0,1/2) and -> 1 only logarithmically
    ok = all(cstar(q) < 1 for q in [1e-8, 1e-4, 0.01, 0.1, 0.25, 0.3,
                                    0.45, 0.499])
    def cstar_exact(q):
        return (h_dec(max(Fraction(1, 2), 1 - 2 * q)) - h_dec(q)) / h_dec(q)
    seq = [cstar_exact(Fraction(1, 10 ** k)) for k in (3, 6, 12, 24, 100)]
    slow = all(seq[i] < seq[i + 1] < 1 for i in range(len(seq) - 1))
    lograte = abs(float(1 - seq[-1])
                  - 1 / (100 * math.log2(10))) < 0.02
    check("048 P2: c*(q) < 1 on (0,1/2), and c*(q) -> 1 only "
          "logarithmically (1 - c*(q) ~ 1/log2(1/q))",
          ok and slow and lograte,
          "c* at q = 1e-3, 1e-6, 1e-12, 1e-24, 1e-100: "
          + ", ".join(f"{v:.5f}" for v in seq))
    # P1: own descent-free sampling floor
    rng = random.Random(4911)
    lo = None
    seen = deficit = 0
    while seen < 400000:
        q = rng.uniform(1e-5, 0.5)
        p1 = rng.uniform(0, 0.5)
        hi = min(0.5, q * (1 - p1) / (1 - q))
        p0 = rng.uniform(0, hi)
        r = L1(q, p0, p1)
        if r is None:
            continue
        seen += 1
        negB, C, t = r
        if negB <= 1e-12:
            continue
        deficit += 1
        ratio = C / negB
        if lo is None or ratio < lo[0]:
            lo = (ratio, q, p0, p1)
    check("048 P1: the ratio C/(-B) stays above 1 on an own 400k sample "
          f"of the branch ({deficit} with a deficit)", lo[0] > 1.0,
          f"own min ratio {lo[0]:.4f} at q={lo[1]:.2e}, p0={lo[2]:.2e}, "
          f"p1={lo[3]:.4f} (checkpoint 1.1492)")
    # P3: which arm binds -- and the count discrepancy
    cnt = {"q": 0, "x-1/2": 0}
    tight = {"q": 0, "x-1/2": 0}
    rng = random.Random(4903)
    for _ in range(600000):
        q = rng.uniform(1e-5, 0.5)
        p1 = rng.uniform(0, 0.5)
        hi = min(0.5, q * (1 - p1) / (1 - q))
        if hi <= 0:
            continue
        p0 = rng.uniform(0, hi)
        r = L1(q, p0, p1)
        if r is None:
            continue
        negB, C, t = r
        x = 1 - q
        k = "q" if (1 - x) <= (x - 0.5) else "x-1/2"
        cnt[k] += 1
        if negB > 1e-12 and C / negB < 1.5:
            tight[k] += 1
    ck = d["P3_arm"]["tight_counts"]
    check("048 P3: every tight (ratio < 1.5) configuration has "
          "t = 1 - x = q; none has t = x - 1/2", tight["x-1/2"] == 0,
          f"own {tight}")
    check("CORRECTION 048 P3: the record quotes \"375 versus 0\" but the "
          "committed checkpoint says 5,698 versus 0 - 375 is the "
          "--fast (40,000-draw) run",
          ck["t=1-x (=q)"] == 5698 and abs(tight["q"] - 5698) <= 5,
          f"checkpoint {ck['t=1-x (=q)']}, own full-run {tight['q']}, "
          "record 375")
    # P4
    rng = random.Random(4912)
    fails = seen = 0
    crude_viol = 0
    while seen < 200000:
        q = rng.uniform(1e-5, 0.5)
        p1 = rng.uniform(0, 0.5)
        hi = min(0.5, q * (1 - p1) / (1 - q))
        p0 = rng.uniform(0, hi)
        r = L1(q, p0, p1)
        if r is None:
            continue
        seen += 1
        negB, C, t = r
        crude = (1 - q) * h(p0) + q * h(p1)
        if negB > crude + 1e-12:
            crude_viol += 1
        if C < crude - 1e-12:
            fails += 1
    check("048 P4: -B <= x h(p0) + (1-x) h(p1) always, and t*Delta "
          "dominates that crude bound on 0% of the branch",
          crude_viol == 0 and fails == seen,
          f"crude-bound violations {crude_viol}, C >= crude on "
          f"{seen - fails}/{seen}")
    # P5: the minimiser table, own grid
    rows = {r["q"]: r for r in d["P5_minimiser"]}
    devs = []
    for q, r in rows.items():
        M = 140
        best = None
        for i in range(M + 1):
            p1 = 0.5 * i / M
            hi = min(0.5, q * (1 - p1) / (1 - q))
            for j in range(M + 1):
                p0 = hi * j / M
                res = L1(q, p0, p1)
                if res is None:
                    continue
                negB, C, t = res
                if negB <= 1e-12:
                    continue
                if best is None or C / negB < best[0]:
                    best = (C / negB, p0, p1)
        at_corner = abs(best[1]) < 1e-9 and abs(best[2] - 0.5) < 1e-9
        devs.append((q, best[0], r["min_ratio"], at_corner, r["at_corner"]))
    ok = all(abs(a - b) < 2e-2 and c == e for _, a, b, c, e in devs)
    check("048 P5: the two-regime minimiser table reproduces on an own "
          "140^2 mesh (corner (0,1/2) for q <= 0.03, interior above)", ok,
          "; ".join(f"q={q}: own {a:.4f} vs {b:.4f} corner={c}"
                    for q, a, b, c, _ in devs))
    check("REPORTING 048 P5: the q = 0.3 minimiser (0.42857, 0.0) is "
          "called \"interior\" though p1 = 0 is on the boundary of the "
          "branch box", rows[0.3]["p1"] == 0.0,
          "cosmetic: 'at_corner' means the (0,1/2) corner only")


def main():
    r1()
    r2()
    r3()
    muQ, n4 = r4()
    r5(muQ, n4)
    r6()
    r7()
    r8()
    print()
    if FAILS:
        print(f"FAILED CHECKS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("All checks pass.")


if __name__ == "__main__":
    main()
