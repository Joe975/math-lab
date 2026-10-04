#!/usr/bin/env python3
"""Skeptic pass on uc_hu_n2_dichotomy.py (attempt 047).  Stance: refute.

Independence: the A/B/C terms are rebuilt here from the RECORD's prose
in nats (own h, own c*, own psi), and -- the check that matters -- the
identity "A + B + C = CR - c*(q) H" is re-verified against a
from-scratch history-recursion HU evaluator that shares no code with
uc_hu_order2 or uc_hu_n2.  If the decomposition were wrong, the whole
dichotomy would be about the wrong quantity.

  S1  A + B + C really is the margin (own evaluator, Case A)
  S2  A >= 0 and C >= 0 attacked directly (the two sign claims the
      dichotomy rests on)
  S3  Lemma N2-ONE-ABOVE attacked with targeted draws
  S4  the dichotomy itself: hunt for a deficit case needing BOTH terms
      or failing outright, including a descent that maximises the
      shortfall
  S5  branch L1 attacked on its own: hunt for q = f0 with B + C < 0

Usage: python uc_hu_n2_dichotomy_skeptic.py   (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
LN2 = math.log(2.0)
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hb(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1 - p) * math.log(1 - p)) / LN2


def zst(p):
    return max(0.5, 1.0 - 2.0 * p)


def cst(p):
    return (hb(zst(p)) - hb(p)) / hb(p)


def sig(p):
    return hb(zst(p)) - hb(p)


def psi(s):
    return hb(min(0.5, s))


def terms(x, p0, p1):
    f0, f1 = 1 - x, x * p0 + (1 - x) * p1
    q = max(f0, f1)
    if not (0 < q < 0.5) or f0 <= 0 or f0 >= 1:
        return None
    A = (cst(f0) - cst(q)) * hb(f0)
    B = (x * (sig(p0) - cst(q) * hb(p0))
         + (1 - x) * (sig(p1) - cst(q) * hb(p1)))
    C = min(x - 0.5, 1 - x) * (2 * psi(p0 + p1) - psi(2 * p0)
                               - psi(2 * p1))
    return A, B, C, q


def hu_margin(x, p0, p1):
    """CR - c*(q) H by an independent history recursion (nats-based h)."""
    u0, u1 = 1 - p0, 1 - p1
    mu = {0b00: x * u0, 0b10: x * (1 - u0),
          0b01: (1 - x) * u1, 0b11: (1 - x) * (1 - u1)}
    mu = {a: w for a, w in mu.items() if w > 1e-15}
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    f0 = sum(w for a, w in mu.items() if a & 1)
    f1 = sum(w for a, w in mu.items() if a & 2)
    q = max(f0, f1)
    if not (0 < q < 0.5):
        return None
    H = -sum(w * math.log(w) for w in mu.values()) / LN2
    cells = [((), (), 1.0)]
    S = 0.0
    for i in range(2):
        nxt = []
        for pa, pb, w in cells:
            def cond(pref):
                sel = [(a, wt) for a, wt in mu.items()
                       if all(((a >> j) & 1) == v
                              for j, v in enumerate(pref))]
                m = sum(wt for _, wt in sel)
                return m, (sum(wt for a, wt in sel
                               if not (a >> i) & 1) / m if m else 0.0)
            ma, xx = cond(pa)
            mb, yy = cond(pb)
            z = min(max(0.5, xx + yy - 1.0), xx, yy)
            S += w * hb(z)
            for va, vb, cw in ((0, 0, z), (0, 1, xx - z),
                               (1, 0, yy - z), (1, 1, 1 - xx - yy + z)):
                if cw > 1e-15:
                    nxt.append((pa + (va,), pb + (vb,), w * cw))
        cells = nxt
    return S - H - cst(q) * H


def main():
    out = json.loads((DATA / "hu_n2_dichotomy.json").read_text())
    rng = random.Random(8123)

    print("S1. A + B + C is the margin (independent evaluator):")
    worst = 0.0
    n = 0
    for _ in range(30000):
        x, p0, p1 = (rng.uniform(0.5, 1), rng.uniform(0, 0.5),
                     rng.uniform(0, 0.5))
        t = terms(x, p0, p1)
        m = hu_margin(x, p0, p1)
        if t is None or m is None:
            continue
        n += 1
        worst = max(worst, abs(t[0] + t[1] + t[2] - m))
    check("decomposition equals the margin", worst < 1e-9,
          f"max |A+B+C - margin| = {worst:.2e} over {n} samples")

    print("S2. The two sign claims:")
    negA = negC = 0
    for _ in range(200000):
        x, p0, p1 = (rng.uniform(0.5, 1), rng.uniform(0, 0.5),
                     rng.uniform(0, 0.5))
        t = terms(x, p0, p1)
        if t is None:
            continue
        if t[0] < -1e-12:
            negA += 1
        if t[2] < -1e-12:
            negC += 1
    check("A >= 0 (c* decreasing, q >= f0)", negA == 0, f"{negA} negatives")
    check("C >= 0 (N2-CONC: psi concave)", negC == 0, f"{negC} negatives")

    print("S3. Lemma N2-ONE-ABOVE, targeted:")
    both = None
    for _ in range(400000):
        x = rng.uniform(0.5, 1)
        p0 = rng.uniform(0.2, 0.5)
        p1 = rng.uniform(0.2, 0.5)
        t = terms(x, p0, p1)
        if t is None:
            continue
        q = t[3]
        if p0 > q + 1e-12 and p1 > q + 1e-12:
            both = (x, p0, p1, q)
            break
    check("no in-regime instance has both conditionals above q",
          both is None, str(both) if both else "as proved")

    print("S4. The dichotomy, attacked (random + descent on shortfall):")
    worstgap = None
    for _ in range(300000):
        x, p0, p1 = (rng.uniform(0.5, 1), rng.uniform(0, 0.5),
                     rng.uniform(0, 0.5))
        t = terms(x, p0, p1)
        if t is None:
            continue
        A, B, C, q = t
        if B >= 0:
            continue
        gap = max(A, C) + B
        if worstgap is None or gap < worstgap[0]:
            worstgap = (gap, x, p0, p1)
    # descent to maximise the shortfall
    g, x, p0, p1 = worstgap
    step = 0.05
    while step > 1e-6:
        improved = False
        for d in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0),
                  (0, 0, 1), (0, 0, -1)):
            c = (min(1, max(0.5, x + step * d[0])),
                 min(0.5, max(0, p0 + step * d[1])),
                 min(0.5, max(0, p1 + step * d[2])))
            t = terms(*c)
            if t is None or t[1] >= 0:
                continue
            v = max(t[0], t[2]) + t[1]
            if v < g - 1e-15:
                x, p0, p1, g = c[0], c[1], c[2], v
                improved = True
        if not improved:
            step /= 2
    check("max(A,C) >= -B survives a descent that maximises the "
          "shortfall", g >= -1e-12,
          f"min max(A,C)+B = {g:+.3e} at x={x:.5f}, p0={p0:.5f}, "
          f"p1={p1:.5f}")

    print("S5. Branch L1 (q = f0) attacked on its own:")
    worstL1 = None
    for _ in range(300000):
        x = rng.uniform(0.5, 1)
        p0 = rng.uniform(0, 0.5)
        p1 = rng.uniform(0, 0.5)
        t = terms(x, p0, p1)
        if t is None:
            continue
        A, B, C, q = t
        if B >= 0 or abs(q - (1 - x)) > 1e-15:
            continue
        v = B + C
        if worstL1 is None or v < worstL1[0]:
            worstL1 = (v, x, p0, p1)
    check("q = f0 branch: the interaction alone always covers",
          worstL1 is None or worstL1[0] >= -1e-12,
          f"min B+C = {worstL1[0]:+.3e} at x={worstL1[1]:.5f}, "
          f"p0={worstL1[2]:.5f}, p1={worstL1[3]:.5f}" if worstL1 else "")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
