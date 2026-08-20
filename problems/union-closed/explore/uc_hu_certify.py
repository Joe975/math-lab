#!/usr/bin/env python3
"""Attempt 033 part B (queue 1(b3)): exact-rational certification of
CR_HU at the two sharpest G-gen floor instances.

Per 029's corrected standard: each certified-log2 kit must certify ALONE
(no intersection as a soundness step; agreement is reported as a
consistency check only).  Kits, run unmodified: kit A =
uc_gap1_skeptic_gemini.certified_log2 (022/024-audited), kit B =
uc_or_agg_probe.log2_enclosure (016/017 standard).

The instances are the cr_attack.json committed endpoints with weights
rationalized by Fraction.limit_denominator(10**7) and renormalized
exactly; every certified statement is about the rationalized instance
(within 1e-7 total variation of the float one).  The HU coupling's cell
tree is finite and every cell quantity (x, y, z, mass) is an exact
Fraction, so CR_HU is a finite sum of c*log2(r) with rational c, r.

Certified statements per instance and per kit, separately:
  C1  CR_HU > 0
  C2  CR_HU - (1/25) H(mu) > 0   (the ratio floor 0.04)
  C3  H(mu) > 1/2, max marginal < 38271/100000 (exact rational checks)
  C4  (the 033 campaign witness, order (0,1,3,2), cap 9/20):
      CR_HU < 0 -- the certified order-dependent violation

Standard library only; deterministic; ~seconds.
Usage: python uc_hu_certify.py
Checkpoint: ../data/hu_certify.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_gap1_skeptic_gemini import certified_log2 as log2_A
from uc_or_agg_probe import log2_enclosure as log2_B

DATA = HERE.parent / "data"
MARG = Fraction(38271, 100000)
HALF = Fraction(1, 2)
OUT = {"instances": []}


def ivl_add(a, b):
    return a[0] + b[0], a[1] + b[1]


def ivl_scale(c, iv):
    return (c * iv[0], c * iv[1]) if c >= 0 else (c * iv[1], c * iv[0])


def xlog2x_ivl(x, log2fn):
    if x == 0:
        return Fraction(0), Fraction(0)
    return ivl_scale(x, log2fn(x))


def h_ivl(z, log2fn):
    """Enclosure of h(z) = -z log2 z - (1-z) log2(1-z)."""
    a = ivl_scale(Fraction(-1), xlog2x_ivl(z, log2fn))
    b = ivl_scale(Fraction(-1), xlog2x_ivl(1 - z, log2fn))
    return ivl_add(a, b)


def hu_cells(n, mu):
    """Exact cell tree of the HU coupling: yields (mass, z) Fractions
    per coordinate/cell, plus returns the exact conditional structure."""
    pref = []
    for i in range(n):
        mask = (1 << i) - 1
        d = {}
        for A, w in mu.items():
            k = A & mask
            e = d.setdefault(k, [Fraction(0), Fraction(0)])
            e[0] += w
            if not (A >> i) & 1:
                e[1] += w
        pref.append(d)
    cells = {(0, 0): Fraction(1)}
    zs = []
    for i in range(n):
        nxt = {}
        for (a, b), w in cells.items():
            da, db = pref[i][a], pref[i][b]
            x = da[1] / da[0]
            y = db[1] / db[0]
            z = min(max(HALF, x + y - 1), x, y)
            zs.append((w, z))
            for (ai, bi), cw in {(0, 0): z, (0, 1): x - z,
                                 (1, 0): y - z,
                                 (1, 1): 1 - x - y + z}.items():
                if cw <= 0:
                    continue
                kk = (a | (0 if ai == 0 else 1 << i),
                      b | (0 if bi == 0 else 1 << i))
                nxt[kk] = nxt.get(kk, Fraction(0)) + w * cw
        cells = nxt
    return zs


def certify_instance(name, n, mu_float):
    mu = {a: Fraction(w).limit_denominator(10 ** 7)
          for a, w in mu_float.items()}
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    marg = max(sum(w for a, w in mu.items() if (a >> i) & 1)
               for i in range(n))
    assert marg < MARG, f"{name}: not in-regime after rationalization"
    zs = hu_cells(n, mu)
    row = {"instance": name, "n": n, "atoms": len(mu),
           "cells": len(zs), "kits": {}}
    for kitname, log2fn in (("A_digit_extraction", log2_A),
                            ("B_atanh_series", log2_B)):
        # H(mu) enclosure
        Hmu = (Fraction(0), Fraction(0))
        for w in mu.values():
            Hmu = ivl_add(Hmu, ivl_scale(Fraction(-1),
                                         xlog2x_ivl(w, log2fn)))
        # sum of cell-mass-weighted h(z)
        S = (Fraction(0), Fraction(0))
        for w, z in zs:
            if 0 < z < 1:
                S = ivl_add(S, ivl_scale(w, h_ivl(z, log2fn)))
        CR = (S[0] - Hmu[1], S[1] - Hmu[0])
        ratio_stat = (CR[0] - Hmu[1] / 25, CR[1] - Hmu[0] / 25)
        ok = CR[0] > 0 and ratio_stat[0] > 0 and Hmu[0] > HALF
        row["kits"][kitname] = {
            "CR_lo": float(CR[0]), "CR_hi": float(CR[1]),
            "H_lo": float(Hmu[0]),
            "ratio_stat_lo": float(ratio_stat[0]), "certifies": ok}
        print(f"  {name} [{kitname}]: CR_HU in "
              f"[{float(CR[0]):+.9e}, {float(CR[1]):+.9e}], "
              f"CR - H/25 > {float(ratio_stat[0]):+.3e}, "
              f"H > {float(Hmu[0]):.4f} "
              f"{'CERTIFIED' if ok else 'NOT CERTIFIED'}", flush=True)
        if not ok:
            raise SystemExit(f"certification failed: {name} {kitname}")
    OUT["instances"].append(row)


# the 033 campaign witness (deterministic descent output, see
# data/hu_attack.json row floor:random0 at cap 0.45), already permuted
# by its violating order (0, 1, 3, 2)
WITNESS = {int(s, 2): w for s, w in {
    "1110": 0.13062344029259873, "1101": 0.13171482140927657,
    "0010": 0.1459089053302761, "1011": 0.14277544270070638,
    "0101": 0.17540846236566682, "0000": 0.2735689279014755}.items()}
W_ORDER = (0, 1, 3, 2)


def certify_violation():
    n = 4
    mu_p = {}
    for a, w in WITNESS.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << W_ORDER[i]
        mu_p[b] = mu_p.get(b, 0.0) + w
    mu = {a: Fraction(w).limit_denominator(10 ** 7)
          for a, w in mu_p.items()}
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    marg = max(sum(w for a, w in mu.items() if (a >> i) & 1)
               for i in range(n))
    assert marg < Fraction(9, 20)
    zs = hu_cells(n, mu)
    row = {"instance": "033-witness-order-0132", "n": n, "kits": {}}
    for kitname, log2fn in (("A_digit_extraction", log2_A),
                            ("B_atanh_series", log2_B)):
        Hmu = (Fraction(0), Fraction(0))
        for w in mu.values():
            Hmu = ivl_add(Hmu, ivl_scale(Fraction(-1),
                                         xlog2x_ivl(w, log2fn)))
        S = (Fraction(0), Fraction(0))
        for w, z in zs:
            if 0 < z < 1:
                S = ivl_add(S, ivl_scale(w, h_ivl(z, log2fn)))
        CR = (S[0] - Hmu[1], S[1] - Hmu[0])
        ok = CR[1] < 0 and Hmu[0] > HALF
        row["kits"][kitname] = {"CR_lo": float(CR[0]),
                                "CR_hi": float(CR[1]),
                                "H_lo": float(Hmu[0]), "certifies": ok}
        print(f"  C4 witness [{kitname}]: CR_HU in "
              f"[{float(CR[0]):+.9e}, {float(CR[1]):+.9e}] "
              f"(marg {float(marg):.6f} < 9/20 exact) "
              f"{'CERTIFIED < 0' if ok else 'NOT CERTIFIED'}", flush=True)
        if not ok:
            raise SystemExit("violation certification failed")
    OUT["instances"].append(row)


def main():
    print("Attempt 033: single-kit certification of CR_HU at the "
          "sharpest floors")
    d = json.loads((DATA / "cr_attack.json").read_text())
    for row in d["rows"]:
        if row["seed"] not in ("mmabskill", "windowkill"):
            continue
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        certify_instance(row["seed"], row["n"], mu)
    certify_violation()
    (DATA / "hu_certify.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    print("checkpoint: data/hu_certify.json")


if __name__ == "__main__":
    main()
