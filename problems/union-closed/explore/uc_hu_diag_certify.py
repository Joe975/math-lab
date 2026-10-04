#!/usr/bin/env python3
"""Attempt 056 (055 lead 3): the n = 6 essential equality endpoint,
certified exactly.

055's constrained anneal converged at n = 6, cap 0.49 to a diagonal
{empty: 1-p, full: p} with p = 0.3826 and margin +6.26e-13 -- the first
time an adversary on this route landed on the equality family in
ESSENTIAL form at n >= 6.  Float evidence of an equality is weak
evidence, so this certifies it.

042's identity for a diagonal at parameter p: every HU cell has
x = y = 1-p, so z = z*(p) = max(1/2, 1-2p) and

    CR_HU = h(z*(p)) - h(p) = c*(p) h(p),    H(mu) = h(p),

independently of n -- the diagonal is an n = 1 instance in disguise as
far as CR is concerned, yet it is genuinely n-dimensional (every
marginal equals p).  So the margin CR/H - c*(p) is EXACTLY 0.

Certified statements, each kit alone (029 standard):
  M1  CR_HU - c*(p) H(mu) = 0 to within the kit enclosures, at the
      rationalized anneal endpoint
  M2  CR_HU > 0 and H > 1/2, exact
  M3  the endpoint is essential: every coordinate marginal equals p
      exactly, in rationals (so this is an n = 6 statement, not an
      embedded smaller one)

Usage: python uc_hu_diag_certify.py
Checkpoint: ../data/hu_diag_certify.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_certify import (hu_cells, ivl_add, ivl_scale, xlog2x_ivl, h_ivl,
                           log2_A, log2_B)

DATA = HERE.parent / "data"
HALF = Fraction(1, 2)


def main():
    d = json.loads((DATA / "hu_anneal_essential.json").read_text())
    r = min(d["A_n6_cap_0.49"]["rows"], key=lambda z: z["floor"])
    n = r["n"]
    muF = {int(s, 2): w for s, w in r["mu"].items() if w > 1e-9}
    full = (1 << n) - 1
    assert set(muF) <= {0, full}, f"not a diagonal: {sorted(muF)}"
    p = Fraction(muF[full]).limit_denominator(10 ** 7)
    mu = {0: 1 - p, full: p}
    marg = [sum(w for a, w in mu.items() if (a >> i) & 1) for i in range(n)]
    print("Attempt 056: the n = 6 essential equality endpoint, certified")
    print(f"  diagonal at p = {float(p):.6f} = {p}")
    print(f"  M3: all {n} marginals equal p exactly: "
          f"{all(m == p for m in marg)}  (essential, so genuinely n = {n})")
    assert all(m == p for m in marg)
    assert Fraction(0) < p < HALF
    zstar = max(HALF, 1 - 2 * p)
    out = {"n": n, "p": float(p), "p_exact": str(p),
           "z_star": float(zstar), "essential": True, "kits": {}}
    for kitname, fn in (("A_digit_extraction", log2_A),
                        ("B_atanh_series", log2_B)):
        Hmu = (Fraction(0), Fraction(0))
        for w in mu.values():
            Hmu = ivl_add(Hmu, ivl_scale(Fraction(-1), xlog2x_ivl(w, fn)))
        S = (Fraction(0), Fraction(0))
        for w, z in hu_cells(n, mu):
            if 0 < z < 1:
                S = ivl_add(S, ivl_scale(w, h_ivl(z, fn)))
        CR = (S[0] - Hmu[1], S[1] - Hmu[0])
        # target: CR = h(z*) - h(p);  margin = CR - c*(p) H = CR - (h(z*) - h(p))
        hz = h_ivl(zstar, fn)
        hp = h_ivl(p, fn)
        tgt = (hz[0] - hp[1], hz[1] - hp[0])
        margin = (CR[0] - tgt[1], CR[1] - tgt[0])
        zero_in = margin[0] <= 0 <= margin[1]
        pos = CR[0] > 0 and Hmu[0] > HALF
        out["kits"][kitname] = {
            "CR_lo": float(CR[0]), "CR_hi": float(CR[1]),
            "target_lo": float(tgt[0]), "target_hi": float(tgt[1]),
            "margin_lo": float(margin[0]), "margin_hi": float(margin[1]),
            "H_lo": float(Hmu[0]), "equality_certified": zero_in,
            "positive_certified": pos}
        print(f"  [{kitname}] M1: CR - c*(p)H in "
              f"[{float(margin[0]):+.3e}, {float(margin[1]):+.3e}] "
              f"{'CONTAINS 0 -> equality certified' if zero_in else 'DOES NOT CONTAIN 0'}")
        print(f"  [{kitname}] M2: CR in [{float(CR[0]):+.9f}, "
              f"{float(CR[1]):+.9f}], H > {float(Hmu[0]):.4f} "
              f"{'CERTIFIED > 0' if pos else 'NOT CERTIFIED'}")
        if not (zero_in and pos):
            raise SystemExit(f"certification failed under {kitname}")
    (DATA / "hu_diag_certify.json").write_text(
        json.dumps(out, indent=1, default=float) + "\n")
    print("checkpoint: data/hu_diag_certify.json")


if __name__ == "__main__":
    main()
