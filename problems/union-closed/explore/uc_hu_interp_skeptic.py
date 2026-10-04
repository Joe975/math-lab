#!/usr/bin/env python3
"""Skeptic pass on uc_hu_interp.py (attempt 041).  Stance: refute.

  S1  small n (<= 8): sampled grid points rebuilt as atom measures and
      re-scored through the independent 037-skeptic stack
      (permute_mu -> half_union_pairs -> cr_eval), margins recomputed
      with an own nats-based c*; marginals checked exact
  S2  all n (incl. 16/32/64): the path endpoints are the two PROVED
      equality families, so margin(t=0) and margin(t=1) must vanish --
      a closed-form check of the DP at every n in the checkpoint
  S3  summary consistency: worst margins, negative counts, refined
      minima recompute from rows / from the DP

Usage: python uc_hu_interp_skeptic.py    (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2_skeptic import cr_indep
from uc_hu_interp import comps_A, comps_B, atoms_of, hu_cr_pmix

DATA = HERE.parent / "data"
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hn(p):
    return 0.0 if p <= 0 or p >= 1 else -(p * math.log(p)
                                           + (1 - p) * math.log(1 - p))


def cstar_nats(p):
    return (hn(max(0.5, 1 - 2 * p)) - hn(p)) / hn(p)


def main():
    out = json.loads((DATA / "hu_interp.json").read_text())
    paths = {"A_convex": comps_A, "B_prodmix": comps_B}

    print("S1. Small-n spot re-evaluation through the independent stack:")
    for pname, cf in paths.items():
        rows = out[pname]["rows"]
        small = [r for r in rows if r["n"] <= 8]
        bad = []
        for r in small[:: max(1, len(small) // 24)]:
            n, p, t = r["n"], r["p"], r["t"]
            mu = atoms_of(n, cf(p, t))
            dev = max(abs(sum(w for a, w in mu.items() if (a >> i) & 1) - p)
                      for i in range(n))
            cr, H = cr_indep(n, mu, tuple(range(n)))
            m = cr / H - cstar_nats(p)
            if abs(m - r["margin"]) > 1e-9 or dev > 1e-12:
                bad.append((n, p, t, m, r["margin"], dev))
        check(f"{pname}: sampled small-n margins reproduce independently "
              f"({len(small[::max(1, len(small)//24)])} points)",
              not bad, f"mismatches {bad[:2]}")

    print("S2. Endpoint identities at every n (the two equality families):")
    for pname in paths:
        rows = out[pname]["rows"]
        bad = [(r["n"], r["p"], r["t"], r["margin"]) for r in rows
               if r["t"] in (0.0, 1.0) and abs(r["margin"]) > 1e-9]
        check(f"{pname}: margin vanishes at t=0 (product) and t=1 "
              "(diagonal) for every n and p", not bad, f"{bad[:3]}")

    print("S3. Summary consistency:")
    for pname, cf in paths.items():
        blk = out[pname]
        rows = blk["rows"]
        wm = min(r["margin"] for r in rows)
        neg = sum(1 for r in rows if r["margin"] < -1e-12)
        ok = (abs(wm - blk["worst_margin"]) < 1e-15
              and neg == blk["negatives"])
        rm = blk["refined_min"]
        cr, H = hu_cr_pmix(rm["n"], cf(rm["p"], rm["t"]))
        ok2 = abs(cr / H - cstar_nats(rm["p"]) - rm["margin"]) < 1e-9
        check(f"{pname}: worst/negatives/refined-min consistent",
              ok and ok2,
              f"worst {wm:+.2e}, negatives {neg}, refined "
              f"{rm['margin']:+.2e}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
