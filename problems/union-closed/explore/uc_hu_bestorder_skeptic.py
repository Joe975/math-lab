#!/usr/bin/env python3
"""Skeptic pass on uc_hu_bestorder.py (attempt 045).  Stance: refute.

Every endpoint's best-order ratio is re-derived by full enumeration
through the independent stack (permute_mu -> half_union_pairs ->
cr_eval; no code shared with the engine's hu_cr_seq), own max marginals
and c* margins recomputed, and each block's global floor / kill list /
sharp-violation list checked against its rows.

Usage: python uc_hu_bestorder_skeptic.py   (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2_skeptic import cr_indep

DATA = HERE.parent / "data"
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hbits(p):
    return 0.0 if p <= 0 or p >= 1 else -(p * math.log2(p)
                                          + (1 - p) * math.log2(1 - p))


def cstar(p):
    return (hbits(max(0.5, 1.0 - 2.0 * p)) - hbits(p)) / hbits(p)


def main():
    out = json.loads((DATA / "hu_bestorder.json").read_text())
    for key, blk in out.items():
        bad = []
        for r in blk["rows"]:
            n = r["n"]
            mu = {int(s, 2): w for s, w in r["mu"].items()}
            tot = sum(mu.values())
            m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
            best, H = None, None
            for s in itertools.permutations(range(n)):
                cr, H = cr_indep(n, m, s)
                if best is None or cr > best:
                    best = cr
            fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
                       for i in range(n))
            if (abs(best / H - r["floor"]) > 1e-8
                    or abs(fmax - r["fmax"]) > 1e-12
                    or abs((best / H - cstar(fmax)) - r["own_margin"])
                    > 1e-8):
                bad.append((r["start"], best / H, r["floor"]))
        floors = [r["floor"] for r in blk["rows"]]
        check(f"{key}: every endpoint (ratio, fmax, own margin) "
              "reproduces independently", not bad, f"mismatches {bad[:2]}")
        check(f"{key}: global floor, kills, sharp violations consistent",
              abs(min(floors) - blk["global_floor"]) < 1e-12
              and len(blk["kills"]) == sum(1 for f in floors
                                           if f < -1e-12)
              and len(blk["sharp_violations"]) == sum(
                  1 for r in blk["rows"]
                  if r["floor"] >= -1e-12 and r["own_margin"] < -1e-12))
    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
