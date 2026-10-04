#!/usr/bin/env python3
"""Attempt 044 (043's standing audit caveat): every stored endpoint of
the 2026-08-20 window, re-scored against its OWN constant.

The window's engines flagged violations by CR_HU/H < 0.  The standing
conjecture (HU-TAX, rollout form) claims more: CR_HU/H >= c*(fmax)
with fmax the instance's own maximum marginal (034's lesson -- compare
against the instance's own constant, not the campaign cap).  043
re-checked 042's endpoints; this file sweeps everything else the
window stored, plus the regenerated 1200-instance census:

  A  every endpoint mu stored in hu_order2.json (C_*, C2_*,
     D_bestorder_*), hu_rollcensus.json (P2_*, P3_*),
     hu_roll_anneal.json (A_*), hu_blocks.json (B_attack, C_crash8,
     C_sat497), and hu_canon.json (B_* rows, re-scored under ROLLOUT
     -- 035's canonical endpoints are the hardest instances of record)
  B  the 038 census, regenerated deterministically (seeds as
     committed), all 1200 instances: margin vs c*(fmax) under rollout

Any negative margin is a (HU-TAX, rollout-form) violation and the
sweep's headline; otherwise the headline is the minimum margin and
where it sits.  D_bestorder rows are scored under rollout too (the
conjecture's object), with the best-order value implied >= it.

Standard library only; deterministic.
Usage: python uc_hu_ownconst.py [--fast]
Checkpoint: ../data/hu_ownconst.json
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import hu_cr_seq, seq_roll, h
from uc_hu_rollcensus import gen_instance
import random

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def own_margin(n, mu):
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-15}
    fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
               for i in range(n))
    if not (0.0 < fmax < 0.5):
        return None
    cr, H = hu_cr_seq(n, m, seq_roll(n, m))
    if H < 0.2:
        return None
    return cr / H - cstar(fmax), fmax, cr / H


def part_A():
    log("A. Stored endpoints vs their own constant (rollout):")
    sources = [
        ("hu_order2.json", ("C_", "C2_", "D_")),
        ("hu_rollcensus.json", ("P2_", "P3_")),
        ("hu_roll_anneal.json", ("A_",)),
        ("hu_blocks.json", ("B_attack", "C_crash8")),
        ("hu_canon.json", ("B_cap",)),
    ]
    rows = []
    for fname, prefixes in sources:
        d = json.loads((DATA / fname).read_text())
        for key, blk in d.items():
            if not any(key.startswith(p) for p in prefixes):
                continue
            for r in blk.get("rows", blk.get("violations", [])):
                if "mu" not in r:
                    continue
                n = r["n"]
                mu = {int(s, 2): w for s, w in r["mu"].items()}
                res = own_margin(n, mu)
                if res is None:
                    continue
                mg, fmax, ratio = res
                rows.append({"src": f"{fname}:{key}:{r.get('start', '?')}",
                             "n": n, "fmax": fmax, "ratio": ratio,
                             "own_margin": mg})
    d = json.loads((DATA / "hu_blocks.json").read_text())
    s = d["C_sat497"]
    mu = {int(k, 2): w for k, w in s["mu"].items()}
    res = own_margin(6, mu)
    if res is not None:
        mg, fmax, ratio = res
        rows.append({"src": "hu_blocks.json:C_sat497", "n": 6,
                     "fmax": fmax, "ratio": ratio, "own_margin": mg})
    neg = [r for r in rows if r["own_margin"] < -1e-12]
    mn = min(rows, key=lambda r: r["own_margin"])
    log(f"  {len(rows)} stored endpoints; {len(neg)} own-constant "
        f"violations; min margin {mn['own_margin']:+.3e} at {mn['src']} "
        f"(fmax {mn['fmax']:.5f})")
    for r in sorted(rows, key=lambda x: x["own_margin"])[:6]:
        log(f"    {r['own_margin']:+.3e}  fmax {r['fmax']:.5f}  {r['src']}")
    OUT["A_endpoints"] = {"rows": rows, "negatives": len(neg),
                          "min_margin": mn["own_margin"],
                          "min_at": mn["src"]}


def part_B():
    log()
    log("B. The 038 census, regenerated, vs own constants (rollout):")
    N = 300 if not FAST else 30
    allrows = []
    for n, cap in ((4, 0.45), (4, 0.49), (5, 0.45), (5, 0.49)):
        rng = random.Random(938000 + n * 1000 + int(cap * 1000))
        got, tried = 0, 0
        while got < N and tried < 20 * N:
            tried += 1
            mu = gen_instance(rng, n, cap)
            H = -sum(w * math.log2(w) for w in mu.values() if w > 0)
            if H < 0.2:
                continue
            got += 1
            res = own_margin(n, mu)
            if res is None:
                continue
            mg, fmax, ratio = res
            allrows.append({"n": n, "cap": cap, "fmax": fmax,
                            "own_margin": mg})
        sub = [r for r in allrows if r["n"] == n and r["cap"] == cap]
        mn = min(sub, key=lambda r: r["own_margin"])
        log(f"  n={n} cap {cap}: {len(sub)} instances, "
            f"{sum(1 for r in sub if r['own_margin'] < -1e-12)} "
            f"violations, min margin {mn['own_margin']:+.4e} "
            f"(fmax {mn['fmax']:.4f})")
    neg = [r for r in allrows if r["own_margin"] < -1e-12]
    mn = min(allrows, key=lambda r: r["own_margin"])
    OUT["B_census"] = {"count": len(allrows), "negatives": len(neg),
                       "min_margin": mn["own_margin"],
                       "min_at": {k: mn[k] for k in ("n", "cap", "fmax")},
                       "rows": allrows}
    log(f"  total {len(allrows)}, {len(neg)} violations, min "
        f"{mn['own_margin']:+.4e}")


def main():
    part_A()
    part_B()
    (DATA / "hu_ownconst.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_ownconst.json")


if __name__ == "__main__":
    main()
