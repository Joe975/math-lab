#!/usr/bin/env python3
"""Skeptic pass on uc_hu_n3_essential.py (attempt 051).  Stance: refute.

Independence: every CR is recomputed with the nats history-recursion
evaluator of the 050 skeptic (own code, no shared evaluation with
uc_hu_order2 or the engine), and the degeneracy claim is re-derived
from the committed 050 checkpoint rather than from the engine's report.

  V1  050's part-D endpoints really are degenerate (recompute the
      coordinate marginals and live-atom counts from the checkpoint)
  V2  every 051 essential endpoint really is essential AND its
      best-order own-constant margin reproduces
  V3  the certificate's float value reproduces independently
  V4  the essentiality constraint is not vacuous: random draws are
      rejected by it at a healthy rate, and the constrained floor is
      strictly above the degenerate one

Usage: python uc_hu_n3_essential_skeptic.py   (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import itertools
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_n3_census_skeptic import cr_of, cstar     # own evaluator
from uc_hu_rollcensus import gen_instance            # generator only

DATA = HERE.parent / "data"
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def norm(mu):
    tot = sum(mu.values())
    return {a: w / tot for a, w in mu.items() if w / tot > 1e-14}


def best_margin(mu):
    m = norm(mu)
    margs = [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(3)]
    fmax = max(margs)
    if not (0 < fmax < 0.5):
        return None
    best, H = None, None
    for s in itertools.permutations(range(3)):
        cr, H = cr_of(3, m, s)
        if best is None or cr > best:
            best = cr
    if H is None or H < 0.2:
        return None
    return best / H - cstar(fmax), margs, sum(1 for w in m.values()
                                              if w >= 1e-3)


def main():
    o = json.loads((DATA / "hu_n3_essential.json").read_text())
    census = json.loads((DATA / "hu_n3_census.json").read_text())

    print("V1. 050's part-D endpoints re-inspected from its checkpoint:")
    for capkey in ("D_cap_0.45", "D_cap_0.49"):
        r = min(census[capkey]["rows"], key=lambda z: z["own_margin"])
        m = norm({int(s, 2): w for s, w in r["mu"].items()})
        margs = [sum(w for a, w in m.items() if (a >> i) & 1)
                 for i in range(3)]
        live = sum(1 for w in m.values() if w >= 1e-3)
        degenerate = min(margs) < 0.05 or live < 4
        check(f"{capkey}: endpoint is degenerate (an n=1 face)",
              degenerate,
              f"marginals {[round(x,4) for x in margs]}, live atoms {live}")

    print("V2. 051 essential endpoints: essential, and margins reproduce:")
    for cap in ("0.45", "0.49"):
        blk = o.get(f"B_cap_{cap}")
        if not blk:
            continue
        bad = []
        for r in blk["rows"]:
            res = best_margin({int(s, 2): w for s, w in r["mu"].items()})
            if res is None:
                bad.append((r["start"], "out of regime"))
                continue
            margin, margs, live = res
            if min(margs) < 0.05 - 1e-9 or live < 4:
                bad.append((r["start"], f"not essential: {margs}, {live}"))
            elif abs(margin - r["floor"]) > 1e-7:
                bad.append((r["start"], f"{margin:+.6f} vs {r['floor']:+.6f}"))
        check(f"cap {cap}: all endpoints essential and reproducing",
              not bad, f"{bad[:2]}")
        floors = [r["floor"] for r in blk["rows"]]
        check(f"cap {cap}: min floor matches and is positive",
              abs(min(floors) - blk["min_margin"]) < 1e-12
              and min(floors) > 0,
              f"min {min(floors):+.6f}")

    print("V3. The certificate's float value, independently:")
    c = o["C_certificate"]
    blk = o[f"B_cap_{c['cap']}"]
    r = min(blk["rows"], key=lambda z: z["floor"])
    res = best_margin({int(s, 2): w for s, w in r["mu"].items()})
    check("certified endpoint reproduces and is positive",
          res is not None and abs(res[0] - c["float_floor"]) < 1e-7
          and all(k["certifies"] for k in c["kits"].values()),
          f"mine {res[0]:+.6e} vs recorded {c['float_floor']:+.6e}")

    print("V4. The essentiality constraint is not vacuous:")
    rng = random.Random(9001)
    rejected = 0
    total = 0
    for _ in range(3000):
        mu = norm(gen_instance(rng, 3, 0.49))
        total += 1
        margs = [sum(w for a, w in mu.items() if (a >> i) & 1)
                 for i in range(3)]
        live = sum(1 for w in mu.values() if w >= 1e-3)
        if min(margs) < 0.05 or live < 4:
            rejected += 1
    check("the constraint rejects a nontrivial share of random draws",
          0 < rejected < total,
          f"{rejected}/{total} rejected")
    degen = min(min(census[k]["rows"], key=lambda z: z["own_margin"])
                ["own_margin"] for k in ("D_cap_0.45", "D_cap_0.49"))
    ess = min(o[f"B_cap_{c}"]["min_margin"] for c in ("0.45", "0.49")
              if f"B_cap_{c}" in o)
    check("the essential floor is strictly above the degenerate one",
          ess > degen, f"essential {ess:+.3e} vs degenerate {degen:+.3e}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
