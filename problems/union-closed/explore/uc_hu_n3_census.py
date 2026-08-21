#!/usr/bin/env python3
"""Attempt 050 (045 lead 1 / 047 lead 3): n = 3 -- the size the HU line
has never examined -- under the own-constant standard.

Everything on this route has been measured at n = 2 (046-048) or at
n >= 4 (033/035/037/038/040/044/045).  n = 3 was skipped, and it is the
first size where a SECOND interaction layer exists, so it is where the
n = 2 structure either generalises or breaks.

Parts:
  A  own-constant census at n = 3: random in-regime measures, FULL
     6-order enumeration, margins CR/H - c*(own fmax) for the
     best order, the rollout order and the identity order
  B  the same at n = 4 (24 orders) for contrast, so the n = 3 numbers
     can be read against a size the line already knows
  C  is the order quantifier vacuous at n = 3?  (046 proved it vacuous
     at n = 2; 033/035/044 proved it essential at n >= 4)  Measured as:
     how often does the identity order alone already meet the bound
  D  descents against the best-order margin at n = 3, caps 0.45/0.49

Violations are flagged against c*(own max marginal) throughout -- the
044 standard, not the weaker CR/H < 0.

Standard library only; deterministic.
Usage: python uc_hu_n3_census.py [--fast]
Checkpoint: ../data/hu_n3_census.json
"""
from __future__ import annotations

import itertools
import json
import math
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import hu_cr_seq, seq_roll, h, descend, best_order_ratio
from uc_hu_rollcensus import gen_instance

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def analyse(n, mu):
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
    if not m:
        return None
    fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
               for i in range(n))
    if not (0 < fmax < 0.5):
        return None
    vals = {}
    H = None
    for s in itertools.permutations(range(n)):
        cr, H = hu_cr_seq(n, m, s)
        vals[s] = cr
    if H is None or H < 0.2:
        return None
    c = cstar(fmax)
    best = max(vals.values())
    ident = vals[tuple(range(n))]
    roll = hu_cr_seq(n, m, seq_roll(n, m))[0]
    return {"fmax": fmax, "H": H, "cstar": c,
            "m_best": best / H - c, "m_ident": ident / H - c,
            "m_roll": roll / H - c,
            "worst_order_cr": min(vals.values()),
            "n_neg_orders": sum(1 for v in vals.values() if v < 0)}


def census(n, cap, count, seed):
    rng = random.Random(seed)
    rows = []
    tried = 0
    while len(rows) < count and tried < 40 * count:
        tried += 1
        r = analyse(n, gen_instance(rng, n, cap))
        if r is not None:
            rows.append(r)
    return rows


def summarise(tag, rows):
    def mn(k):
        return min(r[k] for r in rows)
    viol = {k: sum(1 for r in rows if r[k] < -1e-12)
            for k in ("m_best", "m_ident", "m_roll")}
    log(f"  {tag}: {len(rows)} instances | min own-constant margin: "
        f"best {mn('m_best'):+.5f}, rollout {mn('m_roll'):+.5f}, "
        f"identity {mn('m_ident'):+.5f}")
    log(f"      violations: best {viol['m_best']}, rollout "
        f"{viol['m_roll']}, identity {viol['m_ident']}; instances where "
        f"SOME order has CR < 0: "
        f"{sum(1 for r in rows if r['n_neg_orders'] > 0)}")
    return {"instances": len(rows),
            "min_best": mn("m_best"), "min_roll": mn("m_roll"),
            "min_ident": mn("m_ident"), "violations": viol,
            "some_order_negative":
                sum(1 for r in rows if r["n_neg_orders"] > 0)}


def part_A():
    log("A. n = 3 own-constant census (full 6-order enumeration):")
    OUT["A_n3"] = {}
    for cap in (0.45, 0.49):
        rows = census(3, cap, 3000 if not FAST else 200, 5000 + int(cap * 1000))
        OUT["A_n3"][f"cap_{cap}"] = summarise(f"n=3 cap {cap}", rows)


def part_B():
    log()
    log("B. n = 4 contrast (24 orders):")
    OUT["B_n4"] = {}
    for cap in (0.45, 0.49):
        rows = census(4, cap, 1500 if not FAST else 100, 5100 + int(cap * 1000))
        OUT["B_n4"][f"cap_{cap}"] = summarise(f"n=4 cap {cap}", rows)


def part_C():
    log()
    log("C. Is the order quantifier vacuous at n = 3?")
    for n in (3, 4):
        rows = census(n, 0.49, 2000 if not FAST else 150, 5200 + n)
        ident_ok = sum(1 for r in rows if r["m_ident"] >= -1e-12)
        roll_ok = sum(1 for r in rows if r["m_roll"] >= -1e-12)
        log(f"  n={n}: identity order meets the own-constant bound on "
            f"{ident_ok}/{len(rows)}; rollout on {roll_ok}/{len(rows)}")
        OUT[f"C_n{n}"] = {"instances": len(rows), "identity_ok": ident_ok,
                          "rollout_ok": roll_ok}
    log("  (046 proved the quantifier vacuous at n = 2; 033/035/044 "
        "proved it essential at n >= 4 -- n = 3 is the boundary case)")


def part_D():
    log()
    log("D. Descents against the BEST-ORDER margin at n = 3:")
    for cap in (0.45, 0.49):
        rng = random.Random(5300 + int(cap * 1000))
        rows, gmin = [], 1e9
        for j in range(8 if not FAST else 2):
            mu = gen_instance(rng, 3, cap)
            r = descend(f"n3seed{j}", 3, mu, cap, best_order_ratio,
                        rounds_max=80, step_floor=0.04)
            if r is None:
                continue
            m = {int(s, 2): w for s, w in r["mu"].items()}
            tot = sum(m.values())
            m = {a: w / tot for a, w in m.items()}
            fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
                       for i in range(3))
            r["fmax"] = fmax
            r["own_margin"] = r["floor"] - cstar(fmax)
            rows.append(r)
            gmin = min(gmin, r["own_margin"])
        viol = [r for r in rows if r["own_margin"] < -1e-12]
        log(f"  cap {cap}: {len(rows)} descents, min own-constant margin "
            f"{gmin:+.6f}, {len(viol)} violations")
        OUT[f"D_cap_{cap}"] = {"rows": rows, "min_own_margin": gmin,
                               "violations": viol}


def main():
    part_A()
    part_B()
    part_C()
    part_D()
    (DATA / "hu_n3_census.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_n3_census.json")


if __name__ == "__main__":
    main()
