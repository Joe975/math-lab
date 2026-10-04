#!/usr/bin/env python3
"""Attempt 045 (044 lead 1): the direct attack on BEST-ORDER HU at caps
0.495 / 0.497 / 0.499 -- the line's last unrefuted positivity statement.

After 044, the ladder reads: canonical dies at 0.49 (035), rollout dies
at 0.495 (044), and only the best-order (oracle) form survives every
instance of record.  This file attacks it head-on where the fixed rules
died: weight-descents minimizing the max-over-ALL-orders CR/H, n <= 5
(full enumeration inside the objective), seeded hostile -- the 044 kill
witness and its D_bestorder siblings, the 042 equality family perturbed
off-family, and the standing starts.  Violations are flagged BOTH ways
(044's lesson baked in):

    KILL : best_ratio < 0            -> ends the order-quantifier
                                        program below 1/2 entirely
    SHARP: 0 <= best_ratio < c*(own fmax)
                                     -> refutes the sharp-constant
                                        (HU-TAX, best-order) form

Standard library only; deterministic.
Usage: python uc_hu_bestorder.py [--fast]
Checkpoint: ../data/hu_bestorder.json
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import best_order_ratio, descend, h
from uc_hu_attack import starts
from uc_hu_blocks import block_tensor

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def own_fmax(n, mu):
    tot = sum(mu.values())
    return max(sum(w / tot for a, w in mu.items() if (a >> i) & 1)
               for i in range(n))


def seeds_for(cap):
    out = []
    d2 = json.loads((DATA / "hu_order2.json").read_text())
    kill = [r for r in d2["D_bestorder_cap_0.497"]["rows"]
            if r["start"] == "floor:windowkill"][0]
    out.append(("044kill", kill["n"],
                {int(s, 2): w for s, w in kill["mu"].items()}))
    for key in ("D_bestorder_cap_0.49", "D_bestorder_cap_0.497"):
        for r in sorted(d2[key]["rows"], key=lambda x: x["floor"])[:3]:
            if r["start"] == "floor:windowkill" and key.endswith("497"):
                continue
            out.append((f"D:{key.split('_')[-1]}:{r['start']}", r["n"],
                        {int(s, 2): w for s, w in r["mu"].items()}))
    # equality-family members perturbed OFF-family (043's finding: the
    # adversary stalls near the family -- start it there, off-balance)
    rng = random.Random(944500)
    for blocks in ((("d", 2), ("d", 2)), (("d", 2), ("d", 3)),
                   (("d", 4),)):
        p = cap - 0.002
        n, mu = block_tensor(p, list(blocks))
        pert = {a: w * math.exp(rng.uniform(-0.35, 0.35))
                for a, w in mu.items()}
        out.append((f"family{blocks}-pert", n, pert))
    for name, n, mu in starts(cap):
        if n <= 5:
            out.append((name, n, mu))
    return out


def main():
    caps = (0.495, 0.497, 0.499) if not FAST else (0.497,)
    for cap in caps:
        rows, gmin = [], 1e9
        kills, sharps = [], []
        for name, n, mu in seeds_for(cap):
            r = descend(name, n, mu, cap, best_order_ratio,
                        rounds_max=(80 if not FAST else 10),
                        step_floor=0.04)
            if r is None:
                continue
            m = {int(s, 2): w for s, w in r["mu"].items()}
            fmax = own_fmax(r["n"], m)
            r["fmax"] = fmax
            r["own_margin"] = r["floor"] - cstar(fmax)
            rows.append(r)
            gmin = min(gmin, r["floor"])
            if r["floor"] < -1e-12:
                kills.append(r)
                log(f"    {name:22s}: KILL best_ratio {r['floor']:+.6e}")
            elif r["own_margin"] < -1e-12:
                sharps.append(r)
                log(f"    {name:22s}: SHARP viol margin "
                    f"{r['own_margin']:+.3e} (ratio {r['floor']:+.6e}, "
                    f"fmax {fmax:.5f})")
            else:
                log(f"    {name:22s}: floor {r['floor']:+.6e}, own margin "
                    f"{r['own_margin']:+.3e}")
        log(f"  best-order cap {cap}: {len(rows)} starts, global floor "
            f"{gmin:+.6e}, {len(kills)} kills, {len(sharps)} sharp "
            f"violations [c*(cap) {cstar(cap):+.6e}]")
        OUT[f"cap_{cap}"] = {"rows": rows, "global_floor": gmin,
                             "kills": kills, "sharp_violations": sharps}
    (DATA / "hu_bestorder.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_bestorder.json")


if __name__ == "__main__":
    main()
