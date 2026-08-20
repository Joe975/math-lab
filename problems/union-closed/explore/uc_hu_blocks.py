#!/usr/bin/env python3
"""Attempt 042: the block-tensor equality family of (HU-TAX), attacked
at its own boundary; plus two queued close-outs.

041 observed that CR_HU is additive over independent blocks, so every
tensor product of diagonals diag_m(p) and Bern(p) factors is a further
equality case of the corrected constant c*(p).  Parts:

  A  verify the equality claim numerically: block-tensors at n = 4/6
     (all block partitions of diagonals + mixed diag/Bernoulli
     tensors), margin == 0 under identity, rollout AND random
     interleaved orders (additivity needs no block-contiguity)
  B  attack the NEW equality surface: anneal + polish against rollout
     seeded AT block-tensor equality points (n = 4, 6; caps 0.49 and
     0.497) -- the sharpest falsification surface known after 041
  C  close-outs: (i) crash8 -- the one start never descended under
     rollout (039 correction) -- descended at caps 0.49/0.497;
     (ii) the 0.497 anneal re-seeded at the diagonal so that cap's
     floor is family-saturated (040 lead 3)

Standard library only; deterministic (seeds 951000+).
Usage: python uc_hu_blocks.py [--fast]
Checkpoint: ../data/hu_blocks.json
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
from uc_hu_order2 import hu_cr_seq, seq_roll, rule_ratio, descend, h
from uc_hu_roll_anneal import anneal, project
from uc_hu_attack import starts

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def tensor(mu1, n1, mu2):
    out = {}
    for a, w in mu1.items():
        for b, v in mu2.items():
            out[a | (b << n1)] = w * v
    return out


def block_tensor(p, blocks):
    """blocks: list of ('d', m) diagonal of size m or ('b', 1) Bern(p)."""
    mu = {0: 1.0}
    n = 0
    for kind, m in blocks:
        if kind == "d":
            f = {0: 1.0 - p, (1 << m) - 1: p}
        else:
            f = {0: 1.0 - p, 1: p}
        mu = tensor(mu, n, f) if n else f
        n += m
    return n, mu


def part_A():
    log("A. Block-tensor equality verification (margin == 0, any order):")
    rng = random.Random(951000)
    cases = [
        (0.49, [("d", 2), ("d", 2)]),
        (0.49, [("d", 3), ("d", 3)]),
        (0.49, [("d", 2), ("d", 2), ("d", 2)]),
        (0.49, [("d", 4), ("b", 1), ("b", 1)]),
        (0.45, [("d", 2), ("b", 1), ("d", 3)]),
        (0.38271, [("d", 2), ("d", 2)]),
        (0.30, [("d", 3), ("b", 1), ("d", 2)]),   # sub-1/4? p=0.30 >= 1/4
        (0.20, [("d", 2), ("d", 2)]),             # clamp branch of c*
    ]
    rows, worst = [], 0.0
    for p, blocks in cases:
        n, mu = block_tensor(p, blocks)
        margins = []
        orders = [tuple(range(n)), seq_roll(n, mu)]
        for _ in range(3):
            o = list(range(n))
            rng.shuffle(o)
            orders.append(tuple(o))
        for seq in orders:
            cr, H = hu_cr_seq(n, mu, seq)
            margins.append(cr / H - cstar(p))
        m = max(abs(x) for x in margins)
        worst = max(worst, m)
        rows.append({"p": p, "blocks": blocks, "n": n,
                     "max_abs_margin": m, "orders_tested": len(orders)})
        log(f"  p={p} blocks={blocks}: |margin| <= {m:.2e} over "
            f"{len(orders)} orders")
    OUT["A_equality"] = {"rows": rows, "worst_abs_margin": worst}
    log(f"  worst |margin| across all cases/orders: {worst:.2e}")


def part_B():
    log()
    log("B. Anneal + polish seeded AT block-tensor equality points:")
    caps_blocks = [
        (0.49, 4, [("d", 2), ("d", 2)]),
        (0.49, 6, [("d", 3), ("d", 3)]),
        (0.49, 6, [("d", 2), ("d", 2), ("d", 2)]),
        (0.497, 4, [("d", 2), ("d", 2)]),
        (0.497, 6, [("d", 2), ("d", 2), ("d", 2)]),
    ] if not FAST else [(0.49, 4, [("d", 2), ("d", 2)])]
    steps = 4000 if not FAST else 300
    rows, gmin = [], 1e9
    for cap, n, blocks in caps_blocks:
        p = cap - 0.003
        _, mu0 = block_tensor(p, blocks)
        rng = random.Random(951100 + n + int(cap * 1000))
        r = anneal(f"blocks{blocks}", n, mu0, cap, steps, rng)
        if r is None:
            continue
        av, amu = r
        d = descend(f"blocks{blocks}+polish", n, amu, cap,
                    lambda n_, m_, c_: rule_ratio(n_, m_, c_, seq_roll),
                    rounds_max=80 if not FAST else 15)
        floor = d["floor"] if d is not None and d["floor"] < av else av
        endmu = d["mu"] if d is not None and d["floor"] < av else {
            format(a, f"0{n}b"): w / sum(amu.values())
            for a, w in amu.items()}
        rows.append({"start": str(blocks), "n": n, "cap": cap,
                     "seed_margin": cstar(p) - cstar(cap),
                     "anneal_best": av, "floor": floor, "mu": endmu})
        gmin = min(gmin, floor)
        log(f"  cap {cap} n={n} {blocks}: anneal {av:+.6f}, "
            f"polished {floor:+.6f} [extremal {cstar(cap):+.6f}]")
    viol = [r for r in rows if r["floor"] < 0]
    OUT["B_attack"] = {"rows": rows, "global_floor": gmin,
                       "violations": viol}
    log(f"  global floor {gmin:+.6f}, {len(viol)} violations")


def part_C():
    log()
    log("C. Close-outs:")
    # (i) crash8 under rollout (the start 037/038 never descended)
    rows = []
    for cap in ((0.49, 0.497) if not FAST else (0.49,)):
        for name, n, mu in starts(cap):
            if name != "crash8":
                continue
            r = descend("crash8", n, mu, cap,
                        lambda n_, m_, c_: rule_ratio(n_, m_, c_, seq_roll),
                        rounds_max=40 if not FAST else 10, step_floor=0.05)
            if r is not None:
                rows.append(r)
                log(f"  crash8 cap {cap}: floor {r['floor']:+.6f}")
    OUT["C_crash8"] = {"rows": rows,
                       "violations": [r for r in rows if r["floor"] < 0]}
    # (ii) 0.497 diagonal saturation
    if not FAST:
        cap = 0.497
        p = cap - 0.003
        rng = random.Random(951497)
        r = anneal("diag-seed", 6, {0: 1 - p, 63: p}, cap, 2000, rng)
        av, amu = r
        d = descend("diag-seed+polish", 6, amu, cap,
                    lambda n_, m_, c_: rule_ratio(n_, m_, c_, seq_roll),
                    rounds_max=60)
        floor = d["floor"] if d is not None and d["floor"] < av else av
        endmu = d["mu"] if d is not None and d["floor"] < av else {
            format(a, "06b"): w / sum(amu.values()) for a, w in amu.items()}
        diagval = cstar(p)
        OUT["C_sat497"] = {"floor": floor, "diag_value": diagval,
                           "extremal": cstar(cap), "mu": endmu,
                           "violation": floor < 0}
        log(f"  0.497 diagonal-seeded: floor {floor:+.8f} "
            f"[diag(0.494) value {diagval:+.8f}, extremal "
            f"{cstar(cap):+.8f}]")


def main():
    part_A()
    part_B()
    part_C()
    (DATA / "hu_blocks.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_blocks.json")


if __name__ == "__main__":
    main()
