#!/usr/bin/env python3
"""Skeptic pass on uc_hu_blocks.py (attempt 042).  Stance: refute.

Every value re-derived through the independent 037-skeptic stack
(nats rollout re-implementation; permute_mu -> half_union_pairs ->
cr_eval); block-tensor measures rebuilt from the blocks spec by an own
tensor routine; c* recomputed in nats.

  S1  part A: every block-tensor case rebuilt, margin re-derived under
      identity + rollout + 2 fresh random orders -- |margin| < 1e-12
  S2  parts B/C: every endpoint floor reproduces; violation lists and
      global floors consistent; endpoints in-regime
  S3  C_sat497: the floor equals the diag(cap-0.003) closed form

Usage: python uc_hu_blocks_skeptic.py    (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2_skeptic import cr_indep, my_rule

DATA = HERE.parent / "data"
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hn(p):
    return 0.0 if p <= 0 or p >= 1 else -(p * math.log(p)
                                           + (1 - p) * math.log(1 - p))


def cstar(p):
    return (hn(max(0.5, 1 - 2 * p)) - hn(p)) / hn(p)


def build(p, blocks):
    mu = {0: 1.0}
    n = 0
    for kind, m in blocks:
        f = ({0: 1 - p, (1 << m) - 1: p} if kind == "d"
             else {0: 1 - p, 1: p})
        mu = {a | (b << n): w * v for a, w in mu.items()
              for b, v in f.items()}
        n += m
    return n, mu


def main():
    out = json.loads((DATA / "hu_blocks.json").read_text())
    rng = random.Random(4242)

    print("S1. Block-tensor equality, own build + independent evaluator:")
    bad = []
    for r in out["A_equality"]["rows"]:
        p = r["p"]
        blocks = [(k, m) for k, m in r["blocks"]]
        n, mu = build(p, blocks)
        orders = [tuple(range(n)), my_rule(n, mu, "roll")]
        for _ in range(2):
            o = list(range(n))
            rng.shuffle(o)
            orders.append(tuple(o))
        for seq in orders:
            cr, H = cr_indep(n, mu, seq)
            if abs(cr / H - cstar(p)) > 1e-12:
                bad.append((p, blocks, seq, cr / H - cstar(p)))
    check("every case: margin == c*(p) to 1e-12 under 4 orders each",
          not bad, f"violations {bad[:2]}")

    print("S2. Attack + close-out endpoints:")
    for key in ("B_attack", "C_crash8"):
        blk = out[key]
        bad_v = []
        for r in blk["rows"]:
            n = r["n"]
            mu = {int(s, 2): w for s, w in r["mu"].items()}
            tot = sum(mu.values())
            m = {a: w / tot for a, w in mu.items()}
            fm = max(sum(w for a, w in m.items() if (a >> i) & 1)
                     for i in range(n))
            cr, H = cr_indep(n, m, my_rule(n, m, "roll"))
            if abs(cr / H - r["floor"]) > 1e-8 or fm >= r["cap"]:
                bad_v.append((r["start"], cr / H, r["floor"], fm))
        check(f"{key}: every endpoint reproduces, in-regime", not bad_v,
              f"mismatches {bad_v[:2]}")
        if "global_floor" in blk:
            floors = [r["floor"] for r in blk["rows"]]
            check(f"{key}: global floor + violations consistent",
                  abs(min(floors) - blk["global_floor"]) < 1e-12
                  and sum(1 for f in floors if f < 0)
                  == len(blk["violations"]))

    print("S3. 0.497 saturation:")
    s = out["C_sat497"]
    mu = {int(k, 2): w for k, w in s["mu"].items()}
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items()}
    cr, H = cr_indep(6, m, my_rule(6, m, "roll"))
    check("floor reproduces and equals diag(0.494) closed form",
          abs(cr / H - s["floor"]) < 1e-8
          and abs(s["floor"] - cstar(0.494)) < 1e-7
          and not s["violation"],
          f"floor {s['floor']:+.8f} vs c*(0.494) {cstar(0.494):+.8f}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
