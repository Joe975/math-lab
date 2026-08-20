#!/usr/bin/env python3
"""Skeptic pass on uc_hu_rollcensus.py (attempt 038).  Stance: refute.

Independence: every re-evaluated CR goes through the 037-skeptic stack
(nats-based rule re-implementations; permute_mu -> half_union_pairs ->
cr_eval), which shares no evaluation code with the engine's hu_cr_seq.
The census instances are regenerated deterministically from the seeds
(same generator -- what is independently re-derived is every VALUE).

  S1  every P1 summary field recomputes from its stored rows
  S2  census spot re-evaluation: every 30th instance regenerated and
      roll/canonical/best CR re-derived independently
  S3  every P2/P3 descent endpoint: rollout sequence re-derived by the
      independent rule, floor re-evaluated independently; global
      floors and violation lists consistent

Usage: python uc_hu_rollcensus_skeptic.py   (exit 0 iff nothing refuted)
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
from uc_hu_rollcensus import gen_instance
from uc_hu_order2_skeptic import cr_indep, my_rule

DATA = HERE.parent / "data"
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def main():
    out = json.loads((DATA / "hu_rollcensus.json").read_text())

    print("S1. P1 summaries recompute from stored rows:")
    for key, blk in out.items():
        if not key.startswith("P1_"):
            continue
        rows = blk["rows"]
        ok = (blk["instances"] == len(rows)
              and blk["roll_best"] == sum(1 for r in rows if r["roll_is_best"])
              and blk["worst_roll_rank"] == min(r["roll_rank"] for r in rows)
              and abs(blk["worst_gap_over_H"]
                      - max(r["gap_over_H"] for r in rows)) < 1e-15
              and blk["roll_negative"] == sum(1 for r in rows if r["roll_neg"])
              and blk["any_order_negative"] == sum(1 for r in rows
                                                   if r["any_neg"])
              and blk["canon_below_roll"] == sum(
                  1 for r in rows if r["canon_rank"] < r["roll_rank"]))
        check(f"{key}: all summary fields", ok)

    print("S2. Census spot re-evaluation (every 30th instance, "
          "independent stack):")
    for key, blk in out.items():
        if not key.startswith("P1_"):
            continue
        n = int(key.split("_")[1][1:])
        cap = float(key.split("cap")[1])
        rng = random.Random(938000 + n * 1000 + int(cap * 1000))
        N = blk["instances"]
        regen = []
        tried = 0
        while len(regen) < N and tried < 20 * N:
            tried += 1
            mu = gen_instance(rng, n, cap)
            H = -sum(w * math.log2(w) for w in mu.values() if w > 0)
            if H < 0.2:
                continue
            regen.append(mu)
        bad = []
        for idx in range(0, N, 30):
            mu = regen[idx]
            row = blk["rows"][idx]
            best = max(cr_indep(n, mu, s)[0]
                       for s in itertools.permutations(range(n)))
            rollcr, H = cr_indep(n, mu, my_rule(n, mu, "roll"))
            cancr, _ = cr_indep(n, mu, my_rule(n, mu, "can"))
            if (abs(best - row["best_CR"]) > 1e-9
                    or abs(rollcr - row["roll_CR"]) > 1e-9
                    or abs(H - row["H"]) > 1e-9):
                bad.append((idx, rollcr, row["roll_CR"]))
        check(f"{key}: {len(range(0, N, 30))} spot instances reproduce "
              "(best, roll, H)", not bad, f"mismatches {bad[:2]}")

    print("S3. P2/P3 endpoints, independent rule + evaluator:")
    for key, blk in out.items():
        if not (key.startswith("P2_") or key.startswith("P3_")):
            continue
        bad = []
        for r in blk["rows"]:
            n = r["n"]
            mu = {int(s, 2): w for s, w in r["mu"].items()}
            cr, H = cr_indep(n, mu, my_rule(n, mu, "roll"))
            if abs(cr / H - r["floor"]) > 1e-8:
                bad.append((r["start"], cr / H, r["floor"]))
        floors = [r["floor"] for r in blk["rows"]]
        check(f"{key}: every endpoint floor reproduces", not bad,
              f"mismatches {bad[:2]}")
        check(f"{key}: global floor + violations consistent",
              abs(min(floors) - blk["global_floor"]) < 1e-12
              and sum(1 for f in floors if f < 0) == len(blk["violations"]))

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
