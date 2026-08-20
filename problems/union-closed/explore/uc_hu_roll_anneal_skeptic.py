#!/usr/bin/env python3
"""Skeptic pass on uc_hu_roll_anneal.py (attempt 040).  Stance: refute.

Every endpoint measure in the checkpoint is re-scored with the
037-skeptic independent stack (nats-based rollout re-implementation;
permute_mu -> half_union_pairs -> cr_eval), which shares no evaluation
code with the engine.  Anneal trajectories themselves are not re-run
(the endpoints are the claims); regime membership and summary
consistency are re-checked exactly.

  S1  every endpoint: in-regime (max marginal < cap, H >= 0.2), floor
      value reproduces under the independent rollout + evaluator
  S2  global floors = min over rows; violation lists = negative rows
  S3  the embedded seeds' construction claim: embedding the 035 kill
      measures with independent near-cap coordinates keeps them
      in-regime (checked on the checkpoint rows named embed:*)

Usage: python uc_hu_roll_anneal_skeptic.py  (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import json
import math
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


def main():
    out = json.loads((DATA / "hu_roll_anneal.json").read_text())
    for key, blk in out.items():
        if not key.startswith("A_"):
            continue
        cap = float(key.split("cap_")[1])
        bad_val, bad_reg = [], []
        for r in blk["rows"]:
            n = r["n"]
            mu = {int(s, 2): w for s, w in r["mu"].items()}
            tot = sum(mu.values())
            m = {a: w / tot for a, w in mu.items()}
            fm = max(sum(w for a, w in m.items() if (a >> i) & 1)
                     for i in range(n))
            cr, H = cr_indep(n, m, my_rule(n, m, "roll"))
            if not (fm < cap and H >= 0.2):
                bad_reg.append((r["start"], fm, H))
            if abs(cr / H - r["floor"]) > 1e-8:
                bad_val.append((r["start"], cr / H, r["floor"]))
        floors = [r["floor"] for r in blk["rows"]]
        check(f"{key}: every endpoint in-regime", not bad_reg,
              f"violations {bad_reg[:2]}")
        check(f"{key}: every endpoint floor reproduces independently",
              not bad_val, f"mismatches {bad_val[:2]}")
        check(f"{key}: global floor + violation list consistent",
              abs(min(floors) - blk["global_floor"]) < 1e-12
              and sum(1 for f in floors if f < 0) == len(blk["violations"]))
        emb = [r for r in blk["rows"] if r["start"].startswith("embed:")]
        check(f"{key}: embedded hostile seeds present and in-regime",
              len(emb) >= 1 and all(
                  max(sum(w for a, w in
                          {int(s, 2): w2 for s, w2 in r["mu"].items()}.items()
                          if (a >> i) & 1) for i in range(r["n"])) < cap
                  for r in emb),
              f"{len(emb)} embeds")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
