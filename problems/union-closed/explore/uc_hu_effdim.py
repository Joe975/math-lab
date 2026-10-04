#!/usr/bin/env python3
"""Attempt 052 (051 lead 1): effective-dimension audit of every descent
endpoint this route has recorded.

051 found that 050's n = 3 descents saturated their bound by escaping
to a degenerate face.  That raised an obvious question about every
earlier campaign, none of which was ever checked for support collapse.
This file answers it for all committed checkpoints.

Criterion (and why): a coordinate whose marginal is EXACTLY 0 is
absent, so the instance really lives at a smaller n and any "n = k
floor" quoting it overstates the dimension.  Atom count is NOT a valid
criterion -- a diagonal {empty, S} has two atoms and is genuinely
n-dimensional (it is one of the equality families) -- so this audit
uses marginals only, and reports "min marginal < 0.05" separately as
near-degenerate rather than folding it in.

Parts:
  A  every stored endpoint, by block: how many have an absent
     coordinate, how many are merely near-degenerate
  B  the HEADLINE floors -- the minimum-floor endpoint of each block,
     which is what the records quote -- with claimed n against
     effective n
  C  what this does and does not change

Standard library only; deterministic (pure audit of committed data).
Usage: python uc_hu_effdim.py
Checkpoint: ../data/hu_effdim.json
"""
from __future__ import annotations

import glob
import json
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = {}
NEAR = 0.05


def marginals(n, mu_str):
    mu = {int(s, 2): w for s, w in mu_str.items()}
    tot = sum(mu.values())
    if tot <= 0:
        return None
    m = {a: w / tot for a, w in mu.items()}
    return [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(n)]


def blocks():
    for fname in sorted(glob.glob(str(DATA / "hu_*.json"))):
        try:
            d = json.loads(open(fname).read())
        except Exception:
            continue
        if not isinstance(d, dict):
            continue
        for key, blk in d.items():
            if not isinstance(blk, dict):
                continue
            rows = blk.get("rows")
            if not isinstance(rows, list) or not rows:
                continue
            good = [r for r in rows if isinstance(r, dict)
                    and "mu" in r and "n" in r]
            if good:
                yield os.path.basename(fname), key, good


def part_A():
    print("A. Every stored endpoint, by block:")
    rows = []
    tot = absent = near = 0
    for fname, key, endpoints in blocks():
        a = n_ = 0
        for r in endpoints:
            margs = marginals(r["n"], r["mu"])
            if margs is None:
                continue
            tot += 1
            mm = min(margs)
            if mm <= 1e-12:
                a += 1
            elif mm < NEAR:
                n_ += 1
        absent += a
        near += n_
        rows.append({"source": f"{fname}:{key}", "endpoints": len(endpoints),
                     "absent_coordinate": a, "near_degenerate": n_})
        if a or n_:
            print(f"  {fname}:{key:38s} {len(endpoints):4d} endpoints, "
                  f"{a} with an absent coordinate, {n_} near-degenerate")
    print(f"  TOTAL: {tot} endpoints, {absent} with an absent coordinate "
          f"({100*absent/max(tot,1):.1f}%), {near} near-degenerate")
    OUT["A_by_block"] = {"rows": rows, "total": tot, "absent": absent,
                         "near": near}


def part_B():
    print()
    print("B. The headline floors (minimum-floor endpoint per block):")
    rows = []
    for fname, key, endpoints in blocks():
        cand = [r for r in endpoints
                if "floor" in r or "own_margin" in r]
        if not cand:
            continue
        kf = "own_margin" if "own_margin" in cand[0] else "floor"
        r = min(cand, key=lambda z: z[kf])
        margs = marginals(r["n"], r["mu"])
        if margs is None:
            continue
        eff = sum(1 for x in margs if x > 1e-12)
        rows.append({"source": f"{fname}:{key}", "claimed_n": r["n"],
                     "effective_n": eff, "floor": r[kf],
                     "start": r.get("start", "?"),
                     "overstated": eff < r["n"]})
        mark = "  <-- effectively smaller" if eff < r["n"] else ""
        print(f"  {fname}:{key:38s} claimed n={r['n']}, effective n={eff}, "
              f"floor {r[kf]:+.6f}{mark}")
    bad = [r for r in rows if r["overstated"]]
    print(f"  => {len(bad)} of {len(rows)} headline floors are quoted at a "
          f"higher n than the instance actually uses")
    OUT["B_headlines"] = {"rows": rows, "overstated": len(bad),
                          "total": len(rows)}


def part_C():
    print()
    print("C. What this changes:")
    bad = [r for r in OUT["B_headlines"]["rows"] if r["overstated"]]
    worst = max(bad, key=lambda r: r["claimed_n"] - r["effective_n"])
    print(f"  sharpest overstatement: {worst['source']} -- claimed "
          f"n={worst['claimed_n']}, effective n={worst['effective_n']} "
          f"(start {worst['start']}, floor {worst['floor']:+.6f})")
    print("  positivity and kill claims are UNAFFECTED: every endpoint is "
          "still a valid in-regime instance, and a certified negative CR "
          "at effective n is still a certified negative CR")
    print("  what is affected is dimension labelling: a floor quoted as "
          "'the n = k floor' may be a smaller instance embedded in k, so "
          "it is not evidence about dimension k")
    OUT["C_worst"] = worst


def main():
    part_A()
    part_B()
    part_C()
    (DATA / "hu_effdim.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    print()
    print("checkpoint: data/hu_effdim.json")


if __name__ == "__main__":
    main()
