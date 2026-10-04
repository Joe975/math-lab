#!/usr/bin/env python3
"""Attempt 061 (060 lead 1): one check that would have caught all three
leaked constraints.

Three times this route has recorded endpoints that violate the region
their own campaign claimed to search:

  dimension     (052)  11 of 27 headline floors quoted at a higher n
                       than the instance uses
  essentiality  (051)  descents saturating a bound by escaping to a
                       lower-dimensional face
  the cap       (060)  an "equality at cap 0.499" endpoint whose own
                       max marginal is 0.0360

Each was caught by inspecting endpoints by hand, after the fact.  This
file is the check that finds all three at once, run over every
committed checkpoint: for each stored endpoint, does it satisfy the
region its block name declares?

Rules applied per endpoint:
  R1  no coordinate absent (marginal exactly 0) -- else the effective
      dimension is below the label
  R2  every coordinate marginal >= 0.03 -- the 051 essentiality bar
  R3  max marginal < 1/2 -- in regime at all
  R4  if the block name declares a cap, the endpoint's own max
      marginal is within 0.02 of it -- the 060 bar (loose, so it flags
      only real drift)

Blocks predating a rule are still reported: the point is a map of
which recorded numbers describe the region they claim, not a verdict
on records written before the rule existed.

Standard library only; deterministic (pure audit).
Usage: python uc_hu_region_audit.py
Checkpoint: ../data/hu_region_audit.json
"""
from __future__ import annotations

import glob
import json
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
OUT = {}
CAP_RE = re.compile(r"cap[_-]?(0\.\d+)")


def marginals(n, mu_str):
    mu = {int(s, 2): w for s, w in mu_str.items()}
    tot = sum(mu.values())
    if tot <= 0:
        return None
    m = {a: w / tot for a, w in mu.items()}
    return [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(n)]


def audit_endpoint(n, margs, cap):
    flags = []
    if min(margs) <= 1e-12:
        flags.append("R1_absent_coordinate")
    elif min(margs) < 0.03:
        flags.append("R2_below_essentiality")
    if max(margs) >= 0.5:
        flags.append("R3_out_of_regime")
    if cap is not None and max(margs) < cap - 0.02:
        flags.append("R4_drifted_off_cap")
    return flags


def main():
    rows = []
    tot = flagged = 0
    counts = {}
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
            eps = blk.get("rows")
            if not isinstance(eps, list) or not eps:
                continue
            m = CAP_RE.search(key)
            cap = float(m.group(1)) if m else None
            bad = {}
            n_ok = 0
            for r in eps:
                if not isinstance(r, dict) or "mu" not in r or "n" not in r:
                    continue
                margs = marginals(r["n"], r["mu"])
                if margs is None:
                    continue
                tot += 1
                fl = audit_endpoint(r["n"], margs, cap)
                if fl:
                    flagged += 1
                    for f in fl:
                        bad[f] = bad.get(f, 0) + 1
                        counts[f] = counts.get(f, 0) + 1
                else:
                    n_ok += 1
            if bad:
                src = f"{os.path.basename(fname)}:{key}"
                rows.append({"source": src, "cap": cap, "clean": n_ok,
                             "flags": bad})
                print(f"  {src:44s} clean {n_ok:3d} | " +
                      ", ".join(f"{k.split('_')[0]}×{v}"
                                for k, v in sorted(bad.items())))
    print()
    print(f"  {tot} endpoints audited, {flagged} flagged "
          f"({100*flagged/max(tot,1):.1f}%)")
    for k in sorted(counts):
        print(f"    {k}: {counts[k]}")
    OUT.update({"endpoints": tot, "flagged": flagged, "by_rule": counts,
                "blocks": rows})
    (DATA / "hu_region_audit.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    print()
    print("checkpoint: data/hu_region_audit.json")


if __name__ == "__main__":
    main()
