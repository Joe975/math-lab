#!/usr/bin/env python3
"""Aggregate run_slices.sh output files for the girth>=5 Erdos-Gyarfas search.

Usage: aggregate_slices.py N MOD [--json OUT]

Reads attempts/erdos-gyarfas/data/g5n{N}_s{R}of{MOD}.txt for R = 0..MOD-1,
requires every file to end in a SUMMARY line (i.e. the slice completed),
sums totals / girth histograms / c4-c8-c16 combo histograms, collects all
SURV8 / P8ONLY / CAND lines, and cross-checks the per-slice totals against
the graph counts geng itself reported on stderr (.gengerr files).
"""

import argparse
import json
import os
import re
import sys

DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                    "..", "attempts", "erdos-gyarfas", "data")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("mod", type=int)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()

    total = 0
    girth_hist = {}
    combo = {}
    surv8, p8only, cand = [], [], []
    slice_totals = []
    for r in range(args.mod):
        path = os.path.join(DATA, "g5n%d_s%dof%d.txt" % (args.n, r, args.mod))
        if not os.path.exists(path):
            sys.exit("MISSING slice %d/%d (%s) - run incomplete" % (r, args.mod, path))
        summary = None
        with open(path) as f:
            for line in f:
                line = line.rstrip("\n")
                if line.startswith("SUMMARY "):
                    summary = json.loads(line[8:])
                elif line.startswith("SURV8 "):
                    surv8.append(line)
                elif line.startswith("P8ONLY "):
                    p8only.append(line)
                elif line.startswith("CAND "):
                    cand.append(line)
        if summary is None:
            sys.exit("slice %d/%d has no SUMMARY - incomplete" % (r, args.mod))
        total += summary["total"]
        slice_totals.append(summary["total"])
        for k, v in summary["girth_hist"].items():
            girth_hist[k] = girth_hist.get(k, 0) + v
        for k, v in summary["combo_c4_c8_c16"].items():
            combo[k] = combo.get(k, 0) + v
        # cross-check vs geng's own stderr count
        errpath = path.replace(".txt", ".gengerr")
        if os.path.exists(errpath):
            with open(errpath) as f:
                m = re.search(r">Z\s+(\d+)\s+graphs generated", f.read())
            if m and int(m.group(1)) != summary["total"]:
                sys.exit("slice %d/%d: filter total %d != geng count %s"
                         % (r, args.mod, summary["total"], m.group(1)))

    out = {
        "n": args.n, "mod": args.mod, "total": total,
        "slice_totals": slice_totals,
        "girth_hist": dict(sorted(girth_hist.items(), key=lambda kv: int(kv[0]))),
        "combo_c4_c8_c16": dict(sorted(combo.items())),
        "surv8": surv8, "p8only": p8only, "cand": cand,
    }
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("surv8", "p8only", "cand")}, indent=1))
    print("SURV8 (C8-free) count:", len(surv8))
    for line in surv8:
        print(" ", line)
    print("P8ONLY count:", len(p8only))
    for line in p8only[:50]:
        print(" ", line)
    print("CAND (no C4/C8/C16) count:", len(cand))
    for line in cand:
        print(" ", line)
    if args.json:
        with open(args.json, "w") as f:
            json.dump(out, f, indent=1)
        print("aggregate written to", args.json)


if __name__ == "__main__":
    main()
