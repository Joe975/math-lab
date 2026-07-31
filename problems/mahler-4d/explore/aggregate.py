#!/usr/bin/env python3
"""Aggregate census eval slices into a distribution summary + candidate list.

Usage: aggregate.py KMIN KMAX [--slice 200000] [--near 10.75]
                    [--datadir data] [--json data/summary.json]

Verifies completeness first: for each level k, every expected slice file
data/eval{k}_{start}.txt must exist with the expected line count, and the
masks seen must be exactly the level file's masks (order preserved).
Then reports, per level and overall: status counts, min/max P over `st=ok`
bodies, a histogram of P, all bodies with P < NEAR (written to
data/candidates_k{K}.txt), and every `punt` mask (data/punts.txt) -- punts
carry no float value and MUST be exact-checked separately.
"""
import argparse
import json
import math
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kmin", type=int)
    ap.add_argument("kmax", type=int)
    ap.add_argument("--slice", type=int, default=200000)
    ap.add_argument("--near", type=float, default=10.75)
    ap.add_argument("--datadir", default="data")
    ap.add_argument("--json", default="data/summary.json")
    args = ap.parse_args()

    overall = {"levels": {}, "near_threshold": args.near}
    all_candidates = []
    all_punts = []
    grand_min = (math.inf, None)

    for k in range(args.kmin, args.kmax + 1):
        level_path = os.path.join(args.datadir, f"level{k:02d}.masks")
        with open(level_path) as f:
            masks = [ln.split()[0] for ln in f if ln.strip()]
        n = len(masks)
        stats = {"orbits": n, "ok": 0, "improper": 0, "degen": 0, "punt": 0}
        hist = {}
        kmin_p = (math.inf, None)
        kmax_p = (-math.inf, None)
        cands = []
        pos = 0
        start = 0
        while start < n:
            want = min(args.slice, n - start)
            path = os.path.join(args.datadir, f"eval{k:02d}_{start:08d}.txt")
            if not os.path.exists(path):
                sys.exit(f"MISSING {path} -- eval incomplete")
            with open(path) as f:
                lines = f.read().splitlines()
            if len(lines) != want:
                sys.exit(f"{path}: {len(lines)} lines, expected {want}")
            for ln in lines:
                fields = ln.split()
                if fields[0] != masks[pos]:
                    sys.exit(f"{path}: mask {fields[0]} != level file "
                             f"{masks[pos]} at index {pos}")
                pos += 1
                kv = dict(x.split("=", 1) for x in fields[1:])
                st = kv["st"]
                stats[st] += 1
                if st == "punt":
                    all_punts.append((k, fields[0]))
                if st != "ok":
                    continue
                p = float(kv["P"])
                if p < kmin_p[0]:
                    kmin_p = (p, fields[0])
                if p > kmax_p[0]:
                    kmax_p = (p, fields[0])
                b = f"{math.floor(p * 4) / 4:.2f}"   # bins of width 0.25
                hist[b] = hist.get(b, 0) + 1
                if p < args.near:
                    cands.append((p, fields[0], ln))
            start += args.slice
        if pos != n:
            sys.exit(f"level {k}: consumed {pos} lines, expected {n}")
        cands.sort()
        cpath = os.path.join(args.datadir, f"candidates_k{k:02d}.txt")
        with open(cpath, "w") as f:
            for p, m, ln in cands:
                f.write(ln + "\n")
        all_candidates += [(p, k, m) for p, m, ln in cands]
        stats["min_P"] = kmin_p
        stats["max_P"] = kmax_p
        stats["near_count"] = len(cands)
        stats["hist"] = dict(sorted(hist.items(), key=lambda kv: float(kv[0])))
        overall["levels"][k] = stats
        if kmin_p[0] < grand_min[0]:
            grand_min = (kmin_p[0], (k, kmin_p[1]))
        print(f"level {k:2d}: {stats['orbits']:8d} orbits  ok={stats['ok']:8d} "
              f"improper={stats['improper']:8d} degen={stats['degen']:6d} "
              f"punt={stats['punt']:4d}  minP={kmin_p[0]:.9f} ({kmin_p[1]})  "
              f"near<{args.near}: {len(cands)}")

    ppath = os.path.join(args.datadir, "punts.txt")
    with open(ppath, "w") as f:
        for k, m in all_punts:
            f.write(f"{m} k={k}\n")
    overall["grand_min"] = grand_min
    overall["punt_total"] = len(all_punts)
    all_candidates.sort()
    overall["near_total"] = len(all_candidates)
    print(f"\ngrand minimum float P = {grand_min[0]:.12f} at level "
          f"{grand_min[1][0]} mask {grand_min[1][1]}")
    print(f"32/3 = {32/3:.12f}")
    print(f"bodies with P < {args.near}: {len(all_candidates)}; punts: "
          f"{len(all_punts)} (exact-check both lists)")
    with open(args.json, "w") as f:
        json.dump(overall, f, indent=1)
    print(f"summary written to {args.json}")


if __name__ == "__main__":
    main()
