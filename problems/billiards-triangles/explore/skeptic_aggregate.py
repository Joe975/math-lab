#!/usr/bin/env python3
"""SKEPTIC re-aggregation of the coverage table from the raw census data.

Reads data/cover/map_g32.jsonl and data/cover/anchor_g32.jsonl and recomputes
every headline aggregate with independent bookkeeping:
  * interior/boundary cell counts and areas,
  * area-weighted candidate and anchored fractions by max word length,
  * gamma-band candidate fractions,
  * total exactly-certified area over ALL anchor sources (map + arcs +
    families).
Optionally spot-checks candidate lengths at random cell midpoints with the
skeptic float corridor (own implementation).
"""

import glob
import json
import math
import random
import sys
from fractions import Fraction as F

sys.path.insert(0, "problems/billiards-triangles/explore")
from skeptic_family import width, gamma_deg  # noqa: E402

D = "problems/billiards-triangles/data"


def load_jsonl(path):
    return [json.loads(l) for l in open(path) if l.strip()]


def main():
    cells = load_jsonl(f"{D}/cover/map_g32.jsonl")
    anchors = {(r["i"], r["j"]): r for r in load_jsonl(f"{D}/cover/anchor_g32.jsonl")}

    interior = [r for r in cells if r["status"] != "BOUNDARY"]
    boundary = [r for r in cells if r["status"] == "BOUNDARY"]
    cand = [r for r in interior if r["status"] == "CAND"]
    nocand = [r for r in interior if r["status"] == "NOCAND"]

    def area(rec):
        x0, x1, y0, y1 = (F(v) for v in rec["box"])
        return float((x1 - x0) * (y1 - y0))

    int_area = sum(area(r) for r in interior)
    region_area = math.pi / 16
    print(f"cells: interior {len(interior)} boundary {len(boundary)} "
          f"CAND {len(cand)} NOCAND {len(nocand)} "
          f"anchored {sum(1 for a in anchors.values() if a.get('status')=='ANCHORED')}")
    print(f"interior area {int_area:.4f}  boundary fringe "
          f"{region_area - int_area:.4f}  (region {region_area:.4f})")

    print("\nmax_len  cand%   anchored%")
    for L in (14, 16, 18, 20, 22, 24, 26):
        ca = sum(area(r) for r in cand if r["cand_len"] <= L)
        aa = sum(area(r) for r in interior
                 if anchors.get((r["i"], r["j"]), {}).get("status") == "ANCHORED"
                 and anchors[(r["i"], r["j"])]["len"] <= L)
        print(f"  {L:2d}   {100*ca/int_area:5.1f}   {100*aa/int_area:5.1f}")

    ca26 = sum(area(r) for r in cand)
    print(f"\nno candidate <= 26 at midpoint: {100*(1-ca26/int_area):.1f}% of interior")

    print("\ngamma band   cand%")
    for glo, ghi in ((92, 95), (95, 100), (100, 105), (105, 112.3),
                     (112.3, 120), (120, 135), (135, 180)):
        tot = sum(area(r) for r in interior if glo <= r["gamma_mid"] < ghi)
        cb = sum(area(r) for r in cand if glo <= r["gamma_mid"] < ghi)
        print(f"  {glo}-{ghi}: {100*cb/tot if tot else float('nan'):5.1f}"
              f"  (area {tot:.4f})")

    # total certified area across every anchor source
    tot_cert = 0.0
    n = 0
    for a in anchors.values():
        if a.get("status") == "ANCHORED":
            tot_cert += float((2 * F(a["half_width"])) ** 2)
            n += 1
    for f in glob.glob(f"{D}/growth/arcsweep_*.jsonl"):
        for r in load_jsonl(f):
            if r.get("anchor"):
                tot_cert += float((2 * F(r["anchor"]["half_width"])) ** 2)
                n += 1
    for f in glob.glob(f"{D}/growth/family_W*.json"):
        r = json.load(open(f))
        if r.get("anchor"):
            tot_cert += float((2 * F(r["anchor"]["half_width"])) ** 2)
            n += 1
    print(f"\ncertified area over {n} anchors: {tot_cert:.3g}")

    # spot-check candidate lengths at random midpoints, full universe
    if "--spot" in sys.argv:
        words = {}
        for L in range(6, 27, 2):
            try:
                words[L] = [tuple(int(c) for c in l.strip())
                            for l in open(f"{D}/words/words_L{L:02d}.jsonl")]
            except FileNotFoundError:
                pass
        rng = random.Random(20260730)
        sample = rng.sample(interior, 8)
        print("\nspot-check of cand_len at cell midpoints (own float corridor):")
        bad = 0
        for r in sample:
            mx, my = float(F(r["mid"][0])), float(F(r["mid"][1]))
            mine = None
            for L in sorted(words):
                if any(width(w, mx, my) > 0 for w in words[L]):
                    mine = L
                    break
            ok = mine == r.get("cand_len")
            bad += not ok
            print(f"  cell ({r['i']},{r['j']}) mid {r['mid']}: "
                  f"mine {mine} theirs {r.get('cand_len')} "
                  f"{'OK' if ok else 'MISMATCH'}")
        print("spot-check:", "all agree" if not bad else f"{bad} MISMATCHES")


if __name__ == "__main__":
    main()
