#!/usr/bin/env python3
"""Classify census endpoints and probe whether sub-2 endpoints are genuine
local maxima or optimizer stalls.

Escape probe: at each endpoint with final ratio < 1.97, try NP random
perturbation directions at several scales; record the best ratio found.  An
endpoint from which random perturbation finds ascent is a STALL, not a local
maximum, and is excluded from any 'local maxima all look like X' statement.

Float-only exploration tool; claims are certified separately.
"""

import json
import math
import random
import sys

sys.path.insert(0, "/home/user/work-crouzeix/explore")
from census import objective, unpack, diagnostics, nelder_mead

SCALES = (1e-3, 3e-3, 1e-2, 3e-2)
NPROBE = 60          # directions per scale


def escape_probe(x, v, rng, ntheta=96):
    """Random perturbations at several scales; if none ascend, one full
    Nelder-Mead restart (a fresh simplex sees ridge directions random
    probes miss)."""
    best = v
    n = len(x)
    for s in SCALES:
        for _ in range(NPROBE):
            y = [xi + rng.gauss(0, s) for xi in x]
            fy = objective(y, ntheta)
            if fy > best:
                best = fy
    if best <= v + 5e-4:
        _, vnm = nelder_mead(lambda z: objective(z, 96), x, step=0.08,
                             maxit=600)
        best = max(best, vnm)
    return best


def classify(rec, escaped):
    d = rec["diag"]
    r = d["ratio"]
    if escaped:
        return "stall"
    if r > 1.97:
        # near-2: split by eigenvalue structure
        g = d["eig_gaps"]        # sorted pairwise |ev_i - ev_j| / diam W
        if g[2] < 0.15:
            return "near2-triple"      # all eigenvalues coincide: J3-like
        if g[0] < 0.15:
            return "near2-double"      # two coincide: J2 (+) [c]-like
        return "near2-other"
    if abs(r - 1.0) < 1e-6:
        return "one-exact"
    if r < 1.97:
        return "below2-max"            # survived the probe below 1.97
    return "other"


def main():
    files = sys.argv[1:]
    recs = []
    for fn in files:
        with open(fn) as f:
            for line in f:
                recs.append(json.loads(line))
    print(f"{len(recs)} endpoints")
    rng = random.Random(12345)
    out = []
    for i, rec in enumerate(recs):
        v = rec["diag"]["ratio"]
        escaped = False
        best = v
        if v < 1.97:
            best = escape_probe(rec["x"], objective(rec["x"], 96), rng)
            escaped = best > v + 5e-4
        cls = classify(rec, escaped)
        rec["class"] = cls
        rec["probe_best"] = best
        out.append(rec)
        if (i + 1) % 50 == 0:
            print(f"  probed {i+1}/{len(recs)}", file=sys.stderr)
    with open("/home/user/work-crouzeix/data/classified.jsonl", "w") as f:
        for rec in out:
            f.write(json.dumps(rec) + "\n")

    # summary
    from collections import Counter, defaultdict
    byc = Counter(r["class"] for r in out)
    print("\nclass counts:", dict(byc))
    byfc = defaultdict(Counter)
    for r in out:
        byfc[r["family"]][r["class"]] += 1
    for fam in sorted(byfc):
        print(f"  {fam:11s} {dict(byfc[fam])}")
    print("\nratio distribution of non-stall endpoints:")
    hist = Counter()
    for r in out:
        if r["class"] != "stall":
            hist[round(r["diag"]["ratio"], 1)] += 1
    for k in sorted(hist):
        print(f"  {k:4.1f}: {hist[k]}")
    # the interesting ones: survived probing, not near 2, not exactly 1
    odd = [r for r in out if r["class"] == "below2-max"]
    odd.sort(key=lambda r: -r["diag"]["ratio"])
    print(f"\n{len(odd)} probe-surviving endpoints strictly between 1 and 1.97:")
    for r in odd[:20]:
        d = r["diag"]
        print(f"  {r['family']:11s} seed={r['seed']} ratio={d['ratio']:.5f} "
              f"gaps={[round(g,3) for g in d['eig_gaps']]} "
              f"disc={d['discness']:.3f} peaks={d['peaks']} "
              f"nonnormal={d['nonnormal']:.3f}")


if __name__ == "__main__":
    main()
