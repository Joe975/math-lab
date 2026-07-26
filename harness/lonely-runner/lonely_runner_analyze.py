#!/usr/bin/env python3
"""Verify and analyze a lonely_runner.py scan output JSON.

For each primitive near-tight tuple: checks that every expected scaled copy
(g*v with g*max(v) <= V) was independently hit with identical ML; recomputes
ML three ways (exact sums-only candidate set, exact paranoid larger candidate
set, dense float grid); flags any ML < 1/8 as a counterexample candidate.

Usage: lonely_runner_analyze.py <scan_output.json>
"""
import json, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lonely_runner import ml_exact, ml_exact_paranoid, ml_grid_numeric
from math import gcd
from fractions import Fraction

data = json.load(open(sys.argv[1]))
print(f"V={data['V']} threshold={data['threshold']} scanned={data['n_scanned']} "
      f"raw_hits={data['n_hits_raw']} primitive={data['n_primitive']} time={data['seconds']}s")
prims = {}
for h in data["near_tight_raw"]:
    key = tuple(h["primitive"])
    prims.setdefault(key, []).append(h)

for key in sorted(prims, key=lambda k: (prims[k][0]['ml_float'], k)):
    hs = prims[key]
    mls = {h["ml"] for h in hs}
    scales = sorted(h["gcd"] for h in hs)
    print(f"\nprimitive {key}: ML values across scaled copies {mls}, scales present {scales}")
    # verify scaled-copy completeness: expected scales g with g*max(key) <= V
    expected = list(range(1, data['V'] // max(key) + 1))
    print(f"  expected scales {expected} -> {'COMPLETE' if scales == expected else 'MISSING SOME - BUG?'}")
    # independent re-verification of the primitive
    a, t = ml_exact(key)
    b = ml_exact_paranoid(key)
    g, slack = ml_grid_numeric(key, 200_001)
    status = "OK" if (str(a.numerator)+'/'+str(a.denominator) in mls and a == b and abs(float(a)-g) <= 2*slack) else "MISMATCH"
    print(f"  re-verify: exact={a} paranoid={b} grid={g:.7f}(±{2*slack:.1e}) witness t={t}  {status}")
    below = a < Fraction(1, 8)
    if below:
        print("  *** ML < 1/8: COUNTEREXAMPLE CANDIDATE - UNVERIFIED, ESCALATE ***")
