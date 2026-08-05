#!/usr/bin/env python3
"""Aggressive escape attempt on selected endpoints: high-resolution objective
(512 boundary samples), multiple NM restarts with jitter, 1200 random probes.
An endpoint is a STALL if anything ascends by > 1e-3."""
import json, random, sys
sys.path.insert(0, "/home/user/work-crouzeix/explore")
from census import objective, nelder_mead, compass

def deep(x, v0, rng):
    f = lambda z: objective(z, 512)
    v = f(x)
    best = v
    for s in (1e-3, 3e-3, 1e-2, 3e-2, 1e-1):
        for _ in range(240):
            y = [xi + rng.gauss(0, s) for xi in x]
            best = max(best, objective(y, 128))
    for k in range(3):
        x0 = [xi + rng.gauss(0, 0.02) for xi in x]
        _, vnm = nelder_mead(lambda z: objective(z, 128), x0, step=0.1, maxit=1200)
        best = max(best, vnm)
    xc, vc = compass(f, x, v, step=0.005, maxit=1200)
    best = max(best, vc)
    return v, best

rng = random.Random(999)
targets = [json.loads(l) for l in open("/home/user/work-crouzeix/data/classified.jsonl")
           if json.loads(l)["class"] == "below2-max"]
for r in targets:
    v512, best = deep(r["x"], r["diag"]["ratio"], rng)
    verdict = "STALL" if best > v512 + 1e-3 else "SURVIVES"
    print(f"{r['family']} seed={r['seed']} ratio96={r['diag']['ratio']:.5f} "
          f"ratio512={v512:.5f} best_escape={best:.5f} {verdict}", flush=True)
