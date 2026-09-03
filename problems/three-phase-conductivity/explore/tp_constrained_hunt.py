#!/usr/bin/env python3
"""Hunt for a structure that RESPECTS Cherkaev's constraint (4.26) and still
sits below his B2 bound, below m11.

Attempt 009 found that B2 separates exactly on (4.26) within one laminate
family: every member below B2 violates the constraint, every member respecting
it lies above B2.  A single counterexample here -- a (4.26)-respecting structure
under B2 -- would overturn that reading and reopen the escalation of 006.  So
this is the falsification test for 009, run over MANY topologies rather than
one family.

Search space: all binary laminate trees of a given rank with axis normals and
leaves labelled by the three phases, node fractions optimised by Nelder-Mead
against isotropy and the target volume fractions (a float SCREEN), then the
survivors re-evaluated exactly and classified by (4.26).

Usage:
    python tp_constrained_hunt.py --selftest
    python tp_constrained_hunt.py --m1 3/25 --rank 4 [--restarts 2] [--out X.json]
Standard library only.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import random
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "three-phase-conductivity"))
sys.path.insert(0, HERE)
import laminate as L  # noqa: E402
import tp_fields as TF  # noqa: E402
import tp_cherkaev_bound as CB  # noqa: E402
import tp_search as TS  # noqa: E402

SIGMA = (Fr(1), Fr(2), Fr(5))
PHASES = {"p1": SIGMA[0], "p2": SIGMA[1], "p3": SIGMA[2]}
M2 = Fr(1, 4)
AXES = [(1, 0), (0, 1)]


def slack_4_26(tree, phases=PHASES):
    """min over phases 1,2 of (S - varsigma_N)^2 - D^2, exact. Negative = violates.

    S and D are the rotationally invariant components of the field matrix
    Z = grad u built from the two orthogonal applied fields (Cherkaev section 2.3).
    """
    a = TF.leaf_fields(tree, phases, (Fr(1), Fr(0)))
    b = TF.leaf_fields(tree, phases, (Fr(0), Fr(1)))
    per = {}
    for (na, va, Ea), (nb, vb, Eb) in zip(a, b):
        Z11, Z12, Z21, Z22 = Ea[0], Ea[1], Eb[0], Eb[1]
        per.setdefault(na, []).append((va, Z11 + Z22, Z11 - Z22, Z12 + Z21))
    if "p3" not in per:
        return None
    tot = sum(v for v, _, _, _ in per["p3"])
    sN = sum(v * S for v, S, _, _ in per["p3"]) / tot
    vals = [(S - sN) ** 2 - (Ds * Ds + Dss * Dss)
            for nm in ("p1", "p2") if nm in per for _, S, Ds, Dss in per[nm]]
    return min(vals) if vals else None


def value_of(tree):
    """Effective value as the quadratic form, exact; None if degenerate."""
    e = L.effective(tree, PHASES)
    return (e[0] + e[2]) / 2, e[0] - e[2], e[1]


def hunt(m1: Fr, rank: int, restarts: int, seed: int, tol: float):
    """Screen every axis-normal topology; return exact records for survivors."""
    f = [float(m1), float(M2), float(1 - m1 - M2)]
    sig = [float(x) for x in SIGMA]
    rng = random.Random(seed)
    b2 = CB.B2(SIGMA, m1, Fr(1, 2))
    hs = CB.hs_lo(SIGMA, m1, M2)
    out = []
    required = {0, 1, 2}
    for shape in TS.shapes(rank):
        for labels in itertools.product(range(3), repeat=rank + 1):
            if required - set(labels):
                continue
            for norms in itertools.product(AXES, repeat=rank):
                top = TS.Topology(shape, labels, norms)

                def fn(x, top=top):
                    return TS.objective(top, x, sig, f, 1.0)
                for _ in range(restarts):
                    x0 = [rng.uniform(-2, 2) for _ in range(top.n_params)]
                    x, _ = TS.nelder_mead(fn, x0)
                    x, _ = TS.nelder_mead(fn, x, step=0.05, iters=400)
                    _, mean, aniso, ferr = fn(x)
                    if aniso > tol or ferr > tol:
                        continue
                    params = [TS.sigmoid(v) for v in x]
                    rats = [Fr(p).limit_denominator(10 ** 7) for p in params]
                    if any(not (0 < r < 1) for r in rats):
                        continue
                    tree = top.to_tree([str(r) for r in rats])
                    try:
                        val, dgap, off = value_of(tree)
                        sl = slack_4_26(tree)
                    except (ZeroDivisionError, AssertionError, KeyError):
                        continue
                    if sl is None:
                        continue
                    out.append({
                        "value": val, "slack": sl, "aniso_exact": dgap,
                        "offdiag": off, "labels": labels, "normals": norms,
                        "shape": shape, "params": [str(r) for r in rats],
                    })
    return out, hs, b2


def report(m1: Fr, rank: int, restarts: int, seed: int, tol: float, out_path=None):
    rows, hs, b2 = hunt(m1, rank, restarts, seed, tol)
    below = [r for r in rows if r["value"] < b2]
    resp = [r for r in rows if r["slack"] >= 0]
    both = [r for r in rows if r["value"] < b2 and r["slack"] >= 0]
    print(f"m1={m1} rank={rank}  HS_lo={float(hs):.10f}  B2={float(b2):.10f}")
    print(f"  screened survivors: {len(rows)}")
    print(f"  below B2:                       {len(below)}")
    print(f"  respecting (4.26):              {len(resp)}")
    print(f"  BOTH below B2 and respecting:   {len(both)}   <-- would REOPEN 006")
    if resp:
        m = min(resp, key=lambda r: r["value"])
        print(f"  min value among respecting: {float(m['value']):.10f}"
              f"  (B2 {float(b2):.10f}, margin {float(m['value']-b2):+.3e})")
    if below:
        m = min(below, key=lambda r: r["value"])
        print(f"  min value below B2:         {float(m['value']):.10f}"
              f"  slack {float(m['slack']):+.3e}")
    for r in both[:5]:
        print(f"  COUNTEREXAMPLE: value={float(r['value']):.10f} slack={float(r['slack']):+.3e}"
              f" labels={r['labels']} normals={r['normals']}")
    if out_path:
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump({"m1": str(m1), "rank": rank, "hs": str(hs), "b2": str(b2),
                       "n_survivors": len(rows), "n_below": len(below),
                       "n_respecting": len(resp), "n_both": len(both),
                       "counterexamples": [
                           {k: (str(v) if isinstance(v, Fr) else v) for k, v in r.items()}
                           for r in both[:20]]}, fh, indent=1)
        print(f"  written to {out_path}")
    return both


def selftest() -> None:
    """The known below-B2 structure must be detected as violating (4.26), and
    the known attaining structure as making it exactly active."""
    import tp_attain as A
    _, tree, _ = A.attains(SIGMA, Fr(1, 2), Fr(4, 5))
    assert slack_4_26(tree) == 0, "attaining structure must make (4.26) active"
    path = os.path.join(HERE, "..", "data", "attained", "below-m11-m1-3_25.json")
    rec = json.load(open(path, encoding="utf-8"))
    ph = {k: Fr(v) for k, v in rec["phases"].items()}
    assert slack_4_26(rec["tree"], ph) < 0, "below-B2 structure must violate (4.26)"
    print("selftest: (4.26) slack is 0 at the attaining structure and negative "
          "at the below-B2 structure")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--m1", default="3/25")
    ap.add_argument("--rank", type=int, default=4)
    ap.add_argument("--restarts", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--tol", type=float, default=1e-7)
    ap.add_argument("--out")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    else:
        report(Fr(a.m1), a.rank, a.restarts, a.seed, a.tol, a.out)
