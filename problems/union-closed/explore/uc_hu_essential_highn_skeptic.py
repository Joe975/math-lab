#!/usr/bin/env python3
"""Skeptic pass on uc_hu_essential_highn.py (attempt 054).  Refute.

Independence: CR is recomputed with the nats history-recursion
evaluator of the 050 skeptic (own code, no shared evaluation path with
uc_hu_order2), and the rollout order is re-derived here from the
records' rules rather than imported.

  W1  every recorded endpoint is genuinely essential (no absent
      coordinate, min marginal >= the stated threshold) -- the whole
      point of the rerun
  W2  every recorded floor reproduces under the independent evaluator
  W3  the comparison claim: the essential floors really are ABOVE the
      degenerate ones recorded in 038/040 for the same n and cap
  W4  a direct attack: descents of my own at n = 6 under the same
      constraint, hunting a violation

Usage: python uc_hu_essential_highn_skeptic.py   (exit 0 iff clean)
"""
from __future__ import annotations

import itertools
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_n3_census_skeptic import cr_of, cstar
from uc_hu_rollcensus import gen_instance

DATA = HERE.parent / "data"
MIN_MARG = 0.03
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def norm(mu):
    t = sum(mu.values())
    return {a: w / t for a, w in mu.items() if w / t > 1e-14}


def margs_of(n, m):
    return [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(n)]


def cond_entropy(n, m, i, S):
    d = {}
    for a, w in m.items():
        k = tuple((a >> j) & 1 for j in S)
        e = d.setdefault(k, [0.0, 0.0])
        e[0] += w
        if not (a >> i) & 1:
            e[1] += w
    tot = 0.0
    for mm, z in d.values():
        if mm > 0:
            p = z / mm
            tot += mm * (0.0 if p <= 0 or p >= 1 else
                         -(p * math.log(p) + (1 - p) * math.log(1 - p))
                         / math.log(2))
    return tot


def canon_completion(n, m, chosen):
    seq = list(chosen)
    rem = sorted(set(range(n)) - set(chosen))
    while rem:
        sc = sorted((cond_entropy(n, m, i, seq), i) for i in rem)
        lo = sc[0][0]
        pick = min(i for v, i in sc if v <= lo + 1e-12)
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq)


def roll_seq(n, m):
    """Own rollout: maximise full CR with canonical completion."""
    seq = []
    rem = sorted(range(n))
    while rem:
        best = None
        for i in rem:
            full = canon_completion(n, m, seq + [i])
            cr, _ = cr_of(n, m, full)
            if best is None or cr > best[0] + 1e-12:
                best = (cr, i)
        seq.append(best[1])
        rem.remove(best[1])
    return tuple(seq)


def roll_margin(n, mu):
    m = norm(mu)
    mg = margs_of(n, m)
    if min(mg) < MIN_MARG or max(mg) >= 0.5:
        return None
    cr, H = cr_of(n, m, roll_seq(n, m))
    if H < 0.2:
        return None
    return cr / H - cstar(max(mg))


def main():
    o = json.loads((DATA / "hu_essential_highn.json").read_text())

    print("W1/W2. Endpoints: essential, and floors reproduce:")
    for key, blk in o.items():
        rows = blk.get("rows", [])
        if not rows:
            continue
        n = rows[0]["n"]
        noness, bad = [], []
        for r in rows:
            m = norm({int(s, 2): w for s, w in r["mu"].items()})
            mg = margs_of(n, m)
            if min(mg) <= 1e-12 or min(mg) < MIN_MARG - 1e-9:
                noness.append((r["start"], min(mg)))
            if key.startswith("C_"):
                continue                      # best-order: checked in W2b
            v = roll_margin(n, {int(s, 2): w for s, w in r["mu"].items()})
            if v is None or abs(v - r["floor"]) > 1e-6:
                bad.append((r["start"], v, r["floor"]))
        check(f"{key}: every endpoint keeps all {n} coordinates",
              not noness, f"{noness[:2]}")
        if not key.startswith("C_"):
            check(f"{key}: floors reproduce under the independent "
                  "evaluator", not bad, f"{bad[:2]}")

    print("W2b. The best-order block, by full 720-order enumeration:")
    blk = o.get("C_best_n6_cap_0.49")
    if blk and blk["rows"]:
        bad = []
        for r in blk["rows"]:
            m = norm({int(s, 2): w for s, w in r["mu"].items()})
            mg = margs_of(6, m)
            best, H = None, None
            for s in itertools.permutations(range(6)):
                cr, H = cr_of(6, m, s)
                if best is None or cr > best:
                    best = cr
            v = best / H - cstar(max(mg))
            if abs(v - r["floor"]) > 1e-6:
                bad.append((r["start"], v, r["floor"]))
        check("best-order endpoints reproduce", not bad, f"{bad[:2]}")

    print("W3. Essential floors versus the degenerate ones on record:")
    old = json.loads((DATA / "hu_rollcensus.json").read_text())
    pairs = [("A_n6_cap_0.49", old["P2_cap_0.49"]["global_floor"]),
             ("A_n6_cap_0.497", old["P2_cap_0.497"]["global_floor"]),
             ("B_n7_cap_0.49", old["P3_cap_0.49"]["global_floor"])]
    for key, oldfloor in pairs:
        if key not in o:
            continue
        new = o[key]["floor"]
        check(f"{key}: essential floor is above the recorded one",
              new > oldfloor,
              f"essential {new:+.6f} vs recorded {oldfloor:+.6f}")

    print("W4. My own constrained descents at n = 6, hunting a violation:")
    rng = random.Random(9911)
    worst = None
    for _ in range(4):
        mu = gen_instance(rng, 6, 0.49)
        cur = roll_margin(6, mu)
        if cur is None:
            continue
        step = 0.35
        while step > 0.05:
            improved = False
            for a in sorted(mu):
                for f in (1 + step, 1 / (1 + step)):
                    cand = dict(mu)
                    cand[a] = cand[a] * f
                    v = roll_margin(6, cand)
                    if v is not None and v < cur - 1e-12:
                        mu, cur, improved = cand, v, True
            if not improved:
                step /= 2
        if worst is None or cur < worst:
            worst = cur
    check("no violation found under my own constrained descents",
          worst is None or worst >= -1e-12,
          f"min margin {worst:+.6f}" if worst is not None else "no runs")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
