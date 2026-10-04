#!/usr/bin/env python3
"""Skeptic pass on uc_hu_n3_census.py (attempt 050).  Stance: refute.

Independence: the HU coupling, CR and c* are rebuilt here in NATS by an
explicit history recursion that shares no code with uc_hu_order2 or the
046-048 stack; the census instances are regenerated from the recorded
seeds (the generator is the shared part -- what is independently
recomputed is every VALUE).

  U1  regenerate each census block and recompute best/rollout/identity
      own-constant margins with the independent evaluator; the recorded
      minima and violation counts must reproduce
  U2  attack the n = 3 claim directly: a descent that minimises the
      best-order own-constant margin at n = 3, from many starts
  U3  the part-D descent endpoints: re-evaluate each by full order
      enumeration with the independent evaluator
  U4  sanity: at n = 3 the identity order is not always the best order
      (otherwise "the quantifier is vacuous" would be trivial rather
      than a finding)

Usage: python uc_hu_n3_census_skeptic.py    (exit 0 iff nothing refuted)
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
from uc_hu_rollcensus import gen_instance          # generator only

DATA = HERE.parent / "data"
LN2 = math.log(2.0)
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hb(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1 - p) * math.log(1 - p)) / LN2


def cstar(p):
    return (hb(max(0.5, 1 - 2 * p)) - hb(p)) / hb(p)


def cr_of(n, mu, seq):
    """CR by explicit history recursion (own code, nats-based h)."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items() if w / tot > 1e-14}
    H = -sum(w * math.log(w) for w in mu.values()) / LN2
    cells = [((), (), 1.0)]
    S = 0.0
    for k in range(n):
        i = seq[k]
        nxt = []
        for pa, pb, w in cells:
            def cond(pref):
                sel = [(a, wt) for a, wt in mu.items()
                       if all(((a >> seq[j]) & 1) == v
                              for j, v in enumerate(pref))]
                m = sum(wt for _, wt in sel)
                return m, (sum(wt for a, wt in sel
                               if not (a >> i) & 1) / m if m else 0.0)
            ma, x = cond(pa)
            mb, y = cond(pb)
            z = min(max(0.5, x + y - 1.0), x, y)
            S += w * hb(z)
            for va, vb, cw in ((0, 0, z), (0, 1, x - z),
                               (1, 0, y - z), (1, 1, 1 - x - y + z)):
                if cw > 1e-15:
                    nxt.append((pa + (va,), pb + (vb,), w * cw))
        cells = nxt
    return S - H, H


def margins(n, mu):
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
    if not m:
        return None
    fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
               for i in range(n))
    if not (0 < fmax < 0.5):
        return None
    vals = {}
    H = None
    for s in itertools.permutations(range(n)):
        cr, H = cr_of(n, m, s)
        vals[s] = cr
    if H is None or H < 0.2:
        return None
    c = cstar(fmax)
    return {"best": max(vals.values()) / H - c,
            "ident": vals[tuple(range(n))] / H - c,
            "H": H, "fmax": fmax,
            "best_is_ident": abs(max(vals.values())
                                 - vals[tuple(range(n))]) < 1e-12}


def regen(n, cap, count, seed):
    rng = random.Random(seed)
    out = []
    tried = 0
    while len(out) < count and tried < 40 * count:
        tried += 1
        r = margins(n, gen_instance(rng, n, cap))
        if r is not None:
            out.append(r)
    return out


def main():
    o = json.loads((DATA / "hu_n3_census.json").read_text())

    print("U1. Census blocks recomputed with the independent evaluator:")
    for key, n, seedbase in (("A_n3", 3, 5000), ("B_n4", 4, 5100)):
        for capkey, rec in o[key].items():
            cap = float(capkey.split("_")[1])
            rows = regen(n, cap, rec["instances"], seedbase + int(cap * 1000))
            mb = min(r["best"] for r in rows)
            mi = min(r["ident"] for r in rows)
            vb = sum(1 for r in rows if r["best"] < -1e-12)
            vi = sum(1 for r in rows if r["ident"] < -1e-12)
            ok = (abs(mb - rec["min_best"]) < 1e-8
                  and abs(mi - rec["min_ident"]) < 1e-8
                  and vb == rec["violations"]["m_best"]
                  and vi == rec["violations"]["m_ident"])
            check(f"{key} {capkey}: minima and violation counts", ok,
                  f"best {mb:+.6f} vs {rec['min_best']:+.6f}, "
                  f"identity {mi:+.6f} vs {rec['min_ident']:+.6f}, "
                  f"violations {vb}/{vi}")

    print("U2. n = 3 attacked directly (descent on the best-order margin):")
    rng = random.Random(6001)
    best = None
    for _ in range(60):
        mu = gen_instance(rng, 3, 0.49)
        r = margins(3, mu)
        if r is None:
            continue
        cur = r["best"]
        keys = sorted(mu)
        step = 0.3
        while step > 1e-4:
            improved = False
            for a in keys:
                for f in (1 + step, 1 / (1 + step)):
                    cand = dict(mu)
                    cand[a] = cand[a] * f
                    rr = margins(3, cand)
                    if rr and rr["best"] < cur - 1e-12:
                        mu, cur = cand, rr["best"]
                        improved = True
            if not improved:
                step /= 2
        if best is None or cur < best:
            best = cur
    check("best-order own-constant margin stays >= 0 at n = 3 under "
          "descent", best >= -1e-12, f"min {best:+.6e} over 60 descents")

    print("U3. Part-D endpoints re-evaluated:")
    bad = []
    for capkey in ("D_cap_0.45", "D_cap_0.49"):
        for r in o[capkey]["rows"]:
            mu = {int(s, 2): w for s, w in r["mu"].items()}
            m = margins(3, mu)
            if m is None or abs(m["best"] - r["own_margin"]) > 1e-7:
                bad.append((capkey, r["start"],
                            None if m is None else m["best"],
                            r["own_margin"]))
    check("every part-D endpoint reproduces", not bad, f"{bad[:2]}")

    print("U4. Is the identity order trivially the best at n = 3?")
    rows = regen(3, 0.49, 300, 7007)
    same = sum(1 for r in rows if r["best_is_ident"])
    check("identity is NOT always the best order (so the vacuity "
          "finding has content)", same < len(rows),
          f"identity = best on {same}/{len(rows)}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
