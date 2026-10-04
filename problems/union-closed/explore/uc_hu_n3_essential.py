#!/usr/bin/env python3
"""Attempt 051 (correction to 050 part D): the n = 3 best-order descents
saturate the bound by DEGENERATING, not by finding an n = 3 extremal.

050 reported "min own-constant margin +0.000000, 0 violations" for the
best-order descents at n = 3 and read the saturation as the equality
family binding.  Inspecting the endpoints shows what actually happened:

    cap 0.45 endpoint: mu = {000: 0.5538, 001: 0.4462, 010: 0.0}
    cap 0.49 endpoint: mu = {000: 0.5105, 011: 0.4895, 111: 0.0}

both have one atom driven to weight 0 and only ONE coordinate carrying
mass -- they are n = 1 instances embedded in n = 3, where equality is
automatic.  The descent escaped to a degenerate face; it did not find
an n = 3 near-extremal, and 050's part D therefore says less than it
appears to.

This file re-runs the descent with an ESSENTIALITY constraint -- every
coordinate must keep marginal >= 0.05 and at least four atoms must
keep weight >= 1e-3 -- so the adversary is confined to genuinely
three-dimensional instances, and certifies the resulting floor.

Parts:
  A  the degeneracy, documented at 050's own endpoints
  B  essential-support descents at caps 0.45 / 0.49, own-constant
     flagged (the honest n = 3 best-order floor)
  C  exact certification of the sharpest essential floor, each kit
     alone

Standard library only; deterministic.
Usage: python uc_hu_n3_essential.py [--fast]
Checkpoint: ../data/hu_n3_essential.json
"""
from __future__ import annotations

import itertools
import json
import math
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import hu_cr_seq, h
from uc_hu_rollcensus import gen_instance
from uc_hu_certify import (hu_cells, ivl_add, ivl_scale, xlog2x_ivl, h_ivl,
                           log2_A, log2_B)

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
MIN_MARG = 0.05
MIN_ATOMS = 4
MIN_W = 1e-3


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def essential(n, mu):
    """Genuinely n-dimensional: every coordinate carries mass, and the
    support has not collapsed."""
    tot = sum(mu.values())
    if tot <= 0:
        return None
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-14}
    margs = [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(n)]
    if min(margs) < MIN_MARG or max(margs) >= 0.5:
        return None
    if sum(1 for w in m.values() if w >= MIN_W) < MIN_ATOMS:
        return None
    return m, max(margs)


def best_margin(n, mu):
    e = essential(n, mu)
    if e is None:
        return None
    m, fmax = e
    best, H = None, None
    for s in itertools.permutations(range(n)):
        cr, H = hu_cr_seq(n, m, s)
        if best is None or cr > best:
            best = cr
    if H is None or H < 0.2:
        return None
    return best / H - cstar(fmax)


def part_A():
    log("A. 050's part-D endpoints, inspected:")
    d = json.loads((DATA / "hu_n3_census.json").read_text())
    rows = []
    for capkey in ("D_cap_0.45", "D_cap_0.49"):
        r = min(d[capkey]["rows"], key=lambda z: z["own_margin"])
        mu = {int(s, 2): w for s, w in r["mu"].items()}
        tot = sum(mu.values())
        m = {a: w / tot for a, w in mu.items()}
        margs = [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(3)]
        live = sum(1 for w in m.values() if w >= MIN_W)
        rows.append({"block": capkey, "start": r["start"],
                     "own_margin": r["own_margin"], "marginals": margs,
                     "live_atoms": live,
                     "is_essential": essential(3, mu) is not None})
        log(f"  {capkey} ({r['start']}): margin {r['own_margin']:+.2e}, "
            f"coordinate marginals {[round(x,4) for x in margs]}, "
            f"atoms with weight >= 1e-3: {live} -> "
            f"{'essential' if rows[-1]['is_essential'] else 'DEGENERATE (n=1 face)'}")
    log("  => 050's saturation is the descent escaping to a degenerate "
        "face, not an n = 3 extremal")
    OUT["A_degeneracy"] = rows


def descend_essential(name, mu, cap, rounds=120):
    cur = best_margin(3, mu)
    if cur is None:
        return None
    step = 0.4
    rounds_done = 0
    while step > 0.02 and rounds_done < (rounds if not FAST else 20):
        improved = False
        for a in sorted(mu):
            for f in (1 + step, 1 / (1 + step)):
                cand = dict(mu)
                cand[a] = cand[a] * f
                v = best_margin(3, cand)
                if v is not None and v < cur - 1e-12:
                    mu, cur, improved = cand, v, True
        for extra in (0, 7):
            cand = dict(mu)
            tot = sum(mu.values())
            cand[extra] = cand.get(extra, 0.0) + 0.08 * tot
            v = best_margin(3, cand)
            if v is not None and v < cur - 1e-12:
                mu, cur, improved = cand, v, True
        if not improved:
            step *= 0.5
        rounds_done += 1
    tot = sum(mu.values())
    return {"start": name, "floor": cur,
            "mu": {format(a, "03b"): w / tot for a, w in mu.items()}}


def part_B():
    log()
    log("B. Essential-support descents (every coordinate marginal "
        f">= {MIN_MARG}, >= {MIN_ATOMS} atoms above {MIN_W}):")
    for cap in (0.45, 0.49):
        rng = random.Random(6100 + int(cap * 1000))
        rows, gmin = [], 1e9
        tries = 0
        while len(rows) < (10 if not FAST else 3) and tries < 200:
            tries += 1
            mu = gen_instance(rng, 3, cap)
            if essential(3, mu) is None:
                continue
            r = descend_essential(f"ess{len(rows)}", mu, cap)
            if r is None:
                continue
            rows.append(r)
            gmin = min(gmin, r["floor"])
        viol = [r for r in rows if r["floor"] < -1e-12]
        log(f"  cap {cap}: {len(rows)} descents, min own-constant margin "
            f"{gmin:+.6f}, {len(viol)} violations")
        OUT[f"B_cap_{cap}"] = {"rows": rows, "min_margin": gmin,
                               "violations": viol}


def part_C():
    log()
    log("C. Certifying the sharpest essential floor (each kit alone):")
    best = None
    for cap in (0.45, 0.49):
        blk = OUT.get(f"B_cap_{cap}")
        if not blk or not blk["rows"]:
            continue
        r = min(blk["rows"], key=lambda z: z["floor"])
        if best is None or r["floor"] < best[0]["floor"]:
            best = (r, cap)
    if best is None:
        log("  no essential endpoint to certify")
        return
    r, cap = best
    muF = {int(s, 2): w for s, w in r["mu"].items()}
    muQ = {a: Fraction(w).limit_denominator(10 ** 7) for a, w in muF.items()}
    tot = sum(muQ.values())
    muQ = {a: w / tot for a, w in muQ.items()}
    marg = max(sum(w for a, w in muQ.items() if (a >> i) & 1)
               for i in range(3))
    seq = max(itertools.permutations(range(3)),
              key=lambda s: hu_cr_seq(3, {a: float(w) for a, w in muQ.items()},
                                      s)[0])
    perm = [0] * 3
    for slot, coord in enumerate(seq):
        perm[coord] = slot
    muP = {}
    for a, w in muQ.items():
        b = 0
        for i in range(3):
            if (a >> i) & 1:
                b |= 1 << perm[i]
        muP[b] = muP.get(b, Fraction(0)) + w
    row = {"start": r["start"], "cap": cap, "best_seq": list(seq),
           "max_marginal": float(marg), "float_floor": r["floor"],
           "kits": {}}
    for kitname, fn in (("A_digit_extraction", log2_A),
                        ("B_atanh_series", log2_B)):
        Hmu = (Fraction(0), Fraction(0))
        for w in muP.values():
            Hmu = ivl_add(Hmu, ivl_scale(Fraction(-1), xlog2x_ivl(w, fn)))
        S = (Fraction(0), Fraction(0))
        for w, z in hu_cells(3, muP):
            if 0 < z < 1:
                S = ivl_add(S, ivl_scale(w, h_ivl(z, fn)))
        CR = (S[0] - Hmu[1], S[1] - Hmu[0])
        ok = CR[0] > 0 and Hmu[0] > Fraction(1, 2) and marg < Fraction(1, 2)
        row["kits"][kitname] = {"CR_lo": float(CR[0]), "CR_hi": float(CR[1]),
                                "H_lo": float(Hmu[0]), "certifies": ok}
        log(f"  [{kitname}]: CR_HU in [{float(CR[0]):+.9e}, "
            f"{float(CR[1]):+.9e}], marg {float(marg):.5f} "
            f"{'CERTIFIED > 0' if ok else 'NOT CERTIFIED'}")
        if not ok:
            raise SystemExit("certification failed")
    OUT["C_certificate"] = row


def main():
    part_A()
    part_B()
    part_C()
    (DATA / "hu_n3_essential.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_n3_essential.json")


if __name__ == "__main__":
    main()
