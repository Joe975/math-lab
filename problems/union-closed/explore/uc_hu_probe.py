#!/usr/bin/env python3
"""Attempt 031: the half-union (HU) coupling -- definition and
adversarial sweep.

HU coupling of mu with itself, per history cell (a, b) at coordinate i:
x = mu(A_i=0|a), y = mu(B_i=0|b), both-zero probability
z = min(max(1/2, x + y - 1), x, y) -- the union marginal pushed to 1/2
within the cell's Frechet bounds.  Total (every mu), closed-form (no
Sinkhorn, no lambda).  Per-coordinate idea = Sawin's max-entropy
coupling ingredient / 009 part D's adaptive rule; this file is the
per-history totalization and its adversarial sweep.

Parts:
  A  floor endpoints of cr_attack.json + crash family + product identity
  B  fixed-seed random-support scans at caps 0.38271 / 0.40 / 0.43 /
     0.46 / 0.49 (delta_0-projected into regime), plus deterministic
     pattern-search descents on CR/H from the worst starts
  C  slice DP (exchangeable (w_a, w_b) states) at w0/n up to 0.48,
     n up to 300

Determinism: the only RNG is random.Random with the fixed seeds written
below (911000+, 913000+); everything else is exact iteration.
Usage: python uc_hu_probe.py [--fast]
Checkpoint: ../data/hu_probe.json
"""
from __future__ import annotations

import json
import math
import random
import sys
import time
from math import lgamma, log2
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_mitax import cr_eval, crash_family, h

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def half_union_pairs(n, mu):
    """The HU coupling as an explicit pairs dict."""
    pref = []
    for i in range(n):
        mask = (1 << i) - 1
        d = {}
        for A, w in mu.items():
            k = A & mask
            e = d.setdefault(k, [0.0, 0.0])
            e[0] += w
            if not (A >> i) & 1:
                e[1] += w
        pref.append(d)
    cells = {(0, 0): 1.0}
    for i in range(n):
        nxt = {}
        for (a, b), w in cells.items():
            da, db = pref[i].get(a), pref[i].get(b)
            x = da[1] / da[0]
            y = db[1] / db[0]
            z = min(max(0.5, x + y - 1.0), x, y)
            for (ai, bi), cw in {(0, 0): z, (0, 1): x - z, (1, 0): y - z,
                                 (1, 1): 1 - x - y + z}.items():
                if cw <= 1e-15:
                    continue
                kk = (a | (0 if ai == 0 else 1 << i),
                      b | (0 if bi == 0 else 1 << i))
                nxt[kk] = nxt.get(kk, 0.0) + w * cw
        cells = nxt
    return cells


def hu_eval(n, mu):
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w > 1e-12}
    return cr_eval(n, half_union_pairs(n, m))


def maxmarg(n, mu):
    tot = sum(mu.values())
    return max(sum(w for a, w in mu.items() if (a >> i) & 1)
               for i in range(n)) / tot


# ------------------------------------------------------------------- A
def part_A():
    log("A. Recorded hard instances under HU")
    rows = []
    d = json.loads((DATA / "cr_attack.json").read_text())
    for row in d["rows"]:
        if row["H_mu"] < 0.5:
            continue
        n = row["n"]
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        r = hu_eval(n, mu)
        rows.append({"instance": f"floor:{row['seed']}", "n": n,
                     "prior_CR": row["sup_cr"], "CR_HU": r["CR"],
                     "H": r["H_mu"], "ratio": r["CR"] / r["H_mu"]})
        log(f"  floor:{row['seed']:13s} n={n}: sweep {row['sup_cr']:+.4f}"
            f" -> HU {r['CR']:+.4f} (ratio {r['CR']/r['H_mu']:+.4f})")
        assert r["CR"] > row["sup_cr"] - 1e-9
    for n in (5, 8, 11, 14, 17):
        atoms, logm = crash_family(n)
        mu = {a: 2.0 ** lm for a, lm in zip(atoms, logm)}
        r = hu_eval(n, mu)
        rows.append({"instance": "crash", "n": n, "CR_HU": r["CR"],
                     "H": r["H_mu"], "ratio": r["CR"] / r["H_mu"]})
    log(f"  crash family n=5..17: CR_HU = {rows[-1]['CR_HU']:+.5f} "
        f"(flat in n; ratio {rows[-1]['ratio']:+.4f})")
    # product identity CR = n(1 - h(p))
    for p, n in [(0.3, 5), (0.38271, 6), (0.45, 5)]:
        mu = {a: p ** bin(a).count("1")
              * (1 - p) ** (n - bin(a).count("1"))
              for a in range(1 << n)}
        r = hu_eval(n, mu)
        dev = abs(r["CR"] - n * (1 - h(p)))
        assert dev < 1e-9, (p, n, dev)
    log("  product identity CR = n(1-h(p)): checked to 1e-9 at "
        "p in {0.3, 0.38271, 0.45}")
    OUT["A_instances"] = rows


# ------------------------------------------------------------------- B
def scan_and_descend(cap, nscan, seed_base):
    def ratio_of(n, mu):
        if maxmarg(n, mu) >= cap:
            return None
        r = hu_eval(n, mu)
        if r["H_mu"] < 0.2:
            return None
        return r["CR"] / r["H_mu"], r

    worst = (1e9, None)
    nviol = 0
    pool = []
    for seed in range(nscan):
        rng = random.Random(seed_base + seed)
        n = rng.choice([3, 4, 5, 6])
        k = rng.randint(2, min(14, 1 << n))
        supp = rng.sample(range(1 << n), k)
        mu = {a: math.exp(rng.uniform(-2.5, 2.5)) for a in supp}
        tot = sum(mu.values())
        mu = {a: w / tot for a, w in mu.items()}
        fm = maxmarg(n, mu)
        if fm >= cap - 0.005:
            al = 1.0 - (cap - 0.005) / fm
            mu = {a: w * (1 - al) for a, w in mu.items()}
            mu[0] = mu.get(0, 0.0) + al
        r = ratio_of(n, mu)
        if r is None:
            continue
        if r[1]["CR"] < -1e-12:
            nviol += 1
        if r[0] < worst[0]:
            worst = (r[0], seed)
        pool.append((r[0], seed, n, mu))
    floors = []
    for _, seed, n, mu0 in sorted(pool)[:3]:
        mu = dict(mu0)
        best = ratio_of(n, mu)[0]
        step = 0.5
        while step >= 0.03:
            improved = False
            for a in sorted(mu):
                for f in (1 + step, 1 / (1 + step)):
                    cand = dict(mu)
                    cand[a] = cand[a] * f
                    r = ratio_of(n, cand)
                    if r and r[0] < best - 1e-9:
                        mu, best, improved = cand, r[0], True
            if not improved:
                step *= 0.5
        floors.append(best)
    return {"cap": cap, "scanned": len(pool), "violations": nviol,
            "worst_scan_ratio": worst[0], "descent_floors": floors,
            "product_extremal": (1 - h(cap)) / h(cap)}


def part_B():
    log()
    log("B. Fixed-seed scans + descents by marginal cap:")
    rows = []
    caps = [(0.38271, 300, 911000)] + \
        [(c, 200, 913000) for c in (0.40, 0.43, 0.46, 0.49)]
    if FAST:
        caps = [(0.38271, 40, 911000), (0.46, 30, 913000)]
    for cap, nscan, base in caps:
        r = scan_and_descend(cap, nscan, base)
        rows.append(r)
        log(f"  cap {cap:.5f}: {r['scanned']} cases, "
            f"{r['violations']} violations; worst "
            f"{r['worst_scan_ratio']:+.5f}, descent floors "
            f"{[f'{x:+.5f}' for x in r['descent_floors']]}, "
            f"extremal {r['product_extremal']:+.5f}")
        assert r["violations"] == 0
    OUT["B_caps"] = rows


# ------------------------------------------------------------------- C
def slice_hu_cr(n, w0):
    states = {(0, 0): 1.0}
    sum_hz = 0.0
    for i in range(n):
        rem = n - i
        nxt = {}
        for (wa, wb), P in states.items():
            x = 1.0 - (w0 - wa) / rem
            y = 1.0 - (w0 - wb) / rem
            z = min(max(0.5, x + y - 1.0), x, y)
            if 0.0 < z < 1.0:
                sum_hz += P * h(z)
            for (ai, bi), cw in {(0, 0): z, (0, 1): x - z, (1, 0): y - z,
                                 (1, 1): 1 - x - y + z}.items():
                if cw <= 1e-14:
                    continue
                k = (wa + ai, wb + bi)
                nxt[k] = nxt.get(k, 0.0) + P * cw
        states = nxt
    H = (lgamma(n + 1) - lgamma(w0 + 1) - lgamma(n - w0 + 1)) \
        / math.log(2.0)
    return sum_hz - H, H


def part_C():
    log()
    log("C. Slices (exchangeable-state DP):")
    rows = []
    ns = (50, 100, 200, 300) if not FAST else (50,)
    for frac in (0.36, 0.40, 0.42, 0.45, 0.48):
        line = []
        for n in ns:
            w0 = round(frac * n)
            cr, H = slice_hu_cr(n, w0)
            rows.append({"frac": frac, "n": n, "w0": w0, "CR": cr,
                         "H": H, "ratio": cr / H})
            line.append(f"n={n}: {cr:+8.3f} ({cr/H:+.4f})")
            assert cr > 0
        log(f"  w0/n={frac:.2f} [extremal {(1-h(frac))/h(frac):+.5f}]: "
            + "  ".join(line))
    OUT["C_slices"] = rows


# ------------------------------------------------------------------- D
def part_D():
    """Scope checks: HU vs optimized K-mixing (030 battery, n = 4), and
    coordinate-order dependence at the floor instances."""
    import itertools
    from uc_latentK_probe import latentK_pairs, q_from_t, optimize_t
    log()
    log("D. Scope checks:")
    insts = {
        "lo+mid+hi": [(0.5, 0.01), (0.3, 0.55), (0.2, 0.9)],
        "spread": [(0.6, 0.05), (0.25, 0.7), (0.15, 0.99)],
        "near0+2hi": [(0.62, 0.001), (0.22, 0.9), (0.16, 0.99)],
        "all-low+1mid": [(0.4, 0.2), (0.35, 0.3), (0.25, 0.55)],
        "near-psi": [(0.4, 0.36), (0.3, 0.39), (0.3, 0.40)],
    }
    hh = []
    for tag, comps in insts.items():
        n = 4
        mu = {}
        for P, p in comps:
            if p <= 0.0:
                mu[0] = mu.get(0, 0.0) + P
                continue
            for a in range(1 << n):
                k = bin(a).count("1")
                mu[a] = mu.get(a, 0.0) + P * p ** k * (1 - p) ** (n - k)
        r = hu_eval(n, mu)
        km, _, _ = optimize_t(n, comps)
        hh.append({"tag": tag, "HU": r["CR"], "Kmix": km})
        log(f"  {tag:12s} n=4: HU {r['CR']:+.5f} vs K-mix* {km:+.5f}"
            f"  ({'HU' if r['CR'] >= km else 'K-mix'} wins)")
    OUT["D_headtohead"] = hh

    d = json.loads((DATA / "cr_attack.json").read_text())
    orows = []
    for row in d["rows"]:
        if row["H_mu"] < 0.5:
            continue
        n = row["n"]
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        perms = itertools.permutations(range(n)) if n <= 5 else \
            itertools.islice(itertools.permutations(range(n)), 120)
        vals = []
        for perm in perms:
            mu2 = {}
            for a, w in mu.items():
                b = 0
                for i in range(n):
                    if (a >> i) & 1:
                        b |= 1 << perm[i]
                mu2[b] = mu2.get(b, 0.0) + w
            vals.append(hu_eval(n, mu2)["CR"])
        orows.append({"seed": row["seed"], "n": n, "min": min(vals),
                      "max": max(vals), "identity": hu_eval(n, mu)["CR"]})
        log(f"  order-dep {row['seed']:13s} n={n}: CR in "
            f"[{min(vals):+.5f}, {max(vals):+.5f}] "
            f"(identity {orows[-1]['identity']:+.5f})")
        assert min(vals) > 0
    OUT["D_orders"] = orows


def main():
    part_A()
    part_B()
    part_C()
    part_D()
    (DATA / "hu_probe.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_probe.json")


if __name__ == "__main__":
    main()
