#!/usr/bin/env python3
"""Attempt 037 (035 lead 1 + lead 3): NON-GREEDY / SURPLUS-SCORED order
rules for the half-union coupling, and the best-order ceiling.

035 showed the canonical (greedy conditional-entropy) order rescues cap
0.45 but dies certified at 0.49, picking the WORST of 120 orders on its
n=5 kill: minimising per-step entropy banks predictability early while
CR needs surplus cells late.  This file scores the surplus/deficit
ledger DIRECTLY.  Since H(mu) telescopes along the same revelation
order, CR decomposes per cell as

    CR = sum_k sum_cells w(a,b) * [ h(z) - (h(x) + h(y))/2 ],

so s(x, y) = h(z) - (h(x)+h(y))/2 is the cell's surplus (z the HU
clip).  Rules tested, each a function of (mu, revealed SET) only -- so
each gives a genuine total coupling:

  can    035's canonical rule (baseline; via uc_hu_canon)
  canst  canonical, but exact ties broken by step surplus (036's S5
         note: lowest-index tie-breaks can cost value)
  surp   greedy MAX immediate step surplus (the ledger, one step)
  roll   rollout: pick the coordinate maximising full CR when the
         remainder is completed by the canonical rule (the ledger, to
         the end; non-greedy in the sense 035 lead 1 asks for)

Parts:
  A  the three record witnesses (033 order-kill, both 035 cap-0.49
     kills): every rule vs the full order enumeration
  B  all 51 endpoint measures of 035's canonical descents, re-scored
     under every rule (were those kills rule-specific?)
  C  fresh weight-descent attacks AGAINST surp and roll at caps
     0.49 / 0.497 (order re-derived by the rule at every candidate)
  D  the best-order ceiling (035 lead 3): descend on max-over-ALL-
     orders CR at caps 0.49 / 0.497, n <= 5 starts -- what ANY order
     rule could at best achieve
  C2 hostile seeding (runs after D): attack roll from the canonical and
     surp kill measures and the sharpest C/D endpoints, caps
     0.49 / 0.497 / 0.499

Standard library only; deterministic (fixed seeds in start shapes only).
Usage: python uc_hu_order2.py [--fast]
Checkpoint: ../data/hu_order2.json
"""
from __future__ import annotations

import itertools
import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_canon import canonical_order, cond_entropy, TIE_TOL
from uc_hu_attack import starts

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
HMIN = 0.2
log2 = math.log2


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def h(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * log2(p) + (1 - p) * log2(1 - p))


def zclip(x, y):
    return min(max(0.5, x + y - 1.0), x, y)


# --------------------------------------------------------- HU evaluation
def prefix_tab(atoms, seq_prefix, i):
    """value-tuple on seq_prefix -> [mass, zero-mass of coord i]."""
    d = {}
    for a, w in atoms:
        key = tuple((a >> j) & 1 for j in seq_prefix)
        e = d.setdefault(key, [0.0, 0.0])
        e[0] += w
        if not (a >> i) & 1:
            e[1] += w
    return d


def hu_cr_seq(n, mu, seq):
    """(CR, H) of the HU coupling under revelation sequence seq."""
    atoms = list(mu.items())
    H = -sum(w * log2(w) for _, w in atoms if w > 0.0)
    cells = {((), ()): 1.0}
    Ehz = 0.0
    pre = []
    for k in range(n):
        tab = prefix_tab(atoms, pre, seq[k])
        nxt = {}
        for (pa, pb), w in cells.items():
            ma, za = tab[pa]
            mb, zb = tab[pb]
            x, y = za / ma, zb / mb
            z = zclip(x, y)
            Ehz += w * h(z)
            for va, vb, cw in ((0, 0, z), (0, 1, x - z),
                               (1, 0, y - z), (1, 1, 1.0 - x - y + z)):
                if cw > 1e-15:
                    key = (pa + (va,), pb + (vb,))
                    nxt[key] = nxt.get(key, 0.0) + w * cw
        cells = nxt
        pre.append(seq[k])
    return Ehz - H, H


def step_surplus(cells, tab):
    """Immediate ledger surplus of revealing the tab's coordinate now."""
    s = 0.0
    for (pa, pb), w in cells.items():
        ma, za = tab[pa]
        mb, zb = tab[pb]
        x, y = za / ma, zb / mb
        s += w * (h(zclip(x, y)) - 0.5 * (h(x) + h(y)))
    return s


def advance(cells, tab):
    nxt = {}
    for (pa, pb), w in cells.items():
        ma, za = tab[pa]
        mb, zb = tab[pb]
        x, y = za / ma, zb / mb
        z = zclip(x, y)
        for va, vb, cw in ((0, 0, z), (0, 1, x - z),
                           (1, 0, y - z), (1, 1, 1.0 - x - y + z)):
            if cw > 1e-15:
                key = (pa + (va,), pb + (vb,))
                nxt[key] = nxt.get(key, 0.0) + w * cw
    return nxt


# ------------------------------------------------------------ order rules
def canon_completion(n, mu, chosen):
    """035's greedy rule continued from a given prefix (same TIE_TOL,
    lowest-index break)."""
    seq = list(chosen)
    remaining = sorted(set(range(n)) - set(chosen))
    while remaining:
        scored = sorted((cond_entropy(n, mu, i, seq), i) for i in remaining)
        lo = scored[0][0]
        best = min(i for v, i in scored if v <= lo + TIE_TOL)
        seq.append(best)
        remaining.remove(best)
    return tuple(seq)


def seq_can(n, mu):
    perm = canonical_order(n, mu)
    seq = [0] * n
    for coord, slot in enumerate(perm):
        seq[slot] = coord
    return tuple(seq)


def seq_canst(n, mu):
    """Canonical, exact ties broken by immediate step surplus (then
    lowest index)."""
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items()}
    atoms = list(m.items())
    seq = []
    remaining = sorted(set(range(n)))
    cells = {((), ()): 1.0}
    while remaining:
        scored = sorted((cond_entropy(n, m, i, seq), i) for i in remaining)
        lo = scored[0][0]
        tied = [i for v, i in scored if v <= lo + TIE_TOL]
        if len(tied) > 1:
            tabs = {i: prefix_tab(atoms, seq, i) for i in tied}
            best = max(tied, key=lambda i: (step_surplus(cells, tabs[i]), -i))
        else:
            best = tied[0]
        tab = prefix_tab(atoms, seq, best)
        cells = advance(cells, tab)
        seq.append(best)
        remaining.remove(best)
    return tuple(seq)


def seq_surp(n, mu):
    """Greedy max immediate step surplus; ties (1e-12) -> lowest index."""
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items()}
    atoms = list(m.items())
    seq = []
    remaining = sorted(set(range(n)))
    cells = {((), ()): 1.0}
    while remaining:
        tabs = {i: prefix_tab(atoms, seq, i) for i in remaining}
        scored = sorted(((-step_surplus(cells, tabs[i]), i))
                        for i in remaining)
        lo = scored[0][0]
        best = min(i for v, i in scored if v <= lo + 1e-12)
        cells = advance(cells, tabs[best])
        seq.append(best)
        remaining.remove(best)
    return tuple(seq)


def seq_roll(n, mu):
    """Rollout: at each step pick the coordinate maximising full CR with
    the remainder completed canonically; ties (1e-12) -> lowest index."""
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items()}
    seq = []
    remaining = sorted(set(range(n)))
    while remaining:
        scored = []
        for i in remaining:
            full = canon_completion(n, m, seq + [i])
            cr, _ = hu_cr_seq(n, m, full)
            scored.append((-cr, i))
        scored.sort()
        lo = scored[0][0]
        best = min(i for v, i in scored if v <= lo + 1e-12)
        seq.append(best)
        remaining.remove(best)
    return tuple(seq)


RULES = {"can": seq_can, "canst": seq_canst, "surp": seq_surp,
         "roll": seq_roll}


def rule_ratio(n, mu, cap, rulefn):
    """CR/H under the rule, with 035's regime filters; None if
    out-of-regime."""
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
    fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
               for i in range(n))
    if fmax >= cap:
        return None
    cr, H = hu_cr_seq(n, m, rulefn(n, m))
    if H < HMIN:
        return None
    return cr / H


def best_order_ratio(n, mu, cap):
    """max over ALL n! orders of CR, over H (the oracle ceiling)."""
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
    fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
               for i in range(n))
    if fmax >= cap:
        return None
    H = -sum(w * log2(w) for w in m.values())
    if H < HMIN:
        return None
    best = -1e9
    for s in itertools.permutations(range(n)):
        cr, _ = hu_cr_seq(n, m, s)
        if cr > best:
            best = cr
    return best / H


# ------------------------------------------------------------- descent
def descend(name, n, mu0, cap, objective, rounds_max=200, step_floor=0.03):
    """035's move set, objective(n, mu, cap) -> ratio or None."""
    mu = dict(mu0)
    best = objective(n, mu, cap)
    if best is None:
        return None
    step = 0.5
    rounds = 0
    while step >= step_floor and rounds < (rounds_max if not FAST else 15):
        improved = False
        for a in sorted(mu):
            for f in (1 + step, 1 / (1 + step)):
                cand = dict(mu)
                cand[a] = cand[a] * f
                v = objective(n, cand, cap)
                if v is not None and v < best - 1e-9:
                    mu, best, improved = cand, v, True
        tot = sum(mu.values())
        for extra in (0, (1 << n) - 1):
            cand = dict(mu)
            cand[extra] = cand.get(extra, 0.0) + 0.1 * tot
            v = objective(n, cand, cap)
            if v is not None and v < best - 1e-9:
                mu, best, improved = cand, v, True
        if len(mu) > 3:
            a0 = min(mu, key=mu.get)
            cand = {a: w for a, w in mu.items() if a != a0}
            v = objective(n, cand, cap)
            if v is not None and v < best - 1e-9:
                mu, best, improved = cand, v, True
        if not improved:
            step *= 0.5
        rounds += 1
    tot = sum(mu.values())
    return {"start": name, "n": n, "cap": cap, "floor": best,
            "atoms": len(mu),
            "mu": {format(a, f"0{n}b"): w / tot for a, w in mu.items()}}


# ---------------------------------------------------------------- parts
def part_A():
    log("A. The three record witnesses under every rule:")
    can = json.loads((DATA / "hu_canon.json").read_text())
    att = json.loads((DATA / "hu_attack.json").read_text())
    tests = [("033witness", att["cap_0.45"]["violations"][0], 0.45),
             ("kill_n4", can["B_cap_0.49"]["violations"][0], 0.49),
             ("kill_n5", can["B_cap_0.49"]["violations"][1], 0.49)]
    rows = []
    for tag, v, cap in tests:
        n = v["n"]
        mu = {int(s, 2): w for s, w in v["mu"].items()}
        tot = sum(mu.values())
        m = {a: w / tot for a, w in mu.items()}
        H = -sum(w * log2(w) for w in m.values())
        allcr = sorted(hu_cr_seq(n, m, s)[0]
                       for s in itertools.permutations(range(n)))
        row = {"tag": tag, "n": n, "cap": cap, "H": H,
               "worst_CR": allcr[0], "best_CR": allcr[-1],
               "n_orders": len(allcr),
               "n_negative": sum(1 for c in allcr if c < 0)}
        msg = []
        for rn, fn in RULES.items():
            seq = fn(n, m)
            cr, _ = hu_cr_seq(n, m, seq)
            rank = 1 + sum(1 for c in allcr if c < cr - 1e-12)
            row[rn] = {"seq": list(seq), "CR": cr, "rank": rank}
            msg.append(f"{rn} {cr:+.5f} (r{rank})")
        log(f"  {tag:11s} n={n}: worst {allcr[0]:+.5f}, best {allcr[-1]:+.5f} "
            f"| " + ", ".join(msg))
        rows.append(row)
    OUT["A_witnesses"] = rows


def part_B():
    log()
    log("B. 035's 51 canonical-descent endpoints re-scored under every rule:")
    can = json.loads((DATA / "hu_canon.json").read_text())
    summary = {}
    for capkey in ("B_cap_0.38271", "B_cap_0.45", "B_cap_0.49"):
        cap = float(capkey.split("_")[-1])
        worst = {rn: (1e9, None) for rn in RULES}
        neg = {rn: 0 for rn in RULES}
        for row in can[capkey]["rows"]:
            n = row["n"]
            mu = {int(s, 2): w for s, w in row["mu"].items()}
            for rn, fn in RULES.items():
                v = rule_ratio(n, mu, cap, fn)
                if v is None:
                    continue
                if v < 0:
                    neg[rn] += 1
                if v < worst[rn][0]:
                    worst[rn] = (v, row["start"])
        summary[capkey] = {rn: {"worst": worst[rn][0],
                                "worst_start": worst[rn][1],
                                "negative_rows": neg[rn]}
                           for rn in RULES}
        log(f"  {capkey}: " + ", ".join(
            f"{rn} worst {worst[rn][0]:+.6f} ({neg[rn]} neg)"
            for rn in RULES))
    OUT["B_endpoints"] = summary


def part_C():
    log()
    log("C. Fresh weight-descent attacks against surp and roll:")
    caps = (0.49, 0.497) if not FAST else (0.49,)
    for rn in ("surp", "roll"):
        fn = RULES[rn]
        for cap in caps:
            rows, gmin = [], 1e9
            for name, n, mu in starts(cap):
                if rn == "roll" and n > 5 and not FAST:
                    # rollout descent at n=6 is ~25x an n=4 run; keep the
                    # budget on the caps, log the skip
                    log(f"    (skip {name}: n={n} rollout descent -- "
                        f"budget; covered by surp)")
                    continue
                r = descend(name, n, mu, cap,
                            lambda n_, m_, c_: rule_ratio(n_, m_, c_, fn))
                if r is None:
                    continue
                rows.append(r)
                gmin = min(gmin, r["floor"])
            viol = [r for r in rows if r["floor"] < 0]
            log(f"  {rn} cap {cap}: {len(rows)} starts, global floor "
                f"{gmin:+.6f}, {len(viol)} violations")
            OUT[f"C_{rn}_cap_{cap}"] = {"rows": rows, "global_floor": gmin,
                                        "violations": viol}


def part_C2():
    """Hostile seeding: attack roll from the measures that killed the
    other rules (035's two canonical kills, this file's surp kills), and
    a cap-0.499 probe on the sharpest endpoints."""
    log()
    log("C2. Hostile-seeded descents against roll + cap-0.499 probe:")
    can = json.loads((DATA / "hu_canon.json").read_text())
    seeds = [(f"canonkill:{r['start']}", r["n"],
              {int(s, 2): w for s, w in r["mu"].items()})
             for r in can["B_cap_0.49"]["violations"]]
    for key in ("C_surp_cap_0.49", "C_surp_cap_0.497"):
        if key in OUT:
            for r in OUT[key]["violations"]:
                seeds.append((f"surpkill:{r['start']}@{r['cap']}", r["n"],
                              {int(s, 2): w for s, w in r["mu"].items()}))
    for key in ("C_roll_cap_0.49", "C_roll_cap_0.497",
                "D_bestorder_cap_0.497"):
        if key in OUT:
            rows = sorted(OUT[key]["rows"], key=lambda r: r["floor"])[:2]
            for r in rows:
                seeds.append((f"sharp:{key}:{r['start']}", r["n"],
                              {int(s, 2): w for s, w in r["mu"].items()}))
    for cap in ((0.49, 0.497, 0.499) if not FAST else (0.49,)):
        rows, gmin = [], 1e9
        for name, n, mu in seeds:
            if n > 5:
                continue
            r = descend(name, n, mu, cap,
                        lambda n_, m_, c_: rule_ratio(n_, m_, c_, seq_roll))
            if r is None:
                continue
            rows.append(r)
            gmin = min(gmin, r["floor"])
        viol = [r for r in rows if r["floor"] < 0]
        log(f"  roll cap {cap}: {len(rows)} hostile seeds, global floor "
            f"{gmin:+.6f}, {len(viol)} violations "
            f"[product extremal {(1 - h(cap)) / h(cap):+.6f}]")
        OUT[f"C2_roll_cap_{cap}"] = {"rows": rows, "global_floor": gmin,
                                     "violations": viol}


def part_D():
    log()
    log("D. The best-order ceiling (oracle over all n! orders):")
    caps = (0.49, 0.497) if not FAST else (0.49,)
    for cap in caps:
        rows, gmin = [], 1e9
        for name, n, mu in starts(cap):
            if n > 5:
                continue
            r = descend(name, n, mu, cap, best_order_ratio,
                        rounds_max=60, step_floor=0.05)
            if r is None:
                continue
            rows.append(r)
            gmin = min(gmin, r["floor"])
        viol = [r for r in rows if r["floor"] < 0]
        log(f"  best-order cap {cap}: {len(rows)} starts (n <= 5), global "
            f"floor {gmin:+.6f}, {len(viol)} violations")
        OUT[f"D_bestorder_cap_{cap}"] = {"rows": rows, "global_floor": gmin,
                                         "violations": viol}


def main():
    part_A()
    part_B()
    part_C()
    part_D()
    part_C2()
    (DATA / "hu_order2.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_order2.json")


if __name__ == "__main__":
    main()
