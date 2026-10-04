#!/usr/bin/env python3
"""Skeptic pass on uc_hu_order2.py (attempt 037).  Default stance: refute.

Independence: rule logic re-implemented here in NATS with a different
cell representation (atom-index partitions instead of value-tuple
dicts), and every CR value is recomputed through the pre-existing,
separately-validated stack (permute_mu -> half_union_pairs -> cr_eval),
which shares no code with the engine's hu_cr_seq.

  S1  three witnesses: every rule's recorded sequence and CR reproduce
      (own rule implementations; independent evaluator)
  S2  part C descent endpoints: floor values and rule sequences
      reproduce at every endpoint; global floors and violation lists
      consistent
  S3  part D endpoints: best-order ratio re-derived by full
      enumeration with the independent evaluator
  S4  part B summaries: worst value per rule per cap re-derived from
      the 51 committed 035 endpoints

Usage: python uc_hu_order2_skeptic.py     (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_mitax import cr_eval
from uc_hu_probe import half_union_pairs
from uc_hu_attack import permute_mu

DATA = HERE.parent / "data"
LN2 = math.log(2.0)
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hnat(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1 - p) * math.log(1 - p))


def cr_indep(n, mu, seq):
    """CR/H via the independent stack: seq -> perm, permute, evaluate."""
    perm = [0] * n
    for slot, coord in enumerate(seq):
        perm[coord] = slot
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items()}
    r = cr_eval(n, half_union_pairs(n, permute_mu(n, m, perm)))
    return r["CR"], r["H_mu"]


# ---- own cell machinery: partition atom INDICES by realized prefixes
def cells_after(n, atoms, seq_prefix):
    """List of ((idxA tuple), (idxB tuple), weight) coupling cells after
    revealing seq_prefix, built by direct recursion in nats-free
    arithmetic."""
    cells = [(tuple(range(len(atoms))), tuple(range(len(atoms))), 1.0)]
    for c in seq_prefix:
        nxt = {}
        for ia, ib, w in cells:
            ma = sum(atoms[k][1] for k in ia)
            mb = sum(atoms[k][1] for k in ib)
            a0 = tuple(k for k in ia if not (atoms[k][0] >> c) & 1)
            b0 = tuple(k for k in ib if not (atoms[k][0] >> c) & 1)
            x = sum(atoms[k][1] for k in a0) / ma
            y = sum(atoms[k][1] for k in b0) / mb
            z = min(max(0.5, x + y - 1.0), x, y)
            a1 = tuple(k for k in ia if k not in a0)
            b1 = tuple(k for k in ib if k not in b0)
            for A, B, cw in ((a0, b0, z), (a0, b1, x - z),
                             (a1, b0, y - z), (a1, b1, 1 - x - y + z)):
                if cw > 1e-15 and A and B:
                    key = (A, B)
                    nxt[key] = nxt.get(key, 0.0) + w * cw
        cells = [(A, B, w) for (A, B), w in nxt.items()]
    return cells


def condH_nats(n, atoms, i, S):
    d = {}
    for a, w in atoms:
        key = tuple((a >> j) & 1 for j in S)
        e = d.setdefault(key, [0.0, 0.0])
        e[0] += w
        if not (a >> i) & 1:
            e[1] += w
    return sum(m * hnat(z / m) for m, z in d.values() if m > 0)


def step_surplus_nats(n, atoms, cells, i):
    s = 0.0
    for ia, ib, w in cells:
        ma = sum(atoms[k][1] for k in ia)
        mb = sum(atoms[k][1] for k in ib)
        x = sum(atoms[k][1] for k in ia if not (atoms[k][0] >> i) & 1) / ma
        y = sum(atoms[k][1] for k in ib if not (atoms[k][0] >> i) & 1) / mb
        z = min(max(0.5, x + y - 1.0), x, y)
        s += w * (hnat(z) - 0.5 * (hnat(x) + hnat(y)))
    return s


TIE_NATS = 1e-12 * LN2


def my_canon_completion(n, atoms, chosen):
    seq = list(chosen)
    rem = sorted(set(range(n)) - set(chosen))
    while rem:
        scored = sorted((condH_nats(n, atoms, i, seq), i) for i in rem)
        lo = scored[0][0]
        pick = min(i for v, i in scored if v <= lo + TIE_NATS)
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq)


def my_rule(n, mu, kind):
    tot = sum(mu.values())
    atoms = [(a, w / tot) for a, w in mu.items()]
    if kind == "can":
        return my_canon_completion(n, atoms, [])
    seq, rem = [], sorted(range(n))
    while rem:
        if kind == "canst":
            scored = sorted((condH_nats(n, atoms, i, seq), i) for i in rem)
            lo = scored[0][0]
            tied = [i for v, i in scored if v <= lo + TIE_NATS]
            if len(tied) > 1:
                cells = cells_after(n, atoms, seq)
                pick = max(tied, key=lambda i:
                           (step_surplus_nats(n, atoms, cells, i), -i))
            else:
                pick = tied[0]
        elif kind == "surp":
            cells = cells_after(n, atoms, seq)
            scored = sorted((-step_surplus_nats(n, atoms, cells, i), i)
                            for i in rem)
            lo = scored[0][0]
            pick = min(i for v, i in scored if v <= lo + TIE_NATS)
        elif kind == "roll":
            best = None
            for i in rem:
                full = my_canon_completion(n, atoms, seq + [i])
                cr, _ = cr_indep(n, dict(mu), full)
                if best is None or cr > best[0] + 1e-12:
                    best = (cr, i)
                elif abs(cr - best[0]) <= 1e-12 and i < best[1]:
                    best = (best[0], i)
            pick = best[1]
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq)


def main():
    out = json.loads((DATA / "hu_order2.json").read_text())
    can = json.loads((DATA / "hu_canon.json").read_text())
    att = json.loads((DATA / "hu_attack.json").read_text())

    print("S1. Witness rows: rules + values, independently:")
    src = {"033witness": att["cap_0.45"]["violations"][0],
           "kill_n4": can["B_cap_0.49"]["violations"][0],
           "kill_n5": can["B_cap_0.49"]["violations"][1]}
    for row in out["A_witnesses"]:
        v = src[row["tag"]]
        n = row["n"]
        mu = {int(s, 2): w for s, w in v["mu"].items()}
        allcr = sorted(cr_indep(n, mu, s)[0]
                       for s in itertools.permutations(range(n)))
        check(f"{row['tag']}: worst/best/negative-count reproduce",
              abs(allcr[0] - row["worst_CR"]) < 1e-9
              and abs(allcr[-1] - row["best_CR"]) < 1e-9
              and sum(1 for c in allcr if c < 0) == row["n_negative"],
              f"worst {allcr[0]:+.6f}, best {allcr[-1]:+.6f}")
        for rn in ("can", "canst", "surp", "roll"):
            seq = my_rule(n, mu, rn)
            cr, _ = cr_indep(n, mu, seq)
            check(f"{row['tag']} {rn}: sequence and CR reproduce",
                  list(seq) == row[rn]["seq"]
                  and abs(cr - row[rn]["CR"]) < 1e-9,
                  f"mine {seq} {cr:+.6f} vs {tuple(row[rn]['seq'])} "
                  f"{row[rn]['CR']:+.6f}")

    print("S2. Part C endpoints:")
    for key, blk in out.items():
        if not (key.startswith("C_") or key.startswith("C2_")):
            continue
        _, rn, _, caps = key.split("_")
        cap = float(caps)
        bad = []
        for r in blk["rows"]:
            n = r["n"]
            mu = {int(s, 2): w for s, w in r["mu"].items()}
            seq = my_rule(n, mu, rn)
            cr, H = cr_indep(n, mu, seq)
            if abs(cr / H - r["floor"]) > 1e-8:
                bad.append((r["start"], cr / H, r["floor"]))
        floors = [r["floor"] for r in blk["rows"]]
        check(f"{key}: every endpoint floor reproduces under own rule + "
              "independent evaluator", not bad, f"mismatches {bad[:2]}")
        check(f"{key}: global floor and violations consistent",
              abs(min(floors) - blk["global_floor"]) < 1e-12
              and sum(1 for f in floors if f < 0) == len(blk["violations"]),
              f"min {min(floors):+.6f}")

    print("S3. Part D endpoints (best-order by full enumeration):")
    for key, blk in out.items():
        if not key.startswith("D_"):
            continue
        bad = []
        for r in blk["rows"]:
            n = r["n"]
            mu = {int(s, 2): w for s, w in r["mu"].items()}
            best = max(cr_indep(n, mu, s)[0]
                       for s in itertools.permutations(range(n)))
            _, H = cr_indep(n, mu, tuple(range(n)))
            if abs(best / H - r["floor"]) > 1e-8:
                bad.append((r["start"], best / H, r["floor"]))
        check(f"{key}: every endpoint best-order ratio reproduces",
              not bad, f"mismatches {bad[:2]}")

    print("S4. Part B summaries from the 51 committed 035 endpoints:")
    for capkey in ("B_cap_0.38271", "B_cap_0.45", "B_cap_0.49"):
        cap = float(capkey.split("_")[-1])
        summ = out["B_endpoints"][capkey]
        for rn in ("can", "canst", "surp", "roll"):
            worst, neg = 1e9, 0
            for row in can[capkey]["rows"]:
                n = row["n"]
                mu = {int(s, 2): w for s, w in row["mu"].items()}
                tot = sum(mu.values())
                m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
                fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
                           for i in range(n))
                if fmax >= cap:
                    continue
                seq = my_rule(n, m, rn)
                cr, H = cr_indep(n, m, seq)
                if H < 0.2:
                    continue
                v = cr / H
                worst = min(worst, v)
                neg += 1 if v < 0 else 0
            check(f"{capkey} {rn}: worst and negative-count reproduce",
                  abs(worst - summ[rn]["worst"]) < 1e-8
                  and neg == summ[rn]["negative_rows"],
                  f"worst {worst:+.6f} vs {summ[rn]['worst']:+.6f}, "
                  f"neg {neg}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
