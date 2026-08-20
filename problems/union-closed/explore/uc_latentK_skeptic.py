#!/usr/bin/env python3
"""Attempt 030 skeptic pass.  Default stance: refute.

Independent of the engine: couplings rebuilt by direct (A, B)-mask
double loops from the definition (no product-DP, no shared builder);
CR evaluated with uc_branch_skeptic.cr_chain (natural-log, B-major).
Reads latentK_probe.json only to compare claims.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

HERE2 = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE2))
from uc_branch_skeptic import cr_chain, bern_measure

HERE = Path(__file__).resolve().parent.parent / "data"
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'REFUTED'}] {name}"
          + (f"  ({detail})" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)


def cell(p):
    if p > 0.5:
        return {(0, 0): 1.0 - p, (1, 1): p}
    if p <= 0.0:
        return {(0, 0): 1.0}
    m = 2.0 * p if p <= 0.25 else 0.5
    z = 1.0 - m
    off = 1.0 - p - z
    return {(0, 0): z, (0, 1): off, (1, 0): off, (1, 1): p - off}


def pair_mass(a, b, n, cj):
    w = 1.0
    for i in range(n):
        w *= cj.get(((a >> i) & 1, (b >> i) & 1), 0.0)
        if w == 0.0:
            return 0.0
    return w


def bmeas(n, p):
    if p <= 0.0:
        return {0: 1.0}
    return bern_measure(n, p)


def latentK(n, comps, t):
    """Direct rebuild from the definition, t: dict '{j}{k}' -> mass."""
    K = len(comps)
    P = [c[0] for c in comps]
    Q = [[0.0] * K for _ in range(K)]
    for jk, v in t.items():
        j, k = int(jk[0]), int(jk[1])
        Q[j][k] = Q[k][j] = v
    for j in range(K):
        Q[j][j] = P[j] - sum(Q[j][k] for k in range(K) if k != j)
        assert Q[j][j] > -1e-12
    pairs = {}
    ms = [bmeas(n, p) for _, p in comps]
    for j in range(K):
        for k in range(K):
            if Q[j][k] <= 1e-15:
                continue
            if j == k:
                cj = cell(comps[j][1])
                for a in range(1 << n):
                    for b in range(1 << n):
                        w = Q[j][j] * pair_mass(a, b, n, cj)
                        if w > 0.0:
                            pairs[(a, b)] = pairs.get((a, b), 0.0) + w
            else:
                for a, wa in ms[j].items():
                    for b, wb in ms[k].items():
                        w = Q[j][k] * wa * wb
                        if w > 0.0:
                            pairs[(a, b)] = pairs.get((a, b), 0.0) + w
    # marginal symmetry check
    ma, mb = {}, {}
    for (a, b), w in pairs.items():
        ma[a] = ma.get(a, 0.0) + w
        mb[b] = mb.get(b, 0.0) + w
    mdev = max(abs(ma.get(x, 0.0) - mb.get(x, 0.0))
               for x in set(ma) | set(mb))
    return pairs, mdev


INSTS = {
    "lo+mid+hi": [(0.5, 0.01), (0.3, 0.55), (0.2, 0.9)],
    "spread": [(0.6, 0.05), (0.25, 0.7), (0.15, 0.99)],
    "near0+2hi": [(0.62, 0.001), (0.22, 0.9), (0.16, 0.99)],
    "all-low+1mid": [(0.4, 0.2), (0.35, 0.3), (0.25, 0.55)],
    "near-psi": [(0.4, 0.36), (0.3, 0.39), (0.3, 0.40)],
}


def main():
    ref = json.loads((HERE / "latentK_probe.json").read_text())
    print("Skeptic pass on uc_latentK_probe.py")
    print("=" * 68)
    print("\nS1. Battery rows at the engine's optimal transfers:")
    for row in [r for r in ref["battery"] if r["n"] in (2, 6)]:
        comps = INSTS[row["tag"]]
        pairs, mdev = latentK(row["n"], comps, row["t"])
        r = cr_chain(row["n"], pairs)
        check(f"{row['tag']:12s} n={row['n']}: CR at t*",
              abs(r["CR"] - row["CR_opt"]) < 1e-9 and mdev < 1e-12
              and r["CR"] > 0,
              f"CR={r['CR']:+.6f} vs engine {row['CR_opt']:+.6f}, "
              f"mdev={mdev:.1e}")
        # diagonal (adaptive) value too
        pairs_d, _ = latentK(row["n"], comps,
                             {k: 0.0 for k in row["t"]})
        rd = cr_chain(row["n"], pairs_d)
        check(f"{row['tag']:12s} n={row['n']}: diagonal CR",
              abs(rd["CR"] - row["CR_diag"]) < 1e-9,
              f"{rd['CR']:+.6f} vs {row['CR_diag']:+.6f}")

    print("\nS2. Degeneration rows (n = 4):")
    for row in ref["degeneration"]:
        comps = [(0.62, row["p_light"]), (0.22, 0.9), (0.16, 0.99)]
        # engine's t not stored per-row with keys; re-optimize coarsely:
        # instead verify the CLAIM CR* > em_limit by evaluating the
        # saturated-t01 point the engine reported (t01 = 0.22).
        pairs, mdev = latentK(4, comps, {"01": 0.22, "02": 0.0, "12": 0.0})
        r = cr_chain(4, pairs)
        check(f"p_light={row['p_light']:.0e}: t01=0.22 point beats EM "
              f"limit",
              r["CR"] > row["em_limit"] and mdev < 1e-12,
              f"CR={r['CR']:+.5f} vs EM {row['em_limit']:+.5f}")

    print("\nS3. near-psi CR/H constancy claim (0.0430 at n=2,4,6):")
    vals = [r["cr_over_h"] for r in ref["battery"]
            if r["tag"] == "near-psi"]
    check("CR/H spread < 2e-4 across n", max(vals) - min(vals) < 2e-4,
          f"values {[f'{v:.5f}' for v in vals]}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
    else:
        print("No refutation found.")
    (HERE / "latentK_skeptic_out.txt").write_text(
        "done\n")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
