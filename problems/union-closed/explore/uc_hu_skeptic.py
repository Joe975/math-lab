#!/usr/bin/env python3
"""Attempt 031 skeptic pass (half-union coupling).  Default stance: refute.

Independent paths:
  - HU coupling rebuilt from the definition with a different construction:
    the joint is assembled as an explicit measure over (A, B) mask pairs
    by recursive descent (B-major recursion), not the cell-product loop;
  - CR evaluated with uc_branch_skeptic.cr_chain (natural log, B-major);
  - slice values re-derived by a SECOND DP written state-major with
    natural logs, cross-checked against direct mask enumeration at n = 8;
  - marginal exactness of the HU build checked on every instance.
"""
from __future__ import annotations

import json
import math
import sys
from math import lgamma
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from uc_branch_skeptic import cr_chain, LN2
from uc_mitax import crash_family

HERE = Path(__file__).parent
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'REFUTED'}] {name}"
          + (f"  ({detail})" if detail else ""), flush=True)
    if not ok:
        FAILS.append(name)


def hu_pairs_recursive(n, mu):
    """Independent HU build: recursive descent over coordinates, carrying
    conditional measures explicitly (dict-of-masses per branch)."""
    out = {}

    def rec(i, condA, condB, wa_mask, wb_mask, weight):
        if weight <= 1e-16:
            return
        if i == n:
            out[(wa_mask, wb_mask)] = out.get((wa_mask, wb_mask), 0.0) \
                + weight
            return
        ta = sum(condA.values())
        tb = sum(condB.values())
        a0 = {A: w for A, w in condA.items() if not (A >> i) & 1}
        a1 = {A: w for A, w in condA.items() if (A >> i) & 1}
        b0 = {B: w for B, w in condB.items() if not (B >> i) & 1}
        b1 = {B: w for B, w in condB.items() if (B >> i) & 1}
        x = sum(a0.values()) / ta
        y = sum(b0.values()) / tb
        z = min(max(0.5, x + y - 1.0), x, y)
        for (sa, sb), cw in {(0, 0): z, (0, 1): x - z, (1, 0): y - z,
                             (1, 1): 1 - x - y + z}.items():
            if cw <= 1e-15:
                continue
            ca = a0 if sa == 0 else a1
            cb = b0 if sb == 0 else b1
            if not ca or not cb:
                continue
            rec(i + 1, ca, cb,
                wa_mask | (sa << i), wb_mask | (sb << i), weight * cw)

    rec(0, dict(mu), dict(mu), 0, 0, 1.0)
    return out


def marg_check(pairs, mu):
    ma, mb = {}, {}
    for (a, b), w in pairs.items():
        ma[a] = ma.get(a, 0.0) + w
        mb[b] = mb.get(b, 0.0) + w
    d = 0.0
    for a in set(mu) | set(ma) | set(mb):
        d = max(d, abs(ma.get(a, 0.0) - mu.get(a, 0.0)),
                abs(mb.get(a, 0.0) - mu.get(a, 0.0)))
    return d


def hnat(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log(p) - (1 - p) * math.log(1 - p)


def slice_hu_dp2(n, w0):
    """Second slice DP: state-major, natural logs."""
    states = {(0, 0): 1.0}
    acc = 0.0
    for i in range(n):
        rem = n - i
        nxt = {}
        for (wa, wb), P in sorted(states.items()):
            x = 1.0 - (w0 - wa) / rem
            y = 1.0 - (w0 - wb) / rem
            z = min(max(0.5, x + y - 1.0), x, y)
            acc += P * hnat(z)
            for (sa, sb), cw in {(0, 0): z, (0, 1): x - z, (1, 0): y - z,
                                 (1, 1): 1 - x - y + z}.items():
                if cw <= 1e-14:
                    continue
                k = (wa + sa, wb + sb)
                nxt[k] = nxt.get(k, 0.0) + P * cw
        states = nxt
    H = (lgamma(n + 1) - lgamma(w0 + 1) - lgamma(n - w0 + 1))
    return (acc - H) / LN2


def main():
    print("Skeptic pass: half-union coupling (031)")
    print("=" * 68)

    print("\nS1. Floor endpoints (rebuilt coupling, independent evaluator):")
    d = json.loads((Path(__file__).resolve().parent.parent / 'data' / 'cr_attack.json').read_text())
    claimed = {"windowkill": 0.3090, "mmabskill": 0.0822,
               "mmabskill_n6": 0.3007, "random0": 0.1563,
               "random1": 0.2991}
    for row in d["rows"]:
        if row["H_mu"] < 0.5:
            continue
        n = row["n"]
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        tot = sum(mu.values())
        mu = {a: w / tot for a, w in mu.items()}
        pairs = hu_pairs_recursive(n, mu)
        md = marg_check(pairs, mu)
        r = cr_chain(n, pairs)
        check(f"{row['seed']:14s}: CR_HU and marginals",
              abs(r["CR"] - claimed[row["seed"]]) < 5e-4 and md < 1e-12
              and r["CR"] > row["sup_cr"] - 1e-9,
              f"CR={r['CR']:+.4f} (claimed {claimed[row['seed']]:+.4f}, "
              f"sweep {row['sup_cr']:+.4f}), mdev={md:.0e}")

    print("\nS2. Crash family flat value:")
    for n in (5, 11, 17):
        atoms, logm = crash_family(n)
        mu = {a: 2.0 ** lm for a, lm in zip(atoms, logm)}
        tot = sum(mu.values())
        mu = {a: w / tot for a, w in mu.items()}
        pairs = hu_pairs_recursive(n, mu)
        r = cr_chain(n, pairs)
        check(f"n={n}: CR_HU = +0.19919, n-independent",
              abs(r["CR"] - 0.19919) < 5e-5,
              f"CR={r['CR']:+.5f}")

    print("\nS3. Slice DP: second implementation + n=8 mask enumeration:")
    for n, w0, ref in [(40, 16, 2.262), (100, 40, 4.542),
                       (300, 120, 11.071), (40, 12, 5.680)]:
        v = slice_hu_dp2(n, w0)
        check(f"n={n} w0={w0}: DP2 matches engine", abs(v - ref) < 2e-3,
              f"{v:+.4f} vs {ref:+.4f}")
    # n=8 slice: DP vs full-history mask build
    n, w0 = 8, 3
    from math import comb
    mu = {a: 1.0 / comb(n, w0) for a in range(1 << n)
          if bin(a).count("1") == w0}
    pairs = hu_pairs_recursive(n, mu)
    r = cr_chain(n, pairs)
    v = slice_hu_dp2(n, w0)
    check("n=8 slice: mask enumeration matches DP2",
          abs(r["CR"] - v) < 1e-9, f"{r['CR']:+.6f} vs {v:+.6f}")

    print("\nS4. Product extremal identity CR = n(1-h(p)):")
    import itertools
    for p, n in [(0.3, 5), (0.38271, 6), (0.45, 5)]:
        mu = {}
        for a in range(1 << n):
            k = bin(a).count("1")
            mu[a] = p ** k * (1 - p) ** (n - k)
        pairs = hu_pairs_recursive(n, mu)
        r = cr_chain(n, pairs)
        pred = n * (1 - (-p * math.log2(p) - (1 - p) * math.log2(1 - p)))
        check(f"Bern({p})^{n}: CR = n(1-h(p))",
              abs(r["CR"] - pred) < 1e-9,
              f"{r['CR']:+.6f} vs {pred:+.6f}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
    else:
        print("No refutation found.")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
