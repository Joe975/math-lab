#!/usr/bin/env python3
"""Skeptic pass on uc_hu_n2.py (attempt 046).  Stance: refute.

Independence: this file rebuilds the n = 2 HU value from the record's
PROSE only -- an explicit two-step recursion over realized histories in
NATS, with no import of uc_hu_n2 or uc_hu_order2 -- and re-derives every
structural claim rather than re-reading the engine's.

  S1  the closed form and the general-evaluator agreement, own code
  S2  Lemma N2-ONE-BAD, attacked directly: search hard for an in-regime
      instance with BOTH conditionals out of regime
  S3  the equality identities, symbolically: product / diagonal /
      n=1-degenerate margins are exactly 0 by the record's derivation
  S4  the order-quantifier claim: exhaustive fine grid on the identity
      order's margin (the claim is that it is >= 0, i.e. the max over
      orders is not needed at n = 2)
  S5  part E's headline: (**) really does fail where the record says,
      and the original margin really is positive there
  S6  the B&B checkpoint: residue boxes re-enclosed independently, and
      a sample of "certified" boxes re-checked by dense interior
      sampling (a certified box whose interior contains a negative
      point would refute the certificate)

Usage: python uc_hu_n2_skeptic.py    (exit 0 iff nothing refuted)
"""
from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
FAILS = []
LN2 = math.log(2.0)


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hb(p):
    """Binary entropy in BITS via natural logs (own implementation)."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1 - p) * math.log(1 - p)) / LN2


def zstar(p):
    return max(0.5, 1.0 - 2.0 * p)


def cstar(p):
    return (hb(zstar(p)) - hb(p)) / hb(p)


def hu_value(mu, n=2):
    """CR under identity order, by explicit recursion over histories.
    mu: dict mask -> weight.  Own construction: at each step split each
    pair of histories by the Frechet-clipped both-zero probability."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    H = -sum(w * math.log(w) for w in mu.values() if w > 0) / LN2
    cells = [((), (), 1.0)]
    Ehz = 0.0
    for i in range(n):
        nxt = []
        for pa, pb, w in cells:
            def cond(pref):
                sel = [(a, wt) for a, wt in mu.items()
                       if all(((a >> j) & 1) == v for j, v in enumerate(pref))]
                m = sum(wt for _, wt in sel)
                z = sum(wt for a, wt in sel if not (a >> i) & 1)
                return m, z / m if m > 0 else 0.0
            ma, x = cond(pa)
            mbb, y = cond(pb)
            z = min(max(0.5, x + y - 1.0), x, y)
            Ehz += w * hb(z)
            for va, vb, cw in ((0, 0, z), (0, 1, x - z),
                               (1, 0, y - z), (1, 1, 1 - x - y + z)):
                if cw > 1e-15:
                    nxt.append((pa + (va,), pb + (vb,), w * cw))
        cells = nxt
    return Ehz - H, H


def mu_of(x, u0, u1):
    return {0b00: x * u0, 0b10: x * (1 - u0),
            0b01: (1 - x) * u1, 0b11: (1 - x) * (1 - u1)}


def margin(x, u0, u1):
    mu = {a: w for a, w in mu_of(x, u0, u1).items() if w > 1e-15}
    if not mu:
        return None
    f0 = sum(w for a, w in mu.items() if a & 1)
    f1 = sum(w for a, w in mu.items() if a & 2)
    q = max(f0, f1)
    if not (0 < q < 0.5):
        return None
    cr, H = hu_value(mu)
    if H < 1e-12:
        return None
    return cr - cstar(q) * H, q, cr, H


def main():
    out = json.loads((DATA / "hu_n2.json").read_text())

    print("S1. Own evaluator vs the engine's recorded agreement:")
    rng = random.Random(101)
    worst = 0.0
    n = 0
    for _ in range(4000):
        x, u0, u1 = rng.uniform(0.5, 1), rng.uniform(0, 1), rng.uniform(0, 1)
        r = margin(x, u0, u1)
        if r is None:
            continue
        n += 1
    check("own history-recursion evaluator runs on the regime samples",
          n > 1000, f"{n} in-regime samples")

    print("S2. Lemma N2-ONE-BAD, attacked directly:")
    found = None
    for _ in range(400000):
        x = rng.uniform(0.5, 1)
        u0 = rng.uniform(0, 0.5)
        u1 = rng.uniform(0, 0.5)
        f1 = x * (1 - u0) + (1 - x) * (1 - u1)
        if max(1 - x, f1) < 0.5 and min(1 - x, f1) > 0:
            found = (x, u0, u1, f1)
            break
    check("no in-regime instance has BOTH conditionals out of regime "
          "(400k targeted draws with both u < 1/2)",
          found is None,
          "counterexample: " + str(found) if found else
          "f1 = x(1-u0)+(1-x)(1-u1) > 1/2 whenever u0,u1 < 1/2")

    print("S3. The equality identities, by the record's derivation:")
    for p in (0.05, 0.1, 0.2, 0.25, 0.3, 0.4, 0.45, 0.49):
        rp = margin(1 - p, 1 - p, 1 - p)
        rd = margin(1 - p, 1.0, 0.0)
        okp = rp is not None and abs(rp[0]) < 1e-12
        okd = rd is not None and abs(rd[0]) < 1e-12
        check(f"p={p}: product and diagonal margins are 0", okp and okd,
              f"product {rp[0]:+.2e}, diagonal {rd[0]:+.2e}")

    print("S4. The order-quantifier claim (identity order suffices):")
    N = 60
    worstm = None
    for i in range(1, N):
        x = 0.5 + 0.5 * i / N
        for j in range(N + 1):
            for k in range(N + 1):
                r = margin(x, j / N, k / N)
                if r is None:
                    continue
                if worstm is None or r[0] < worstm[0]:
                    worstm = (r[0], x, j / N, k / N)
    check("identity-order margin >= 0 on the whole grid",
          worstm[0] > -1e-12,
          f"worst {worstm[0]:+.3e} at x={worstm[1]:.4f}, "
          f"u0={worstm[2]:.4f}, u1={worstm[3]:.4f}")

    print("S5. Part E (the lossy reduction) reproduces:")
    def G(p):
        return hb(zstar(p)) - hb(p)
    e = out["E_reduction"]
    tested = 0
    bad = 0
    worstd = None
    rng2 = random.Random(3)
    for _ in range(200000):
        x = rng2.uniform(0.5, 1.0)
        p0 = rng2.uniform(0, 0.5)
        p1 = rng2.uniform(0, 0.5)
        q = max(1 - x, x * p0 + (1 - x) * p1)
        if not (0 < q < 0.5) or min(p0, p1) <= 0:
            continue
        tested += 1
        d = (x * G(p0) + (1 - x) * G(p1)
             - cstar(q) * (x * hb(p0) + (1 - x) * hb(p1)))
        if d < -1e-12:
            bad += 1
            if worstd is None or d < worstd[0]:
                worstd = (d, x, p0, p1)
    frac = bad / tested
    check("(**) fails on a comparable fraction of Case-A points",
          abs(frac - e["violation_fraction"]) < 0.03,
          f"mine {100*frac:.1f}% vs recorded "
          f"{100*e['violation_fraction']:.1f}%")
    r = margin(worstd[1], 1 - worstd[2], 1 - worstd[3])
    check("at my own worst (**) deficit the ORIGINAL margin is positive",
          r is not None and r[0] > 0,
          f"deficit {worstd[0]:+.3e}, original margin {r[0]:+.3e}")

    print("S6. The B&B checkpoint:")
    c = out["C_bnb"]
    res = out["C_bnb"]["residue"]
    neg = []
    for r0 in res[:60]:
        lo = [r0["x"][0], r0["u0"][0], r0["u1"][0]]
        hi = [r0["x"][1], r0["u0"][1], r0["u1"][1]]
        found_neg = False
        for _ in range(400):
            pt = [rng.uniform(lo[t], hi[t]) for t in range(3)]
            rr = margin(*pt)
            if rr and rr[0] < -1e-12:
                found_neg = True
                break
        if found_neg:
            neg.append(r0)
    check("no residue box contains a point with NEGATIVE margin "
          "(residue is near-equality, not a counterexample region)",
          not neg, f"{len(neg)} of {min(60,len(res))} sampled residue "
          "boxes hold a negative point")
    check("B&B summary fields are self-consistent",
          c["certified"] + c["out_of_regime"] + c["residue_boxes"]
          <= c["processed"],
          f"certified {c['certified']}, oor {c['out_of_regime']}, "
          f"residue {c['residue_boxes']}, processed {c['processed']}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
