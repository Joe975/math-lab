#!/usr/bin/env python3
"""Attempt 041 (040 lead 1): interpolate between (HU-TAX)'s two equality
families and look for a dip below the constant.

040 identified the two-point diagonal (1-p)delta_empty + p delta_full
as a second equality family of the corrected (HU-TAX) constant
c*(p) = (h(max(1/2, 1-2p)) - h(p))/h(p), at the opposite correlation
extreme from products.  If the conjecture is false in the mixture
world, the natural place to look is BETWEEN the two equality cases.
Two paths, both exchangeable mixtures of product measures:

  A  convex mixture   mu_rho = (1-rho) Bern(p)^n + rho * diag(p)
                      (= mixture of Bern(p), Bern(0), Bern(1) products)
  B  product-mixture  mu_t = (1-lam) Bern(a)^n + lam Bern(b)^n with
     a = p(1-t), b = p + (1-p)t, lam = (p-a)/(b-a): t=0 the product,
     t=1 exactly the diagonal -- the path stays inside the
     mixture-of-products genre the recipe's latent-mixing branch owns.

Because both paths are K <= 3 product mixtures, the HU coupling
collapses exactly to a DP over (ones-in-a-prefix, ones-in-b-prefix)
count pairs -- (k+1)^2 cells at step k -- so the sweep runs to n = 64.
The DP is cross-checked against the general evaluator (uc_hu_order2
hu_cr_seq) at n <= 6, and order-independence (exchangeability) is
spot-checked against the rollout order.

For every grid point: margin(mu) = CR_HU/H - c*(p).  p is the exact
marginal of every coordinate along both paths (checked at the DP
level by construction, at the atom level in the crosscheck), so the
margin compares against the instance's OWN constant (034's lesson).
Any negative margin refutes (HU-TAX) in its re-posed form.

Standard library only; deterministic (no RNG).
Usage: python uc_hu_interp.py [--fast]
Checkpoint: ../data/hu_interp.json
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import hu_cr_seq, seq_roll, h

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def zclip(x, y):
    return min(max(0.5, x + y - 1.0), x, y)


def H_pmix(n, comps):
    """Exact H(mu) in bits for a mixture of product measures (counts)."""
    tot = 0.0
    for m in range(n + 1):
        pm = sum(lam * q ** m * (1 - q) ** (n - m) for q, lam in comps)
        if pm > 0.0:
            tot -= math.comb(n, m) * pm * math.log2(pm)
    return tot


def hu_cr_pmix(n, comps):
    """(CR, H) of the HU coupling for an exchangeable product mixture,
    by exact DP over (ones-in-a-prefix, ones-in-b-prefix) cells.  The
    posterior over components given a prefix depends only on its ones
    count, so per-cell conditionals are exact."""
    def cond_zero(k, j):
        # P(next coordinate = 0 | prefix of length k with j ones)
        num = den = 0.0
        for q, lam in comps:
            w = lam * q ** j * (1 - q) ** (k - j)
            den += w
            num += w * (1 - q)
        return num / den if den > 0 else None

    cells = {(0, 0): 1.0}
    Ehz = 0.0
    for k in range(n):
        xz = {}
        nxt = {}
        for (ja, jb), w in cells.items():
            if ja not in xz:
                xz[ja] = cond_zero(k, ja)
            if jb not in xz:
                xz[jb] = cond_zero(k, jb)
            x, y = xz[ja], xz[jb]
            z = zclip(x, y)
            Ehz += w * h(z)
            for da, db, cw in ((0, 0, z), (0, 1, x - z),
                               (1, 0, y - z), (1, 1, 1.0 - x - y + z)):
                if cw > 1e-18:
                    key = (ja + da, jb + db)
                    nxt[key] = nxt.get(key, 0.0) + w * cw
        cells = nxt
    H = H_pmix(n, comps)
    return Ehz - H, H


def comps_A(p, rho):
    return [(p, 1.0 - rho), (0.0, rho * (1.0 - p)), (1.0, rho * p)]


def comps_B(p, t):
    if t <= 0.0:
        return [(p, 1.0)]
    a = p * (1 - t)
    b = p + (1 - p) * t
    lam = (p - a) / (b - a)
    return [(a, 1.0 - lam), (b, lam)]


def atoms_of(n, comps):
    mu = {}
    for a in range(1 << n):
        m = bin(a).count("1")
        w = sum(lam * q ** m * (1 - q) ** (n - m) for q, lam in comps)
        if w > 0.0:
            mu[a] = w
    return mu


def sweep(pathname, compfn, ns, ps, grid):
    rows = []
    worst = (1e9, None)
    for n in ns:
        for p in ps:
            for g in grid:
                cr, H = hu_cr_pmix(n, compfn(p, g))
                m = cr / H - cstar(p)
                rows.append({"n": n, "p": p, "t": g, "margin": m,
                             "CR": cr, "H": H})
                if m < worst[0]:
                    worst = (m, (n, p, g))
            sub = [r for r in rows if r["n"] == n and r["p"] == p]
            mn = min(sub, key=lambda r: r["margin"])
            log(f"  {pathname} n={n} p={p}: min margin {mn['margin']:+.6f} "
                f"at t={mn['t']:.2f}")
    neg = [r for r in rows if r["margin"] < -1e-12]
    OUT[pathname] = {"rows": rows, "worst_margin": worst[0],
                     "worst_at": worst[1], "negatives": len(neg)}
    log(f"  {pathname}: worst margin {worst[0]:+.8f} at {worst[1]}, "
        f"{len(neg)} negatives")


def part_C_checks():
    """DP vs the general evaluator (n <= 6, identity order), and
    exchangeability (identity vs rollout order) on the atom measures."""
    log("C. Cross-checks: DP vs general evaluator, and order freedom:")
    bad_dp, bad_ord = [], []
    for n, p, t, cf in ((2, 0.49, 0.5, comps_A), (4, 0.45, 0.3, comps_A),
                        (6, 0.49, 0.85, comps_A), (4, 0.49, 0.5, comps_B),
                        (6, 0.38271, 0.7, comps_B), (5, 0.30, 0.15, comps_B)):
        comps = cf(p, t)
        crd, Hd = hu_cr_pmix(n, comps)
        mu = atoms_of(n, comps)
        dev = max(abs(sum(w for a, w in mu.items() if (a >> i) & 1) - p)
                  for i in range(n))
        ci, Hi = hu_cr_seq(n, mu, tuple(range(n)))
        cr, _ = hu_cr_seq(n, mu, seq_roll(n, mu))
        if abs(crd - ci) > 1e-10 or abs(Hd - Hi) > 1e-10 or dev > 1e-12:
            bad_dp.append((cf.__name__, n, p, t, crd, ci, dev))
        if abs(ci - cr) > 1e-12:
            bad_ord.append((cf.__name__, n, p, t, ci, cr))
    log(f"  DP vs evaluator: {'OK (<=1e-10, marginals exact)' if not bad_dp else 'MISMATCH ' + str(bad_dp)}")
    log(f"  identity vs rollout: {'OK (<=1e-12)' if not bad_ord else 'MISMATCH ' + str(bad_ord)}")
    OUT["C_crosschecks"] = {"bad_dp": bad_dp, "bad_order": bad_ord}


def main():
    ns = (2, 4, 8, 16, 32, 64) if not FAST else (2, 4)
    ps = (0.30, 0.38271, 0.45, 0.49) if not FAST else (0.49,)
    grid = [k / 40 for k in range(41)] if not FAST else [k / 10 for k in range(11)]
    log("A. Convex mixture path product -> diagonal:")
    sweep("A_convex", comps_A, ns, ps, grid)
    log()
    log("B. Product-mixture path (within the mixture genre):")
    sweep("B_prodmix", comps_B, ns, ps, grid)
    log()
    part_C_checks()
    for pathname, cf in (("A_convex", comps_A), ("B_prodmix", comps_B)):
        n, p, g0 = OUT[pathname]["worst_at"]
        fine = [max(0.0, min(1.0, g0 + k / 400)) for k in range(-10, 11)]
        best = (1e9, None)
        for g in fine:
            cr, H = hu_cr_pmix(n, cf(p, g))
            m = cr / H - cstar(p)
            if m < best[0]:
                best = (m, g)
        OUT[pathname]["refined_min"] = {"margin": best[0], "t": best[1],
                                        "n": n, "p": p}
        log(f"  {pathname} refined min: {best[0]:+.8f} at t={best[1]:.4f} "
            f"(n={n}, p={p})")
    (DATA / "hu_interp.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_interp.json")


if __name__ == "__main__":
    main()
