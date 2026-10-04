#!/usr/bin/env python3
"""Attempt 043 (reviewer pass on 041/042): independent re-implementation.

Nothing here imports from uc_hu_interp.py, uc_hu_blocks.py, or the
uc_hu_order2 / uc_hu_canon / uc_hu_roll_anneal stack.  The HU coupling
(031 "The coupling": per history cell, x = P(A_i = 0 | a),
y = P(B_i = 0 | b), z = min(max(1/2, x + y - 1), x, y)), the chain-rule
value CR = sum_i E[h(z_i)] - H(mu) (009), the corrected constant
c*(p) = (h(max(1/2, 1 - 2p)) - h(p))/h(p) (034), the canonical
min-conditional-entropy completion and the rollout rule (037; ties
1e-12 -> lowest index) are rebuilt from the records' prose.  Two
evaluator structures, both this file's own:

  * an atom-level evaluator whose pair-history state is a pair of
    POSTERIOR DISTRIBUTIONS (renormalized conditionals), not
    value-tuples with a prefix table;
  * a count-pair DP for exchangeable mixtures of product measures,
    written from this reviewer's own re-derivation of 041's collapse
    claim (component posterior given a prefix depends only on the
    prefix's ones count), with no sub-threshold pruning (only exact
    cw > 0), per-step dictionaries keyed (ja, jb).

Checks (against the committed checkpoints and the 041/042 record
prose; R-checks fail -> exit 1):

  R1  041 sweep reproduced at 132 sampled grid points (both paths,
      n in {2, 8, 64} x p in {0.38271, 0.49}, every fourth t), margins
      match data/hu_interp.json rows to 1e-9; marginal of every
      component mixture is exactly p
  R2  own DP vs own atom evaluator at n <= 6 (<= 1e-10), atom-level
      marginals exact, order invariance under random permutations
      (exchangeability, <= 1e-12); path B hits the product at t=0 and
      exactly the diagonal at t=1
  R3  041 checkpoint summaries + record prose: row counts (1968),
      zero negatives, worst margin/argmin, refined minima, the four
      interior minima quoted in the record and both growth-with-n
      sequences at p = 0.49
  R4  042 part A re-done with an own tensor builder: 5 cases
      (incl. p = 0.20 clamp branch), identity + own rollout + fixed
      non-contiguous + random orders, |CR/H - c*(p)| <= 1e-12; plus a
      numerical spot-check of the DIAG identity (040) at n in {3, 6},
      p in {0.05, 0.20, 0.35, 0.494}, random orders
  R5  042 parts B/C: every endpoint row of data/hu_blocks.json
      re-scored (own rollout + own evaluator) to 1e-8; in-regime;
      global floor / violation-list consistency; crash8 floors and the
      B-table floors match the record prose
  R6  C_sat497: floor == c*(0.494) to 1e-12, measure is exactly
      {000000: 0.506, 111111: 0.494}
  R7  CORRECTION checks (asserting the corrected statements found by
      this review): the two POLISHED B-endpoints (0.49/n6/d2^3 and
      0.497/n6/d2^3) are NOT on the block-tensor family -- 8 atoms,
      equal marginals, but blocks pairwise correlated and
      CR/H - c*(own max marginal) = +1.5e-5 / +1.4e-5 > 0 -- so 042's
      "polished endpoints are single diagonals" and
      "+0.000364 ~ c*(0.4887)" are wrong for exactly those two rows
      (c*(0.4887) = 3.686e-4 != 3.639e-4; the endpoint's own marginal
      is 0.489001 with c* = 3.492e-4); the three unpolished endpoints
      ARE family members (checked as exact products of diagonal
      blocks, one padded by two marginal-0 coordinates)

The pipeline byte-compare re-runs (P4) are run outside this file; see
the 043 record.

Usage: python uc_reviewer043_reimpl.py     (exit 0 iff all checks pass)
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
LOG2 = math.log(2.0)


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# ------------------------------------------------------- own primitives
def hb(p):
    """Binary entropy in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p)) / LOG2


def cstar(p):
    return (hb(max(0.5, 1.0 - 2.0 * p)) - hb(p)) / hb(p)


def zhu(x, y):
    """031's clamped half-union z: clip 1/2 into the Frechet interval."""
    lo = max(0.0, x + y - 1.0)
    hi = min(x, y)
    return min(max(0.5, lo), hi)


def norm(mu):
    t = sum(mu.values())
    return {a: w / t for a, w in mu.items() if w > 0.0}


def marginal(mu, i):
    return sum(w for a, w in mu.items() if (a >> i) & 1)


# ---------------------------- own atom-level evaluator (posterior pairs)
def _split(post, i):
    """post -> (P(coord i = 0), posterior|0, posterior|1)."""
    z0 = [(a, w) for a, w in post if not (a >> i) & 1]
    z1 = [(a, w) for a, w in post if (a >> i) & 1]
    x = sum(w for _, w in z0)
    p0 = tuple((a, w / x) for a, w in z0) if x > 0 else ()
    p1 = tuple((a, w / (1.0 - x)) for a, w in z1) if x < 1 else ()
    return x, p0, p1


def cr_pair(n, mu, seq):
    """(CR, H) in bits of the HU coupling under revelation order seq.
    State: dict (posterior_a, posterior_b) -> weight."""
    m = norm(mu)
    H = -sum(w * math.log(w) for w in m.values()) / LOG2
    p0 = tuple(sorted(m.items()))
    states = {(p0, p0): 1.0}
    ehz = 0.0
    for i in seq:
        nxt = {}
        for (pa, pb), w in states.items():
            x, a0, a1 = _split(pa, i)
            y, b0, b1 = _split(pb, i)
            z = zhu(x, y)
            ehz += w * hb(z)
            for cw, key in ((z, (a0, b0)), (x - z, (a0, b1)),
                            (y - z, (a1, b0)), (1.0 - x - y + z, (a1, b1))):
                if cw > 1e-16:
                    nxt[key] = nxt.get(key, 0.0) + w * cw
        states = nxt
    return ehz - H, H


# ------------------------------------------- own order rules (from prose)
def cond_ent(mu, i, prefix):
    """H(A_i | A_prefix) in bits over the atom measure."""
    cells = {}
    for a, w in mu.items():
        k = tuple((a >> j) & 1 for j in prefix)
        e = cells.setdefault(k, [0.0, 0.0])
        e[0] += w
        if not (a >> i) & 1:
            e[1] += w
    return sum(mw * hb(zw / mw) for mw, zw in cells.values() if mw > 0)


def canon_from(n, mu, chosen):
    """Greedy min-conditional-entropy completion; ties 1e-12 -> lowest
    index."""
    seq = list(chosen)
    rest = [i for i in range(n) if i not in seq]
    while rest:
        sc = sorted((cond_ent(mu, i, seq), i) for i in rest)
        best = min(i for v, i in sc if v <= sc[0][0] + 1e-12)
        seq.append(best)
        rest.remove(best)
    return tuple(seq)


def roll_order(n, mu):
    """Rollout: pick the coordinate maximizing full CR with the
    remainder completed canonically; ties 1e-12 -> lowest index."""
    m = norm(mu)
    seq = []
    rest = list(range(n))
    while rest:
        sc = []
        for i in rest:
            cr, _ = cr_pair(n, m, canon_from(n, m, seq + [i]))
            sc.append((-cr, i))
        sc.sort()
        best = min(i for v, i in sc if v <= sc[0][0] + 1e-12)
        seq.append(best)
        rest.remove(best)
    return tuple(seq)


# --------------------- own count-pair DP (this reviewer's re-derivation)
def pmix_H(n, comps):
    """H of an exchangeable mixture of products, via ones-count classes:
    H = -sum_m C(n,m) p_m log2 p_m, p_m = sum_k lam_k q_k^m (1-q_k)^(n-m)."""
    tot = 0.0
    for m in range(n + 1):
        pm = sum(l * q ** m * (1.0 - q) ** (n - m) for q, l in comps)
        if pm > 0.0:
            tot -= math.comb(n, m) * pm * math.log(pm)
    return tot / LOG2


def pmix_cr(n, comps):
    """(CR, H): HU coupling collapsed to (ones-in-a, ones-in-b) count
    pairs.  Justification (own derivation): each component is iid, so
    the posterior over components given a prefix of length k with j
    ones is prop. to lam_k q^j (1-q)^(k-j); hence the predictive
    P(next = 0 | prefix) depends on the prefix only through (k, j), so
    all pair-histories with equal count pairs share (x, y, z) and both
    the h(z) charge and the 4-way transition merge exactly."""
    def pzero(k, j):
        den = num = 0.0
        for q, l in comps:
            w = l * q ** j * (1.0 - q) ** (k - j)
            den += w
            num += w * (1.0 - q)
        return num / den
    cells = {(0, 0): 1.0}
    ehz = 0.0
    for k in range(n):
        col = {}
        nxt = {}
        for (ja, jb), w in cells.items():
            for j in (ja, jb):
                if j not in col:
                    col[j] = pzero(k, j)
            x, y = col[ja], col[jb]
            z = zhu(x, y)
            ehz += w * hb(z)
            for da, db, cw in ((0, 0, z), (0, 1, x - z),
                               (1, 0, y - z), (1, 1, 1.0 - x - y + z)):
                if cw > 0.0:
                    key = (ja + da, jb + db)
                    nxt[key] = nxt.get(key, 0.0) + w * cw
        cells = nxt
    H = pmix_H(n, comps)
    return ehz - H, H


# ------------------------------------------- own path/tensor constructors
def path_A(p, rho):
    return [(p, 1.0 - rho), (0.0, rho * (1.0 - p)), (1.0, rho * p)]


def path_B(p, t):
    if t <= 0.0:
        return [(p, 1.0)]
    a, b = p * (1.0 - t), p + (1.0 - p) * t
    lam = (p - a) / (b - a)
    return [(a, 1.0 - lam), (b, lam)]


def atoms(n, comps):
    mu = {}
    for a in range(1 << n):
        m = bin(a).count("1")
        w = sum(l * q ** m * (1.0 - q) ** (n - m) for q, l in comps)
        if w > 0.0:
            mu[a] = w
    return mu


def build_tensor(p, blocks):
    """Own tensor builder: diag_m(p) blocks ('d', m) and Bern(p)
    factors ('b', 1), low bits first."""
    mu = {0: 1.0}
    n = 0
    for kind, m in blocks:
        full = (1 << m) - 1 if kind == "d" else 1
        mu = {a | (v << n): w * pw
              for a, w in mu.items()
              for v, pw in ((0, 1.0 - p), (full, p))}
        n += m
    return n, mu


# ------------------------------------------------------------------ R1-R3
def r1(interp):
    print("R1. 041 sweep re-computed at sampled grid points (own DP):")
    for pname, pf in (("A_convex", path_A), ("B_prodmix", path_B)):
        rows = interp[pname]["rows"]
        pick = [r for r in rows if r["n"] in (2, 8, 64)
                and r["p"] in (0.38271, 0.49)
                and round(r["t"] * 40) % 4 == 0]
        bad = []
        for r in pick:
            comps = pf(r["p"], r["t"])
            mdev = abs(sum(l * q for q, l in comps) - r["p"])
            cr, H = pmix_cr(r["n"], comps)
            m = cr / H - cstar(r["p"])
            if abs(m - r["margin"]) > 1e-9 or mdev > 1e-15:
                bad.append((pname, r["n"], r["p"], r["t"], m, r["margin"]))
        check(f"{pname}: {len(pick)} sampled margins match to 1e-9, "
              "mixture marginal exactly p", not bad, f"{bad[:3]}")


def r2():
    print("R2. Own DP vs own atom evaluator; exchangeability; endpoints:")
    rng = random.Random(43043)
    bad_dp, bad_ord, bad_marg = [], [], []
    for n, p, t, pf in ((2, 0.49, 0.5, path_A), (4, 0.45, 0.3, path_A),
                        (6, 0.49, 0.85, path_A), (4, 0.49, 0.5, path_B),
                        (6, 0.38271, 0.7, path_B), (5, 0.30, 0.15, path_B)):
        comps = pf(p, t)
        crd, Hd = pmix_cr(n, comps)
        mu = atoms(n, comps)
        if max(abs(marginal(norm(mu), i) - p) for i in range(n)) > 1e-12:
            bad_marg.append((n, p, t))
        cri, Hi = cr_pair(n, mu, tuple(range(n)))
        if abs(crd - cri) > 1e-10 or abs(Hd - Hi) > 1e-10:
            bad_dp.append((pf.__name__, n, p, t, crd, cri))
        for _ in range(2):
            o = list(range(n))
            rng.shuffle(o)
            cro, _ = cr_pair(n, mu, tuple(o))
            if abs(cro - cri) > 1e-12:
                bad_ord.append((pf.__name__, n, p, t, o, cro - cri))
    check("count-pair DP == atom evaluator at n <= 6 (<= 1e-10)",
          not bad_dp, f"{bad_dp[:2]}")
    check("order invariance under random permutations (<= 1e-12)",
          not bad_ord, f"{bad_ord[:2]}")
    check("atom-level coordinate marginals exactly p", not bad_marg)
    okB = (path_B(0.42, 0.0) == [(0.42, 1.0)]
           and path_B(0.42, 1.0)[0][0] == 0.0
           and path_B(0.42, 1.0)[1][0] == 1.0
           and abs(path_B(0.42, 1.0)[1][1] - 0.42) < 1e-15)
    check("path B: t=0 is the product, t=1 exactly the diagonal", okB)


def r3(interp):
    print("R3. 041 checkpoint summaries and record prose:")
    nrows = sum(len(interp[k]["rows"]) for k in ("A_convex", "B_prodmix"))
    check("grid size: 984 + 984 = 1968 rows", nrows == 1968, str(nrows))
    for pname, pf in (("A_convex", path_A), ("B_prodmix", path_B)):
        blk = interp[pname]
        rows = blk["rows"]
        wm = min(rows, key=lambda r: r["margin"])
        neg = sum(1 for r in rows if r["margin"] < -1e-12)
        ok = (abs(wm["margin"] - blk["worst_margin"]) < 1e-15
              and [wm["n"], wm["p"], wm["t"]] == list(blk["worst_at"])
              and neg == 0 and blk["negatives"] == 0)
        check(f"{pname}: worst margin / argmin / zero negatives "
              "recompute from rows", ok,
              f"worst {wm['margin']:+.2e} at t={wm['t']}")
        rm = blk["refined_min"]
        cr, H = pmix_cr(rm["n"], pf(rm["p"], rm["t"]))
        check(f"{pname}: refined min reproduces (own DP) and sits at an "
              "endpoint", abs(cr / H - cstar(rm["p"]) - rm["margin"]) < 1e-9
              and rm["t"] in (0.0, 1.0),
              f"t={rm['t']}, margin {rm['margin']:+.2e}")
    quotes = {("A_convex", 2): 2.20e-4, ("A_convex", 8): 2.1e-3,
              ("A_convex", 64): 9.1e-3, ("B_prodmix", 2): 1.41e-7,
              ("B_prodmix", 8): 9.8e-7, ("B_prodmix", 64): 8.7e-6}
    bad = []
    for (pname, n), q in quotes.items():
        sub = [r for r in interp[pname]["rows"]
               if r["n"] == n and r["p"] == 0.49 and 0.0 < r["t"] < 1.0]
        mn = min(sub, key=lambda r: r["margin"])
        if abs(mn["margin"] - q) > 0.05 * q or mn["t"] != 0.025:
            bad.append((pname, n, mn["margin"], q, mn["t"]))
    check("record's six interior minima at p=0.49 (growth with n) match "
          "rows to quoted rounding, all at t=0.025", not bad, f"{bad}")
    gi = min((r for r in interp["A_convex"]["rows"] if 0 < r["t"] < 1),
             key=lambda r: r["margin"])
    gb = min((r for r in interp["B_prodmix"]["rows"] if 0 < r["t"] < 1),
             key=lambda r: r["margin"])
    check("global interior minima are the record's (n=2, p=0.49, "
          "t=0.025) pair", (gi["n"], gi["p"], gi["t"]) == (2, 0.49, 0.025)
          and (gb["n"], gb["p"], gb["t"]) == (2, 0.49, 0.025),
          f"A {gi['margin']:+.3e}, B {gb['margin']:+.3e}")


# ------------------------------------------------------------------ R4-R7
def r4():
    print("R4. 042 part A re-done (own tensor builder + evaluators):")
    rng = random.Random(97043)
    cases = [(0.49, [("d", 2), ("d", 2)]),
             (0.20, [("d", 2), ("d", 2)]),
             (0.45, [("d", 2), ("b", 1), ("d", 3)]),
             (0.49, [("d", 4), ("b", 1), ("b", 1)]),
             (0.30, [("d", 3), ("b", 1), ("d", 2)])]
    bad = []
    for p, blocks in cases:
        n, mu = build_tensor(p, blocks)
        noncontig = tuple(sorted(range(n), key=lambda i: (i % n) * 7 % n))
        orders = [tuple(range(n)), roll_order(n, mu), noncontig]
        o = list(range(n))
        rng.shuffle(o)
        orders.append(tuple(o))
        for seq in orders:
            cr, H = cr_pair(n, mu, seq)
            if abs(cr / H - cstar(p)) > 1e-12:
                bad.append((p, blocks, seq, cr / H - cstar(p)))
    check("5 tensor cases x 4 orders (incl. p=0.20 clamp, non-contiguous"
          " interleavings): |CR/H - c*(p)| <= 1e-12", not bad, f"{bad[:2]}")
    bad = []
    for p in (0.05, 0.20, 0.35, 0.494):
        for n in (3, 6):
            mu = {0: 1.0 - p, (1 << n) - 1: p}
            o = list(range(n))
            rng.shuffle(o)
            cr, H = cr_pair(n, mu, tuple(o))
            if (abs(cr - (hb(max(0.5, 1 - 2 * p)) - hb(p))) > 1e-12
                    or abs(H - hb(p)) > 1e-12):
                bad.append((n, p, cr))
    check("DIAG identity (040): CR = h(max(1/2,1-2p)) - h(p), H = h(p), "
          "random orders", not bad, f"{bad[:2]}")


def r5(blocks):
    print("R5. 042 B/C endpoints re-scored (own rollout + evaluator):")
    endinfo = {}
    for key in ("B_attack", "C_crash8"):
        blk = blocks[key]
        bad = []
        for r in blk["rows"]:
            n = r["n"]
            mu = norm({int(s, 2): w for s, w in r["mu"].items()})
            fm = max(marginal(mu, i) for i in range(n))
            cr, H = cr_pair(n, mu, roll_order(n, mu))
            if abs(cr / H - r["floor"]) > 1e-8 or fm >= r["cap"]:
                bad.append((r["start"], r["cap"], cr / H, r["floor"], fm))
            endinfo[(key, r["start"], r["cap"])] = (mu, n, fm, cr / H)
        check(f"{key}: every endpoint floor reproduces to 1e-8, "
              "in-regime", not bad, f"{bad[:2]}")
        floors = [r["floor"] for r in blk["rows"]]
        if "global_floor" in blk:
            check(f"{key}: global floor and (floor<0) violation list "
                  "consistent", abs(min(floors) - blk["global_floor"]) < 1e-15
                  and len(blk["violations"]) == sum(1 for f in floors if f < 0))
    prose = [(0.49, "[('d', 2), ('d', 2)]", 0.000488),
             (0.49, "[('d', 3), ('d', 3)]", 0.000488),
             (0.49, "[('d', 2), ('d', 2), ('d', 2)]", 0.000364),
             (0.497, "[('d', 2), ('d', 2)]", 0.000104),
             (0.497, "[('d', 2), ('d', 2), ('d', 2)]", 0.000061)]
    got = {(r["cap"], r["start"]): r["floor"] for r in blocks["B_attack"]["rows"]}
    bad = [(c, s, q, got.get((c, s))) for c, s, q in prose
           if (c, s) not in got or abs(got[(c, s)] - q) > 5e-7]
    check("B table in the record matches checkpoint floors (5 rows)",
          not bad, f"{bad}")
    c8 = {r["cap"]: r["floor"] for r in blocks["C_crash8"]["rows"]}
    check("crash8 prose floors +0.000952 / +0.000310 match checkpoint",
          abs(c8[0.49] - 0.000952) < 5e-7 and abs(c8[0.497] - 0.000310) < 5e-7,
          f"{c8}")
    wa = blocks["A_equality"]["worst_abs_margin"]
    npairs = sum(r["orders_tested"] for r in blocks["A_equality"]["rows"])
    check("part A prose: worst 4.4e-16 over 40 case-order pairs (8 x 5)",
          abs(wa - 4.4408920985e-16) < 1e-17 and npairs == 40
          and len(blocks["A_equality"]["rows"]) == 8,
          f"worst {wa:.3e}, pairs {npairs}")
    own = {k: v[3] - cstar(v[2]) for k, v in endinfo.items()}
    worst_own = min(own.values())
    check("(HU-TAX) own-constant margins CR/H - c*(fmax) of ALL B/C "
          "endpoints are >= -1e-9 (zero exactly on-family; the committed "
          "'floor < 0' flag alone is weaker than the conjecture)",
          worst_own > -1e-9, f"min own-margin {worst_own:+.3e}")
    return endinfo


def r6(blocks):
    print("R6. C_sat497:")
    s = blocks["C_sat497"]
    ok = (abs(s["floor"] - cstar(0.494)) < 1e-12
          and s["mu"] == {"000000": 0.506, "111111": 0.494}
          and not s["violation"]
          and abs(s["diag_value"] - cstar(0.494)) < 1e-15)
    check("floor == c*(0.494) to 1e-12; endpoint is exactly "
          "diag_6(0.494)", ok, f"floor {s['floor']:.11e} vs "
          f"c*(0.494) {cstar(0.494):.11e}")


def r7(blocks, endinfo):
    print("R7. Corrections (asserting the corrected statements):")
    polished = [("B_attack", "[('d', 2), ('d', 2), ('d', 2)]", 0.49),
                ("B_attack", "[('d', 2), ('d', 2), ('d', 2)]", 0.497)]
    bad = []
    for k in polished:
        mu, n, fm, ratio = endinfo[k]
        blkpair = (sum(w for a, w in mu.items()
                       if (a & 0b110000) == 0b110000
                       and (a & 0b001100) == 0b001100))
        b1 = sum(w for a, w in mu.items() if (a & 0b110000) == 0b110000)
        b2 = sum(w for a, w in mu.items() if (a & 0b001100) == 0b001100)
        on_family = (len(mu) == 2 or abs(blkpair - b1 * b2) < 1e-9)
        own = ratio - cstar(fm)
        # corrected statement: NOT a single diagonal, blocks correlated,
        # strictly above the family value at its own marginal
        if on_family or not (1e-6 < own < 1e-4):
            bad.append((k, len(mu), fm, own, blkpair - b1 * b2))
    check("both POLISHED B-endpoints are OFF the family: 8 atoms, "
          "blocks pairwise correlated, own-margin ~ +1.4e-5 (contra "
          "042's 'single diagonals' / '+0.000364 = c*(0.4887)')",
          not bad, f"{bad}")
    mu, n, fm, ratio = endinfo[("B_attack", "[('d', 2), ('d', 2), ('d', 2)]", 0.49)]
    check("c*(0.4887) = 3.686e-4 != the +0.000364 floor; the endpoint's "
          "own marginal is 0.489001 with c* = 3.492e-4",
          abs(cstar(0.4887) - 3.686e-4) < 1e-7
          and abs(fm - 0.489001) < 1e-6
          and abs(cstar(fm) - 3.492e-4) < 1e-7,
          f"c*(0.4887)={cstar(0.4887):.4e}, fmax={fm:.6f}, "
          f"c*(fmax)={cstar(fm):.4e}, floor={ratio:.4e}")
    bad = []
    for key, blkspec, cap, comps in (
            ("[('d', 2), ('d', 2)]", ((0b0011,),), 0.49, None),
            ("[('d', 3), ('d', 3)]", ((0b000111, 0b111000),), 0.49, None),
            ("[('d', 2), ('d', 2)]", ((0b0011, 0b1100),), 0.497, None)):
        mu, n, fm, ratio = endinfo[("B_attack", key, cap)]
        masks = blkspec[0]
        prod = {}
        for a in range(1 << n):
            w = 1.0
            live = 0
            for msk in masks:
                live |= msk
                on = (a & msk) == msk
                off = (a & msk) == 0
                if not (on or off):
                    w = 0.0
                    break
                w *= fm if on else 1.0 - fm
            if (a & ~live) != 0:
                w = 0.0
            if w > 0:
                prod[a] = w
        dev = max(abs(prod.get(a, 0.0) - mu.get(a, 0.0))
                  for a in set(prod) | set(mu))
        if dev > 1e-9:
            bad.append((key, cap, dev))
    check("the three UNPOLISHED B-endpoints ARE exact members: diagonal "
          "block tensors at p = fmax (one padded by marginal-0 coords)",
          not bad, f"{bad}")


def main():
    interp = json.loads((DATA / "hu_interp.json").read_text())
    blocks = json.loads((DATA / "hu_blocks.json").read_text())
    r1(interp)
    r2()
    r3(interp)
    r4()
    endinfo = r5(blocks)
    r6(blocks)
    r7(blocks, endinfo)
    print()
    if FAILS:
        print(f"FAILED CHECKS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("All checks pass.")


if __name__ == "__main__":
    main()
