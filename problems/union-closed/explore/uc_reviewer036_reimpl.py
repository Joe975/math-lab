#!/usr/bin/env python3
"""Attempt 036 (reviewer pass on 035): independent re-implementation.

Nothing here imports from the 031/033/035 evaluator stack: the HU
coupling, the chain-rule value CR, and the canonical (greedy
most-predictable-first) order rule are rebuilt from the records' prose
alone (031 "The coupling", 009 definitions, 035 "Approach").  Weights
and z-values are exact Fractions; every entropy is evaluated in 60-digit
Decimal arithmetic (Decimal.ln), so agreement with the committed floats
is a genuine third-path check, and ties in the order rule are detected
two ways -- exact conditional-profile equality (sufficient), and a
1e-40 Decimal gap (far below 035's TIE_TOL = 1e-12).

The one deliberate exception: check R4 re-derives the exact enclosure
endpoints with 035's own kit code (uc_hu_certify) in order to test
CONTAINMENT of this file's independent 60-digit values -- the value
being checked is ours, the enclosure being checked is theirs.

Checks (all against committed checkpoints / record prose):
  R1  the 033 witness: canonical order, CR/H, 24-order enumeration
  R2  the n=4 cap-0.49 kill: order, CR, 16/24 negative, tie-free trace
  R3  the n=5 cap-0.49 kill: order, CR, 60/120 negative, canonical is
      worst; the step-0 tie is EXACT (two zero-marginal coordinates) and
      CR is invariant under their swap and under all 120 relabels
  R4  the three certificates: 60-digit values inside both kits' exact
      enclosures, on the same limit_denominator(1e7) rationalizations;
      canonical order of the RATIONALIZED measure re-derived (the
      certify script derives it from the float measure)
  R5  all 51 descent-endpoint rows: canonical order + floor value +
      global floor + violation lists + product-extremal constants
  R6  equivariance on all 24 relabels of the tie-free n=4 kill
  R7  sweep worst-seed values as recorded (C blocks)

Usage: python uc_reviewer036_reimpl.py       (exit 0 iff all checks pass)
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
getcontext().prec = 60
LN2 = Decimal(2).ln()
HALF = Fraction(1, 2)

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# ---------------------------------------------------------------- entropies
def dlog2(q: Fraction) -> Decimal:
    return (Decimal(q.numerator).ln() - Decimal(q.denominator).ln()) / LN2


def hbits(q: Fraction) -> Decimal:
    """Binary entropy of a Fraction, in bits, 60-digit Decimal."""
    if q == 0 or q == 1:
        return Decimal(0)
    return -(Decimal(float(0)) + _xlog2(q) + _xlog2(1 - q))


def _xlog2(q: Fraction) -> Decimal:
    if q == 0:
        return Decimal(0)
    return (Decimal(q.numerator) / Decimal(q.denominator)) * dlog2(q)


def H_of(mu) -> Decimal:
    return -sum((_xlog2(w) for w in mu.values() if w > 0), Decimal(0))


# ------------------------------------------------------- the HU coupling/CR
def hu_cr(n, mu, seq):
    """CR(mu, HU) in bits under revelation sequence seq (seq[k] = the
    coordinate revealed at step k).  mu: dict int->Fraction, normalized.
    Own construction: cells keyed by the two revealed-prefix value
    tuples; per cell x = mu(A_i=0 | a), y = mu(B_i=0 | b),
    z = min(max(1/2, x+y-1), x, y); cell splits (z, x-z, y-z, 1-x-y+z).
    Returns (CR, H, sum_Ehz) as Decimals."""
    atoms = list(mu.items())
    # conditional tables per step: values-on-seq[:k] -> (mass, zero-mass of seq[k])
    tabs = []
    for k in range(n):
        d = {}
        for a, w in atoms:
            key = tuple((a >> seq[j]) & 1 for j in range(k))
            e = d.setdefault(key, [Fraction(0), Fraction(0)])
            e[0] += w
            if not (a >> seq[k]) & 1:
                e[1] += w
        tabs.append(d)
    cells = {((), ()): Fraction(1)}
    Ehz = Decimal(0)
    for k in range(n):
        nxt = {}
        for (pa, pb), w in cells.items():
            ma, za = tabs[k][pa]
            mb, zb = tabs[k][pb]
            x, y = za / ma, zb / mb
            z = min(max(HALF, x + y - 1), x, y)
            Ehz += (Decimal(w.numerator) / Decimal(w.denominator)) * hb(z)
            for va, vb, cw in ((0, 0, z), (0, 1, x - z),
                               (1, 0, y - z), (1, 1, 1 - x - y + z)):
                if cw > 0:
                    key = (pa + (va,), pb + (vb,))
                    nxt[key] = nxt.get(key, Fraction(0)) + w * cw
        cells = nxt
    H = H_of(mu)
    return Ehz - H, H, Ehz


_hb_cache = {}


def hb(z: Fraction) -> Decimal:
    if z in _hb_cache:
        return _hb_cache[z]
    v = Decimal(0) if z in (0, 1) else -(_xlog2(z) + _xlog2(1 - z))
    _hb_cache[z] = v
    return v


# -------------------------------------------------------- the order rule
TIE_DEC = Decimal("1e-40")


def cond_profile(n, mu, i, S):
    """Multiset of (group mass, zero-mass) pairs for H(A_i | A_S)."""
    d = {}
    for a, w in mu.items():
        key = tuple((a >> j) & 1 for j in S)
        e = d.setdefault(key, [Fraction(0), Fraction(0)])
        e[0] += w
        if not (a >> i) & 1:
            e[1] += w
    return sorted((m, z) for m, z in d.values())


def cond_H_dec(profile) -> Decimal:
    tot = Decimal(0)
    for m, z in profile:
        if m > 0:
            tot += (Decimal(m.numerator) / Decimal(m.denominator)) * hb(z / m)
    return tot


def canon_seq(n, mu):
    """Own greedy rule.  Returns (seq, tie_trace) where tie_trace[k] is
    the number of coordinates tied at the step-k minimum (exact-profile
    OR <= 1e-40 Decimal gap).  Ties break lowest index."""
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    seq, trace = [], []
    remaining = sorted(set(range(n)))
    while remaining:
        profs = {i: cond_profile(n, mu, i, seq) for i in remaining}
        vals = {i: cond_H_dec(profs[i]) for i in remaining}
        lo = min(vals.values())
        tied = [i for i in remaining
                if vals[i] - lo <= TIE_DEC or profs[i] == profs[min(
                    j for j in remaining if vals[j] == lo)]]
        trace.append(len(tied))
        pick = min(tied)
        seq.append(pick)
        remaining.remove(pick)
    return tuple(seq), trace


def perm_to_seq(perm):
    """035 stores perm[coord] = slot; convert to seq[slot] = coord."""
    n = len(perm)
    seq = [None] * n
    for coord, slot in enumerate(perm):
        seq[slot] = coord
    return tuple(seq)


def load_mu(d):
    return {int(s, 2): Fraction(w) for s, w in d.items()}


def normalized(mu):
    tot = sum(mu.values())
    return {a: w / tot for a, w in mu.items()}


def max_marginal(n, mu):
    return max(sum(w for a, w in mu.items() if (a >> i) & 1) for i in range(n))


def enumerate_orders(n, mu):
    """CR under every revelation sequence; float-precision via Decimal."""
    out = {}
    for seq in itertools.permutations(range(n)):
        cr, H, _ = hu_cr(n, mu, seq)
        out[seq] = (cr, H)
    return out


def relabel(n, mu, sigma):
    """sigma[i] = new position of coordinate i."""
    out = {}
    for a, w in mu.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << sigma[i]
        out[b] = out.get(b, Fraction(0)) + w
    return out


def main():
    can = json.loads((DATA / "hu_canon.json").read_text())
    cert = json.loads((DATA / "hu_canon_certify.json").read_text())
    att = json.loads((DATA / "hu_attack.json").read_text())

    # ---------------------------------------------------------------- R1
    print("R1. The 033 witness under the canonical rule:")
    wrec = can["A_witness"][0]
    w033 = att["cap_0.45"]["violations"][0]
    n = w033["n"]
    mu = normalized(load_mu(w033["mu"]))
    seq, trace = canon_seq(n, mu)
    check("canonical order matches",
          seq == perm_to_seq(wrec["canon_order"]),
          f"mine {seq}, 035 {perm_to_seq(wrec['canon_order'])}, ties {trace}")
    cr, H, _ = hu_cr(n, mu, seq)
    ratio = float(cr / H)
    check("canonical CR/H reproduces", abs(ratio - wrec["canon_ratio"]) < 1e-9,
          f"{ratio:+.9f} vs {wrec['canon_ratio']:+.9f}")
    allv = enumerate_orders(n, mu)
    ratios = sorted(float(c / h) for c, h in allv.values())
    negs = sum(1 for r in ratios if r < 0)
    # the checkpoint's canon_rank is 1-based (engine stores rank+1)
    rank = 1 + sum(1 for r in ratios if r < ratio - 1e-12)
    check("24-order enumeration: 1 negative, canonical rank 23 (1-based), "
          "worst/best match",
          negs == 1 and rank == wrec["canon_rank"]
          and abs(ratios[0] - wrec["worst_ratio"]) < 1e-9
          and abs(ratios[-1] - wrec["best_ratio"]) < 1e-9,
          f"negs {negs}, rank {rank}/{len(ratios)}, worst {ratios[0]:+.6f}, "
          f"best {ratios[-1]:+.6f}")

    # ---------------------------------------------------------------- R2
    print("R2. The n=4 cap-0.49 kill:")
    v4 = can["B_cap_0.49"]["violations"][0]
    n4, mu4 = v4["n"], normalized(load_mu(v4["mu"]))
    seq4, trace4 = canon_seq(n4, mu4)
    check("canonical order matches", seq4 == perm_to_seq(v4["canon_order"]),
          f"mine {seq4}, ties {trace4}")
    check("tie trace is [1,1,1,1] (tie-free at 1e-40 resolution)",
          trace4 == [1, 1, 1, 1], f"{trace4}")
    cr4, H4, _ = hu_cr(n4, mu4, seq4)
    check("CR/H equals recorded floor",
          abs(float(cr4 / H4) - v4["floor"]) < 1e-9,
          f"{float(cr4 / H4):+.9f} vs {v4['floor']:+.9f}")
    o4 = enumerate_orders(n4, mu4)
    negs4 = sum(1 for c, _ in o4.values() if c < 0)
    best4 = max(float(c) for c, _ in o4.values())
    check("record part D: 16 of 24 orders negative, best +0.1507",
          negs4 == 16 and abs(best4 - 0.1507) < 5e-5,
          f"negs {negs4}/24, best CR {best4:+.5f}")

    # ---------------------------------------------------------------- R3
    print("R3. The n=5 cap-0.49 kill:")
    v5 = can["B_cap_0.49"]["violations"][1]
    n5, mu5 = v5["n"], normalized(load_mu(v5["mu"]))
    seq5, trace5 = canon_seq(n5, mu5)
    check("canonical order matches", seq5 == perm_to_seq(v5["canon_order"]),
          f"mine {seq5}, ties {trace5}")
    marg = [float(sum(w for a, w in mu5.items() if (a >> i) & 1))
            for i in range(n5)]
    zeros = [i for i in range(n5) if marg[i] == 0.0]
    check("step-0 tie is EXACT: two coordinates with marginal exactly 0",
          trace5[0] == 2 and zeros == [1, 3],
          f"tie trace {trace5}, zero-marginal coords {zeros}")
    cr5, H5, _ = hu_cr(n5, mu5, seq5)
    check("CR/H equals recorded floor",
          abs(float(cr5 / H5) - v5["floor"]) < 1e-9,
          f"{float(cr5 / H5):+.9f} vs {v5['floor']:+.9f}")
    o5 = enumerate_orders(n5, mu5)
    vals5 = sorted(float(c) for c, _ in o5.values())
    negs5 = sum(1 for v in vals5 if v < 0)
    check("record part D: 60 of 120 negative, canonical IS worst, "
          "best +0.0350",
          negs5 == 60 and abs(vals5[0] - float(cr5)) < 1e-12
          and abs(vals5[-1] - 0.0350) < 5e-5,
          f"negs {negs5}/120, worst {vals5[0]:+.6f} (=canonical "
          f"{float(cr5):+.6f}), best {vals5[-1]:+.5f}")
    # the tie is between two DETERMINISTIC coordinates: swapping them in
    # the revelation sequence cannot change the coupling
    seq5b = tuple(3 if c == 1 else 1 if c == 3 else c for c in seq5)
    cr5b, _, _ = hu_cr(n5, mu5, seq5b)
    check("CR invariant under swapping the tied (deterministic) pair "
          "in the order", abs(cr5b - cr5) < Decimal("1e-50"),
          f"delta {float(cr5b - cr5):.1e}")
    # full labelling-independence: canonical CR identical on all 120 relabels
    worst_dev = Decimal(0)
    for sigma in itertools.permutations(range(n5)):
        mus = normalized(relabel(n5, mu5, sigma))
        sq, _ = canon_seq(n5, mus)
        crs, _, _ = hu_cr(n5, mus, sq)
        worst_dev = max(worst_dev, abs(crs - cr5))
    check("canonical CR identical on ALL 120 relabels -> the n=5 kill "
          "is labelling-independent in fact (stronger than 035 claims)",
          worst_dev < Decimal("1e-50"), f"max |dev| {float(worst_dev):.1e}")

    # ---------------------------------------------------------------- R4
    print("R4. The three certificates (rationalized measures, 60-digit "
          "values vs exact kit enclosures):")
    sys.path.insert(0, str(HERE))
    from uc_hu_certify import hu_cells, ivl_add, ivl_scale, xlog2x_ivl, \
        h_ivl, log2_A, log2_B          # kit code: containment leg only

    def rationalize(mu_float):
        m = {a: Fraction(w).limit_denominator(10 ** 7)
             for a, w in mu_float.items()}
        tot = sum(m.values())
        return {a: w / tot for a, w in m.items()}

    for crow, (src, cap) in zip(cert["witnesses"],
                                [(v4, Fraction(49, 100)),
                                 (v5, Fraction(49, 100)),
                                 (w033, Fraction(9, 20))]):
        tag = crow["witness"]
        n_ = src["n"]
        muQ = rationalize({int(s, 2): w for s, w in src["mu"].items()})
        seqQ, traceQ = canon_seq(n_, muQ)
        check(f"{tag}: canonical order of the RATIONALIZED measure equals "
              "the certified order",
              seqQ == perm_to_seq(crow["canon_order"]),
              f"mine {seqQ}, certified {perm_to_seq(crow['canon_order'])}, "
              f"ties {traceQ}")
        crQ, HQ, _ = hu_cr(n_, muQ, seqQ)
        margQ = max_marginal(n_, muQ)
        check(f"{tag}: exact marginal < cap and H > 1/2",
              margQ < cap and HQ > Decimal("0.5"),
              f"marg {float(margQ):.6f} vs cap {float(cap)}, H {float(HQ):.4f}")
        # rebuild each kit's exact enclosure and test containment of ours
        perm = crow["canon_order"]
        muP = {}
        for a, w in muQ.items():
            b = 0
            for i in range(n_):
                if (a >> i) & 1:
                    b |= 1 << perm[i]
            muP[b] = muP.get(b, Fraction(0)) + w
        zs = hu_cells(n_, muP)
        for kitname, fn in (("A_digit_extraction", log2_A),
                            ("B_atanh_series", log2_B)):
            Hiv = (Fraction(0), Fraction(0))
            for w in muP.values():
                Hiv = ivl_add(Hiv, ivl_scale(Fraction(-1), xlog2x_ivl(w, fn)))
            S = (Fraction(0), Fraction(0))
            for w, z in zs:
                if 0 < z < 1:
                    S = ivl_add(S, ivl_scale(w, h_ivl(z, fn)))
            lo, hi = S[0] - Hiv[1], S[1] - Hiv[0]
            mine = Fraction(str(crQ))
            inside = lo <= mine <= hi
            sign_ok = (hi < 0) if crow["want"] == "neg" else (lo > 0)
            check(f"{tag} [{kitname}]: my 60-digit CR inside the exact "
                  "enclosure; sign certifies",
                  inside and sign_ok,
                  f"CR {float(mine):+.9e} in [{float(lo):+.9e}, "
                  f"{float(hi):+.9e}]")

    # ---------------------------------------------------------------- R5
    print("R5. All 51 descent-endpoint rows (order + floor + summary):")
    hb_cap = lambda p: -(p * math.log2(p) + (1 - p) * math.log2(1 - p))
    for capkey, capv in (("B_cap_0.38271", 0.38271), ("B_cap_0.45", 0.45),
                         ("B_cap_0.49", 0.49)):
        blk = can[capkey]
        bad_order, bad_val = [], []
        floors = []
        for row in blk["rows"]:
            n_ = row["n"]
            m = normalized(load_mu(row["mu"]))
            sq, _ = canon_seq(n_, m)
            if sq != perm_to_seq(row["canon_order"]):
                bad_order.append(row["start"])
            cr_, H_, _ = hu_cr(n_, m, sq)
            val = float(cr_ / H_)
            floors.append(val)
            if abs(val - row["floor"]) > 1e-9:
                bad_val.append((row["start"], val, row["floor"]))
        check(f"{capkey}: every row's canonical order matches",
              not bad_order, f"mismatches: {bad_order}")
        check(f"{capkey}: every row's floor value reproduces (<=1e-9)",
              not bad_val, f"mismatches: {bad_val[:2]}")
        check(f"{capkey}: global floor = min over rows; violations = "
              "negative rows",
              abs(min(floors) - blk["global_floor"]) < 1e-12
              and sum(1 for f in floors if f < 0) == len(blk["violations"]),
              f"min {min(floors):+.9f} vs {blk['global_floor']:+.9f}, "
              f"negs {sum(1 for f in floors if f < 0)}")
        pe = (1 - hb_cap(capv)) / hb_cap(capv)
        check(f"{capkey}: product-extremal constant",
              abs(pe - blk["product_extremal"]) < 1e-12, f"{pe:.9f}")

    # ---------------------------------------------------------------- R6
    print("R6. Equivariance on all 24 relabels of the tie-free n=4 kill:")
    bad = 0
    for sigma in itertools.permutations(range(n4)):
        mus = normalized(relabel(n4, mu4, sigma))
        sq, tr = canon_seq(n4, mus)
        expect = tuple(sigma[c] for c in seq4)
        if sq != expect or tr != [1, 1, 1, 1]:
            bad += 1
    check("canonical order commutes with every relabel, all tie-free",
          bad == 0, f"failures {bad}/24")

    # ---------------------------------------------------------------- R7
    print("R7. Sweep-block summaries as recorded:")
    for capkey, wc, mg in (("C_cap_0.45", 0.040165171987841884,
                            0.07623552767743069),
                           ("C_cap_0.49", 0.019985769184519192,
                            0.07939672024859677)):
        c = can[capkey]
        check(f"{capkey}: 120 tested, 0 canonical / 0 any-order violations, "
              "worst + mean-gap match record",
              c["tested"] == 120 and c["canon_violations"] == 0
              and c["any_order_violations"] == 0
              and abs(c["worst_canon"] - wc) < 1e-15
              and abs(c["mean_gap"] - mg) < 1e-15,
              f"worst {c['worst_canon']:+.6f}, gap {c['mean_gap']:+.5f}")

    print()
    if FAILS:
        print(f"REFUTATIONS/FAILURES: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("All checks pass; no refutation found.")


if __name__ == "__main__":
    main()
