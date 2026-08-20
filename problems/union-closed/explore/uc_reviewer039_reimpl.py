#!/usr/bin/env python3
"""Attempt 039 (reviewer pass on 036/037/038): independent re-implementation.

Nothing in the value-producing path imports from the 031/033/035/037
evaluator stack or from uc_reviewer036_reimpl: the HU coupling, CR, and
all four order rules (can / canst / surp / roll) are rebuilt from the
records' prose alone (031 "The coupling", 009 CR definition, 035
"Approach", 037 "Approach").  Structure is deliberately different from
all three committed implementations: CR is a memoized depth-first
RECURSION over history trees keyed by atom-index subsets (the engine
uses iterative value-tuple cell dicts, the 037 skeptic iterative
index-partition lists, 036 iterative Fraction cell dicts); entropies
are accumulated in NATS via math.log and converted to bits once at the
end; the certificate leg re-runs the same recursion with exact Fraction
weights and 60-digit Decimal logs.

Two deliberate exceptions, both non-value-producing:
  * check V3 imports 035/037's own kit code (uc_hu_certify) to rebuild
    the exact enclosures whose CONTAINMENT of this file's independent
    60-digit values is under test (the 036 R4 precedent) -- the value
    is ours, the enclosure is theirs;
  * check V5 property-tests the COMMITTED rule functions themselves
    (uc_hu_canon.cond_entropy, uc_hu_order2.canon_completion/seq_can/
    seq_roll), because what is audited there is the committed code's
    own consistency (the (ROLL-DOM) recursion hinge), not a value.
The census generator gen_instance is replicated verbatim (fixed-seed
instance IDENTITY is shared by design; every VALUE is re-derived here).

Checks:
  V1  037 part A: the three record witnesses -- all four rules'
      sequences + CR, full 24/120-order enumerations, worst/best,
      record-prose key values (roll = best on all three; surp = worst
      of 120 on kill_n5)
  V2  every committed descent endpoint: hu_order2.json C_*/C2_*
      (own rule + own evaluator), D_* (own full enumeration),
      B_endpoints from the 51 committed 035 rows; hu_rollcensus.json
      P2_*/P3_*; global floors = min over rows, violation lists
  V3  the five certificates: 60-digit CR on the same
      limit_denominator(1e7) rationalizations inside both kits' exact
      enclosures, sign positive, exact side conditions; rollout order
      of each rationalized measure re-derived at 60 digits WITH
      60-digit canonical completions (stronger than the committed
      fixed-point check, whose completions are float-derived)
  V4  census: 12 spot rows per P1 block (every 25th) regenerated and
      fully re-scored (roll/canon CR, best/worst, H, ranks, flags);
      all summary fields recomputed from stored rows; canon_rank <=
      roll_rank on all 1200 rows ((ROLL-DOM) numerically)
  V5  (ROLL-DOM) hinges on the committed code: cond_entropy depends on
      the revealed SET only; seq_can == canon_completion from the
      empty prefix; the completion recursion identity; and a float
      tie-tolerance stress test (perturbed exchangeable measures,
      CR_roll >= CR_canon - n*TIE - float slack)
  V6  the 0.499 fixed-point story: the four step-0 rollout scores of
      mu_Q are EXACTLY equal at 60 digits and the unrationalized float
      measure's step-0 gap is ~1.6e-10 with roll(mu_float) !=
      roll(mu_Q) -- 037 part E verified -- but mu_Q is NOT
      exactly-exchangeable as 037 lead 4 claims: the 8.1e-11 atom
      rationalizes to 0 and the 2-atom remainder has stabilizer S3 on
      {0,2,3} only (a CORRECTION, codified here)
  V7  reporting audit: every number quoted in the 037/038 record prose
      re-checked against the committed checkpoints, plus the corrected
      count of the part-C roll-skipped starts (two n=6 and ONE n=8 --
      not "three n=6" -- and P2 covers only the two n=6 ones)

Usage: python uc_reviewer039_reimpl.py       (exit 0 iff all checks pass)
"""
from __future__ import annotations

import itertools
import json
import math
import random
import sys
from decimal import Decimal, getcontext
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
getcontext().prec = 60
LN2 = math.log(2.0)
DLN2 = Decimal(2).ln()
HALF = Fraction(1, 2)

FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}"
          + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


# ===================================================== float path (nats)
def hn(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -(p * math.log(p) + (1.0 - p) * math.log(1.0 - p))


def norm_atoms(mu):
    tot = sum(mu.values())
    return [(a, w / tot) for a, w in sorted(mu.items())]


def hu_cr(n, atoms, seq):
    """(CR, H) in bits.  Memoized DFS over history trees: state =
    (depth, alive A-atom indices, alive B-atom indices); the recursion
    returns the conditional expected future sum of h(z) in nats per
    unit cell mass."""
    memo = {}
    masses = [w for _, w in atoms]

    def cond(k, ia, ib):
        if k == n:
            return 0.0
        key = (k, ia, ib)
        v = memo.get(key)
        if v is not None:
            return v
        c = seq[k]
        ma = sum(masses[t] for t in ia)
        mb = sum(masses[t] for t in ib)
        ia0 = tuple(t for t in ia if not (atoms[t][0] >> c) & 1)
        ib0 = tuple(t for t in ib if not (atoms[t][0] >> c) & 1)
        x = sum(masses[t] for t in ia0) / ma
        y = sum(masses[t] for t in ib0) / mb
        z = min(max(0.5, x + y - 1.0), x, y)
        val = hn(z)
        ia1 = tuple(t for t in ia if (atoms[t][0] >> c) & 1)
        ib1 = tuple(t for t in ib if (atoms[t][0] >> c) & 1)
        for A, B, cw in ((ia0, ib0, z), (ia0, ib1, x - z),
                         (ia1, ib0, y - z), (ia1, ib1, 1.0 - x - y + z)):
            if cw > 0.0 and A and B:
                val += cw * cond(k + 1, A, B)
        memo[key] = val
        return val

    all_idx = tuple(range(len(atoms)))
    Ehz = cond(0, all_idx, all_idx)
    Hn = -sum(w * math.log(w) for w in masses if w > 0.0)
    return (Ehz - Hn) / LN2, Hn / LN2


def cells_after(n, atoms, prefix):
    """[(ia, ib, w)] coupling cells after revealing prefix (DFS)."""
    masses = [w for _, w in atoms]
    out = []

    def rec(k, ia, ib, w):
        if k == len(prefix):
            out.append((ia, ib, w))
            return
        c = prefix[k]
        ma = sum(masses[t] for t in ia)
        mb = sum(masses[t] for t in ib)
        ia0 = tuple(t for t in ia if not (atoms[t][0] >> c) & 1)
        ib0 = tuple(t for t in ib if not (atoms[t][0] >> c) & 1)
        x = sum(masses[t] for t in ia0) / ma
        y = sum(masses[t] for t in ib0) / mb
        z = min(max(0.5, x + y - 1.0), x, y)
        ia1 = tuple(t for t in ia if (atoms[t][0] >> c) & 1)
        ib1 = tuple(t for t in ib if (atoms[t][0] >> c) & 1)
        for A, B, cw in ((ia0, ib0, z), (ia0, ib1, x - z),
                         (ia1, ib0, y - z), (ia1, ib1, 1.0 - x - y + z)):
            if cw > 0.0 and A and B:
                rec(k + 1, A, B, w * cw)

    all_idx = tuple(range(len(atoms)))
    rec(0, all_idx, all_idx, 1.0)
    return out


def condH_bits(n, atoms, i, S):
    """H(A_i | A_S) in bits; groups by the integer projection onto S."""
    groups = {}
    for a, w in atoms:
        p = 0
        for t, j in enumerate(S):
            if (a >> j) & 1:
                p |= 1 << t
        g = groups.get(p)
        if g is None:
            groups[p] = g = [0.0, 0.0]
        g[0] += w
        if not (a >> i) & 1:
            g[1] += w
    return sum(m * hn(zm / m) for m, zm in groups.values() if m > 0.0) / LN2


TIE = 1e-12  # bits; the records' stated tolerance


def rule_can(n, atoms, prefix=()):
    """035's greedy min-conditional-entropy rule from a prefix; ties
    within 1e-12 -> lowest index."""
    seq = list(prefix)
    rem = [i for i in range(n) if i not in seq]
    while rem:
        vals = [(condH_bits(n, atoms, i, seq), i) for i in rem]
        lo = min(v for v, _ in vals)
        pick = min(i for v, i in vals if v <= lo + TIE)
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq)


def step_surplus(n, atoms, cells, i):
    """Sigma_cells w*(h(z) - (h(x)+h(y))/2), bits."""
    masses = [w for _, w in atoms]
    s = 0.0
    for ia, ib, w in cells:
        ma = sum(masses[t] for t in ia)
        mb = sum(masses[t] for t in ib)
        x = sum(masses[t] for t in ia if not (atoms[t][0] >> i) & 1) / ma
        y = sum(masses[t] for t in ib if not (atoms[t][0] >> i) & 1) / mb
        z = min(max(0.5, x + y - 1.0), x, y)
        s += w * (hn(z) - 0.5 * (hn(x) + hn(y)))
    return s / LN2


def rule_surp(n, atoms):
    seq, rem = [], list(range(n))
    while rem:
        cells = cells_after(n, atoms, seq)
        vals = [(step_surplus(n, atoms, cells, i), i) for i in rem]
        hi = max(v for v, _ in vals)
        pick = min(i for v, i in vals if v >= hi - TIE)
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq)


def rule_canst(n, atoms):
    seq, rem = [], list(range(n))
    while rem:
        vals = [(condH_bits(n, atoms, i, seq), i) for i in rem]
        lo = min(v for v, _ in vals)
        tied = [i for v, i in vals if v <= lo + TIE]
        if len(tied) > 1:
            cells = cells_after(n, atoms, seq)
            sc = [(step_surplus(n, atoms, cells, i), i) for i in tied]
            hi = max(v for v, _ in sc)
            pick = min(i for v, i in sc if v >= hi - TIE)
        else:
            pick = tied[0]
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq)


def rule_roll(n, atoms):
    """Rollout: maximise full CR with canonical completion; ties within
    1e-12 -> lowest index."""
    seq, rem = [], list(range(n))
    while rem:
        scores = {}
        for i in rem:
            full = rule_can(n, atoms, tuple(seq) + (i,))
            scores[i], _ = hu_cr(n, atoms, full)
        hi = max(scores.values())
        pick = min(i for i in rem if scores[i] >= hi - TIE)
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq)


MY_RULES = {"can": lambda n, at: rule_can(n, at), "canst": rule_canst,
            "surp": rule_surp, "roll": rule_roll}


def ratio_under_rule(n, mu, cap, rulename):
    """CR/H with the records' regime filters; None if out of regime."""
    tot = sum(mu.values())
    m = {a: w / tot for a, w in mu.items() if w / tot > 1e-12}
    fmax = max(sum(w for a, w in m.items() if (a >> i) & 1)
               for i in range(n))
    if fmax >= cap:
        return None
    atoms = norm_atoms(m)
    cr, H = hu_cr(n, atoms, MY_RULES[rulename](n, atoms))
    if H < 0.2:
        return None
    return cr / H


def load_mu(d):
    return {int(s, 2): w for s, w in d.items()}


# ================================================ 60-digit Decimal path
_dln_cache = {}


def dln(q: Fraction) -> Decimal:
    v = _dln_cache.get(q)
    if v is None:
        v = Decimal(q.numerator).ln() - Decimal(q.denominator).ln()
        _dln_cache[q] = v
    return v


def hd(q: Fraction) -> Decimal:
    """Binary entropy of a Fraction in NATS, 60-digit Decimal."""
    if q == 0 or q == 1:
        return Decimal(0)
    qd = Decimal(q.numerator) / Decimal(q.denominator)
    rd = Decimal((1 - q).numerator) / Decimal((1 - q).denominator)
    return -(qd * dln(q) + rd * dln(1 - q))


def frac_dec(q: Fraction) -> Decimal:
    return Decimal(q.numerator) / Decimal(q.denominator)


def hu_cr_dec(n, muQ, seq):
    """(CR, H) in bits as 60-digit Decimals; exact Fraction cells, same
    DFS structure as hu_cr."""
    atoms = sorted(muQ.items())
    masses = [w for _, w in atoms]
    memo = {}

    def cond(k, ia, ib):
        if k == n:
            return Decimal(0)
        key = (k, ia, ib)
        v = memo.get(key)
        if v is not None:
            return v
        c = seq[k]
        ma = sum(masses[t] for t in ia)
        mb = sum(masses[t] for t in ib)
        ia0 = tuple(t for t in ia if not (atoms[t][0] >> c) & 1)
        ib0 = tuple(t for t in ib if not (atoms[t][0] >> c) & 1)
        x = sum(masses[t] for t in ia0) / ma
        y = sum(masses[t] for t in ib0) / mb
        z = min(max(HALF, x + y - 1), x, y)
        val = hd(z)
        ia1 = tuple(t for t in ia if (atoms[t][0] >> c) & 1)
        ib1 = tuple(t for t in ib if (atoms[t][0] >> c) & 1)
        for A, B, cw in ((ia0, ib0, z), (ia0, ib1, x - z),
                         (ia1, ib0, y - z), (ia1, ib1, 1 - x - y + z)):
            if cw > 0 and A and B:
                val += frac_dec(cw) * cond(k + 1, A, B)
        memo[key] = val
        return val

    all_idx = tuple(range(len(atoms)))
    Ehz = cond(0, all_idx, all_idx)
    Hn = -sum((frac_dec(w) * dln(w) for w in masses if w > 0), Decimal(0))
    return (Ehz - Hn) / DLN2, Hn / DLN2


def cond_profile(n, muQ, i, S):
    d = {}
    for a, w in muQ.items():
        p = 0
        for t, j in enumerate(S):
            if (a >> j) & 1:
                p |= 1 << t
        e = d.setdefault(p, [Fraction(0), Fraction(0)])
        e[0] += w
        if not (a >> i) & 1:
            e[1] += w
    return sorted((m, zm) for m, zm in d.values())


def rule_can_dec(n, muQ, prefix=()):
    """Canonical rule on an exact-Fraction measure: 60-digit conditional
    entropies, ties by exact profile equality or a <=1e-40 Decimal gap,
    broken lowest-index."""
    seq = list(prefix)
    rem = [i for i in range(n) if i not in seq]
    while rem:
        profs = {i: cond_profile(n, muQ, i, seq) for i in rem}
        vals = {i: sum((frac_dec(m) * hd(zm / m) for m, zm in profs[i]
                        if m > 0), Decimal(0)) for i in rem}
        lo = min(vals.values())
        lo_i = min(i for i in rem if vals[i] == lo)
        tied = [i for i in rem if vals[i] - lo <= Decimal("1e-40")
                or profs[i] == profs[lo_i]]
        pick = min(tied)
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq)


def rule_roll_dec(n, muQ):
    """Rollout on an exact-Fraction measure, scores AND canonical
    completions at 60 digits (stronger than the committed fixed-point
    check, which derives completions from the float projection).
    Returns (seq, per-step score gap diagnostics)."""
    seq, rem = [], list(range(n))
    gaps = []
    while rem:
        scores = {}
        for i in rem:
            full = rule_can_dec(n, muQ, tuple(seq) + (i,))
            scores[i], _ = hu_cr_dec(n, muQ, full)
        hi = max(scores.values())
        gaps.append(float(hi - min(scores.values())))
        pick = min(i for i in rem if hi - scores[i] <= Decimal("1e-12"))
        seq.append(pick)
        rem.remove(pick)
    return tuple(seq), gaps


def rationalize(mu_float):
    m = {a: Fraction(w).limit_denominator(10 ** 7)
         for a, w in mu_float.items()}
    tot = sum(m.values())
    return {a: w / tot for a, w in m.items()}


def relabel_frac(n, muQ, sigma):
    out = {}
    for a, w in muQ.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << sigma[i]
        out[b] = out.get(b, Fraction(0)) + w
    return out


# ====================================================== census generator
def gen_instance(rng, n, cap):
    """Replicated verbatim from uc_hu_rollcensus.py (instance identity
    shared by design; every value re-derived independently)."""
    k = rng.randint(3, min(14, 1 << n))
    supp = rng.sample(range(1 << n), k)
    mu = {a: math.exp(rng.uniform(-2.5, 2.5)) for a in supp}
    tot = sum(mu.values())
    mu = {a: w / tot for a, w in mu.items()}
    fm = max(sum(w for a, w in mu.items() if (a >> i) & 1)
             for i in range(n))
    if fm >= cap - 0.005:
        al = 1.0 - (cap - 0.005) / fm
        mu = {a: w * (1 - al) for a, w in mu.items()}
        mu[0] = mu.get(0, 0.0) + al
    return mu


# =============================================================== checks
def v1(ord2, can, att):
    print("V1. 037 part A: the three witnesses under all four rules:")
    src = {"033witness": att["cap_0.45"]["violations"][0],
           "kill_n4": can["B_cap_0.49"]["violations"][0],
           "kill_n5": can["B_cap_0.49"]["violations"][1]}
    for row in ord2["A_witnesses"]:
        v = src[row["tag"]]
        n = v["n"]
        atoms = norm_atoms(load_mu(v["mu"]))
        allcr = sorted(hu_cr(n, atoms, s)[0]
                       for s in itertools.permutations(range(n)))
        check(f"{row['tag']}: worst/best/negative-count reproduce",
              abs(allcr[0] - row["worst_CR"]) < 1e-9
              and abs(allcr[-1] - row["best_CR"]) < 1e-9
              and sum(1 for c in allcr if c < 0) == row["n_negative"],
              f"worst {allcr[0]:+.6f}, best {allcr[-1]:+.6f}, "
              f"neg {sum(1 for c in allcr if c < 0)}/{len(allcr)}")
        for rn in ("can", "canst", "surp", "roll"):
            seq = MY_RULES[rn](n, atoms)
            cr, _ = hu_cr(n, atoms, seq)
            rank = 1 + sum(1 for c in allcr if c < cr - 1e-12)
            check(f"{row['tag']} {rn}: sequence, CR and rank reproduce",
                  list(seq) == row[rn]["seq"]
                  and abs(cr - row[rn]["CR"]) < 1e-9
                  and rank == row[rn]["rank"],
                  f"mine {seq} {cr:+.6f} r{rank}")
        # the record-prose key claims
        if row["tag"] == "kill_n5":
            surp_cr = hu_cr(n, atoms, MY_RULES["surp"](n, atoms))[0]
            check("record: surp picks the worst of 120 on kill_n5 "
                  "(-0.02889)",
                  abs(surp_cr - allcr[0]) < 1e-12
                  and abs(surp_cr - (-0.02889)) < 5e-6,
                  f"surp {surp_cr:+.6f} vs worst {allcr[0]:+.6f}")
        roll_cr = hu_cr(n, atoms, MY_RULES["roll"](n, atoms))[0]
        check(f"record: roll lands on the best order ({row['tag']})",
              abs(roll_cr - allcr[-1]) < 1e-12,
              f"roll {roll_cr:+.6f} = best {allcr[-1]:+.6f}")


def v2(ord2, can, roc):
    print("V2. Every committed descent endpoint, own rule + evaluator:")
    for key, blk in ord2.items():
        if key.startswith(("C_", "C2_")):
            rn = key.split("_")[1] if key.startswith("C_") else "roll"
            cap = float(key.rsplit("_", 1)[-1])
            bad = []
            for r in blk["rows"]:
                v = ratio_under_rule(r["n"], load_mu(r["mu"]), cap, rn)
                if v is None or abs(v - r["floor"]) > 1e-8:
                    bad.append((r["start"], v, r["floor"]))
            floors = [r["floor"] for r in blk["rows"]]
            check(f"{key}: every endpoint floor reproduces", not bad,
                  f"mismatches {bad[:2]}")
            check(f"{key}: global floor = min, violations = negatives",
                  abs(min(floors) - blk["global_floor"]) < 1e-12
                  and sum(1 for f in floors if f < 0)
                  == len(blk["violations"]))
        elif key.startswith("D_"):
            bad = []
            for r in blk["rows"]:
                n = r["n"]
                atoms = norm_atoms(load_mu(r["mu"]))
                best = max(hu_cr(n, atoms, s)[0]
                           for s in itertools.permutations(range(n)))
                H = hu_cr(n, atoms, tuple(range(n)))[1]
                if abs(best / H - r["floor"]) > 1e-8:
                    bad.append((r["start"], best / H, r["floor"]))
            floors = [r["floor"] for r in blk["rows"]]
            check(f"{key}: every best-order endpoint reproduces", not bad,
                  f"mismatches {bad[:2]}")
            check(f"{key}: global floor = min, violations = negatives",
                  abs(min(floors) - blk["global_floor"]) < 1e-12
                  and sum(1 for f in floors if f < 0)
                  == len(blk["violations"]))
    for capkey in ("B_cap_0.38271", "B_cap_0.45", "B_cap_0.49"):
        cap = float(capkey.split("_")[-1])
        summ = ord2["B_endpoints"][capkey]
        for rn in ("can", "canst", "surp", "roll"):
            worst, neg = 1e9, 0
            for row in can[capkey]["rows"]:
                v = ratio_under_rule(row["n"], load_mu(row["mu"]), cap, rn)
                if v is None:
                    continue
                worst = min(worst, v)
                neg += 1 if v < 0 else 0
            check(f"{capkey} {rn}: worst + negative count reproduce",
                  abs(worst - summ[rn]["worst"]) < 1e-8
                  and neg == summ[rn]["negative_rows"],
                  f"worst {worst:+.6f}, neg {neg}")
    for key, blk in roc.items():
        if not key.startswith(("P2_", "P3_")):
            continue
        cap = float(key.split("_")[-1])
        bad = []
        for r in blk["rows"]:
            v = ratio_under_rule(r["n"], load_mu(r["mu"]), cap, "roll")
            if v is None or abs(v - r["floor"]) > 1e-8:
                bad.append((r["start"], v, r["floor"]))
        floors = [r["floor"] for r in blk["rows"]]
        check(f"{key}: every endpoint floor reproduces (own rollout, "
              f"n={blk['rows'][0]['n']})", not bad, f"mismatches {bad[:2]}")
        check(f"{key}: global floor = min, violations = negatives",
              abs(min(floors) - blk["global_floor"]) < 1e-12
              and sum(1 for f in floors if f < 0) == len(blk["violations"]))


def v3(ord2, can, cert):
    print("V3. The five certificates (60-digit values, exact enclosures, "
          "rationalized rollout orders):")
    sys.path.insert(0, str(HERE))
    from uc_hu_certify import (hu_cells, ivl_add, ivl_scale, xlog2x_ivl,
                               h_ivl, log2_A, log2_B)  # enclosure leg only
    srcs = []
    for r in can["B_cap_0.49"]["violations"]:
        srcs.append((load_mu(r["mu"]), r["n"], Fraction(49, 100)))
    for key, capf in (("C_roll_cap_0.49", Fraction(49, 100)),
                      ("C_roll_cap_0.497", Fraction(497, 1000)),
                      ("C2_roll_cap_0.499", Fraction(499, 1000))):
        row = min(ord2[key]["rows"], key=lambda r: r["floor"])
        srcs.append((load_mu(row["mu"]), row["n"], capf))
    for crow, (muF, n, cap) in zip(cert["witnesses"], srcs):
        tag = crow["witness"]
        muQ = rationalize(muF)
        seq_mine, gaps = rule_roll_dec(n, muQ)
        check(f"{tag}: rollout order of the RATIONALIZED measure "
              "(60-digit scores AND completions) equals roll_seq",
              list(seq_mine) == crow["roll_seq"],
              f"mine {seq_mine}, certified {crow['roll_seq']}")
        crD, HD = hu_cr_dec(n, muQ, tuple(crow["roll_seq"]))
        margQ = max(sum(w for a, w in muQ.items() if (a >> i) & 1)
                    for i in range(n))
        check(f"{tag}: exact side conditions (marginal < cap, H > 1/2)",
              margQ < cap and HD > Decimal("0.5"),
              f"marg {float(margQ):.6f} < {float(cap)}, H {float(HD):.4f}")
        check(f"{tag}: checkpoint max_marginal matches",
              abs(float(margQ) - crow["max_marginal"]) < 1e-15)
        # their exact enclosures, my value inside, sign certifies
        perm = [0] * n
        for slot, coord in enumerate(crow["roll_seq"]):
            perm[coord] = slot
        muP = {}
        for a, w in muQ.items():
            b = 0
            for i in range(n):
                if (a >> i) & 1:
                    b |= 1 << perm[i]
            muP[b] = muP.get(b, Fraction(0)) + w
        zs = hu_cells(n, muP)
        mine = Fraction(str(crD))
        for kitname, fn in (("A_digit_extraction", log2_A),
                            ("B_atanh_series", log2_B)):
            Hiv = (Fraction(0), Fraction(0))
            for w in muP.values():
                Hiv = ivl_add(Hiv, ivl_scale(Fraction(-1),
                                             xlog2x_ivl(w, fn)))
            S = (Fraction(0), Fraction(0))
            for w, z in zs:
                if 0 < z < 1:
                    S = ivl_add(S, ivl_scale(w, h_ivl(z, fn)))
            lo, hi = S[0] - Hiv[1], S[1] - Hiv[0]
            kc = crow["kits"][kitname]
            check(f"{tag} [{kitname}]: my 60-digit CR inside the exact "
                  "enclosure; CR_lo > 0; checkpoint endpoints match",
                  lo <= mine <= hi and lo > 0
                  and float(lo) == kc["CR_lo"] and float(hi) == kc["CR_hi"],
                  f"CR {float(mine):+.9e} in [{float(lo):+.9e}, "
                  f"{float(hi):+.9e}]")


def v4(roc):
    print("V4. Census: spot re-scores, summary fields, canon_rank <= "
          "roll_rank on all 1200 rows:")
    total_rows = 0
    dom_bad = 0
    for key, blk in roc.items():
        if not key.startswith("P1_"):
            continue
        rows = blk["rows"]
        total_rows += len(rows)
        dom_bad += sum(1 for r in rows if r["canon_rank"] > r["roll_rank"])
        ok = (blk["instances"] == len(rows)
              and blk["roll_best"] == sum(1 for r in rows
                                          if r["roll_is_best"])
              and blk["worst_roll_rank"] == min(r["roll_rank"]
                                                for r in rows)
              and abs(blk["worst_gap_over_H"]
                      - max(r["gap_over_H"] for r in rows)) < 1e-15
              and blk["roll_negative"] == sum(1 for r in rows
                                              if r["roll_neg"])
              and blk["any_order_negative"] == sum(1 for r in rows
                                                   if r["any_neg"])
              and blk["canon_below_roll"] == sum(
                  1 for r in rows if r["canon_rank"] < r["roll_rank"]))
        check(f"{key}: all summary fields recompute from rows", ok)
        n = int(key.split("_")[1][1:])
        cap = float(key.split("cap")[1])
        rng = random.Random(938000 + n * 1000 + int(cap * 1000))
        regen, tried = [], 0
        while len(regen) < blk["instances"] and tried < 20 * blk["instances"]:
            tried += 1
            mu = gen_instance(rng, n, cap)
            H = -sum(w * math.log2(w) for w in mu.values() if w > 0)
            if H < 0.2:
                continue
            regen.append(mu)
        bad = []
        dom_val_bad = []
        for idx in range(0, blk["instances"], 25):
            mu, row = regen[idx], rows[idx]
            atoms = norm_atoms(mu)
            allcr = sorted(hu_cr(n, atoms, s)[0]
                           for s in itertools.permutations(range(n)))
            rollcr, H = hu_cr(n, atoms, rule_roll(n, atoms))
            cancr, _ = hu_cr(n, atoms, rule_can(n, atoms))
            rank = 1 + sum(1 for c in allcr if c < rollcr - 1e-12)
            crank = 1 + sum(1 for c in allcr if c < cancr - 1e-12)
            if not (abs(rollcr - row["roll_CR"]) < 1e-9
                    and abs(allcr[-1] - row["best_CR"]) < 1e-9
                    and abs(allcr[0] - row["worst_CR"]) < 1e-9
                    and abs(H - row["H"]) < 1e-9
                    and rank == row["roll_rank"]
                    and crank == row["canon_rank"]
                    and row["roll_is_best"] == (allcr[-1] - rollcr <= 1e-12)
                    and abs(row["gap_over_H"]
                            - (allcr[-1] - rollcr) / H) < 1e-9
                    and row["any_neg"] == (allcr[0] < 0)
                    and row["roll_neg"] == (rollcr < 0)):
                bad.append(idx)
            if cancr > rollcr + 1e-9:
                dom_val_bad.append(idx)
        check(f"{key}: 12 spot rows fully re-scored (roll/canon/best/"
              "worst/H/ranks/flags)", not bad, f"mismatched idx {bad}")
        check(f"{key}: CR_canon <= CR_roll on all spot rows "
              "((ROLL-DOM) by value)", not dom_val_bad,
              f"violations {dom_val_bad}")
    check("all 1200 rows: canon_rank <= roll_rank ((ROLL-DOM) by rank)",
          total_rows == 1200 and dom_bad == 0,
          f"rows {total_rows}, rank violations {dom_bad}")


def v5(can, att):
    print("V5. (ROLL-DOM) hinges on the committed code:")
    sys.path.insert(0, str(HERE))
    from uc_hu_canon import cond_entropy
    from uc_hu_order2 import canon_completion, seq_can, seq_roll, hu_cr_seq
    rng = random.Random(939000)
    insts = []
    for v in (att["cap_0.45"]["violations"][0],
              can["B_cap_0.49"]["violations"][0],
              can["B_cap_0.49"]["violations"][1]):
        mu = load_mu(v["mu"])
        tot = sum(mu.values())
        insts.append((v["n"], {a: w / tot for a, w in mu.items()}))
    for _ in range(200):
        n = rng.choice([3, 4, 5])
        mu = gen_instance(rng, n, 0.49)
        insts.append((n, mu))
    bad_set, bad_v0, bad_rec = [], [], []
    for n, m in insts:
        S = list(rng.sample(range(n), rng.randint(1, n - 1)))
        i = rng.choice([j for j in range(n) if j not in S])
        S2 = list(S)
        rng.shuffle(S2)
        if cond_entropy(n, m, i, S) != cond_entropy(n, m, i, S2):
            bad_set.append((n, S, S2))
        if seq_can(n, m) != canon_completion(n, m, []):
            bad_v0.append(n)
        pref = list(seq_can(n, m)[:rng.randint(0, n - 1)])
        comp = canon_completion(n, m, pref)
        f0 = comp[len(pref)]
        if canon_completion(n, m, pref + [f0]) != comp:
            bad_rec.append((n, pref))
    check("cond_entropy depends on the revealed SET only "
          "(order-shuffle invariant, exact, 203 instances)", not bad_set,
          f"{bad_set[:2]}")
    check("seq_can == canon_completion(empty prefix) (the V0 >= "
          "CR(canon) hinge)", not bad_v0, f"{bad_v0[:3]}")
    check("completion recursion: completion(P) == completion(P+[f(P)])",
          not bad_rec, f"{bad_rec[:2]}")
    # tie-tolerance stress: perturbed exchangeable measures
    worst_margin = 1e9
    tested = 0
    for trial in range(150):
        r2 = random.Random(940000 + trial)
        n = r2.choice([3, 4])
        p = r2.uniform(0.30, 0.485)
        mu = {a: p ** bin(a).count("1") * (1 - p) ** (n - bin(a).count("1"))
              for a in range(1 << n)}
        eps = 10.0 ** r2.uniform(-14, -9)
        for a in list(mu):
            mu[a] *= 1.0 + eps * r2.uniform(-1, 1)
        tot = sum(mu.values())
        m = {a: w / tot for a, w in mu.items()}
        cr_r = hu_cr_seq(n, m, seq_roll(n, m))[0]
        cr_c = hu_cr_seq(n, m, seq_can(n, m))[0]
        worst_margin = min(worst_margin, cr_r - cr_c)
        tested += 1
    check(f"tie-tolerance stress ({tested} perturbed-exchangeable "
          "instances): CR_roll >= CR_canon - n*TIE - float slack",
          worst_margin > -5 * 1e-12,
          f"worst CR_roll - CR_canon = {worst_margin:+.2e}")


def v6(ord2):
    print("V6. The 0.499 fixed-point story:")
    row = min(ord2["C2_roll_cap_0.499"]["rows"], key=lambda r: r["floor"])
    n = row["n"]
    muF = load_mu(row["mu"])
    muQ = rationalize(muF)
    # CORRECTION basis (037 lead 4 says "rationalizes to an
    # exactly-exchangeable measure"): in fact the 8.1e-11 atom
    # rationalizes to weight exactly 0, leaving TWO atoms {0000, 1101};
    # coordinate 1 has marginal exactly 0 and the stabilizer is only S3
    # on {0,2,3} (order 6 of 24) -- NOT exchangeable, and perfectly
    # correlated rather than (near-)product.  The 4-way step-0 tie is
    # real, but it comes from this symmetry + the inert coordinate.
    stab = sum(1 for sigma in itertools.permutations(range(n))
               if relabel_frac(n, muQ, sigma) == muQ)
    support = sorted(a for a, w in muQ.items() if w > 0)
    marg1 = sum(w for a, w in muQ.items() if (a >> 1) & 1)
    check("CORRECTION basis: mu_Q is NOT exchangeable -- 2-atom support "
          "{0000, 1101} (the 8.1e-11 atom rationalizes to 0), coord-1 "
          "marginal exactly 0, stabilizer order 6 of 24",
          stab == 6 and support == [0, 0b1101] and marg1 == 0,
          f"stabilizer {stab}/24, support {support}, marg1 {marg1}")
    scores = []
    for i in range(n):
        full = rule_can_dec(n, muQ, (i,))
        scores.append(hu_cr_dec(n, muQ, full)[0])
    check("four step-0 rollout scores of mu_Q EXACTLY equal at 60 digits",
          max(scores) == min(scores),
          f"spread {float(max(scores) - min(scores)):.1e}")
    # unrationalized float measure: exact-Fraction lift of the floats
    muE = {a: Fraction(w) for a, w in muF.items()}
    tot = sum(muE.values())
    muE = {a: w / tot for a, w in muE.items()}
    sc = []
    for i in range(n):
        full = rule_can_dec(n, muE, (i,))
        sc.append(hu_cr_dec(n, muE, full)[0])
    gap = float(max(sc) - min(sc))
    check("unrationalized step-0 score gap ~1.6e-10 (record's number)",
          1e-11 < gap < 1e-8, f"gap {gap:.3e}")
    atomsF = norm_atoms(muF)
    seqF = rule_roll(n, atomsF)
    check("roll(mu_float) != roll(mu_Q): rationalization changes the "
          "derived order (the instance the fixed-point check caught)",
          list(seqF) != [0, 1, 2, 3] or gap > 1e-12,
          f"roll(mu_float) = {seqF}")


def v7(ord2, roc, cert):
    print("V7. Reporting audit (record prose vs checkpoints):")
    be = ord2["B_endpoints"]
    check("037 B: roll worst +0.041983 / +0.007920 / +0.000679, 0 neg",
          abs(be["B_cap_0.38271"]["roll"]["worst"] - 0.041983) < 5e-7
          and abs(be["B_cap_0.45"]["roll"]["worst"] - 0.007920) < 5e-7
          and abs(be["B_cap_0.49"]["roll"]["worst"] - 0.000679) < 5e-7
          and all(be[k]["roll"]["negative_rows"] == 0 for k in be))
    check("037 C: surp floors -0.012831 / -0.020503 (2 violations each); "
          "roll floors +0.000679 / +0.000037 (0 violations)",
          abs(ord2["C_surp_cap_0.49"]["global_floor"] + 0.012831) < 5e-7
          and abs(ord2["C_surp_cap_0.497"]["global_floor"] + 0.020503) < 5e-7
          and len(ord2["C_surp_cap_0.49"]["violations"]) == 2
          and len(ord2["C_surp_cap_0.497"]["violations"]) == 2
          and abs(ord2["C_roll_cap_0.49"]["global_floor"] - 0.000679) < 5e-7
          and abs(ord2["C_roll_cap_0.497"]["global_floor"] - 0.000037) < 5e-7
          and not ord2["C_roll_cap_0.49"]["violations"]
          and not ord2["C_roll_cap_0.497"]["violations"])
    check("037 C2: floors +0.000679 / +0.000037 / +0.0000030, 0 violations",
          abs(ord2["C2_roll_cap_0.49"]["global_floor"] - 0.000679) < 5e-7
          and abs(ord2["C2_roll_cap_0.497"]["global_floor"] - 0.000037) < 5e-7
          and abs(ord2["C2_roll_cap_0.499"]["global_floor"] - 0.0000030) < 5e-8
          and all(not ord2[k]["violations"]
                  for k in ("C2_roll_cap_0.49", "C2_roll_cap_0.497",
                            "C2_roll_cap_0.499")))
    check("037 D: best-order floors +0.000694 / +0.000219, 0 violations",
          abs(ord2["D_bestorder_cap_0.49"]["global_floor"] - 0.000694) < 5e-7
          and abs(ord2["D_bestorder_cap_0.497"]["global_floor"]
                  - 0.000219) < 5e-7)
    quoted = {"rollrescue-floor:random0-n4": 1.507474125e-1,
              "rollrescue-seed4-n5": 3.495449146e-2,
              "rollfloor-0.49-product4": 2.714371934e-3,
              "rollfloor-0.497-seed2": 3.693561440e-5,
              "rollfloor-0.499-sharp:C_roll_cap_0.497:seed0": 3.003306123e-6}
    ok = True
    for w in cert["witnesses"]:
        lo = w["kits"]["A_digit_extraction"]["CR_lo"]
        if abs(lo - quoted[w["witness"]]) > 5e-10 * max(1, abs(lo) * 1e3):
            ok = False
    check("037 E: all five quoted certificate values match kit-A CR_lo",
          ok and all(w["fixed_point_60dig"] for w in cert["witnesses"]))
    p1 = [("P1_n4_cap0.45", 181, 7, 0.101), ("P1_n4_cap0.49", 171, 9, 0.207),
          ("P1_n5_cap0.45", 117, 61, 0.153), ("P1_n5_cap0.49", 158, 41, 0.091)]
    ok = all(roc[k]["instances"] == 300 and roc[k]["roll_best"] == rb
             and roc[k]["worst_roll_rank"] == wr
             and abs(roc[k]["worst_gap_over_H"] - wg) < 5e-4
             and roc[k]["roll_negative"] == 0 for k, rb, wr, wg in p1)
    anyneg = sum(roc[k]["any_order_negative"] for k, _, _, _ in p1)
    eqs = [sum(1 for r in roc[k]["rows"]
               if r["canon_rank"] == r["roll_rank"]) for k, _, _, _ in p1]
    check("038 P1: table (60.3/57.0/39.0/52.7%, ranks 7/9/61/41, gaps "
          ".101/.207/.153/.091), 0 roll-neg, 5 any-neg, equal on "
          "52/52/24/28", ok and anyneg == 5 and eqs == [52, 52, 24, 28],
          f"anyneg {anyneg}, equal {eqs}")
    check("038 P2/P3: floors +0.030430 / +0.053740 / +0.099692, 0 viol",
          abs(roc["P2_cap_0.49"]["global_floor"] - 0.030430) < 5e-7
          and abs(roc["P2_cap_0.497"]["global_floor"] - 0.053740) < 5e-7
          and abs(roc["P3_cap_0.49"]["global_floor"] - 0.099692) < 5e-7
          and all(not roc[k]["violations"]
                  for k in ("P2_cap_0.49", "P2_cap_0.497", "P3_cap_0.49")))
    # the corrected skipped-starts story (037 "three n=6 starts" / 038
    # "the three starts 037 skipped"): in fact two n=6 and ONE n=8, and
    # P2 covers only the two n=6 ones -- crash8 has no roll descent.
    from uc_hu_attack import starts
    sk = [(name, n) for name, n, _ in starts(0.49) if n > 5]
    p2names = {r["start"] for r in roc["P2_cap_0.49"]["rows"]}
    check("CORRECTION basis: part-C skip set is {mmabskill_n6 (n=6), "
          "crash8 (n=8), seed7 (n=6)}; P2 contains only the two n=6 "
          "starts (crash8 never roll-descended)",
          sorted(sk) == [("crash8", 8), ("floor:mmabskill_n6", 6),
                         ("seed7", 6)]
          and {"floor:mmabskill_n6", "seed7"} <= p2names
          and "crash8" not in p2names,
          f"skipped {sk}, P2 starts {sorted(p2names)}")


def main():
    ord2 = json.loads((DATA / "hu_order2.json").read_text())
    can = json.loads((DATA / "hu_canon.json").read_text())
    att = json.loads((DATA / "hu_attack.json").read_text())
    cert = json.loads((DATA / "hu_order2_certify.json").read_text())
    roc = json.loads((DATA / "hu_rollcensus.json").read_text())
    v1(ord2, can, att)
    v2(ord2, can, roc)
    v3(ord2, can, cert)
    v4(roc)
    v5(can, att)
    v6(ord2)
    v7(ord2, roc, cert)
    print()
    if FAILS:
        print(f"REFUTATIONS/FAILURES: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("All checks pass; no refutation found.")


if __name__ == "__main__":
    main()
