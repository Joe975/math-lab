#!/usr/bin/env python3
"""Attempt 037: exact certification of the ROLLOUT-order HU positives.

Certified statements, per instance and per kit (single-kit per 029):
  E1  CR_HU > 0 under the rollout order at both 035 cap-0.49 kill
      measures (the rescues)
  E2  CR_HU > 0 at the sharpest rollout descent floors, caps 49/100,
      497/1000 and 499/1000
  E3  side conditions: max marginal < cap and H(mu) > 1/2, exact

Order provenance: the rollout order is derived from the RATIONALIZED
measure (the certified object) via its float projection, then applied
exactly.  Going beyond the 035 precedent, per 036's review, a 60-digit
Decimal FIXED-POINT check runs on the rationalized measure itself: at
every step of the certified sequence, every alternative coordinate's
rollout score (full CR with canonical completion) is recomputed at 60
digits and the certified choice is confirmed to be the argmax (ties ->
lowest index).  This caught a real instance: at the 0.499 floor the
limit_denominator(1e7) rationalization makes the four step-0 scores
EXACTLY equal (gap ~1.6e-10 on the original float measure), so
roll(mu_Q) != roll(mu_float) there -- the certified order must come
from mu_Q.

Usage: python uc_hu_order2_certify.py
Checkpoint: ../data/hu_order2_certify.json
"""
from __future__ import annotations

import json
import sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_certify import (hu_cells, ivl_add, ivl_scale, xlog2x_ivl, h_ivl,
                           log2_A, log2_B)
from uc_hu_order2 import seq_roll, canon_completion
from uc_reviewer036_reimpl import hu_cr as hu_cr_dec

DATA = HERE.parent / "data"
HALF = Fraction(1, 2)
OUT = {"witnesses": []}


def rationalize(mu_float):
    mu = {a: Fraction(w).limit_denominator(10 ** 7)
          for a, w in mu_float.items()}
    tot = sum(mu.values())
    return {a: w / tot for a, w in mu.items()}


def fixed_point_check(n, muQ, seq):
    """60-digit rollout fixed-point on the rationalized measure: at each
    step the certified coordinate maximizes full CR with canonical
    completion (candidate completions derived by the float rule on the
    float projection of muQ; scores at 60 digits)."""
    muF = {a: float(w) for a, w in muQ.items()}
    chosen = []
    for k, pick in enumerate(seq):
        rem = sorted(set(range(n)) - set(chosen))
        scores = {}
        for i in rem:
            full = canon_completion(n, muF, chosen + [i])
            cr, _, _ = hu_cr_dec(n, muQ, full)
            scores[i] = cr
        hi = max(scores.values())
        tied = [i for i in rem if hi - scores[i] <= Decimal("1e-12")]
        if min(tied) != pick:
            return False, f"step {k}: argmax {min(tied)} != certified {pick}"
        chosen.append(pick)
    return True, "fixed point at 60 digits"


def certify(tag, n, mu_float, cap_frac):
    muQ = rationalize(mu_float)
    # the certified object is muQ, so the rule is applied to muQ (via its
    # float projection; the 60-digit fixed-point check below then confirms
    # the choice on muQ itself)
    seq = seq_roll(n, {a: float(w) for a, w in muQ.items()})
    perm = [0] * n
    for slot, coord in enumerate(seq):
        perm[coord] = slot
    ok_fp, fp_msg = fixed_point_check(n, muQ, seq)
    muP = {}
    for a, w in muQ.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << perm[i]
        muP[b] = muP.get(b, Fraction(0)) + w
    marg = max(sum(w for a, w in muP.items() if (a >> i) & 1)
               for i in range(n))
    assert marg < cap_frac, f"{tag}: marginal {float(marg)} >= cap"
    zs = hu_cells(n, muP)
    row = {"witness": tag, "n": n, "cap": str(cap_frac),
           "roll_seq": list(seq), "max_marginal": float(marg),
           "fixed_point_60dig": ok_fp, "kits": {}}
    for kitname, log2fn in (("A_digit_extraction", log2_A),
                            ("B_atanh_series", log2_B)):
        Hmu = (Fraction(0), Fraction(0))
        for w in muP.values():
            Hmu = ivl_add(Hmu, ivl_scale(Fraction(-1),
                                         xlog2x_ivl(w, log2fn)))
        S = (Fraction(0), Fraction(0))
        for w, z in zs:
            if 0 < z < 1:
                S = ivl_add(S, ivl_scale(w, h_ivl(z, log2fn)))
        CR = (S[0] - Hmu[1], S[1] - Hmu[0])
        ok = CR[0] > 0 and Hmu[0] > HALF
        row["kits"][kitname] = {"CR_lo": float(CR[0]), "CR_hi": float(CR[1]),
                                "H_lo": float(Hmu[0]), "certifies": ok}
        print(f"  {tag} [{kitname}]: CR_HU in [{float(CR[0]):+.9e}, "
              f"{float(CR[1]):+.9e}], H > {float(Hmu[0]):.4f}, "
              f"marg {float(marg):.5f} "
              f"{'CERTIFIED > 0' if ok else 'NOT CERTIFIED'}", flush=True)
        if not ok:
            raise SystemExit(f"certification failed: {tag} {kitname}")
    print(f"  {tag} [order]: {fp_msg}", flush=True)
    if not ok_fp:
        raise SystemExit(f"fixed-point check failed: {tag}")
    OUT["witnesses"].append(row)


def main():
    print("Attempt 037: certification of rollout-order HU rescues and "
          "floors")
    can = json.loads((DATA / "hu_canon.json").read_text())
    ord2 = json.loads((DATA / "hu_order2.json").read_text())
    for r in can["B_cap_0.49"]["violations"]:
        mu = {int(s, 2): w for s, w in r["mu"].items()}
        certify(f"rollrescue-{r['start']}-n{r['n']}", r["n"], mu,
                Fraction(49, 100))
    for key, capf in (("C_roll_cap_0.49", Fraction(49, 100)),
                      ("C_roll_cap_0.497", Fraction(497, 1000)),
                      ("C2_roll_cap_0.499", Fraction(499, 1000))):
        row = min(ord2[key]["rows"], key=lambda r: r["floor"])
        mu = {int(s, 2): w for s, w in row["mu"].items()}
        certify(f"rollfloor-{key.split('_')[-1]}-{row['start']}",
                row["n"], mu, capf)
    (DATA / "hu_order2_certify.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    print("checkpoint: data/hu_order2_certify.json")


if __name__ == "__main__":
    main()
