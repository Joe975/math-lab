#!/usr/bin/env python3
"""Attempt 044: exact certification of the rollout-order kill found by
the own-constant sweep.

The witness: the D_bestorder_cap_0.497 endpoint bred from
floor:windowkill (hu_order2.json), a 9-atom n=4 measure with max
marginal 0.49500 -- CR under the ROLLOUT order is negative (both float
stacks agree at -9.597e-4), while the best order stays positive
(+0.0161).  Statements, per kit (single-kit per 029), on the
limit_denominator(1e7) rationalization:

  K1  CR_HU under the rollout order < 0
  K2  max marginal < 497/1000 (and < 1/2), H > 1/2, exact
  K3  the certified order is the rollout order of the rationalized
      measure (60-digit fixed point, 037/043 pattern)
  K4  positive control: CR_HU under the best enumerated order > 0
      (so the kill is order-specific, not a dead measure)

Usage: python uc_hu_ownconst_certify.py
Checkpoint: ../data/hu_ownconst_certify.json
"""
from __future__ import annotations

import itertools
import json
import sys
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_certify import (hu_cells, ivl_add, ivl_scale, xlog2x_ivl, h_ivl,
                           log2_A, log2_B)
from uc_hu_order2 import seq_roll, canon_completion, hu_cr_seq
from uc_reviewer036_reimpl import hu_cr as hu_cr_dec

DATA = HERE.parent / "data"
HALF = Fraction(1, 2)
OUT = {}


def kit_enclosure(n, muP, log2fn):
    Hmu = (Fraction(0), Fraction(0))
    for w in muP.values():
        Hmu = ivl_add(Hmu, ivl_scale(Fraction(-1), xlog2x_ivl(w, log2fn)))
    S = (Fraction(0), Fraction(0))
    for w, z in hu_cells(n, muP):
        if 0 < z < 1:
            S = ivl_add(S, ivl_scale(w, h_ivl(z, log2fn)))
    return (S[0] - Hmu[1], S[1] - Hmu[0]), Hmu


def permute(n, mu, perm):
    out = {}
    for a, w in mu.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << perm[i]
        out[b] = out.get(b, Fraction(0)) + w
    return out


def main():
    d = json.loads((DATA / "hu_order2.json").read_text())
    row = [r for r in d["D_bestorder_cap_0.497"]["rows"]
           if r["start"] == "floor:windowkill"][0]
    n = row["n"]
    muF = {int(s, 2): w for s, w in row["mu"].items()}
    muQ = {a: Fraction(w).limit_denominator(10 ** 7)
           for a, w in muF.items()}
    tot = sum(muQ.values())
    muQ = {a: w / tot for a, w in muQ.items()}

    seq = seq_roll(n, {a: float(w) for a, w in muQ.items()})
    # K3: 60-digit rollout fixed point on muQ
    chosen = []
    fp_ok = True
    for k, pick in enumerate(seq):
        rem = sorted(set(range(n)) - set(chosen))
        scores = {}
        for i in rem:
            full = canon_completion(n, {a: float(w) for a, w in muQ.items()},
                                    chosen + [i])
            cr, _, _ = hu_cr_dec(n, muQ, full)
            scores[i] = cr
        hi = max(scores.values())
        tied = [i for i in rem if hi - scores[i] <= Decimal("1e-12")]
        if min(tied) != pick:
            fp_ok = False
            break
        chosen.append(pick)
    print(f"K3 fixed point at 60 digits: {'PASS' if fp_ok else 'FAIL'} "
          f"(roll seq {seq})")
    if not fp_ok:
        raise SystemExit("fixed-point check failed")

    # K4 positive control: best order by float enumeration, certified > 0
    best_seq = max(itertools.permutations(range(n)),
                   key=lambda s: hu_cr_seq(
                       n, {a: float(w) for a, w in muQ.items()}, s)[0])

    marg = max(sum(w for a, w in muQ.items() if (a >> i) & 1)
               for i in range(n))
    assert marg < Fraction(497, 1000) < HALF
    OUT.update({"witness": "ownconst-kill-D497-windowkill", "n": n,
                "roll_seq": list(seq), "best_seq": list(best_seq),
                "max_marginal": float(marg),
                "fixed_point_60dig": fp_ok, "kits": {}})
    for tag, s, want in (("K1_roll", seq, "neg"),
                         ("K4_best", best_seq, "pos")):
        perm = [0] * n
        for slot, coord in enumerate(s):
            perm[coord] = slot
        muP = permute(n, muQ, perm)
        for kitname, fn in (("A_digit_extraction", log2_A),
                            ("B_atanh_series", log2_B)):
            CR, Hmu = kit_enclosure(n, muP, fn)
            ok = ((CR[1] < 0 if want == "neg" else CR[0] > 0)
                  and Hmu[0] > HALF)
            OUT["kits"][f"{tag}:{kitname}"] = {
                "CR_lo": float(CR[0]), "CR_hi": float(CR[1]),
                "H_lo": float(Hmu[0]), "certifies": ok}
            print(f"  {tag} [{kitname}]: CR_HU in [{float(CR[0]):+.9e}, "
                  f"{float(CR[1]):+.9e}], marg {float(marg):.5f} "
                  f"{'CERTIFIED' if ok else 'NOT CERTIFIED'} "
                  f"{'< 0' if want == 'neg' else '> 0'}", flush=True)
            if not ok:
                raise SystemExit(f"certification failed: {tag} {kitname}")
    (DATA / "hu_ownconst_certify.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    print("checkpoint: data/hu_ownconst_certify.json")


if __name__ == "__main__":
    main()
