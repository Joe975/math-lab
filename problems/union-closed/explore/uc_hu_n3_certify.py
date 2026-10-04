#!/usr/bin/env python3
"""Attempt 050: exact certification of the n = 3 identity-order failure.

050 found that the order quantifier, vacuous at n = 2 (046), is already
essential at n = 3: the identity order fails the own-constant bound on
census instances at cap 0.49.  The sharpest of them is stronger than
"below c*" -- its CR is NEGATIVE:

  mu (4 atoms, n = 3, max marginal 0.485):
      000: 0.515, 011: 0.046779089, 101: 0.317092622, 111: 0.121128289

Statements, per kit alone (029 standard):
  N1  CR_HU < 0 under the IDENTITY order          (the failure)
  N2  CR_HU > 0 under the best enumerated order   (so it is the order,
                                                   not the measure)
  N3  max marginal < 1/2 and H(mu) > 1/2, exact rational checks

Usage: python uc_hu_n3_certify.py
Checkpoint: ../data/hu_n3_certify.json
"""
from __future__ import annotations

import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_certify import (hu_cells, ivl_add, ivl_scale, xlog2x_ivl, h_ivl,
                           log2_A, log2_B)
from uc_hu_order2 import hu_cr_seq

DATA = HERE.parent / "data"
HALF = Fraction(1, 2)

# masks as integers (the record displays them as format(a, "03b"),
# most-significant-bit first): 000 -> 0, 011 -> 3, 101 -> 5, 111 -> 7
WITNESS = {0: 0.515, 3: 0.046779089, 5: 0.317092622, 7: 0.121128289}


def permute(n, mu, perm):
    out = {}
    for a, w in mu.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << perm[i]
        out[b] = out.get(b, Fraction(0)) + w
    return out


def enclose(n, muP, fn):
    Hmu = (Fraction(0), Fraction(0))
    for w in muP.values():
        Hmu = ivl_add(Hmu, ivl_scale(Fraction(-1), xlog2x_ivl(w, fn)))
    S = (Fraction(0), Fraction(0))
    for w, z in hu_cells(n, muP):
        if 0 < z < 1:
            S = ivl_add(S, ivl_scale(w, h_ivl(z, fn)))
    return (S[0] - Hmu[1], S[1] - Hmu[0]), Hmu


def main():
    n = 3
    muQ = {a: Fraction(w).limit_denominator(10 ** 7)
           for a, w in WITNESS.items()}
    tot = sum(muQ.values())
    muQ = {a: w / tot for a, w in muQ.items()}
    marg = max(sum(w for a, w in muQ.items() if (a >> i) & 1)
               for i in range(n))
    muF = {a: float(w) for a, w in muQ.items()}
    best_seq = max(itertools.permutations(range(n)),
                   key=lambda s: hu_cr_seq(n, muF, s)[0])
    print("Attempt 050: the n = 3 identity-order failure, certified")
    print(f"  exact max marginal {float(marg):.6f} < 1/2; best order "
          f"{best_seq}")
    out = {"n": n, "max_marginal": float(marg),
           "best_seq": list(best_seq), "kits": {}}
    for tag, seq, want in (("N1_identity", tuple(range(n)), "neg"),
                           ("N2_best", best_seq, "pos")):
        perm = [0] * n
        for slot, coord in enumerate(seq):
            perm[coord] = slot
        muP = permute(n, muQ, perm)
        for kitname, fn in (("A_digit_extraction", log2_A),
                            ("B_atanh_series", log2_B)):
            CR, Hmu = enclose(n, muP, fn)
            ok = ((CR[1] < 0 if want == "neg" else CR[0] > 0)
                  and Hmu[0] > HALF and marg < HALF)
            out["kits"][f"{tag}:{kitname}"] = {
                "CR_lo": float(CR[0]), "CR_hi": float(CR[1]),
                "H_lo": float(Hmu[0]), "certifies": ok}
            print(f"  {tag} [{kitname}]: CR_HU in "
                  f"[{float(CR[0]):+.9e}, {float(CR[1]):+.9e}], "
                  f"H > {float(Hmu[0]):.4f} "
                  f"{'CERTIFIED ' + ('< 0' if want == 'neg' else '> 0') if ok else 'NOT CERTIFIED'}")
            if not ok:
                raise SystemExit(f"certification failed: {tag} {kitname}")
    (DATA / "hu_n3_certify.json").write_text(
        json.dumps(out, indent=1, default=float) + "\n")
    print("checkpoint: data/hu_n3_certify.json")


if __name__ == "__main__":
    main()
