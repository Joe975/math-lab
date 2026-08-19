#!/usr/bin/env python3
"""Attempt 035: exact certification of the CANONICAL-order HU failures at
marginal cap 49/100 (and of the canonical rescue at cap 9/20).

Single-kit per 029: each certified-log2 kit must certify ALONE; agreement
is reported as a consistency check only.  Interval helpers and the exact
HU cell tree are reused from uc_hu_certify (033, reviewer-checked at
float level by uc_hu_canon_skeptic before certification here).

Statements, per witness and per kit:
  D1  canonical-order CR_HU < 0 at cap 49/100 (the two 035 witnesses)
  D2  max marginal < 49/100 and H(mu) > 1/2, exact rational checks
  D3  the 033 witness under the CANONICAL order: CR_HU > 0 at cap 9/20
      (the rescue, certified)

Usage: python uc_hu_canon_certify.py
Checkpoint: ../data/hu_canon_certify.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_certify import (hu_cells, ivl_add, ivl_scale, xlog2x_ivl, h_ivl,
                           log2_A, log2_B)
from uc_hu_canon import canonical_order

DATA = HERE.parent / "data"
HALF = Fraction(1, 2)
OUT = {"witnesses": []}


def rationalize(mu_float):
    mu = {a: Fraction(w).limit_denominator(10 ** 7)
          for a, w in mu_float.items()}
    tot = sum(mu.values())
    return {a: w / tot for a, w in mu.items()}


def permute(n, mu, perm):
    out = {}
    for a, w in mu.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << perm[i]
        out[b] = out.get(b, Fraction(0)) + w
    return out


def certify(tag, n, mu_float, cap_frac, want):
    """want = 'neg' or 'pos'; the canonical order is derived from the
    FLOAT measure (the rule's own definition) and then applied exactly."""
    perm = canonical_order(n, {a: float(w) for a, w in mu_float.items()})
    mu = permute(n, rationalize(mu_float), perm)
    marg = max(sum(w for a, w in mu.items() if (a >> i) & 1)
               for i in range(n))
    assert marg < cap_frac, f"{tag}: marginal {float(marg)} >= cap"
    zs = hu_cells(n, mu)
    row = {"witness": tag, "n": n, "cap": str(cap_frac),
           "canon_order": list(perm), "max_marginal": float(marg),
           "want": want, "kits": {}}
    for kitname, log2fn in (("A_digit_extraction", log2_A),
                            ("B_atanh_series", log2_B)):
        Hmu = (Fraction(0), Fraction(0))
        for w in mu.values():
            Hmu = ivl_add(Hmu, ivl_scale(Fraction(-1),
                                         xlog2x_ivl(w, log2fn)))
        S = (Fraction(0), Fraction(0))
        for w, z in zs:
            if 0 < z < 1:
                S = ivl_add(S, ivl_scale(w, h_ivl(z, log2fn)))
        CR = (S[0] - Hmu[1], S[1] - Hmu[0])
        ok = (CR[1] < 0 if want == "neg" else CR[0] > 0) and Hmu[0] > HALF
        row["kits"][kitname] = {"CR_lo": float(CR[0]),
                                "CR_hi": float(CR[1]),
                                "H_lo": float(Hmu[0]), "certifies": ok}
        print(f"  {tag} [{kitname}]: CR_HU in "
              f"[{float(CR[0]):+.9e}, {float(CR[1]):+.9e}], "
              f"H > {float(Hmu[0]):.4f}, marg {float(marg):.5f} "
              f"{'CERTIFIED ' + ('< 0' if want == 'neg' else '> 0') if ok else 'NOT CERTIFIED'}",
              flush=True)
        if not ok:
            raise SystemExit(f"certification failed: {tag} {kitname}")
    OUT["witnesses"].append(row)


def main():
    print("Attempt 035: certification of canonical-order HU at the "
          "0.49 failures and the 0.45 rescue")
    can = json.loads((DATA / "hu_canon.json").read_text())
    for r in can["B_cap_0.49"]["violations"]:
        mu = {int(s, 2): w for s, w in r["mu"].items()}
        certify(f"canon-fail-{r['start']}-n{r['n']}", r["n"], mu,
                Fraction(49, 100), "neg")
    att = json.loads((DATA / "hu_attack.json").read_text())
    w = att["cap_0.45"]["violations"][0]
    mu = {int(s, 2): w2 for s, w2 in w["mu"].items()}
    certify("canon-rescue-033witness", w["n"], mu, Fraction(9, 20), "pos")
    (DATA / "hu_canon_certify.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    print("checkpoint: data/hu_canon_certify.json")


if __name__ == "__main__":
    main()
