#!/usr/bin/env python3
"""Build and save the below-m11 laminates that sit under Cherkaev's B2 bound.

These are the structures of the open escalation recorded in attempt 006.  They
come from the same family as `tp_attain.py` but with f1 pushed BELOW
m11 = 2*Theta*sqrt(m2)*(1-sqrt(m2)), where that closed form no longer applies:
a1 is re-solved from the volume constraint and a5 is re-solved for isotropy by
exact bisection, leaving a structure that is isotropic to ~1e-60 rather than
exactly.  That residual is many orders of magnitude below the ~1e-4 gap to B2,
so it cannot explain the discrepancy.

An adjudication pass flagged that only one of the three tabulated rows had a
saved, auditable record.  This module builds and writes ALL of them, so every
number in 006's table is reproducible from the repo.

Usage:
    python tp_below_m11.py --selftest
    python tp_below_m11.py --write        writes data/attained/below-m11-*.json
Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "three-phase-conductivity"))
import laminate as L  # noqa: E402
import verify_laminate as V  # noqa: E402

sys.path.insert(0, HERE)
import tp_cherkaev_bound as CB  # noqa: E402

SIGMA = (Fr(1), Fr(2), Fr(5))
PHASES = {"p1": SIGMA[0], "p2": SIGMA[1], "p3": SIGMA[2]}
M2 = Fr(1, 4)

# (m1, a0, a3) -- a0 and a3 are the two free parameters, taken from the coarse
# minimisation recorded in 006; they are not claimed optimal.
ROWS = [
    (Fr(31, 250), Fr(1, 2), Fr(11, 20)),
    (Fr(3, 25), Fr(101, 200), Fr(1, 2)),
    (Fr(11, 100), Fr(51, 100), Fr(1, 2)),
]


def build(m1: Fr, a0: Fr, a3: Fr, a5: Fr):
    Q = (1 - a0) - M2
    if Q <= 0:
        return None
    a4 = (M2 / (1 - a0) - (1 - a3)) / a3
    a1 = (m1 - Q * (1 - a5)) / a0
    if not all(0 < x < 1 for x in (a0, a1, a3, a4, a5)):
        return None
    A = {"fraction": str(a1), "normal": [0, 1],
         "layers": [{"phase": "p1"}, {"phase": "p3"}]}
    D = {"fraction": str(a5), "normal": [1, 0],
         "layers": [{"phase": "p3"}, {"phase": "p1"}]}
    C = {"fraction": str(a4), "normal": [0, 1], "layers": [{"phase": "p2"}, D]}
    B = {"fraction": str(a3), "normal": [0, 1], "layers": [C, {"phase": "p2"}]}
    return {"fraction": str(a0), "normal": [1, 0], "layers": [A, B]}


def isotropic_member(m1: Fr, a0: Fr, a3: Fr, iters: int = 220):
    """Bisect a5 for sigma11 == sigma22. Exact Fractions throughout."""
    Q = (1 - a0) - M2

    def aniso(x):
        t = build(m1, a0, a3, x)
        if t is None:
            return None
        e = L.effective(t, PHASES)
        return e[0] - e[2]

    lo = max(Fr(0), 1 - m1 / Q) + Fr(1, 10 ** 9)
    hi = 1 - Fr(1, 10 ** 9)
    glo, ghi = aniso(lo), aniso(hi)
    if glo is None or ghi is None or (glo > 0) == (ghi > 0):
        return None
    for _ in range(iters):
        mid = (lo + hi) / 2
        g = aniso(mid)
        if g is None:
            return None
        if (g > 0) == (glo > 0):
            lo = mid
        else:
            hi = mid
    return build(m1, a0, a3, (lo + hi) / 2)


def record_for(tree):
    e = L.effective(tree, PHASES)
    fr = L.fractions_of(tree)
    names = sorted(fr)
    fracs = [fr[n] for n in names]
    sigs = [PHASES[n] for n in names]
    return {
        "tree": tree,
        "phases": {k: str(v) for k, v in PHASES.items()},
        "tensor": [str(x) for x in e],
        "fractions": {k: str(v) for k, v in fr.items()},
        "wiener": [str(x) for x in L.wiener_bounds(fracs, sigs)],
        "hs": [str(x) for x in L.hs_bounds(fracs, sigs)],
    }


def summarise(m1: Fr, a0: Fr, a3: Fr):
    tree = isotropic_member(m1, a0, a3)
    assert tree is not None, m1
    e = L.effective(tree, PHASES)
    fr = L.fractions_of(tree)
    assert fr["p1"] == m1 and fr["p2"] == M2, (fr, m1)
    assert e == V.effective(tree, PHASES), "harness routes disagree"
    assert L.keller_check(tree, PHASES), "Keller-Dykhne fails"
    assert e[1] == 0, "off-diagonal must vanish"
    hs = CB.hs_lo(SIGMA, m1, M2)
    b2 = CB.B2(SIGMA, m1, Fr(1, 2))
    val = (e[0] + e[2]) / 2
    return tree, {"m1": str(m1), "value": val, "hs": hs, "b2": b2,
                  "iso_residual": e[0] - e[2]}


def selftest() -> None:
    for m1, a0, a3 in ROWS:
        _, s = summarise(m1, a0, a3)
        assert s["value"] > s["hs"], "must exceed HS"
        assert s["value"] < s["b2"], "the open discrepancy: below B2"
        assert abs(s["iso_residual"]) < Fr(1, 10 ** 40), s["iso_residual"]
    print("selftest: all three below-m11 rows rebuild, are isotropic to <1e-40, "
          "exceed HS_lo and sit below the transcribed B2")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--write", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    if a.write:
        out = os.path.join(HERE, "..", "data", "attained")
        os.makedirs(out, exist_ok=True)
        for m1, a0, a3 in ROWS:
            tree, s = summarise(m1, a0, a3)
            name = f"below-m11-m1-{m1.numerator}_{m1.denominator}.json"
            with open(os.path.join(out, name), "w", encoding="utf-8") as fh:
                json.dump(record_for(tree), fh, indent=1)
            print(f"{name}: value {float(s['value']):.10f}  HS {float(s['hs']):.10f}"
                  f"  B2 {float(s['b2']):.10f}  ours-B2 {float(s['value']-s['b2']):+.3e}")


if __name__ == "__main__":
    main()
