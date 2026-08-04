#!/usr/bin/env python3
"""Independent re-computation of hierarchical-laminate records (2D conductivity).

Same objects as laminate.py, deliberately different algebra: instead of the
frame-free projection formula, each lamination is computed by

  1. rotating both tensors into the layer frame (normal along e1).  For an
     integer direction n = (u, v) the rotation conjugation is exactly rational:
     with S = [[u, -v], [v, u]],  R = S / sqrt(u^2 + v^2)  and
     R^T A R = S^T A S / (u^2 + v^2);
  2. averaging in the layer frame with the classical interface conditions
     (tangential E and normal J continuous).  For tensors [[a, b], [b, c]]
     with weights w summing to 1:
         s11 = <1/a>^{-1}
         s12 = s11 * <b/a>
         s22 = <c - b^2/a> + s11 * <b/a>^2
  3. rotating back:  sigma* = S sigma*_frame S^T / (u^2 + v^2).

Given a JSON record produced by laminate.py (tree, phases, tensor, fractions,
bounds), this tool recomputes every field from scratch -- effective tensor by
the route above, volume fractions, Wiener and Hashin-Shtrikman bounds from
independently coded formulas, and the exact Keller-Dykhne duality identity --
and reports any disagreement.  Agreement of the two algebraic routes on the
same tree is the two-implementation bar computational claims must clear.

Usage:
    python verify_laminate.py RECORD.json
    python verify_laminate.py --selftest

Standard library only.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction

Tensor = tuple[Fraction, Fraction, Fraction]  # [[a, b], [b, c]]


def iso(s: Fraction) -> Tensor:
    return (s, Fraction(0), s)


def conj_S(t: Tensor, u: int, v: int, inverse: bool) -> Tensor:
    """S^T t S / (u^2+v^2) with S = [[u, -v], [v, u]]; inverse swaps S and S^T."""
    a, b, c = t
    if inverse:
        u, v = u, -v  # (S^T)^T t S^T pattern via sign flip
    # M = S^T A S, entries expanded
    m00 = u * (a * u + b * v) + v * (b * u + c * v)
    m01 = u * (-a * v + b * u) + v * (-b * v + c * u)
    m11 = -v * (-a * v + b * u) + u * (-b * v + c * u)
    d = Fraction(u * u + v * v)
    return (m00 / d, m01 / d, m11 / d)


def backus_pair(wA: Fraction, A: Tensor, B: Tensor) -> Tensor:
    """Layer-frame average of A (weight wA) and B (1 - wA), normal along e1."""
    wB = 1 - wA
    inv_a = wA / A[0] + wB / B[0]
    s11 = 1 / inv_a
    b_over_a = wA * A[1] / A[0] + wB * B[1] / B[0]
    schur = wA * (A[2] - A[1] * A[1] / A[0]) + wB * (B[2] - B[1] * B[1] / B[0])
    s12 = s11 * b_over_a
    s22 = schur + s11 * b_over_a * b_over_a
    return (s11, s12, s22)


def laminate_frame(A: Tensor, B: Tensor, m: Fraction, n: tuple[int, int]) -> Tensor:
    u, v = int(n[0]), int(n[1])
    assert (u, v) != (0, 0) and 0 < m < 1
    Af = conj_S(A, u, v, inverse=False)
    Bf = conj_S(B, u, v, inverse=False)
    Ef = backus_pair(m, Af, Bf)
    return conj_S(Ef, u, v, inverse=True)


def effective(tree: dict, phases: dict[str, Fraction]) -> Tensor:
    if "phase" in tree:
        return iso(phases[tree["phase"]])
    return laminate_frame(
        effective(tree["layers"][0], phases),
        effective(tree["layers"][1], phases),
        Fraction(tree["fraction"]),
        tuple(tree["normal"]),
    )


def fractions_of(tree: dict) -> dict[str, Fraction]:
    if "phase" in tree:
        return {tree["phase"]: Fraction(1)}
    m = Fraction(tree["fraction"])
    out: dict[str, Fraction] = {}
    for sub, w in ((tree["layers"][0], m), (tree["layers"][1], 1 - m)):
        for k, v in fractions_of(sub).items():
            out[k] = out.get(k, Fraction(0)) + w * v
    return out


def t_inv(t: Tensor) -> Tensor:
    a, b, c = t
    det = a * c - b * b
    return (c / det, -b / det, a / det)


def keller_identity(tree: dict, phases: dict[str, Fraction]) -> bool:
    sigma = effective(tree, phases)
    dual = effective(tree, {k: 1 / v for k, v in phases.items()})
    a, b, c = t_inv(sigma)
    return dual == (c, -b, a)  # R (.) R^T for R = 90-degree rotation


def wiener(fracs: list[Fraction], sigs: list[Fraction]) -> tuple[Fraction, Fraction]:
    return 1 / sum(f / s for f, s in zip(fracs, sigs)), sum(f * s for f, s in zip(fracs, sigs))


def hashin_shtrikman(fracs: list[Fraction], sigs: list[Fraction]) -> tuple[Fraction, Fraction]:
    lo_c, hi_c = min(sigs), max(sigs)
    return (
        1 / sum(f / (s + lo_c) for f, s in zip(fracs, sigs)) - lo_c,
        1 / sum(f / (s + hi_c) for f, s in zip(fracs, sigs)) - hi_c,
    )


# ---------------------------------------------------------------------------


def verify_record(path: str) -> list[str]:
    with open(path, encoding="utf-8") as fh:
        rec = json.load(fh)
    phases = {k: Fraction(v) for k, v in rec["phases"].items()}
    tree = rec["tree"]
    problems: list[str] = []

    t = effective(tree, phases)
    claimed = tuple(Fraction(x) for x in rec["tensor"])
    if t != claimed:
        problems.append(f"tensor mismatch: recomputed {t}, claimed {claimed}")

    fr = fractions_of(tree)
    claimed_fr = {k: Fraction(v) for k, v in rec["fractions"].items()}
    if fr != claimed_fr:
        problems.append(f"fraction mismatch: recomputed {fr}, claimed {claimed_fr}")
    if sum(fr.values()) != 1:
        problems.append("fractions do not sum to 1")

    names = sorted(fr)
    fracs = [fr[n] for n in names]
    sigs = [phases[n] for n in names]
    if [str(x) for x in wiener(fracs, sigs)] != rec["wiener"]:
        problems.append("Wiener bounds mismatch")
    if [str(x) for x in hashin_shtrikman(fracs, sigs)] != rec["hs"]:
        problems.append("Hashin-Shtrikman bounds mismatch")

    if not keller_identity(tree, phases):
        problems.append("Keller-Dykhne identity fails (implementation or record bug)")
    return problems


def selftest() -> None:
    F = Fraction
    s1, s2 = F(1), F(3)
    # rank-1, normal e1: diag(harmonic, arithmetic)
    m = F(2, 5)
    t = laminate_frame(iso(s1), iso(s2), m, (1, 0))
    harm = 1 / (m / s1 + (1 - m) / s2)
    arith = m * s1 + (1 - m) * s2
    assert t == (harm, F(0), arith), t
    # rank-1, normal e2: the transpose arrangement
    t = laminate_frame(iso(s1), iso(s2), m, (0, 1))
    assert t == (arith, F(0), harm), t
    # direction scale-invariance and duality on a nested tree
    tree = {
        "fraction": "1/3",
        "normal": [2, -1],
        "layers": [
            {"phase": "p1"},
            {"fraction": "4/9", "normal": [1, 5], "layers": [{"phase": "p2"}, {"phase": "p3"}]},
        ],
    }
    phases = {"p1": F(1), "p2": F(3), "p3": F(7)}
    assert laminate_frame(iso(s1), iso(s2), m, (3, 4)) == laminate_frame(iso(s1), iso(s2), m, (-6, -8))
    assert keller_identity(tree, phases)
    # laminating a material with itself changes nothing
    assert laminate_frame(iso(s2), iso(s2), F(1, 7), (5, 2)) == iso(s2)
    print("selftest: all checks passed")


def main() -> None:
    ap = argparse.ArgumentParser(description="Independent laminate-record verification (2D)")
    ap.add_argument("record", nargs="?", help="JSON record produced by the reference tool")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest()
        return
    if not args.record:
        ap.error("need RECORD (or --selftest)")
    problems = verify_record(args.record)
    if problems:
        print("VERIFICATION FAILED:")
        for p in problems:
            print(f"  - {p}")
        sys.exit(1)
    print("record verified: tensor, fractions, bounds and duality all reproduce")


if __name__ == "__main__":
    main()
