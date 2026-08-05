#!/usr/bin/env python3
"""Exact effective tensors of hierarchical laminates in 2D conductivity.

Objects: 2x2 symmetric positive-definite conductivity tensors over Q,
represented as (a, b, c) for [[a, b], [b, c]].  A *hierarchical laminate* is a
binary tree whose leaves are (isotropic, rational) phases and whose internal
nodes laminate their two children in a rational direction with a rational
volume fraction.  Every such structure has an exactly rational effective
tensor: the lamination formula

    sigma* = m A + (1-m) B - m(1-m) (dA n)(dA n)^T / (n . ((1-m)A + mB) n),

with dA = A - B and unnormalized integer direction n (the formula is
scale-invariant in n), is a rational function of its inputs.  So attainment
claims can be exact equalities in Q, which is the verification currency of
this problem.

Also computed, all in Q:

  * per-phase volume fractions of a tree,
  * Wiener (series/parallel) bounds and the 2D multiphase Hashin-Shtrikman
    bounds in their comparison-medium form
        sum_i f_i / (sigma_i + sigma_c) = 1 / (sigma_HS + sigma_c),
    with sigma_c the smallest (lower bound) or largest (upper bound) phase
    conductivity,
  * a Keller-Dykhne duality check: homogenization commutes with pointwise
    inversion-plus-90deg-rotation, so the same tree with phases 1/sigma_i must
    have effective tensor R (sigma*)^{-1} R^T -- an exact identity that any
    correct implementation must reproduce on every tree.

Tree JSON:  {"phase": "p1"}  or
            {"fraction": "2/5", "normal": [3, 1], "layers": [T1, T2]}
(fraction refers to layers[0]); phases are given as {"p1": "1", "p2": "5/2"}.

Usage:
    python laminate.py TREE.json --phases PHASES.json [--json OUT]
    python laminate.py --selftest

Standard library only.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction

Tensor = tuple[Fraction, Fraction, Fraction]  # [[a, b], [b, c]]


def iso(s: Fraction) -> Tensor:
    return (s, Fraction(0), s)


def is_spd(t: Tensor) -> bool:
    a, b, c = t
    return a > 0 and a * c - b * b > 0


def t_add(x: Tensor, y: Tensor) -> Tensor:
    return (x[0] + y[0], x[1] + y[1], x[2] + y[2])


def t_scale(x: Tensor, s: Fraction) -> Tensor:
    return (x[0] * s, x[1] * s, x[2] * s)


def t_inv(x: Tensor) -> Tensor:
    a, b, c = x
    det = a * c - b * b
    return (c / det, -b / det, a / det)


def t_rot90_conj(x: Tensor) -> Tensor:
    """R x R^T for the 90-degree rotation R = [[0, -1], [1, 0]]."""
    a, b, c = x
    return (c, -b, a)


def laminate(A: Tensor, B: Tensor, m: Fraction, n: tuple[int, int]) -> Tensor:
    """Effective tensor of layering A (fraction m) with B (1 - m), normal n."""
    assert 0 < m < 1, "degenerate fraction"
    u, v = n
    assert (u, v) != (0, 0)
    dA = (A[0] - B[0], A[1] - B[1], A[2] - B[2])
    w = (dA[0] * u + dA[1] * v, dA[1] * u + dA[2] * v)  # (A - B) n
    comp = t_add(t_scale(A, 1 - m), t_scale(B, m))
    denom = comp[0] * u * u + 2 * comp[1] * u * v + comp[2] * v * v  # n . comp n
    assert denom != 0
    coeff = m * (1 - m) / denom
    mean = t_add(t_scale(A, m), t_scale(B, 1 - m))
    return (
        mean[0] - coeff * w[0] * w[0],
        mean[1] - coeff * w[0] * w[1],
        mean[2] - coeff * w[1] * w[1],
    )


# ---------------------------------------------------------------------------
# trees
# ---------------------------------------------------------------------------


def effective(tree: dict, phases: dict[str, Fraction]) -> Tensor:
    if "phase" in tree:
        return iso(phases[tree["phase"]])
    m = Fraction(tree["fraction"])
    u, v = tree["normal"]
    A = effective(tree["layers"][0], phases)
    B = effective(tree["layers"][1], phases)
    return laminate(A, B, m, (int(u), int(v)))


def fractions_of(tree: dict) -> dict[str, Fraction]:
    if "phase" in tree:
        return {tree["phase"]: Fraction(1)}
    m = Fraction(tree["fraction"])
    out: dict[str, Fraction] = {}
    for sub, w in ((tree["layers"][0], m), (tree["layers"][1], 1 - m)):
        for k, v in fractions_of(sub).items():
            out[k] = out.get(k, Fraction(0)) + w * v
    return out


def keller_check(tree: dict, phases: dict[str, Fraction]) -> bool:
    """Exact Keller-Dykhne identity: dual phases give R (sigma*)^{-1} R^T."""
    sigma = effective(tree, phases)
    dual_phases = {k: 1 / v for k, v in phases.items()}
    dual_sigma = effective(tree, dual_phases)
    return dual_sigma == t_rot90_conj(t_inv(sigma))


# ---------------------------------------------------------------------------
# bounds (isotropic phases, prescribed fractions)
# ---------------------------------------------------------------------------


def wiener_bounds(fracs: list[Fraction], sigs: list[Fraction]) -> tuple[Fraction, Fraction]:
    harm = 1 / sum(f / s for f, s in zip(fracs, sigs))
    arith = sum(f * s for f, s in zip(fracs, sigs))
    return harm, arith


def hs_bounds(fracs: list[Fraction], sigs: list[Fraction]) -> tuple[Fraction, Fraction]:
    """2D multiphase Hashin-Shtrikman bounds, comparison-medium form."""
    lo_c, hi_c = min(sigs), max(sigs)
    lo = 1 / sum(f / (s + lo_c) for f, s in zip(fracs, sigs)) - lo_c
    hi = 1 / sum(f / (s + hi_c) for f, s in zip(fracs, sigs)) - hi_c
    return lo, hi


def eig_range_within(t: Tensor, lo: Fraction, hi: Fraction) -> bool:
    """Both eigenvalues of t in [lo, hi], decided exactly (PSD shift tests)."""
    a, b, c = t

    def psd(a_: Fraction, c_: Fraction) -> bool:
        return a_ >= 0 and c_ >= 0 and a_ * c_ - b * b >= 0

    return psd(a - lo, c - lo) and psd(hi - a, hi - c)


# ---------------------------------------------------------------------------
# self-test (exact throughout)
# ---------------------------------------------------------------------------


def selftest() -> None:
    F = Fraction
    s1, s2, s3 = F(1), F(3), F(7)

    # rank-1 closed form: normal e1 gives diag(harmonic, arithmetic)
    m = F(2, 5)
    t = laminate(iso(s1), iso(s2), m, (1, 0))
    harm = 1 / (m / s1 + (1 - m) / s2)
    arith = m * s1 + (1 - m) * s2
    assert t == (harm, F(0), arith), t

    # scale-invariance in the direction vector
    assert laminate(iso(s1), iso(s2), m, (2, 3)) == laminate(iso(s1), iso(s2), m, (-4, -6))

    # Keller-Dykhne duality on nested three-phase trees, several shapes
    trees = [
        {
            "fraction": "1/3",
            "normal": [1, 0],
            "layers": [
                {"phase": "p1"},
                {
                    "fraction": "3/7",
                    "normal": [1, 2],
                    "layers": [{"phase": "p2"}, {"phase": "p3"}],
                },
            ],
        },
        {
            "fraction": "2/7",
            "normal": [0, 1],
            "layers": [
                {
                    "fraction": "1/2",
                    "normal": [3, -1],
                    "layers": [{"phase": "p3"}, {"phase": "p1"}],
                },
                {
                    "fraction": "5/9",
                    "normal": [2, 5],
                    "layers": [{"phase": "p2"}, {"phase": "p1"}],
                },
            ],
        },
    ]
    phases = {"p1": s1, "p2": s2, "p3": s3}
    for tree in trees:
        assert keller_check(tree, phases), "Keller-Dykhne identity failed"
        fr = fractions_of(tree)
        assert sum(fr.values()) == 1

    # two-phase HS lower bound: comparison-medium form == textbook closed form
    for f1 in (F(1, 4), F(1, 2), F(9, 10)):
        lo, _ = hs_bounds([f1, 1 - f1], [s1, s2])
        closed = s1 + (1 - f1) / (1 / (s2 - s1) + f1 / (2 * s1))
        assert lo == closed, (lo, closed)

    # every laminate obeys the Wiener box (eigenvalue test is exact)
    for tree in trees:
        fr = fractions_of(tree)
        names = sorted(fr)
        w_lo, w_hi = wiener_bounds([fr[n] for n in names], [phases[n] for n in names])
        t = effective(tree, phases)
        assert is_spd(t)
        assert eig_range_within(t, w_lo, w_hi), "Wiener bounds violated"

    # three-phase HS bounds bracket a crude isotropic-ish laminate (sanity)
    fr3 = [F(1, 3), F(1, 3), F(1, 3)]
    hs_lo, hs_hi = hs_bounds(fr3, [s1, s2, s3])
    w_lo, w_hi = wiener_bounds(fr3, [s1, s2, s3])
    assert w_lo < hs_lo < hs_hi < w_hi, "HS must tighten Wiener strictly here"

    print("selftest: all checks passed")


# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser(description="Exact hierarchical-laminate effective tensors (2D)")
    ap.add_argument("tree", nargs="?", help="JSON file with the laminate tree")
    ap.add_argument("--phases", help="JSON file mapping phase names to rational conductivities")
    ap.add_argument("--json", metavar="PATH", help="write the full exact record here")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    if args.selftest:
        selftest()
        return
    if not args.tree or not args.phases:
        ap.error("need TREE and --phases (or --selftest)")

    with open(args.tree, encoding="utf-8") as fh:
        tree = json.load(fh)
    with open(args.phases, encoding="utf-8") as fh:
        phases = {k: Fraction(v) for k, v in json.load(fh).items()}

    t = effective(tree, phases)
    fr = fractions_of(tree)
    names = sorted(fr)
    fracs = [fr[n] for n in names]
    sigs = [phases[n] for n in names]
    w = wiener_bounds(fracs, sigs)
    hs = hs_bounds(fracs, sigs)
    duality = keller_check(tree, phases)

    print(f"effective tensor  [[{t[0]}, {t[1]}], [{t[1]}, {t[2]}]]")
    print("fractions        " + ", ".join(f"{n}={fr[n]}" for n in names))
    print(f"Wiener bounds     [{w[0]}, {w[1]}]")
    print(f"HS bounds         [{hs[0]}, {hs[1]}]  (isotropic-composite bounds)")
    print(f"Keller duality    {'ok' if duality else 'FAILED'}")
    assert duality, "Keller-Dykhne identity failed on this tree"

    if args.json:
        record = {
            "tree": tree,
            "phases": {k: str(v) for k, v in phases.items()},
            "tensor": [str(x) for x in t],
            "fractions": {k: str(v) for k, v in fr.items()},
            "wiener": [str(x) for x in w],
            "hs": [str(x) for x in hs],
        }
        with open(args.json, "w", encoding="utf-8") as fh:
            json.dump(record, fh, indent=1)
        print(f"record written to {args.json}")


if __name__ == "__main__":
    main()
