#!/usr/bin/env python3
"""Exact coated-laminate constructions for 2D conductivity (two and three phases).

Everything here is over Q.  The two building blocks:

  * rank-1 coating identity.  Laminating a core tensor T (fraction c) with an
    isotropic matrix s I in unit normal n gives
        (sigma' - sI)^{-1} = (1/c) (T - sI)^{-1} + ((1-c)/(c s)) n n^T.
    Checked exactly against harness laminate() on every use (check_coating).
  * iterating it in the two orthogonal normals e1, e2 with core fractions
    c1, c2 (total core fraction F = c1 c2) gives
        (sigma* - sI)^{-1} = (1/F)(T - sI)^{-1} + ((1-F)/(F s)) diag(rho, 1-rho),
    rho = (1-c1)/(1-F).  Isotropy is LINEAR in rho, so the isotropic rank-2
    coated laminate around ANY diagonal core is a rational object.

Constructions:
  two_phase(s_core, s_mat, f_core)     rank-2 coated laminate, isotropic
  three_phase_lower(f, sigma)          core = rank-1 laminate of phases 2|3,
                                       coated by phase 1 (the worst) -> a
                                       candidate for the HS lower bound
  three_phase_upper(f, sigma)          the same with phase 3 as the coat

Each returns (tree, phases) ready for harness/three-phase-conductivity.

Usage:
    python tp_coated.py --selftest
    python tp_coated.py --scan          prints the attainability boundary
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

E1, E2 = (1, 0), (0, 1)


def leaf(p: str) -> dict:
    return {"phase": p}


def node(frac: Fr, normal, a: dict, b: dict) -> dict:
    return {"fraction": str(frac), "normal": list(normal), "layers": [a, b]}


def check_coating(T, s: Fr, c: Fr, n) -> None:
    """Exact check of the rank-1 coating identity against the harness."""
    got = L.laminate(T, L.iso(s), c, n)
    lhs = L.t_inv((got[0] - s, got[1], got[2] - s))
    u, v = n
    d = Fr(u * u + v * v)
    nn = (Fr(u * u) / d, Fr(u * v) / d, Fr(v * v) / d)
    base = L.t_scale(L.t_inv((T[0] - s, T[1], T[2] - s)), 1 / c)
    rhs = L.t_add(base, L.t_scale(nn, (1 - c) / (c * s)))
    assert lhs == rhs, ("coating identity fails", lhs, rhs)


def iso_coat(core_tree: dict, core_T, mat: str, s: Fr, F: Fr):
    """Coat a DIAGONAL core (fraction F) with matrix phase `mat` (conductivity s)
    in normals e1 then e2 so that the result is isotropic.  Returns (tree, rho)
    or None if the required rho is outside [0, 1] (construction impossible)."""
    h, b, r = core_T
    assert b == 0, "core must be diagonal"
    assert h != s and r != s
    rho = (s * (1 / (r - s) - 1 / (h - s)) + (1 - F)) / (2 * (1 - F))
    if not (0 <= rho <= 1):
        return None, rho
    c1 = 1 - rho * (1 - F)
    c2 = F / c1
    t1 = core_tree if c1 == 1 else node(c1, E1, core_tree, leaf(mat))
    t2 = t1 if c2 == 1 else node(c2, E2, t1, leaf(mat))
    return t2, rho


def two_phase(s_core: Fr, s_mat: Fr, f_core: Fr):
    phases = {"p1": s_core, "p2": s_mat}
    tree, rho = iso_coat(leaf("p1"), L.iso(s_core), "p2", s_mat, f_core)
    assert tree is not None and rho == Fr(1, 2)
    return tree, phases


def three_phase(f, sigma, coat: int):
    """coat = 0 -> phase 1 is the matrix (lower-bound candidate);
    coat = 2 -> phase 3 is the matrix (upper-bound candidate)."""
    f1, f2, f3 = f
    s1, s2, s3 = sigma
    phases = {"p1": s1, "p2": s2, "p3": s3}
    names = ["p1", "p2", "p3"]
    mat = names[coat]
    core_names = [n for n in names if n != mat]
    fc = [f[names.index(n)] for n in core_names]
    F = fc[0] + fc[1]
    if F == 0 or f[coat] == 0:
        return None, None, None
    m = fc[0] / F
    core = node(m, E1, leaf(core_names[0]), leaf(core_names[1]))
    core_T = L.effective(core, phases)
    tree, rho = iso_coat(core, core_T, mat, sigma[coat], F)
    return tree, phases, rho


def record(tree: dict, phases: dict) -> dict:
    t = L.effective(tree, phases)
    fr = L.fractions_of(tree)
    names = sorted(fr)
    fracs = [fr[n] for n in names]
    sigs = [phases[n] for n in names]
    assert L.keller_check(tree, phases)
    return {
        "tree": tree,
        "phases": {k: str(v) for k, v in phases.items()},
        "tensor": [str(x) for x in t],
        "fractions": {k: str(v) for k, v in fr.items()},
        "wiener": [str(x) for x in L.wiener_bounds(fracs, sigs)],
        "hs": [str(x) for x in L.hs_bounds(fracs, sigs)],
    }


def is_iso(t) -> bool:
    return t[1] == 0 and t[0] == t[2]


def selftest() -> None:
    # coating identity, random-ish exact instances
    for T, s, c, n in [
        ((Fr(3), Fr(0), Fr(5)), Fr(1), Fr(2, 7), (1, 0)),
        ((Fr(3), Fr(1, 2), Fr(5)), Fr(7), Fr(4, 9), (2, -3)),
        ((Fr(1, 3), Fr(-1, 5), Fr(2)), Fr(5, 2), Fr(1, 2), (1, 1)),
    ]:
        check_coating(T, s, c, n)

    # two-phase: coated laminate hits the HS bound exactly, both sides
    for s_a, s_b in [(Fr(1), Fr(3)), (Fr(2, 7), Fr(11)), (Fr(1), Fr(1000))]:
        for f in [Fr(1, 10), Fr(1, 3), Fr(1, 2), Fr(7, 9), Fr(99, 100)]:
            # worse phase as matrix -> lower bound
            tree, ph = two_phase(s_b, s_a, 1 - f)  # core = better phase, fraction 1-f
            t = L.effective(tree, ph)
            lo, hi = L.hs_bounds([f, 1 - f], [s_a, s_b])
            assert is_iso(t) and t[0] == lo, (t, lo)
            # better phase as matrix -> upper bound
            tree, ph = two_phase(s_a, s_b, f)
            t = L.effective(tree, ph)
            assert is_iso(t) and t[0] == hi, (t, hi)
            assert L.keller_check(tree, ph)
    print("selftest: coating identity + two-phase HS attainment exact")


def scan(sigma, N: int) -> None:
    s1, s2, s3 = sigma
    print(f"sigma = {sigma}; grid f_i = k/{N}")
    print("f1 f2 f3 | lower: rho, attains HS_lo | upper: rho, attains HS_hi")
    for i in range(1, N):
        for j in range(1, N - i):
            k = N - i - j
            f = [Fr(i, N), Fr(j, N), Fr(k, N)]
            row = f"{i:2d} {j:2d} {k:2d} |"
            for coat, idx in ((0, 0), (2, 1)):
                tree, ph, rho = three_phase(f, sigma, coat)
                hs = L.hs_bounds(f, list(sigma))[idx]
                if tree is None:
                    row += f" rho={float(rho):+.3f} --      |"
                else:
                    t = L.effective(tree, ph)
                    ok = is_iso(t) and t[0] == hs
                    row += f" rho={float(rho):+.3f} {'HS' if ok else 'no'} {float(t[0]):.4f}/{float(hs):.4f} |"
            print(row)



# ---------------------------------------------------------------------------
# Milton-type assemblage: each core phase coated separately by the matrix,
# with core fractions phi_i chosen so every coated laminate has the SAME
# scalar conductivity v.  Then any lamination of them is trivially v.
# v solves 1/(v + s_m) = sum_i f_i/(s_i + s_m), i.e. v is the multiphase HS
# value with s_m as comparison medium; feasibility is phi_i in (0, 1].
# ---------------------------------------------------------------------------


def hs_value(f, sigma, s_m: Fr) -> Fr:
    return 1 / sum(fi / (si + s_m) for fi, si in zip(f, sigma)) - s_m


def milton(f, sigma, coat: int):
    """Returns (tree, phases, info).  tree is None when infeasible."""
    names = ["p1", "p2", "p3"]
    phases = {n: s for n, s in zip(names, sigma)}
    s_m = sigma[coat]
    v = hs_value(f, sigma, s_m)
    parts = []
    info = {"v": v, "phi": {}, "w": {}}
    for i in range(3):
        if i == coat or f[i] == 0:
            continue
        s_i = sigma[i]
        # two-phase HS with matrix s_m, core fraction phi:
        #   1/(v+s_m) = (1-phi)/(2 s_m) + phi/(s_i+s_m)
        phi = (1 / (v + s_m) - 1 / (2 * s_m)) / (1 / (s_i + s_m) - 1 / (2 * s_m))
        info["phi"][names[i]] = phi
        if not (0 < phi <= 1):
            return None, phases, info
        w = f[i] / phi
        info["w"][names[i]] = w
        if phi == 1:
            sub = leaf(names[i])
        else:
            sub, rho = iso_coat(leaf(names[i]), L.iso(s_i), names[coat], s_m, phi)
            assert sub is not None
        parts.append((w, sub))
    assert sum(w for w, _ in parts) == 1, parts
    if len(parts) == 1:
        tree = parts[0][1]
    else:
        tree = node(parts[0][0], (1, 1), parts[0][1], parts[1][1])
    t = L.effective(tree, phases)
    assert is_iso(t) and t[0] == v, (t, v)
    return tree, phases, info


def scan_milton(sigma, N: int) -> None:
    s1, s2, s3 = sigma
    print(f"Milton assemblage, sigma = {tuple(map(str, sigma))}; grid k/{N}")
    print("f1 f2 f3 | HS_lo   attained? | HS_hi   attained?")
    for i in range(1, N):
        for j in range(1, N - i):
            k = N - i - j
            f = [Fr(i, N), Fr(j, N), Fr(k, N)]
            lo, hi = L.hs_bounds(f, list(sigma))
            tl, _, il = milton(f, sigma, 0)
            tu, _, iu = milton(f, sigma, 2)
            print(f"{i:2d} {j:2d} {k:2d} | {float(lo):.4f} {'YES' if tl else 'no '} (v<=s2: {lo <= s2})"
                  f" | {float(hi):.4f} {'YES' if tu else 'no '} (v>=s2: {hi >= s2})")
            # the feasibility condition must coincide with v vs s2 exactly
            assert (tl is not None) == (lo <= s2)
            assert (tu is not None) == (hi >= s2)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--sigma", default="1,2,5")
    ap.add_argument("--N", type=int, default=8)
    a = ap.parse_args()
    sigma = tuple(Fr(x) for x in a.sigma.split(","))
    if a.selftest:
        selftest()
    if a.scan:
        scan(sigma, a.N)
        scan_milton(sigma, a.N)
    if a.scan:
        scan_milton(sigma, a.N)
