#!/usr/bin/env python3
"""Float SCREEN: best isotropic effective conductivity of bounded-rank
hierarchical laminates at fixed phase fractions.  Nothing here certifies
anything; tp_certify.py turns a screened tree into an exact record.

Search space: binary trees with R internal nodes (rank R), leaves labelled by
phases (repeats allowed), normals from a fixed integer set, fractions
continuous.  For each topology the fractions are optimised by Nelder-Mead on
    objective = (+/-) trace(sigma*)/2  + W * anisotropy + W * fraction error
where anisotropy = (lambda_max - lambda_min) / trace and fraction error is the
L1 distance of the tree's per-phase fractions from the target.

Deterministic (seeded).  Standard library only.

Usage:
    python tp_search.py --f 1/8,1/8,6/8 --sigma 1,2,5 --rank 3 [--side lower|upper]
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from fractions import Fraction as Fr

NORMALS_ALL = [(1, 0), (0, 1), (1, 1), (1, -1)]


# ---- float laminate algebra (mirrors harness laminate.py) ------------------

def lam(A, B, m, n):
    u, v = n
    dA = (A[0] - B[0], A[1] - B[1], A[2] - B[2])
    w = (dA[0] * u + dA[1] * v, dA[1] * u + dA[2] * v)
    comp = ((1 - m) * A[0] + m * B[0], (1 - m) * A[1] + m * B[1], (1 - m) * A[2] + m * B[2])
    den = comp[0] * u * u + 2 * comp[1] * u * v + comp[2] * v * v
    k = m * (1 - m) / den
    return (m * A[0] + (1 - m) * B[0] - k * w[0] * w[0],
            m * A[1] + (1 - m) * B[1] - k * w[0] * w[1],
            m * A[2] + (1 - m) * B[2] - k * w[1] * w[1])


# ---- topologies -------------------------------------------------------------
# A shape is a nested tuple: leaf index (int) or (left, right).

def shapes(n_internal: int):
    """All binary tree shapes with n_internal internal nodes, leaves numbered."""
    def build(k):
        if k == 0:
            yield "L"
            return
        for i in range(k):
            for a in build(i):
                for b in build(k - 1 - i):
                    yield (a, b)
    out = []
    for s in build(n_internal):
        counter = itertools.count()

        def number(t):
            if t == "L":
                return next(counter)
            return (number(t[0]), number(t[1]))
        out.append(number(s))
    return out


def internal_nodes(shape):
    if isinstance(shape, int):
        return []
    return [shape] + internal_nodes(shape[0]) + internal_nodes(shape[1])


class Topology:
    def __init__(self, shape, labels, normals):
        self.shape, self.labels, self.normals = shape, labels, normals
        self.nodes = internal_nodes(shape)
        self.n_params = len(self.nodes)

    def evaluate(self, params, sig):
        """returns (tensor, fractions) in float."""
        it = iter(params)

        def rec(t):
            if isinstance(t, int):
                s = sig[self.labels[t]]
                fr = [0.0, 0.0, 0.0]
                fr[self.labels[t]] = 1.0
                return (s, 0.0, s), fr
            m = next(it)
            idx = self.nodes.index(t)
            A, fa = rec(t[0])
            B, fb = rec(t[1])
            fr = [m * x + (1 - m) * y for x, y in zip(fa, fb)]
            return lam(A, B, m, self.normals[idx]), fr
        return rec(self.shape)

    def to_tree(self, params_frac):
        it = iter(params_frac)

        def rec(t):
            if isinstance(t, int):
                return {"phase": f"p{self.labels[t] + 1}"}
            m = next(it)
            idx = self.nodes.index(t)
            return {"fraction": str(m), "normal": list(self.normals[idx]),
                    "layers": [rec(t[0]), rec(t[1])]}
        return rec(self.shape)


# ---- optimisation -----------------------------------------------------------

def sigmoid(x):
    return 1 / (1 + math.exp(-x)) if x > -700 else 0.0


def objective(top, x, sig, f_target, sign, W=200.0):
    params = [sigmoid(v) for v in x]
    t, fr = top.evaluate(params, sig)
    a, b, c = t
    tr = a + c
    disc = math.sqrt(max((a - c) ** 2 + 4 * b * b, 0.0))
    aniso = disc / tr
    ferr = sum(abs(p - q) for p, q in zip(fr, f_target))
    return sign * tr / 2 + W * (aniso + ferr), tr / 2, aniso, ferr


def nelder_mead(fn, x0, step=1.0, iters=400, tol=1e-12):
    n = len(x0)
    pts = [list(x0)]
    for i in range(n):
        p = list(x0)
        p[i] += step
        pts.append(p)
    vals = [fn(p)[0] for p in pts]
    for _ in range(iters):
        order = sorted(range(n + 1), key=lambda i: vals[i])
        pts = [pts[i] for i in order]
        vals = [vals[i] for i in order]
        if abs(vals[-1] - vals[0]) < tol:
            break
        cen = [sum(p[j] for p in pts[:-1]) / n for j in range(n)]
        worst = pts[-1]
        refl = [cen[j] + (cen[j] - worst[j]) for j in range(n)]
        fr_ = fn(refl)[0]
        if fr_ < vals[0]:
            exp = [cen[j] + 2 * (cen[j] - worst[j]) for j in range(n)]
            fe = fn(exp)[0]
            if fe < fr_:
                pts[-1], vals[-1] = exp, fe
            else:
                pts[-1], vals[-1] = refl, fr_
        elif fr_ < vals[-2]:
            pts[-1], vals[-1] = refl, fr_
        else:
            con = [cen[j] + 0.5 * (worst[j] - cen[j]) for j in range(n)]
            fc = fn(con)[0]
            if fc < vals[-1]:
                pts[-1], vals[-1] = con, fc
            else:
                best = pts[0]
                for i in range(1, n + 1):
                    pts[i] = [best[j] + 0.5 * (pts[i][j] - best[j]) for j in range(n)]
                    vals[i] = fn(pts[i])[0]
    i = min(range(n + 1), key=lambda i: vals[i])
    return pts[i], vals[i]


def search(f_target, sig, rank, side, normals, restarts, seed, tol_report):
    rng = random.Random(seed)
    sign = 1.0 if side == "lower" else -1.0
    best = []  # (value, aniso, ferr, topology, params)
    n_leaf = rank + 1
    count = 0
    required = set(i for i, f in enumerate(f_target) if f > 0)
    for shape in shapes(rank):
        for labels in itertools.product(range(3), repeat=n_leaf):
            if required - set(labels):
                continue  # a required phase is missing
            for norms in itertools.product(normals, repeat=rank):
                top = Topology(shape, labels, norms)
                count += 1

                def fn(x, top=top):
                    return objective(top, x, sig, f_target, sign)
                for _ in range(restarts):
                    x0 = [rng.uniform(-2, 2) for _ in range(top.n_params)]
                    x, val = nelder_mead(fn, x0)
                    # polish
                    x, val = nelder_mead(fn, x, step=0.05, iters=400)
                    _, mean, aniso, ferr = fn(x)
                    if aniso < tol_report and ferr < tol_report:
                        best.append((mean, aniso, ferr, top, [sigmoid(v) for v in x]))
    best.sort(key=lambda r: sign * r[0])
    return best, count


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--f", default="1/8,1/8,6/8")
    ap.add_argument("--sigma", default="1,2,5")
    ap.add_argument("--rank", type=int, default=3)
    ap.add_argument("--side", choices=["lower", "upper"], default="lower")
    ap.add_argument("--normals", default="all", help="all or axes")
    ap.add_argument("--restarts", type=int, default=3)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--top", type=int, default=5)
    ap.add_argument("--tol", type=float, default=1e-6)
    ap.add_argument("--out", help="write the best screened trees as JSON")
    a = ap.parse_args()
    f = [float(Fr(x)) for x in a.f.split(",")]
    sig = [float(Fr(x)) for x in a.sigma.split(",")]
    normals = NORMALS_ALL if a.normals == "all" else NORMALS_ALL[:2]
    best, count = search(f, sig, a.rank, a.side, normals, a.restarts, a.seed, a.tol)
    F = [Fr(x) for x in a.f.split(",")]
    S = [Fr(x) for x in a.sigma.split(",")]
    lo_c, hi_c = min(S), max(S)
    hs_lo = 1 / sum(p / (s + lo_c) for p, s in zip(F, S)) - lo_c
    hs_hi = 1 / sum(p / (s + hi_c) for p, s in zip(F, S)) - hi_c
    print(f"f={a.f} sigma={a.sigma} rank={a.rank} side={a.side} topologies={count}")
    print(f"HS_lo={float(hs_lo):.6f} HS_hi={float(hs_hi):.6f}")
    seen = set()
    shown = []
    for mean, aniso, ferr, top, params in best:
        key = (top.shape, top.labels, top.normals)
        if key in seen:
            continue
        seen.add(key)
        shown.append((mean, aniso, ferr, top, params))
        print(f"  {mean:.7f}  aniso={aniso:.1e} ferr={ferr:.1e}  labels={top.labels} normals={top.normals} shape={top.shape}")
        print("     fractions=" + ", ".join(f"{p:.6f}" for p in params))
        if len(shown) >= a.top:
            break
    if a.out:
        with open(a.out, "w", encoding="utf-8") as fh:
            json.dump([{"mean": m, "aniso": an, "ferr": fe,
                        "tree": t.to_tree([repr(p) for p in pr])}
                       for m, an, fe, t, pr in shown], fh, indent=1)


if __name__ == "__main__":
    main()
