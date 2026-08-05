#!/usr/bin/env python3
"""Independent escape probes for the skeptic review of 001.

Two things distinguish this from 001's classify.py / deep_probe.py probes:

1. **Consistent resolution.** 001's deep probe compared candidate escape
   values computed at 128 boundary samples against a baseline computed at
   512.  A coarser boundary sampling UNDER-estimates max_{dW}|p| and so
   OVER-estimates the ratio, which can manufacture apparent ascent.  Here
   every accepted ascent is re-evaluated at 512 samples against a 512-sample
   baseline before it counts.

2. **Different scheme, different RNG.**  Sparse coordinate-subset Gaussian
   perturbations with log-uniform scales (ridges move few coordinates), plus
   jittered Nelder-Mead restarts, seed 20260804 (001 used 999 / 12345 with
   dense perturbations at fixed scales).

Also re-evaluates each claimed ascent with a second, independently written
float ratio evaluator (Jacobi eigenvalues rather than the census's trig-cubic
closed form) to guard against an evaluator artifact.

Usage:
  skeptic_probe.py stalls SEED [SEED ...]     probe the 4 deep-probed stalls
  skeptic_probe.py maxima SEED [SEED ...]     attack claimed ratio-1 maxima

Float-only exploration tool; verdicts here feed the review record, claims are
certified elsewhere.  Standard library only.
"""

import cmath
import json
import math
import os
import random
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from census import objective, unpack, nelder_mead, compass  # noqa: E402

DATA = os.path.join(HERE, "..", "data")


# ----------------------------------------------------------------------
# independent float ratio evaluator: Jacobi rotations, not the trig cubic
# ----------------------------------------------------------------------

def jacobi_eigs(h, sweeps=30):
    """Eigenvalues of a Hermitian 3x3 (python complex rows) by cyclic complex
    Jacobi rotations.  Shares nothing with census.herm_eigs."""
    a = [row[:] for row in h]
    n = 3
    v = None
    for _ in range(sweeps):
        off = sum(abs(a[i][j]) ** 2 for i in range(n) for j in range(n) if i != j)
        if off < 1e-26:
            break
        for p in range(n):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-18:
                    continue
                app = a[p][p].real
                aqq = a[q][q].real
                # complex Jacobi: phase then angle
                phase = apq / abs(apq)
                tau = (aqq - app) / (2 * abs(apq))
                t = (1 if tau >= 0 else -1) / (abs(tau) + math.sqrt(1 + tau * tau))
                c = 1 / math.sqrt(1 + t * t)
                s = t * c
                for k in range(n):
                    akp = a[k][p]
                    akq = a[k][q]
                    a[k][p] = c * akp - s * (phase.conjugate() * akq)
                    a[k][q] = s * (phase * akp) + c * akq
                for k in range(n):
                    apk = a[p][k]
                    aqk = a[q][k]
                    a[p][k] = c * apk - s * (phase * aqk)
                    a[q][k] = s * (phase.conjugate() * apk) + c * aqk
    return sorted(a[i][i].real for i in range(3))


def ratio_indep(x, ntheta=512):
    """Crouzeix ratio via Jacobi eigenvalues throughout."""
    a, p = unpack(x)
    np_ = math.sqrt(sum(abs(c) ** 2 for c in p))
    if np_ < 1e-12:
        return 0.0
    p = [c / np_ for c in p]
    # numerator
    m = [[0j] * 3 for _ in range(3)]
    power = [[1 + 0j if i == j else 0j for j in range(3)] for i in range(3)]
    for c in p:
        for i in range(3):
            for j in range(3):
                m[i][j] += c * power[i][j]
        power = [[sum(power[i][k] * a[k][j] for k in range(3))
                  for j in range(3)] for i in range(3)]
    g = [[sum(m[k][i].conjugate() * m[k][j] for k in range(3))
          for j in range(3)] for i in range(3)]
    num = math.sqrt(max(0.0, jacobi_eigs(g)[2]))
    # denominator: boundary support points via Jacobi top eigenvalue + a
    # fine tangential |p| sweep on the support values themselves is not
    # available here, so use the same boundary-point construction idea but
    # with Jacobi: for each angle take the top eigenvector by inverse power
    # iteration on (lam*I - H + eps).
    den = 0.0
    for k in range(ntheta):
        th = 2 * math.pi * k / ntheta
        e = cmath.exp(-1j * th)
        h = [[(e * a[i][j] + (e * a[j][i]).conjugate()) / 2 for j in range(3)]
             for i in range(3)]
        lam = jacobi_eigs(h)[2]
        vec = top_vec(h, lam)
        z = sum(vec[i].conjugate() * a[i][j] * vec[j]
                for i in range(3) for j in range(3))
        val = 0j
        for c in reversed(p):
            val = val * z + c
        den = max(den, abs(val))
    if den < 1e-14:
        return 0.0
    return num / den


def top_vec(h, lam, iters=3):
    """Inverse iteration on (H - lam I + eps I)^-1 via 3x3 solve."""
    eps = 1e-9 * (1 + abs(lam))
    m = [[h[i][j] - ((lam + eps) if i == j else 0) for j in range(3)]
         for i in range(3)]
    v = [1.0 + 0j, 0.7 + 0.2j, 0.4 - 0.5j]
    for _ in range(iters):
        v = solve3(m, v)
        nv = math.sqrt(sum(abs(x) ** 2 for x in v))
        if nv < 1e-300:
            break
        v = [x / nv for x in v]
    return v


def solve3(m, b):
    a = [row[:] + [bi] for row, bi in zip(m, b)]
    for c in range(3):
        piv = max(range(c, 3), key=lambda r: abs(a[r][c]))
        a[c], a[piv] = a[piv], a[c]
        if abs(a[c][c]) < 1e-300:
            a[c][c] = 1e-300
        for r in range(3):
            if r != c:
                f = a[r][c] / a[c][c]
                a[r] = [x - f * y for x, y in zip(a[r], a[c])]
    return [a[i][3] / a[i][i] for i in range(3)]


# ----------------------------------------------------------------------
# the probe
# ----------------------------------------------------------------------

def probe(x, rng, nprobe=900, nm_restarts=3, label=""):
    """Search for ascent near x; every candidate is confirmed at 512-vs-512
    resolution before it counts.  Returns (v512, best512, best_x)."""
    v512 = objective(x, 512)
    n = len(x)
    cands = []
    for _ in range(nprobe):
        s = 10 ** rng.uniform(-3.5, -0.7)
        k = rng.choice((1, 2, 3, 5, 8, n))
        idx = rng.sample(range(n), k)
        y = list(x)
        for i in idx:
            y[i] += rng.gauss(0, s)
        fy = objective(y, 128)
        if fy > v512 - 5e-4:
            cands.append((fy, y))
    for _ in range(nm_restarts):
        y0 = [xi + rng.gauss(0, 0.02) for xi in x]
        y, fy = nelder_mead(lambda z: objective(z, 128), y0, step=0.1,
                            maxit=1200)
        cands.append((fy, y))
    yc, vc = compass(lambda z: objective(z, 256), x, objective(x, 256),
                     step=0.004, maxit=1200)
    cands.append((vc, yc))
    cands.sort(key=lambda t: -t[0])
    best512, best_x = v512, x
    for fy, y in cands[:25]:
        f512 = objective(y, 512)
        if f512 > best512:
            best512, best_x = f512, y
    return v512, best512, best_x


def run(seeds, kind):
    recs = {r["seed"]: r for r in
            (json.loads(l) for l in open(os.path.join(DATA, "classified.jsonl")))}
    rng = random.Random(20260804)
    out_path = os.path.join(DATA, "skeptic_probe.jsonl")
    results = []
    with open(out_path, "a") as out:
        for seed in seeds:
            rec = recs[int(seed)]
            x = rec["x"]
            v512, best, bx = probe(x, rng, label=str(seed))
            ascent = best - v512
            # confirm any ascent with the independent evaluator
            ind_base = ind_best = None
            if ascent > 1e-4:
                ind_base = ratio_indep(x, 512)
                ind_best = ratio_indep(bx, 512)
            verdict = ("ASCENT" if ascent > 1e-3 else
                       "MARGINAL" if ascent > 1e-4 else "NO-ASCENT")
            row = {"kind": kind, "seed": rec["seed"], "family": rec["family"],
                   "class": rec["class"], "ratio512": v512,
                   "best512": best, "ascent": ascent,
                   "indep_base": ind_base, "indep_best": ind_best,
                   "verdict": verdict}
            results.append(row)
            out.write(json.dumps(row) + "\n")
            out.flush()
            print(f"{kind} seed={rec['seed']} {rec['class']:11s} "
                  f"v512={v512:.6f} best512={best:.6f} ascent={ascent:+.6f} "
                  f"{verdict}"
                  + (f" indep {ind_base:.6f}->{ind_best:.6f}"
                     if ind_base is not None else ""), flush=True)
    return results


if __name__ == "__main__":
    run(sys.argv[2:], sys.argv[1])
