#!/usr/bin/env python3
"""Geometric detail for probe-surviving census endpoints (float, exploration).

For each endpoint of a chosen class: eigenvalues of A and roots of p placed
relative to W(A) (inside / on the boundary / outside, with signed margin as a
fraction of diam W), |p|-peak angles on the boundary, and the active set:
which boundary points attain max |p| and whether ||p(A)|| is attained with a
multiple top singular value of p(A).
"""

import cmath
import json
import math
import sys

sys.path.insert(0, "/home/user/work-crouzeix/explore")
from census import (unpack, boundary_points, peval, herm_part, mat_scale,
                    herm_eigs, poly_mat, mat_mul, adjoint, fro)


def roots_dk(coeffs):
    """Durand-Kerner for a complex polynomial, ascending coefficients."""
    cs = list(coeffs)
    while cs and abs(cs[-1]) < 1e-10:
        cs.pop()
    n = len(cs) - 1
    if n < 1:
        return []
    lead = cs[-1]
    cs = [c / lead for c in cs]
    rs = [complex(0.4, 0.9) ** k for k in range(n)]
    for _ in range(300):
        new = []
        for i, r in enumerate(rs):
            num = 0j
            v = 0j
            for c in reversed(cs):
                v = v * r + c
            den = 1.0 + 0j
            for j, s in enumerate(rs):
                if j != i:
                    den *= (r - s)
            new.append(r - v / den if abs(den) > 1e-30 else r)
        shift = max(abs(a - b) for a, b in zip(rs, new))
        rs = new
        if shift < 1e-13:
            break
    return rs


def margin(z, a, ntheta=256):
    """max_theta [Re(e^{-i th} z) - h(th)] / diamW: negative inside W(A)."""
    worst = -1e18
    for k in range(ntheta):
        th = 2 * math.pi * k / ntheta
        e = cmath.exp(-1j * th)
        h = herm_part(mat_scale(e, a))
        lam = herm_eigs(h)[2]
        worst = max(worst, (e * z).real - lam)
    return worst


def detail(rec):
    a, p = unpack(rec["x"])
    np_ = math.sqrt(sum(abs(c) ** 2 for c in p))
    p = [c / np_ for c in p]
    pts = boundary_points(a, 512)
    diam = max(abs(u - v) for u in pts for v in pts)
    evs = [a[0][0], a[1][1], a[2][2]]
    rts = roots_dk(p)
    ev_m = [margin(z, a) / diam for z in evs]
    rt_m = [margin(z, a) / diam for z in rts]
    vals = [abs(peval(p, z)) for z in pts]
    vmax = max(vals)
    peak_idx = [k for k in range(len(pts))
                if vals[k] >= 0.999 * vmax
                and vals[k] >= vals[k - 1] and vals[k] >= vals[(k + 1) % len(pts)]]
    peak_angles = [round(360.0 * k / len(pts), 1) for k in peak_idx]
    # top singular-value multiplicity of p(A)
    g = mat_mul(adjoint(poly_mat(p, a)), poly_mat(p, a))
    s = herm_eigs(g)
    sv_gap = (s[2] - s[1]) / max(s[2], 1e-15)
    return {
        "family": rec["family"], "seed": rec["seed"],
        "ratio": rec["diag"]["ratio"],
        "deg_p": len(rts),
        "eig_margins": [round(m, 4) for m in ev_m],
        "root_margins": [round(m, 4) for m in sorted(rt_m)],
        "n_peaks": len(peak_idx), "peak_angles": peak_angles,
        "sv_gap": round(sv_gap, 6),
    }


def main():
    cls = sys.argv[2] if len(sys.argv) > 2 else "below2-max"
    recs = [json.loads(line) for line in open(sys.argv[1])]
    recs = [r for r in recs if r.get("class") == cls]
    recs.sort(key=lambda r: -r["diag"]["ratio"])
    for r in recs:
        print(json.dumps(detail(r)))


if __name__ == "__main__":
    main()
