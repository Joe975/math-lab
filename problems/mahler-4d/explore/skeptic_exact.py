#!/usr/bin/env python3
"""SKEPTIC exact volume-product implementation, written from scratch for the
adversarial review.  Shares no code with census.c or the tier-0 harness.

Differences from the harness (polytope.py):
  - normals by 4D/3D generalized cross product (cofactor expansion), not
    rref nullspace;
  - integer arithmetic throughout (points scaled to Z^d), Fractions only in
    the final accumulations;
  - 3D volumes by cone-from-centroid over supporting planes (harness uses
    origin-free signed divergence);
  - 2D base case by centroid fan over supporting lines (harness recurses
    to 1D max-min);
  - polar facets EITHER by brute-force hull of the polar vertex set
    ("bf" mode) or by the vertex<->facet duality ("dual" mode, checked);
    both modes cross-checked against each other on demand.

Usage:
  skeptic_exact.py mask HEX [HEX...]
  skeptic_exact.py file LIST.txt [OUT.txt]         (first field = hex mask)
  skeptic_exact.py sample EVALFILE N SEED [OUT]    (recompute st=ok lines)
  skeptic_exact.py crosscheck LIST.txt N SEED      (bf vs dual polar modes)
"""
import random
import sys
from fractions import Fraction
from itertools import combinations
from math import gcd

# ---------------------------------------------------------------- pairs

PAIRS = []
for c0 in (-1, 0, 1):
    for c1 in (-1, 0, 1):
        for c2 in (-1, 0, 1):
            for c3 in (-1, 0, 1):
                p = (c0, c1, c2, c3)
                nz = [x for x in p if x]
                if nz and nz[0] == 1:
                    PAIRS.append(p)
assert len(PAIRS) == 40


def mask_points(mask):
    pts = []
    for i in range(40):
        if mask >> i & 1:
            pts.append(PAIRS[i])
            pts.append(tuple(-x for x in PAIRS[i]))
    return pts


# ---------------------------------------------------- integer linear alg

def det3(a, b, c):
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def cross4(u, v, w):
    """Integer vector orthogonal to u, v, w in Z^4 (generalized cross
    product; zero iff u,v,w linearly dependent)."""
    m = (u, v, w)
    out = []
    for j in range(4):
        cols = [c for c in range(4) if c != j]
        d = det3([m[0][c] for c in cols], [m[1][c] for c in cols],
                 [m[2][c] for c in cols])
        out.append(d if j % 2 == 0 else -d)
    return tuple(out)


def cross3(u, v):
    return (u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0])


def ivec_primitive(a):
    g = 0
    for x in a:
        g = gcd(g, abs(x))
    return a if g <= 1 else tuple(x // g for x in a)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def sub(a, b):
    return tuple(x - y for x, y in zip(a, b))


# ---------------------------------------------------------------- hulls
# A "facet list" is [(a, b, on)] with a integer outward normal (primitive),
# b integer/Fraction offset, on = tuple of point indices with <a,p> == b.

def facets_bruteforce(pts, d):
    """All supporting hyperplanes spanned by d affinely independent points;
    exact, integer, complete by construction.  pts: integer tuples."""
    n = len(pts)
    seen = {}
    crossf = cross4 if d == 4 else (lambda u, v, w=None: cross3(u, v))
    for combo in combinations(range(n), d):
        base = pts[combo[0]]
        diffs = [sub(pts[i], base) for i in combo[1:]]
        a = cross4(*diffs) if d == 4 else cross3(*diffs)
        if not any(a):
            continue
        b = dot(a, base)
        pos = neg = False
        for p in pts:
            s = dot(a, p) - b
            if s > 0:
                pos = True
            elif s < 0:
                neg = True
            if pos and neg:
                break
        if pos and neg:
            continue
        if pos:
            a = tuple(-x for x in a)
            b = -b
        a = ivec_primitive(a)
        b = dot(a, base)
        if (a, b) in seen:
            continue
        on = tuple(i for i, p in enumerate(pts) if dot(a, p) == b)
        seen[(a, b)] = (a, b, on)
    return list(seen.values())


def area2(pts):
    """Area of conv(pts), pts integer tuples in Z^2, via centroid fan over
    supporting lines.  Returns Fraction."""
    n = len(pts)
    if n < 3:
        return Fraction(0)
    # centroid (of the point set), scaled by n to stay integer
    cx = sum(p[0] for p in pts)
    cy = sum(p[1] for p in pts)
    # supporting lines through pairs of points
    seen = set()
    total = Fraction(0)
    for i, j in combinations(range(n), 2):
        u, v = pts[i], pts[j]
        a = (v[1] - u[1], u[0] - v[0])
        if a == (0, 0):
            continue
        b = dot(a, u)
        pos = neg = False
        for p in pts:
            s = dot(a, p) - b
            if s > 0:
                pos = True
            elif s < 0:
                neg = True
            if pos and neg:
                break
        if pos and neg:
            continue
        if pos:
            a, b = (-a[0], -a[1]), -b
        a = ivec_primitive(a)
        b_ = dot(a, pts[i])
        if (a, b_) in seen:
            continue
        seen.add((a, b_))
        on = [p for p in pts if dot(a, p) == b_]
        # edge endpoints: extremes along the line direction
        dvec = (-a[1], a[0])
        on.sort(key=lambda p: dot(dvec, p))
        e0, e1 = on[0], on[-1]
        # triangle (centroid, e0, e1): centroid = (cx/n, cy/n)
        tw = (e0[0] * n - cx) * (e1[1] * n - cy) \
            - (e0[1] * n - cy) * (e1[0] * n - cx)
        total += Fraction(abs(tw), 2 * n * n)
    return total


def vol3(pts):
    """Volume of conv(pts), pts integer tuples in Z^3, cone from centroid
    over the supporting planes."""
    n = len(pts)
    facs = facets_bruteforce(pts, 3)
    cx = tuple(sum(p[c] for p in pts) for c in range(3))  # centroid * n
    total = Fraction(0)
    for a, b, on in facs:
        j = next(c for c in range(3) if a[c])
        face2 = [tuple(x for c, x in enumerate(pts[i]) if c != j)
                 for i in on]
        ar = area2(face2)          # area of projected face
        if ar == 0:
            continue
        # height factor: |b - <a, centroid>| ; area scaling |a_j|
        h = abs(Fraction(b * n - dot(a, cx), n))
        total += h * ar / abs(a[j])
    return total / 3


def vol4(pts, facs=None):
    """Volume of conv(pts), pts integer tuples in Z^4, by divergence over
    facets with 3D volumes from vol3 (centroid-cone)."""
    if facs is None:
        facs = facets_bruteforce(pts, 4)
    total = Fraction(0)
    for a, b, on in facs:
        j = next(c for c in range(4) if a[c])
        face3 = [tuple(x for c, x in enumerate(pts[i]) if c != j)
                 for i in on]
        v = vol3(face3)
        total += Fraction(b) * v / abs(a[j])
    return total / 4


def rank_int(rows, d):
    """Rank over Q of integer rows (fraction-free Gaussian elimination)."""
    m = [list(r) for r in rows]
    r = 0
    for c in range(d):
        piv = next((i for i in range(r, len(m)) if m[i][c]), None)
        if piv is None:
            continue
        m[r], m[piv] = m[piv], m[r]
        for i in range(r + 1, len(m)):
            if m[i][c]:
                f = m[i][c]
                g = m[r][c]
                m[i] = [g * x - f * y for x, y in zip(m[i], m[r])]
        r += 1
    return r


def lcm(a, b):
    return a * b // gcd(a, b)


def exact_product(mask, polar_mode="auto"):
    """Returns dict with exact vol, pvol, P (Fractions), nv, nf, status."""
    pts = mask_points(mask)
    if rank_int(pts, 4) < 4:
        return {"status": "degen"}
    facs = facets_bruteforce(pts, 4)
    assert facs, "no facets?!"
    for a, b, on in facs:
        assert b > 0, f"origin not interior for {mask:x}"
    # vertex test: point i is a vertex iff its tight normals span R^4
    tight = {i: [] for i in range(len(pts))}
    for a, b, on in facs:
        for i in on:
            tight[i].append(a)
    nvert = sum(1 for i in range(len(pts)) if rank_int(tight[i], 4) == 4)
    status = "ok" if nvert == len(pts) else "improper"
    vol = vol4(pts, facs)
    # polar vertex set: a_F / b_F, scaled by L = lcm(b_F)
    L = 1
    for a, b, on in facs:
        L = lcm(L, b)
    ppts = [tuple(x * (L // b) for x in a) for a, b, on in facs]
    assert len(set(ppts)) == len(ppts)
    # guard: every polar point satisfies <v,y> <= 1 for ALL v (i.e. <= L)
    for y in ppts:
        for p in pts:
            assert dot(p, y) <= L, f"polar point outside K-polar, {mask:x}"
    if polar_mode == "auto":
        polar_mode = "bf" if len(ppts) <= 26 else "dual"
    if polar_mode == "bf":
        pfacs = facets_bruteforce(ppts, 4)
    else:  # dual: one facet of K-polar per VERTEX of K, offset L
        pfacs = []
        for i, v in enumerate(pts):
            if rank_int(tight[i], 4) < 4:
                continue  # not a vertex of K: redundant constraint, no facet
            on = tuple(k for k, y in enumerate(ppts) if dot(v, y) == L)
            assert len(on) >= 4, f"dual facet too small, {mask:x}"
            assert rank_int([sub(ppts[k], ppts[on[0]]) for k in on[1:]], 4) \
                == 3, f"dual facet not 3-dimensional, {mask:x}"
            pfacs.append((v, L, on))
    pvol_scaled = vol4(ppts, pfacs)
    pvol = pvol_scaled / Fraction(L) ** 4
    return {"status": status, "nv": nvert, "nf": len(facs),
            "vol": vol, "pvol": pvol, "P": vol * pvol,
            "npolarpts": len(ppts)}


def fmt(q):
    return f"{q.numerator}/{q.denominator}" if q.denominator != 1 \
        else str(q.numerator)


def run_file(path, out=None):
    bound = Fraction(32, 3)
    lines_out = []
    n_eq = n_above = n_below = 0
    min_above = None
    with open(path) as f:
        rows = [ln.split() for ln in f if ln.strip()]
    for i, fs in enumerate(rows):
        mask = int(fs[0], 16)
        r = exact_product(mask)
        if r["status"] == "degen":
            line = f"{fs[0]} DEGENERATE"
        else:
            P = r["P"]
            rel = "EQ" if P == bound else ("ABOVE" if P > bound else "BELOW")
            if rel == "EQ":
                n_eq += 1
            elif rel == "ABOVE":
                n_above += 1
                if min_above is None or P < min_above[0]:
                    min_above = (P, fs[0])
            else:
                n_below += 1
            line = (f"{fs[0]} st={r['status']} nv={r['nv']} nf={r['nf']} "
                    f"vol={fmt(r['vol'])} pvol={fmt(r['pvol'])} "
                    f"P={fmt(P)} {rel}")
        lines_out.append(line)
        if out and (i + 1) % 50 == 0:
            with open(out, "w") as g:
                g.write("\n".join(lines_out) + "\n")
            print(f"...{i + 1}/{len(rows)}", file=sys.stderr, flush=True)
    summ = (f"# SUMMARY {path}: {len(rows)} masks, {n_eq} == 32/3, "
            f"{n_above} above, {n_below} BELOW; min above = "
            f"{fmt(min_above[0]) if min_above else '-'} "
            f"(mask {min_above[1] if min_above else '-'})")
    lines_out.append(summ)
    if out:
        with open(out, "w") as g:
            g.write("\n".join(lines_out) + "\n")
    print(summ)


def run_sample(evalfile, n, seed, out=None):
    rng = random.Random(seed)
    with open(evalfile) as f:
        ok = [ln.split() for ln in f if " st=ok " in ln]
    chosen = ok if n >= len(ok) else rng.sample(ok, n)
    worst = Fraction(0)
    worst_mask = "-"
    lines = []
    for fs in chosen:
        mask = int(fs[0], 16)
        pf = float(fs[7][2:])  # P=...
        r = exact_product(mask)
        assert r["status"] == "ok", f"{fs[0]}: my status {r['status']} != ok"
        assert r["nv"] == int(fs[3][3:]), f"{fs[0]} nv mismatch"
        assert r["nf"] == int(fs[4][3:]), f"{fs[0]} nf mismatch"
        dev = abs(Fraction(pf) - r["P"])
        if dev > worst:
            worst, worst_mask = dev, fs[0]
        lines.append(f"{fs[0]} floatP={pf!r} exactP={fmt(r['P'])} "
                     f"dev={float(dev):.3e}")
    lines.append(f"# {evalfile}: {len(chosen)} sampled (seed={seed}), "
                 f"max |float-exact| = {float(worst):.3e} (mask {worst_mask})")
    print(lines[-1])
    if out:
        with open(out, "a") as g:
            g.write("\n".join(lines) + "\n")


def run_crosscheck(path, n, seed):
    rng = random.Random(seed)
    with open(path) as f:
        masks = [ln.split()[0] for ln in f if ln.strip()]
    chosen = masks if n >= len(masks) else rng.sample(masks, n)
    for h in chosen:
        m = int(h, 16)
        a = exact_product(m, polar_mode="bf")
        b = exact_product(m, polar_mode="dual")
        tag = "AGREE" if (a["P"] == b["P"] and a["pvol"] == b["pvol"]) \
            else "DISAGREE"
        print(f"{h}: bf P={fmt(a['P'])} dual P={fmt(b['P'])} {tag}")
        assert tag == "AGREE"
    print(f"# crosscheck {path}: {len(chosen)} masks, bf==dual everywhere")


if __name__ == "__main__":
    cmd = sys.argv[1]
    if cmd == "mask":
        for h in sys.argv[2:]:
            r = exact_product(int(h, 16))
            if r["status"] == "degen":
                print(f"{h}: DEGENERATE")
            else:
                print(f"{h}: st={r['status']} nv={r['nv']} nf={r['nf']} "
                      f"vol={fmt(r['vol'])} pvol={fmt(r['pvol'])} "
                      f"P={fmt(r['P'])} = {float(r['P']):.10f}")
    elif cmd == "file":
        run_file(sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    elif cmd == "sample":
        run_sample(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]),
                   sys.argv[5] if len(sys.argv) > 5 else None)
    elif cmd == "crosscheck":
        run_crosscheck(sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
    else:
        print(__doc__)
        sys.exit(2)
