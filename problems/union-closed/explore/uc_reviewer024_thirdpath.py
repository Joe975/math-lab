#!/usr/bin/env python3
"""Attempt 024 (fresh-session reviewer pass): third-path check of the 022 certificates.

Independent recomputation of the 022 spec's functionals, written from
swarm022_skeptic_brief.md alone. Exact Fractions for the coupling, census,
marginals and the Plackett root (bisection to 2^-70); float log2 only at
the very end. Values should land inside the certified enclosures (which
have width < 1e-14) to ~1e-12 absolute.
"""
import json, math, sys
from fractions import Fraction

def plackett_root(x, y, t):
    # unique root of F(z) = z(1-x-y+z) - t(x-z)(y-z) in (xy, min(x,y))
    F = lambda z: z * (1 - x - y + z) - t * (x - z) * (y - z)
    lo, hi = x * y, min(x, y)
    assert F(lo) < 0 < F(hi), "bracket invalid"
    for _ in range(70):
        mid = (lo + hi) / 2
        if F(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2

def check(w):
    n = w["n"]
    t = Fraction(*w["t"])
    u = {int(k, 2): Fraction(*v) for k, v in w["u"].items()}
    atoms = sorted(u)
    # coupling
    pi = {}
    Z = Fraction(0)
    for A in atoms:
        for B in atoms:
            p = u[A] * u[B] * t ** bin(A & B).count("1")
            pi[(A, B)] = p
            Z += p
    # measure + marginals (exact)
    mu = {A: sum(pi[(A, B)] for B in atoms) / Z for A in atoms}
    margs = [sum(m for A, m in mu.items() if A >> j & 1) for j in range(n)]
    in_regime = all(m < Fraction(38271, 100000) for m in margs)
    # census over i = 1..n-1
    num = 0.0
    agg = 0.0
    n_tables = 0
    mixed_zero = 0
    all_rplus = True
    for i in range(1, n):
        m_i = (1 << (i - 1)) - 1
        c_i = 1 << (i - 1)
        groups = {}
        for (A, B), p in pi.items():
            key = (A & m_i, B & m_i)
            cell = (1 if A & c_i else 0, 1 if B & c_i else 0)
            g = groups.setdefault(key, {})
            g[cell] = g.get(cell, Fraction(0)) + p / Z
        for g in groups.values():
            cells = [g.get(c, Fraction(0)) for c in ((0,0),(0,1),(1,0),(1,1))]
            pos = [c > 0 for c in cells]
            if all(pos):
                p00, p01, p10, p11 = cells
                mass = sum(cells)
                x = (p00 + p01) / mass
                y = (p00 + p10) / mass
                OR = (p00 * p11) / (p01 * p10)
                dev = math.log2(OR / t)
                z = plackett_root(x, y, t)
                if z >= Fraction(1, 2):
                    all_rplus = False
                hp = math.log2((1 - z) / (z))
                dz = float((x - z) * (y - z) / ((1 - x - y + 2 * z) + t * (x + y - 2 * z)))
                wgt = abs(hp * dz)
                num += float(mass) * wgt * dev
                agg += float(mass) * dev
                n_tables += 1
            elif any(pos):
                # degenerate table must have a full zero row or column
                r0 = cells[0] == 0 and cells[1] == 0
                r1 = cells[2] == 0 and cells[3] == 0
                c0 = cells[0] == 0 and cells[2] == 0
                c1 = cells[1] == 0 and cells[3] == 0
                if not (r0 or r1 or c0 or c1):
                    mixed_zero += 1
    return dict(name=w["name"], num=num, agg=agg, in_regime=in_regime,
                max_marginal=float(max(margs)), n_tables=n_tables,
                mixed_zero=mixed_zero,
                all_rplus=all_rplus if w["kind"] == "mm_abs_rplus" else None)

def main():
    data = json.load(open(sys.argv[1]))
    cert = {}
    for line in open(sys.argv[2]):
        line = line.strip()
        if line.startswith("{"):
            r = json.loads(line)
            cert[r["name"]] = r
    all_ok = True
    for w in data["witnesses"]:
        r = check(w)
        c = cert[r["name"]]
        clo = Fraction(c["num_lo"] or c["agg_lo"])
        chi = Fraction(c["num_hi"] or c["agg_hi"])
        mine = r["num"] if c["num_lo"] else r["agg"]
        inside = float(clo) - 1e-11 <= mine <= float(chi) + 1e-11
        match = (r["in_regime"] == c["in_regime"]
                 and r["n_tables"] == c["n_tables"]
                 and r["mixed_zero"] == c["mixed_zero_tables"]
                 and r["all_rplus"] == c["all_rplus"]
                 and abs(r["max_marginal"] - c["max_marginal"]) < 1e-12)
        ok = inside and match and mine < 0
        all_ok &= ok
        print(f"{r['name']:28s} mine={mine:+.10e} cert=[{float(clo):+.6e}] "
              f"inside={inside} meta_match={match} => {'OK' if ok else 'MISMATCH'}")
    print("ALL OK" if all_ok else "FAILURES PRESENT")

main()
