#!/usr/bin/env python3
"""SKEPTIC check of the Hanner identification (claim 3).

Own Hanner-family generator (segments under l1/linf sums), own B4 canonical
form, own certificate verification: for each line of
data/hanner_match_attainers.txt either
  - "B4-equivalent ... mask M":   my_canon(mask) == my_canon(M), and M is in
    MY generated Hanner family;
  - "GL-equivalent ... M via T=[...]": T is invertible (det != 0) and maps
    the mask's 2k points bijectively onto MY vertex set of M.

Usage: skeptic_hanner.py CERTFILE [SAMPLE] [SEED]
"""
import itertools
import random
import re
import sys
from fractions import Fraction

from skeptic_burnside import PAIRS, PIDX, norm_pair, group_pair_perms
from skeptic_exact import exact_product


def canon(mask, perms):
    best = mask
    bits = [i for i in range(40) if mask >> i & 1]
    for g in perms:
        m = 0
        for i in bits:
            m |= 1 << g[i]
        if m < best:
            best = m
    return best


def mask_points(mask):
    pts = []
    for i in range(40):
        if mask >> i & 1:
            pts.append(PAIRS[i])
            pts.append(tuple(-x for x in PAIRS[i]))
    return pts


def points_mask(pts):
    m = 0
    for p in pts:
        m |= 1 << PIDX[norm_pair(tuple(int(x) for x in p))]
    return m


# --- my own Hanner generator -------------------------------------------

def hanner(n):
    """All Hanner vertex sets in R^n from the definition, as frozensets of
    integer tuples."""
    if n == 1:
        return {frozenset({(1,), (-1,)})}
    out = set()
    for k in range(1, n):
        for A in hanner(k):
            for B in hanner(n - k):
                # linf sum: cartesian product
                out.add(frozenset(a + b for a in A for b in B))
                # l1 sum: embed each in complementary coordinates
                z1, z2 = (0,) * (n - k), (0,) * k
                out.add(frozenset(a + z1 for a in A)
                        | frozenset(z2 + b for b in B))
    return out


def main():
    certfile = sys.argv[1]
    sample = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 777
    perms = sorted(set(group_pair_perms()))

    fam = sorted(hanner(4), key=lambda s: (len(s), sorted(s)))
    fam_masks = {points_mask(s): len(s) for s in fam}
    orbits = {}
    for m, nv in fam_masks.items():
        orbits.setdefault(canon(m, perms), []).append(m)
    print(f"my Hanner family in R^4: {len(fam)} vertex sets, "
          f"{len(orbits)} B4-orbits, vertex counts "
          f"{sorted(set(fam_masks.values()))}")
    # each is a genuine attainer of 32/3 by my own exact code
    for m in fam_masks:
        r = exact_product(m)
        assert r["P"] == Fraction(32, 3), f"Hanner mask {m:x} P={r['P']}!?"
    print("all my Hanner masks have exact P = 32/3 (my code)")

    glre = re.compile(
        r"^([0-9a-f]+): GL-equivalent to harness Hanner mask ([0-9a-f]+) "
        r"via T=\[(.*)\]$")
    b4re = re.compile(
        r"^([0-9a-f]+): B4-equivalent to harness Hanner "
        r"\((\d+) vertices, mask ([0-9a-f]+)\)$")

    with open(certfile) as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    gl_lines, b4_lines = [], []
    for ln in lines:
        if glre.match(ln):
            gl_lines.append(ln)
        elif b4re.match(ln):
            b4_lines.append(ln)
        else:
            print(f"UNPARSEABLE: {ln}")
            sys.exit(1)
    print(f"{len(lines)} certificates: {len(gl_lines)} GL, {len(b4_lines)} B4")

    # all B4 lines, verified with my canon
    for ln in b4_lines:
        h, nv, hm = b4re.match(ln).groups()
        mask, hmask = int(h, 16), int(hm, 16)
        assert hmask in fam_masks, f"{hm} not in MY Hanner family"
        assert fam_masks[hmask] == int(nv)
        assert canon(mask, perms) == canon(hmask, perms), f"B4 cert {h} FAILS"
    print(f"all {len(b4_lines)} B4 certificates verified (my canon)")

    chosen = gl_lines if not sample or sample >= len(gl_lines) \
        else random.Random(seed).sample(gl_lines, sample)
    targets = set()
    for ln in chosen:
        h, hm, tstr = glre.match(ln).groups()
        mask, hmask = int(h, 16), int(hm, 16)
        assert hmask in fam_masks, f"target {hm} not in MY Hanner family"
        targets.add(hmask)
        T = [[Fraction(x) for x in row.split(",")]
             for row in tstr.split("; ")]
        assert len(T) == 4 and all(len(r) == 4 for r in T)
        # det via Laplace on Fractions
        def det(m):
            if len(m) == 1:
                return m[0][0]
            return sum((-1) ** j * m[0][j]
                       * det([r[:j] + r[j + 1:] for r in m[1:]])
                       for j in range(len(m)))
        d = det(T)
        assert d != 0, f"cert {h}: T SINGULAR"
        V = mask_points(mask)
        W = set(mask_points(hmask))
        img = {tuple(sum(T[r][c] * p[c] for c in range(4)) for r in range(4))
               for p in V}
        assert len(img) == len(V), f"cert {h}: T not injective on V"
        assert img == W, f"cert {h}: T(V) != W"
    scope = f"all {len(chosen)}" if len(chosen) == len(gl_lines) \
        else f"sample {len(chosen)}/{len(gl_lines)} seed={seed}"
    print(f"GL certificates verified ({scope}): T invertible, "
          f"T(V) = W bijectively; {len(targets)} distinct targets")
    print("ALL CERTIFICATE CHECKS PASSED")


if __name__ == "__main__":
    main()
