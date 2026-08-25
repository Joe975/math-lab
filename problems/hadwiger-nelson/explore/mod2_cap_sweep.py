"""Where does the mod-2 capping argument break?

For a group G = sum_j r_j E of rotated triangular lattices, define the colouring

    psi(sum_j r_j Z_j) = sum_j mu_j Z_j   in  E/2E,

for a fixed mu in (E/2E)^s.  E/2E has four elements, so psi is a 4-colouring,
and it is proper exactly when no unit vector of G lies in ker psi.  If such a mu
exists, *every* finite unit-distance graph on G is 4-colourable and no
construction inside G can ever reach five colours.

Two things have to be checked before that argument means anything:

* **G must be a direct sum.** If two coefficient tuples give the same point,
  psi is not a function of the point at all.  Injectivity is checked directly on
  the box used.
* **The unit vectors must be enumerated completely** for the box in question,
  which the exhaustive field search does by construction.

This script sweeps rotation sets and reports, for each, whether a capping mu
exists.  A group with no capping mu is a candidate arena for a 5-chromatic
graph -- it is where to search next, and the only place searching can pay.

Run:  python problems/hadwiger-nelson/explore/mod2_cap_sweep.py --max-n 27
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def _harness_dir():
    d = HERE
    while True:
        cand = os.path.join(d, "harness", "hadwiger-nelson")
        if os.path.isdir(cand):
            return cand
        parent = os.path.dirname(d)
        if parent == d:
            raise RuntimeError("could not locate harness/hadwiger-nelson")
        d = parent


sys.path.insert(0, _harness_dir())
sys.path.insert(0, HERE)

from homomorphism_colouring import eis_mul  # noqa: E402
from lattice import eisenstein_patch, norms_up_to  # noqa: E402
from multi_rot_group import MultiRotGroup  # noqa: E402

OUT = os.environ.get("MATHLAB_OUT", os.path.join(os.getcwd(), "out"))

F4 = [(0, 0), (1, 0), (0, 1), (1, 1)]


def mulf(p, q):
    r = eis_mul(p, q)
    return (r[0] % 2, r[1] % 2)


def addf(p, q):
    return ((p[0] + q[0]) % 2, (p[1] + q[1]) % 2)


def capping_functional(unit_vectors, s):
    """A mu in (E/2E)^s killing no unit vector, or None."""
    residues = {tuple((z[0] % 2, z[1] % 2) for z in v) for v in unit_vectors}
    for mu in itertools.product(F4, repeat=s):
        if all(mu_j == (0, 0) for mu_j in mu):
            continue
        good = True
        for v in residues:
            acc = (0, 0)
            for m, z in zip(mu, v):
                acc = addf(acc, mulf(m, z))
            if acc == (0, 0):
                good = False
                break
        if good:
            return mu, len(residues)
    return None, len(residues)


def check_direct_sum(group, patch_norm):
    """Is the coefficient tuple -> point map injective on this box?

    A collision means the group has a relation, and then the colouring is not
    well defined on points at all.
    """
    patch = eisenstein_patch(patch_norm)
    seen = {}
    collisions = 0

    def rec(prefix, j):
        nonlocal collisions
        if j == group.s:
            p = group.point(prefix)
            if p in seen:
                collisions += 1
            else:
                seen[p] = prefix
            return
        for ab in patch:
            rec(prefix + (ab,), j + 1)

    rec((), 0)
    return collisions, len(seen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-n", type=int, default=27)
    ap.add_argument("--patch", type=int, default=9)
    ap.add_argument("--injectivity-patch", type=int, default=3)
    ap.add_argument("--out", default=os.path.join(OUT, "mod2_cap_sweep.json"))
    args = ap.parse_args()

    norms = [n for n in norms_up_to(args.max_n) if n > 1]
    rows = []
    hdr = (f"{'rotations':>14} {'s':>3} {'|U|':>5} {'residues':>9} "
           f"{'inject':>7} {'capped by mu':>28} {'s':>7}")
    print(hdr)
    print("-" * len(hdr))
    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    specs = [([n], [(0,), (1,)]) for n in norms]
    specs += [([a, b], [(0, 0), (1, 0), (0, 1)])
              for a, b in itertools.combinations(norms, 2)
              if a <= 12 and b <= 21]

    for rot_norms, words in specs:
        t0 = time.time()
        try:
            group = MultiRotGroup(rot_norms, words)
        except Exception as exc:  # coprime-generator failures etc
            print(f"{str(rot_norms):>14}  skipped: {exc}")
            continue
        coll, distinct = check_direct_sum(group, args.injectivity_patch)
        uv = group.unit_vectors(args.patch)
        mu, nres = capping_functional(uv, group.s)
        row = {
            "rotations": rot_norms, "words": [list(w) for w in words],
            "s": group.s, "field": repr(group.field),
            "unit_vectors": len(uv), "mod2_residues": nres,
            "direct_sum_collisions": coll,
            "capping_mu": [list(m) for m in mu] if mu else None,
            "capped": mu is not None,
            "seconds": round(time.time() - t0, 1),
        }
        rows.append(row)
        print(f"{str(rot_norms):>14} {group.s:>3} {len(uv):>5} {nres:>9} "
              f"{'ok' if coll == 0 else f'{coll} COLL':>7} "
              f"{str(row['capping_mu']) if mu else 'NONE -- CANDIDATE ARENA':>28} "
              f"{row['seconds']:>7}", flush=True)
        with open(args.out, "w") as fh:
            json.dump(rows, fh, indent=2)

    uncapped = [r for r in rows if not r["capped"]]
    print(f"\n{len(rows)} groups tested, {len(uncapped)} with no capping "
          f"functional")
    for r in uncapped:
        print(f"  candidate arena: rotations {r['rotations']} "
              f"({r['unit_vectors']} unit vectors)")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
