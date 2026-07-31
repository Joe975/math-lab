#!/usr/bin/env python3
"""Exact (Fraction) volume products for census masks, via the tier-0 harness.

The census kernel (census.c) filters in double precision with exact integer
combinatorics.  Nothing from it is claimable; this script recomputes bodies
EXACTLY in Q with harness/mahler-4d/polytope.py, which is an entirely
independent implementation (different language, different code paths for
facets/polarity/volume).

Pair indexing must match census.c: iterate (c0,c1,c2,c3), ci in {-1,0,1}
ascending, c0 outermost, keep points whose first nonzero coordinate is +1,
index in discovery order.  Verify with `census points` if in doubt.

Usage:
  exact_check.py mask HEX [HEX...]        exact P for the given masks
  exact_check.py file LIST.txt            exact P for hex masks listed in a file
                                          (first whitespace field of each line)
  exact_check.py sample EVAL.txt N        exact-recompute every Nth `st=ok` line
                                          of an eval file, compare with the
                                          kernel's float P (tolerance 1e-9)
  exact_check.py claim HEX OUT.json       write a claim record for
                                          verify_product.py --check

Output for mask/file: one line per mask:
  <hex> k=<k> P=<exact fraction> float=<...> <BOUND|ABOVE> delta=<P-32/3 exact>
Exit status 1 if any body came out BELOW 32/3 (escalation required).
"""
import itertools
import json
import os
import sys
from fractions import Fraction

HARNESS = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "..", "..", "..", "harness", "mahler-4d")
sys.path.insert(0, HARNESS)
import polytope
from polytope import Q

BOUND = Fraction(32, 3)


def pair_reps():
    out = []
    for p in itertools.product((-1, 0, 1), repeat=4):
        if p == (0, 0, 0, 0):
            continue
        if next(x for x in p if x) == 1:
            out.append(p)
    assert len(out) == 40
    return out


REPS = pair_reps()


def mask_points(mask):
    pts = []
    for i in range(40):
        if mask >> i & 1:
            p = tuple(Q(x) for x in REPS[i])
            pts.append(p)
            pts.append(tuple(-x for x in p))
    return pts


def exact_product(mask):
    pts = mask_points(mask)
    assert polytope.is_symmetric(pts)
    vol = polytope.volume(pts)
    pvol = polytope.volume(polytope.polar(pts))
    return vol, pvol, vol * pvol


def show(mask):
    k = bin(mask).count("1")
    vol, pvol, p = exact_product(mask)
    rel = "BOUND" if p == BOUND else ("BELOW-ESCALATE" if p < BOUND else "ABOVE")
    print(f"{mask:x} k={k} P={p} float={float(p):.12f} {rel} "
          f"delta={p - BOUND}")
    return p


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    cmd = sys.argv[1]
    below = False
    if cmd == "mask":
        for h in sys.argv[2:]:
            below |= show(int(h, 16)) < BOUND
    elif cmd == "file":
        with open(sys.argv[2]) as f:
            for line in f:
                line = line.split()
                if line and not line[0].startswith("#"):
                    below |= show(int(line[0], 16)) < BOUND
    elif cmd == "sample":
        path, step = sys.argv[2], int(sys.argv[3])
        n = bad = 0
        with open(path) as f:
            oks = [ln for ln in f if " st=ok " in ln]
        for ln in oks[::step]:
            fields = dict(kv.split("=", 1) for kv in ln.split()[1:])
            mask = int(ln.split()[0], 16)
            _, _, p = exact_product(mask)
            err = abs(float(p) - float(fields["P"]))
            n += 1
            if err > 1e-9:
                bad += 1
                print(f"MISMATCH {mask:x} exact={float(p)} kernel={fields['P']}")
            below |= p < BOUND
        print(f"sample: {n} bodies exact-checked from {path}, "
              f"{bad} float mismatches > 1e-9")
        if bad:
            return 1
    elif cmd == "claim":
        mask = int(sys.argv[2], 16)
        vol, pvol, p = exact_product(mask)
        pts = mask_points(mask)
        rec = {
            "mask": f"{mask:x}",
            "vertices": [[str(x) for x in v] for v in pts],
            "volume": f"{vol.numerator}/{vol.denominator}",
            "volume_polar": f"{pvol.numerator}/{pvol.denominator}",
            "product": f"{p.numerator}/{p.denominator}",
        }
        with open(sys.argv[3], "w") as f:
            json.dump(rec, f, indent=1)
        print(f"claim for {mask:x} (P={p}) written to {sys.argv[3]}")
        below = p < BOUND
    else:
        print(__doc__)
        return 2
    if below:
        print("*** AT LEAST ONE BODY BELOW 32/3 -- ESCALATE, DO NOT RECORD "
              "WITHOUT INDEPENDENT RE-VERIFICATION ***")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
