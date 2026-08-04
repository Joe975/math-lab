#!/usr/bin/env python3
"""Certified enclosures for census endpoints, via the tier-0 harness.

Each endpoint (A upper triangular 3x3 complex, p of degree <= 3) is rounded
entrywise to denominator DEN, and harness/crouzeix/ratio.py computes a
certified enclosure [ratio_lo, ratio_hi] of the Crouzeix ratio *of the rounded
point*, together with the exact-rational refutation test
||p(A)||^2 > 4 (max_{W} |p|)^2.  The claim certified is therefore about the
rounded rational configurations, which lie within 2^-1/DEN of the float
endpoints entrywise; the float values are exploration only.

Top endpoints additionally go through verify_ratio.verify, the independent
route (Faddeev-LeVerrier + Sturm), per the verification contract.

Usage:
  certify_endpoints.py classified.jsonl --den 256 --out certified.jsonl \
      [--top-verify 8] [--limit N]
"""

import argparse
import json
import sys
import time
from fractions import Fraction as Q

sys.path.insert(0, "/home/user/work-crouzeix/harness/crouzeix")
import ratio as R
import verify_ratio as V
sys.path.insert(0, "/home/user/work-crouzeix/explore")
from census import unpack


def to_rational(x, den):
    a, p = unpack(x)
    A = [[R.Cx(Q(round(a[i][j].real * den), den),
               Q(round(a[i][j].imag * den), den)) for j in range(3)]
         for i in range(3)]
    # normalize p as the census did (unit norm) BEFORE rounding, so the
    # rounded coefficients are O(1)
    import math
    nrm = math.sqrt(sum(abs(c) ** 2 for c in p))
    P = [R.Cx(Q(round(c.real / nrm * den), den),
              Q(round(c.imag / nrm * den), den)) for c in p]
    return A, P


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("infile")
    ap.add_argument("--den", type=int, default=256)
    ap.add_argument("--out", required=True)
    ap.add_argument("--top-verify", type=int, default=0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--tol", type=int, default=6, help="tol = 10^-this")
    ap.add_argument("--dirs", type=int, default=3)
    args = ap.parse_args()

    recs = [json.loads(line) for line in open(args.infile)]
    if args.limit:
        recs = recs[:args.limit]
    tol = Q(1, 10 ** args.tol)

    t_all = time.time()
    out = open(args.out, "w")
    n_ok = n_ref = 0
    for i, rec in enumerate(recs):
        A, P = to_rational(rec["x"], args.den)
        t0 = time.time()
        r = R.crouzeix_ratio(A, P, dir_bound=args.dirs, tol=tol)
        dt = time.time() - t0
        if r["degenerate"]:
            row = {"seed": rec["seed"], "family": rec["family"],
                   "class": rec.get("class"), "degenerate": True}
        else:
            if r["refutes"]:
                n_ref += 1
            # floor the exact rational at 12 decimals so the stored decimal
            # is still a certified lower bound
            lo_scaled = int(r["ratio_lo"] * 10 ** 12)   # exact floor
            lo_dec = f"{lo_scaled // 10**12}.{lo_scaled % 10**12:012d}"
            row = {"seed": rec["seed"], "family": rec["family"],
                   "class": rec.get("class"),
                   "float_ratio": rec["diag"]["ratio"],
                   "ratio_lo": lo_dec,
                   "ratio_hi_f": float(r["ratio_hi"]),
                   "refutes": r["refutes"], "den": args.den,
                   "dirs": args.dirs, "tol_exp": args.tol, "sec": round(dt, 2)}
            n_ok += 1
        out.write(json.dumps(row) + "\n")
        out.flush()
        if (i + 1) % 25 == 0:
            print(f"  {i+1}/{len(recs)} certified, "
                  f"{time.time()-t_all:.0f}s elapsed", file=sys.stderr)
    out.close()
    print(f"certified {n_ok} points, refutations: {n_ref}, "
          f"total {time.time()-t_all:.0f}s")

    if args.top_verify:
        recs.sort(key=lambda r: -r["diag"]["ratio"])
        print(f"\nindependent verification (verify_ratio.py) of the top "
              f"{args.top_verify} endpoints by float ratio:")
        for rec in recs[:args.top_verify]:
            A, P = to_rational(rec["x"], args.den)
            print(f"[{rec['family']} seed={rec['seed']} "
                  f"float={rec['diag']['ratio']:.5f}]")
            ok = V.verify(A, P, dir_bound=args.dirs, tol=tol)
            print("  verify:", "AGREE" if ok else "MISMATCH")


if __name__ == "__main__":
    main()
