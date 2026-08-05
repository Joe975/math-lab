#!/usr/bin/env python3
"""Independent structural checks on an equilibria.py certificate (skeptic pass).

Written for the skeptic review of attempt 001 (abk-witness). Shares no code
with harness/maxwell-equilibria/equilibria.py or verify_equilibria.py; every
check below is re-implemented from scratch in exact rational arithmetic
(fractions.Fraction + math.isqrt only).

Checks, and why each is stronger than (or independent of) what the two
harness engines already did:

  TREE   The verifier's tiling check is prefix-freeness + the Kraft equality
         sum 2^-depth = 1. That is NOT sufficient for a leaf set whose paths
         carry per-node axis choices: the set {"01", "11"} (axis-0 upper
         half, axis-1 upper half) is prefix-free with Kraft sum 1 and tiles
         nothing -- the two boxes overlap and jointly miss a quadrant. The
         check here reconstructs the subdivision tree from the paths and
         verifies AXIS CONSISTENCY at every internal node: each node is
         either a leaf (path ends there, no continuations) or ALL paths
         through it agree on the split axis and both children subtrees are
         present and well-formed. By induction this proves the leaf boxes
         partition the root box exactly. Kraft and prefix-freeness are then
         theorems, but are recomputed anyway as cross-checks.

  DISJ   The 24 isolation enclosures are pairwise disjoint, by exact
         rational comparison (some coordinate strictly separates each pair).

  BOX    Every isolated leaf's recorded box equals the box rebuilt from its
         split path (own midpoint-splitting code), its enclosure is contained
         in it, and count == number of isolated leaves.

  INDEX  All det signs are +-1 and sum to 1 - n (Poincare-Hopf identity for
         all-positive configurations), recomputed from the leaf data.

  REGION R0^2 = max |x_i|^2 recomputed from the config; the region must
         contain the closed localization ball (certified rational sqrt upper
         bound, own isqrt code).

  RHO    The domination inequality |q_i|/rho_i^2 > sum_j |q_j|/(d_ij-rho_i)^2
         is re-instantiated exactly (own certified lower bounds for d_ij).

  BALL   Every excluded-ball leaf's box lies in the closed ball of its
         domination radius (exact); every excluded-outside leaf's box has
         min |x|^2 >= R0^2 (exact).

  TALLY  Leaf-kind tallies match the certificate's summary fields.

Usage:  skeptic_cert_checks.py CERT.json [--centroid] [--expect-count N]
        --centroid additionally checks that (1/3,1/3,1/3) lies in exactly
        one isolation enclosure (the abk-witness family has a symmetric
        rational equilibrium there).

Exit status 0 iff every check passes. Standard library only.
"""

from __future__ import annotations

import argparse
import json
import sys
from fractions import Fraction
from math import isqrt

FAILURES = []


def check(name: str, cond: bool, detail: str = "") -> None:
    line = ("ok   " if cond else "FAIL ") + name + (f"  [{detail}]" if detail else "")
    print(line)
    if not cond:
        FAILURES.append(name)


def frac(s) -> Fraction:
    return Fraction(s)


def sqrt_up(f: Fraction, bits: int = 128) -> Fraction:
    """Certified upper bound for sqrt(f), f >= 0. Own implementation."""
    s = 1 << bits
    t = f.numerator * f.denominator * s * s
    r = isqrt(t)
    if r * r != t:
        r += 1
    return Fraction(r, f.denominator * s)


def sqrt_dn(f: Fraction, bits: int = 128) -> Fraction:
    s = 1 << bits
    return Fraction(isqrt(f.numerator * f.denominator * s * s), f.denominator * s)


# --------------------------------------------------------------------------
# TREE: full subdivision-tree reconstruction with axis consistency
# --------------------------------------------------------------------------

def check_tree(paths: list[str]) -> None:
    """Prove the leaf paths are exactly the leaves of one well-formed binary
    subdivision tree (every internal node splits along a single axis into
    both halves). This implies the leaf boxes tile the root box."""
    TERM = "$"
    root: dict = {}
    dup = badchar = 0
    for p in paths:
        if len(p) % 2:
            badchar += 1
            continue
        node = root
        ok = True
        for i in range(0, len(p), 2):
            pair = p[i:i + 2]
            if pair[0] not in "012" or pair[1] not in "01":
                ok = False
                break
            node = node.setdefault(pair, {})
        if not ok:
            badchar += 1
            continue
        if TERM in node:
            dup += 1
        node[TERM] = True
    check("TREE: all paths well-formed (even length, axis in 0-2, side in 0-1)",
          badchar == 0, f"{badchar} malformed")
    check("TREE: no duplicate leaf paths", dup == 0, f"{dup} duplicates")

    bad_mixed = bad_onechild = bad_termint = empty = 0
    n_leaves = 0
    max_depth = 0
    stack = [(root, 0)]
    while stack:
        node, depth = stack.pop()
        kids = [k for k in node if k != TERM]
        term = TERM in node
        if term and kids:
            bad_termint += 1      # a leaf path is a strict prefix of another
        if not term and not kids:
            empty += 1            # cannot happen structurally, defensive
        if term and not kids:
            n_leaves += 1
            max_depth = max(max_depth, depth)
            continue
        axes = {k[0] for k in kids}
        if len(axes) != 1:
            bad_mixed += 1        # two different split axes at one node
        else:
            a = axes.pop()
            if set(kids) != {a + "0", a + "1"}:
                bad_onechild += 1  # only one half present -> a hole
        for k in kids:
            stack.append((node[k], depth + 1))
    check("TREE: no leaf path is a prefix of another", bad_termint == 0,
          f"{bad_termint} violations")
    check("TREE: single split axis at every internal node", bad_mixed == 0,
          f"{bad_mixed} mixed-axis nodes")
    check("TREE: both children present at every internal node", bad_onechild == 0,
          f"{bad_onechild} one-child nodes")
    check("TREE: leaf count from tree equals path count",
          n_leaves == len(paths), f"{n_leaves} vs {len(paths)}")
    # cross-checks (implied by the above, but cheap)
    kraft = sum(Fraction(1, 1 << (len(p) // 2)) for p in paths)
    check("TREE: Kraft equality sum 2^-depth == 1", kraft == 1, str(kraft))
    print(f"       leaves={n_leaves} max_depth={max_depth} splits")


# --------------------------------------------------------------------------
# box rebuild (own code)
# --------------------------------------------------------------------------

def rebuild(root_box, path: str):
    box = [list(iv) for iv in root_box]
    for i in range(0, len(path), 2):
        a, s = int(path[i]), int(path[i + 1])
        lo, hi = box[a]
        mid = (lo + hi) / 2
        box[a] = [lo, mid] if s == 0 else [mid, hi]
    return [tuple(iv) for iv in box]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("cert")
    ap.add_argument("--centroid", action="store_true")
    ap.add_argument("--expect-count", type=int, default=None)
    args = ap.parse_args()

    with open(args.cert, encoding="utf-8") as fh:
        cert = json.load(fh)

    charges = [(frac(c["q"]), tuple(frac(v) for v in c["pos"]))
               for c in cert["config"]["charges"]]
    n = len(charges)
    root_box = [(frac(lo), frac(hi)) for lo, hi in cert["region"]]
    leaves = cert["leaves"]
    isolated = cert["isolated"]
    print(f"cert: n={n} charges, {len(leaves)} leaves, count={cert['count']}, "
          f"unresolved={cert['unresolved']}, complete={cert.get('complete')}")

    # ---- TREE ----
    check_tree([lf["path"] for lf in leaves])

    # ---- TALLY ----
    tally: dict = {}
    for lf in leaves:
        tally[lf["kind"]] = tally.get(lf["kind"], 0) + 1
    check("TALLY: isolated leaves == count field",
          tally.get("isolated", 0) == cert["count"],
          f"{tally.get('isolated', 0)} vs {cert['count']}")
    check("TALLY: unresolved leaves == unresolved field",
          tally.get("unresolved", 0) == cert["unresolved"],
          f"{tally.get('unresolved', 0)} vs {cert['unresolved']}")
    for kind, cnt in cert["excluded"].items():
        check(f"TALLY: {kind} == summary", tally.get(kind, 0) == cnt,
              f"{tally.get(kind, 0)} vs {cnt}")
    known = {"isolated", "unresolved"} | set(cert["excluded"])
    check("TALLY: no unknown leaf kinds", set(tally) <= known, str(set(tally) - known))
    if args.expect_count is not None:
        check(f"TALLY: count == expected {args.expect_count}",
              cert["count"] == args.expect_count, str(cert["count"]))

    # ---- BOX ----
    iso_paths = {lf["path"] for lf in leaves if lf["kind"] == "isolated"}
    check("BOX: isolated records match isolated leaves 1:1",
          iso_paths == {it["path"] for it in isolated}
          and len(isolated) == len(iso_paths))
    bad_box = bad_enc = 0
    encs = []
    for it in isolated:
        want = rebuild(root_box, it["path"])
        got = [tuple(frac(v) for v in iv) for iv in it["box"]]
        if want != got:
            bad_box += 1
        enc = [tuple(frac(v) for v in iv) for iv in it["enclosure"]]
        for k in range(3):
            if not (got[k][0] <= enc[k][0] <= enc[k][1] <= got[k][1]):
                bad_enc += 1
                break
        encs.append(enc)
    check("BOX: every recorded isolated box == box rebuilt from its path",
          bad_box == 0, f"{bad_box} mismatches")
    check("BOX: every enclosure inside its leaf box", bad_enc == 0,
          f"{bad_enc} violations")

    # ---- DISJ ----
    overlap = 0
    min_gap = None
    for i in range(len(encs)):
        for j in range(i + 1, len(encs)):
            seps = []
            for k in range(3):
                if encs[i][k][1] < encs[j][k][0]:
                    seps.append(encs[j][k][0] - encs[i][k][1])
                elif encs[j][k][1] < encs[i][k][0]:
                    seps.append(encs[i][k][0] - encs[j][k][1])
            if not seps:
                overlap += 1
            else:
                g = max(seps)
                min_gap = g if min_gap is None else min(min_gap, g)
    check("DISJ: isolation enclosures pairwise disjoint (exact)", overlap == 0,
          f"{overlap} overlapping pairs of {len(encs) * (len(encs) - 1) // 2}")
    if min_gap is not None:
        print(f"       minimum pairwise separating gap ~ {float(min_gap):.3e}")

    # ---- INDEX ----
    signs = [it["det_sign"] for it in isolated]
    check("INDEX: every det sign is +1 or -1", all(s in (1, -1) for s in signs))
    check(f"INDEX: sum of det signs == 1 - n == {1 - n}",
          sum(signs) == 1 - n,
          f"sum={sum(signs)}, +1s={signs.count(1)}, -1s={signs.count(-1)}")

    # ---- REGION ----
    r0_sq = max(p[0] ** 2 + p[1] ** 2 + p[2] ** 2 for _, p in charges)
    if all(q > 0 for q, _ in charges):
        r_up = sqrt_up(r0_sq)
        ok = all(root_box[k][0] <= -r_up and r_up <= root_box[k][1] for k in range(3))
        check("REGION: region contains the closed localization ball |x| <= R0", ok,
              f"R0^2 = {r0_sq}")
    else:
        print("REGION: charges not all positive; localization not applicable")

    # ---- RHO ----
    rhos = [frac(r) for r in cert["domination_rhos"]]
    ok = len(rhos) == n
    if ok:
        for i, (qi, pi) in enumerate(charges):
            rho = rhos[i]
            if rho <= 0:
                ok = False
                break
            rhs = Fraction(0)
            for j, (qj, pj) in enumerate(charges):
                if j == i:
                    continue
                d_lo = sqrt_dn(sum((pi[k] - pj[k]) ** 2 for k in range(3)))
                if d_lo <= rho:
                    ok = False
                    break
                rhs += abs(qj) / (d_lo - rho) ** 2
            if not ok or not abs(qi) / rho ** 2 > rhs:
                ok = False
                break
    check("RHO: domination inequality re-instantiated exactly for every charge", ok)

    # ---- BALL / OUTSIDE ----
    bad_ball = bad_out = 0
    for lf in leaves:
        if lf["kind"] == "excluded-ball":
            box = rebuild(root_box, lf["path"])
            i = lf["ball"]
            if not (isinstance(i, int) and 0 <= i < n):
                bad_ball += 1
                continue
            c = charges[i][1]
            s = sum(max(abs(box[k][0] - c[k]), abs(box[k][1] - c[k])) ** 2
                    for k in range(3))
            if s > rhos[i] ** 2:
                bad_ball += 1
        elif lf["kind"] == "excluded-outside":
            box = rebuild(root_box, lf["path"])
            s = Fraction(0)
            for k in range(3):
                if box[k][0] > 0:
                    s += box[k][0] ** 2
                elif box[k][1] < 0:
                    s += box[k][1] ** 2
            if s < r0_sq:
                bad_out += 1
    check("BALL: every excluded-ball box inside its domination ball (exact)",
          bad_ball == 0, f"{bad_ball} violations")
    check("BALL: every excluded-outside box outside the localization ball (exact)",
          bad_out == 0, f"{bad_out} violations")

    # ---- CENTROID (optional) ----
    if args.centroid:
        c = Fraction(1, 3)
        hits = sum(1 for enc in encs
                   if all(enc[k][0] <= c <= enc[k][1] for k in range(3)))
        check("CENTROID: (1/3,1/3,1/3) lies in exactly one enclosure", hits == 1,
              f"{hits} hits")

    print("FAILURES:", len(FAILURES))
    return 1 if FAILURES else 0


if __name__ == "__main__":
    sys.exit(main())
