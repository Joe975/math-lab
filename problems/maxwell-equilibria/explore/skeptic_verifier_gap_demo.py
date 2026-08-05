#!/usr/bin/env python3
"""Demonstration: verify_equilibria.py's tiling check is not sound by itself.

The independent checker validates the leaf decomposition of a certificate by
(i) prefix-freeness of the split paths and (ii) the Kraft equality
sum 2^-depth = 1. For paths that carry a PER-NODE AXIS CHOICE those two
conditions do not imply a tiling: {"01", "11"} (axis-0 upper half and axis-1
upper half of the root) is prefix-free with Kraft sum 1, yet the two boxes
overlap and jointly miss the (lo, lo) quadrant of the root.

This script shows the consequence at both severities:

  DEMO 1  A synthetic certificate whose two leaves overlap and leave a hole
          (no equilibrium in the region, every leaf genuinely sign-excluded):
          verify() returns PASS on a non-tiling leaf list.

  DEMO 2  The real driver's certificate for two unit charges at (+-1, 0, 0)
          (exactly one equilibrium, at the origin) is tampered with: the
          isolated leaf is replaced by two mixed-axis exclusion leaves that
          each avoid the origin but jointly leave a hole containing it, and
          the count is edited from 1 to 0. If verify() PASSes, the checker
          has endorsed a "complete census" that misses a real equilibrium.

  Both demos are also run through the skeptic's own tree reconstruction
  (skeptic_cert_checks.check_tree), which rejects them at the mixed-axis
  node -- i.e. axis consistency is exactly the missing condition.

The driver itself (equilibria.py) always emits both children of a split with
the same axis digit, so certificates it actually produces are axis-consistent
and DO tile; the gap is in what the independent checker is able to catch,
not a wrong count in any existing certificate. See attempt 002.

Standard library only. Deterministic. Runtime ~1 minute (DEMO 2 runs the
two-charge census at 64 bits).
"""

from __future__ import annotations

import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "harness", "maxwell-equilibria"))
sys.path.insert(0, HERE)

import equilibria as eq                      # noqa: E402  (the driver)
import verify_equilibria as ve               # noqa: E402  (the checker under test)
import skeptic_cert_checks as scc            # noqa: E402  (this review's tree check)


def own_tree_check(paths) -> bool:
    """Run the skeptic tree reconstruction; True iff it accepts."""
    scc.FAILURES.clear()
    scc.check_tree(list(paths))
    ok = not scc.FAILURES
    scc.FAILURES.clear()
    return ok


def demo1() -> bool:
    print("=" * 72)
    print("DEMO 1: overlapping leaves + hole, no zero involved")
    print("=" * 72)
    cfg2 = {"charges": [{"q": "1", "pos": ["1", "0", "0"]},
                        {"q": "1", "pos": ["-1", "0", "0"]}]}
    # Region from the checker's own validation suite: F is sign-definite on
    # all of it, so ANY sub-box is genuinely sign-excludable.
    bad = {
        "config": cfg2,
        "region": [["2", "9/4"], ["0", "1/4"], ["0", "1/4"]],
        "complete": False,
        "domination_rhos": ["1/4", "1/4"],
        "count": 0, "unresolved": 0, "isolated": [],
        "leaves": [{"path": "01", "kind": "excluded-sign", "ball": None},
                   {"path": "11", "kind": "excluded-sign", "ball": None}],
    }
    # the two boxes: axis-0 upper half and axis-1 upper half of the region --
    # they overlap on (x hi, y hi) and miss (x lo, y lo) entirely.
    res = ve.verify(bad)
    print("verify_equilibria.verify() verdict:", res["verdict"], res["failures"])
    mine = own_tree_check(["01", "11"])
    print("skeptic tree reconstruction accepts:", mine)
    exploited = res["verdict"] == "PASS" and not mine
    print("GAP DEMONSTRATED:" if exploited else "gap NOT demonstrated:",
          "checker PASSes a non-tiling leaf list that the axis-consistency "
          "check rejects" if exploited else res)
    return exploited


def demo2() -> bool:
    print("=" * 72)
    print("DEMO 2: hiding the real origin equilibrium of two unit charges")
    print("=" * 72)
    eq.set_precision(64)
    cfg = eq.Config([(Fraction(1), (Fraction(1), Fraction(0), Fraction(0))),
                     (Fraction(1), (Fraction(-1), Fraction(0), Fraction(0)))])
    cert = eq.count_equilibria(cfg)
    print(f"honest driver run: count={cert['count']} "
          f"unresolved={cert['unresolved']} leaves={len(cert['leaves'])}")
    assert cert["count"] == 1 and cert["unresolved"] == 0
    iso = cert["isolated"][0]
    path = iso["path"]
    box = [(Fraction(lo), Fraction(hi)) for lo, hi in iso["box"]]
    # the unique equilibrium is the origin (by symmetry); it is interior here
    assert all(lo < 0 < hi for lo, hi in box)

    # Replace the isolated leaf with two mixed-axis half-leaves that each
    # avoid the origin; their union misses the quadrant of the origin.
    # Use the two axes whose midpoints are farthest from 0 (best exclusion
    # margin), and on each keep the half NOT containing the origin.
    mids = [(box[a][0] + box[a][1]) / 2 for a in range(3)]
    assert all(m != 0 for m in mids), "face through the origin"
    ax1, ax2 = sorted(range(3), key=lambda a: abs(mids[a]), reverse=True)[:2]
    fakes = []
    for axis in sorted((ax1, ax2)):
        mid = mids[axis]
        side = "1" if mid > 0 else "0"   # origin is in the other half
        fakes.append(f"{path}{axis}{side}")
    print("isolated leaf", repr(path), "replaced by mixed-axis leaves",
          [repr(f) for f in fakes])

    tampered = dict(cert)
    tampered["count"] = 0
    tampered["isolated"] = []
    tampered["index_sum"] = None
    tampered["index_sum_ok"] = None
    tampered["leaves"] = (
        [lf for lf in cert["leaves"] if lf["path"] != path]
        + [{"path": f, "kind": "excluded-newton", "ball": None} for f in fakes])
    tampered = json.loads(json.dumps(tampered))   # same round-trip as a file

    res = ve.verify(tampered)
    print("verify_equilibria.verify() verdict:", res["verdict"],
          "count:", res["count"], res["failures"][:3])
    mine = own_tree_check([lf["path"] for lf in tampered["leaves"]])
    print("skeptic tree reconstruction accepts:", mine)
    exploited = res["verdict"] == "PASS" and res["count"] == 0 and not mine
    print("GAP DEMONSTRATED:" if exploited else "gap NOT demonstrated:",
          "checker endorses a census claiming 0 equilibria for a region "
          "that provably contains one (the origin)" if exploited else res)
    return exploited


def main() -> int:
    d1 = demo1()
    d2 = demo2()
    print("=" * 72)
    print(f"DEMO 1 (overlap+hole PASSes): {d1};  "
          f"DEMO 2 (hidden real zero PASSes): {d2}")
    return 0 if d1 else 1   # demo 1 is the load-bearing demonstration


if __name__ == "__main__":
    sys.exit(main())
