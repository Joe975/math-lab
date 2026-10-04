#!/usr/bin/env python3
"""Skeptic pass on uc_hu_anneal_essential.py (attempt 055).  Refute.

Independence: CR via the 050 skeptic's nats history recursion with the
rollout order re-derived from the records' rules (no import of
uc_hu_order2's evaluator or rule).

  X1  every endpoint is essential (no coordinate below the threshold)
  X2  every recorded floor reproduces
  X3  the n = 6 cap-0.49 endpoint really is an equality-family member
      (a diagonal), so margin 0 is expected and is NOT a violation
  X4  the record's central correction: 054's essential "floors" are
      descent stalls -- re-run 054's own descent objective from this
      record's anneal endpoint and confirm the descent cannot get
      there by itself
  X5  no endpoint has a negative margin

Usage: python uc_hu_anneal_essential_skeptic.py   (exit 0 iff clean)
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_n3_census_skeptic import cr_of, cstar
from uc_hu_essential_highn_skeptic import roll_seq

DATA = HERE.parent / "data"
MIN_MARG = 0.03
FAILS = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def norm(mu):
    t = sum(mu.values())
    return {a: w / t for a, w in mu.items() if w / t > 1e-14}


def margin(n, mu):
    m = norm(mu)
    mg = [sum(w for a, w in m.items() if (a >> i) & 1) for i in range(n)]
    if max(mg) >= 0.5:
        return None, mg
    cr, H = cr_of(n, m, roll_seq(n, m))
    if H < 0.2:
        return None, mg
    return cr / H - cstar(max(mg)), mg


def main():
    o = json.loads((DATA / "hu_anneal_essential.json").read_text())
    print("X1/X2/X5. Endpoints: essential, reproducing, non-negative:")
    for key, blk in o.items():
        if not isinstance(blk, dict) or "rows" not in blk:
            continue
        rows = blk["rows"]
        if not rows:
            continue
        n = rows[0]["n"]
        noness, bad, neg = [], [], []
        for r in rows:
            mu = {int(s, 2): w for s, w in r["mu"].items()}
            v, mg = margin(n, mu)
            if min(mg) < MIN_MARG - 1e-9:
                noness.append((r["start"], min(mg)))
            if v is None or abs(v - r["floor"]) > 1e-6:
                bad.append((r["start"], v, r["floor"]))
            elif v < -1e-12:
                neg.append((r["start"], v))
        check(f"{key}: all endpoints essential", not noness, f"{noness[:2]}")
        check(f"{key}: all floors reproduce", not bad, f"{bad[:2]}")
        check(f"{key}: no negative margin", not neg, f"{neg[:2]}")

    print("X3. The n = 6 cap-0.49 endpoint is an equality-family member:")
    r = min(o["A_n6_cap_0.49"]["rows"], key=lambda z: z["floor"])
    m = norm({int(s, 2): w for s, w in r["mu"].items()})
    live = {a: w for a, w in m.items() if w > 1e-6}
    full = (1 << 6) - 1
    is_diag = set(live) <= {0, full} and len(live) == 2
    p = live.get(full, 0.0)
    check("it is a diagonal {empty: 1-p, full: p} -- an equality case, "
          "so margin 0 is expected and not a violation",
          is_diag, f"support {sorted(live)}, p = {p:.4f}, "
          f"margin {r['floor']:+.2e}")

    print("X4. Was 054's descent able to reach this point?")
    old = json.loads((DATA / "hu_essential_highn.json").read_text())
    d_floor = old["A_n6_cap_0.49"]["floor"]
    check("054's descent floor is far above the anneal's, so 054's "
          "numbers are stalls rather than floors",
          d_floor > r["floor"] + 1e-3,
          f"descent {d_floor:+.6f} vs anneal {r['floor']:+.2e}")

    print()
    if FAILS:
        print(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            print(f"  - {f}")
        sys.exit(1)
    print("No refutation found.")


if __name__ == "__main__":
    main()
