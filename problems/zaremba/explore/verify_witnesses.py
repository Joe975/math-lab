#!/usr/bin/env python3
"""Integrity + spot-check of a zscan witness file (tier 1, explore).

Checks, against the independent stdlib implementation in
harness/zaremba/zaremba.py (not the C kernel that wrote the file):

1. the file has exactly one line per n = 1..N, in order;
2. every witness in a fixed-seed random sample of `samples` lines is a
   valid Zaremba certificate (coprime, all canonical quotients <= 5);
3. every explicitly listed n (e.g. the per-slice worst cases from the scan
   log) verifies too.

Usage:
    python verify_witnesses.py <witness_file> [samples] [n1 n2 ...]

Exit 0 iff everything checks. Deterministic: the sample uses seed 20260727.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "harness" / "zaremba"))
from zaremba import cf  # noqa: E402


def main() -> int:
    path = sys.argv[1]
    samples = int(sys.argv[2]) if len(sys.argv) > 2 else 2000
    extra = {int(a) for a in sys.argv[3:]}

    witness: dict[int, int] = {}
    expect = 1
    for line in open(path, encoding="utf-8"):
        n, m = map(int, line.split())
        if n != expect:
            print(f"FAIL sequence: expected n={expect}, got n={n}")
            return 1
        witness[n] = m
        expect += 1
    top = expect - 1
    print(f"sequence OK: one line per n = 1..{top}, in order")

    rng = random.Random(20260727)
    picks = sorted(rng.sample(range(2, top + 1), samples) if top > samples else
                   range(2, top + 1))
    picks = sorted(set(picks) | (extra & set(witness)))
    missing = extra - set(witness)
    if missing:
        print(f"FAIL: requested n not in file: {sorted(missing)}")
        return 1

    bad = 0
    for n in picks:
        quots = cf(n, witness[n])  # raises if not coprime / out of range
        if max(quots) > 5:
            print(f"FAIL n={n} m={witness[n]}: K={max(quots)} > 5")
            bad += 1
    print(f"verified {len(picks)} witnesses "
          f"({len(extra & set(witness))} explicitly requested): {bad} bad")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
