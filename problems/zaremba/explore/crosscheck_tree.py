#!/usr/bin/env python3
"""Independent z(n) computation by forward generation (tier 1, explore).

Adversarial cross-check for the Euclid-based kernel in harness/zaremba/.
This computes z(n) by the OPPOSITE route: instead of expanding every m/n and
taking a minimum, it generates every canonical bounded-quotient continued
fraction [0; a1..ak] (ai <= 5, ak >= 2) with denominator <= N via the
continuant recurrence q_i = a_i*q_{i-1} + q_{i-2}, and records for each
denominator the smallest maximum digit along any generating word. Continuants
of a word are automatically coprime, so no gcd is ever computed; no Euclidean
division happens anywhere. Agreement with the kernel is therefore evidence
about correctness, not about determinism.

Usage:
    python crosscheck_tree.py N [table_file]

Prints the z histogram for 2..N and every n with z(n) >= 4; verifies
coverage (every n in 2..N reached with some digit <= 5). With table_file,
writes "n z" per line for diffing against `zscan exact 2 N <file>`.

Standard library only; exact integer arithmetic.
"""

from __future__ import annotations

import sys
from collections import Counter


def z_table(limit: int) -> bytearray:
    """best[n] = min over canonical words with continuant n of the max digit.

    0 means unreached. Iterative DFS over (q_prev, q_cur, max_digit_so_far);
    a state is a word a1..ak with k >= 1; it is a canonical word iff its last
    digit is >= 2, in which case q_cur is a denominator it witnesses.
    """
    best = bytearray(limit + 1)
    # State: (q_prev, q_cur, max digit so far, last digit). Roots are the
    # length-1 words [0; a] = 1/a with q_0 = 1, q_1 = a. A word is canonical
    # iff its LAST digit is >= 2 -- tracked explicitly, because clever
    # continuant-ratio tests break on the word [1, 1].
    stack = [(1, a, a, a) for a in range(1, 6) if a <= limit]
    while stack:
        q_prev, q_cur, worst, last = stack.pop()
        if last >= 2 and (best[q_cur] == 0 or worst < best[q_cur]):
            best[q_cur] = worst
        for a in range(1, 6):
            q_next = a * q_cur + q_prev
            if q_next > limit:
                break
            stack.append((q_cur, q_next, max(worst, a), a))
    best[1] = 1  # z(1) = 1 by convention
    return best


def main() -> int:
    limit = int(sys.argv[1])
    table_path = sys.argv[2] if len(sys.argv) > 2 else None

    best = z_table(limit)

    unreached = [n for n in range(2, limit + 1) if best[n] == 0]
    for n in unreached:
        print(f"FAIL {n}: no witness with digits <= 5")

    hist = Counter(best[n] for n in range(2, limit + 1) if best[n])
    for n in range(2, limit + 1):
        if best[n] >= 4:
            print(f"z {n} = {best[n]}")
    print(f"histogram [2,{limit}]:")
    for v in sorted(hist):
        print(f"  z={v}: {hist[v]}")

    if table_path:
        with open(table_path, "w", encoding="utf-8") as fh:
            for n in range(2, limit + 1):
                fh.write(f"{n} {best[n]}\n")

    return 1 if unreached else 0


if __name__ == "__main__":
    raise SystemExit(main())
