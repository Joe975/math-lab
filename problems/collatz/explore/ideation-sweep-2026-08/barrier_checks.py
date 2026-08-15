#!/usr/bin/env python3
"""Verify the two barrier claims from the 2026-08 ideation sweep (001).

Both are classical facts rediscovered by sweep workers; this re-derives them
exactly in range rather than trusting worker prose.

Check 1 — the -1 tower kills finite-memory pointwise descent certificates.
    For n = 2^L - 1, the first L-1 steps of the accelerated odd map
    A(n) = (3n+1)/2^{v2(3n+1)} all have valuation 1, so the block multiplier
    is exactly (3/2)^{L-1} and n stays ≡ -1 mod shrinking powers of 2.
    Consequence: for EVERY modulus 2^M the residue automaton has the state
    -1 mod 2^M as a self-loop with log-multiplier log(3/2) > 0, so any
    difference-constraint / LP system demanding pointwise descent with
    bounded lookahead over residue states is infeasible at every M.
    (Classical root: unbounded stopping times on 2^L - 1, Terras 1976.)

Check 2 — parity blocks carry zero entropy signal.
    For the SHORTENED map T(n) = n/2 (even), (3n+1)/2 (odd), the map
    (n mod 2^k) -> first k parity bits is a bijection (Terras 1976), so the
    parity block of a uniform residue is EXACTLY uniform on {0,1}^k:
    H = k bits, KL to fair coin = 0. Entropy/KL functionals of parity
    statistics cannot see the conjecture at any finite depth.
    NB: sweep worker 016 asserted this bijection for the UNSHORTENED map
    (odd -> 3n+1), where it is FALSE — residues 1 and 3 share parity word
    "10" at k = 2, because an odd step forces an even successor. The
    barrier statement only holds in the shortened form checked here.

Usage: python3 barrier_checks.py   (exact integer arithmetic, ~seconds)
"""

L_MAX = 200
K_MAX = 16


def v2(n: int) -> int:
    return (n & -n).bit_length() - 1


def accelerated_valuations(n: int, steps: int) -> list[int]:
    out = []
    for _ in range(steps):
        m = 3 * n + 1
        a = v2(m)
        out.append(a)
        n = m >> a
    return out


def check_tower() -> None:
    for L in range(2, L_MAX + 1):
        n = (1 << L) - 1
        vals = accelerated_valuations(n, L - 1)
        assert vals == [1] * (L - 1), (L, vals[:10])
        # and the orbit stays in the -1 residue class at every scale:
        m = n
        for k in range(L - 1, 0, -1):
            assert m % (1 << k) == (1 << k) - 1, (L, k)
            m = (3 * m + 1) >> 1
    print(f"check 1 OK: n = 2^L - 1 has L-1 valuation-1 steps for all L <= {L_MAX}")
    print("  => every residue automaton mod 2^M contains the -1 self-loop "
          "with multiplier 3/2 > 1")


def parity_word(n: int, k: int) -> int:
    """First k parity bits under the shortened map T."""
    w = 0
    for i in range(k):
        b = n & 1
        w |= b << i
        n = (3 * n + 1) >> 1 if b else n >> 1
    return w


def check_bijection() -> None:
    for k in range(1, K_MAX + 1):
        words = {parity_word(r, k) for r in range(1 << k)}
        assert len(words) == 1 << k, k
    print(f"check 2 OK: residue -> parity-word map is a bijection mod 2^k "
          f"for all k <= {K_MAX}")
    print("  => parity blocks of a uniform residue are exactly uniform: "
          "H = k bits, KL = 0")


if __name__ == "__main__":
    check_tower()
    check_bijection()
