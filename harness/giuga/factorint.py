"""Deterministic integer factorisation, with an explicit give-up.

Used where a *complete* divisor list is needed, so the interesting part of the
interface is the failure mode: ``factorize_full`` returns the cofactor it could
not split rather than pretending, and callers are expected to record that
branch as unresolved instead of silently dropping it.

Deterministic: Pollard-Brent is seeded with c = 1, 2, 3, ... and y = 2, never
with a random start, so the same input always costs the same and the same
branches fail on every run.

Standard library only.
"""

from math import gcd

from .conditions import is_prime

_SMALL = None


def _small_primes(limit=1_000_000):
    global _SMALL
    if _SMALL is None:
        s = bytearray([1]) * (limit + 1)
        s[0:2] = b"\x00\x00"
        i = 2
        while i * i <= limit:
            if s[i]:
                s[i * i::i] = bytearray(len(s[i * i::i]))
            i += 1
        _SMALL = [i for i in range(limit + 1) if s[i]]
    return _SMALL


def _brent(n, c, budget):
    """A nontrivial factor of composite n, or None if the budget runs out."""
    if n % 2 == 0:
        return 2
    y, m, steps = 2, 128, 0
    g = r = q = 1
    x = ys = y
    while g == 1:
        x = y
        for _ in range(r):
            y = (y * y + c) % n
        k = 0
        while k < r and g == 1:
            ys = y
            for _ in range(min(m, r - k)):
                y = (y * y + c) % n
                q = q * abs(x - y) % n
            steps += min(m, r - k)
            if steps > budget:
                return None
            g = gcd(q, n)
            k += m
        r *= 2
    if g == n:
        g = 1
        while g == 1:
            ys = (ys * ys + c) % n
            g = gcd(abs(x - ys), n)
            steps += 1
            if steps > budget:
                return None
    return g if 1 < g < n else None


def factorize_full(n, small_limit=1_000_000, budget=400_000):
    """Return (factors, cofactor).

    ``factors`` maps prime -> exponent for the part that was resolved;
    ``cofactor`` is 1 on success, otherwise the composite that was not split.
    A caller that needs every divisor of n must treat cofactor != 1 as
    "this branch is unresolved", not as "n has no other factors".
    """
    factors = {}
    if n <= 1:
        return factors, 1
    for p in _small_primes(small_limit):
        if p * p > n:
            break
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
    if n == 1:
        return factors, 1
    if is_prime(n):
        factors[n] = factors.get(n, 0) + 1
        return factors, 1

    stack = [n]
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            factors[m] = factors.get(m, 0) + 1
            continue
        d = None
        for c in range(1, 12):
            d = _brent(m, c, budget)
            if d:
                break
        if not d:
            return factors, m
        stack.append(d)
        stack.append(m // d)
    return factors, 1


def divisors(factors):
    """Every divisor of a fully factored number, unsorted."""
    out = [1]
    for p, e in factors.items():
        out = [d * p ** i for d in out for i in range(e + 1)]
    return out
