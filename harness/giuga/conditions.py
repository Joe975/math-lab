"""Condition checkers for Giuga's conjecture.

Three independent faces of the same test, kept separate on purpose so that a
number produced by one can be re-derived by another:

  1. ``power_sum_mod``          -- literal sum  sum_{k=1}^{n-1} k^(n-1) mod n.
  2. ``power_sum_mod_structural`` -- the same residue rebuilt from the prime
     factorisation via CRT, using Fermat/"sum of powers mod p" per factor.
  3. ``agoh_holds``             -- n*B_{n-1} == -1 (mod n) with exact rational
     Bernoulli numbers.  Shares no arithmetic with 1 or 2.

Plus the structural conditions a composite solution has to satisfy:
``is_giuga_number``, ``is_carmichael``, ``is_counterexample``.

Exact integer / rational arithmetic only.  Standard library only.

Definitions used (all standard):

* A *Giuga number* is a composite n with  p | (n/p - 1)  for every prime p | n.
* A *Carmichael number* is a squarefree composite n with  (p-1) | (n-1)  for
  every prime p | n.
* Giuga's congruence is  sum_{k=1}^{n-1} k^(n-1) == -1 (mod n).
"""

from fractions import Fraction
from math import comb, gcd

# --------------------------------------------------------------------------
# factorisation / primality
# --------------------------------------------------------------------------

_SMALL_PRIME_BASES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)


def is_prime(n: int) -> bool:
    """Deterministic Miller-Rabin (correct for all n < 3.3e24)."""
    if n < 2:
        return False
    for p in _SMALL_PRIME_BASES:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _SMALL_PRIME_BASES:
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def factorize(n: int) -> dict:
    """Trial division.  Intended for n small enough that this is fine."""
    if n < 1:
        raise ValueError("factorize needs n >= 1")
    f = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            f[d] = f.get(d, 0) + 1
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def prime_factors(n: int) -> list:
    return sorted(factorize(n))


def is_squarefree(n: int) -> bool:
    return all(e == 1 for e in factorize(n).values())


# --------------------------------------------------------------------------
# structural conditions
# --------------------------------------------------------------------------

def is_giuga_number(n: int, factors=None) -> bool:
    """Composite n with p | (n/p - 1) for every prime p | n.

    Note this forces squarefreeness: if p^2 | n then p | n/p, so
    p | (n/p - 1) would give p | 1.
    """
    if n < 2:
        return False
    f = factorize(n) if factors is None else factors
    if len(f) == 1 and next(iter(f.values())) == 1:
        return False  # prime
    for p, e in f.items():
        if e != 1:
            return False
        if (n // p - 1) % p != 0:
            return False
    return True


def is_carmichael(n: int, factors=None) -> bool:
    """Korselt: squarefree composite n with (p-1) | (n-1) for every p | n."""
    if n < 2:
        return False
    f = factorize(n) if factors is None else factors
    if len(f) == 1 and next(iter(f.values())) == 1:
        return False  # prime
    if any(e != 1 for e in f.values()):
        return False
    return all((n - 1) % (p - 1) == 0 for p in f)


def is_counterexample(n: int) -> bool:
    """Composite n satisfying Giuga's congruence, via the structural criterion.

    A composite n satisfies the congruence iff it is simultaneously a Giuga
    number and a Carmichael number, i.e. iff for every prime p | n both
    p | (n/p - 1) and (p-1) | (n/p - 1).
    """
    f = factorize(n)
    return is_giuga_number(n, f) and is_carmichael(n, f)


def counterexample_conditions(primes) -> bool:
    """The pointwise form: for each p, p | (n/p - 1) and (p-1) | (n/p - 1).

    ``primes`` is an iterable of distinct primes; n is their product.  This is
    written directly from the pointwise statement (no is_giuga/is_carmichael
    call) so it can disagree with ``is_counterexample`` if either is wrong.
    """
    ps = list(primes)
    if len(ps) < 2 or len(set(ps)) != len(ps):
        return False
    n = 1
    for p in ps:
        n *= p
    for p in ps:
        q = n // p - 1
        if q % p != 0 or q % (p - 1) != 0:
            return False
    return True


def giuga_parity_ok(m: int, k: int, odd: bool) -> bool:
    """Necessary parity relation between factor count and k for odd n.

    If n is odd and squarefree with m prime factors, every n/p_i is odd, so
    sum_i n/p_i == m (mod 2).  The Giuga condition is sum_i n/p_i - 1 = k*n, and
    n odd gives k*n == k (mod 2).  Hence

        m == 1 + k   (mod 2),

    so for k = 1 an odd Giuga number has an EVEN number of prime factors.  The
    relation says nothing when n is even, since then n/p_i need not be odd.
    """
    if not odd:
        return True
    return (m - 1 - k) % 2 == 0


def giuga_sum_condition(primes) -> bool:
    """sum 1/p - 1/prod(p) is a positive integer (Giuga condition, rational form)."""
    ps = list(primes)
    n = 1
    for p in ps:
        n *= p
    v = sum((Fraction(1, p) for p in ps), Fraction(0)) - Fraction(1, n)
    return v.denominator == 1 and v > 0


# --------------------------------------------------------------------------
# face 1: the literal power sum
# --------------------------------------------------------------------------

def power_sum_mod(n: int) -> int:
    """sum_{k=1}^{n-1} k^(n-1) mod n, by direct summation.  O(n log n) bigint ops."""
    e = n - 1
    return sum(pow(k, e, n) for k in range(1, n)) % n


def satisfies_congruence(n: int) -> bool:
    return power_sum_mod(n) == (n - 1) % n


# --------------------------------------------------------------------------
# face 2: the same residue rebuilt from the factorisation
# --------------------------------------------------------------------------

def _crt(residues) -> int:
    """residues: list of (r, m) with pairwise coprime m."""
    r0, m0 = 0, 1
    for r, m in residues:
        g = gcd(m0, m)
        assert g == 1, "moduli must be coprime"
        # solve x == r0 (mod m0), x == r (mod m)
        inv = pow(m0, -1, m)
        t = ((r - r0) * inv) % m
        r0 = r0 + m0 * t
        m0 *= m
    return r0 % m0


def power_sum_mod_structural(n: int):
    """sum_{k=1}^{n-1} k^(n-1) mod n for squarefree n, via CRT.

    For each prime p | n the residues k mod p run over 0..p-1 exactly n/p times,
    and sum_{r=1}^{p-1} r^t == -1 (mod p) when (p-1) | t and 0 otherwise.  So

        S(n) == -(n/p)  (mod p)   if (p-1) | (n-1)
        S(n) == 0       (mod p)   otherwise.

    Returns None for non-squarefree n (this routine does not cover it).
    """
    f = factorize(n)
    if any(e != 1 for e in f.values()):
        return None
    residues = []
    for p in f:
        if (n - 1) % (p - 1) == 0:
            residues.append(((-(n // p)) % p, p))
        else:
            residues.append((0, p))
    return _crt(residues)


def satisfies_congruence_structural(n: int) -> bool:
    """Whether S(n) == -1 (mod n), decided structurally, all n >= 2.

    If p^2 | n then p | n/p, so S(n) == (n/p) * sum_{r} r^(n-1) == 0 (mod p),
    while -1 is not 0 mod p; hence no non-squarefree n satisfies the congruence.
    """
    if n < 2:
        return False
    f = factorize(n)
    if any(e != 1 for e in f.values()):
        return False
    for p in f:
        if (n - 1) % (p - 1) != 0:
            return False          # S(n) == 0 (mod p), but -1 is not
        if (n // p - 1) % p != 0:
            return False          # S(n) == -(n/p) (mod p), need n/p == 1
    return True


# --------------------------------------------------------------------------
# face 3: Agoh's Bernoulli-number form
# --------------------------------------------------------------------------

def bernoulli(m: int) -> list:
    """[B_0 .. B_m] as exact Fractions, convention B_1 = -1/2."""
    b = [Fraction(0)] * (m + 1)
    b[0] = Fraction(1)
    for k in range(1, m + 1):
        s = sum(comb(k + 1, j) * b[j] for j in range(k))
        b[k] = -s / (k + 1)
    return b


def agoh_holds(n: int, b=None) -> bool:
    """n * B_{n-1} == -1 (mod n), as a congruence of n-integral rationals.

    x == -1 (mod n) means x + 1 = a/b in lowest terms with gcd(b, n) == 1 and
    n | a.  Bernoulli numbers with odd index >= 3 vanish, so every even n >= 4
    fails here immediately -- an independent route to "no even solution".
    """
    if n < 2:
        return False
    if b is None:
        b = bernoulli(n - 1)
    x = Fraction(n) * b[n - 1] + 1
    if gcd(x.denominator, n) != 1:
        return False
    return x.numerator % n == 0
