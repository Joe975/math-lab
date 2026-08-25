"""Giuga sequences: the same condition over integers instead of primes.

A *Giuga sequence* is a list of distinct integers 2 <= n_1 < ... < n_m with

    sum_i 1/n_i - 1/(n_1...n_m)  =  k,  an integer.

A Giuga sequence whose entries are all prime is exactly the factorisation of a
Giuga number, so this is a strictly wider search: it can reach candidates the
prime search cannot see, and any all-prime member of it is a Giuga number.
That matters for the odd case, where the object is open -- odd entries may be
composite (3, 9, 15, 21, ...) so the reciprocal sum reaches 1 with fewer terms
than odd *primes* allow.

Two facts make the last two entries cheap, and neither uses primality:

  B = N.  Writing rem = k - sum_{i<=j} 1/n_i over the common denominator
      N = n_1...n_j, its reduced denominator b divides N.  Now suppose the
      branch completes with any r >= 1 further entries x_1..x_r, product P:

          sum_i 1/x_i - 1/(N*P) = a/b
          =>  b*N*sum_i (P/x_i) - b = a*N*P
          =>  N * (a*P - b*sum_i (P/x_i)) = -b,   so N | b.

      With b | N that forces b = N at EVERY node of a surviving branch, not
      just near the leaves.  And b = N holds exactly when the chosen entries are
      pairwise coprime: if a prime l divides two of them, every N/n_i still
      carries l so the numerator k*N - sum_i N/n_i is divisible by l; if l
      divides only n_j, the numerator is -N/n_j mod l, which is coprime to l.

      So: **all but the last entry of a Giuga sequence are pairwise coprime.**
      The last entry is unconstrained -- (2,3,7,43,1805) is a Giuga sequence and
      1805 = 5 * 19^2 is not even squarefree. This is the rule that makes the
      integer search tractable: candidates sharing a factor with N are dead.

  RECT.  With b = N that equation is axy - Nx - Ny + 1 = 0, and multiplying by
      a completes the rectangle:

          (a*x - N) * (a*y - N)  =  N^2 - a.

      So the last two entries come from the divisor pairs of N^2 - a, and the
      last one alone (r == 1) is forced to x = (N-1)/a.

Every other rule is a necessary condition on the remaining reciprocal budget,
so nothing that could extend to a solution is discarded:

  S1  k >= 1.                     S2  k <= floor(largest available sum).
  S3  rem = k - s > 0.            S4  rem < the largest sum r more entries can
                                      contribute (the r smallest allowed
                                      integers above ``last``).
  S5  the next entry x satisfies r/x >= (remaining sum) > rem, so x < r/rem.

Completeness depends on N^2 - a factoring; branches where it does not are
reported, never dropped.

Standard library only, exact rational arithmetic.
"""

from fractions import Fraction
from math import gcd

from .conditions import is_prime
from .factorint import factorize_full
from .lasttwo import completions


def _allowed(x, odd_only):
    return (x % 2 == 1) if odd_only else True


def _next_allowed(last, count, odd_only):
    """The ``count`` smallest allowed integers strictly greater than ``last``."""
    step = 2 if odd_only else 1
    start = last + 1
    if odd_only and start % 2 == 0:
        start += 1
    return [start + i * step for i in range(count)]


def max_tail_sum(last, r, odd_only):
    return sum((Fraction(1, x) for x in _next_allowed(last, r, odd_only)), Fraction(0))


def k_upper_bound(m, min_entry, odd_only):
    return int(sum((Fraction(1, x) for x in _next_allowed(min_entry - 1, m, odd_only)),
                   Fraction(0)))


def enumerate_sequences(m, min_entry=2, odd_only=False, k=None, on_node=None,
                        unresolved=None, rho_budget=400_000,
                        window_cap=2_000_000):
    """All Giuga sequences of length m with entries >= min_entry.

    Yields sorted tuples.  ``unresolved`` collects branches whose N^2 - a would
    not factor; a non-empty list means the enumeration is complete only outside
    those branches.
    """
    ks = ([k] if k is not None
          else list(range(1, k_upper_bound(m, min_entry, odd_only) + 1)))
    for kk in ks:
        yield from _rec([], Fraction(0), 1, m, kk, min_entry - 1, odd_only,
                        on_node, unresolved, rho_budget, window_cap)


def _rec(chosen, s, N, m, k, last, odd_only, on_node, unresolved, rho_budget,
         window_cap=2_000_000):
    if on_node is not None:
        on_node()
    r = m - len(chosen)
    rem = k - s
    if rem <= 0:                                        # S3
        return
    a, b = rem.numerator, rem.denominator

    if b != N:                                          # B = N, at every level
        return
    if r == 1:
        if (N - 1) % a == 0:
            x = (N - 1) // a
            if x > last and _allowed(x, odd_only):
                yield tuple(chosen + [x])
        return
    if r == 2:
        pairs, complete = completions(last, N, a, window_cap=window_cap,
                                      rho_budget=rho_budget, odd_only=odd_only,
                                      require_prime=False)
        if not complete:
            if unresolved is not None:
                unresolved.append({"chosen": list(chosen), "N": N, "a": a,
                                   "M": N * N - a,
                                   "cofactor": factorize_full(N * N - a,
                                                              budget=rho_budget)[1]})
            return
        for x, y in pairs:
            # x joins the prefix, so it must be coprime to N; y is the last
            # entry and is unconstrained.
            if gcd(x, N) == 1:
                yield tuple(chosen + [x, y])
        return

    if max_tail_sum(last, r, odd_only) <= rem:          # S4
        return
    hi = Fraction(r) / rem                              # S5
    limit = hi.numerator // hi.denominator
    step = 2 if odd_only else 1
    start = last + 1
    if odd_only and start % 2 == 0:
        start += 1
    for x in range(start, limit + 1, step):
        if Fraction(x) >= hi:
            break
        if gcd(x, N) != 1:      # every entry but the last joins the prefix
            continue
        yield from _rec(chosen + [x], s + Fraction(1, x), N * x, m, k, x,
                        odd_only, on_node, unresolved, rho_budget, window_cap)


def all_prime(seq):
    return all(is_prime(x) for x in seq)


def verify_sequence(seq):
    """Independent re-check from the definition, with exact rationals."""
    seq = list(seq)
    if len(seq) != len(set(seq)) or seq != sorted(seq) or any(x < 2 for x in seq):
        return False
    N = 1
    for x in seq:
        N *= x
    v = sum((Fraction(1, x) for x in seq), Fraction(0)) - Fraction(1, N)
    return v.denominator == 1 and v >= 1
