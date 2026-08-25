"""Two independent enumerators for Giuga numbers.

``scan_giuga_numbers``  -- brute force over the integers below a limit, using a
smallest-prime-factor sieve.  Obviously correct, hopeless past ~1e8.

``enumerate_by_prime_sets`` -- recursion over sets of primes using the rational
form of the Giuga condition:

    n = p_1 ... p_m squarefree,   p_i | (n/p_i - 1) for all i
      <=>  sum_i 1/p_i - 1/n = k  for some integer k >= 1.

(Proof of the equivalence: multiplying by n, the condition n | sum_i n/p_i - 1
splits by CRT into, for each j, p_j | (sum_i n/p_i - 1); every term with i != j
is divisible by p_j, leaving p_j | (n/p_j - 1).)

The two share no arithmetic: one tests integers, the other builds prime sets.

Pruning rules used by ``enumerate_by_prime_sets``, with their justification --
each is a *necessary* condition for a completion to exist, so nothing that
could extend to a solution is discarded:

  P1  k >= 1.  sum 1/p_i - 1/n is positive (m >= 2) and an integer.
  P2  k <= floor(sum of 1/q over the first m primes).  The sum over any m
      distinct primes is at most that.
  P3  With j primes chosen and sum s, the remaining r = m - j primes must
      contribute exactly k - s + 1/n > k - s, so k - s > 0.
  P4  The remaining r primes are all > last and distinct, so their reciprocal
      sum is at most that of the r smallest primes exceeding last; if that is
      <= k - s, no completion exists.
  P5  The next prime p is the smallest of the r remaining, so
      r/p >= (their sum) > k - s, i.e. p < r/(k - s).
  P6  At r == 1 the last prime is forced: 1/p - 1/(Np) = k - s = a/b gives
      p = b(N-1)/(aN); accept only if that is an integer, prime, and > last.
  P7  At r == 2 the larger of the two is forced by the smaller:
      1/p + 1/q - 1/(Npq) = a/b gives q = b(Np-1)/(N(ap-b)), so only p is
      searched, over the window b/a < p < 2b/a implied by q > p.
  P8  At r == 2, a closed form that replaces that scan with a factorisation.
      The remainder k - s is always a/N in lowest terms: writing it over N, the
      numerator is kN - sum_i N/p_i, which mod p_j is -N/p_j, and N/p_j is
      coprime to p_j because N is squarefree -- so the numerator is coprime to
      every prime factor of N.  With b = N the last-two equation
      Nq + Np - 1 = apq becomes apq - Np - Nq + 1 = 0, and multiplying by a
      completes the rectangle:

          (a*p - N) * (a*q - N) = N^2 - a.

      So p and q come from the divisor pairs of N^2 - a, and the search cost
      moves from scanning a window of width ~N/a to factoring one integer.
      Both signs are tried, since the two factors need only share a sign.
      P8 is complete exactly when N^2 - a factors; when it does not, the branch
      is appended to ``unresolved`` rather than dropped, because silently
      skipping it would turn a partial enumeration into a false "none exist".

Standard library only, exact rational arithmetic (fractions.Fraction).
"""

from fractions import Fraction

from .conditions import is_prime, is_giuga_number
from .factorint import factorize_full
from .lasttwo import completions


# --------------------------------------------------------------------------
# implementation 1: sieve over the integers
# --------------------------------------------------------------------------

def scan_giuga_numbers(limit: int):
    """All Giuga numbers n <= limit, by direct test of every n."""
    spf = list(range(limit + 1))
    i = 2
    while i * i <= limit:
        if spf[i] == i:
            for j in range(i * i, limit + 1, i):
                if spf[j] == j:
                    spf[j] = i
        i += 1

    out = []
    for n in range(2, limit + 1):
        m, ps, ok = n, [], True
        while m > 1:
            p = spf[m]
            m //= p
            if m % p == 0:
                ok = False
                break
            ps.append(p)
        if not ok or len(ps) < 2:
            continue
        if all((n // p - 1) % p == 0 for p in ps):
            out.append(n)
    return out


# --------------------------------------------------------------------------
# implementation 2: recursion over prime sets
# --------------------------------------------------------------------------

class SieveTooLarge(RuntimeError):
    """The branch needed candidate primes beyond any sievable range."""


class _Primes:
    """Growable sieved list of primes with a 'next prime after x' lookup."""

    def __init__(self):
        self._limit = 1 << 12
        self._ps = self._sieve(self._limit)

    @staticmethod
    def _sieve(n):
        s = bytearray([1]) * (n + 1)
        s[0:2] = b"\x00\x00"
        i = 2
        while i * i <= n:
            if s[i]:
                s[i * i::i] = bytearray(len(s[i * i::i]))
            i += 1
        return [i for i in range(n + 1) if s[i]]

    MAX_SIEVE = 1 << 27          # ~134M; beyond this the enumeration is not viable

    def _grow(self, upto):
        if upto <= self._limit:
            return
        if upto > self.MAX_SIEVE:
            raise SieveTooLarge(
                f"needed primes up to {upto}, cap is {self.MAX_SIEVE}")
        while self._limit < upto:
            self._limit *= 2
        self._ps = self._sieve(self._limit)

    def _index_above(self, x, ps):
        lo, hi = 0, len(ps)
        while lo < hi:                       # bisect without importing
            mid = (lo + hi) // 2
            if ps[mid] <= x:
                lo = mid + 1
            else:
                hi = mid
        return lo

    def upto(self, limit):
        """Every prime <= limit, as a fresh list."""
        self._grow(limit)
        return self._ps[:self._index_above(limit, self._ps)]

    def between(self, low, high):
        """Primes p with low < p <= high, without copying the list.

        ``upto`` copies its prefix on every call, which at a limit of 8e6 costs
        ~2.5 ms -- more than everything else a search node does. Growing the
        sieve first and then capturing the list means a deeper call that grows
        it again rebinds ``self._ps`` without disturbing this iteration, and the
        captured list already covers everything up to ``high``.
        """
        self._grow(high)
        ps = self._ps
        i = self._index_above(low, ps)
        n = self._index_above(high, ps)
        while i < n:
            yield ps[i]
            i += 1

    def next_after(self, x, count):
        """The ``count`` smallest primes strictly greater than x."""
        out = []
        n = x
        while len(out) < count:
            n += 1
            if is_prime(n):
                out.append(n)
        return out


PRIMES = _Primes()


def max_tail_sum(last: int, r: int) -> Fraction:
    """Largest possible sum of 1/p over r distinct primes all > last."""
    return sum((Fraction(1, p) for p in PRIMES.next_after(last, r)), Fraction(0))


def k_upper_bound(m: int, min_prime: int = 2) -> int:
    """P2: floor of the largest reciprocal sum available to m distinct primes >= min_prime."""
    ps = PRIMES.next_after(min_prime - 1, m)
    return int(sum((Fraction(1, p) for p in ps), Fraction(0)))


def enumerate_by_prime_sets(m: int, min_prime: int = 2, k=None, on_node=None,
                            use_p4: bool = True, use_p7: bool = True,
                            use_p8: bool = False, unresolved=None,
                            rho_budget: int = 400_000,
                            window_cap: int = 2_000_000):
    """All sets of m distinct primes >= min_prime whose product is a Giuga number.

    Yields sorted tuples of primes.  ``k`` restricts to a single value of the
    integer sum 1/p_i - 1/n; None means all admissible k (P1, P2).

    ``use_p4`` / ``use_p7`` switch off individual pruning rules.  Turning one
    off must not change the output -- only the cost.  That is the check that
    the rule was necessary rather than merely convenient.

    ``use_p8`` replaces the last-two-primes scan with the factorisation closed
    form, which is what makes larger factor counts reachable.  Pass a list as
    ``unresolved`` to collect the branches whose N^2 - a did not factor; if that
    list is non-empty the enumeration is complete only outside those branches,
    and saying so is the caller's job.
    """
    ks = [k] if k is not None else list(range(1, k_upper_bound(m, min_prime) + 1))
    for kk in ks:
        yield from _rec([], Fraction(0), 1, m, kk, min_prime - 1, on_node,
                        use_p4, use_p7, use_p8, unresolved, rho_budget, window_cap)


def _rec(chosen, s, N, m, k, last, on_node, use_p4=True, use_p7=True,
         use_p8=False, unresolved=None, rho_budget=400_000,
         window_cap=2_000_000):
    if on_node is not None:
        on_node()
    r = m - len(chosen)
    rem = k - s
    if rem <= 0:                                   # P3
        return
    if r == 1:                                     # P6: last prime is forced
        a, b = rem.numerator, rem.denominator
        num, den = b * (N - 1), a * N
        if num % den == 0:
            p = num // den
            if p > last and is_prime(p):
                yield tuple(chosen + [p])
        return
    if r == 2 and use_p8:                          # P8: last two, closed form
        a, b = rem.numerator, rem.denominator
        assert b == N, ("denominator lemma violated", chosen, rem)
        pairs, complete = completions(last, N, a, window_cap=window_cap,
                                      rho_budget=rho_budget, require_prime=True)
        if not complete:
            if unresolved is not None:
                unresolved.append({
                    "chosen": list(chosen), "N": N, "a": a, "M": N * N - a,
                    "cofactor": factorize_full(N * N - a, budget=rho_budget)[1]})
            return
        for p, q in pairs:
            yield tuple(chosen + [p, q])
        return
    if r == 2 and use_p7:                          # P7: last two, q forced by p
        a, b = rem.numerator, rem.denominator
        # 1/p + 1/q - 1/(Npq) = a/b  =>  q = b(Np - 1) / (N(ap - b)),
        # so p > b/a; and q > p forces 2/p > a/b, so p < 2b/a.
        lo = b // a
        hi2 = (2 * b) // a
        for p in PRIMES.between(max(last, lo), hi2):
            den = N * (a * p - b)
            num = b * (N * p - 1)
            if den <= 0 or num % den:
                continue
            q = num // den
            if q > p and is_prime(q):
                yield tuple(chosen + [p, q])
        return
    if use_p4 and max_tail_sum(last, r) <= rem:    # P4
        return
    hi = Fraction(r) / rem                         # P5: next prime < r/rem
    limit = hi.numerator // hi.denominator
    for p in PRIMES.between(last, limit):
        if Fraction(p) >= hi:
            break
        yield from _rec(chosen + [p], s + Fraction(1, p), N * p, m, k, p,
                        on_node, use_p4, use_p7, use_p8, unresolved, rho_budget,
                        window_cap)


def product(ps):
    n = 1
    for p in ps:
        n *= p
    return n


def verify_set(ps) -> bool:
    """Independent re-check of a returned set: the pointwise divisibility form."""
    n = product(ps)
    return is_giuga_number(n, {p: 1 for p in ps})
