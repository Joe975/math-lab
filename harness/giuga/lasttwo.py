"""The last two entries of a Giuga sequence, by whichever route is cheaper.

Given a prefix with product N and remainder a/N, the final pair (x, y) satisfies

    a*x*y - N*x - N*y + 1 = 0,

which can be solved two ways:

  WINDOW   y = (N*x - 1) / (a*x - N), and q > x forces N/a < x <= 2N/a.  Walking
           that window settles the prefix with no factorisation at all.  Cost is
           the window width N/a.

  RECT     (a*x - N)(a*y - N) = N^2 - a, so the pair comes from a divisor pair
           of N^2 - a.  Cost is one factorisation, independent of the window.

Neither dominates: a prefix whose reciprocal sum is far from 1 has a narrow
window and is best walked, while one very close to 1 has a huge window but its
N^2 - a may factor instantly.  ``completions`` picks per prefix and reports
whether the answer is complete, so a caller can never mistake "gave up" for
"none exist".
"""

from .conditions import is_prime
from .factorint import divisors, factorize_full


def completions(last, N, a, window_cap=2_000_000, rho_budget=400_000,
                odd_only=False, require_prime=True):
    """(pairs, complete).

    ``pairs`` are the (x, y) with last < x < y completing the prefix.
    ``complete`` is False only when both routes were refused, in which case
    ``pairs`` is empty and the branch is undecided rather than empty.
    """
    lo = max(N // a, last)
    hi = (2 * N) // a + 1
    if hi - lo <= window_cap:
        out = []
        for x in range(lo + 1, hi + 1):
            d = a * x - N
            if d <= 0:
                continue
            num = N * x - 1
            if num % d:
                continue
            y = num // d
            if y <= x:
                continue
            if odd_only and (x % 2 == 0 or y % 2 == 0):
                continue
            if require_prime and not (is_prime(x) and is_prime(y)):
                continue
            out.append((x, y))
        return out, True

    M = N * N - a
    fac, cofactor = factorize_full(M, budget=rho_budget)
    if cofactor != 1:
        return [], False
    out, seen = [], set()
    for d in divisors(fac):
        for sd in (d, -d):
            num = sd + N
            if num <= 0 or num % a:
                continue
            x = num // a
            if x <= last:
                continue
            num2 = (M // d if sd > 0 else -(M // d)) + N
            if num2 <= 0 or num2 % a:
                continue
            y = num2 // a
            if y <= x or (x, y) in seen:
                continue
            if odd_only and (x % 2 == 0 or y % 2 == 0):
                continue
            if require_prime and not (is_prime(x) and is_prime(y)):
                continue
            seen.add((x, y))
            out.append((x, y))
    return out, True
