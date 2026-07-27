"""Lower bounds on the number of prime factors of a composite solution.

A composite n satisfying Giuga's congruence has a prime set S = {p : p | n} that
is *admissible* in the following sense (each clause is proved in
``docs``-free form here; see the module tests for the mechanised checks):

  A1  n is squarefree, so S determines n.
  A2  every p in S is odd.
  A3  for any two distinct p, q in S,  q does not divide p - 1.
  A4  sum_{p in S} 1/p > 1.

Call a finite set of primes satisfying A2-A4 an *admissible set*.  Then

    (number of prime factors of a solution) >= min { |S| : S admissible }
    (a solution n)                          >= min { prod S : S admissible }

and both minima can be bounded from below by exhaustive search, because A3
thins the supply of usable primes so aggressively that A4 needs a great many of
them.

This module computes certified lower bounds for those two minima.  It never
uses floating point: reciprocals are represented as integers scaled by
``SCALE`` with *directed* rounding, so every comparison it makes is a
conservative one.

Search shape
------------
Fix a branch limit X.  Every admissible set S splits as S = T u B with
T = S n [3, X] and B = S n (X, inf).  The search enumerates every candidate T
(as an increasing sequence of primes <= X satisfying A3) and, for each, bounds
the contribution of B without enumerating it: B's members are all > X and all
T-admissible, so their reciprocal sum is at most that of the |B| smallest
T-admissible primes above X.  That gives the smallest |B| that could possibly
close the gap to 1, hence a lower bound |T| + |B|_min.  The minimum of that
over all T is a lower bound on |S|, and it is valid for *every* admissible set,
including sets whose members above X are not the ones the bound used.

Pruning (all of it sound, i.e. it only ever discards branches whose every
completion already has a bound at least as large as one already recorded):

  R1  A prime q may be appended to T only if no p in T divides q - 1 (A3
      itself, not an optimisation).
  R2  At a node (T, last) the whole subtree has bound >= |T| + t where t is the
      least number of T-admissible primes above ``last`` whose reciprocals sum
      to more than 1 - sum(1/p for p in T).  Every completion of T draws its
      remaining primes from exactly that pool, so no completion beats this.
      If |T| + t >= best-so-far, the subtree cannot lower the minimum: prune.
  R3  In the loop over the next prime q, the *relaxed* child bound (filtering
      only by T, not by T u {q}) is non-increasing in the achievable sum and so
      non-decreasing in q; once it reaches best-so-far, every later q is also
      hopeless: break.

R2 and R3 lower-bound the subtree, and the quantity being minimised is itself a
lower bound, so pruning cannot make the reported bound larger than the truth --
only smaller (i.e. weaker).  That is the safe direction, and it is checked
directly by ``search(..., prune=False)``.
"""

from fractions import Fraction

SCALE = 10 ** 24


def primes_upto(n: int):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    i = 2
    while i * i <= n:
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
        i += 1
    return [i for i in range(n + 1) if sieve[i]]


class PoolExhausted(RuntimeError):
    """The prime pool was too small to certify a bound.  Enlarge it."""


class Pool:
    """Odd primes below ``limit``, with reciprocals scaled to integers.

    ``hi[i]`` >= SCALE/p and ``lo[i]`` <= SCALE/p, so sums of ``hi`` over-state
    and sums of ``lo`` under-state the true reciprocal sums.

    ``mask[i]`` has bit j set iff ``branch[j]`` divides ``p_i - 1``, where
    ``branch`` is the list of odd primes <= X.  A set T of branchable primes is
    a bitmask ``tm``; prime i is T-admissible iff ``mask[i] & tm == 0``.
    """

    def __init__(self, limit: int, branch_limit: int):
        self.limit = limit
        self.branch_limit = branch_limit
        self.p = [q for q in primes_upto(limit) if q > 2]
        self.hi = [-(-SCALE // q) for q in self.p]
        self.lo = [SCALE // q for q in self.p]
        self.branch = [q for q in self.p if q <= branch_limit]
        self.nbranch = len(self.branch)
        self.first_above_X = self.nbranch  # index of first pool prime > X
        self.cap_hits = 0
        mask = []
        for q in self.p:
            m = 0
            qm1 = q - 1
            for j, b in enumerate(self.branch):
                if qm1 % b == 0:
                    m |= 1 << j
            mask.append(m)
        self.mask = mask

    def t_min(self, tm: int, s_hi: int, start: int, budget: int):
        """Least t <= budget with s_hi + (sum of hi over the t smallest
        T-admissible pool primes at index >= start) > SCALE.

        Returns None if no such t <= budget exists (so the bound is > budget).
        Raises PoolExhausted if it runs off the end of the pool first.
        """
        if budget <= 0:
            return None
        acc = s_hi
        t = 0
        mask, hi = self.mask, self.hi
        n = len(self.p)
        i = start
        while i < n:
            if mask[i] & tm == 0:
                acc += hi[i]
                t += 1
                if acc > SCALE:
                    return t
                if t >= budget:
                    return None
            i += 1
        # Pool exhausted.  The t primes just used are the largest reciprocals
        # available above ``start``, and they still do not reach 1, so any
        # completion needs strictly more than t of them: the bound is >= t + 1.
        # That is enough to answer "> budget" whenever t + 1 >= budget; if not,
        # the pool really is too small to decide and we must say so.
        if t + 1 >= budget:
            return None
        raise PoolExhausted(
            f"pool of odd primes < {self.limit} exhausted at t={t} "
            f"with budget {budget}; sum still {acc}/{SCALE}")

    def prod_min(self, tm: int, s_hi: int, start: int, cap: int):
        """(t, product) for the least t as in ``t_min``, product = the product
        of those t primes -- a lower bound on the product of any completion.

        ``cap`` bounds t so a hopeless branch does not run forever.  A capped
        branch returns (None, partial product) and increments ``cap_hits``; the
        partial product is still a valid lower bound, so the caller can prune
        on it, but a nonzero ``cap_hits`` means the reported bound depended on
        the cap and should be re-run with a larger one.
        """
        acc = s_hi
        t = 0
        prod = 1
        mask, hi, ps = self.mask, self.hi, self.p
        n = len(self.p)
        i = start
        while i < n:
            if mask[i] & tm == 0:
                acc += hi[i]
                prod *= ps[i]
                t += 1
                if acc > SCALE:
                    return t, prod
                if t >= cap:
                    self.cap_hits += 1
                    return None, prod
            i += 1
        # exhausted: >= t + 1 primes are needed, so prod * (next prime) is a
        # lower bound; use prod, which is weaker and still valid.
        self.cap_hits += 1
        return None, prod


def search(pool: Pool, best: int = 10 ** 9, prune: bool = True, node_cap=None):
    """Lower bound on |S| over admissible sets, by the split search above.

    Returns (bound, witness_T, nodes).  ``witness_T`` is the small-prime set of
    the leaf that realised the minimum -- the configuration the bound is
    weakest at, i.e. the one an improvement has to attack.
    """
    stats = {"nodes": 0}
    result = {"best": best, "witness": None}
    fx = pool.first_above_X

    def rec(T, tm, s_hi, li):
        stats["nodes"] += 1
        if node_cap is not None and stats["nodes"] > node_cap:
            raise TimeoutError("node cap reached")
        k = len(T)
        if prune:
            t = pool.t_min(tm, s_hi, li + 1, result["best"] - k)
            if t is None or k + t >= result["best"]:
                return
        # leaf: no further prime <= X, everything else lies above X
        tleaf = pool.t_min(tm, s_hi, fx, result["best"] - k)
        if tleaf is not None and k + tleaf < result["best"]:
            result["best"] = k + tleaf
            result["witness"] = (list(T), tleaf)
        for i in range(li + 1, pool.nbranch):
            if pool.mask[i] & tm:
                continue
            child_hi = s_hi + pool.hi[i]
            if prune:
                # R3: relaxed (filter by T only), monotone in i
                tr = pool.t_min(tm, child_hi, i + 1, result["best"] - k - 1)
                if tr is None or k + 1 + tr >= result["best"]:
                    break
            rec(T + [pool.p[i]], tm | (1 << i), child_hi, i)

    rec([], 0, 0, -1)
    return result["best"], result["witness"], stats["nodes"]


def search_product(pool: Pool, cap: int, best=None, node_cap=None):
    """Lower bound on prod(S) over admissible sets.  ``cap`` limits how many
    primes a leaf bound may use before the branch is abandoned as hopeless."""
    stats = {"nodes": 0}
    result = {"best": best, "witness": None}
    fx = pool.first_above_X

    def rec(T, tm, s_hi, li, prod):
        stats["nodes"] += 1
        if node_cap is not None and stats["nodes"] > node_cap:
            raise TimeoutError("node cap reached")
        t, tail = pool.prod_min(tm, s_hi, li + 1, cap)
        if result["best"] is not None and prod * tail >= result["best"]:
            return
        # A capped leaf still yields a valid (weaker) lower bound, so it must be
        # offered to the minimum -- skipping it could report a bound that is too
        # large, which is the unsound direction.
        tl, taill = pool.prod_min(tm, s_hi, fx, cap)
        cand = prod * taill
        if result["best"] is None or cand < result["best"]:
            result["best"] = cand
            result["witness"] = (list(T), tl)
        for i in range(li + 1, pool.nbranch):
            if pool.mask[i] & tm:
                continue
            child_hi = s_hi + pool.hi[i]
            tr, tailr = pool.prod_min(tm, child_hi, i + 1, cap)
            if result["best"] is not None and prod * pool.p[i] * tailr >= result["best"]:
                break
            rec(T + [pool.p[i]], tm | (1 << i), child_hi, i, prod * pool.p[i])

    rec([], 0, 0, -1, 1)
    return result["best"], result["witness"], stats["nodes"]


# --------------------------------------------------------------------------
# exact re-check of a single leaf, with Fractions rather than scaled integers
# --------------------------------------------------------------------------

def leaf_bound_exact(T, branch_limit, pool_limit):
    """Recompute |T| + t for one leaf using exact Fractions and plain division.

    Written to share nothing with Pool: no bitmask, no scaling, no sieve reuse.
    """
    ps = [q for q in primes_upto(pool_limit) if q > 2]
    s = sum((Fraction(1, p) for p in T), Fraction(0))
    t = 0
    for q in ps:
        if q <= branch_limit:
            continue
        if any((q - 1) % p == 0 for p in T):
            continue
        s += Fraction(1, q)
        t += 1
        if s > 1:
            return len(T) + t
    raise PoolExhausted("exact recheck ran out of primes")
