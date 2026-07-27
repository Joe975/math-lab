"""Strongly regular graph tools in exact integer arithmetic.

A strongly regular graph SRG(n, k, lam, mu) is a k-regular graph on n vertices
in which adjacent vertices have exactly `lam` common neighbours and
non-adjacent vertices have exactly `mu`.

Everything here is integer/rational.  No floating point is used anywhere in a
decision: square roots go through `exact_isqrt`, which returns None unless the
argument is a perfect square, and every inequality is between integers.

Two independent things live here:

  * `feasibility(n, k, lam, mu)` -- the standard admissibility conditions on a
    *parameter set*, decided exactly.
  * `check_srg(adj, ...)` -- a direct check that an *explicit graph* really is
    strongly regular with the given parameters, by counting common
    neighbourhoods for every pair.

Graphs are represented as a list of int bitmasks: `adj[v]` has bit u set iff
u ~ v.  `n` is implicit in `len(adj)`.
"""

from itertools import combinations


# --------------------------------------------------------------------------
# exact square root
# --------------------------------------------------------------------------

def exact_isqrt(m):
    """Return t >= 0 with t*t == m, or None if m is not a perfect square."""
    if m < 0:
        return None
    t = 0
    # integer sqrt, then confirm.  math.isqrt is exact for ints.
    import math
    t = math.isqrt(m)
    return t if t * t == m else None


# --------------------------------------------------------------------------
# parameter feasibility
# --------------------------------------------------------------------------

def is_sum_of_two_squares(m):
    """Exact: is m = a^2 + b^2 for integers a, b?

    Decided by the classical criterion -- every prime factor congruent to
    3 mod 4 occurs to an even power -- computed by trial division.
    """
    if m < 0:
        return False
    d = 2
    while d * d <= m:
        if m % d == 0:
            e = 0
            while m % d == 0:
                m //= d
                e += 1
            if d % 4 == 3 and e % 2:
                return False
        d += 1
    return not (m % 4 == 3)


def basic_identity(n, k, lam, mu):
    """k(k - lam - 1) == (n - k - 1) mu, the standard counting identity."""
    return k * (k - lam - 1) == (n - k - 1) * mu


def spectrum(n, k, lam, mu):
    """Exact spectrum of SRG(n, k, lam, mu).

    Returns a dict with r, s (the restricted eigenvalues) and f, g (their
    multiplicities), or a dict with 'ok': False and a 'reason'.

    r, s are the roots of x^2 - (lam - mu) x - (k - mu) = 0.  The
    multiplicities are pinned by two exact linear equations:
        f + g       = n - 1          (dimension count)
        k + f r + g s = 0            (trace of the adjacency matrix is 0)
    We solve that system in integers rather than quoting the closed form, so
    the answer is a re-derivation rather than a transcription.
    """
    disc = (lam - mu) ** 2 + 4 * (k - mu)
    t = exact_isqrt(disc)
    conference = (2 * k + (n - 1) * (lam - mu) == 0)

    if conference:
        # r, s are irrational conjugates; f = g = (n-1)/2 is forced.
        if (n - 1) % 2:
            return {"ok": False, "reason": "conference case but n-1 odd"}
        return {
            "ok": True, "conference": True, "disc": disc,
            "r": None, "s": None, "f": (n - 1) // 2, "g": (n - 1) // 2,
        }

    if t is None:
        return {"ok": False, "reason": f"discriminant {disc} is not a perfect square",
                "disc": disc, "conference": False}

    # roots of x^2 - (lam-mu)x - (k-mu): ((lam-mu) +- t)/2
    num_r = (lam - mu) + t
    num_s = (lam - mu) - t
    if num_r % 2 or num_s % 2:
        return {"ok": False, "reason": "restricted eigenvalues are not integers",
                "disc": disc, "conference": False}
    r, s = num_r // 2, num_s // 2

    # f + g = n-1 ; k + f r + g s = 0  =>  f (r - s) = -k - (n-1) s
    denom = r - s
    num = -k - (n - 1) * s
    if denom == 0 or num % denom:
        return {"ok": False, "reason": "multiplicities are not integers",
                "disc": disc, "r": r, "s": s, "conference": False}
    f = num // denom
    g = (n - 1) - f
    if f < 0 or g < 0:
        return {"ok": False, "reason": "negative multiplicity",
                "disc": disc, "r": r, "s": s, "f": f, "g": g, "conference": False}
    # re-check the trace equation exactly rather than trusting the algebra
    assert f + g == n - 1
    assert k + f * r + g * s == 0
    return {"ok": True, "conference": False, "disc": disc,
            "r": r, "s": s, "f": f, "g": g}


def krein_conditions(k, r, s):
    """The two Krein inequalities, as exact integer comparisons.

        (r+1)(k + r + 2rs) <= (k + r)(s + 1)^2
        (s+1)(k + s + 2rs) <= (k + s)(r + 1)^2
    """
    c1 = (r + 1) * (k + r + 2 * r * s) <= (k + r) * (s + 1) ** 2
    c2 = (s + 1) * (k + s + 2 * r * s) <= (k + s) * (r + 1) ** 2
    return c1, c2


def absolute_bound(n, f, g):
    """n <= f(f+3)/2 and n <= g(g+3)/2 (multiplicities of the two eigenspaces)."""
    return 2 * n <= f * (f + 3), 2 * n <= g * (g + 3)


def feasibility(n, k, lam, mu):
    """Decide the standard feasibility conditions exactly.

    Returns (ok, report).  `report` lists every condition with its verdict, so
    a caller can say *which* condition killed a parameter set.
    """
    conds = []

    def add(name, ok, detail=""):
        conds.append({"name": name, "ok": bool(ok), "detail": detail})

    add("ranges", 0 < k < n and 0 <= lam <= k - 1 and 0 <= mu <= k,
        f"n={n} k={k} lam={lam} mu={mu}")
    add("counting-identity", basic_identity(n, k, lam, mu),
        f"k(k-lam-1)={k*(k-lam-1)} vs (n-k-1)mu={(n-k-1)*mu}")

    spec = spectrum(n, k, lam, mu)
    add("integral-spectrum", spec["ok"], spec.get("reason", ""))

    if spec["ok"] and not spec.get("conference"):
        r, s, f, g = spec["r"], spec["s"], spec["f"], spec["g"]
        c1, c2 = krein_conditions(k, r, s)
        add("krein-1", c1, f"r={r} s={s}")
        add("krein-2", c2, f"r={r} s={s}")
        a1, a2 = absolute_bound(n, f, g)
        add("absolute-bound-f", a1, f"n={n} f={f} f(f+3)/2={f*(f+3)//2}")
        add("absolute-bound-g", a2, f"n={n} g={g} g(g+3)/2={g*(g+3)//2}")

    if spec["ok"] and spec.get("conference"):
        # A conference graph on n vertices yields a symmetric conference matrix
        # of order n+1; Belevitch / Bruck-Ryser-Chowla then force n = 1 mod 4
        # and n to be a sum of two integer squares.
        add("conference-n-mod-4", n % 4 == 1, f"n={n}")
        add("conference-sum-two-squares", is_sum_of_two_squares(n), f"n={n}")

    if lam == 1:
        # a vertex neighbourhood is a graph on k vertices that is 1-regular,
        # i.e. a perfect matching, so k must be even.
        add("lam1-even-degree", k % 2 == 0, f"k={k}")

    ok = all(c["ok"] for c in conds)
    return ok, {"params": (n, k, lam, mu), "spectrum": spec, "conditions": conds}


def failed_conditions(report):
    return [c["name"] for c in report["conditions"] if not c["ok"]]


# --------------------------------------------------------------------------
# enumerating a parameter family
# --------------------------------------------------------------------------

def enumerate_feasible(lam, mu, n_max):
    """All feasible parameter sets with the given lam, mu and n <= n_max.

    Iterates over k and solves the counting identity for n exactly, rather
    than scanning (n, k) pairs.
    """
    out = []
    k = 1
    while True:
        # k(k-lam-1) = (n-k-1) mu  =>  n = k(k-lam-1)/mu + k + 1
        num = k * (k - lam - 1)
        if mu == 0:
            break
        if num % mu == 0:
            n = num // mu + k + 1
            if n > n_max:
                break
            ok, rep = feasibility(n, k, lam, mu)
            if ok:
                out.append((n, k, lam, mu, rep))
        else:
            # n would not be an integer; still need a stopping rule
            n_lo = num // mu + k + 1
            if n_lo > n_max:
                break
        k += 1
        if k > 4 * n_max:
            break
    return out


# --------------------------------------------------------------------------
# checking an explicit graph
# --------------------------------------------------------------------------

def check_srg(adj, n=None, k=None, lam=None, mu=None):
    """Verify an explicit graph is strongly regular; return its parameters.

    Counts |N(u) & N(v)| for every pair by direct bitmask intersection, i.e.
    from the definition, not from the spectrum.  If k/lam/mu are given they
    are asserted rather than inferred.
    """
    N = len(adj)
    if n is not None and n != N:
        return False, f"n mismatch: {N} != {n}"
    # simple + irreflexive + symmetric
    for v in range(N):
        if adj[v] >> v & 1:
            return False, f"vertex {v} has a loop"
        for u in range(N):
            if (adj[v] >> u & 1) != (adj[u] >> v & 1):
                return False, f"asymmetric at {u},{v}"
    degs = {bin(adj[v]).count("1") for v in range(N)}
    if len(degs) != 1:
        return False, f"not regular: degrees {sorted(degs)}"
    kk = degs.pop()
    lams, mus = set(), set()
    for u, v in combinations(range(N), 2):
        c = bin(adj[u] & adj[v]).count("1")
        if adj[u] >> v & 1:
            lams.add(c)
        else:
            mus.add(c)
    if len(lams) != 1:
        return False, f"adjacent common-neighbour counts vary: {sorted(lams)}"
    if len(mus) != 1:
        return False, f"non-adjacent common-neighbour counts vary: {sorted(mus)}"
    ll, mm = lams.pop(), mus.pop()
    got = (N, kk, ll, mm)
    if k is not None and (kk, ll, mm) != (k, lam, mu):
        return False, f"parameters {got} != {(n, k, lam, mu)}"
    return True, got


def neighbourhood_components(adj, v):
    """Multiset of connected-component sizes of the subgraph induced on N(v).

    For lam = 1 this must be {2: k/2}: the neighbourhood is a perfect matching.
    """
    nb = [u for u in range(len(adj)) if adj[v] >> u & 1]
    idx = {u: i for i, u in enumerate(nb)}
    seen = set()
    sizes = []
    for u in nb:
        if u in seen:
            continue
        stack, comp = [u], []
        seen.add(u)
        while stack:
            w = stack.pop()
            comp.append(w)
            for x in nb:
                if x not in seen and (adj[w] >> x & 1):
                    seen.add(x)
                    stack.append(x)
        sizes.append(len(comp))
    return sorted(sizes)


def degree_sequence_in(adj, vs):
    """Degrees of the subgraph induced on the vertex set `vs`."""
    mask = 0
    for v in vs:
        mask |= 1 << v
    return sorted(bin(adj[v] & mask).count("1") for v in vs)
