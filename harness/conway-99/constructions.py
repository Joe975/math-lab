"""Explicit strongly regular graphs, built from scratch, for calibration.

Every graph here is generated from its definition (a field, a code, a set
system) rather than read from a table, so that a checker run against them is
testing the checker rather than testing a transcription.

Of particular interest are the two realised members of the lam=1, mu=2 family:

    paley9()   -- SRG(9, 4, 1, 2)
    bvls243()  -- SRG(243, 22, 1, 2), the coset graph of the perfect ternary
                  Golay code [11, 6, 5]_3

Graphs are lists of int bitmasks.
"""


def _from_pairs(n, pairs):
    adj = [0] * n
    for u, v in pairs:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return adj


# --------------------------------------------------------------------------
# Paley graphs
# --------------------------------------------------------------------------

def paley_prime(q):
    """Paley graph on GF(q), q prime, q = 1 mod 4.  SRG(q, (q-1)/2, (q-5)/4, (q-1)/4)."""
    assert q % 4 == 1
    sq = {(x * x) % q for x in range(1, q)}
    return _from_pairs(q, [(u, v) for u in range(q) for v in range(u + 1, q)
                           if (v - u) % q in sq])


def paley9():
    """SRG(9, 4, 1, 2) on GF(9) = GF(3)[i], i^2 = -1.

    Built as a Cayley graph on the additive group of GF(9) with connection set
    the four nonzero squares.
    """
    # elements a + b*i, a,b in GF(3); index = 3a + b
    def mul(x, y):
        a, b = divmod(x, 3)
        c, d = divmod(y, 3)
        # (a+bi)(c+di) = (ac - bd) + (ad + bc) i
        return ((a * c - b * d) % 3) * 3 + (a * d + b * c) % 3

    def sub(x, y):
        a, b = divmod(x, 3)
        c, d = divmod(y, 3)
        return ((a - c) % 3) * 3 + (b - d) % 3

    squares = {mul(x, x) for x in range(1, 9)}
    return _from_pairs(9, [(u, v) for u in range(9) for v in range(u + 1, 9)
                           if sub(u, v) in squares])


# --------------------------------------------------------------------------
# small classical graphs
# --------------------------------------------------------------------------

def petersen():
    """SRG(10, 3, 0, 1): Kneser graph on 2-subsets of a 5-set."""
    from itertools import combinations
    vs = list(combinations(range(5), 2))
    idx = {v: i for i, v in enumerate(vs)}
    return _from_pairs(10, [(idx[a], idx[b]) for a in vs for b in vs
                            if a < b and not set(a) & set(b)])


def triangular(m):
    """T(m) = Johnson/line graph of K_m: SRG(m(m-1)/2, 2(m-2), m-2, 4)."""
    from itertools import combinations
    vs = list(combinations(range(m), 2))
    idx = {v: i for i, v in enumerate(vs)}
    return _from_pairs(len(vs), [(idx[a], idx[b]) for a in vs for b in vs
                                 if a < b and set(a) & set(b)])


def rook(m):
    """m x m rook's graph (lattice graph L2(m)): SRG(m^2, 2(m-1), m-2, 2)."""
    n = m * m
    pairs = []
    for u in range(n):
        for v in range(u + 1, n):
            if (u // m == v // m) or (u % m == v % m):
                pairs.append((u, v))
    return _from_pairs(n, pairs)


def shrikhande():
    """SRG(16, 6, 2, 2): Cayley graph on Z4 x Z4 with connection set
    {+-(1,0), +-(0,1), +-(1,1)}.  Same parameters as rook(4) but not isomorphic
    -- the standard test pair for an isomorphism check."""
    conn = [(1, 0), (3, 0), (0, 1), (0, 3), (1, 1), (3, 3)]
    idx = lambda a, b: 4 * (a % 4) + (b % 4)
    pairs = set()
    for a in range(4):
        for b in range(4):
            for (c, d) in conn:
                u, v = idx(a, b), idx(a + c, b + d)
                pairs.add((min(u, v), max(u, v)))
    return _from_pairs(16, sorted(pairs))


def clebsch():
    """SRG(16, 5, 0, 2): folded 5-cube.  Vertices = even-weight subsets of a
    5-set modulo complement; simpler model: GF(2)^4 with connection set the
    four unit vectors and the all-ones vector."""
    conn = [1, 2, 4, 8, 15]
    pairs = set()
    for x in range(16):
        for c in conn:
            u, v = x, x ^ c
            pairs.add((min(u, v), max(u, v)))
    return _from_pairs(16, sorted(pairs))


# --------------------------------------------------------------------------
# the perfect ternary Golay code and its coset graph
# --------------------------------------------------------------------------

def _poly_divides_xn_minus_1(g, n, p=3):
    """Does monic g (list of coeffs, low to high) divide x^n - 1 over GF(p)?"""
    # long division of x^n - 1 by g
    rem = [0] * n + [1]
    rem[0] = (-1) % p
    dg = len(g) - 1
    for i in range(len(rem) - 1, dg - 1, -1):
        c = rem[i]
        if c:
            for j in range(dg + 1):
                rem[i - dg + j] = (rem[i - dg + j] - c * g[j]) % p
    return all(c == 0 for c in rem)


def ternary_golay_generator():
    """Find a degree-5 monic divisor of x^11 - 1 over GF(3) by exhaustive search.

    x^11 - 1 factors over GF(3) as (x-1) times two irreducible quintics
    (because 3 has multiplicative order 5 mod 11); either quintic generates a
    perfect [11, 6, 5]_3 code.  We locate one by brute force over all 3^5
    monic quintics rather than quoting a generator polynomial.
    """
    for code in range(3 ** 5):
        c, coeffs = code, []
        for _ in range(5):
            coeffs.append(c % 3)
            c //= 3
        g = coeffs + [1]          # monic, degree 5
        if _poly_divides_xn_minus_1(g, 11):
            return g
    raise RuntimeError("no quintic divisor found")


def ternary_golay_code():
    """All 3^6 = 729 codewords of a perfect [11, 6, 5]_3 code, as tuples."""
    g = ternary_golay_generator()
    n, k = 11, 6
    # generator matrix: cyclic shifts of g
    rows = []
    for i in range(k):
        row = [0] * n
        for j, c in enumerate(g):
            row[i + j] = c
        rows.append(row)
    code = set()
    for m in range(3 ** k):
        x, coef = m, []
        for _ in range(k):
            coef.append(x % 3)
            x //= 3
        w = [0] * n
        for ci, row in zip(coef, rows):
            if ci:
                for j in range(n):
                    w[j] = (w[j] + ci * row[j]) % 3
        code.add(tuple(w))
    return code


def _weight(v):
    return sum(1 for c in v if c)


def bvls243():
    """SRG(243, 22, 1, 2): the coset graph of the perfect ternary Golay code.

    Vertices are the 3^5 = 243 cosets of the [11, 6, 5]_3 code in GF(3)^11.
    Since the code is perfect with covering radius 2, every coset has a unique
    representative of weight <= 2, and we use those as vertex names.  Two
    cosets are adjacent iff they differ by a vector of weight 1; there are
    11 * 2 = 22 such vectors, giving degree 22.

    This is the Berlekamp-van Lint-Seidel graph.  Returns (adj, reps).
    """
    from itertools import combinations, product
    code = ternary_golay_code()
    n = 11
    reps = [tuple([0] * n)]
    for i in range(n):
        for a in (1, 2):
            v = [0] * n
            v[i] = a
            reps.append(tuple(v))
    for i, j in combinations(range(n), 2):
        for a, b in product((1, 2), repeat=2):
            v = [0] * n
            v[i], v[j] = a, b
            reps.append(tuple(v))
    assert len(reps) == 243, len(reps)
    idx = {}
    for r in reps:
        for c in code:
            idx[tuple((ri + ci) % 3 for ri, ci in zip(r, c))] = r
    assert len(idx) == 3 ** n, (len(idx), 3 ** n)   # confirms the code is perfect

    rnum = {r: i for i, r in enumerate(reps)}
    pairs = set()
    for r in reps:
        for i in range(n):
            for a in (1, 2):
                w = list(r)
                w[i] = (w[i] + a) % 3
                other = idx[tuple(w)]
                u, v = rnum[r], rnum[other]
                if u != v:
                    pairs.add((min(u, v), max(u, v)))
    return _from_pairs(243, sorted(pairs)), reps


CATALOGUE = {
    "paley9": (lambda: paley9(), (9, 4, 1, 2)),
    "paley13": (lambda: paley_prime(13), (13, 6, 2, 3)),
    "paley17": (lambda: paley_prime(17), (17, 8, 3, 4)),
    "petersen": (lambda: petersen(), (10, 3, 0, 1)),
    "triangular6": (lambda: triangular(6), (15, 8, 4, 4)),
    "rook3": (lambda: rook(3), (9, 4, 1, 2)),
    "rook4": (lambda: rook(4), (16, 6, 2, 2)),
    "shrikhande": (lambda: shrikhande(), (16, 6, 2, 2)),
    "clebsch": (lambda: clebsch(), (16, 5, 0, 2)),
    "bvls243": (lambda: bvls243()[0], (243, 22, 1, 2)),
}
