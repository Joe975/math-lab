"""Exact check of the character-block multiplicity lemma on graphs that exist.

Let tau be an automorphism of prime order p of a graph with adjacency A,
with t orbits of size p and f fixed points. C^n splits into the isotypic
components of <tau>: the trivial one (dimension t + f, on which A acts as the
quotient / orbit matrix R) and, for each non-trivial character chi = zeta^u,
one of dimension t on which A acts as

    M^(u)_ij = sum over m in Z_p with o_i ~ tau^m(o_j) of zeta^(u m)

(o_i a representative of orbit i). Hence for every eigenvalue theta

    mult_A(theta) = mult_R(theta) + sum_u mult_{M^(u)}(theta),

and since M^(u) is the Galois conjugate of M^(1), all p-1 non-trivial terms
are equal:  mult_A(theta) = a_theta + (p-1) r_theta.  That identity is what
attempt 003 uses on the order-7 case of the 99-graph (p = 7, f = 1, t = 14)
to restrict the trace of the orbit matrix. This script verifies it, in exact
arithmetic over Q(zeta_p), on four graphs with a Z_3 automorphism -- two
semiregular and two with exactly one fixed point, the shape the 99 case has.

    python problems/conway-99/explore/char_blocks.py
"""

import os
import sys
from fractions import Fraction
from itertools import combinations

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))

import constructions as C  # noqa: E402
import srg  # noqa: E402


# --- Q(zeta_p) as vectors of length p-1 over Q ------------------------------

class Cyclo:
    """Arithmetic in Q(zeta_p), basis 1, zeta, ..., zeta^(p-2)."""

    def __init__(self, p):
        self.p = p
        self.d = p - 1

    def zero(self):
        return [Fraction(0)] * self.d

    def one(self):
        return self.const(1)

    def const(self, c):
        v = self.zero()
        v[0] = Fraction(c)
        return v

    def zeta_pow(self, m):
        m %= self.p
        v = self.zero()
        if m < self.d:
            v[m] = Fraction(1)
        else:  # zeta^(p-1) = -(1 + zeta + ... + zeta^(p-2))
            for i in range(self.d):
                v[i] = Fraction(-1)
        return v

    def add(self, a, b):
        return [x + y for x, y in zip(a, b)]

    def sub(self, a, b):
        return [x - y for x, y in zip(a, b)]

    def scale(self, a, c):
        return [x * c for x in a]

    def mul(self, a, b):
        out = self.zero()
        for i, x in enumerate(a):
            if x == 0:
                continue
            for j, y in enumerate(b):
                if y == 0:
                    continue
                z = self.zeta_pow(i + j)
                for k in range(self.d):
                    out[k] += x * y * z[k]
        return out

    def is_zero(self, a):
        return all(x == 0 for x in a)

    def inv(self, a):
        """Solve a * x = 1 as a linear system over Q."""
        d = self.d
        cols = [self.mul(a, self.zeta_pow(j)) for j in range(d)]  # a * basis_j
        # matrix M with M[k][j] = cols[j][k]; solve M x = e_0
        M = [[cols[j][k] for j in range(d)] + [Fraction(1 if k == 0 else 0)]
             for k in range(d)]
        for c in range(d):
            piv = next(r for r in range(c, d) if M[r][c] != 0)
            M[c], M[piv] = M[piv], M[c]
            pv = M[c][c]
            M[c] = [x / pv for x in M[c]]
            for r in range(d):
                if r != c and M[r][c] != 0:
                    f = M[r][c]
                    M[r] = [x - f * y for x, y in zip(M[r], M[c])]
        return [M[k][d] for k in range(d)]

    def rank(self, rows):
        """Rank of a matrix with entries in Q(zeta_p)."""
        A = [row[:] for row in rows]
        n, m = len(A), len(A[0]) if A else 0
        r = 0
        for c in range(m):
            piv = next((i for i in range(r, n) if not self.is_zero(A[i][c])), None)
            if piv is None:
                continue
            A[r], A[piv] = A[piv], A[r]
            inv = self.inv(A[r][c])
            A[r] = [self.mul(x, inv) for x in A[r]]
            for i in range(n):
                if i != r and not self.is_zero(A[i][c]):
                    f = A[i][c]
                    A[i] = [self.sub(x, self.mul(f, y)) for x, y in zip(A[i], A[r])]
            r += 1
            if r == n:
                break
        return r


def rank_q(rows):
    A = [[Fraction(x) for x in row] for row in rows]
    n, m = len(A), len(A[0])
    r = 0
    for c in range(m):
        piv = next((i for i in range(r, n) if A[i][c] != 0), None)
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        pv = A[r][c]
        A[r] = [x / pv for x in A[r]]
        for i in range(n):
            if i != r and A[i][c] != 0:
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        r += 1
        if r == n:
            break
    return r


# --- graphs and automorphisms ---------------------------------------------

def orbits_of(perm, n):
    seen, orbs = [False] * n, []
    for v in range(n):
        if seen[v]:
            continue
        o, w = [], v
        while not seen[w]:
            seen[w] = True
            o.append(w)
            w = perm[w]
        orbs.append(o)
    return orbs


def is_automorphism(adj, perm):
    n = len(adj)
    return all((adj[u] >> v & 1) == (adj[perm[u]] >> perm[v] & 1)
               for u in range(n) for v in range(n))


def orbit_matrix(adj, orbs):
    t = len(orbs)
    R = [[0] * t for _ in range(t)]
    for i, Oi in enumerate(orbs):
        for j, Oj in enumerate(orbs):
            vals = {sum(1 for w in Oj if adj[u] >> w & 1) for u in Oi}
            assert len(vals) == 1, "not equitable"
            R[i][j] = vals.pop()
    return R


def character_block(adj, orbs, p, u, K):
    """M^(u) for the size-p orbits, over Q(zeta_p)."""
    big = [o for o in orbs if len(o) == p]
    t = len(big)
    M = [[K.zero() for _ in range(t)] for _ in range(t)]
    for i, Oi in enumerate(big):
        for j, Oj in enumerate(big):
            for m in range(p):
                if adj[Oi[0]] >> Oj[m] & 1:  # Oj[m] = tau^m(o_j)
                    M[i][j] = K.add(M[i][j], K.zeta_pow(u * m))
    return M


def multiplicities(adj, perm, p, thetas):
    n = len(adj)
    assert is_automorphism(adj, perm)
    orbs = orbits_of(perm, n)
    assert {len(o) for o in orbs} <= {1, p}
    R = orbit_matrix(adj, orbs)
    K = Cyclo(p)
    blocks = [character_block(adj, orbs, p, u, K) for u in range(1, p)]
    t = len(blocks[0])
    out = {}
    for theta in thetas:
        a = len(R) - rank_q([[R[i][j] - (theta if i == j else 0) for j in range(len(R))]
                             for i in range(len(R))])
        rs = []
        for M in blocks:
            shifted = [[K.sub(M[i][j], K.const(theta)) if i == j else M[i][j]
                        for j in range(t)] for i in range(t)]
            rs.append(t - K.rank(shifted))
        out[theta] = (a, rs)
    return R, out


def rook(m):
    adj = [0] * (m * m)
    for a in range(m):
        for b in range(m):
            for c in range(m):
                if c != b:
                    adj[a * m + b] |= 1 << (a * m + c)
                if c != a:
                    adj[a * m + b] |= 1 << (c * m + b)
    return adj


def petersen():
    vs = list(combinations(range(5), 2))
    idx = {v: i for i, v in enumerate(vs)}
    adj = [0] * 10
    for x in vs:
        for y in vs:
            if not set(x) & set(y):
                adj[idx[x]] |= 1 << idx[y]
    return adj, vs, idx


def cases():
    yield "paley9 / Z3 (semiregular)", C.paley9(), [(v + 3) % 9 for v in range(9)], 3
    yield "rook3 / Z3 (semiregular)", rook(3), [((v // 3 + 1) % 3) * 3 + v % 3 for v in range(9)], 3
    sig = {0: 1, 1: 2, 2: 0, 3: 3}
    yield ("rook4 / Z3 (one fixed point)", rook(4),
           [sig[v // 4] * 4 + sig[v % 4] for v in range(16)], 3)
    adj, vs, idx = petersen()
    yield ("petersen / Z3 (one fixed point)", adj,
           [idx[tuple(sorted(sig[x] if x < 3 else x for x in v))] for v in vs], 3)


def main():
    ok_all = True
    for name, adj, perm, p in cases():
        good, params = srg.check_srg(adj)
        assert good, params
        n, k, lam, mu = params
        spec = srg.spectrum(n, k, lam, mu)
        disc = (lam - mu) ** 2 + 4 * (k - mu)
        sq = int(round(disc ** 0.5))
        assert sq * sq == disc and (lam - mu + sq) % 2 == 0, "integral spectrum needed"
        r, s_ = (lam - mu + sq) // 2, (lam - mu - sq) // 2
        thetas = [k, r, s_]
        mults = {k: 1, r: spec["f"], s_: spec["g"]}
        R, out = multiplicities(adj, perm, p, thetas)
        fixed = sum(1 for o in orbits_of(perm, n) if len(o) == 1)
        print(f"{name}: SRG{params}, fixed points {fixed}, quotient {len(R)}x{len(R)}")
        for theta, (a, rs) in out.items():
            lhs = mults[theta]
            rhs = a + sum(rs)
            same = len(set(rs)) == 1
            flag = "ok" if (lhs == rhs and same) else "MISMATCH"
            ok_all &= flag == "ok"
            print(f"  theta={theta}: mult {lhs} = a {a} + sum r {rs} -> {rhs}  [{flag}]")
    print("ALL OK" if ok_all else "FAILURE")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
