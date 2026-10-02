#!/usr/bin/env python3
"""Polish a numerical MUB-defect critical point to high precision (stdlib).

    python problems/mub-six/explore/hiprec_polish.py problems/mub-six/data/polished_d6k4.json [digits]

The input must carry the float Hessian spectrum written by hessian.py; the
float Hessian is recomputed here and used as a fixed preconditioner (inverse on
the gauge complement), while the gradient, unitarity and the step are done in
decimal arithmetic at `digits` precision. Converges linearly with ratio
~ (float Hessian error) per iteration. Prints the distinct overlap values
|<a_i|b_j>|^2 and L at full precision.
"""
from __future__ import annotations

import json
import sys
from decimal import Decimal as D, getcontext

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import hessian as H  # float machinery: generators, gauge, complement, Hessian


class C:
    """Minimal complex number over Decimal."""
    __slots__ = ("re", "im")

    def __init__(self, re, im=D(0)):
        self.re, self.im = re, im

    def __add__(self, o): return C(self.re + o.re, self.im + o.im)
    def __sub__(self, o): return C(self.re - o.re, self.im - o.im)
    def __mul__(self, o):
        if isinstance(o, C):
            return C(self.re * o.re - self.im * o.im, self.re * o.im + self.im * o.re)
        return C(self.re * o, self.im * o)
    def conj(self): return C(self.re, -self.im)
    def abs2(self): return self.re * self.re + self.im * self.im


Z = lambda: C(D(0), D(0))


def mm(A, B):
    n, k, m = len(A), len(B), len(B[0])
    out = [[Z() for _ in range(m)] for _ in range(n)]
    for i in range(n):
        for t in range(k):
            a = A[i][t]
            for j in range(m):
                out[i][j] = out[i][j] + a * B[t][j]
    return out


def dag(A):
    return [[A[j][i].conj() for j in range(len(A))] for i in range(len(A[0]))]


def gram_schmidt(U):
    d = len(U)
    cols = [[U[r][c] for r in range(d)] for c in range(d)]
    out = []
    for v in cols:
        for q in out:
            ip = Z()
            for a, b in zip(q, v):
                ip = ip + a.conj() * b
            v = [b - a * ip for a, b in zip(q, v)]
        n = sum((x.abs2() for x in v), D(0)).sqrt()
        out.append([C(x.re / n, x.im / n) for x in v])
    return [[out[c][r] for c in range(d)] for r in range(d)]


def loss_grad(U, d, gens):
    inv = D(1) / d; L = D(0); K = len(U)
    Ge = [[[Z() for _ in range(d)] for _ in range(d)] for _ in range(K)]
    for a in range(K):
        for b in range(a + 1, K):
            M = mm(dag(U[a]), U[b])
            R = [[None] * d for _ in range(d)]
            for i in range(d):
                for j in range(d):
                    e = M[i][j].abs2() - inv
                    L += e * e
                    R[i][j] = M[i][j] * (2 * e)
            UaR = mm(U[a], R); UbRh = mm(U[b], dag(R))
            for r in range(d):
                for c in range(d):
                    Ge[b][r][c] = Ge[b][r][c] + UaR[r][c]
                    Ge[a][r][c] = Ge[a][r][c] + UbRh[r][c]
    g = []
    for b in range(1, K):
        W = mm(dag(Ge[b]), U[b])
        for G in gens:
            s = D(0)
            for i in range(d):
                for j in range(d):
                    gji = G[j][i]
                    if gji != 0:
                        s += W[i][j].re * D(gji.real) - W[i][j].im * D(gji.imag)
            g.append(2 * s)
    return L, g


def move(U, x, d, gens):
    """U_b exp(A_b) with A tiny: exp by Taylor to 6th order, then re-orthonormalise."""
    n = len(gens); V = [U[0]]
    for b in range(1, len(U)):
        xs = x[(b - 1) * n:b * n]
        A = [[Z() for _ in range(d)] for _ in range(d)]
        for k, G in enumerate(gens):
            if xs[k] == 0:
                continue
            for i in range(d):
                for j in range(d):
                    if G[i][j] != 0:
                        A[i][j] = A[i][j] + C(D(G[i][j].real), D(G[i][j].imag)) * xs[k]
        E = [[C(D(int(i == j))) for j in range(d)] for i in range(d)]
        T = [row[:] for row in E]
        for k in range(1, 7):
            T = [[x * (D(1) / k) for x in row] for row in mm(T, A)]
            E = [[E[i][j] + T[i][j] for j in range(d)] for i in range(d)]
        V.append(gram_schmidt(mm(U[b], E)))
    return V


def overlap_values(U, d):
    vals = []
    for a in range(len(U)):
        for b in range(a + 1, len(U)):
            M = mm(dag(U[a]), U[b])
            vals += [M[i][j].abs2() for i in range(d) for j in range(d)]
    return vals


def main():
    digits = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    getcontext().prec = digits + 10
    data = json.load(open(sys.argv[1]))
    d = len(data["bases"][0])
    Uf = [[[complex(*data["bases"][b][c][r]) for c in range(d)] for r in range(d)] for b in range(len(data["bases"]))]
    gens = H.generators(d)
    N = len(gens) * (len(Uf) - 1)
    # float preconditioner on the gauge complement
    G = H.orthonormalize(H.gauge_vectors(Uf, d, gens))
    Q = H.complement(G, N)
    lam, V = H.jacobi_eig(H.project(H.hessian(Uf, d, gens), Q))
    m = len(Q)
    U = [gram_schmidt([[C(D(z.real), D(z.imag)) for z in row] for row in M]) for M in Uf]
    for it in range(12):
        L, g = loss_grad(U, d, gens)
        gn = max(abs(x) for x in g)
        print(f"iter {it}: |grad|_inf = {float(gn):.3e}  L = {L}")
        if gn < D(10) ** (-digits + 2):
            break
        gq = [sum((D(q[k]) * g[k] for k in range(N)), D(0)) for q in Q]
        y = [D(0)] * m
        for i in range(m):
            vi = [D(V[k][i]) for k in range(m)]
            coef = sum((a * b for a, b in zip(vi, gq)), D(0)) / D(lam[i])
            y = [a - coef * b for a, b in zip(y, vi)]
        step = [sum((y[r] * D(Q[r][k]) for r in range(m)), D(0)) for k in range(N)]
        U = move(U, step, d, gens)
    vals = sorted(overlap_values(U, d))
    distinct = []
    for v in vals:
        if not distinct or abs(v - distinct[-1][0]) > D(10) ** (-digits // 2):
            distinct.append([v, 1])
        else:
            distinct[-1][1] += 1
    getcontext().prec = digits
    print(f"L = {+L}")
    for v, n in distinct:
        print(f"{+v}  (x{n})")
    return U, L, distinct


if __name__ == "__main__":
    main()
