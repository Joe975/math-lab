#!/usr/bin/env python3
"""Second-order check of a numerical MUB-defect minimum.

    python problems/mub-six/explore/hessian.py problems/mub-six/data/best_d6k4_20k.json [out.json]

Coordinates: basis b >= 1 moves as U_b exp(A_b), A_b anti-Hermitian
(36 real coordinates per basis; basis 0 is fixed). At the input point it

1. Newton-polishes the configuration (pseudo-inverse on the non-gauge
   subspace) until the gradient is at roundoff;
2. builds the Hessian by central differences of the analytic gradient;
3. builds the gauge subspace explicitly (column phases of each moving basis,
   and a common left diagonal unitary, which leaves every |overlap| fixed)
   and reports its rank;
4. reports the spectrum of the Hessian on the orthogonal complement of the
   gauge subspace.

A critical point whose complement spectrum is strictly positive (by a margin
far above the finite-difference error) is a nondegenerate local minimum
modulo gauge -- numerically. Floats throughout: EVIDENCE, not a certificate.
Standard library only.
"""
from __future__ import annotations

import json
import math
import sys

def load(path):
    data = json.load(open(path))
    bases = data["bases"]
    d = len(bases[0])
    # U[b][r][c] = component r of vector c
    U = [[[complex(*bases[b][c][r]) for c in range(d)] for r in range(d)] for b in range(len(bases))]
    return d, U

def mm(A, B):
    n = len(A); m = len(B[0]); k = len(B)
    return [[sum(A[i][t] * B[t][j] for t in range(k)) for j in range(m)] for i in range(n)]

def dag(A):
    return [[A[j][i].conjugate() for j in range(len(A))] for i in range(len(A[0]))]

def expm(A):
    """exp of a small anti-Hermitian matrix: scaling and squaring + Taylor-18."""
    n = len(A)
    nrm = max(sum(abs(x) for x in row) for row in A)
    s = max(0, int(math.ceil(math.log2(nrm / 0.25))) if nrm > 0.25 else 0)
    B = [[x / (2 ** s) for x in row] for row in A]
    E = [[complex(i == j) for j in range(n)] for i in range(n)]
    T = [row[:] for row in E]
    for k in range(1, 19):
        T = [[x / k for x in row] for row in mm(T, B)]
        E = [[E[i][j] + T[i][j] for j in range(n)] for i in range(n)]
    for _ in range(s):
        E = mm(E, E)
    return E

def generators(d):
    """Real basis of anti-Hermitian d x d matrices (d^2 of them)."""
    gens = []
    for i in range(d):
        G = [[0j] * d for _ in range(d)]; G[i][i] = 1j; gens.append(G)
    for i in range(d):
        for j in range(i + 1, d):
            G = [[0j] * d for _ in range(d)]; G[i][j] = 1; G[j][i] = -1; gens.append(G)
            G = [[0j] * d for _ in range(d)]; G[i][j] = 1j; G[j][i] = 1j; gens.append(G)
    return gens

def loss_and_egrad(U, d):
    K = len(U); inv = 1.0 / d; L = 0.0
    Ge = [[[0j] * d for _ in range(d)] for _ in range(K)]
    for a in range(K):
        for b in range(a + 1, K):
            M = mm(dag(U[a]), U[b])
            R = [[0j] * d for _ in range(d)]
            for i in range(d):
                for j in range(d):
                    e = abs(M[i][j]) ** 2 - inv
                    L += e * e
                    R[i][j] = 2 * e * M[i][j]
            UaR = mm(U[a], R); UbRh = mm(U[b], dag(R))
            for r in range(d):
                for c in range(d):
                    Ge[b][r][c] += UaR[r][c]
                    Ge[a][r][c] += UbRh[r][c]
    return L, Ge

def coord_grad(U, d, gens):
    """dL/dx_k for U_b -> U_b exp(sum x_k G_k), b >= 1."""
    L, Ge = loss_and_egrad(U, d)
    g = []
    for b in range(1, len(U)):
        W = mm(dag(Ge[b]), U[b])   # dL = 2 Re tr(Ge^H U E)
        for G in gens:
            g.append(2 * sum((W[i][j] * G[j][i]).real for i in range(d) for j in range(d)))
    return L, g

def move(U, x, d, gens):
    n = len(gens); V = [U[0]]
    for b in range(1, len(U)):
        xs = x[(b - 1) * n:b * n]
        A = [[sum(xs[k] * gens[k][i][j] for k in range(n)) for j in range(d)] for i in range(d)]
        V.append(mm(U[b], expm(A)))
    return V

def gauge_vectors(U, d, gens):
    """Tangent vectors (in x-coordinates) of the gauge orbit at U."""
    n = len(gens); K = len(U); vecs = []
    def coords_of(A):          # anti-Hermitian A -> coefficients on gens
        out = [A[i][i].imag for i in range(d)]
        for i in range(d):
            for j in range(i + 1, d):
                out.append(A[i][j].real); out.append(A[i][j].imag)
        return out
    # column phases of each moving basis: A = i e_c e_c^T on basis b
    for b in range(1, K):
        for c in range(d):
            v = [0.0] * (n * (K - 1)); v[(b - 1) * n + c] = 1.0; vecs.append(v)
    # common left diagonal phase D = exp(i diag(phi)): U_b -> D U_b = U_b exp(U_b^H i diag U_b)
    for r in range(d):
        v = []
        for b in range(1, K):
            P = [[0j] * d for _ in range(d)]; P[r][r] = 1j
            v += coords_of(mm(dag(U[b]), mm(P, U[b])))
        vecs.append(v)
    return vecs

def orthonormalize(vecs, tol=1e-9):
    basis = []
    for v in vecs:
        w = v[:]
        for _ in range(2):
            for q in basis:
                p = sum(a * b for a, b in zip(w, q)); w = [a - p * b for a, b in zip(w, q)]
        nrm = math.sqrt(sum(a * a for a in w))
        if nrm > tol:
            basis.append([a / nrm for a in w])
    return basis

def complement(basis, N):
    """Orthonormal basis of the orthogonal complement of span(basis) in R^N."""
    comp = orthonormalize(basis + [[float(i == j) for j in range(N)] for i in range(N)])
    return comp[len(basis):]

def jacobi_eig(S, tol=1e-15, sweeps=60):
    n = len(S); A = [row[:] for row in S]
    V = [[float(i == j) for j in range(n)] for i in range(n)]
    for _ in range(sweeps):
        off = sum(A[i][j] ** 2 for i in range(n) for j in range(n) if i != j)
        if off < tol ** 2:
            break
        for p in range(n):
            for q in range(p + 1, n):
                if abs(A[p][q]) < 1e-300:
                    continue
                th = (A[q][q] - A[p][p]) / (2 * A[p][q])
                t = (1 if th >= 0 else -1) / (abs(th) + math.sqrt(th * th + 1))
                c = 1 / math.sqrt(t * t + 1); s = t * c
                for k in range(n):
                    akp, akq = A[k][p], A[k][q]
                    A[k][p] = c * akp - s * akq; A[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = A[p][k], A[q][k]
                    A[p][k] = c * apk - s * aqk; A[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p] = c * vkp - s * vkq; V[k][q] = s * vkp + c * vkq
    return [A[i][i] for i in range(n)], V

def hessian(U, d, gens, h=1e-5):
    N = len(gens) * (len(U) - 1); H = []
    for j in range(N):
        e = [0.0] * N; e[j] = h
        _, gp = coord_grad(move(U, e, d, gens), d, gens)
        e[j] = -h
        _, gm = coord_grad(move(U, e, d, gens), d, gens)
        H.append([(a - b) / (2 * h) for a, b in zip(gp, gm)])
    return [[0.5 * (H[i][j] + H[j][i]) for j in range(N)] for i in range(N)]

def project(H, Q):
    HQ = [[sum(H[i][k] * q[k] for k in range(len(q))) for i in range(len(H))] for q in Q]
    return [[sum(a * b for a, b in zip(Q[r], HQ[c])) for c in range(len(Q))] for r in range(len(Q))]

def main():
    d, U = load(sys.argv[1])
    gens = generators(d)
    N = len(gens) * (len(U) - 1)
    L, g = coord_grad(U, d, gens)
    print(f"input: L = {L:.17g}, |grad| = {math.sqrt(sum(x*x for x in g)):.3e}, coordinates N = {N}")

    for it in range(4):                      # Newton polish on the gauge complement
        G = orthonormalize(gauge_vectors(U, d, gens))
        Q = complement(G, N)
        H = hessian(U, d, gens)
        Hc = project(H, Q)
        lam, V = jacobi_eig(Hc)
        gq = [sum(a * b for a, b in zip(q, g)) for q in Q]
        # step = -sum_i (v_i . gq / lam_i) v_i, in Q-coords, then back to x
        m = len(Q); y = [0.0] * m
        for i in range(m):
            vi = [V[k][i] for k in range(m)]
            coef = sum(a * b for a, b in zip(vi, gq)) / lam[i]
            y = [a - coef * b for a, b in zip(y, vi)]
        step = [sum(y[r] * Q[r][k] for r in range(m)) for k in range(N)]
        U = move(U, step, d, gens)
        L, g = coord_grad(U, d, gens)
        gn = math.sqrt(sum(x * x for x in g))
        print(f"newton {it}: L = {L:.17g}  |grad| = {gn:.3e}  gauge rank = {len(G)}  "
              f"min/max complement eig = {min(lam):.6g} / {max(lam):.6g}")
        if gn < 1e-14:
            break

    G = orthonormalize(gauge_vectors(U, d, gens))
    Q = complement(G, N)
    H = hessian(U, d, gens)
    # how flat is H along the gauge directions? (should be ~ FD error)
    gauge_leak = max(math.sqrt(sum(sum(H[i][k] * q[k] for k in range(N)) ** 2 for i in range(N))) for q in G)
    lam_all, _ = jacobi_eig(H)
    lam, _ = jacobi_eig(project(H, Q))
    lam.sort(); lam_all.sort()
    near0 = sum(1 for x in lam_all if abs(x) < 1e-6)
    print(f"gauge dimension = {len(G)}; max |H v| over gauge v = {gauge_leak:.2e}")
    print(f"full Hessian: {near0} eigenvalues with |lambda| < 1e-6; smallest 30: "
          + " ".join(f"{x:.3g}" for x in lam_all[:30]))
    print(f"complement ({len(Q)} dims): min eig = {lam[0]:.6g}, max eig = {lam[-1]:.6g}")
    print("complement spectrum (sorted): " + " ".join(f"{x:.5g}" for x in lam))
    if len(sys.argv) > 2:
        out = {"d": d, "k": len(U), "L": L,
               "bases": [[[[U[b][r][c].real, U[b][r][c].imag] for r in range(d)] for c in range(d)] for b in range(len(U))],
               "hessian_complement_eigs": lam, "gauge_dimension": len(G),
               "full_hessian_eigs": lam_all,
               "convention": "bases[b][c] is basis vector c of basis b, entries [re, im]"}
        json.dump(out, open(sys.argv[2], "w"))

if __name__ == "__main__":
    main()
