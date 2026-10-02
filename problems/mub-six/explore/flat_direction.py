#!/usr/bin/env python3
"""Is an extra Hessian null direction a true valley or a higher-order bump?

    python problems/mub-six/explore/flat_direction.py <config.json>

At a polished critical point, take the Hessian null vectors that are NOT gauge
directions. For each step t along such a vector v, move to x = t v and then
re-minimise L over the gauge complement minus v (Newton with the base-point
Hessian as preconditioner, gradient projected off v). Prints
relaxed L(t) - L*. A true one-parameter valley gives ~roundoff for all t;
a degenerate-but-isolated minimum gives growth like t^4 (or t^3 for a saddle).
Standard library only.
"""
from __future__ import annotations

import math
import sys

sys.path.insert(0, __file__.rsplit("/", 1)[0])
import hessian as H


def main():
    d, U = H.load(sys.argv[1])
    gens = H.generators(d); N = len(gens) * (len(U) - 1)

    def polish(U, fixed=(), steps=6):
        for _ in range(steps):
            G = H.orthonormalize(H.gauge_vectors(U, d, gens) + [list(f) for f in fixed])
            Q = H.complement(G, N)
            Hc = H.project(H.hessian(U, d, gens), Q)
            lam, V = H.jacobi_eig(Hc)
            L, g = H.coord_grad(U, d, gens)
            gq = [sum(a * b for a, b in zip(q, g)) for q in Q]
            m = len(Q); y = [0.0] * m
            for i in range(m):
                if abs(lam[i]) < 1e-6:
                    continue
                vi = [V[k][i] for k in range(m)]
                coef = sum(a * b for a, b in zip(vi, gq)) / lam[i]
                y = [a - coef * b for a, b in zip(y, vi)]
            U = H.move(U, [sum(y[r] * Q[r][k] for r in range(m)) for k in range(N)], d, gens)
        L, g = H.coord_grad(U, d, gens)
        return U, L

    U, Lstar = polish(U, steps=3)
    G = H.orthonormalize(H.gauge_vectors(U, d, gens))
    Q = H.complement(G, N)
    lam, V = H.jacobi_eig(H.project(H.hessian(U, d, gens), Q))
    null = [i for i in range(len(Q)) if abs(lam[i]) < 1e-6]
    print(f"L* = {Lstar:.17g}; gauge rank {len(G)}; non-gauge null directions: {len(null)}")
    for i in null:
        v = [sum(V[r][i] * Q[r][k] for r in range(len(Q))) for k in range(N)]
        for t in (1e-3, 1e-2, 3e-2, 1e-1, 3e-1):
            Ut = H.move(U, [t * x for x in v], d, gens)
            raw = H.coord_grad(Ut, d, gens)[0] - Lstar
            # relax: keep the component along v fixed by excluding the
            # parallel-transported direction (approximated by v itself)
            Ur, Lr = polish(Ut, fixed=[v], steps=4)
            print(f"  t = {t:<6g} raw L - L* = {raw:.3e}   relaxed L - L* = {Lr - Lstar:.3e}")


if __name__ == "__main__" and "--bargmann" not in sys.argv:
    main()


def invariants(U, d):
    """Sorted |overlap|^2 values per pair: unchanged by every gauge move."""
    out = []
    for a in range(len(U)):
        for b in range(a + 1, len(U)):
            M = H.mm(H.dag(U[a]), U[b])
            out.append(sorted(round(abs(M[i][j]) ** 2, 10) for i in range(d) for j in range(d)))
    return out


def bargmann(U, d):
    """Sorted Re and Im of all triple products <a_i|b_j><b_j|c_k><c_k|a_i>
    over the first three bases. Invariant under column phases and a common
    left unitary, so a change here means a genuinely different configuration
    (up to relabelling, which only permutes the multiset)."""
    A, B, Cm = U[0], U[1], U[2]
    ab = H.mm(H.dag(A), B); bc = H.mm(H.dag(B), Cm); ca = H.mm(H.dag(Cm), A)
    t = [ab[i][j] * bc[j][k] * ca[k][i] for i in range(d) for j in range(d) for k in range(d)]
    return sorted(round(z.real, 9) for z in t), sorted(round(z.imag, 9) for z in t)


def relax_along(U, v, t, d, gens, iters=4000):
    """Step t along v, then projected gradient descent off gauge + v."""
    Ut = H.move(U, [t * x for x in v], d, gens)
    L, g = H.coord_grad(Ut, d, gens); eta = 0.05; Gs = None
    for it in range(iters):
        if it % 200 == 0:
            Gs = H.orthonormalize(H.gauge_vectors(Ut, d, gens) + [v])
        gp = g[:]
        for q in Gs:
            p = sum(a * b for a, b in zip(gp, q)); gp = [a - p * b for a, b in zip(gp, q)]
        Un = H.move(Ut, [-eta * x for x in gp], d, gens); Ln, gn = H.coord_grad(Un, d, gens)
        if Ln < L:
            Ut, L, g = Un, Ln, gn; eta *= 1.2
        else:
            eta *= 0.5
        if math.sqrt(sum(x * x for x in gp)) < 1e-13:
            break
    return Ut, L


if __name__ == "__main__" and len(sys.argv) > 2 and sys.argv[2] == "--bargmann":
    d, U = H.load(sys.argv[1]); gens = H.generators(d); N = len(gens) * (len(U) - 1)
    G = H.orthonormalize(H.gauge_vectors(U, d, gens)); Q = H.complement(G, N)
    lam, V = H.jacobi_eig(H.project(H.hessian(U, d, gens), Q))
    i = min(range(len(Q)), key=lambda j: abs(lam[j]))
    v = [sum(V[r][i] * Q[r][k] for r in range(len(Q))) for k in range(N)]
    re0, im0 = bargmann(U, d)
    print(f"null eigenvalue used: {lam[i]:.2e}")
    for t in (0.1, 0.3, 0.6):
        Ut, L = relax_along(U, v, t, d, gens)
        re1, im1 = bargmann(Ut, d)
        dre = max(abs(a - b) for a, b in zip(re0, re1)); dim = max(abs(a - b) for a, b in zip(im0, im1))
        print(f"t = {t}: relaxed L - L* = {L - 0.0512492189962838278:.2e}; "
              f"max change in sorted triple products: Re {dre:.3e}, Im {dim:.3e}")
