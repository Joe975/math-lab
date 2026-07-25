#!/usr/bin/env python3
"""uc_weighted_kl.py -- evaluation of the weighted-KL ladder (Idea B of
attempts/union-closed/001-entropy-barrier-map.md; verdict recorded in
attempts/union-closed/002-weighted-kl-ladder.md).

Idea B proposed the family of distributional statements, for A,B i.i.d. from a
distribution mu on subsets of [n] and U = A u B:

  (S_c) at threshold p:  all marginals Pr[i in A] < p and H(A) > 0
                         ==>  H(U) + c*D(U || mu) > H(A).

For mu = Unif(F), F exactly union-closed, the identity
H(X) + D(X||Unif(F)) = log|F| (any X supported on F) gives
H(U) + c*D = log|F| - (1-c)*D <= log|F| = H(A), so (S_c) at threshold p
implies Frankl with constant p.  Since f_mu(c) := H(U) + c*D - H(A) is LINEAR
in c with slope D >= 0, each violating mu is summarized by one number

  c*(mu) = (H(A) - H(U)) / D(U||mu)     (when H(U) < H(A), 0 < D < inf):

mu violates (S_c) exactly for c <= c*(mu), at every threshold p > max marginal.

This script computes:
  A. the Bernoulli-product obstruction curve p(c) (the "ladder"), verifying
     the table in attempt 001, and the weight c needed to beat Liu's 0.38271;
  B. the c-profile of Ellis's n=2 counterexample (arXiv:2211.12401) to
     Gilmer's Conjecture 1 -- refutes only c ~ 1 at marginals ~ 1/2;
  C. EXACT finite-n evaluation of Sawin's Proposition 6 distribution
     (arXiv:2211.11504): k ~ Geometric(theta) truncated, then
     A|k ~ Bernoulli(1-(1-ubar)^(k+1))^n.  The component class is closed
     under union (k' = kA + kB + 1), so U is a mixture of the SAME products;
     both laws are exchangeable and everything is computed exactly from
     weight profiles.  Result: c*(mu) >> 1 at marginals arbitrarily close to
     psi = (3-sqrt5)/2 -- every rung c of the ladder is refuted at every
     threshold p > psi.
  D. an independent cross-check construction (slice at w0 ~ p*n smoothed by
     an exponentially small Bernoulli(q) product component), with CERTIFIED
     closed-form bounds F(c) >= f(c), c in [0,1], no exact law needed.

Conclusion (see part E summary): the true ceiling of the ladder is
P(c) = psi for every c >= 0.  Idea B is dead.

Usage: python3 uc_weighted_kl.py [--big]   (--big: larger n in part C/D)
Stdlib only, deterministic, ~seconds (default) / ~1 min (--big).
"""

import argparse
import math
from math import lgamma, log, log2

LN2 = log(2.0)
PSI = (3.0 - math.sqrt(5.0)) / 2.0        # 0.381966...
RECORD = 0.38271                          # Liu, arXiv:2306.08824 (ISIT 2024)


# ---------------------------------------------------------------- primitives

def h(p):
    """Binary entropy in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p)


def dkl(a, b):
    """KL divergence D(Bern(a)||Bern(b)) in bits."""
    if a <= 0.0:
        return -log2(1.0 - b) if b < 1.0 else float("inf")
    if a >= 1.0:
        return -log2(b) if b > 0.0 else float("inf")
    if b <= 0.0 or b >= 1.0:
        return float("inf")
    return a * log2(a / b) + (1.0 - a) * log2((1.0 - a) / (1.0 - b))


def log2C(n, w):
    """log2 of binomial coefficient C(n, w)."""
    return (lgamma(n + 1) - lgamma(w + 1) - lgamma(n - w + 1)) / LN2


def logsumexp2(xs):
    m = max(xs)
    if m == float("-inf"):
        return m
    return m + log2(sum(2.0 ** (x - m) for x in xs))


# ------------------------------------------- A. product obstruction curve p(c)

def c_star_product(p):
    """c*(Bern(p)^n) = (h(p)-h(q))/d(q||p), q = 2p-p^2.  Positive iff p>psi.
    This is the exact weight c at which the Bernoulli(p) product switches from
    supporting (S_c) to violating it."""
    q = 2.0 * p - p * p
    return (h(p) - h(q)) / dkl(q, p)


def p_of_c(c):
    """Product-obstruction threshold: the marginal level where Bernoulli
    products stop violating (S_c).  c_star_product is increasing on
    (psi, 1/2), from 0 to 1, so bisect."""
    if c <= 0.0:
        return PSI
    if c >= 1.0:
        return 0.5
    lo, hi = PSI + 1e-15, 0.5 - 1e-15
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if c_star_product(mid) < c:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def part_A():
    print("=" * 78)
    print("A. Bernoulli-product obstruction curve p(c)  [the 'ladder' of 001]")
    print("=" * 78)
    print("p(c) = first root of  h(2p-p^2) + c*d(2p-p^2||p) = h(p)")
    print(f"{'c':>8}   {'p(c)':>10}   residual g(c,p(c))")
    table_001 = {0.0: 0.38197, 0.1: 0.39264, 0.2: 0.40356, 0.3: 0.41473,
                 0.5: 0.43782, 0.7: 0.46192, 0.9: 0.48705}
    for c in [0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.99, 1.0]:
        p = p_of_c(c)
        q = 2 * p - p * p
        res = h(q) + c * dkl(q, p) - h(p)
        ref = table_001.get(c)
        note = ""
        if ref is not None:
            note = ("  == 001 table" if abs(p - ref) < 6e-6
                    else f"  MISMATCH vs 001 value {ref}")
        print(f"{c:8.3f}   {p:10.6f}   {res:+.2e}{note}")
    c_needed = c_star_product(RECORD)
    print(f"\nWeight needed for the product curve to beat Liu's record "
          f"{RECORD}:\n  c_needed = c*_prod({RECORD}) = {c_needed:.6f}")
    print("So even c ~ 0.007, IF (S_c) held at p(c), would beat the record.")
    print("Parts B-D show (S_c) is false above psi for EVERY c, so p(c) is a")
    print("ceiling for products only -- a mirage, not a ladder.")
    return c_needed


# --------------------------------- B. finite-support exact evaluation (Ellis)

def dist_eval(mu, n, c_list):
    """Exact H(A), H(U), D(U||A), marginals, c* for a distribution given as
    {bitmask: prob} on subsets of [n].  U-law by exact pair convolution."""
    nu = {}
    for s, ps in mu.items():
        for t, pt in mu.items():
            u = s | t
            nu[u] = nu.get(u, 0.0) + ps * pt
    if any(u not in mu or mu[u] <= 0 for u in nu):
        return None  # D infinite (support escape)
    HA = -sum(p * log2(p) for p in mu.values() if p > 0)
    HU = -sum(p * log2(p) for p in nu.values() if p > 0)
    D = sum(p * log2(p / mu[u]) for u, p in nu.items() if p > 0)
    marg = [sum(p for s, p in mu.items() if s >> i & 1) for i in range(n)]
    cstar = (HA - HU) / D if D > 0 else float("nan")
    fs = {c: HU + c * D - HA for c in c_list}
    return dict(HA=HA, HU=HU, D=D, marg=marg, cstar=cstar, f=fs)


def part_B(c_list):
    print("\n" + "=" * 78)
    print("B. Ellis's n=2 counterexample to Gilmer's Conjecture 1 "
          "(arXiv:2211.12401)")
    print("=" * 78)
    print("mu = {emptyset:0.3, {1}:0.2, {2}:0.2, {1,2}:0.3} and perturbations")
    print("moving delta of mass {1,2} -> emptyset (marginals 1/2 - delta).")
    print(f"{'delta':>7} {'max marg':>9} {'H(A)':>8} {'H(U)':>8} {'D':>8} "
          f"{'f(1)':>8} {'c*':>7}")
    for delta in [0.0, 0.01, 0.02, 0.05, 0.10]:
        mu = {0: 0.3 + delta, 1: 0.2, 2: 0.2, 3: 0.3 - delta}
        r = dist_eval(mu, 2, c_list)
        print(f"{delta:7.2f} {max(r['marg']):9.4f} {r['HA']:8.4f} "
              f"{r['HU']:8.4f} {r['D']:8.4f} {r['f'][1.0]:+8.4f} "
              f"{r['cstar']:7.4f}")
    print("Reading: Ellis's example refutes c <= c* ~ 1.15 but only at")
    print("marginals ~ 1/2; by delta = 0.05 it no longer refutes c = 1.")
    print("Small-n examples of this type leave the small-c rungs untouched --")
    print("that gap is closed by Sawin's construction (part C).")


# ------------------------- C. Sawin Proposition 6, exact finite-n evaluation

def sawin_eval(n, ubar, theta, c_list, tail=1e-28):
    """Exact evaluation of the (truncated-geometric) Sawin Prop 6 distribution.

    A: draw k in {0..K} w.p. Pk ~ (1-theta)theta^k (renormalized), then
       A | k ~ Bernoulli(p_k)^n with p_k = 1-(1-ubar)^(k+1).
    U = A u B: mixture over k' = kA+kB+1 of Bernoulli(p_k')^n  (exact, since
       (1-p_a)(1-p_b) = (1-ubar)^(a+b+2) = 1 - p_{a+b+1}).
    Both laws are exchangeable: atom prob depends on weight w only, so
    H(A) = -sum_w C(n,w) a_w log a_w   etc., computed in log space.
    All reported quantities are exact (float rounding only) for the truncated
    model, which is itself a legitimate probability distribution."""
    K = 0
    while theta ** (K + 1) > tail:
        K += 1
    raw = [(1 - theta) * theta ** k for k in range(K + 1)]
    Z = sum(raw)
    Pk = [r / Z for r in raw]
    l2u = log2(1.0 - ubar)                      # < 0
    pk = [1.0 - (1.0 - ubar) ** (k + 1) for k in range(K + 1)]
    margA = sum(Pk[k] * pk[k] for k in range(K + 1))

    P2 = [0.0] * (2 * K + 2)                    # index = k' in 1..2K+1
    for a in range(K + 1):
        for b in range(K + 1):
            P2[a + b + 1] += Pk[a] * Pk[b]

    def comps(weights, kprimes):
        """Per component: (base, slope) with log2 atom = base + w*slope."""
        out = []
        for wgt, kp in zip(weights, kprimes):
            if wgt <= 0.0:
                continue
            l1mp = (kp + 1) * l2u               # log2(1-p_{kp}), exact
            lp = log2(1.0 - (1.0 - ubar) ** (kp + 1))
            out.append((log2(wgt) + n * l1mp, lp - l1mp))
        return out

    cA = comps(Pk, range(K + 1))
    cU = comps(P2[1:], range(1, 2 * K + 2))

    HA = HU = CE = sumA = sumU = 0.0
    for w in range(n + 1):
        lc = log2C(n, w)
        la = logsumexp2([b + w * s for b, s in cA])
        lu = logsumexp2([b + w * s for b, s in cU])
        mA = 2.0 ** (lc + la)
        mU = 2.0 ** (lc + lu)
        sumA += mA
        sumU += mU
        HA -= mA * la
        HU -= mU * lu
        CE -= mU * la
    assert abs(sumA - 1.0) < 1e-6 and abs(sumU - 1.0) < 1e-6, (sumA, sumU)
    D = CE - HU
    cstar = (HA - HU) / D
    return dict(n=n, ubar=ubar, theta=theta, K=K, marg=margA, HA=HA, HU=HU,
                D=D, cstar=cstar, f={c: HU + c * D - HA for c in c_list},
                ideal=h(2 * ubar - ubar ** 2) - h(ubar))


def part_C(c_list, big):
    print("\n" + "=" * 78)
    print("C. Sawin's Proposition 6 distribution, EXACT finite-n evaluation")
    print("=" * 78)
    print("mu: k ~ Geom(theta) (truncated), A|k ~ Bern(1-(1-ubar)^(k+1))^n.")
    print("Sawin (arXiv:2211.11504): H(U) <= d*H(A) with d ~ h(2u-u^2)/h(u)")
    print("< 1 for ubar > psi, while D(U||A) = O(1) as n -> inf.")
    print("Hence f(c) = H(U)+cD-H(A) < 0 for every fixed c: c*(mu) -> inf.\n")
    cases = [(2000, 0.390, 0.05),
             (20000, 0.386, 0.02),
             ((150000 if big else 60000), 0.3823, 0.001)]
    for n, ubar, theta in cases:
        r = sawin_eval(n, ubar, theta, c_list)
        print(f"n={n}, ubar={ubar}, theta={theta} (K={r['K']}):")
        print(f"  max marginal    = {r['marg']:.6f}"
              f"{'  < 0.38271 = record!' if r['marg'] < RECORD else ''}")
        print(f"  H(A) = {r['HA']:.2f} bits   H(U) = {r['HU']:.2f} bits   "
              f"D(U||A) = {r['D']:.3f} bits")
        print(f"  (H(U)-H(A))/n = {(r['HU'] - r['HA']) / n:+.6f}   "
              f"[theta->0 ideal: {r['ideal']:+.6f}]")
        fstr = "   ".join(f"f({c:g})={v:+.2f}" for c, v in r["f"].items())
        print(f"  {fstr}")
        print(f"  c*(mu) = {r['cstar']:.1f}   ==> refutes (S_c) for ALL "
              f"c <= {r['cstar']:.1f} at thresholds p > {r['marg']:.5f}\n")
    print("Reading: D stays ~ log2(1/theta)+O(1) bits while H(A)-H(U) grows")
    print("linearly in n, so c* is unbounded.  With ubar -> psi+ the marginal")
    print("comes arbitrarily close to psi: every c is refuted at every p>psi.")


# ------------------- D. cross-check: smoothed-slice construction (certified)

def slice_cert(n, p0, eps):
    """Independent construction: mu = (1-t) Unif(slice w0=round(p0 n))
    + t Bern(q)^n, t = 2^(-eps n), q = 1-(1-m)^2 with m the marginal
    (fixed point).  CERTIFIED bounds, valid for all c in [0,1]:
      H(U)  <= n h(q)                    (subadditivity; Pr[i in U] = q)
      CE(U,mu) <= log2(1/t) + n h(q)     (mu >= t Bern(q)^n pointwise;
                                          cross-entropy vs a product measure
                                          depends only on marginals of U)
      H(mu) >= (1-t) log2 C(n,w0) + t n h(q)     (concavity)
    hence f(c) <= F(c) := n h(q) + c*eps*n - Hlb, linear increasing in c:
    F(1) < 0 certifies violation of (S_c) for every c in [0,1]."""
    w0 = round(p0 * n)
    p0 = w0 / n
    t = 2.0 ** (-eps * n)
    m = p0
    for _ in range(200):
        m = (1 - t) * p0 + t * (1.0 - (1.0 - m) ** 2)
    q = 1.0 - (1.0 - m) ** 2
    Hlb = (1 - t) * log2C(n, w0) + t * n * h(q)
    F0 = n * h(q) - Hlb
    F1 = n * h(q) + eps * n - Hlb
    return m, F0, F1


def part_D(big):
    print("\n" + "=" * 78)
    print("D. Cross-check: slice + 2^(-eps n) Bernoulli smoothing "
          "(certified bounds)")
    print("=" * 78)
    print("A 'regularized Chase-Lovett' distribution: full support (so D is")
    print("finite), union mass caught by a planted Bern(q) component of mass")
    print("t = 2^(-eps n), costing only eps*n bits of KL vs a Theta(n)")
    print("entropy drop.  F(1) < 0 certifies f(c) < 0 for all c in [0,1].\n")
    print(f"{'n':>9} {'p0':>7} {'eps':>7} {'max marg':>9} "
          f"{'F(0)':>10} {'F(1)':>10}")
    cases = [(200000, 0.3927, 0.004),
             (2000000, 0.3835, 0.0004),
             ((40000000 if big else 20000000), 0.3822, 0.00005)]
    for n, p0, eps in cases:
        m, F0, F1 = slice_cert(n, p0, eps)
        note = "  (marg < record)" if m < RECORD else ""
        print(f"{n:>9} {p0:7.4f} {eps:7.5f} {m:9.5f} {F0:+10.1f} "
              f"{F1:+10.1f}{note}")
    print("\nAll F(1) < 0: independent confirmation of part C.")


# ----------------------------------------------------------------- summary

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--big", action="store_true",
                    help="larger n in parts C/D (slower)")
    args = ap.parse_args()

    c_needed = part_A()
    c_list = [round(c_needed, 4), 0.1, 0.5, 1.0, 2.0]
    part_B(c_list)
    part_C(c_list, args.big)
    part_D(args.big)

    print("\n" + "=" * 78)
    print("E. VERDICT")
    print("=" * 78)
    print(f"psi = {PSI:.6f}.  For every c >= 0:")
    print("  * (S_c) HOLDS at threshold psi    (implied by the c=0 theorem of")
    print("    AHS/Chase-Lovett/Sawin, since f_mu(c) >= f_mu(0));")
    print("  * (S_c) FAILS at every threshold p > psi  (Sawin Prop 6 with")
    print("    ubar in (psi, p): H(A)-H(U) = Theta(n), D = O(1), c* -> inf).")
    print("The weighted-KL ladder is exactly flat: P(c) = psi for all c.")
    print("No weight c beats psi = 0.381966, let alone the 0.38271 record.")


if __name__ == "__main__":
    main()
