#!/usr/bin/env python3
"""uc_mitax_skeptic.py -- adversarial re-verification of 009 (Gap 2 probe,
attempts/009-mutual-information-tax.md), for the skeptic record 011.

Everything here is re-implemented FROM SCRATCH -- no imports from
uc_mitax.py, no shared code paths.  Deliberate implementation differences:

  * natural-log internals everywhere (uc_mitax.py works in log2);
  * exchangeable couplings evaluated by a 4-cell RECURSION over one future
    coordinate (T(m,wa,wb) = sum_{ab} rho^{ab} T(m-1,wa+a,wb+b)) plus
    forward filtering of the prefix state, instead of 009's closed-form
    binomial sums S(m,ra,rb) and prefix*future mass products;
  * asymmetric alternating Sinkhorn (row scale, then column scale;
    symmetry of the result CHECKED, not imposed) instead of 009's damped
    symmetric iteration;
  * the taxes are cross-checked by DIRECT conditional-mutual-information
    computation (T_A = sum_i I(A_i; B_{<i} | A_{<i}) from the joint law)
    and the U-side chain rule (H(U) = sum_i H(U_i|U_{<i})), rather than
    only through 009's entropy-difference identities.

Parts:
  A. evaluator validation: product+Plackett closed forms; pure-slice tilt
     n=8 vs the checked-in 003_partB2.json gain; identity cross-checks
     (chain rule for U; direct-MI tax vs identity tax); exchangeable
     recursion engine vs the atom engine on a slice AND on a
     Sinkhorn-fitted product mixture.
  B. crash family (005 Prop 3): reproduce 009 part B, verify the exact
     crash-OR identity log2 OR = lam(3-n) from my own fitted couplings,
     extend the lambda grid past both edges, extend n to 14/17/20, and
     bisect the positive-CR window edge lam_max(n) to test the scaling
     lam_max ~ c/(n-3) suggested by the crash exponent.
  C. half-mixing genre (004 R1): own construction; ST = 0 (mini-lemma HM)
     at machine precision plus a DIRECT sigma-field check (conditional
     P(U_i=1|a,b) constant across (a,b) refining the same u); gain vs
     004's closed form; two parameter points 009 never ran.
  D. Sawin certificates: own profile arithmetic (own truncation) for
     H(mu), max marginal, H(P), the block-coupling bound CR_lb =
     n sum_k P_k h(m_k) - H(mu); legality check p_k >= 1-1/sqrt(2)
     wherever m_k = 1/2 is used; independent H(U_block) to verify
     0 <= gain_block - CR_lb <= H(P) and 003_partD's gains.
  E. Chase-Lovett pure slice: reproduce all six 009 part-E sweeps with
     the recursion engine; extend lambda past the grid edge (3.5..6);
     extend p past 0.425 (w0 = 35, 36 at n = 80); ST-scaling campaign at
     fixed p with n up to 240 and lambda in {1.5, 3}; least-squares fits
     ST = a + b log2 n vs ST = c n^alpha.
  F. THE UNTESTED CELL: tilt-recipe CR on Sawin gadgets at large n
     (profile Sinkhorn + recursion engine), instances of 003 part (a)
     ((200, .40, .05) and (300, .3823, .02)) and n-sweeps of both shapes;
     block-recipe CR_lb at the same instances for comparison.

Usage:  python3 problems/union-closed/explore/uc_mitax_skeptic.py [PARTS]
        PARTS defaults to ABCDEF.  Deterministic, stdlib only.
Checkpoints: problems/union-closed/data/mitaxsk_part[A-F].json
"""

import json
import os
import sys
import time
from array import array
from math import exp, inf, lgamma, log, sqrt

LN2 = log(2.0)
NEG = -inf
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def hb(p):
    """binary entropy in BITS"""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return (-p * log(p) - (1.0 - p) * log(1.0 - p)) / LN2


def lnC(n, k):
    if k < 0 or k > n:
        return NEG
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def lse(xs):
    m = max(xs)
    if m == NEG:
        return NEG
    return m + log(sum(exp(x - m) for x in xs))


def save(name, obj):
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, name)
    if os.path.exists(path) and not name.startswith("mitaxsk_"):
        raise RuntimeError(f"refusing to overwrite {path}")
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    print(f"  [checkpoint: {os.path.relpath(path)}]")


# =============================================================== atom engine
def eval_atoms(n, pairs):
    """Full-history Gap-2 accounting for an explicit coupling, from scratch.

    pairs: dict (Amask, Bmask) -> prob.  Returns everything in BITS, plus
    internal identity checks:
      - uchain_err: |H(U) - sum_i H(U_i | U_{<i})|  (chain rule, U side)
      - taxA_mi_err: |tax_A(identity) - sum_i I(A_i; B_{<i}|A_{<i})| where
        the MI is computed directly from the joint conditional tables.
    """
    ma, mb, ulaw = {}, {}, {}
    tot = 0.0
    for (A, B), w in pairs.items():
        ma[A] = ma.get(A, 0.0) + w
        mb[B] = mb.get(B, 0.0) + w
        ulaw[A | B] = ulaw.get(A | B, 0.0) + w
        tot += w
    marg_dev = max(abs(ma.get(k, 0.0) - mb.get(k, 0.0))
                   for k in set(ma) | set(mb))
    Hmu = -sum(w * log(w) for w in ma.values() if w > 0.0) / LN2
    HU = -sum(w * log(w) for w in ulaw.values() if w > 0.0) / LN2
    sum_hz = sum_hx = sum_hy = 0.0
    sum_hu_uhist = 0.0          # sum_i H(U_i | U_{<i})
    taxA_mi = 0.0               # sum_i I(A_i; B_{<i} | A_{<i}) directly
    for i in range(n):
        mask = (1 << i) - 1
        gj = {}                 # (a,b)-history cells
        gu = {}                 # u-history cells
        ga = {}                 # a-history marginal cells: [mass, A_i=0 mass]
        gab = {}                # (a, b) -> [mass, A_i=0 mass]  (for MI)
        for (A, B), w in pairs.items():
            a, b = A & mask, B & mask
            ai, bi = (A >> i) & 1, (B >> i) & 1
            cell = gj.setdefault((a, b), [0.0, 0.0, 0.0, 0.0])
            cell[0] += w
            if ai == 0 and bi == 0:
                cell[1] += w
            if ai == 0:
                cell[2] += w
            if bi == 0:
                cell[3] += w
            u = a | b
            uc = gu.setdefault(u, [0.0, 0.0])
            uc[0] += w
            if (ai | bi) == 0:
                uc[1] += w
            ac = ga.setdefault(a, [0.0, 0.0])
            ac[0] += w
            if ai == 0:
                ac[1] += w
        for m_, z_, x_, y_ in gj.values():
            if m_ > 0.0:
                sum_hz += m_ * hb(z_ / m_)
                sum_hx += m_ * hb(x_ / m_)
                sum_hy += m_ * hb(y_ / m_)
        for m_, z_ in gu.values():
            if m_ > 0.0:
                sum_hu_uhist += m_ * hb(z_ / m_)
        # direct conditional MI: I(A_i; B_{<i} | A_{<i})
        #   = sum_a P(a) [ h(x(a)) - sum_b P(b|a) h(x(a,b)) ]
        mi = 0.0
        for a, (m_, x_) in ga.items():
            if m_ > 0.0:
                mi += m_ * hb(x_ / m_)
        for (a, b), cell in gj.items():
            if cell[0] > 0.0:
                mi -= cell[0] * hb(cell[2] / cell[0])
        taxA_mi += mi
    CR = sum_hz - Hmu
    gain = HU - Hmu
    taxA = Hmu - sum_hx
    # fact (iii) check material: I(A;B) = H(A) + H(B) - H(A,B)
    HB_ = -sum(w * log(w) for w in mb.values() if w > 0.0) / LN2
    Hj = -sum(w * log(w) for w in pairs.values() if w > 0.0) / LN2
    return {"n": n, "H_mu": Hmu, "H_U": HU, "gain": gain, "CR": CR,
            "tax_A": taxA, "tax_B": Hmu - sum_hy,
            "second_tax": gain - CR, "mass": tot, "marg_dev": marg_dev,
            "I_AB": Hmu + HB_ - Hj,
            "uchain_err": abs(HU - sum_hu_uhist),
            "taxA_mi_err": abs(taxA - taxA_mi)}


def sinkhorn_atoms(atoms, lnmu, lam2, tol=1e-13, itmax=60000):
    """Alternating (asymmetric) scaling of K(A,B)=2^{lam2|A&B|} to margins
    exp(lnmu) on both sides.  Natural-log domain.  Returns (pairs, resid,
    sym_dev): resid = max abs ln-error of BOTH margins, sym_dev = max
    |r - c - const| deviation of the two scaling vectors."""
    L = lam2 * LN2
    N = len(atoms)
    K = [[L * bin(a & b).count("1") for b in atoms] for a in atoms]
    r = list(lnmu)
    c = [0.0] * N
    resid = inf
    for it in range(itmax):
        r = [lnmu[i] - lse([c[j] + K[i][j] for j in range(N)])
             for i in range(N)]
        c = [lnmu[j] - lse([r[i] + K[i][j] for i in range(N)])
             for j in range(N)]
        # row-margin residual (col margins exact after the c-update)
        resid = max(abs(lse([c[j] + K[i][j] for j in range(N)])
                        + r[i] - lnmu[i]) for i in range(N))
        if resid < tol:
            break
    d = [r[i] - c[i] for i in range(N)]
    sym_dev = max(d) - min(d)
    pairs = {}
    for i, a in enumerate(atoms):
        for j, b in enumerate(atoms):
            pairs[(a, b)] = exp(r[i] + c[j] + K[i][j])
    s = sum(pairs.values())
    for k in pairs:
        pairs[k] /= s
    return pairs, resid, sym_dev


def plackett_z(x, y, rho):
    """P(both zero) for the 2x2 coupling with zero-margins x, y and odds
    ratio rho.  Own derivation: z(1-x-y+z) = rho (x-z)(y-z)  =>
    (rho-1) z^2 - (1+(rho-1)(x+y)) z + rho x y = 0, smaller root."""
    if x <= 0.0 or y <= 0.0:
        return 0.0
    if x >= 1.0:
        return y
    if y >= 1.0:
        return x
    if abs(rho - 1.0) < 1e-14:
        return x * y
    a = rho - 1.0
    S = 1.0 + a * (x + y)
    z = (S - sqrt(S * S - 4.0 * a * rho * x * y)) / (2.0 * a)
    return min(max(z, max(0.0, x + y - 1.0)), min(x, y))


# ==================================================== exchangeable recursion
def exch_cr(n, lam2, logf, prune=1e-16):
    """Gap-2 accounting for pi(A,B) ∝ f(|A|) f(|B|) 2^(lam2 |A∩B|).

    Method (independent of 009's): T(m, wa, wb) = total weight of futures
    of length m given prefix weights (wa, wb), by the one-coordinate
    recursion T(m,wa,wb) = sum_{al,be} rho^{al be} T(m-1, wa+al, wb+be),
    T(0,wa,wb) = f(wa) f(wb).  Prefix-state law by forward filtering with
    the conditional 4-cell laws  q(al,be|l,wa,wb) ∝ rho^{al be}
    T(n-l-1, wa+al, wb+be) / T(n-l, wa, wb).  H(U) from the exchangeable
    U-profile.  All outputs in bits."""
    L = lam2 * LN2
    # ---- T tables, flattened per level
    T = [None] * (n + 1)
    lvl0 = array("d", [0.0]) * 0
    lvl0 = array("d", [logf[wa] + logf[wb]
                       for wa in range(n + 1) for wb in range(n + 1)])
    T[0] = lvl0
    for m in range(1, n + 1):
        sz = n - m + 1          # wa, wb in 0..n-m  -> sz values
        psz = sz + 1            # previous level row length
        prev = T[m - 1]
        cur = array("d", bytes(8 * sz * sz))
        idx = 0
        for wa in range(sz):
            base = wa * psz
            for wb in range(sz):
                t00 = prev[base + wb]
                t01 = prev[base + wb + 1]
                t10 = prev[base + psz + wb]
                t11 = prev[base + psz + wb + 1] + L
                mx = t00
                if t01 > mx:
                    mx = t01
                if t10 > mx:
                    mx = t10
                if t11 > mx:
                    mx = t11
                if mx == NEG:
                    cur[idx] = NEG
                else:
                    cur[idx] = mx + log(exp(t00 - mx) + exp(t01 - mx)
                                        + exp(t10 - mx) + exp(t11 - mx))
                idx += 1
        T[m] = cur
    lnZ = T[n][0]
    # ---- marginal profile of the coupling: P(A) ∝ f(wa) * M(wa),
    #      M(wa) = sum_B f(|B|) rho^{|A∩B|}; from T with the identity
    #      M(wa) = T-like sum: use W-recursion instead (independent):
    #      W(j, m) = ln sum_u C(m,u) f(j+u); V(w)=lse_j[lnC(w,j)+L j+W(j,n-w)]
    W = [[NEG] * (n + 1) for _ in range(n + 1)]   # W[j][m]
    for j in range(n + 1):
        W[j][0] = logf[j]
    for m in range(1, n + 1):
        for j in range(n + 1 - m):
            W[j][m] = lse([W[j][m - 1], W[j + 1][m - 1]])
    lnmarg = []                  # per-atom ln prob of the A-marginal
    for w in range(n + 1):
        V = lse([lnC(w, j) + L * j + W[j][n - w] for j in range(w + 1)])
        lnmarg.append(logf[w] + V - lnZ)
    mass_marg = sum(exp(lnC(n, w) + lm) for w, lm in enumerate(lnmarg)
                    if lm > NEG)
    Hmu = -sum(exp(lnC(n, w) + lm) * lm for w, lm in enumerate(lnmarg)
               if lm > NEG) / LN2
    max_freq = 0.0
    # element frequency: P(elem in A) = sum_w (w/n) profile_mass(w)
    max_freq = sum(exp(lnC(n, w) + lm) * w / n
                   for w, lm in enumerate(lnmarg) if lm > NEG)
    # ---- H(U) via the exchangeable U-profile
    HU = 0.0
    mass_u = 0.0
    for wu in range(n + 1):
        terms = []
        for wa in range(wu + 1):
            for j in range(wa + 1):
                # |A|=wa, |A∩B|=j, |B| = j + wu - wa; partition of u into
                # A∩B (j), A-only (wa-j), B-only (wu-wa)
                lmult = (lgamma(wu + 1) - lgamma(j + 1)
                         - lgamma(wa - j + 1) - lgamma(wu - wa + 1))
                t = lmult + logf[wa] + logf[j + wu - wa] + L * j
                if t > NEG:
                    terms.append(t)
        if not terms:
            continue
        lq = lse(terms) - lnZ
        m_ = exp(lnC(n, wu) + lq)
        mass_u += m_
        HU -= m_ * lq / LN2
    # ---- forward filtering for the chain-rule sums
    P = {(0, 0): 1.0}
    sum_hz = sum_hx = sum_hy = 0.0
    mass_min = inf
    lost = 0.0
    for l in range(n):
        m1 = n - l - 1
        sz1 = n - m1            # row length of level m1 is (n-m1+1)
        row1 = m1 and 0
        Tm1 = T[m1]
        rl1 = n - m1 + 1        # row length in level m1
        newP = {}
        mass_here = sum(P.values())
        mass_min = min(mass_min, mass_here)
        for (wa, wb), pw in P.items():
            if pw < prune:
                lost += pw
                continue
            t00 = Tm1[wa * rl1 + wb]
            t01 = Tm1[wa * rl1 + wb + 1] if wb + 1 < rl1 else NEG
            t10 = Tm1[(wa + 1) * rl1 + wb] if wa + 1 < rl1 else NEG
            t11 = (Tm1[(wa + 1) * rl1 + wb + 1] + L
                   if wa + 1 < rl1 and wb + 1 < rl1 else NEG)
            mx = max(t00, t01, t10, t11)
            if mx == NEG:
                continue
            e00, e01 = exp(t00 - mx), exp(t01 - mx)
            e10, e11 = exp(t10 - mx), exp(t11 - mx)
            s = e00 + e01 + e10 + e11
            z = e00 / s
            x = (e00 + e01) / s
            y = (e00 + e10) / s
            sum_hz += pw * hb(z)
            sum_hx += pw * hb(x)
            sum_hy += pw * hb(y)
            for (al, be), ew in (((0, 0), e00), ((0, 1), e01),
                                 ((1, 0), e10), ((1, 1), e11)):
                if ew > 0.0:
                    k = (wa + al, wb + be)
                    newP[k] = newP.get(k, 0.0) + pw * ew / s
        P = newP
    CR = sum_hz - Hmu
    gain = HU - Hmu
    return {"n": n, "lam2": lam2, "H_mu": Hmu, "H_U": HU, "gain": gain,
            "CR": CR, "tax_A": Hmu - sum_hx, "tax_B": Hmu - sum_hy,
            "second_tax": gain - CR, "mass_marg": mass_marg,
            "mass_u": mass_u, "mass_min": mass_min, "lost": lost,
            "max_freq": max_freq}


def sinkhorn_profile(n, lam2, lnmu_atom, tol=1e-12, itmax=20000):
    """Fit exchangeable ln-potential g(w) so that the tilt coupling's
    A-marginal has per-atom ln-probability lnmu_atom(w).  Fixed point
    g(w) = lnmu(w) - ln V_g(w), V_g(w) = sum_j C(w,j) rho^j W_g(j, n-w),
    W_g(j,m) = sum_u C(m,u) e^{g(j+u)}.  Damped half-steps.  Returns
    (g, resid) with resid = max_w |ln-marginal - lnmu|."""
    L = lam2 * LN2
    g = [x * 0.5 for x in lnmu_atom]
    resid = inf
    for it in range(itmax):
        W = [[NEG] * (n + 1) for _ in range(n + 1)]
        for j in range(n + 1):
            W[j][0] = g[j]
        for m in range(1, n + 1):
            for j in range(n + 1 - m):
                W[j][m] = lse([W[j][m - 1], W[j + 1][m - 1]])
        lnV = [lse([lnC(w, j) + L * j + W[j][n - w] for j in range(w + 1)])
               for w in range(n + 1)]
        gnew = [lnmu_atom[w] - lnV[w] for w in range(n + 1)]
        # residual uses the CURRENT g's marginal: ln-marg(w) = g + lnV - lnZ
        # normalize via total mass
        lnmass = lse([lnC(n, w) + g[w] + lnV[w] for w in range(n + 1)])
        resid = max(abs(g[w] + lnV[w] - lnmass - lnmu_atom[w])
                    for w in range(n + 1))
        if resid < tol:
            break
        g = [(a + b) * 0.5 for a, b in zip(g, gnew)]
    # center g so that pairs are O(1): fold -lnZ/2 into g
    lnV_half = lnmass / 2.0
    g = [x - lnV_half for x in g]
    return g, resid


# ------------------------------------------------------------------- part A
def part_A():
    print("=" * 78)
    print("SK-A. Evaluator validation + identity cross-checks")
    print("=" * 78)
    out = {}
    # A1: product + Plackett closed forms
    n, p = 8, 0.40
    x = 1.0 - p
    rows = []
    for rho in (1.0, 1.5, 3.0):
        z = plackett_z(x, x, rho)
        # verify the cross-ratio of the 2x2 table equals rho
        w11 = 1.0 - 2.0 * x + z
        cr_ratio = z * w11 / ((x - z) ** 2) if rho != 1.0 else 1.0
        cell = {(0, 0): z, (0, 1): x - z, (1, 0): x - z, (1, 1): w11}
        pairs = {}
        for A in range(1 << n):
            for B in range(1 << n):
                w = 1.0
                for i in range(n):
                    w *= cell[((A >> i) & 1, (B >> i) & 1)]
                if w > 0.0:
                    pairs[(A, B)] = w
        r = eval_atoms(n, pairs)
        closed = n * (hb(z) - hb(p))
        ok = (abs(r["CR"] - closed) < 1e-9 and abs(r["tax_A"]) < 1e-9
              and abs(r["second_tax"]) < 1e-9
              and r["uchain_err"] < 1e-9 and r["taxA_mi_err"] < 1e-9)
        rows.append({"rho": rho, "CR": r["CR"], "closed": closed,
                     "tax_A": r["tax_A"], "ST": r["second_tax"],
                     "or_err": abs(cr_ratio - rho),
                     "uchain_err": r["uchain_err"],
                     "taxA_mi_err": r["taxA_mi_err"], "ok": ok})
        print(f"  A1 product Bern({p})^{n} rho={rho:4g}: CR {r['CR']:+.9f} "
              f"closed {closed:+.9f}  |OR-rho| {abs(cr_ratio-rho):.1e} "
              f"uchain {r['uchain_err']:.1e} miTA {r['taxA_mi_err']:.1e} "
              f"{'OK' if ok else 'FAIL'}")
    out["A1_product_plackett"] = rows

    # A2: pure slice tilt n=8 w0=3 lam=0.7, atom engine
    n, w0, lam2 = 8, 3, 0.7
    rho = 2.0 ** lam2
    atoms = [m for m in range(1 << n) if bin(m).count("1") == w0]
    pairs = {}
    for a in atoms:
        for b in atoms:
            pairs[(a, b)] = rho ** bin(a & b).count("1")
    s = sum(pairs.values())
    for k in pairs:
        pairs[k] /= s
    r = eval_atoms(n, pairs)
    with open(os.path.join(DATA, "003_partB2.json")) as f:
        ref_gain = json.load(f)["gain_direct"]
    with open(os.path.join(DATA, "mitax_partA.json")) as f:
        ref_a = json.load(f)["slice_n8"]
    out["A2_slice_atom"] = {
        "gain": r["gain"], "gain_ref_003": ref_gain,
        "dgain_003": abs(r["gain"] - ref_gain),
        "CR": r["CR"], "dCR_009": abs(r["CR"] - ref_a["CR"]),
        "tax_A": r["tax_A"], "ST": r["second_tax"],
        "dST_009": abs(r["second_tax"] - ref_a["second_tax"]),
        "uchain_err": r["uchain_err"], "taxA_mi_err": r["taxA_mi_err"]}
    print(f"  A2 slice n=8 w0=3 lam=0.7 (atoms): gain {r['gain']:+.9f} "
          f"(003: {ref_gain:+.9f}, d={abs(r['gain']-ref_gain):.1e}); "
          f"CR d={abs(r['CR']-ref_a['CR']):.1e} vs 009; "
          f"tax_A {r['tax_A']:.2e}")

    # A3: recursion engine vs atom engine, pure slice
    logf = [NEG] * (n + 1)
    logf[w0] = 0.0
    re_ = exch_cr(n, lam2, logf)
    d = {k: abs(re_[k] - r[k]) for k in ("CR", "gain", "tax_A",
                                         "second_tax", "H_mu", "H_U")}
    out["A3_recursion_vs_atoms_slice"] = d
    print(f"  A3 recursion vs atoms (slice): "
          + " ".join(f"d{k}={v:.1e}" for k, v in d.items()))

    # A4: recursion+profile-Sinkhorn vs atoms+atom-Sinkhorn on a mixture
    n = 8
    comp = [(0.6, 0.35), (0.4, 0.15)]
    lnmu = []
    for A in range(1 << n):
        w = bin(A).count("1")
        lnmu.append(log(sum(P * q ** w * (1 - q) ** (n - w)
                            for P, q in comp)))
    lam2 = 0.8
    pairs, resid, sym_dev = sinkhorn_atoms(list(range(1 << n)), lnmu, lam2)
    ra = eval_atoms(n, pairs)
    lnmu_atom = [log(sum(P * q ** w * (1 - q) ** (n - w) for P, q in comp))
                 for w in range(n + 1)]
    g, presid = sinkhorn_profile(n, lam2, lnmu_atom)
    rp = exch_cr(n, lam2, g)
    d = {k: abs(rp[k] - ra[k]) for k in ("CR", "gain", "tax_A",
                                         "second_tax", "H_mu", "H_U")}
    out["A4_mixture"] = {"atom_resid": resid, "sym_dev": sym_dev,
                         "prof_resid": presid, **{"d" + k: v
                                                  for k, v in d.items()},
                         "CR": ra["CR"], "tax_A": ra["tax_A"],
                         "ST": ra["second_tax"], "I_AB": ra["I_AB"],
                         "fact_iii_slack": ra["I_AB"] - ra["tax_A"]
                         - ra["tax_B"]}
    print(f"     fact (iii): tax_A = {ra['tax_A']:.6f} <= I(A;B) = "
          f"{ra['I_AB']:.6f} "
          f"({'OK' if ra['tax_A'] <= ra['I_AB'] + 1e-12 else 'FAIL'})")
    print(f"  A4 mixture n=8 lam=0.8: atom-resid {resid:.1e} "
          f"sym {sym_dev:.1e} prof-resid {presid:.1e}; "
          + " ".join(f"d{k}={v:.1e}" for k, v in d.items()))
    print(f"     (CR {ra['CR']:+.6f} tax_A {ra['tax_A']:.6f} "
          f"ST {ra['second_tax']:.6f}; identities: uchain "
          f"{ra['uchain_err']:.1e}, miTA {ra['taxA_mi_err']:.1e})")
    save("mitaxsk_partA.json", out)
    return out


# ------------------------------------------------------------------- part B
def crash_atoms(n):
    S1 = 1 << 1
    S2 = ((1 << n) - 1) ^ 0b11
    S3 = (1 << n) - 1
    S4 = 1 << 0
    return [S1, S2, S3, S4], [log(0.33), log(0.33), log(0.04), log(0.30)]


def crash_cr(n, lam2):
    atoms, lnmu = crash_atoms(n)
    pairs, resid, sym_dev = sinkhorn_atoms(atoms, lnmu, lam2)
    r = eval_atoms(n, pairs)
    # crash OR at coordinate 2 (bit 1), history (A_1=0, B_1=1) (bit 0)
    cells = {}
    for (A, B), w in pairs.items():
        if (A & 1) == 0 and (B & 1) == 1:
            key = ((A >> 1) & 1, (B >> 1) & 1)
            cells[key] = cells.get(key, 0.0) + w
    or_log = None
    if all(cells.get(c, 0.0) > 0.0 for c in ((0, 0), (0, 1), (1, 0), (1, 1))):
        or_log = (log(cells[(1, 1)]) + log(cells[(0, 0)])
                  - log(cells[(1, 0)]) - log(cells[(0, 1)])) / LN2
    r["resid"] = resid
    r["sym_dev"] = sym_dev
    r["log2_OR_crash"] = or_log
    return r


def crash_iid_closed():
    """CR at lam=0 by the 2-coordinate reduction, closed form (no engine)."""
    m = [0.33, 0.33, 0.04, 0.30]     # S1={2},S2={3..n},S3=[n],S4={1}
    # coordinate 1 (element 1): in S3, S4
    pA1 = m[2] + m[3]
    z1 = (1 - pA1) ** 2
    # coordinate 2 (element 2): in S1, S3
    # histories (A1, B1); conditional zero-prob of coord 2:
    #  A1=0 -> atoms {S1,S2}: P(A2=0|A1=0) = m[1]/(m[0]+m[1])
    #  A1=1 -> atoms {S3,S4}: P(A2=0|A1=1) = m[3]/(m[2]+m[3])
    x0 = m[1] / (m[0] + m[1])
    x1 = m[3] / (m[2] + m[3])
    e2 = 0.0
    for a1, xa in ((0, x0), (1, x1)):
        for b1, xb in ((0, x0), (1, x1)):
            pa = (m[0] + m[1]) if a1 == 0 else (m[2] + m[3])
            pb = (m[0] + m[1]) if b1 == 0 else (m[2] + m[3])
            e2 += pa * pb * hb(xa * xb)
    Hmu = -sum(x * log(x) for x in m) / LN2
    return hb(z1) + e2 - Hmu


def part_B():
    print("\n" + "=" * 78)
    print("SK-B. Crash family: reproduce, extend past grid edges, window law")
    print("=" * 78)
    out = {"iid_CR_closed": crash_iid_closed()}
    print(f"  closed-form iid CR (n-independent) = "
          f"{out['iid_CR_closed']:+.6f}  (009: +0.1114)")
    ref = {}
    with open(os.path.join(DATA, "mitax_partB.json")) as f:
        for blk in json.load(f):
            for row in blk["rows"]:
                ref[(blk["n"], row["lam2"])] = row
    LAMS = [-2.0, -1.5, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5,
            2.0, 3.0, 3.5, 4.0, 5.0]
    sweeps = []
    maxdiff = 0.0
    for n in (5, 8, 11, 14, 17):
        rows = []
        best = (-inf, None)
        for lam2 in LAMS:
            r = crash_cr(n, lam2)
            row = {"lam2": lam2, "CR": r["CR"], "gain": r["gain"],
                   "tax_A": r["tax_A"], "ST": r["second_tax"],
                   "resid": r["resid"], "log2_OR": r["log2_OR_crash"],
                   "OR_pred_err": (abs(r["log2_OR_crash"] - lam2 * (3 - n))
                                   if r["log2_OR_crash"] is not None
                                   and lam2 != 0.0 else None)}
            rows.append(row)
            if r["CR"] > best[0]:
                best = (r["CR"], lam2)
            key = (n, lam2)
            if key in ref:
                d = max(abs(row["CR"] - ref[key]["CR"]),
                        abs(row["gain"] - ref[key]["gain"]))
                maxdiff = max(maxdiff, d)
        edge_note = (f"CR(-2)={rows[0]['CR']:+.4f} "
                     f"CR(5)={rows[-1]['CR']:+.4f}")
        orerrs = [row["OR_pred_err"] for row in rows
                  if row["OR_pred_err"] is not None]
        print(f"  n={n:2d}: sup CR {best[0]:+.5f} at lam {best[1]}; "
              f"{edge_note}; max|log2OR - lam(3-n)| = {max(orerrs):.1e}")
        sweeps.append({"n": n, "rows": rows, "best_CR": best[0],
                       "best_lam2": best[1], "max_or_err": max(orerrs)})
    out["max_diff_vs_009_partB"] = maxdiff
    print(f"  max |CR/gain - 009 partB| over shared grid = {maxdiff:.2e}")
    # window edges by bisection: lam_max (upper zero of CR), lam_min (lower)
    win = []
    for n in (5, 8, 11, 14, 17, 20):
        def cr(l):
            return crash_cr(n, l)["CR"]
        # find upper crossing: scan up from 0.05
        lo, hi = None, None
        lam = 0.05
        prev = cr(lam)
        while lam < 6.0:
            nxt = lam + 0.1
            v = cr(nxt)
            if prev > 0.0 and v <= 0.0:
                lo, hi = lam, nxt
                break
            prev, lam = v, nxt
        lmax = None
        if lo is not None:
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if cr(mid) > 0.0:
                    lo = mid
                else:
                    hi = mid
            lmax = 0.5 * (lo + hi)
        # lower crossing (negative lam side)
        lo2, hi2 = None, None
        lam = -0.05
        prev = cr(lam)
        while lam > -3.0:
            nxt = lam - 0.1
            v = cr(nxt)
            if prev > 0.0 and v <= 0.0:
                lo2, hi2 = lam, nxt
                break
            prev, lam = v, nxt
        lmin = None
        if lo2 is not None:
            for _ in range(40):
                mid = 0.5 * (lo2 + hi2)
                if cr(mid) > 0.0:
                    lo2 = mid
                else:
                    hi2 = mid
            lmin = 0.5 * (lo2 + hi2)
        win.append({"n": n, "lam_max": lmax, "lam_min": lmin,
                    "lam_max_times_nm3": (lmax * (n - 3)
                                          if lmax is not None else None)})
        print(f"  window n={n:2d}: lam_max = {lmax:.5f}  lam_min = "
              f"{lmin:.5f}  lam_max*(n-3) = {lmax*(n-3):.4f}")
    out["window"] = win
    out["sweeps"] = sweeps
    save("mitaxsk_partB.json", out)
    return out


# ------------------------------------------------------------------- part C
def half_mixing(n, ph, P1):
    """Own construction of the 004-R1 half-mixing coupling (needs P1<=1/2).
    Cells per drawn s: (0,0),(0,s),(s,0) each P1/2; (s,s) 1-1.5*P1."""
    assert 0.0 < P1 <= 0.5
    pairs = {}
    for s in range(1 << n):
        k = bin(s).count("1")
        ps = ph ** k * (1.0 - ph) ** (n - k)
        for key, w in (((0, 0), P1 / 2), ((0, s), P1 / 2),
                       ((s, 0), P1 / 2), ((s, s), 1.0 - 1.5 * P1)):
            pairs[key] = pairs.get(key, 0.0) + ps * w
    return pairs


def hm_sigma_check(n, pairs, ph):
    """Direct check of mini-lemma HM: over histories (a,b) with u=a|b != 0,
    P(U_i=1|a,b) must equal ph exactly; and for u=0 the only history is
    (0,0).  Returns (max |P-ph| over u!=0, max #distinct (a,b) per u=0)."""
    worst = 0.0
    for i in range(n):
        mask = (1 << i) - 1
        cells = {}
        for (A, B), w in pairs.items():
            a, b = A & mask, B & mask
            c = cells.setdefault((a, b), [0.0, 0.0])
            c[0] += w
            if ((A | B) >> i) & 1:
                c[1] += w
        for (a, b), (m_, u1) in cells.items():
            if (a | b) != 0 and m_ > 1e-300:
                worst = max(worst, abs(u1 / m_ - ph))
    return worst


def part_C():
    print("\n" + "=" * 78)
    print("SK-C. Half-mixing genre: ST=0 identity (HM), gain closed form")
    print("=" * 78)
    out = []
    for ph, P1, tag in ((0.60, 0.50, "004 mu_n"),
                        (0.51, 0.163, "recipe-killer"),
                        (0.75, 0.40, "new point A (009 never ran)"),
                        (0.51, 0.50, "new point B (P1 at its cap)")):
        rows = []
        for n in (6, 10, 14, 16):
            pairs = half_mixing(n, ph, P1)
            r = eval_atoms(n, pairs)
            closed = P1 / 2.0 * n * hb(ph) - (hb(P1) - hb(P1 / 2.0))
            sig = hm_sigma_check(n, pairs, ph)
            rows.append({"n": n, "CR": r["CR"], "gain": r["gain"],
                         "gain_004_closed": closed,
                         "dgain_closed": abs(r["gain"] - closed),
                         "ST": r["second_tax"], "tax_A": r["tax_A"],
                         "sigma_check": sig, "marg_dev": r["marg_dev"],
                         "uchain_err": r["uchain_err"]})
            print(f"  ph={ph} P1={P1} n={n:2d}: ST = {r['second_tax']:.2e}  "
                  f"|gain-closed| = {abs(r['gain']-closed):.2e}  "
                  f"tax_A = {r['tax_A']:.6f}  sigma-check = {sig:.2e}")
        out.append({"ph": ph, "P1": P1, "tag": tag, "rows": rows})
    save("mitaxsk_partC.json", out)
    return out


# ------------------------------------------------------------------- part D
def sawin_profile(n, ubar, theta, tail=1e-30):
    """Own truncation + own lse.  Returns (lnP list, p list) truncated."""
    K = 0
    while theta ** (K + 1) >= tail:
        K += 1
    raw = [(1.0 - theta) * theta ** k for k in range(K + 1)]
    Z = sum(raw)
    P = [x / Z for x in raw]
    p = [1.0 - (1.0 - ubar) ** (k + 1) for k in range(K + 1)]
    return P, p


def sawin_H_profile(n, P, p):
    """H (bits) of  mu = sum_k P_k Bern(p_k)^n  by profile arithmetic."""
    lgp = [log(x) for x in P]
    lpp = [(log(q), log(1.0 - q)) for q in p]
    H = 0.0
    tot = 0.0
    for w in range(n + 1):
        lc = lnC(n, w)
        la = lse([g + w * a + (n - w) * b for g, (a, b) in zip(lgp, lpp)])
        m = exp(lc + la)
        tot += m
        H -= m * la
    return H / LN2, tot


def part_D():
    print("\n" + "=" * 78)
    print("SK-D. Sawin certificates: H(mu), block CR bound, gain sandwich")
    print("=" * 78)
    ref002 = {2000: 1927.70, 20000: 19239.59, 60000: 57578.99}
    with open(os.path.join(DATA, "003_partD.json")) as f:
        g003 = {row["n"]: row.get("gain_coupled") for row in json.load(f)
                if "ubar" in row}
    with open(os.path.join(DATA, "mitax_partD.json")) as f:
        r009 = {row["n"]: row for row in json.load(f)}
    out = []
    PMIN = 1.0 - 1.0 / sqrt(2.0)
    for n, ubar, theta in ((2000, 0.390, 0.05), (20000, 0.386, 0.02),
                           (60000, 0.3823, 0.001)):
        P, p = sawin_profile(n, ubar, theta)
        H, tot = sawin_H_profile(n, P, p)
        marg = sum(Pk * pk for Pk, pk in zip(P, p))
        HP = -sum(Pk * log(Pk) for Pk in P) / LN2
        # legality of m_k = 1/2 for p_k <= 1/2 (needs p_k >= 1-1/sqrt2)
        legal = all(pk >= PMIN for pk in p if pk <= 0.5)
        mk = [0.5 if pk <= 0.5 else pk for pk in p]
        crlb = n * sum(Pk * hb(m) for Pk, m in zip(P, mk)) - H
        # independent H(U_block): law(U) = sum_k P_k Bern(m_k)^n
        HU, totU = sawin_H_profile(n, P, mk)
        gain_blk = HU - H
        row = {"n": n, "ubar": ubar, "theta": theta, "H_mu": H,
               "H_002": ref002[n], "dH_009": abs(H - r009[n]["H_mu"]),
               "marg": marg, "dmarg_009": abs(marg - r009[n]["marg"]),
               "H_P": HP, "CR_lb": crlb,
               "dCRlb_009": abs(crlb - r009[n]["CR_lower_bound"]),
               "mk_legal": legal, "gain_block": gain_blk,
               "gain_003": g003.get(n),
               "sandwich_lo": gain_blk - crlb,       # must be in [0, H_P]
               "mass_dev": abs(tot - 1.0), "massU_dev": abs(totU - 1.0)}
        out.append(row)
        print(f"  n={n:6d}: H(mu) = {H:.4f} (002: {ref002[n]}; "
              f"d009 = {row['dH_009']:.2e}); marg = {marg:.6f}; "
              f"CR_lb = {crlb:+.4f} (d009 = {row['dCRlb_009']:.2e})")
        print(f"       m_k=1/2 legality (p_k >= {PMIN:.4f}): {legal}; "
              f"gain_block = {gain_blk:+.4f} (003: {g003.get(n)}); "
              f"gain-CR_lb = {row['sandwich_lo']:.6f} vs H(P) = {HP:.6f}")
    save("mitaxsk_partD.json", out)
    return out


# ------------------------------------------------------------------- part E
def slice_run(n, w0, lam2):
    logf = [NEG] * (n + 1)
    logf[w0] = 0.0
    return exch_cr(n, lam2, logf)


def lsq_fit_log(ns, ys):
    """least squares y = a + b*log2 n; returns a, b, rms residual."""
    xs = [log(x) / LN2 for x in ns]
    N = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    b = (N * sxy - sx * sy) / (N * sxx - sx * sx)
    a = (sy - b * sx) / N
    rms = sqrt(sum((a + b * x - y) ** 2 for x, y in zip(xs, ys)) / N)
    return a, b, rms


def lsq_fit_pow(ns, ys):
    """least squares ln y = ln c + alpha ln n; returns c, alpha, rms of y."""
    xs = [log(x) for x in ns]
    ls = [log(y) for y in ys]
    N = len(xs)
    sx, sy = sum(xs), sum(ls)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ls))
    al = (N * sxy - sx * sy) / (N * sxx - sx * sx)
    c = exp((sy - al * sx) / N)
    rms = sqrt(sum((c * n ** al - y) ** 2 for n, y in zip(ns, ys)) / N)
    return c, al, rms


def part_E():
    print("\n" + "=" * 78)
    print("SK-E. Pure slice: reproduce 009 part E, extend lambda, p, and n")
    print("=" * 78)
    out = {}
    # E1: reproduce all six sweeps
    with open(os.path.join(DATA, "mitax_partE.json")) as f:
        ref = json.load(f)["sweeps"]
    maxd = 0.0
    for sweep in ref:
        n, w0 = sweep["n"], sweep["w0"]
        for row in sweep["rows"]:
            r = slice_run(n, w0, row["lam2"])
            d = max(abs(r["CR"] - row["CR"]), abs(r["gain"] - row["gain"]),
                    abs(r["second_tax"] - row["second_tax"]))
            maxd = max(maxd, d)
    out["E1_max_diff_vs_009_partE"] = maxd
    print(f"  E1 reproduce all 42 (n,w0,lam) points: max |diff| = {maxd:.2e}")
    # E2: lambda edge  (grid stopped at 3.0)
    rows = []
    for n, w0 in ((24, 10), (80, 34)):
        seq = []
        for lam2 in (3.0, 3.5, 4.0, 5.0, 6.0, 8.0):
            r = slice_run(n, w0, lam2)
            seq.append({"lam2": lam2, "CR": r["CR"], "gain": r["gain"],
                        "ST": r["second_tax"]})
        rows.append({"n": n, "w0": w0, "seq": seq})
        s = " ".join(f"{q['lam2']:g}:{q['CR']:+.4f}" for q in seq)
        print(f"  E2 lambda edge n={n} w0={w0}: CR at lam {s}")
    out["E2_lambda_edge"] = rows
    # E3: p edge at n=80  (grid stopped at w0=34, p=0.425)
    rows = []
    for w0 in (34, 35, 36, 38):
        best = (-inf, None)
        for lam2 in (0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0):
            r = slice_run(80, w0, lam2)
            if r["CR"] > best[0]:
                best = (r["CR"], lam2)
        rows.append({"w0": w0, "p": w0 / 80.0, "best_CR": best[0],
                     "best_lam2": best[1]})
        print(f"  E3 p edge n=80 w0={w0} (p={w0/80:.4f}): sup CR = "
              f"{best[0]:+.4f} at lam {best[1]}")
    out["E3_p_edge"] = rows
    # E4: ST scaling campaign
    campaigns = [
        ("p=5/12, lam=3", [(24, 10), (48, 20), (72, 30), (96, 40),
                           (120, 50), (168, 70), (240, 100)], 3.0),
        ("p=5/12, lam=1.5", [(24, 10), (48, 20), (72, 30), (96, 40),
                             (120, 50)], 1.5),
        ("p=0.425, lam=3", [(40, 17), (80, 34), (120, 51), (160, 68),
                            (200, 85)], 3.0),
        ("p=0.3833, lam=3", [(60, 23), (120, 46), (180, 69), (240, 92)],
         3.0),
    ]
    camps = []
    for tag, pts, lam2 in campaigns:
        ns, sts, crs = [], [], []
        for n, w0 in pts:
            t0 = time.time()
            r = slice_run(n, w0, lam2)
            ns.append(n)
            sts.append(r["second_tax"])
            crs.append(r["CR"])
            print(f"  E4 {tag}: n={n:3d} w0={w0:3d}  ST = "
                  f"{r['second_tax']:.4f}  CR = {r['CR']:+.4f}  "
                  f"tax_A = {r['tax_A']:.1e}  [{time.time()-t0:.1f}s]")
        a, b, rms_l = lsq_fit_log(ns, sts)
        c, al, rms_p = lsq_fit_pow(ns, sts)
        camps.append({"tag": tag, "lam2": lam2, "ns": ns, "ST": sts,
                      "CR": crs, "fit_log": {"a": a, "b": b, "rms": rms_l},
                      "fit_pow": {"c": c, "alpha": al, "rms": rms_p}})
        print(f"     fit ST = {a:+.4f} + {b:.4f}*log2(n)  rms {rms_l:.4f}"
              f"   |   ST = {c:.4f} * n^{al:.4f}  rms {rms_p:.4f}")
    out["E4_scaling"] = camps
    save("mitaxsk_partE.json", out)
    return out


# ------------------------------------------------------------------- part F
def sawin_atom_profile(n, ubar, theta):
    P, p = sawin_profile(n, ubar, theta)
    lgp = [log(x) for x in P]
    lpp = [(log(q), log(1.0 - q)) for q in p]
    return [lse([g + w * a + (n - w) * b for g, (a, b) in zip(lgp, lpp)])
            for w in range(n + 1)], P, p


def part_F():
    print("\n" + "=" * 78)
    print("SK-F. UNTESTED CELL: tilt-recipe CR on Sawin gadgets, large n")
    print("=" * 78)
    out = []
    jobs = [
        # (ubar, theta, n, lambda grid)
        (0.40, 0.05, 60, (0.1, 0.25, 0.5, 1.0, 1.5, 2.5)),
        (0.40, 0.05, 100, (0.1, 0.25, 0.5, 1.0, 1.5, 2.5)),
        (0.40, 0.05, 200, (0.1, 0.25, 0.5, 1.0, 1.5, 2.5, 3.5)),
        (0.3823, 0.02, 100, (0.1, 0.25, 0.5, 1.0, 1.5, 2.5)),
        (0.3823, 0.02, 200, (0.25, 1.0, 2.5)),
        (0.3823, 0.02, 300, (0.25, 1.0, 2.5, 3.5)),
    ]
    for ubar, theta, n, lams in jobs:
        lnmu_atom, P, p = sawin_atom_profile(n, ubar, theta)
        H, _ = sawin_H_profile(n, P, p)
        HP = -sum(Pk * log(Pk) for Pk in P) / LN2
        mk = [0.5 if pk <= 0.5 else pk for pk in p]
        crlb_block = n * sum(Pk * hb(m) for Pk, m in zip(P, mk)) - H
        marg = sum(Pk * pk for Pk, pk in zip(P, p))
        rows = []
        best = (-inf, None)
        for lam2 in lams:
            t0 = time.time()
            g, resid = sinkhorn_profile(n, lam2, lnmu_atom)
            r = exch_cr(n, lam2, g)
            rows.append({"lam2": lam2, "CR": r["CR"], "gain": r["gain"],
                         "tax_A": r["tax_A"], "ST": r["second_tax"],
                         "sink_resid": resid,
                         "dH_check": abs(r["H_mu"] - H)})
            if r["CR"] > best[0]:
                best = (r["CR"], lam2)
            print(f"  ({n},{ubar},{theta}) lam={lam2:4g}: CR = "
                  f"{r['CR']:+9.4f}  gain = {r['gain']:+9.4f}  tax_A = "
                  f"{r['tax_A']:7.4f}  ST = {r['second_tax']:8.4f}  "
                  f"resid {resid:.1e}  [{time.time()-t0:.0f}s]")
        verdict = "SURVIVES" if best[0] > 0 else "TILT-CR NEGATIVE"
        print(f"  ==> n={n} marg={marg:.4f}: sup_lam CR_tilt = "
              f"{best[0]:+.4f} at lam={best[1]}  ({verdict}); "
              f"block CR_lb = {crlb_block:+.4f}\n")
        out.append({"n": n, "ubar": ubar, "theta": theta, "marg": marg,
                    "H_mu": H, "H_P": HP, "rows": rows,
                    "best_CR_tilt": best[0], "best_lam2": best[1],
                    "block_CR_lb": crlb_block, "verdict": verdict})
    save("mitaxsk_partF.json", out)
    return out


def main():
    parts = sys.argv[1] if len(sys.argv) > 1 else "ABCDEF"
    t0 = time.time()
    if "A" in parts:
        part_A()
    if "B" in parts:
        part_B()
    if "C" in parts:
        part_C()
    if "D" in parts:
        part_D()
    if "E" in parts:
        part_E()
    if "F" in parts:
        part_F()
    print(f"\n[total time: {time.time()-t0:.1f}s]")


if __name__ == "__main__":
    main()
