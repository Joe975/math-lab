#!/usr/bin/env python3
"""uc_mitax.py -- Gap 2 probe: the mutual-information tax of the
dependent-couplings route (attempts/009-mutual-information-tax.md).

Gap 2 as recorded in 003 ("Why it failed", Gap 2) is one vague sentence:
the tax  T_A(pi) = sum_i I(A_i; B_{<i} | A_{<i})  "must be beaten by the
per-coordinate gains; no analogue of the (KEY) inequality with tax is
known".  This probe pins the statement down and runs its first
falsifiable tests.

Precise objects (fixed coordinate order 1..n; histories a = A_{<i},
b = B_{<i}; pi a coupling with BOTH marginals mu):

  z_i(a,b)   = Pr_pi[A_i = 0, B_i = 0 | a, b]     (union bit is 1 - z)
  xt_i(a,b)  = Pr_pi[A_i = 0 | a, b]              (A bit, joint history)
  CR(mu,pi)  = sum_i E_pi[h(z_i)] - H(mu)         (chain-rule assembly value)
  T_A(pi)    = H(mu) - sum_i E_pi[h(xt_i)]        (the mutual-information tax
             = sum_i I(A_i; B_{<i} | A_{<i}), an identity, not a bound)
  ST(pi)     = Gain - CR                          (second tax:
             = sum_i I(U_i; (A_{<i},B_{<i}) | U_{<i}) >= 0)

Rigorous facts used (one line each):
  (i)  H_pi(U) >= sum_i E_pi[h(z_i)]: chain rule for U plus "conditioning
       on the finer sigma-field (a,b) >= U_{<i} reduces entropy".
       Hence CR <= Gain and CR > 0 certifies Gain > 0.
  (ii) CR = sum_i E_pi[h(z_i) - h(xt_i)] - T_A(pi): per-coordinate gains
       minus tax -- Gap 2 says exactly "CR > 0".
  (iii) T_A(pi) <= I_pi(A;B)  (chain rule:
       I(A_i;B|A_{<i}) >= I(A_i;B_{<i}|A_{<i}), summed).

Candidate statement under test (constructed here; 003/004 leave Gap 2
without a precise inequality):

  (TAX at p): for every mu with H(mu) > 0 and all element marginals < p,
  the route's recipe coupling pi(mu) has CR(mu, pi(mu)) > 0.

(TAX at p) implies (S-coup at p) implies Frankl at p (licensing lemma,
003/004).  "Recipe" is per-genre, as in 003: Sinkhorn/pure overlap tilt
with a lambda sweep on generic/slice mu; block-adaptive shared-component
couplings on product mixtures; half-mixing on the {0}u[1/2,1) genre.

Parts:
  A. validation of the generic full-history CR evaluator on products
     with per-coordinate Plackett(rho) couplings (closed forms: CR =
     n(h(z_rho)-h(p)), T_A = 0, ST = 0) and on the pure-slice tilt at
     n=8 against the checked-in 003_partB2.json gain.
  B. 4-atom crash family of 005 Prop 3 (the Gap-1 killer), Sinkhorn
     tilt, lambda sweep: is sup_lambda CR > 0?  Cross-check: the crash
     odds ratio 2^{lam(3-n)} is reproduced from the fitted coupling.
  C. half-mixing genre of 004 R1 (delta_empty (+) Bern(p_h)^n, the
     genre where 003's mini-theorem (e) coupling proves nothing):
     exact CR of the half-mixing coupling, n sweep -> Theta(n)?
  D. Sawin Prop-6 geometric mixtures, exactly the 002/003 certificate
     instances (n = 2000/20000/60000): certified lower bound
     CR_block >= n*sum_k P_k h(m_k) - H(mu) for the 003 part-D block
     coupling (rigorous: E[h(z_i)] = H(U_i|a,b) >= H(U_i|a,b,k)
     = sum_k P_k h(m_k) since (A_i,B_i) is independent of histories
     given k).  This is the 002-no-go smoothing-sensitivity check: the
     tax-paid functional must respond at Theta(n) on the no-go's own
     certificate instances.
  E. Chase-Lovett pure slice under the pure overlap tilt: exact CR via
     a (prefix-weights) dynamic program, n up to 80: sign of
     sup_lambda CR, and the scaling of the losses (T_A, ST) in n.

Usage: python3 uc_mitax.py            (all parts, ~1 min, deterministic)
Standard library only.  Checkpoints: ../data/mitax_part[A-E].json.
"""

import json
import os
import time
from math import inf, lgamma, log, log2, sqrt

LN2 = log(2.0)
PSI = (3.0 - sqrt(5.0)) / 2.0
RECORD = 0.38271
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data")


def h(p):
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * log2(p) - (1.0 - p) * log2(1.0 - p)


def log2C(n, k):
    if k < 0 or k > n:
        return -inf
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / LN2


def lse2(xs):
    m = max(xs)
    if m == -inf:
        return m
    return m + log2(sum(2.0 ** (x - m) for x in xs))


def save(name, obj):
    os.makedirs(DATA, exist_ok=True)
    path = os.path.join(DATA, name)
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    print(f"  [checkpoint: {os.path.relpath(path)}]")


# ------------------------------------------------------------------ core
def cr_eval(n, pairs):
    """Full-history chain-rule accounting for an explicit coupling.

    pairs: dict (Amask, Bmask) -> weight, weights summing to ~1, both
    marginals equal (checked).  Returns every Gap-2 quantity exactly.
    """
    ma, mb = {}, {}
    tot = 0.0
    for (A, B), w in pairs.items():
        ma[A] = ma.get(A, 0.0) + w
        mb[B] = mb.get(B, 0.0) + w
        tot += w
    keys = set(ma) | set(mb)
    marg_dev = max(abs(ma.get(k, 0.0) - mb.get(k, 0.0)) for k in keys)
    Hmu = -sum(w * log2(w) for w in ma.values() if w > 0.0)
    ulaw = {}
    for (A, B), w in pairs.items():
        ulaw[A | B] = ulaw.get(A | B, 0.0) + w
    HU = -sum(w * log2(w) for w in ulaw.values() if w > 0.0)
    freqs = [sum(w for A, w in ma.items() if (A >> i) & 1)
             for i in range(n)]
    Ehz, Ehx, Ehy = [], [], []
    for i in range(n):
        mask = (1 << i) - 1
        groups = {}
        for (A, B), w in pairs.items():
            key = (A & mask, B & mask)
            g = groups.get(key)
            if g is None:
                g = groups[key] = [0.0, 0.0, 0.0, 0.0]
            g[0] += w
            ai, bi = (A >> i) & 1, (B >> i) & 1
            if not ai and not bi:
                g[1] += w
            if not ai:
                g[2] += w
            if not bi:
                g[3] += w
        ez = ex = ey = 0.0
        for m_, z_, x_, y_ in groups.values():
            if m_ > 0.0:
                ez += m_ * h(z_ / m_)
                ex += m_ * h(x_ / m_)
                ey += m_ * h(y_ / m_)
        Ehz.append(ez)
        Ehx.append(ex)
        Ehy.append(ey)
    CR = sum(Ehz) - Hmu
    gain = HU - Hmu
    return {"n": n, "H_mu": Hmu, "H_U": HU, "gain": gain, "CR": CR,
            "tax_A": Hmu - sum(Ehx), "tax_B": Hmu - sum(Ehy),
            "second_tax": gain - CR, "mass": tot, "marg_dev": marg_dev,
            "max_freq": max(freqs), "Ehz": Ehz}


def atom_sinkhorn(atoms, la, lam2, tol=1e-12, itmax=20000):
    """Fit g so pi(A,B) = 2^(g_A + g_B + lam2 |A&B|) has marginals 2^la.
    Pure python; fine for small supports.  Returns pairs dict + resid."""
    K = [[lam2 * bin(a & b).count("1") for b in atoms] for a in atoms]
    g = [x / 2.0 for x in la]
    resid = inf
    for _ in range(itmax):
        gnew = []
        for r, lar in zip(K, la):
            gnew.append(lar - lse2([g[j] + r[j] for j in range(len(atoms))]))
        resid = max(abs(a - b) for a, b in zip(g, gnew))
        g = [(a + b) / 2.0 for a, b in zip(g, gnew)]
        if resid < tol:
            break
    pairs = {}
    for ia, a in enumerate(atoms):
        for ib, b in enumerate(atoms):
            pairs[(a, b)] = 2.0 ** (g[ia] + g[ib] + K[ia][ib])
    tot = sum(pairs.values())
    for k in pairs:
        pairs[k] /= tot
    return pairs, resid


def z_rho(x, y, rho):
    """Both-zero probability of the Plackett coupling, odds ratio rho."""
    if x <= 0 or y <= 0:
        return 0.0
    if x >= 1:
        return y
    if y >= 1:
        return x
    if abs(rho - 1.0) < 1e-12:
        return x * y
    a = rho - 1.0
    S = 1.0 + a * (x + y)
    disc = S * S - 4.0 * a * rho * x * y
    z = (S - sqrt(disc)) / (2.0 * a)
    return min(max(z, max(0.0, x + y - 1.0)), min(x, y))


# ---------------------------------------------------------------- part A
def part_A():
    print("=" * 78)
    print("A. Evaluator validation (products + Plackett; slice tilt vs 003)")
    print("=" * 78)
    out = {}
    n, p = 8, 0.40
    x = 1.0 - p
    checks = []
    for rho in [1.0, 1.5, 3.0]:
        z = z_rho(x, x, rho)
        cell = {(0, 0): z, (0, 1): x - z, (1, 0): x - z,
                (1, 1): 1.0 - 2.0 * x + z}
        pairs = {}
        for A in range(1 << n):
            for B in range(1 << n):
                w = 1.0
                for i in range(n):
                    w *= cell[((A >> i) & 1, (B >> i) & 1)]
                if w > 0.0:
                    pairs[(A, B)] = w
        r = cr_eval(n, pairs)
        cr_pred = n * (h(z) - h(p))
        ok = (abs(r["CR"] - cr_pred) < 1e-9 and abs(r["tax_A"]) < 1e-9
              and abs(r["second_tax"]) < 1e-9 and r["marg_dev"] < 1e-12)
        checks.append({"rho": rho, "CR": r["CR"], "CR_closed": cr_pred,
                       "tax_A": r["tax_A"], "second_tax": r["second_tax"],
                       "ok": ok})
        print(f"  product Bern({p})^{n}, Plackett rho={rho:4g}: "
              f"CR = {r['CR']:+.9f} vs closed {cr_pred:+.9f}, "
              f"tax_A = {r['tax_A']:.2e}, ST = {r['second_tax']:.2e}  "
              f"{'OK' if ok else 'FAIL'}")
    out["product_plackett"] = checks

    # pure slice tilt at n=8, w0=3, lam2=0.7: gain must match 003_partB2
    n, w0, lam2 = 8, 3, 0.7
    rho = 2.0 ** lam2
    atoms = [m for m in range(1 << n) if bin(m).count("1") == w0]
    pairs = {}
    for a in atoms:
        for b in atoms:
            pairs[(a, b)] = rho ** bin(a & b).count("1")
    tot = sum(pairs.values())
    for k in pairs:
        pairs[k] /= tot
    r = cr_eval(n, pairs)
    ref = 1.5857909077473256          # data/003_partB2.json gain_direct
    okg = abs(r["gain"] - ref) < 1e-9
    print(f"  slice n=8 w0=3 lam2=0.7: gain = {r['gain']:+.9f} "
          f"(003_partB2: {ref:+.9f}) {'OK' if okg else 'FAIL'}; "
          f"CR = {r['CR']:+.6f}, tax_A = {r['tax_A']:.6f}, "
          f"ST = {r['second_tax']:.6f}")
    out["slice_n8"] = {"gain": r["gain"], "gain_003_partB2": ref,
                       "ok": okg, "CR": r["CR"], "tax_A": r["tax_A"],
                       "second_tax": r["second_tax"]}
    out["slice_n8"]["dp_crosscheck"] = None  # filled by part E
    save("mitax_partA.json", out)
    return out


# ---------------------------------------------------------------- part B
def crash_family(n):
    """005 Prop 3: S1={2}, S2={3..n}, S3=[n], S4={1} (1-indexed coords;
    bit i-1 = coordinate i), masses (0.33, 0.33, 0.04, 0.30)."""
    S1 = 1 << 1
    S2 = ((1 << n) - 1) ^ 0b11
    S3 = (1 << n) - 1
    S4 = 1 << 0
    atoms = [S1, S2, S3, S4]
    mass = [0.33, 0.33, 0.04, 0.30]
    return atoms, [log2(m) for m in mass]


def part_B():
    print("\n" + "=" * 78)
    print("B. 4-atom crash family (005 Prop 3), Sinkhorn tilt, lambda sweep")
    print("=" * 78)
    print("Marginals <= 0.37 < record; H(mu) = 1.76 bits; the family that")
    print("crashed the pointwise OR control.  Question: sup_lam CR > 0 ?\n")
    LAMS = [-1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 1.5, 2.0, 3.0]
    out = []
    for n in [5, 8, 11]:
        atoms, la = crash_family(n)
        rows = []
        best = (-inf, None)
        for lam2 in LAMS:
            pairs, resid = atom_sinkhorn(atoms, la, lam2)
            r = cr_eval(n, pairs)
            # crash OR cross-check at coordinate 2, history (A_1=0, B_1=1)
            or_log = None
            if lam2 != 0.0:
                cells = {}
                for (A, B), w in pairs.items():
                    if (A & 1) == 0 and (B & 1) == 1:
                        cells[((A >> 1) & 1, (B >> 1) & 1)] = \
                            cells.get(((A >> 1) & 1, (B >> 1) & 1), 0.0) + w
                if all(cells.get(c, 0.0) > 0 for c in
                       [(0, 0), (0, 1), (1, 0), (1, 1)]):
                    or_log = (log2(cells[(1, 1)]) + log2(cells[(0, 0)])
                              - log2(cells[(1, 0)]) - log2(cells[(0, 1)]))
            rows.append({"lam2": lam2, "CR": r["CR"], "gain": r["gain"],
                         "tax_A": r["tax_A"], "second_tax": r["second_tax"],
                         "resid": resid, "marg_dev": r["marg_dev"],
                         "log2_OR_crash": or_log,
                         "log2_OR_pred": lam2 * (3 - n) if or_log is not None
                         else None})
            if r["CR"] > best[0]:
                best = (r["CR"], lam2)
            oc = (f"  OR: {or_log:+8.3f} (pred {lam2*(3-n):+8.3f})"
                  if or_log is not None else "")
            print(f"  n={n:2d} lam2={lam2:5.2f}: CR = {r['CR']:+8.4f}  "
                  f"gain = {r['gain']:+8.4f}  tax_A = {r['tax_A']:6.4f}  "
                  f"ST = {r['second_tax']:6.4f}{oc}")
        print(f"  ==> n={n}: sup_lam CR = {best[0]:+.4f} at lam2 = {best[1]}"
              f"  ({'TAX candidate SURVIVES' if best[0] > 0 else 'KILLED'})\n")
        out.append({"n": n, "rows": rows, "best_CR": best[0],
                    "best_lam2": best[1]})
    save("mitax_partB.json", out)
    return out


# ---------------------------------------------------------------- part C
def half_mixing_pairs(n, ph, P1):
    """004 R1 half-mixing coupling for mu = P1 delta_empty + (1-P1)
    Bern(ph)^n (needs P1 <= 1/2): draw s ~ Bern(ph)^n; with prob 2*P1
    (A,B) iid uniform on {empty, s}; else A = B = s."""
    pairs = {}

    def add(k, w):
        pairs[k] = pairs.get(k, 0.0) + w

    for s in range(1 << n):
        ps = ph ** bin(s).count("1") * (1 - ph) ** (n - bin(s).count("1"))
        add((0, 0), ps * P1 / 2.0)
        add((0, s), ps * P1 / 2.0)
        add((s, 0), ps * P1 / 2.0)
        add((s, s), ps * (1.0 - 1.5 * P1))
    return pairs


def part_C():
    print("\n" + "=" * 78)
    print("C. Half-mixing genre (004 R1): delta_empty (+) Bern(p_h)^n")
    print("=" * 78)
    print("The genre 003's mini-theorem-(e) coupling provably cannot touch")
    print("(its gain is 0); 004's half-mixing repair has gain Theta(n).")
    print("Question: does chain-rule certification survive, i.e. CR = ")
    print("Theta(n) too, or does the tax eat the repair?\n")
    out = []
    for ph, P1, tag in [(0.60, 0.50, "004's mu_n (marginals 0.30)"),
                        (0.51, 0.163, "recipe-killer shape, cap 0.427")]:
        print(f"  {tag}: p_h = {ph}, P1 = {P1}")
        theory = P1 / 2.0 * h(ph)     # per-coordinate gain slope (004 R1)
        rows = []
        for n in [6, 8, 10, 12, 14]:
            pairs = half_mixing_pairs(n, ph, P1)
            r = cr_eval(n, pairs)
            gain_pred = (n * theory
                         - (h(P1) - h(P1 / 2.0) * 1.0)
                         - (P1 / 2.0) * 0.0)
            # closed form from 004 R1: (P1/2) n h(ph) - [H2(P1)-H2(P1/2)]
            gain_pred = P1 / 2.0 * n * h(ph) - (h(P1) - h(P1 / 2.0))
            rows.append({"n": n, "CR": r["CR"], "gain": r["gain"],
                         "gain_004_closed": gain_pred,
                         "tax_A": r["tax_A"],
                         "second_tax": r["second_tax"],
                         "max_freq": r["max_freq"],
                         "marg_dev": r["marg_dev"]})
            print(f"    n={n:2d}: CR = {r['CR']:+8.4f}  gain = "
                  f"{r['gain']:+8.4f} (004 closed form {gain_pred:+8.4f})  "
                  f"tax_A = {r['tax_A']:6.4f}  ST = {r['second_tax']:6.4f}  "
                  f"maxfreq = {r['max_freq']:.3f}")
        dCR = rows[-1]["CR"] - rows[-2]["CR"]
        print(f"    slope dCR/dn (last step) = {dCR/2:+.4f} vs per-coord "
              f"gain slope (P1/2)h(p_h) = {theory:+.4f}\n")
        out.append({"ph": ph, "P1": P1, "tag": tag, "rows": rows,
                    "slope_last": dCR / 2.0, "slope_theory": theory})
    save("mitax_partC.json", out)
    return out


# ---------------------------------------------------------------- part D
def part_D():
    print("\n" + "=" * 78)
    print("D. Sawin certificates: certified CR lower bound, block coupling")
    print("   (= the 002-no-go smoothing-sensitivity check)")
    print("=" * 78)
    print("Rigorous: E[h(z_i)] = H(U_i|a,b) >= H(U_i|a,b,k) = sum_k P_k")
    print("h(m_k) for the shared-k block coupling of 003 part D (given k,")
    print("(A_i,B_i) is independent of histories).  So CR >= n sum_k P_k")
    print("h(m_k) - H(mu) exactly -- computed below on 002's certificate")
    print("instances, whose H(A) values are re-derived and matched.\n")
    ref_HA = {2000: 1927.70, 20000: 19239.59, 60000: 57578.99}  # 002 table
    out = []
    print(f"  {'n':>6} {'ubar':>7} {'theta':>6} {'maxmarg':>9} "
          f"{'H(mu)':>10} {'H(P)':>6} {'CR_lb':>10}")
    for n, ubar, theta in [(2000, 0.390, 0.05), (20000, 0.386, 0.02),
                           (60000, 0.3823, 0.001)]:
        K = 0
        while theta ** (K + 1) > 1e-28:
            K += 1
        raw = [(1 - theta) * theta ** k for k in range(K + 1)]
        Z = sum(raw)
        Pk = [x / Z for x in raw]
        pk = [1.0 - (1.0 - ubar) ** (k + 1) for k in range(K + 1)]
        lg = [log2(P) for P in Pk]
        lp = [(log2(p), log2(1.0 - p)) for p in pk]
        H = 0.0
        tot = 0.0
        for w in range(n + 1):
            lc = log2C(n, w)
            la = lse2([g + w * a + (n - w) * b
                       for g, (a, b) in zip(lg, lp)])
            mass = 2.0 ** (lc + la)
            tot += mass
            H -= mass * la
        assert abs(tot - 1.0) < 1e-6, tot
        mk = [0.5 if p <= 0.5 else p for p in pk]
        HP = -sum(P * log2(P) for P in Pk)
        crlb = n * sum(P * h(m) for P, m in zip(Pk, mk)) - H
        marg = max(pk[0], 0) and sum(P * p for P, p in zip(Pk, pk))
        okH = abs(H - ref_HA[n]) < 0.5
        print(f"  {n:6d} {ubar:7.4f} {theta:6.3f} {marg:9.6f} "
              f"{H:10.2f} {HP:6.3f} {crlb:+10.2f}   "
              f"H(mu) {'matches 002' if okH else 'MISMATCH vs 002!'}")
        out.append({"n": n, "ubar": ubar, "theta": theta, "marg": marg,
                    "H_mu": H, "H_mu_002": ref_HA[n], "H_P": HP,
                    "CR_lower_bound": crlb, "H_matches_002": okH})
    print("\n  Reading: the tax-paid chain-rule functional is Theta(n)")
    print("  POSITIVE on exactly the instances that certify the 002 no-go")
    print("  (where every smoothing-insensitive functional is Theta(n)")
    print("  negative) -- the functional responds at probability scale,")
    print("  i.e. it is smoothing-sensitive; the no-go does not apply.")
    save("mitax_partD.json", out)
    return out


# ---------------------------------------------------------------- part E
def make_logS(n, lam2):
    """log2 S(m, ra, rb) = log2 sum_j C(m,j) C(m-j, ra-j) C(m-ra, rb-j)
    * 2^(lam2 j): tilted count of pairs (x, y) in ([m] choose ra) x
    ([m] choose rb) weighted by 2^(lam2 |x n y|).  Memoized."""
    cache = {}

    def logS(m, ra, rb):
        if ra < 0 or rb < 0 or ra > m or rb > m:
            return -inf
        key = (m, ra, rb)
        v = cache.get(key)
        if v is None:
            terms = [log2C(m, j) + log2C(m - j, ra - j)
                     + log2C(m - ra, rb - j) + lam2 * j
                     for j in range(max(0, ra + rb - m), min(ra, rb) + 1)]
            v = lse2(terms) if terms else -inf
            cache[key] = v
        return v
    return logS


def slice_gain(n, w0, lam2):
    """H_pi(U) - H(mu) for the pure slice tilt (independent
    reimplementation of the 003 closed form)."""
    js = range(max(0, 2 * w0 - n), w0 + 1)
    lw = [log2C(w0, j) + log2C(n - w0, w0 - j) + lam2 * j for j in js]
    Z = lse2(lw)
    HU = 0.0
    for j, l in zip(js, lw):
        pj = 2.0 ** (l - Z)
        if pj > 0.0:
            HU += pj * (log2C(n, 2 * w0 - j) - (l - Z))
    return HU - log2C(n, w0)


def slice_cr(n, w0, lam2, prune=60.0):
    """Exact CR / tax accounting for pi ∝ 2^(lam2 |A n B|) on
    slice(w0) x slice(w0), via the sufficient prefix state (wa, wb):
    conditionals depend on the prefix pair only through its two weights
    (the 2^(lam2 wj) prefix factor cancels)."""
    logS = make_logS(n, lam2)
    lZ = logS(n, w0, w0)
    sum_hz = sum_hx = 0.0
    mass_min = inf
    cons_err = 0.0
    for l in range(n):           # prefix length; coordinate i = l + 1
        m1 = n - l - 1
        mass = 0.0
        for wa in range(0, min(l, w0) + 1):
            lpa = logS(l, wa, 0) if False else None
            for wb in range(0, min(l, w0) + 1):
                lp = (logS(l, wa, wb)
                      + logS(n - l, w0 - wa, w0 - wb) - lZ)
                if lp < -prune:
                    continue
                P = 2.0 ** lp
                lN = {}
                for al in (0, 1):
                    for be in (0, 1):
                        lN[(al, be)] = (lam2 * al * be
                                        + logS(m1, w0 - wa - al,
                                               w0 - wb - be))
                ltot = lse2(list(lN.values()))
                cons_err = max(cons_err,
                               abs(ltot - logS(n - l, w0 - wa, w0 - wb)))
                z = 2.0 ** (lse2([lN[(0, 0)]]) - ltot)
                x = 2.0 ** (lse2([lN[(0, 0)], lN[(0, 1)]]) - ltot)
                sum_hz += P * h(z)
                sum_hx += P * h(x)
                mass += P
        mass_min = min(mass_min, mass)
    Hmu = log2C(n, w0)
    CR = sum_hz - Hmu
    gain = slice_gain(n, w0, lam2)
    return {"n": n, "w0": w0, "lam2": lam2, "CR": CR, "gain": gain,
            "tax_A": Hmu - sum_hx, "second_tax": gain - CR,
            "mass_min": mass_min, "consistency_err": cons_err}


def part_E(partA_out):
    print("\n" + "=" * 78)
    print("E. Chase-Lovett pure slice, pure overlap tilt: exact CR (DP)")
    print("=" * 78)
    print("The sufficient prefix state is (wa, wb); DP validated against")
    print("the atom-level evaluator at n=8 and the identity")
    print("sum_cells N = S(future).  Sweep lambda, track losses vs n.\n")

    # cross-validation against part A's atom-level slice run
    r_dp = slice_cr(8, 3, 0.7)
    ra = partA_out["slice_n8"]
    ok = (abs(r_dp["CR"] - ra["CR"]) < 1e-9
          and abs(r_dp["gain"] - ra["gain"]) < 1e-9
          and abs(r_dp["tax_A"] - ra["tax_A"]) < 1e-9)
    print(f"  DP vs atom evaluator (n=8, w0=3, lam2=0.7): "
          f"|dCR| = {abs(r_dp['CR']-ra['CR']):.2e}, "
          f"|dgain| = {abs(r_dp['gain']-ra['gain']):.2e}  "
          f"{'OK' if ok else 'FAIL'}\n")
    out = {"dp_crosscheck_ok": ok, "sweeps": []}

    LAMS = [0.0, 0.3, 0.6, 1.0, 1.5, 2.0, 3.0]
    for n, p in [(24, 0.42), (40, 0.42), (60, 0.42), (80, 0.42),
                 (60, 0.3823), (80, 0.3823)]:
        w0 = round(p * n)
        rows = []
        best = (-inf, None)
        for lam2 in LAMS:
            r = slice_cr(n, w0, lam2)
            rows.append(r)
            if r["CR"] > best[0]:
                best = (r["CR"], lam2)
            print(f"  n={n:3d} w0={w0:3d} lam2={lam2:4.1f}: "
                  f"CR = {r['CR']:+9.4f}  gain = {r['gain']:+9.4f}  "
                  f"tax_A = {r['tax_A']:7.4f}  ST = {r['second_tax']:7.4f}")
        tag = "SURVIVES" if best[0] > 0 else "KILLED at this n"
        print(f"  ==> n={n}, p={w0/n:.4f}: sup_lam CR = {best[0]:+.4f} "
              f"at lam2 = {best[1]}   ({tag})\n")
        out["sweeps"].append({"n": n, "p": w0 / n, "w0": w0, "rows": rows,
                              "best_CR": best[0], "best_lam2": best[1]})
    save("mitax_partE.json", out)
    return out


def main():
    t0 = time.time()
    a = part_A()
    part_B()
    part_C()
    part_D()
    part_E(a)
    print(f"\n[total time: {time.time()-t0:.1f}s]")


if __name__ == "__main__":
    main()
