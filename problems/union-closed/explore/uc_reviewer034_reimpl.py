#!/usr/bin/env python3
"""Attempt 034 (reviewer pass on 030/031/033): independent
re-implementation.  Default stance: REFUTE.

Nothing here imports 030/031/033's builders, evaluators or optimizers
(`uc_hu_probe`, `uc_hu_attack`, `uc_hu_certify`, `uc_latentK_probe`,
`uc_mitax.cr_eval`, `uc_branch_skeptic.cr_chain`, `uc_branch_attack`).
Every coupling and every functional below is written from the prose
definitions in the three records.  The only things read from the repo
are *instance data* (`data/cr_attack.json` weights, `data/hu_attack.json`
witness weights, `data/latentK_probe.json` claimed values) and the
published crash family (005 Prop 3), which is restated here.

Third construction path for HU: neither the probe's prefix-cell forward
product nor the skeptic's recursive descent, but a direct closed-form
product over atom pairs, pi(A, B) = prod_i c_i(A_i, B_i | A_<i, B_<i).
Evaluator: own prefix grouping, natural logs.

Parts:
  P1  HU on the five cr_attack floor endpoints, the crash family, and
      products: CR, Lemma HU-notax (T_A = T_B = 0), the per-cell ledger
      identity, and CR = n(1-h(p)).
  P2  own slice DP (2-D array, natural logs) at w0/n = 0.40, n in
      {100, 300}; cross-checked at n = 8 against P1's mask builder.
  P3  own K-component latent-mixing builder + own optimizer: two K = 3
      battery rows, and the degeneration claim vs an own re-derivation
      of the fused 2-block EM limit.
  P4  the 033 order-kill: witness rebuilt from data/hu_attack.json, all
      24 orders enumerated, per-coordinate ledger dissection, and an
      own min-over-all-24-orders descent for the cap trace / crossing.
  P5  high-precision (decimal, 60 digits) exact-rational recomputation
      of the three 033 certificates, checked against the committed
      enclosures; plus the rationalization total-variation claim.

Stdlib only, deterministic.  Usage: python uc_reviewer034_reimpl.py
Runtime ~2 min.  Output: ../data/uc_reviewer034_reimpl_out.txt
"""
from __future__ import annotations

import itertools
import json
import math
from decimal import Decimal, getcontext
from fractions import Fraction
from math import lgamma, log
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
LN2 = log(2.0)
FAILS = []
LINES = []


def say(m=""):
    print(m, flush=True)
    LINES.append(m)


def check(name, ok, detail=""):
    say(f"  [{'ok' if ok else 'REFUTED'}] {name}" + (f"  ({detail})" if detail else ""))
    if not ok:
        FAILS.append(name)


def hb(p):
    """binary entropy in bits, natural-log implementation"""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return (-p * log(p) - (1.0 - p) * log(1.0 - p)) / LN2


# ===================================================================== P1
def prefix_cond(n, mu):
    """cond[i][prefix] = (mass, mass with A_i = 0)."""
    out = []
    for i in range(n):
        mask = (1 << i) - 1
        d = {}
        for a, w in mu.items():
            e = d.setdefault(a & mask, [0.0, 0.0])
            e[0] += w
            if not (a >> i) & 1:
                e[1] += w
        out.append(d)
    return out


def hu_joint(n, mu):
    """HU coupling as an explicit measure on atom pairs, built by the
    direct product formula (no DP, no recursion)."""
    cond = prefix_cond(n, mu)
    supp = sorted(mu)
    pi = {}
    for A in supp:
        for B in supp:
            w = 1.0
            for i in range(n):
                mask = (1 << i) - 1
                ea, eb = cond[i].get(A & mask), cond[i].get(B & mask)
                if ea is None or eb is None or ea[0] <= 0.0 or eb[0] <= 0.0:
                    w = 0.0
                    break
                x = ea[1] / ea[0]
                y = eb[1] / eb[0]
                z = min(max(0.5, x + y - 1.0), x, y)
                ai, bi = (A >> i) & 1, (B >> i) & 1
                c = (z if bi == 0 else x - z) if ai == 0 else \
                    (y - z if bi == 0 else 1.0 - x - y + z)
                if c <= 0.0:
                    w = 0.0
                    break
                w *= c
            if w > 0.0:
                pi[(A, B)] = pi.get((A, B), 0.0) + w
    return pi


def evaluate(n, pi):
    """Own chain-rule accounting: CR, taxes, gain, ledger."""
    ma, mb, ul = {}, {}, {}
    tot = 0.0
    for (A, B), w in pi.items():
        ma[A] = ma.get(A, 0.0) + w
        mb[B] = mb.get(B, 0.0) + w
        ul[A | B] = ul.get(A | B, 0.0) + w
        tot += w
    H = -sum(w * log(w) for w in ma.values() if w > 0.0) / LN2
    HU = -sum(w * log(w) for w in ul.values() if w > 0.0) / LN2
    ledger = []
    Sz = Sx = Sy = 0.0
    for i in range(n):
        mask = (1 << i) - 1
        g = {}
        for (A, B), w in pi.items():
            k = (A & mask, B & mask)
            e = g.setdefault(k, [0.0, 0.0, 0.0, 0.0])
            e[0] += w
            ai, bi = (A >> i) & 1, (B >> i) & 1
            if not ai and not bi:
                e[1] += w
            if not ai:
                e[2] += w
            if not bi:
                e[3] += w
        ez = ex = ey = 0.0
        for m_, z_, x_, y_ in g.values():
            if m_ > 0.0:
                ez += m_ * hb(z_ / m_)
                ex += m_ * hb(x_ / m_)
                ey += m_ * hb(y_ / m_)
        ledger.append(ez - 0.5 * (ex + ey))
        Sz += ez
        Sx += ex
        Sy += ey
    mdev = max(abs(ma.get(k, 0.0) - mb.get(k, 0.0))
               for k in set(ma) | set(mb))
    return {"H": H, "H_U": HU, "CR": Sz - H, "T_A": H - Sx, "T_B": H - Sy,
            "gain": HU - H, "ST": (HU - H) - (Sz - H), "mass": tot,
            "marg_dev": mdev, "ledger": ledger}


def norm(mu):
    t = sum(mu.values())
    return {a: w / t for a, w in mu.items() if w / t > 0.0}


def crash_family_local(n):
    """005 Prop 3 as restated in the 031/033 records: S1 = {2},
    S2 = {3..n}, S3 = [n], S4 = {1}; masses .33/.33/.04/.30."""
    return {1 << 1: 0.33, ((1 << n) - 1) ^ 0b11: 0.33,
            (1 << n) - 1: 0.04, 1 << 0: 0.30}


def maxmarg(n, mu):
    t = sum(mu.values())
    return max(sum(w for a, w in mu.items() if (a >> i) & 1)
               for i in range(n)) / t


def floor_instances():
    d = json.loads((DATA / "cr_attack.json").read_text())
    out = []
    for row in d["rows"]:
        if row["H_mu"] < 0.5:
            continue
        out.append((row["seed"], row["n"],
                    norm({int(s, 2): w for s, w in row["mu"].items()}),
                    row["sup_cr"]))
    return out


def part1():
    say("P1. HU rebuilt (direct pair-product), own evaluator")
    claimed = {"windowkill": 0.3090, "mmabskill": 0.0822,
               "mmabskill_n6": 0.3007, "random0": 0.1563,
               "random1": 0.2991}
    for name, n, mu, prior in floor_instances():
        pi = hu_joint(n, mu)
        r = evaluate(n, pi)
        led = sum(r["ledger"])
        check(f"floor:{name:13s} n={n}: CR_HU = {claimed[name]:+.4f}",
              abs(r["CR"] - claimed[name]) < 5e-5
              and r["marg_dev"] < 1e-12 and abs(r["mass"] - 1) < 1e-12
              and r["CR"] > prior,
              f"CR={r['CR']:+.6f} (031 sweep {prior:+.4f}), mdev={r['marg_dev']:.0e}")
        check(f"floor:{name:13s} HU-notax T_A = T_B = 0 and ledger",
              abs(r["T_A"]) < 1e-12 and abs(r["T_B"]) < 1e-12
              and abs(led - r["CR"]) < 1e-12,
              f"T_A={r['T_A']:+.2e} T_B={r['T_B']:+.2e} "
              f"ledger-CR={led - r['CR']:+.2e}")
    say()
    vals = []
    for n in (5, 8, 11, 14, 17):
        mu = norm(crash_family_local(n))
        r = evaluate(n, hu_joint(n, mu))
        vals.append(r["CR"])
        if n == 5:
            say(f"   crash n=5: CR={r['CR']:.10f} H={r['H']:.6f} "
                f"ratio={r['CR'] / r['H']:.6f} ST={r['ST']:+.2e}")
    check("crash family n=5..17: CR_HU = +0.19919 flat in n",
          all(abs(v - 0.19919) < 5e-6 for v in vals)
          and max(vals) - min(vals) < 1e-12,
          f"spread {max(vals) - min(vals):.1e}, value {vals[0]:.8f}")
    say()
    for p, n in ((0.3, 5), (0.38271, 6), (0.45, 5), (0.25, 4), (0.2, 4)):
        mu = {a: p ** bin(a).count("1") * (1 - p) ** (n - bin(a).count("1"))
              for a in range(1 << n)}
        r = evaluate(n, hu_joint(n, mu))
        pred = n * (1 - hb(p))
        check(f"Bern({p})^{n}: CR_HU = n(1-h(p)) = {pred:.6f}",
              abs(r["CR"] - pred) < 1e-9, f"CR={r['CR']:.9f}")
    c = (1 - hb(0.38271)) / hb(0.38271)
    check("constant (1-h(0.38271))/h(0.38271) = 0.041739",
          abs(c - 0.041739) < 5e-7, f"{c:.9f}")
    c381 = (1 - hb(0.381)) / hb(0.381)
    say(f"   (1-h(0.381))/h(0.381) = {c381:.8f}  [030 near-psi ratio anchor]")
    say()


# ===================================================================== P2
def slice_dp(n, w0):
    """Own slice DP: uniform measure on the weight-w0 slice, states are
    (ones so far in A, ones so far in B), 2-D lists, natural logs."""
    cur = [[0.0] * (w0 + 1) for _ in range(w0 + 1)]
    cur[0][0] = 1.0
    acc = 0.0
    for i in range(n):
        rem = n - i
        nxt = [[0.0] * (w0 + 1) for _ in range(w0 + 1)]
        for wa in range(w0 + 1):
            row = cur[wa]
            for wb in range(w0 + 1):
                P = row[wb]
                if P <= 0.0:
                    continue
                x = 1.0 - (w0 - wa) / rem
                y = 1.0 - (w0 - wb) / rem
                z = min(max(0.5, x + y - 1.0), x, y)
                acc += P * hb(z)
                for (ai, bi), c in (((0, 0), z), ((0, 1), x - z),
                                    ((1, 0), y - z),
                                    ((1, 1), 1.0 - x - y + z)):
                    if c <= 0.0:
                        continue
                    na, nb = wa + ai, wb + bi
                    if na <= w0 and nb <= w0:
                        nxt[na][nb] += P * c
        cur = nxt
    H = (lgamma(n + 1) - lgamma(w0 + 1) - lgamma(n - w0 + 1)) / LN2
    return acc - H, H


def part2():
    say("P2. Slice DP (own 2-D implementation)")
    n, w0 = 8, 3
    mu = norm({a: 1.0 for a in range(1 << n) if bin(a).count("1") == w0})
    r = evaluate(n, hu_joint(n, mu))
    cr, H = slice_dp(n, w0)
    check("n=8 w0=3: own DP == own mask build", abs(cr - r["CR"]) < 1e-11,
          f"DP {cr:+.9f} vs masks {r['CR']:+.9f}")
    for n, frac, ref in ((100, 0.40, 4.542), (300, 0.40, 11.071)):
        w0 = round(frac * n)
        cr, H = slice_dp(n, w0)
        check(f"n={n} w0/n={frac}: 031 claims {ref:+.3f}",
              abs(cr - ref) < 5e-4,
              f"own DP {cr:+.5f}, CR/H = {cr / H:.5f}, "
              f"per-coord {cr / n:.5f}")
    cr1, _ = slice_dp(200, 80)
    cr2, _ = slice_dp(300, 120)
    check("0.40-series slope ~ 0.032 >= product slope 1-h(0.40) = 0.029",
          abs((cr2 - cr1) / 100 - 0.032) < 1e-3
          and (cr2 - cr1) / 100 > 1 - hb(0.40),
          f"slope {(cr2 - cr1) / 100:.5f}")
    say()


# ===================================================================== P3
def adaptive_cell(p):
    """009 part D / 031 per-coordinate rule for a Bern(p) product:
    both-zero probability pushed to 1/2 inside the Frechet window
    [max(0, 2(1-p)-1), 1-p]."""
    x = 1.0 - p
    z = min(max(0.5, 2 * x - 1.0), x)
    return {(0, 0): z, (0, 1): x - z, (1, 0): x - z, (1, 1): 1 - 2 * x + z}


def bern(n, p):
    if p <= 0.0:
        return {0: 1.0}
    return {a: p ** bin(a).count("1") * (1 - p) ** (n - bin(a).count("1"))
            for a in range(1 << n)}


def latentK_joint(n, comps, Q):
    """Label-coupled mixture: diagonal (j,j) = the adaptive product
    joint for p_j, off-diagonal (j,k) = independent product of the two
    component laws."""
    K = len(comps)
    pi = {}
    laws = [bern(n, p) for _, p in comps]
    for j in range(K):
        for k in range(K):
            w = Q[j][k]
            if w <= 0.0:
                continue
            if j == k:
                cell = adaptive_cell(comps[j][1])
                for A in range(1 << n):
                    for B in range(1 << n):
                        c = w
                        for i in range(n):
                            c *= cell[((A >> i) & 1, (B >> i) & 1)]
                            if c == 0.0:
                                break
                        if c > 0.0:
                            pi[(A, B)] = pi.get((A, B), 0.0) + c
            else:
                for A, wa in laws[j].items():
                    for B, wb in laws[k].items():
                        c = w * wa * wb
                        if c > 0.0:
                            pi[(A, B)] = pi.get((A, B), 0.0) + c
    return pi


def qmat(P, t):
    K = len(P)
    Q = [[0.0] * K for _ in range(K)]
    for (j, k), v in t.items():
        Q[j][k] = Q[k][j] = v
    for j in range(K):
        d = P[j] - sum(Q[j][k] for k in range(K) if k != j)
        if d < -1e-12:
            return None
        Q[j][j] = max(d, 0.0)
    return Q


def optimize(n, comps, rounds=4, steps=6):
    """Own optimizer: nested bisection-style shrinking box on the
    symmetric transfers (different search from 030's grid+refine)."""
    P = [c[0] for c in comps]
    idx = [(j, k) for j in range(len(P)) for k in range(len(P)) if j < k]
    box = {jk: (0.0, min(P[jk[0]], P[jk[1]])) for jk in idx}
    best = (None, None)
    for _ in range(rounds):
        grids = [[lo + (hi - lo) * i / steps for i in range(steps + 1)]
                 for lo, hi in (box[jk] for jk in idx)]
        for combo in itertools.product(*grids):
            t = dict(zip(idx, combo))
            Q = qmat(P, t)
            if Q is None:
                continue
            cr = evaluate(n, latentK_joint(n, comps, Q))["CR"]
            if best[0] is None or cr > best[0]:
                best = (cr, t)
        nb = {}
        for gi, jk in enumerate(idx):
            c = best[1][jk]
            w = (box[jk][1] - box[jk][0]) / steps
            nb[jk] = (max(0.0, c - w), min(min(P[jk[0]], P[jk[1]]), c + w))
        box = nb
    return best


def em_fused_limit(n, scomps, P1):
    """Own re-derivation of the fused 2-block EM value: A = eps_A s,
    B = eps_B s with s ~ the fused s-law, joint eps with P(0,0) = q;
    ST = 0 (025 EM), so CR(q) = H(q d_0 + (1-q) S) - H(P1 d_0 + (1-P1) S).
    Measures built atom-wise; no log-sum-exp."""
    def Hmix(w):
        mu = {0: w}
        for S, p in scomps:
            for a in range(1 << n):
                k = bin(a).count("1")
                mu[a] = mu.get(a, 0.0) + (1 - w) * S * p ** k \
                    * (1 - p) ** (n - k)
        return -sum(m * log(m) for m in mu.values() if m > 0.0) / LN2
    base = Hmix(P1)
    lo = max(0.0, 2 * P1 - 1.0)
    return max(Hmix(lo + (P1 - lo) * j / 400) - base for j in range(401))


def part3():
    say("P3. K-component latent-mixing (own builder + own optimizer)")
    ref = json.loads((DATA / "latentK_probe.json").read_text())
    insts = {
        "lo+mid+hi": [(0.5, 0.01), (0.3, 0.55), (0.2, 0.9)],
        "near-psi": [(0.4, 0.36), (0.3, 0.39), (0.3, 0.40)],
        "near0+2hi": [(0.62, 0.001), (0.22, 0.9), (0.16, 0.99)],
    }
    for tag, n in (("lo+mid+hi", 4), ("near-psi", 6)):
        row = [r for r in ref["battery"] if r["tag"] == tag
               and r["n"] == n][0]
        cr, t = optimize(n, insts[tag])
        rr = evaluate(n, latentK_joint(
            n, insts[tag], qmat([c[0] for c in insts[tag]], t)))
        check(f"{tag} n={n}: 030 CR* = {row['CR_opt']:+.5f}",
              abs(cr - row["CR_opt"]) < 2e-3 and cr > 0,
              f"own optimum {cr:+.6f} (t={ {f'{j}{k}': round(v, 4) for (j, k), v in t.items()} }), "
              f"H={rr['H']:.4f} CR/H={cr / rr['H']:.5f}")
    # near-psi: HU-mix identity CR = n - H
    for n in (2, 4, 6):
        comps = insts["near-psi"]
        Q = qmat([c[0] for c in comps], {})
        r = evaluate(n, latentK_joint(n, comps, Q))
        check(f"HU-mix identity at near-psi n={n}: CR = n - H(mu)",
              abs(r["CR"] - (n - r["H"])) < 1e-12,
              f"CR={r['CR']:.9f} n-H={n - r['H']:.9f} "
              f"CR/H={r['CR'] / r['H']:.8f} H_U={r['H_U']:.9f}(=n?)")
    # degeneration
    scomps = [(0.22 / 0.38, 0.9), (0.16 / 0.38, 0.99)]
    em = em_fused_limit(4, scomps, 0.62)
    check("EM fused 2-block limit = +0.35233 (own re-derivation)",
          abs(em - 0.3523253313258097) < 1e-6, f"{em:+.7f}")
    for pl in (1e-2, 1e-3, 1e-4):
        comps = [(0.62, pl), (0.22, 0.9), (0.16, 0.99)]
        cr, t = optimize(4, comps)
        rowref = [r for r in ref["degeneration"]
                  if abs(r["p_light"] - pl) < 1e-12][0]
        check(f"degeneration p_light={pl:.0e}: K=3 optimum beats EM",
              cr > em + 1e-3,
              f"own {cr:+.5f} vs EM {em:+.5f}; 030 checkpoint CR_opt "
              f"{rowref['CR_opt']:+.5f}, 030 prose "
              f"{'+0.4037/+0.3970/+0.3950'}")
    say()


# ===================================================================== P4
def permute(n, mu, perm):
    out = {}
    for a, w in mu.items():
        b = 0
        for i in range(n):
            if (a >> i) & 1:
                b |= 1 << perm[i]
        out[b] = out.get(b, 0.0) + w
    return out


def witness_from_checkpoint():
    d = json.loads((DATA / "hu_attack.json").read_text())
    v = d["cap_0.45"]["violations"]
    assert len(v) == 1, len(v)
    row = v[0]
    return row, norm({int(s, 2): w for s, w in row["mu"].items()})


def part4():
    say("P4. The 033 order-kill, rebuilt from data/hu_attack.json")
    row, mu = witness_from_checkpoint()
    n = row["n"]
    say(f"   start={row['start']} n={n} atoms={len(mu)} "
        f"worst_order={row['worst_order']} floor={row['floor']:+.6f} "
        f"orders_tried={row['orders_tried']}")
    vals = {}
    for perm in itertools.permutations(range(n)):
        r = evaluate(n, hu_joint(n, permute(n, mu, perm)))
        vals[perm] = (r["CR"], r["H"], r["ledger"])
    neg = sorted(p for p, v in vals.items() if v[0] < 0)
    check("exactly 1 of 24 orders gives CR_HU < 0",
          len(neg) == 24 - len([1 for v in vals.values() if v[0] >= 0])
          and len(neg) == 1, f"negative orders: {neg}")
    check("the negative order is (0, 1, 3, 2)", neg == [(0, 1, 3, 2)],
          f"{neg}")
    cr_bad = vals[(0, 1, 3, 2)][0]
    cr_id = vals[(0, 1, 2, 3)][0]
    cr_best = max(v[0] for v in vals.values())
    H = vals[(0, 1, 2, 3)][1]
    check("bad-order CR_HU = -0.00805", abs(cr_bad + 0.00805) < 5e-6,
          f"{cr_bad:+.8f} (ratio {cr_bad / H:+.6f})")
    check("identity CR = +0.0148, best order = +0.2567 (ratio +0.1016)",
          abs(cr_id - 0.0148) < 5e-5 and abs(cr_best - 0.2567) < 5e-5
          and abs(cr_best / H - 0.1016) < 5e-5,
          f"identity {cr_id:+.6f}, best {cr_best:+.6f} "
          f"(ratio {cr_best / H:+.6f})")
    mm = maxmarg(n, mu)
    check("max marginal < 9/20 (= 0.44990) and H = 2.5270",
          mm < 0.45 and abs(mm - 0.44990) < 5e-6 and abs(H - 2.5270) < 5e-5,
          f"maxmarg {mm:.6f}, H {H:.6f}")
    for perm, tag in (((0, 1, 2, 3), "identity"), ((0, 1, 3, 2), "(0,1,3,2)")):
        led = vals[perm][2]
        say(f"   ledger {tag:10s}: "
            + "  ".join(f"{x:+.4f}" for x in led)
            + f"   -> {sum(led):+.4f}")
    # 033 part E: is bit 2 deterministic given the others?
    det = {}
    for a in mu:
        det.setdefault(a & ~(1 << 2), set()).add((a >> 2) & 1)
    check("033 part E: coordinate 2 is a function of the other three",
          all(len(s) == 1 for s in det.values()),
          f"{len(det)} co-histories, max branching "
          f"{max(len(s) for s in det.values())}")
    # own min-over-all-24-orders descent: the cap trace / crossing
    say("   own all-24-order descent from the random0 endpoint:")
    start = [x for x in floor_instances() if x[0] == "random0"][0]
    trace = {}
    for cap in (0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45):
        trace[cap] = own_descent(start[1], start[2], cap)
        say(f"     cap {cap:.2f}: floor min-over-24 CR/H = {trace[cap]:+.6f}")
    check("crossing of the all-orders floor lies in (0.44, 0.45]",
          all(trace[c] > 0 for c in (0.39, 0.40, 0.41, 0.42, 0.43, 0.44))
          and trace[0.45] < 0,
          f"0.44 -> {trace[0.44]:+.6f}, 0.45 -> {trace[0.45]:+.6f}")
    say()
    return mu


def own_descent(n, mu0, cap, rounds=120):
    """Deterministic multiplicative pattern search on min over ALL n!
    orders of CR_HU/H, in-regime, H >= 0.2."""
    perms = list(itertools.permutations(range(n)))

    def obj(mu):
        m = norm(mu)
        if maxmarg(n, m) >= cap:
            return None
        vs = []
        for perm in perms:
            r = evaluate(n, hu_joint(n, permute(n, m, perm)))
            if r["H"] < 0.2:
                return None
            vs.append(r["CR"] / r["H"])
        return min(vs)

    mu = dict(mu0)
    best = obj(mu)
    if best is None:
        return None
    step, it = 0.5, 0
    while step >= 0.03 and it < rounds:
        improved = False
        for a in sorted(mu):
            for f in (1 + step, 1 / (1 + step)):
                cand = dict(mu)
                cand[a] *= f
                v = obj(cand)
                if v is not None and v < best - 1e-9:
                    mu, best, improved = cand, v, True
        tot = sum(mu.values())
        for extra in (0, (1 << n) - 1):
            cand = dict(mu)
            cand[extra] = cand.get(extra, 0.0) + 0.1 * tot
            v = obj(cand)
            if v is not None and v < best - 1e-9:
                mu, best, improved = cand, v, True
        if len(mu) > 3:
            a0 = min(mu, key=mu.get)
            cand = {a: w for a, w in mu.items() if a != a0}
            v = obj(cand)
            if v is not None and v < best - 1e-9:
                mu, best, improved = cand, v, True
        if not improved:
            step *= 0.5
        it += 1
    return best


# ===================================================================== P5
getcontext().prec = 60


def dlog2(fr):
    return Decimal(fr.numerator).ln() / Decimal(2).ln() \
        - Decimal(fr.denominator).ln() / Decimal(2).ln()


def dh(z):
    if z <= 0 or z >= 1:
        return Decimal(0)
    return -Decimal(z.numerator) / Decimal(z.denominator) * dlog2(z) \
        - Decimal((1 - z).numerator) / Decimal((1 - z).denominator) \
        * dlog2(1 - z)


def exact_hu_cr(n, mu):
    """CR_HU in exact Fractions for the cell tree, evaluated to 60
    decimal digits (Decimal.ln); an algorithm independent of both
    certified-log2 kits."""
    pref = []
    for i in range(n):
        mask = (1 << i) - 1
        d = {}
        for A, w in mu.items():
            e = d.setdefault(A & mask, [Fraction(0), Fraction(0)])
            e[0] += w
            if not (A >> i) & 1:
                e[1] += w
        pref.append(d)
    cells = {(0, 0): Fraction(1)}
    S = Decimal(0)
    for i in range(n):
        nxt = {}
        for (a, b), w in cells.items():
            da, db = pref[i][a], pref[i][b]
            x, y = da[1] / da[0], db[1] / db[0]
            z = min(max(Fraction(1, 2), x + y - 1), x, y)
            S += Decimal(w.numerator) / Decimal(w.denominator) * dh(z)
            for (ai, bi), c in (((0, 0), z), ((0, 1), x - z),
                                ((1, 0), y - z),
                                ((1, 1), 1 - x - y + z)):
                if c <= 0:
                    continue
                k = (a | (ai << i), b | (bi << i))
                nxt[k] = nxt.get(k, Fraction(0)) + w * c
        cells = nxt
    H = Decimal(0)
    for w in mu.values():
        H -= Decimal(w.numerator) / Decimal(w.denominator) * dlog2(w)
    return S - H, H


def rationalize(mu):
    m = {a: Fraction(w).limit_denominator(10 ** 7) for a, w in mu.items()}
    t = sum(m.values())
    return {a: w / t for a, w in m.items()}


def part5():
    say("P5. Certificates recomputed at 60 digits (Decimal ln), against "
        "data/hu_certify.json")
    cert = json.loads((DATA / "hu_certify.json").read_text())
    byname = {r["instance"]: r for r in cert["instances"]}
    for name, n, mu, _ in floor_instances():
        if name not in ("windowkill", "mmabskill"):
            continue
        rat = rationalize(mu)
        tv = sum(abs(float(rat[a]) - mu[a]) for a in mu) / 2
        cr, H = exact_hu_cr(n, rat)
        ok = True
        for kit, kd in byname[name]["kits"].items():
            ok &= Decimal(repr(kd["CR_lo"])) - Decimal("1e-9") <= cr \
                <= Decimal(repr(kd["CR_hi"])) + Decimal("1e-9")
            ok &= kd["certifies"]
        ok &= cr > H / 25 and cr > 0
        check(f"{name}: CR_HU inside both kits' enclosures, > H/25",
              ok, f"CR={cr:.12f} H={H:.9f} H/25={H / 25:.9f}, "
              f"rationalization TV={tv:.2e}")
    row, mu = witness_from_checkpoint()
    n = row["n"]
    mup = permute(n, mu, (0, 1, 3, 2))
    rat = rationalize(mup)
    tv = sum(abs(float(rat[a]) - mup[a]) for a in mup) / 2
    cr, H = exact_hu_cr(n, rat)
    kd = byname["033-witness-order-0132"]["kits"]
    ok = cr < 0 and all(k["certifies"] and k["CR_hi"] < 0
                        for k in kd.values())
    ok &= all(Decimal(repr(k["CR_lo"])) - Decimal("1e-9") <= cr
              <= Decimal(repr(k["CR_hi"])) + Decimal("1e-9")
              for k in kd.values())
    mm = max(sum(w for a, w in rat.items() if (a >> i) & 1)
             for i in range(n))
    check("033 witness: CR_HU < 0 at 60 digits, inside both enclosures, "
          "marginals < 9/20 exactly",
          ok and mm < Fraction(9, 20),
          f"CR={cr:.12f}, kit hi={kd['A_digit_extraction']['CR_hi']:.9e}, "
          f"maxmarg={float(mm):.6f}, TV={tv:.2e}")
    say()


def main():
    say("Attempt 034 reviewer pass: independent re-implementation "
        "(030 / 031 / 033)")
    say("=" * 70)
    part1()
    part2()
    part3()
    part4()
    part5()
    say("=" * 70)
    if FAILS:
        say(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            say(f"  - {f}")
    else:
        say("No refutation found.")
    (DATA / "uc_reviewer034_reimpl_out.txt").write_text("\n".join(LINES) + "\n")
    raise SystemExit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
