#!/usr/bin/env python3
"""Attempt 017: skeptic re-computation of 016's aggregated-control kill.

INDEPENDENCE CONTRACT.  This file imports nothing from uc_or_avg.py,
uc_or_avg_skeptic.py, uc_agg_ctrl_probe.py or uc_agg_ctrl_probe2.py.  The
quantity is re-derived from the PROSE of 007 section 1 and 007 section 5 /
013 part C:

    tilt coupling  pi(A,B) prop u(A)u(B) t^{|A&B|},  t = 2^lambda, both
    marginals mu (any positive u is its own marginal's potential, 005);
    for coordinate i (1..n; coordinate j is bit j-1), histories are the
    (i-1)-bit prefixes; N_i = prefixes of supp(mu) with BOTH continuations
    reachable; for (a,b) in N_i^2 the conditional 2x2 table over the
    response bits has four positive cells and OR_i(a,b) is its odds ratio;
    m_ab = pi-mass of the history pair; M_i = sum m_ab log2 OR / sum m_ab;
    w_i = sum_{N_i^2} m_ab (defined history mass);

    A(mu, lambda) = sum_{i<n, N_i nonempty} w_i (M_i - lambda) / sum w_i.

At rational t the sign of A is the sign of  sum_i sum_rows m * log2(OR/t)
with every OR/t exactly rational, certified here by MY OWN log2 enclosure:
repeated interval squaring with directed dyadic rounding (bit-extraction),
NOT 013's atanh series.  ln 2 bounds come from the exact series
ln 2 = sum_{k>=1} 1/(k 2^k) with its elementary geometric tail bound.

All couplings, tables, odds ratios and marginals are computed in scaled
INTEGER arithmetic (every input weight is rational; multiply through by the
lcm of denominators and by den(t)^n), so no rounding exists anywhere except
inside the log2 enclosure, which is directed.

Standard library only.  Deterministic.  Instances are read as DATA from
016's checkpoints (theta parameters / atom lists); the family construction
MU(n, r) is re-implemented from the prose of 016 section 1 and
cross-checked atom-by-atom in scripts/ -- not imported.

Usage:
    python problems/union-closed/explore/uc_agg_ctrl_skeptic.py [--quick]
        [--n128] [--n160]
Checkpoints: problems/union-closed/data/aggsk15_*.json
"""

from __future__ import annotations

import json
import math
import sys
import time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / "data"

RECORD_NUM, RECORD_DEN = 38271, 100000        # threshold 0.38271 (Liu)
F0 = Fraction(0)

# ======================================================================
# certified log2 enclosure -- interval squaring with directed rounding
# ======================================================================

# ln 2 = sum_{k>=1} 1/(k 2^k); partial sum is a strict lower bound and the
# tail after K terms is < sum_{k>K} 1/((K+1) 2^k) = 1/((K+1) 2^K).
_LN2_TERMS = 64
LN2_LO = sum(Fraction(1, k * (1 << k)) for k in range(1, _LN2_TERMS + 1))
LN2_HI = LN2_LO + Fraction(1, (_LN2_TERMS + 1) * (1 << _LN2_TERMS))
assert LN2_LO < LN2_HI


def log2_enclosure(q: Fraction, K: int = 64, prec: int = 192
                   ) -> tuple[Fraction, Fraction]:
    """[lo, hi] with lo <= log2(q) <= hi, q > 0 rational.

    Method: write q = 2^e * x, x in [1, 2).  Repeat K times: square the
    dyadic interval enclosing x (lower endpoint rounded down, upper rounded
    up, at 2^-prec resolution) and renormalize by powers of two, tracking
    the accumulated exponent S = sum_k s_k 2^{K-k}.  Then
        log2 x = (S + log2 x_K) / 2^K
    and x_K lies in a known tiny interval well inside (1/4, 4), where the
    elementary bounds 1 - 1/v <= ln v <= v - 1 (all v > 0) and the LN2
    bounds above give a crude but certified enclosure of log2 x_K whose
    error is divided by 2^K.  Total width ~ 2^-60 or better.
    """
    num, den = q.numerator, q.denominator
    if num <= 0:
        raise ValueError("log2 of nonpositive")
    if num & (num - 1) == 0 and den & (den - 1) == 0:
        e = (num.bit_length() - 1) - (den.bit_length() - 1)
        return Fraction(e), Fraction(e)
    e = num.bit_length() - den.bit_length()
    if e >= 0:
        nn, dd = num, den << e
    else:
        nn, dd = num << (-e), den
    if nn < dd:                          # x < 1: double it
        nn <<= 1
        e -= 1
    # now dd <= nn < 2 dd, i.e. x = nn/dd in [1, 2)
    scale = 1 << prec
    two_scale = scale << 1
    xlo = (nn << prec) // dd             # floor
    xhi = -((-(nn << prec)) // dd)       # ceil
    S = 0
    for _ in range(K):
        xlo = (xlo * xlo) >> prec                    # floor(x^2 * scale)
        xhi = -((-(xhi * xhi)) >> prec)              # ceil
        s = 0
        while xlo >= two_scale:                      # whole interval >= 2
            xlo >>= 1
            xhi = -((-xhi) >> 1)
            s += 1
        while xhi < scale:                           # whole interval < 1
            xlo <<= 1
            xhi <<= 1
            s -= 1
        S = (S << 1) + s
    vlo = Fraction(xlo, scale)
    vhi = Fraction(xhi, scale)
    assert 0 < vlo <= vhi
    # certified crude bounds on log2 v over v in [vlo, vhi]
    ln_lo = 1 - Fraction(1) / vlo        # <= ln vlo
    ln_hi = vhi - 1                      # >= ln vhi
    l2_lo = ln_lo / (LN2_HI if ln_lo >= 0 else LN2_LO)
    l2_hi = ln_hi / (LN2_LO if ln_hi >= 0 else LN2_HI)
    pw = 1 << K
    return e + Fraction(S + l2_lo, 1) / pw, e + Fraction(S + l2_hi, 1) / pw


def _selftest_log2() -> dict:
    """Certified-enclosure self-test on known values (attack item g)."""
    out = {"cases": [], "ok": True}
    cases = [
        (Fraction(1), 0.0), (Fraction(2), 1.0), (Fraction(1024), 10.0),
        (Fraction(1, 8), -3.0), (Fraction(3), math.log2(3)),
        (Fraction(10), math.log2(10)), (Fraction(1, 3), -math.log2(3)),
        (Fraction(181, 16), math.log2(181 / 16)),
        (Fraction(7, 5) ** 10, 10 * math.log2(1.4)),
        (Fraction(10 ** 40 + 1, 10 ** 40), math.log2(1 + 1e-40)),
    ]
    for q, ref in cases:
        lo, hi = log2_enclosure(q)
        w = float(hi - lo)
        ok = lo <= hi and w < 1e-17 and lo - 1e-12 <= ref <= hi + 1e-12
        out["cases"].append({"q": str(q)[:60], "ref": ref, "lo": float(lo),
                             "hi": float(hi), "width": w, "ok": ok})
        out["ok"] = out["ok"] and ok
    # additivity consistency: log2(ab) inside sum of enclosures
    a, b = Fraction(355, 113), Fraction(1234567, 7654321)
    la, ha = log2_enclosure(a)
    lb, hb = log2_enclosure(b)
    lab, hab = log2_enclosure(a * b)
    add_ok = la + lb <= hab and lab <= ha + hb
    out["additivity_ok"] = bool(add_ok)
    out["ok"] = out["ok"] and add_ok
    # exactness on powers of two
    lo, hi = log2_enclosure(Fraction(1, 2 ** 100))
    out["pow2_exact"] = bool(lo == hi == -100)
    out["ok"] = out["ok"] and out["pow2_exact"]
    return out


# ======================================================================
# exact aggregate (scaled-integer census)
# ======================================================================

def exact_aggregate(u: dict[int, Fraction], n: int, t: Fraction,
                    progress: bool = False) -> dict:
    """Certified enclosure of A(mu, lambda) at rational tilt t for the
    measure whose tilt potential is u (lambda = log2 t).  Everything up to
    the log2 calls is exact integer arithmetic."""
    t_start = time.monotonic()
    supp = sorted(u)
    assert all(w > 0 for w in u.values()), "zero/negative weight atom"
    # scale weights to integers
    L = 1
    for w in u.values():
        L = L * w.denominator // math.gcd(L, w.denominator)
    U = {a: int(u[a] * L) for a in supp}
    tn, td = t.numerator, t.denominator
    # kernel k(A,B) = U_A U_B tn^p td^(n-p), p = |A & B|  (global scale
    # L^2 td^n; constant across pairs and coordinates, cancels everywhere)
    pw_tn = [1] * (n + 1)
    pw_td = [1] * (n + 1)
    for j in range(1, n + 1):
        pw_tn[j] = pw_tn[j - 1] * tn
        pw_td[j] = pw_td[j - 1] * td
    m_atoms = len(supp)
    kern_rows: list[list[int]] = []
    for A in supp:
        UA = U[A]
        row = []
        for B in supp:
            p = bin(A & B).count("1")
            row.append(UA * U[B] * pw_tn[p] * pw_td[n - p])
        kern_rows.append(row)
    Z = sum(sum(r) for r in kern_rows)          # total coupling mass, scaled

    # exact elementwise marginals of mu: mu(A) = sum_B k(A,B)
    mu_row = [sum(r) for r in kern_rows]
    marg_num = [0] * n
    for idx, A in enumerate(supp):
        v = mu_row[idx]
        for j in range(n):
            if A >> j & 1:
                marg_num[j] += v
    in_regime = all(mn * RECORD_DEN < RECORD_NUM * Z for mn in marg_num)
    jmax = max(range(n), key=lambda j: marg_num[j])
    max_marg_float = marg_num[jmax] / Z
    # entropy of mu (float, reporting only)
    H = -sum((v / Z) * math.log2(v / Z) for v in mu_row if v)

    # census
    num_terms: dict[Fraction, int] = {}     # q = OR/t  ->  sum of masses
    W_tot = 0                               # sum_{i<n} w_i, scaled
    w_last = 0                              # w_n (for the inclusion variant)
    dich_ok = True
    n_rows = 0
    per_i = []
    for i in range(1, n + 1):
        pmask = (1 << (i - 1)) - 1
        bit = 1 << (i - 1)
        has0: dict[int, bool] = {}
        has1: dict[int, bool] = {}
        pref = []
        resp = []
        for A in supp:
            c = A & pmask
            r = 1 if A & bit else 0
            pref.append(c)
            resp.append(r)
            if r:
                has1[c] = True
            else:
                has0[c] = True
        active = {c for c in has0 if c in has1}
        tables: dict[tuple[int, int], list[int]] = {}
        for ia in range(m_atoms):
            ca, ra, row = pref[ia], resp[ia], kern_rows[ia]
            for ib in range(m_atoms):
                key = (ca, pref[ib])
                tab = tables.get(key)
                if tab is None:
                    tab = tables[key] = [0, 0, 0, 0]
                tab[(ra << 1) | resp[ib]] += row[ib]
        w_i = 0
        rows_i = []
        for (a, b), (t00, t01, t10, t11) in tables.items():
            nondeg = t00 > 0 and t01 > 0 and t10 > 0 and t11 > 0
            if nondeg != (a in active and b in active):
                dich_ok = False
            if nondeg:
                m = t00 + t01 + t10 + t11
                q = Fraction(t00 * t11, t01 * t10) / t
                rows_i.append((m, q))
                w_i += m
        if i < n:
            for m, q in rows_i:
                num_terms[q] = num_terms.get(q, 0) + m
            W_tot += w_i
            n_rows += len(rows_i)
        else:
            w_last = w_i
            # free anchor (005 Prop 1): OR = t exactly at every
            # nondegenerate history at i = n, i.e. q = OR/t = 1
            prop1_ok = all(q == 1 for _m, q in rows_i)
        per_i.append({"i": i, "n_active": len(active),
                      "n_rows": len(rows_i), "w_i_over_Z": w_i / Z if Z else 0})
        if progress and i % 16 == 0:
            print(f"    ... census i={i}/{n} "
                  f"({time.monotonic() - t_start:.0f}s)", flush=True)

    # certified numerator enclosure (deduplicated log2 calls)
    N_lo = N_hi = F0
    for q, m in num_terms.items():
        lo, hi = log2_enclosure(q)
        N_lo += m * lo
        N_hi += m * hi
    agg_lo = N_lo / W_tot
    agg_hi = N_hi / W_tot
    # alternative bookkeeping (013): every degenerate positive-mass history
    # scored at exactly lambda -> per-i defined mass becomes the full
    # history mass Z, numerator unchanged: A_alt = N / ((n-1) Z).
    alt_lo = N_lo / ((n - 1) * Z)
    alt_hi = N_hi / ((n - 1) * Z)
    # i = n inclusion variant: M_n = lambda identically, so only the
    # denominator grows.
    incl_lo = N_lo / (W_tot + w_last)
    incl_hi = N_hi / (W_tot + w_last)
    # union-closedness of the support (informational: the control statement
    # quantifies over ALL mu, closure is not a hypothesis -- checked so the
    # record can say so explicitly)
    sset = set(supp)
    uc = all((A | B) in sset for A in supp for B in supp)
    return {
        "n": n, "k_atoms": m_atoms, "t": t,
        "agg_lo": agg_lo, "agg_hi": agg_hi,
        "alt_lo": alt_lo, "alt_hi": alt_hi,
        "incl_lo": incl_lo, "incl_hi": incl_hi,
        "width": agg_hi - agg_lo,
        "W_over_Z": W_tot / Z, "dich_ok": dich_ok, "prop1_ok": prop1_ok,
        "in_regime": in_regime, "max_marg_float": max_marg_float,
        "argmax_marg_coord": jmax + 1, "H_bits": H,
        "n_log_terms": len(num_terms), "n_rows": n_rows,
        "union_closed_support": uc,
        "per_i": per_i,
        "secs": round(time.monotonic() - t_start, 1),
    }


# ======================================================================
# instance builders (re-implemented from 016's prose; data from
# checkpoints)
# ======================================================================

def mu_ladder_from_theta(n: int, r: int, th: dict[str, float]
                         ) -> dict[int, Fraction]:
    """MU(n, r) from the prose of 016 section 1.  Coordinate j = bit j-1.
    Block a: {1,n-1}:A1, {1,n}:A2, light slices {1,c_j,n}:light;
    block b: {}:B0, {n-1}:B1, {n}:B2, responders {c_j,n-1}:resp;
    dilution {2},{3},{4}:dil; response coordinates c_j = 5..4+r.
    Weights taken as EXACT dyadic rationals of the checkpoint floats."""
    assert 1 <= r <= n - 6
    hi = 1 << (n - 1)          # coordinate n
    lo = 1 << (n - 2)          # coordinate n-1
    F = Fraction
    u = {
        1 | lo: F(th["A1"]), 1 | hi: F(th["A2"]), 0: F(th["B0"]),
        lo: F(th["B1"]), hi: F(th["B2"]),
    }
    for j in range(r):
        c = 1 << (4 + j)       # coordinate 5+j
        u[1 | c | hi] = F(th["light"])
        u[c | lo] = F(th["resp"])
    for d in (1, 2, 3):        # coordinates 2, 3, 4
        u[1 << d] = F(th["dil"])
    assert len(u) == 2 * r + 8
    return u


def tidy(u: dict[int, Fraction], den: int = 10 ** 4) -> dict[int, Fraction]:
    """016's tidy-rational convention: limit_denominator(1e4), zeros
    dropped after rationalization."""
    out = {a: w.limit_denominator(den) for a, w in u.items()}
    return {a: w for a, w in out.items() if w > 0}


# 007's 10-atom witness (weights are DATA copied from 007's record /
# checkpoint; n = 7, response coordinate 5, futures 6 and 7)
WITNESS_THETA = {
    "A1": 3.802705744653898, "A2": 8.094727225702503,
    "light": 0.004792968160237178, "B0": 0.2526660836353666,
    "B1": 27.28641199464732, "B2": 12.529780963546393,
    "resp": 0.3609234435898425, "dil": 20.0,
}


def fmt(r: dict) -> str:
    return (f"A in [{float(r['agg_lo']):+.12f}, {float(r['agg_hi']):+.12f}] "
            f"(width {float(r['width']):.1e}); "
            f"{'KILL' if r['agg_hi'] < 0 else 'POSITIVE' if r['agg_lo'] > 0 else 'STRADDLE'}; "
            f"in-regime {r['in_regime']} (max marg {r['max_marg_float']:.5f} "
            f"at coord {r['argmax_marg_coord']}); dich {r['dich_ok']}; "
            f"prop1 {r['prop1_ok']}; "
            f"H {r['H_bits']:.3f}; W/Z {float(r['W_over_Z']):.4f}; "
            f"alt-conv [{float(r['alt_lo']):+.6f}, {float(r['alt_hi']):+.6f}]; "
            f"union-closed {r['union_closed_support']}; {r['secs']}s")


def save(name: str, obj) -> None:
    DATA.mkdir(exist_ok=True)

    def default(o):
        if isinstance(o, Fraction):
            return str(o) if o.denominator < 10 ** 40 else repr(float(o))
        return str(o)
    (DATA / name).write_text(json.dumps(obj, indent=1, sort_keys=True,
                                        default=default) + "\n",
                             encoding="utf-8")


def public(r: dict) -> dict:
    keep = ("n", "k_atoms", "in_regime", "max_marg_float",
            "argmax_marg_coord", "H_bits", "dich_ok", "prop1_ok",
            "n_log_terms",
            "n_rows", "union_closed_support", "secs")
    out = {k: r[k] for k in keep}
    out["t"] = str(r["t"])
    for k in ("agg_lo", "agg_hi", "alt_lo", "alt_hi", "incl_lo", "incl_hi",
              "width", "W_over_Z"):
        out[k + "_float"] = float(r[k])
    out["verdict"] = ("CERTIFIED NEGATIVE (KILL)" if r["agg_hi"] < 0 else
                      "CERTIFIED POSITIVE" if r["agg_lo"] > 0 else
                      "STRADDLES 0")
    return out


def main() -> None:
    quick = "--quick" in sys.argv
    t0 = time.monotonic()

    print("[S0] log2 enclosure self-test (attack g)")
    st = _selftest_log2()
    for c in st["cases"]:
        print(f"    q={c['q'][:40]:40s} ref={c['ref']:+.12f} "
              f"[{c['lo']:+.12f},{c['hi']:+.12f}] w={c['width']:.1e} "
              f"{'ok' if c['ok'] else 'FAIL'}")
    print(f"    additivity {st['additivity_ok']}, pow2 exact "
          f"{st['pow2_exact']}, ALL {'OK' if st['ok'] else 'FAIL'}")
    save("aggsk15_selftest.json", st)
    assert st["ok"], "log2 enclosure self-test failed"

    results = {}

    # ---- anchor: 007's witness at t = 181/16 (013 part C: +1.844669005)
    uw = mu_ladder_from_theta(7, 1, WITNESS_THETA)
    r = exact_aggregate(uw, 7, Fraction(181, 16))
    print("[S1] witness n=7 t=181/16 (expect +1.844669005):")
    print("    " + fmt(r))
    results["witness_t181_16"] = public(r)

    # ---- the load-bearing kill: theta*(2.0) ladder at n = 96, t = 4
    d = json.loads((DATA / "aggprobe2_partP2.json").read_text())
    th96 = [row for row in d["rows"] if row["stage"] == "transfer"
            and row["lambda"] == 2.0 and row["n"] == 96][0]["theta"]
    u96 = mu_ladder_from_theta(96, 90, th96)
    print(f"[S2] theta* ladder n=96 r=90 ({len(u96)} atoms) t=4 "
          f"(016 cert: [-0.000759865185844166087, ...086]):")
    r96 = exact_aggregate(u96, 96, Fraction(4), progress=True)
    print("    " + fmt(r96))
    results["theta2_n96_t4"] = public(r96)
    save("aggsk15_results.json", results)

    # ---- positive control at n = 32: P3-best, tidy zero-dropped, t = 4
    if not quick:
        d3 = json.loads((DATA / "aggprobe2_partP3.json").read_text())
        b = d3["best"]
        u32 = tidy({int(k, 2): Fraction(v) for k, v in b["u"].items()})
        print(f"[S3] P3-best n={b['n']} tidy ({len(u32)} atoms) t=4 "
              f"(016 cert: +0.297679573):")
        r32 = exact_aggregate(u32, b["n"], Fraction(4))
        print("    " + fmt(r32))
        results["p3best_n32_tidy_t4"] = public(r32)
        save("aggsk15_results.json", results)

    # ---- surviving-positive-space check: RAW witness-weight ladder at
    #      n = 96 with 016's P6 adversarial dilution weight, t = 4
    if not quick:
        d6 = json.loads((DATA / "aggprobe2_partP6.json").read_text())
        row6 = [x for x in d6["rows"] if x["n"] == 96
                and x["lambda"] == 2.0][0]
        thraw = dict(WITNESS_THETA)
        thraw["dil"] = row6["wdil"]
        uraw = mu_ladder_from_theta(96, 90, thraw)
        print(f"[S4] RAW witness ladder n=96 wdil={row6['wdil']:.6f} t=4 "
              f"(016 P6 float: {row6['agg']:+.6f}):")
        rr = exact_aggregate(uraw, 96, Fraction(4), progress=True)
        print("    " + fmt(rr))
        results["raw_ladder_n96_t4"] = public(rr)
        save("aggsk15_results.json", results)

    # ---- optional deeper kills
    if "--n128" in sys.argv:
        th128 = [row for row in d["rows"] if row["stage"] == "transfer"
                 and row["lambda"] == 2.0 and row["n"] == 128][0]["theta"]
        u128 = tidy(mu_ladder_from_theta(128, 122, th128))
        print(f"[S5] theta* ladder n=128 tidy ({len(u128)} atoms) t=4 "
              f"(016 cert: -0.015405727):")
        r128 = exact_aggregate(u128, 128, Fraction(4), progress=True)
        print("    " + fmt(r128))
        results["theta2_n128_tidy_t4"] = public(r128)
        save("aggsk15_results.json", results)
    if "--n160" in sys.argv:
        d7 = json.loads((DATA / "aggprobe2_partP7.json").read_text())
        u160 = {int(k, 2): Fraction(float(v))
                for k, v in d7["u_target"].items()}
        print(f"[S6] theta* ladder n=160 dyadic ({len(u160)} atoms) t=4 "
              f"(016 cert: -0.024225909):")
        r160 = exact_aggregate(u160, 160, Fraction(4), progress=True)
        print("    " + fmt(r160))
        results["theta2_n160_dyadic_t4"] = public(r160)
        save("aggsk15_results.json", results)

    print(f"done in {time.monotonic() - t0:.0f}s")


if __name__ == "__main__":
    main()
