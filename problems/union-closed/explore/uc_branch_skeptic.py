#!/usr/bin/env python3
"""Attempt 025 skeptic pass: independent re-implementation of every
load-bearing number in uc_branch_attack.py.  Default stance: refute.

Deliberately different code paths:
  - natural-log arithmetic throughout, converted to bits only on output;
  - entropies by direct mask enumeration over explicit measures (no
    closed forms, no profile arithmetic);
  - chain-rule CR by a B-major nested-dict grouping (vs uc_mitax's
    (a,b)-tuple keyed flat dict);
  - couplings built from 004's *verbal* description ("with prob 2*P1,
    (A,B) iid uniform on {empty, s}; else A = B = s"), not from the
    q-parametrized builder under test;
  - rescue existence checked with the plain iid coupling mu x mu
    (no Sinkhorn dependency).

Shares NOTHING with uc_branch_attack.py / uc_mitax.py but the stdlib.
Reads data/branch_attack.json only to compare against its claims.

Usage: python uc_branch_skeptic.py
Output: ../data/branch_skeptic_out.txt (tee of stdout), plus exit 1 on
any refutation.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
LN2 = math.log(2.0)
MARG = 0.38271
FAILS = []
LINES = []


def say(m=""):
    print(m, flush=True)
    LINES.append(m)


def check(name, ok, detail=""):
    say(f"  [{'ok' if ok else 'REFUTED'}] {name}" + (f"  ({detail})" if detail
                                                     else ""))
    if not ok:
        FAILS.append(name)


def hnat(p):
    """binary entropy of p, in NATS."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log(p) - (1.0 - p) * math.log(1.0 - p)


def H_of(measure):
    """Shannon entropy (nats) of a dict mass function."""
    return -sum(w * math.log(w) for w in measure.values() if w > 0.0)


def bern_measure(n, p):
    out = {}
    for a in range(1 << n):
        k = bin(a).count("1")
        w = (p ** k) * ((1.0 - p) ** (n - k))
        if w > 0.0:
            out[a] = w
    return out


def mix(measures_weights):
    out = {}
    for meas, wt in measures_weights:
        for a, w in meas.items():
            out[a] = out.get(a, 0.0) + wt * w
    return out


def cr_chain(n, pairs):
    """CR, gain, ST (bits) by an independent accounting: B-major nested
    dicts, natural logs, coordinate i revealed in increasing bit order
    (the convention of record)."""
    A_marg, U_law = {}, {}
    for (A, B), w in pairs.items():
        A_marg[A] = A_marg.get(A, 0.0) + w
        U_law[A | B] = U_law.get(A | B, 0.0) + w
    Hmu = H_of(A_marg)
    HU = H_of(U_law)
    sum_hz = 0.0
    for i in range(n):
        lowmask = (1 << i) - 1
        by_b = {}
        for (A, B), w in pairs.items():
            d = by_b.setdefault(B & lowmask, {})
            g = d.setdefault(A & lowmask, [0.0, 0.0])
            g[0] += w
            if not ((A >> i) & 1 or (B >> i) & 1):
                g[1] += w
        for d in by_b.values():
            for m_, z_ in d.values():
                if m_ > 0.0:
                    zz = z_ / m_
                    if 0.0 < zz < 1.0:
                        sum_hz += m_ * (-zz * math.log(zz)
                                        - (1 - zz) * math.log(1 - zz))
    CR = (sum_hz - Hmu) / LN2
    gain = (HU - Hmu) / LN2
    return {"CR": CR, "gain": gain, "ST": gain - CR, "H_mu": Hmu / LN2}


def halfmix_pairs_004(n, ph, P1):
    """004's verbal recipe, built from the (branch, s) enumeration:
    draw s; with prob 2*P1 pick (A,B) iid uniform on {0, s}; else A=B=s."""
    pairs = {}
    for s, ps in bern_measure(n, ph).items():
        for A in (0, s):
            for B in (0, s):
                w = ps * 2.0 * P1 * 0.25
                pairs[(A, B)] = pairs.get((A, B), 0.0) + w
        pairs[(s, s)] = pairs.get((s, s), 0.0) + ps * (1.0 - 2.0 * P1)
    return pairs


def eps_mixing_pairs(n, smeasure, P1, q):
    """Generalized empty-mixing coupling by explicit (eps_A, eps_B, s)
    triple enumeration."""
    eps_law = {(0, 0): q, (0, 1): P1 - q, (1, 0): P1 - q,
               (1, 1): 1.0 - 2.0 * P1 + q}
    assert min(eps_law.values()) > -1e-12
    pairs = {}
    for (ea, eb), we in eps_law.items():
        if we <= 0.0:
            continue
        for s, ps in smeasure.items():
            k = (s if ea else 0, s if eb else 0)
            pairs[k] = pairs.get(k, 0.0) + we * ps
    return pairs


def main():
    say("Skeptic pass on uc_branch_attack.py (independent implementation)")
    say("=" * 70)
    ref = json.loads((DATA / "branch_attack.json").read_text())

    # ---------------------------------------------------------------- S1
    say("\nS1. Headline half-mixing refutation instances (004 verbal build):")
    ph16 = next(r["ph"] for r in ref["A2_family"] if r["n"] == 16)
    for n, ph, P1 in [(2, 0.9935, 0.64), (6, 0.9985, 0.64),
                      (16, ph16, 0.64)]:
        pairs = halfmix_pairs_004(n, ph, P1)
        r = cr_chain(n, pairs)
        mu = mix([(bern_measure(n, ph), 1.0 - P1), ({0: 1.0}, P1)])
        marg = max(sum(w for a, w in mu.items() if (a >> i) & 1)
                   for i in range(n))
        Hmu_bits = H_of(mu) / LN2
        check(f"CR_halfmix(n={n}, ph={ph}, P1={P1}) < 0 in-regime",
              r["CR"] < 0 and marg < MARG and Hmu_bits > 0.9,
              f"CR={r['CR']:+.6f}, maxmarg={marg:.4f}, H={Hmu_bits:.4f}")
        # against the engine's A2/A5 records where available
        for row in ref.get("A2_family", []) + [
                {"n": x["n"], "ph": x["ph"], "P1": x["P1"],
                 "CR_halfmix": x["CR_halfmix"]} for x in ref["A5_rescue"]]:
            if row["n"] == n and abs(row["ph"] - ph) < 5e-7 \
                    and abs(row["P1"] - P1) < 1e-12:
                check(f"  matches engine value at n={n}",
                      abs(row["CR_halfmix"] - r["CR"]) < 1e-9,
                      f"delta={row['CR_halfmix'] - r['CR']:.2e}")

    # ---------------------------------------------------------------- S2
    say("\nS2. Extended-HM lemma (ST = 0) on independently built eps-mixing")
    say("    couplings, incl. a mixture s-law; and iid-coupling rescue:")
    worst_st = 0.0
    for n, smeas, P1, q in [
            (4, bern_measure(4, 0.997), 0.64, 0.5),
            (5, bern_measure(5, 0.6), 0.5, 0.17),
            (6, mix([(bern_measure(6, 0.90), 0.75),
                     (bern_measure(6, 0.99), 0.25)]), 0.60, 0.2),
            (3, mix([(bern_measure(3, 0.55), 0.5),
                     (bern_measure(3, 0.95), 0.5)]), 0.30, 0.11)]:
        r = cr_chain(n, eps_mixing_pairs(n, smeas, P1, q))
        worst_st = max(worst_st, abs(r["ST"]))
    check("ST = 0 for arbitrary (eps,s) mixing couplings",
          worst_st < 1e-10, f"max |ST| = {worst_st:.2e}")

    for n, ph, P1 in [(2, 0.9935, 0.64), (6, 0.9985, 0.64)]:
        mu = mix([(bern_measure(n, ph), 1.0 - P1), ({0: 1.0}, P1)])
        pairs = {(a, b): wa * wb for a, wa in mu.items()
                 for b, wb in mu.items()}
        r = cr_chain(n, pairs)
        check(f"iid coupling rescues (n={n}) without Sinkhorn",
              r["CR"] > 0, f"CR_iid={r['CR']:+.5f}")

    # 3block empty-mixing rescue at the engine's q*
    smeas = mix([(bern_measure(6, 0.90), 0.75),
                 (bern_measure(6, 0.99), 0.25)])
    q_star = ref["B3_rescue"]["empty_mixing_q"]
    r = cr_chain(6, eps_mixing_pairs(6, smeas, 0.60, q_star))
    check("3block empty-mixing rescue value reproduces",
          abs(r["CR"] - ref["B3_rescue"]["empty_mixing_CR"]) < 1e-9
          and r["CR"] > 0.5,
          f"CR={r['CR']:+.6f} vs engine "
          f"{ref['B3_rescue']['empty_mixing_CR']:+.6f}")

    # ---------------------------------------------------------------- S3
    say("\nS3. Block-branch gains by direct mask enumeration (no profile")
    say("    arithmetic), n = 8 and 12:")
    insts = {"2block": [(0.62, 0.0), (0.38, 0.95)],
             "3block": [(0.60, 0.0), (0.30, 0.90), (0.10, 0.99)],
             "hi+lo": [(0.50, 0.10), (0.50, 0.60)]}
    for tag, comps in insts.items():
        engine_rows = {r["n"]: r for r in
                       next(x for x in ref["B1_gain"] if x["tag"] == tag)
                       ["rows"]}
        for n in (8, 12):
            mu = mix([(bern_measure(n, p), P) for P, p in comps])
            uu = mix([(bern_measure(n, 1.0 - (1.0 - p) ** 2), P)
                      for P, p in comps])
            gain = (H_of(uu) - H_of(mu)) / LN2
            ok = abs(gain - engine_rows[n]["gain_block"]) < 1e-8
            check(f"{tag} gain(n={n}) reproduces and is negative"
                  if gain < 0 else f"{tag} gain(n={n}) reproduces",
                  ok and (gain < 0 or n == 2),
                  f"gain={gain:+.5f} vs engine "
                  f"{engine_rows[n]['gain_block']:+.5f}")

    # comonotone coupling CR = 0 identity, independent build
    n = 6
    mu = mix([(bern_measure(n, p), P) for P, p in insts["3block"]])
    r = cr_chain(n, {(a, a): w for a, w in mu.items()})
    check("comonotone CR = 0 exactly", abs(r["CR"]) < 1e-12,
          f"CR={r['CR']:+.2e}")

    # ---------------------------------------------------------------- S4
    say("\nS4. Coverage-lemma boundary: sign of F'(P1) by central")
    say("    difference of brute-force entropy, n in {1, 2}:")
    for n in (1, 2):
        bad = 0
        pts = 1001
        eps = 1e-6
        for j in range(pts):
            ph = 0.5 + (2 * MARG - 0.5) * j / (pts - 1)
            P1 = max(1e-9, 1.0 - MARG / ph + 1e-9)
            s = bern_measure(n, ph)

            def F(w):
                return H_of(mix([({0: 1.0}, w), (s, 1.0 - w)])) / LN2

            d = (F(P1 + eps) - F(P1 - eps)) / (2 * eps)
            if d >= 0.0:
                bad += 1
        check(f"F'(P1) < 0 on the whole n={n} boundary grid", bad == 0,
              f"{pts - bad}/{pts} negative")
    # analytic constants of case (3)
    h2m = hnat(2 * MARG) / LN2
    r_min = (1 - 2 * MARG) / (2 * MARG)
    xmax = max(-x * math.log(x) / LN2 for x in
               [j / 100000.0 for j in range(1, 100000)])
    margin3 = 3 * h2m - xmax - (-math.log(r_min) / LN2)
    check("case-(3) constants and n>=3 margin > 0",
          abs(h2m - 0.785910) < 1e-5 and margin3 > 0.12,
          f"h(2M)={h2m:.6f}, max(-xlogx)={xmax:.6f}, margin={margin3:+.4f}")

    # ---------------------------------------------------------------- S5
    say("\nS5. psi sign flip: root of h((1-p)^2) = h(p):")
    lo, hi = 0.3, 0.45
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if hnat((1.0 - mid) ** 2) - hnat(mid) > 0:
            lo = mid
        else:
            hi = mid
    root = 0.5 * (lo + hi)
    psi = (3.0 - math.sqrt(5.0)) / 2.0
    check("root equals (3-sqrt5)/2", abs(root - psi) < 1e-12,
          f"root={root:.15f}")

    # ---------------------------------------------------------------- S6
    say("\nS6. Floor-ratio table: recompute CR/H from cr_attack.json:")
    d = json.loads((DATA / "cr_attack.json").read_text())
    worst = None
    for row in d["rows"]:
        if row["H_mu"] > 0.5:  # in-regime, H-bounded endpoints
            ratio = row["sup_cr"] / row["H_mu"]
            eng = next(r for r in ref["C_ratios"]
                       if r["seed"] == row["seed"])
            if abs(eng["cr_over_h"] - ratio) > 1e-12:
                check(f"ratio mismatch on {row['seed']}", False)
            if worst is None or ratio < worst[1]:
                worst = (row["seed"], ratio)
    check("min CR/H over H-bounded endpoints matches engine table",
          worst is not None and abs(worst[1] - 0.0387) < 5e-4,
          f"min = {worst[1]:.4f} at {worst[0]}")

    say("\n" + "=" * 70)
    if FAILS:
        say(f"REFUTATIONS: {len(FAILS)}")
        for f in FAILS:
            say(f"  - {f}")
    else:
        say("No refutation found: every load-bearing number reproduces "
            "under an independent implementation.")
    (DATA / "branch_skeptic_out.txt").write_text("\n".join(LINES) + "\n")
    sys.exit(1 if FAILS else 0)


if __name__ == "__main__":
    main()
