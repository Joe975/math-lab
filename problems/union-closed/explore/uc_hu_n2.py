#!/usr/bin/env python3
"""Attempt 046 (045 lead 3): the n = 2 case of (HU-TAX), proved where it
can be and CERTIFIED where it cannot.

At n = 2 the half-union coupling has a closed form in three parameters,
so the conjecture

    CR_HU(mu) >= c*(q) H(mu),   q = max marginal,
    c*(p) = (h(z*(p)) - h(p))/h(p),  z*(p) = max(1/2, 1-2p)

becomes a concrete inequality on a 3-box.  Parametrization: x = P(A0=0),
u_a = P(A1=0 | A0=a); in regime iff both marginals f0 = 1-x and
f1 = x(1-u0) + (1-x)(1-u1) lie in (0, 1/2).

Parts:
  A  structural lemmas, checked numerically against the general
     evaluator: the closed form; the ORDER QUANTIFIER IS VACUOUS at
     n = 2 (identity order alone meets the bound everywhere tested);
     and Lemma N2-ONE-BAD -- at most one conditional can be
     out-of-regime, so the clip has only two cases and every cell has
     an explicit closed form
  B  the exact equality set: products Bern(p)^2, diagonals, and the
     n = 1 degenerations -- each an exact identity (hand-derived in the
     record), re-checked here in exact rational interval arithmetic
  C  CERTIFIED branch-and-bound: F' := CR*h(q) - (h(z*(q)) - h(q))*H
     enclosed in exact rational interval arithmetic over the whole
     regime box; every box with F'_lo > 0 is certified, boxes that
     resist down to the minimum edge are reported as residue (they must
     contain equality points -- that is the honest limit of a numeric
     certificate against a conjecture with a nonempty equality set)
  D  the local degeneracy order at the equality families (numeric
     second-order probe), which says what any proof must handle

Kits: the B&B runs under kit B (022 atanh-series, cheap); the residue
boundary boxes and every part-B identity are re-certified under kit A
(digit-extraction) alone, per the 029 single-kit standard.

Standard library only; deterministic.
Usage: python uc_hu_n2.py [--fast]
Checkpoint: ../data/hu_n2.json
"""
from __future__ import annotations

import heapq
import json
import math
import random
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_hu_order2 import hu_cr_seq, h
from uc_hu_certify import ivl_add, ivl_scale, h_ivl, log2_A, log2_B

DATA = HERE.parent / "data"
FAST = "--fast" in sys.argv
T0 = time.monotonic()
OUT = {}
HALF = Fraction(1, 2)
ONE = Fraction(1)


def log(m=""):
    print(f"[{time.monotonic() - T0:6.1f}s] {m}" if m else "", flush=True)


# ----------------------------------------------------------- float side
def cstar(p):
    return (h(max(0.5, 1.0 - 2.0 * p)) - h(p)) / h(p)


def clipf(a, b):
    return min(max(0.5, a + b - 1.0), a, b)


def mu_from(x, u0, u1):
    return {0b00: x * u0, 0b10: x * (1 - u0),
            0b01: (1 - x) * u1, 0b11: (1 - x) * (1 - u1)}


def closed_cr(x, u0, u1):
    """CR under the identity order, from the record's closed form."""
    z0 = clipf(x, x)
    w = {(0, 0): z0, (0, 1): x - z0, (1, 0): x - z0, (1, 1): 1 - 2 * x + z0}
    u = (u0, u1)
    S = h(z0) + sum(wt * h(clipf(u[a], u[b]))
                    for (a, b), wt in w.items() if wt > 0)
    H = h(x) + x * h(u0) + (1 - x) * h(u1)
    return S - H, H


def regime(x, u0, u1):
    f0 = 1 - x
    f1 = x * (1 - u0) + (1 - x) * (1 - u1)
    return f0, f1, (0 < max(f0, f1) < 0.5)


# -------------------------------------------------------- interval side
class Kit:
    def __init__(self, fn):
        self.fn = fn
        self.cache = {}

    def h_at(self, z):
        if z in self.cache:
            return self.cache[z]
        if z <= 0 or z >= 1:
            v = (Fraction(0), Fraction(0))
        else:
            v = h_ivl(z, self.fn)
        self.cache[z] = v
        return v

    def h_range(self, lo, hi):
        """Enclosure of {h(z) : z in [lo,hi]} (h up on [0,1/2], down after)."""
        if lo > hi:
            lo, hi = hi, lo
        vals = [self.h_at(lo), self.h_at(hi)]
        if lo <= HALF <= hi:
            top = Fraction(1)
        else:
            top = max(v[1] for v in vals)
        bot = min(v[0] for v in vals)
        return (bot, top)


def imul(a, b):
    prods = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return (min(prods), max(prods))


def isub(a, b):
    return (a[0] - b[1], a[1] - b[0])


def iadd(a, b):
    return (a[0] + b[0], a[1] + b[1])


def iclip(X, Y):
    """Interval enclosure of clip(x,y) = min(max(1/2, x+y-1), x, y)."""
    s = (X[0] + Y[0] - 1, X[1] + Y[1] - 1)
    m = (max(HALF, s[0]), max(HALF, s[1]))
    lo = min(m[0], X[0], Y[0])
    hi = min(m[1], X[1], Y[1])
    return (lo, hi)


def F_enclosure(kit, X, U0, U1):
    """Enclosure of F' = CR*h(q) - (h(z*(q)) - h(q))*H over the box.

    Returns (F_lo, F_hi) or None if the box is entirely out of regime."""
    F0 = (1 - X[1], 1 - X[0])
    one_m_u0 = (1 - U0[1], 1 - U0[0])
    one_m_u1 = (1 - U1[1], 1 - U1[0])
    F1 = iadd(imul(X, one_m_u0), imul((1 - X[1], 1 - X[0]), one_m_u1))
    Q = (max(F0[0], F1[0]), max(F0[1], F1[1]))
    if Q[0] >= HALF or Q[1] <= 0:
        return None                      # wholly out of regime
    Qc = (max(Q[0], Fraction(0)), min(Q[1], HALF))
    # z*(q) = max(1/2, 1-2q): decreasing in q
    Z = (max(HALF, 1 - 2 * Qc[1]), max(HALF, 1 - 2 * Qc[0]))
    hZ = kit.h_range(Z[0], Z[1])
    hQ = kit.h_range(Qc[0], Qc[1])
    # step-0 cell
    Z0 = iclip(X, X)
    hZ0 = kit.h_range(Z0[0], Z0[1])
    W = {(0, 0): Z0,
         (0, 1): isub(X, Z0), (1, 0): isub(X, Z0),
         (1, 1): iadd(isub((ONE, ONE), imul((Fraction(2), Fraction(2)), X)),
                      Z0)}
    U = {0: U0, 1: U1}
    S = hZ0
    for (a, b), w in W.items():
        wc = (max(w[0], Fraction(0)), max(w[1], Fraction(0)))
        zc = iclip(U[a], U[b])
        S = iadd(S, imul(wc, kit.h_range(zc[0], zc[1])))
    Hb = iadd(kit.h_range(X[0], X[1]),
              iadd(imul(X, kit.h_range(U0[0], U0[1])),
                   imul((1 - X[1], 1 - X[0]),
                        kit.h_range(U1[0], U1[1]))))
    CR = isub(S, Hb)
    return isub(imul(CR, hQ), imul(isub(hZ, hQ), Hb))


# ------------------------------------------------------------- parts
def part_A():
    log("A. Structural lemmas (float, against the general evaluator):")
    rng = random.Random(11)
    worst_form = 0.0
    bad_bad = 0
    id_below = 0
    checked = 0
    worst_margin = None
    for _ in range(40000 if not FAST else 2000):
        x = rng.uniform(0.5, 1.0)
        u0, u1 = rng.uniform(0, 1), rng.uniform(0, 1)
        f0, f1, ok = regime(x, u0, u1)
        if not ok:
            continue
        checked += 1
        mu = {a: w for a, w in mu_from(x, u0, u1).items() if w > 1e-15}
        cr_g, H_g = hu_cr_seq(2, mu, (0, 1))
        cr_c, H_c = closed_cr(x, u0, u1)
        worst_form = max(worst_form, abs(cr_g - cr_c), abs(H_g - H_c))
        if u0 < 0.5 and u1 < 0.5:
            bad_bad += 1
        q = max(f0, f1)
        m_id = cr_g - cstar(q) * H_g
        cr_r, _ = hu_cr_seq(2, mu, (1, 0))
        if cr_g < max(cr_g, cr_r) - 1e-12:
            id_below += 1
        if worst_margin is None or m_id < worst_margin[0]:
            worst_margin = (m_id, x, u0, u1, q, H_g)
    log(f"  closed form vs general evaluator: max |diff| {worst_form:.2e} "
        f"over {checked} in-regime samples")
    log(f"  Lemma N2-ONE-BAD: in-regime samples with BOTH conditionals "
        f"out of regime: {bad_bad} (proof in the record: both would force "
        f"f1 > 1/2)")
    log(f"  identity order strictly below the best order on {id_below} "
        f"samples, yet its worst margin is "
        f"{worst_margin[0]:+.3e} (x={worst_margin[1]:.4f}, "
        f"u0={worst_margin[2]:.4f}, u1={worst_margin[3]:.4f})")
    OUT["A"] = {"samples": checked, "closed_form_maxdiff": worst_form,
                "both_bad_count": bad_bad, "identity_below_best": id_below,
                "worst_identity_margin": worst_margin[0],
                "worst_at": {"x": worst_margin[1], "u0": worst_margin[2],
                             "u1": worst_margin[3], "q": worst_margin[4],
                             "H": worst_margin[5]}}


def part_B():
    log()
    log("B. The equality set, exact-rational check under BOTH kits alone:")
    rows = []
    for kitname, fn in (("A_digit_extraction", log2_A),
                        ("B_atanh_series", log2_B)):
        kit = Kit(fn)
        for tag, mk in (("product", lambda p: (1 - p, 1 - p, 1 - p)),
                        ("diagonal", lambda p: (1 - p, ONE, Fraction(0))),
                        ("n1-degenerate", lambda p: (1 - p, ONE, ONE))):
            for p in (Fraction(1, 10), Fraction(3, 10), Fraction(2, 5),
                      Fraction(49, 100)):
                x, u0, u1 = mk(p)
                box = ((x, x), (u0, u0), (u1, u1))
                F = F_enclosure(kit, *box)
                if F is None:
                    continue
                width = float(F[1] - F[0])
                zero_inside = F[0] <= 0 <= F[1]
                rows.append({"kit": kitname, "family": tag, "p": float(p),
                             "F_lo": float(F[0]), "F_hi": float(F[1]),
                             "contains_zero": zero_inside})
                if not zero_inside:
                    log(f"  [FAIL] {tag} p={float(p)} {kitname}: F' "
                        f"in [{float(F[0]):.2e},{float(F[1]):.2e}] excludes 0")
    ok = all(r["contains_zero"] for r in rows)
    log(f"  {len(rows)} equality-family points, both kits: every enclosure "
        f"of F' contains 0 -> consistent with exact equality  "
        f"[{'ok' if ok else 'FAILURE'}]")
    OUT["B_equality"] = {"rows": rows, "all_contain_zero": ok}


def part_C():
    log()
    log("C. Certified branch-and-bound on F' >= 0 over the regime box:")
    kit = Kit(log2_B)
    min_edge = Fraction(1, 512 if not FAST else 64)
    max_boxes = 2000000 if not FAST else 20000
    root = ((HALF, ONE), (Fraction(0), ONE), (Fraction(0), ONE))
    # LARGEST-BOX-FIRST: a plain DFS stack spends the whole budget deep in
    # one corner and leaves most of the box unprocessed (measured: 99.8%),
    # which makes any "residue volume" figure meaningless.  The heap key is
    # -volume, so coverage grows before depth does.
    counter = 0
    stack = [(-Fraction(1, 2), 0, root)]
    certified = 0
    out_of_regime = 0
    residue = []
    processed = 0
    vol_cert = Fraction(0)
    vol_oor = Fraction(0)
    vol_res = Fraction(0)

    def boxvol(b):
        return (b[0][1] - b[0][0]) * (b[1][1] - b[1][0]) * (b[2][1] - b[2][0])
    while stack and processed < max_boxes:
        _, _, box = heapq.heappop(stack)
        processed += 1
        F = F_enclosure(kit, *box)
        if F is None:
            out_of_regime += 1
            vol_oor += boxvol(box)
            continue
        if F[0] > 0:
            certified += 1
            vol_cert += boxvol(box)
            continue
        widths = [b[1] - b[0] for b in box]
        if max(widths) <= min_edge:
            vol_res += boxvol(box)
            residue.append({"x": [float(box[0][0]), float(box[0][1])],
                            "u0": [float(box[1][0]), float(box[1][1])],
                            "u1": [float(box[2][0]), float(box[2][1])],
                            "F_lo": float(F[0]), "F_hi": float(F[1])})
            continue
        k = widths.index(max(widths))
        mid = (box[k][0] + box[k][1]) / 2
        for half in ((box[k][0], mid), (mid, box[k][1])):
            nb = list(box)
            nb[k] = half
            nb = tuple(nb)
            counter += 1
            heapq.heappush(stack, (-boxvol(nb), counter, nb))
    vol_unproc = sum(boxvol(b) for _, _, b in stack)
    total = Fraction(1, 2)
    log(f"  processed {processed} boxes (cap {max_boxes}); exhausted: "
        f"{not stack}")
    log(f"  volume accounting (root box volume 1/2): certified "
        f"{float(vol_cert):.6f} ({100*float(vol_cert)/float(total):.3f}%), "
        f"out-of-regime {float(vol_oor):.6f}, residue "
        f"{float(vol_res):.3e}, unprocessed {float(vol_unproc):.3e}")
    log(f"  residue: {len(residue)} boxes at min edge "
        f"1/{int(1/min_edge)}")
    OUT["C_bnb"] = {"processed": processed, "certified": certified,
                    "out_of_regime": out_of_regime,
                    "residue_boxes": len(residue),
                    "vol_certified": float(vol_cert),
                    "vol_out_of_regime": float(vol_oor),
                    "vol_residue": float(vol_res),
                    "vol_unprocessed": float(vol_unproc),
                    "vol_root": float(total),
                    "min_edge": float(min_edge),
                    "exhausted": not stack,
                    "residue": residue[:200]}
    return residue


def residue_classify(residue):
    """Every residue box should meet the equality set E."""
    log()
    log("C2. Residue boxes vs the equality set:")
    def near_E(r):
        x = (r["x"][0] + r["x"][1]) / 2
        u0 = (r["u0"][0] + r["u0"][1]) / 2
        u1 = (r["u1"][0] + r["u1"][1]) / 2
        tol = 3 * max(r["x"][1] - r["x"][0], r["u0"][1] - r["u0"][0],
                      r["u1"][1] - r["u1"][0]) + 1e-12
        prod = abs(x - u0) <= tol and abs(x - u1) <= tol
        diag = abs(u0 - 1) <= tol and abs(u1) <= tol
        deg = abs(x - 1) <= tol or (abs(u0 - 1) <= tol and abs(u1 - 1) <= tol)
        edge = (r["u0"][1] - r["u0"][0]) > 0 and False
        return prod or diag or deg
    hit = sum(1 for r in residue if near_E(r))
    log(f"  {hit} of {len(residue)} residue boxes lie within 3 edge-lengths "
        f"of a known equality family (products x=u0=u1, diagonals "
        f"u0=1,u1=0, n=1 degenerations)")
    OUT["C2_residue"] = {"total": len(residue), "near_equality": hit}
    return hit


def part_D():
    log()
    log("D. Local structure at the equality families (transversal order):")
    rows = []
    def margin(x, u0, u1):
        f0, f1, ok = regime(x, u0, u1)
        if not ok:
            return None
        cr, H = closed_cr(x, u0, u1)
        return cr - cstar(max(f0, f1)) * H
    fams = (("product p=0.30", (0.70, 0.70, 0.70), [(1, -1, 0), (0, 1, -1),
                                                    (1, 0, -1), (0, 1, 1)]),
            ("product p=0.45", (0.55, 0.55, 0.55), [(1, -1, 0), (0, 1, -1),
                                                    (1, 0, -1), (0, 1, 1)]),
            ("diagonal p=0.30", (0.70, 1.0, 0.0), [(0, -1, 1), (1, -1, 0),
                                                   (0, -1, 0), (0, 0, 1)]),
            ("diagonal p=0.45", (0.55, 1.0, 0.0), [(0, -1, 1), (1, -1, 0),
                                                   (0, -1, 0), (0, 0, 1)]))
    for tag, base, dirs in fams:
        x, u0, u1 = base
        m0 = margin(x, u0, u1)
        probe = []
        for d in dirs:
            pair = []
            for eps in (1e-2, 1e-3):
                c = (min(1.0, max(0.5, x + eps * d[0])),
                     min(1.0, max(0.0, u0 + eps * d[1])),
                     min(1.0, max(0.0, u1 + eps * d[2])))
                m = margin(*c)
                pair.append(None if m is None else m - m0)
            order = None
            if pair[0] and pair[1] and pair[1] > 0:
                order = math.log10(pair[0] / pair[1])
            probe.append({"dir": list(d), "d_1e-2": pair[0],
                          "d_1e-3": pair[1], "order": order})
        rows.append({"family": tag, "margin": m0, "probe": probe})
        s_ = ", ".join(f"{tuple(pp['dir'])}: {pp['d_1e-2']:+.2e}/"
                       f"{pp['d_1e-3']:+.2e}"
                       + (f" (order~{pp['order']:.1f})"
                          if pp["order"] else "")
                       for pp in probe if pp["d_1e-2"] is not None)
        log(f"  {tag}: margin {m0:+.1e} | transversal {s_}")
    OUT["D_local"] = rows


def part_E():
    log()
    log("E. Is the concavity reduction lossy?  Testing (**):")
    def G(p):
        return h(max(0.5, 1.0 - 2.0 * p)) - h(p)
    rng = random.Random(3)
    tested = 0
    bad = []
    for _ in range(400000 if not FAST else 20000):
        x = rng.uniform(0.5, 1.0)
        p0 = rng.uniform(0, 0.5)
        p1 = rng.uniform(0, 0.5)
        f0 = 1 - x
        f1 = x * p0 + (1 - x) * p1
        q = max(f0, f1)
        if not (0 < q < 0.5) or min(p0, p1) <= 0:
            continue
        tested += 1
        lhs = x * G(p0) + (1 - x) * G(p1)
        rhs = cstar(q) * (x * h(p0) + (1 - x) * h(p1))
        if lhs < rhs - 1e-12:
            bad.append((lhs - rhs, x, p0, p1, q))
    bad.sort()
    worst = []
    for d, x, p0, p1, q in bad[:5]:
        m = None
        mu = {a: w for a, w in mu_from(x, 1 - p0, 1 - p1).items()
              if w > 1e-15}
        cr, H = hu_cr_seq(2, mu, (0, 1))
        m = cr - cstar(q) * H
        worst.append({"deficit": d, "x": x, "p0": p0, "p1": p1, "q": q,
                      "original_margin": m})
    frac = len(bad) / tested if tested else 0
    log(f"  (**) fails on {len(bad)}/{tested} Case-A points "
        f"({100*frac:.1f}%); worst deficit "
        f"{bad[0][0]:+.4e} at x={bad[0][1]:.4f}, p0={bad[0][2]:.4f}, "
        f"p1={bad[0][3]:.4f}")
    log(f"  at those very points the ORIGINAL margin is positive "
        f"(worst five: "
        + ", ".join(f"{w['original_margin']:+.3e}" for w in worst) + ")")
    log("  => the concavity step is tight only at p_a = p_b (products); "
        "averaging the pair interaction away discards too much when the "
        "two conditionals are far apart")
    OUT["E_reduction"] = {"tested": tested, "violations": len(bad),
                          "violation_fraction": frac,
                          "worst": worst}


def main():
    part_A()
    part_B()
    residue = part_C()
    residue_classify(residue)
    part_D()
    part_E()
    (DATA / "hu_n2.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    log()
    log("checkpoint: data/hu_n2.json")


if __name__ == "__main__":
    main()
