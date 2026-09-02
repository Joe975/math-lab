#!/usr/bin/env python3
"""Exact per-phase fields inside a hierarchical laminate, and the
Hashin-Shtrikman attainment fields.

Why this exists.  Attaining the multiphase HS lower bound (comparison medium
sigma1) forces the field to be UNIFORM inside every phase i != 1, with the value

    E_i = (sigma_HS + sigma1) / (sigma_i + sigma1) * E0                 (*)

(sum_i f_i E_i = E0 recovers the HS identity 1/(s_HS+s1) = sum f_i/(s_i+s1)).
Since sigma_HS > sigma1 always, (*) constrains only the phases with nonvanishing polarization, i.e. i != 1
(phase 1 IS the comparison medium, so its polarization vanishes identically
and its field is free pointwise -- only its volume average is pinned).  And

    E_2 <= E0   <==>   sigma_HS <= sigma2,

which is exactly the feasibility condition of the coated (Milton) assemblage.
So Milton's attainability threshold IS a pointwise field constraint: the
middle phase must not be required to carry more than the applied field.  That
makes "what field does phase 2 actually carry?" a measurable diagnostic on any
candidate structure.

What this computes.  For a laminate tree and an applied average field E0, the
exact field in every leaf.  At an internal node with integer normal n = (u,v),
put S = [[u,-v],[v,u]] and work with the unnormalised frame components
w = S^T E (legitimate because every equation below is homogeneous of degree 1
in the fields, and D_frame = sigma_frame w with the same sigma_frame that
verify_laminate.py uses).  In that frame, with children A (fraction m) and B:

    tangential E continuous:  w^A_2 = w^B_2 = w0_2
    normal D continuous:      (sigma^A w^A)_1 = (sigma^B w^B)_1
    average:                  m w^A + (1-m) w^B = w0

two linear equations for the two unknown normal components, solved in Q.
Recurse with E^child = S w^child / (u^2+v^2).

Self-checks (all exact, in --selftest):
  * energy:      E0 . sigma* E0  ==  sum over leaves of vol * E . sigma_leaf E
  * average:     sum over leaves of vol * E  ==  E0
  * flux:        sum over leaves of vol * sigma_leaf E  ==  sigma* E0
  * the Milton assemblage reproduces (*) exactly in every phase.

Usage:
    python tp_fields.py --selftest
    python tp_fields.py --report        the diagnostic table for attempt 001
Standard library only.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "three-phase-conductivity"))
import laminate as L  # noqa: E402

sys.path.insert(0, HERE)
import tp_coated as TC  # noqa: E402


def conj_S(t, u: int, v: int):
    """S^T t S / (u^2+v^2) with S = [[u,-v],[v,u]] -- the layer-frame tensor."""
    a, b, c = t
    m00 = u * (a * u + b * v) + v * (b * u + c * v)
    m01 = u * (-a * v + b * u) + v * (-b * v + c * u)
    m11 = -v * (-a * v + b * u) + u * (-b * v + c * u)
    d = Fr(u * u + v * v)
    return (m00 / d, m01 / d, m11 / d)


def mv(t, w):
    """tensor times vector."""
    return (t[0] * w[0] + t[1] * w[1], t[1] * w[0] + t[2] * w[1])


def leaf_fields(tree: dict, phases: dict[str, Fr], E0, vol=Fr(1)):
    """[(phase name, volume fraction, field E)] for every leaf, exact."""
    if "phase" in tree:
        return [(tree["phase"], vol, E0)]
    m = Fr(tree["fraction"])
    u, v = (int(x) for x in tree["normal"])
    A = L.effective(tree["layers"][0], phases)
    B = L.effective(tree["layers"][1], phases)
    Af, Bf = conj_S(A, u, v), conj_S(B, u, v)
    # unnormalised frame components of the average field
    w0 = (u * E0[0] + v * E0[1], -v * E0[0] + u * E0[1])
    t = w0[1]  # tangential component, shared by both layers
    # (Af w^A)_1 = (Bf w^B)_1  and  m w^A_1 + (1-m) w^B_1 = w0_1
    #   a = (Bf00 b + (Bf01 - Af01) t) / Af00
    denom = m * Bf[0] / Af[0] + (1 - m)
    b = (w0[0] - m * (Bf[1] - Af[1]) * t / Af[0]) / denom
    a = (Bf[0] * b + (Bf[1] - Af[1]) * t) / Af[0]
    d = Fr(u * u + v * v)
    EA = ((u * a - v * t) / d, (v * a + u * t) / d)
    EB = ((u * b - v * t) / d, (v * b + u * t) / d)
    return (leaf_fields(tree["layers"][0], phases, EA, vol * m)
            + leaf_fields(tree["layers"][1], phases, EB, vol * (1 - m)))


def per_phase(tree, phases, E0):
    """{phase: [(vol, E), ...]} grouped."""
    out: dict[str, list] = {}
    for name, vol, E in leaf_fields(tree, phases, E0):
        out.setdefault(name, []).append((vol, E))
    return out


def hs_target_fields(fracs, sigs, E0):
    """The fields (*) that attaining the HS lower bound would force."""
    s1 = min(sigs)
    hs_lo = 1 / sum(f / (s + s1) for f, s in zip(fracs, sigs)) - s1
    return hs_lo, [tuple((hs_lo + s1) / (s + s1) * c for c in E0) for s in sigs]


def check_identities(tree, phases, E0) -> None:
    """Energy, average-field and average-flux identities -- all exact."""
    sig = L.effective(tree, phases)
    lf = leaf_fields(tree, phases, E0)
    avgE = (sum(v * E[0] for _, v, E in lf), sum(v * E[1] for _, v, E in lf))
    assert avgE == tuple(E0), ("field average", avgE, E0)
    avgD = (sum(v * mv(L.iso(phases[n]), E)[0] for n, v, E in lf),
            sum(v * mv(L.iso(phases[n]), E)[1] for n, v, E in lf))
    assert avgD == mv(sig, E0), ("flux average", avgD, mv(sig, E0))
    en = sum(v * (phases[n] * (E[0] * E[0] + E[1] * E[1])) for n, v, E in lf)
    macro = E0[0] * mv(sig, E0)[0] + E0[1] * mv(sig, E0)[1]
    assert en == macro, ("energy", en, macro)


def field_ratio_range(entries, E0):
    """min and max of |E|/|E0| over the sublayers of one phase (squared, exact)."""
    n0 = E0[0] * E0[0] + E0[1] * E0[1]
    rs = [(E[0] * E[0] + E[1] * E[1]) / n0 for _, E in entries]
    return min(rs), max(rs)


def selftest() -> None:
    F = Fr
    sigma = (F(1), F(2), F(5))
    phases = {"p1": sigma[0], "p2": sigma[1], "p3": sigma[2]}
    E0 = (F(1), F(1))

    # identities on assorted trees, isotropic and not
    trees = [
        {"fraction": "1/3", "normal": [1, 0],
         "layers": [{"phase": "p1"},
                    {"fraction": "3/7", "normal": [1, 2],
                     "layers": [{"phase": "p2"}, {"phase": "p3"}]}]},
        {"fraction": "2/7", "normal": [0, 1],
         "layers": [{"fraction": "1/2", "normal": [3, -1],
                     "layers": [{"phase": "p3"}, {"phase": "p1"}]},
                    {"fraction": "5/9", "normal": [2, 5],
                     "layers": [{"phase": "p2"}, {"phase": "p1"}]}]},
    ]
    for tr in trees:
        for e in [(F(1), F(0)), (F(0), F(1)), (F(1), F(1)), (F(3), F(-2))]:
            check_identities(tr, phases, e)

    # the Milton assemblage realises the HS attainment fields (*) exactly
    for f in [[F(3, 8), F(2, 8), F(3, 8)], [F(5, 8), F(2, 8), F(1, 8)],
              [F(6, 8), F(1, 8), F(1, 8)]]:
        tree, ph, info = TC.milton(f, sigma, 0)
        assert tree is not None, "expected a feasible assemblage"
        check_identities(tree, ph, E0)
        hs_lo, targets = hs_target_fields(f, list(sigma), E0)
        got = per_phase(tree, ph, E0)
        # phases 2 and 3 (nonzero polarization) must hit (*) exactly;
        # phase 1 is the comparison medium, its polarization vanishes and its
        # field is free -- only its volume average is pinned.
        for i, name in [(1, "p2"), (2, "p3")]:
            for _, E in got[name]:
                assert E == targets[i], (name, E, targets[i])
    print("selftest: field identities exact; Milton assemblage matches (*)")


def report() -> None:
    F = Fr
    sigma = [F(1), F(2), F(5)]
    phases = {"p1": sigma[0], "p2": sigma[1], "p3": sigma[2]}
    E0 = (F(1), F(1))
    cert = os.path.join(HERE, "..", "data", "certified")
    print("Per-phase field ratios |E|^2/|E0|^2 in certified structures")
    print("target = the value HS attainment would force, (s_HS+s1)^2/(s_i+s1)^2")
    print()
    for name in sorted(os.listdir(cert)):
        if not name.endswith(".json") or name.endswith(".meta.json"):
            continue
        if "lower" not in name:
            continue
        rec = json.load(open(os.path.join(cert, name), encoding="utf-8"))
        ph = {k: Fr(v) for k, v in rec["phases"].items()}
        if sorted(ph.values()) != sigma:
            continue
        tree = rec["tree"]
        check_identities(tree, ph, E0)
        fr = L.fractions_of(tree)
        names = sorted(fr)
        fracs = [fr[n] for n in names]
        sigs = [ph[n] for n in names]
        hs_lo, targets = hs_target_fields(fracs, sigs, E0)
        got = per_phase(tree, ph, E0)
        n0 = 2
        line = [f"{name[:-5]:28s}"]
        for i, nm in enumerate(["p1", "p2", "p3"]):
            lo, hi = field_ratio_range(got[nm], E0)
            tg = (targets[i][0] ** 2 + targets[i][1] ** 2) / n0
            line.append(f"{nm}: [{float(lo):.4f},{float(hi):.4f}] tgt {float(tg):.4f}")
        print("  ".join(line))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
    if a.report:
        report()


# ---------------------------------------------------------------------------
# Exact decomposition of the HS gap
# ---------------------------------------------------------------------------


def required_variance(fracs, sigs):
    """The phase-1 field variance that attaining HS_lo forces, per |E0|^2.

    With c = 1/sum_i f_i/(s_i+s1) and S2 = sum_i f_i/(s_i+s1)^2,

        Var_1 = (c^2 S2 - 1) / f_1,

    which is STRICTLY POSITIVE unless all phases coincide: Cauchy-Schwarz on
    (sum f_i/(s_i+s1))^2 <= (sum f_i)(sum f_i/(s_i+s1)^2) gives c^2 S2 >= 1.
    So attaining the bound does not merely permit field fluctuation in the
    comparison-medium phase, it requires a precisely determined amount of it.
    """
    s1 = min(sigs)
    c = 1 / sum(f / (s + s1) for f, s in zip(fracs, sigs))
    S2 = sum(f / (s + s1) ** 2 for f, s in zip(fracs, sigs))
    f1 = fracs[list(sigs).index(s1)]
    return (c * c * S2 - 1) / f1


def moments(entries, E0):
    """(volume, mean field, variance per |E0|^2) of one phase, exact."""
    n0 = E0[0] * E0[0] + E0[1] * E0[1]
    tot = sum(v for v, _ in entries)
    mean = tuple(sum(v * E[k] for v, E in entries) / tot for k in (0, 1))
    var = sum(v * ((E[0] - mean[0]) ** 2 + (E[1] - mean[1]) ** 2)
              for v, E in entries) / tot / n0
    return tot, mean, var


def gap_decomposition(tree, phases, E0):
    """Exact split of sigma* - HS_lo into mean-deviation and variance parts.

        E0.sigma* E0/|E0|^2 - HS_lo = sum_i f_i s_i (|mean_i|^2 - |E_i^tgt|^2)/|E0|^2
                         + sum_i f_i s_i Var_i/|E0|^2
                         - f_1 s_1 Var_req

    Both sides are computed independently here and asserted equal, so this is
    a checked identity rather than an assertion.
    """
    n0 = E0[0] * E0[0] + E0[1] * E0[1]
    fr = L.fractions_of(tree)
    names = sorted(fr)
    fracs = [fr[n] for n in names]
    sigs = [phases[n] for n in names]
    hs_lo, targets = hs_target_fields(fracs, sigs, E0)
    got = per_phase(tree, phases, E0)
    s1 = min(sigs)
    mean_term = Fr(0)
    var_term = Fr(0)
    for i, nm in enumerate(names):
        vol, mean, var = moments(got[nm], E0)
        tg = targets[i]
        mean_term += fracs[i] * sigs[i] * (
            (mean[0] ** 2 + mean[1] ** 2) - (tg[0] ** 2 + tg[1] ** 2)) / n0
        var_term += fracs[i] * sigs[i] * var
    f1 = fracs[list(sigs).index(s1)]
    total = mean_term + var_term - f1 * s1 * required_variance(fracs, sigs)
    sig = L.effective(tree, phases)
    # use the quadratic form E0.sigma* E0/|E0|^2, exact for any tensor, so the
    # identity does not need the tree to be exactly isotropic (certified
    # records carry a tiny rationalisation residual).
    eff = (E0[0] * mv(sig, E0)[0] + E0[1] * mv(sig, E0)[1]) / n0
    assert total == eff - hs_lo, (total, eff - hs_lo)
    return {"gap": eff - hs_lo, "mean_term": mean_term,
            "var_term": var_term, "var_required": f1 * s1 * required_variance(fracs, sigs)}
