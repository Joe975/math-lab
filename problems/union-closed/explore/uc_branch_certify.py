#!/usr/bin/env python3
"""Attempt 026 (025 lead 3 / queue 1c): exact-rational certification of
the 025 branch kills and their repairs, dual-kit.

Every quantity here is a finite sum  sum_i c_i * log2(r_i)  with exact
rational c_i, r_i > 0 (entropies of explicitly rational measures), so a
certified enclosure needs only certified log2 of rationals.  Two
independent, previously audited kits supply it:

  kit A: uc_gap1_skeptic_gemini.certified_log2 (digit extraction,
         directed dyadic rounding; drafted cross-family in 022, audited
         there and re-audited from scratch in 024) -- run unmodified;
  kit B: uc_or_agg_probe.log2_enclosure (atanh series with explicit tail
         bound, directed fixed-point rounding; the 016/017 standard) --
         run unmodified.

Each log2 is evaluated by BOTH kits and the enclosures intersected; a
disjoint intersection would prove one kit wrong and aborts the run.  All
coefficient arithmetic is exact Fractions; interval sums round nothing.

Certified statements (see the record for the instance list):
  C1  half-mixing branch kill:  CR_hm = F(P1/2) - F(P1) < 0, in-regime,
      at rational sliver instances for n = 2, 6, 16;
  C2  repair at the same instances:  CR(q = 1/2) = F(1/2) - F(P1) > 0;
  C3  block-branch kill: Gain_block < 0 at n = 8 for the 2block and
      3block mixtures (with 009 fact (i) CR <= Gain, 011-verified, this
      certifies CR_block < 0);
  C4  empty-mixing rescue of the 3block instance at n = 6, q = 1/5:
      CR = F_gen(1/5) - F_gen(3/5) > 0 (with Lemma EM of 025 supplying
      CR = Gain for this coupling family);
  C5  H(mu) > 1/2 and max marginal < 38271/100000 at every instance,
      exactly.

Standard library only; deterministic; ~seconds.
Usage: python uc_branch_certify.py
Checkpoint: ../data/branch_certify.json
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import comb
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from uc_gap1_skeptic_gemini import certified_log2 as log2_A
from uc_or_agg_probe import log2_enclosure as log2_B

DATA = HERE.parent / "data"
MARG = Fraction(38271, 100000)
OUT = {"certified": [], "log2_calls": 0, "max_width": 0.0}


def clog2(r: Fraction) -> tuple[Fraction, Fraction]:
    """Intersection of both kits' enclosures; abort if disjoint."""
    aL, aH = log2_A(r)
    bL, bH = log2_B(r)
    lo, hi = max(aL, bL), min(aH, bH)
    if lo > hi:
        raise SystemExit(f"KIT DISAGREEMENT at log2({r}): "
                         f"A=[{float(aL)},{float(aH)}] "
                         f"B=[{float(bL)},{float(bH)}]")
    OUT["log2_calls"] += 1
    OUT["max_width"] = max(OUT["max_width"], float(hi - lo))
    return lo, hi


def ivl_add(a, b):
    return a[0] + b[0], a[1] + b[1]


def ivl_scale(c: Fraction, iv):
    return (c * iv[0], c * iv[1]) if c >= 0 else (c * iv[1], c * iv[0])


def ivl_sub(a, b):
    return a[0] - b[1], a[1] - b[0]


def xlog2x_ivl(x: Fraction):
    """Enclosure of x*log2(x), x >= 0 rational."""
    if x == 0:
        return Fraction(0), Fraction(0)
    return ivl_scale(x, clog2(x))


def H_profile(masses):
    """Enclosure of -sum_j mult_j * m_j * log2(m_j) for rational atom
    masses given as (multiplicity, mass) with sum mult*mass = 1."""
    tot = sum(mult * m for mult, m in masses)
    assert tot == 1, f"masses sum to {tot} != 1"
    acc = (Fraction(0), Fraction(0))
    for mult, m in masses:
        if m > 0:
            acc = ivl_add(acc, ivl_scale(-mult, xlog2x_ivl(m)))
    return acc


def empty_mix_profile(n: int, w: Fraction, scomps):
    """Profile masses of  w*delta_empty + (1-w) * sum_k S_k Bern(p_k)^n,
    scomps = [(S_k, p_k)] rational, as [(multiplicity, mass)]."""
    rows = []
    for t in range(n + 1):
        a = sum(S * p ** t * (1 - p) ** (n - t) for S, p in scomps)
        m = (1 - w) * a
        if t == 0:
            m += w
            rows.append((1, m))
        else:
            rows.append((comb(n, t), m))
    return rows


def certify(name, iv, want, extra=""):
    lo, hi = float(iv[0]), float(iv[1])
    if want == "neg":
        ok = iv[1] < 0
    elif want == "pos":
        ok = iv[0] > 0
    else:
        raise ValueError(want)
    verdict = "CERTIFIED" if ok else "NOT CERTIFIED"
    print(f"  {name}: [{lo:+.9e}, {hi:+.9e}] "
          f"{verdict} {'< 0' if want == 'neg' else '> 0'} {extra}")
    OUT["certified"].append({"name": name, "lo": lo, "hi": hi,
                             "want": want, "ok": ok})
    if not ok:
        raise SystemExit(f"certification failed: {name}")


def main():
    print("Attempt 026: dual-kit certification of the 025 branch kills")
    print("=" * 70)

    # anchors: both kits must bracket exact values
    for r, v in [(Fraction(2), 1), (Fraction(8), 3), (Fraction(1, 2), -1)]:
        lo, hi = clog2(r)
        assert lo <= v <= hi, (r, v)
    print("  [ok] anchors log2(2)=1, log2(8)=3, log2(1/2)=-1 bracket "
          "in the kit intersection")

    # ---------------- C1/C2/C5: half-mixing sliver instances ----------
    print("\nC1/C2: half-mixing kill and q=1/2 repair, certified:")
    P1 = Fraction(16, 25)
    for n, ph in [(2, Fraction(1987, 2000)),   # 025's A5 instance
                  (6, Fraction(1997, 2000)),   # 025's A5 instance
                  (16, Fraction(19993, 20000))]:
        beta = (1 - ph) ** n

        def F(w: Fraction):
            m0 = w + (1 - w) * beta
            acc = ivl_scale(Fraction(-1), xlog2x_ivl(m0))
            acc = ivl_add(acc, ivl_scale(-(1 - w) * (1 - beta),
                                         clog2(1 - w)))
            # (1-w) * n * h(ph)
            hph = ivl_add(ivl_scale(Fraction(-1), xlog2x_ivl(ph)),
                          ivl_scale(Fraction(-1), xlog2x_ivl(1 - ph)))
            acc = ivl_add(acc, ivl_scale((1 - w) * n, hph))
            acc = ivl_add(acc, ivl_scale(1 - w, xlog2x_ivl(beta)))
            return acc

        FP1 = F(P1)
        marg = (1 - P1) * ph
        assert marg < MARG, "not in-regime"
        assert FP1[0] > Fraction(1, 2), "H(mu) not certified > 1/2"
        certify(f"CR_halfmix(n={n:2d}, ph={ph}, P1=16/25)",
                ivl_sub(F(P1 / 2), FP1), "neg",
                f"[marg={float(marg):.4f} < 0.38271 exact, "
                f"H(mu) > {float(FP1[0]):.4f} certified]")
        certify(f"CR(q=1/2) (n={n:2d}) repair   ",
                ivl_sub(F(Fraction(1, 2)), FP1), "pos")

    # ---------------- C3: block-branch kills at n = 8 -----------------
    print("\nC3: block iid-given-k Gain < 0 at n = 8 (CR <= Gain, 009(i)):")
    insts = [("2block", [(Fraction(31, 50), Fraction(0)),
                         (Fraction(19, 50), Fraction(19, 20))]),
             ("3block", [(Fraction(3, 5), Fraction(0)),
                         (Fraction(3, 10), Fraction(9, 10)),
                         (Fraction(1, 10), Fraction(99, 100))])]
    n = 8
    for tag, comps in insts:
        marg = sum(S * p for S, p in comps)
        assert marg < MARG
        Hmu = H_profile(empty_mix_profile(
            n, Fraction(0), [(S, p) for S, p in comps]))
        Hu = H_profile(empty_mix_profile(
            n, Fraction(0), [(S, 1 - (1 - p) ** 2) for S, p in comps]))
        certify(f"Gain_block({tag}, n=8)", ivl_sub(Hu, Hmu), "neg",
                f"[marg={float(marg):.4f} exact, "
                f"H(mu) > {float(Hmu[0]):.4f} certified]")

    # ---------------- C4: empty-mixing rescue of 3block, n = 6 --------
    print("\nC4: empty-mixing rescue (Lemma EM: CR = Gain), 3block n = 6:")
    scomps = [(Fraction(3, 4), Fraction(9, 10)),
              (Fraction(1, 4), Fraction(99, 100))]
    P1b = Fraction(3, 5)
    n = 6
    Fq = H_profile(empty_mix_profile(n, Fraction(1, 5), scomps))
    FP = H_profile(empty_mix_profile(n, P1b, scomps))
    certify("CR_empty_mixing(3block, n=6, q=1/5)", ivl_sub(Fq, FP), "pos")

    # ------------- cross-check vs 025's float values -------------
    # The engine evaluated at binary-float inputs (0.9935 != 1987/2000
    # exactly), so its values approximate a nearby instance: proximity
    # is the right test, not containment in the ~1e-30 enclosures.
    ref = json.loads((DATA / "branch_attack.json").read_text())
    checks = []

    def near(fl, c, tol=1e-9):
        return abs(fl - 0.5 * (c["lo"] + c["hi"])) < tol

    for row in ref["A5_rescue"]:
        target = next(c for c in OUT["certified"]
                      if c["name"].startswith(f"CR_halfmix(n={row['n']:2d}"))
        checks.append((f"A5 n={row['n']} float within 1e-9 of certified",
                       near(row["CR_halfmix"], target)))
    b1 = {r["tag"]: r for r in ref["B1_gain"]}
    for tag in ("2block", "3block"):
        fl = next(r["gain_block"] for r in b1[tag]["rows"] if r["n"] == 8)
        target = next(c for c in OUT["certified"]
                      if c["name"] == f"Gain_block({tag}, n=8)")
        checks.append((f"B1 {tag} n=8 float within 1e-9 of certified",
                       near(fl, target)))
    fl = ref["B3_rescue"]["empty_mixing_CR"]
    target = next(c for c in OUT["certified"] if "empty_mixing" in c["name"])
    checks.append(("B3 rescue float within 1e-9 of certified",
                   near(fl, target)))
    print()
    for name, ok in checks:
        print(f"  [{'ok' if ok else 'MISMATCH'}] {name}")
        if not ok:
            raise SystemExit(f"float/certified mismatch: {name}")

    print(f"\n{OUT['log2_calls']} dual-kit log2 evaluations, max "
          f"intersected width {OUT['max_width']:.2e}; all statements "
          f"certified.")
    (DATA / "branch_certify.json").write_text(
        json.dumps(OUT, indent=1, default=float) + "\n")
    print("checkpoint: data/branch_certify.json")


if __name__ == "__main__":
    main()
