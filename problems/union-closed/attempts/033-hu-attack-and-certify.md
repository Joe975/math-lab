# 033 — HU attacked with the order as a weapon: the for-all-orders form dies at cap 0.45 (certified), the working-cap form survives and its floor is certified

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** adversarial campaign (queue 1(b1)) + exact certification
  (queue 1(b3)) against/of the half-union coupling of 031. Persistence
  window (round 2), same author as 031; joins the next reviewer batch.
- **Tools:** `explore/uc_hu_attack.py` (new; joint pattern-search descent
  on CR_HU/H over measure weights AND coordinate orders, in-regime,
  H ≥ 0.2; deterministic — RNG only in fixed-seed start shapes, seeds
  917000+; witnesses saved per row; checkpoint `data/hu_attack.json`);
  `explore/uc_hu_certify.py` (new; exact-rational certification, each
  kit passing ALONE per 029's standard; checkpoint
  `data/hu_certify.json`). Reproduce: run both.
- **Sources:** none.

## Approach

031 posed (HU-TAX) for every coordinate order but attacked mostly the
identity order and left its floors float-level. Here the adversary gets
the order as a degree of freedom: the objective is min over an order
pool (identity + greedy-adjacent-swap descents; full 24-order
enumeration at n = 4 in the follow-up traces) of CR_HU/H, minimized
jointly with pattern-search weight moves, atom additions (δ∅, full set)
and deletions, from every recorded hard instance plus near-cap products
and eight fixed-seed hostile shapes. Caps 0.38271 (working) and 0.45
(probing past ψ). Certification upgrades the sharpest floors — and, as
it turned out, the kill — from float to certified, single-kit.

## What was done

**A. Campaign at the working cap (0.38271): no violation.** Seventeen
starts, all floors positive; global floor CR_HU/H = **+0.04198**
(at a fixed-seed start, n = 4), sitting just above the product-extremal
constant (1−h(0.38271))/h(0.38271) = +0.04174. The order adversary never
produced a meaningfully lower value than the weight adversary alone
(order pools stayed at 1–3 orders; greedy swaps rarely beat identity).
Per-start floors +0.042…+0.215 in `data/hu_attack.json`.

**B. Campaign at cap 0.45: an order-dependent kill.** The descent from
the random0 endpoint reached CR_HU/H = **−0.00319** — and the witness is
clean: n = 4, six atoms

    μ = {1110: .13062, 1101: .13171, 0010: .14591,
         1011: .14278, 0101: .17541, 0000: .27357}
    max marginal 0.44990, H(μ) = 2.5270

with, over the **full 24-order enumeration** (independent evaluator,
025-skeptic stack): exactly **1/24 orders negative** — order (0,1,3,2)
gives CR = −0.00805 — while the identity order gives +0.0148 and the
best order +0.2567 (ratio +0.1016). So at cap 0.45 the *for-all-orders*
positivity of HU is false, while the identity-order and best-order forms
survive on the same witness by wide margins.

**C. The crossing.** Full-order-min descents from the same start at
intermediate caps: floors +0.0426 / +0.0278 / +0.0282 / +0.0299 /
+0.0049 / +0.0009 / −0.0032 at caps 0.39…0.45 — the for-all-orders
floor crosses zero in (0.44, 0.45] on this attack family (floors are
upper bounds on the true inf; the true crossing could sit lower).
Notably far above both ψ = 0.38197 and the working cap.

**D. Certification (single-kit, each kit alone; kit A = 022
digit-extraction, kit B = 016 atanh-series; instances rationalized by
limit_denominator(10⁷), ≤ 1e-7 total variation from the committed
floats; the HU cell tree is exact Fractions, so CR_HU is a finite
Σ c·log₂(r)):

    windowkill (n=4): CR_HU ∈ [+3.089728272e-1, ·]  certified > 0,
                      CR_HU − H/25 certified > +2.02e-1,  H > 2.6751
    mmabskill  (n=5): CR_HU ∈ [+8.223788945e-2, ·]  certified > 0,
                      CR_HU − H/25 certified > +5.527e-3, H > 1.9178
    033 witness, order (0,1,3,2), cap 9/20:
                      CR_HU ∈ [·, −8.048433087e-3]  certified < 0,
                      max marginal < 9/20 exact,     H > 2.5269

  every line under BOTH kits separately (agreement is reported as a
  consistency check, not used for soundness, per 029). The sharpest
  G-gen floor of record is now **certified ≥ H/25 = 0.04·H** at its
  instance, and the order-kill is **certified**, not float.

**E. Witness dissection (lead 1, executed).** Per-coordinate ledgers of
h(z) − (h(x)+h(y))/2 (the HU-notax decomposition, 031):

    identity:  +0.0073  +0.0443  +0.0072  −0.0440   → +0.0148
    (0,1,3,2): +0.0073  +0.0443  −0.0596   0.0000   → −0.0080

The witness has a hidden functional dependency — bit 2 is deterministic
given the other three coordinates — and the bad order reveals it LAST,
where all 30 of its cells contribute exactly zero (a deterministic
conditional has h(z) = (h(x)+h(y))/2 = 0 in every cell): the position is
surplus-sterile. Meanwhile the deficit-prone coordinate (bit 3, which is
genuinely ambiguous given the rest) moves to slot 2 where its
conditional spread — and hence its deficit — is maximal (−0.088 vs
−0.076 in identity). The order kill is exactly: waste a slot on a
determined coordinate, and pay the ambiguous coordinate's deficit
without the banked surplus. Canonical-order heuristic this suggests:
reveal (near-)deterministic coordinates first.

## Outcome

- **REFUTED (certified, at the stated instance): the for-all-orders
  form of HU positivity — hence of (HU-TAX) — at marginal cap 9/20.**
  One explicit witness, one bad order out of 24, CR certified in
  [·, −8.048e-3].
- **EVIDENCE (survival at the working cap):** zero violations across
  the seventeen-start joint (μ, order) campaign at cap 0.38271; global
  floor +0.04198 tracking the product extremal +0.04174 from above.
- **VERIFIED (certified): the working-cap floor values** CR_HU > 0 and
  CR_HU ≥ H/25 at the two sharpest rationalized floor instances.
- **Not claimed:** anything about the true inf below the descent floors;
  identity-order or best-order violations anywhere (none found, incl. at
  cap 0.45); the crossing's exact location (bracketed (0.44, 0.45] by
  one family only); any proof. (HU-TAX) remains SPECULATION and now
  **must be re-posed with an order convention** — for-all-orders is dead
  above ≈ 0.44, best-order and identity-order forms are unrefuted
  everywhere tested.

## Why it failed / what survived

The for-all-orders quantifier was the weakest joint in 031's statement,
and the adversary found it — but only well past ψ: the order degree of
freedom becomes load-bearing somewhere in (0.44, 0.45), while at the
working cap every order of every attacked instance stayed positive.
Structure worth keeping: (i) the kill needed BOTH weapons — the witness
μ is a weight-descent product AND only 1 of its 24 orders fails — so
order-robustness degrades gradually, not catastrophically; (ii) the
identity/best-order forms survive at 0.45 with margins ~30× the kill's
magnitude, so the recipe-facing conjecture (the recipe picks its order)
is untouched; (iii) certification of a *negative* CR value via the exact
cell tree worked unchanged — the machinery is sign-agnostic.

Consequence for the queue's (b2)/(b): the proof target becomes
"CR_HU ≥ 0 (then ≥ c·H) in-regime **for every order**" at p̄ ≤ 0.38271
— now known to be the strongest order-quantifier that can possibly hold
past 0.44, and false there — or the best-order form at higher p̄. The
averaging inequality of 031 must, above ≈ 0.44, use something the bad
orders lack: the deficit-cell mass balance is order-dependent, and the
witness gives the first concrete example where one order's deficit cells
outweigh its surplus cells. Dissecting exactly which prefix tree does it
is the sharpest available probe into the averaging problem.

## Leads generated

1. **Dissect the witness**: per-coordinate, per-cell surplus/deficit
   tables for the bad order vs the identity order on the same μ — the
   minimal concrete instance of the averaging problem (031 §Why).
2. **(b2) restated**: prove all-orders CR_HU ≥ 0 for marginals ≤
   0.38271 (consistent with everything known); separately, best-order
   positivity at caps up to 0.49 (unrefuted).
3. **Map the crossing properly**: more start families and n ≥ 5 full
   enumerations to tighten (0.44, 0.45] and test whether the threshold
   is n-dependent or converges.
4. Canonical-order rule for the recipe (031 lead 3) — now mandatory
   above ≈ 0.44 rather than optional.

## References

- This repo: 031 (the coupling, (HU-TAX), cell taxonomy), 029
  (single-kit standard), 026/022/016 (kits), 023 (floor instances),
  `data/hu_attack.json`, `data/hu_certify.json`.
- No external sources.
