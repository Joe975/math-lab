# 037 — Rollout-order HU: non-greedy passes the discriminating test, survives to cap 0.499, certified

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** construction + adversarial attack + exact certification
  (035 leads 1 and 3; the 036 tie-break design note folded in).
- **Tools:** `explore/uc_hu_order2.py` (new; the rules, the witness and
  endpoint batteries, fresh descents, the best-order ceiling;
  deterministic — RNG only in start shapes via `uc_hu_attack.starts`;
  checkpoint `data/hu_order2.json`); `explore/uc_hu_order2_skeptic.py`
  (new; nats-based rule re-implementations with atom-index-partition
  cells, every CR through the pre-existing independent
  `permute_mu → half_union_pairs → cr_eval` stack; exit 0);
  `explore/uc_hu_order2_certify.py` (new; each kit alone per 029, plus
  a 60-digit rollout fixed-point check on the rationalized measure —
  new to the certification pattern; checkpoint
  `data/hu_order2_certify.json`). Reproduce: run the three in that
  order.
- **Sources:** none.

## Approach

035 killed the canonical (greedy conditional-entropy) order at cap
0.49 and diagnosed the mechanism: greedy banks predictability early,
CR needs surplus cells late. Its lead 1 asks for rules that score the
surplus/deficit ledger directly, with the n=5 kill — where canonical
picks the worst of 120 orders — as the discriminating test.

Since H(μ) telescopes along the same revelation order, CR decomposes
per cell as

    CR = Σ_k Σ_cells w(a,b) · [ h(z) − (h(x) + h(y))/2 ],

so s(x,y) = h(z) − (h(x)+h(y))/2 is the cell's ledger entry (z the HU
clip). Rules tested — each reads only (μ, revealed set), so each is a
genuine total coupling:

- **can** — 035's canonical rule (baseline);
- **canst** — canonical, exact ties broken by step surplus (036's S5
  note: lowest-index tie-breaks can cost value);
- **surp** — greedy maximum immediate step surplus: the ledger, one
  step at a time;
- **roll** — rollout: pick the coordinate maximising **full CR** with
  the remainder completed canonically — the ledger scored to the end;
  non-greedy in exactly the sense 035 lead 1 asks for. O(n²) full
  evaluations per order derivation, same complexity class as one CR
  evaluation, so still recipe-usable.

Why these rather than a learned or exhaustively-optimized order: the
point is to test 035's mechanism claim — if the failure is *greed*,
surp (greedy in the right currency) should still die and roll should
survive; if the failure is the *currency* (entropy vs surplus), surp
should survive. The best-order oracle (035 lead 3) is measured
alongside as the ceiling any rule can hope for.

## What was done

**A. The three record witnesses** (033 order-kill; both 035 cap-0.49
kills), all n! orders enumerated alongside:

    witness      worst CR   best CR   can        canst      surp       roll
    033witness   −0.00805   +0.25670  +0.25628   +0.25670   +0.25670   +0.25670
    kill_n4      −0.03959   +0.15075  −0.03795   −0.03795   −0.03096   +0.15075
    kill_n5      −0.02889   +0.03495  −0.02889   −0.02889   −0.02889   +0.03495

  **The discriminating test splits exactly as 035's mechanism claim
  predicts: surp picks the worst of 120 orders on kill_n5, same as
  canonical — greed is the failure, not the currency — while roll
  lands on the best order of all three witnesses.** (canst fixes the
  033-witness tie flagged in 036, +0.25628 → +0.25670, but is
  identical to can on the tie-free kills, as it must be.)

**B. All 51 endpoints of 035's canonical descents, re-scored:** roll
is positive on every one (worst +0.041983 / +0.007920 / **+0.000679**
at caps 0.38271 / 0.45 / 0.49); surp still goes negative on the two
0.49 kills.

**C. Fresh weight-descents against the rules themselves** (17 starts
per cap, order re-derived by the rule at every candidate; n = 6 starts
skipped for roll, budget):

    surp  cap 0.49 : floor −0.012831, 2 violations   → surp is DEAD
    surp  cap 0.497: floor −0.020503, 2 violations
    roll  cap 0.49 : floor +0.000679, 0 violations
    roll  cap 0.497: floor +0.000037, 0 violations

**C2. Hostile seeding against roll** — descents seeded at the
canonical kills, the surp kills, and the sharpest C/D endpoints:

    cap 0.49 : floor +0.000679   [product extremal +0.000289]
    cap 0.497: floor +0.000037   [product extremal +0.000026]
    cap 0.499: floor +0.0000030  [product extremal +0.0000029]

  Zero violations anywhere; the attacked floor tracks the
  product-extremal constant (1−h(cap))/h(cap) from above at every cap,
  pinching onto it as cap → 1/2 — the (HU-TAX) shape, now visible in
  the region where canonical is dead.

**D. The best-order ceiling** (035 lead 3; oracle over all n! orders,
n ≤ 5 starts): descents against the oracle stall at +0.000694 (cap
0.49) and +0.000219 (cap 0.497), zero violations — at this scale no
ceiling below 1/2 is visible for what ANY order rule could achieve.

**E. Certification** (each kit alone; rationalizations
`limit_denominator(1e7)`): CR_HU under the rollout order certified
**positive** at five instances — both 035 kill measures (the rescues,
+1.507474125e-1 and +3.495449146e-2), and the sharpest descent floors
at caps 49/100, 497/1000, 499/1000 (+2.714371934e-3, +3.693561440e-5,
**+3.003306123e-6** — the last with max marginal 0.49898). New to the
certification pattern, per 036's review: a **60-digit rollout
fixed-point check on the rationalized measure** — at every step of the
certified sequence, every alternative's rollout score is recomputed at
60 digits and the certified choice confirmed as the argmax (ties →
lowest index). The check caught a real instance before it shipped: at
the 0.499 floor, rationalization makes the four step-0 scores
**exactly** equal (the gap on the unrationalized measure is ~1.6e-10),
so roll(μ_rationalized) ≠ roll(μ_float) there — the certified order
must be derived from the certified measure. 035's canonical
certificates were confirmed safe on this point by 036 after the fact;
this pattern makes it a construction-time check.

## Outcome

- **EVIDENCE / LIVE: rollout-order HU survives everywhere canonical
  died.** Zero violations across the witness battery, all 51 canonical
  endpoints, fresh and hostile-seeded descents at caps 0.49, 0.497 and
  0.499, with certified-positive floors at all three caps and the
  floor pinching onto the product extremal from above. The
  [0.45, 1/2) interval that 035 said "needs a different object" now
  has a candidate: total, deterministic, closed-form per step,
  polynomial order derivation.
- **REFUTED: the greedy-surplus rule** (fresh descents kill it at
  0.49 and 0.497; it picks the worst order on kill_n5 exactly like
  canonical). 035's mechanism statement is sharpened: **non-greediness
  is the load-bearing property, not the scoring currency.**
- **EVIDENCE: the best-order ceiling stays positive** at caps
  0.49/0.497 under descents against the oracle itself (n ≤ 5).
- **Not claimed:** any proof; any n > 6 statement; that roll always
  finds the best order (it did on all three witnesses — see lead 3);
  that the 0.499 survival extends to caps beyond those tested. Floors
  are upper bounds on the true inf. Roll descents skipped the three
  n = 6 starts (covered against surp only).

## Why it failed / what survived

Nothing failed that was supposed to survive; the informative kill is
that **surp dies exactly where canonical does**: at the n=5 kill both
pick the worst of 120 orders while optimizing different one-step
scores. The witness's geometry punishes *any* rule that spends its
early slots on immediate reward — the deficit is structural to the
prefix, not to the score. That is why rollout works: its step-k score
already contains the cost of every later slot. The open question it
sharpens: rollout's completions are canonical, so a proof about roll
inherits the canonical rule's analysis inside a max — awkward; lead 2
asks whether the simpler *best-order* form is provable instead, since
the oracle data says best-order positivity is what actually holds.

## Leads generated

1. **(HU-TAX, rollout form):** re-pose 034's corrected conjecture —
   CR_HU ≥ c*(p̄)·H(μ), p̄ ∈ [1/4, 1/2), c*(p̄) =
   (h(min(2p̄,1))−h(p̄))/h(p̄) — with the rollout order. First
   falsifiable step: n = 6–7 descents with larger supports (this
   record's roll descents are n ≤ 5).
2. **Best-order positivity as the proof target:** the oracle floor is
   positive at every cap tested; "for every in-regime μ SOME order has
   CR_HU ≥ 0" is weaker than any fixed-rule statement, is what part D
   actually measures, and roll is its constructive certificate-finder.
   A proof would give the recipe an existence license with roll as the
   implementation.
3. **Is roll = best-order at small n?** It found the exact best order
   on all three witnesses. Census: random in-regime μ at n = 4–5,
   roll rank vs the full enumeration. If roll is near-best with high
   probability, lead 2's gap between rule and oracle is small in
   practice; a counterexample family would show where one-step
   lookahead with canonical completion still fails.
4. **The descent attractor at high caps is (near-)product:** the 0.499
   endpoint rationalizes to an exactly-exchangeable measure (the
   step-0 scores tie 4-ways). Characterize the binding adversary as
   cap → 1/2; if it is always the product boundary, the equality
   structure of (HU-TAX) is doing the resisting, which is evidence the
   conjectured extremal is right.

## References

- This repo: 035 (the discriminating test and both kill witnesses),
  036 (tie-break design note; the fixed-point-check idea), 033/031
  (HU, the ledger framing), 030 (HU-mix), 034 (the corrected
  conjecture constant), 029 (single-kit standard), 026/022/016 (kits).
  `data/hu_order2.json`, `data/hu_order2_certify.json`.
- No external sources.
