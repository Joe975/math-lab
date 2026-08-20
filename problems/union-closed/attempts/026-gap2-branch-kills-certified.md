# 026 — The 025 branch kills and repairs, certified in exact arithmetic

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** exact-rational certification (025 lead 3 / queue 1(c)) of the
  branch kills and their repairs, dual-kit; same persistence window as
  025, same author — implementation-level verification, not a
  reviewer-independence pass.
- **Tools:** `explore/uc_branch_certify.py` (new; stdlib only;
  deterministic; ~5 s; checkpoint `data/branch_certify.json`). Imports,
  unmodified, the two previously audited certified-log₂ kits:
  `uc_gap1_skeptic_gemini.certified_log2` (digit extraction; 022-drafted,
  024-audited) and `uc_or_agg_probe.log2_enclosure` (atanh series with
  explicit tail bound; the 016/017 standard). Reproduce:
  `python3 problems/union-closed/explore/uc_branch_certify.py`.
- **Sources:** none.

## Approach

Every load-bearing number in 025 is a finite sum Σ cᵢ·log₂(rᵢ) with
rational cᵢ, rᵢ (entropies of explicitly rational measures under the EM
closed form or profile arithmetic), so certification needs nothing beyond
certified log₂ of rationals — machinery the repo already owns twice over.
Rather than trust either kit alone, every log₂ is evaluated by BOTH kits
and the enclosures intersected: a disjoint intersection would prove one
kit unsound and aborts the run, so the result is certified under either
kit's audit. Coefficient arithmetic is exact `Fraction`s; interval sums
round nothing. Why this rather than a fresh certified-log₂: a third
implementation adds nothing over the intersection of two audited ones,
and 025 lead 1 (fresh-session reviewer) is the real independence step —
this attempt deliberately does not claim it.

## What was done

Nine statements certified, at rational instances chosen inside 025's
failure regions (n = 2 and 6 match 025's A5 float instances exactly as
rationals; the n = 16 sliver point and the block/rescue instances are the
025 constructions with rational parameters). All enclosures have width
≤ 2.5e-32; 98 dual-kit log₂ evaluations, zero kit disagreements.

**C1 — half-mixing branch kill** (P₁ = 16/25, so q = P₁/2 = 8/25):
CR_hm = F(P₁/2) − F(P₁) certified negative, in-regime, H(μ) certified
> 0.968:

    n= 2, ph = 1987/2000 : CR_hm ∈ [−2.267927906e-3, ·]   marg 0.3577
    n= 6, ph = 1997/2000 : CR_hm ∈ [−7.133120856e-3, ·]   marg 0.3595
    n=16, ph = 19993/20000: CR_hm ∈ [−1.514407394e-2, ·]  marg 0.3599

  (Marginal comparisons (1−P₁)·ph < 38271/100000 are exact rational
  inequalities, not enclosures.)

**C2 — repair at the same instances**: CR(q = 1/2) = F(1/2) − F(P₁)
certified positive: +7.307503323e-2 (n = 2), +7.095307791e-2 (n = 6),
+6.744828594e-2 (n = 16).

**C3 — block iid-given-k kill at n = 8**: Gain_block certified negative
on the 025 mixtures 2block {p=0: 31/50, p=19/20: 19/50} and 3block
{p=0: 3/5, p=9/10: 3/10, p=99/100: 1/10}:

    Gain(2block, n=8) ∈ [−7.940026768e-1, ·]   marg 0.3610 exact
    Gain(3block, n=8) ∈ [−1.050431681,   ·]   marg 0.3690 exact

  With 009 fact (i) (CR ≤ Gain for every coupling; 011-verified), this
  certifies CR_block < 0 at these instances. The dependence on fact (i)
  is a proved lemma, not a computation, and is stated here explicitly.

**C4 — ∅-mixing rescue of the 3block instance** (n = 6, q = 1/5,
s ~ (3/4)Bern(9/10) + (1/4)Bern(99/100)): CR = F_gen(1/5) − F_gen(3/5)
certified positive, +6.929688470e-1. This leg leans on **Lemma EM (025;
same-session skeptic-passed, reviewer pass pending)** for CR = Gain; the
certified quantity itself is the entropy difference.

**C5 — regime and nondegeneracy**: at every instance, max marginal
< 38271/100000 exactly and H(μ) certified above 0.968 (half-mixing) /
1.82 (blocks).

**Cross-checks:** all five of 025's float values for these instances lie
within 1e-9 of the certified midpoints (containment is the wrong test —
the engine's binary-float inputs 0.9935 ≠ 1987/2000 sit ~1e-16 off the
rationals, far outside 1e-32 enclosures; the certify script documents
this); the two kits' enclosures intersected at every one of the 98 calls.

## Outcome

**REFUTED (now certified, at the stated instances):** the half-mixing
branch (q = P₁/2) of the 009 recipe has certified CR < 0 at in-regime
genre instances with n = 2, 6, 16 and H(μ) > 0.968; the iid-given-k block
coupling has certified Gain < 0 (hence CR < 0 by 009(i)) at n = 8 on
in-regime mixtures with above-ψ components. 025's kills are no longer
float-level at these witnesses.

**VERIFIED (certified positive, same instances):** the q = 1/2 ∅-mixing
repair on all three sliver instances, and the q = 1/5 ∅-mixing rescue of
the 3block mixture (the latter modulo Lemma EM for the CR = Gain
identification).

**Not claimed:** anything at instances other than the six certified ones
(025's grids and families remain float-level evidence); reviewer-level
independence (same session, same author as 025 — the 025 lead-1 pass
covers both records and remains queued); certification of the EM /
EM-coverage lemmas themselves (they are proofs to review, not numbers to
enclose); any statement about Gap 2's candidate beyond what 025 already
records.

## Why it failed / what survived

Nothing failed. Worth recording: the dual-kit intersection pattern cost
three lines and turns "trust this kit's audit" into "trust either audit"
— it caught nothing today (the kits agreed to 1e-32 across all 98 calls,
which is itself fresh mutual evidence for both), but it is the cheapest
hardening available for any future certification and should be the
default whenever two audited kits exist. Reusable: `uc_branch_certify.py`
(interval sum/scale helpers over Fraction endpoints, the profile-entropy
enclosure for exchangeable mixtures with an ∅ atom, and the
float-vs-rational-input proximity cross-check pattern).

## Leads generated

1. Unchanged from 025: the fresh-session reviewer pass (025 lead 1) now
   covers 025 + 026 together; queue (a') already says this.
2. If a future pass wants the *grids* certified rather than spot
   instances: the G-criterion of EM-coverage is also a rational-log₂ sum,
   so the same machinery certifies G > 0 on rational grid points — only
   worth doing if the reviewer pass finds reason to doubt the float
   sweeps.

## References

- This repo: 025 (the kills, Lemma EM, instance constructions), 022/024
  (kit A and its audits), 016/017 (kit B and its audit), 009/011
  (fact (i), used as a lemma in C3), `data/branch_certify.json`.
- No external sources.
