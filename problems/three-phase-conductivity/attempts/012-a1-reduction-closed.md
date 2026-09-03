# 012 — The last reduction: a₁ = rΘ closed by hand, so the construction is a theorem

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (closes 011 lead 1)
- **Type:** derivation
- **Tools:** `explore/tp_attain.py`. Tests: `tests/test_three_phase_attain.py`.
- **Sources:** this repo's 002, 004, 010, 011.

## Approach

011 derived a₁ = (1 − ρ₃)/(y − ρ₃) = σ₁(σ₃+σ₁−c)/(c(σ₃−σ₁)) from tangential
continuity, but its reduction to rΘ was still verified instance-wise. The whole
reduction is one identity about c = σ\* + σ₁ at the attaining fractions, and it
turns out to collapse in four lines.

## What was done

**Claim.** At f₁ = m₁₁ = 2Θ·r(1−r), f₂ = r², f₃ = 1 − f₁ − f₂, with
c = HS_lo + σ₁ and K = (σ₃−σ₂)/(σ₂+σ₁),

    (σ₃+σ₁)/c = 1 + rK.

*Proof.* Expand (σ₃+σ₁)/c = (σ₃+σ₁)·Σ_i f_i/(σ_i+σ₁) and pair the terms with
the pieces of f₃ = 1 − f₁ − f₂:

- The f₁ term against −f₁ gives, using m₁₁ = 2Θr(1−r) and
  Θ = σ₁K/(σ₃−σ₁),

      K·r(1−r)·[(σ₃+σ₁) − 2σ₁]/(σ₃−σ₁) = K·r(1−r),

  since (σ₃+σ₁) − 2σ₁ = σ₃−σ₁.
- The f₂ term against −f₂ gives

      r²·[(σ₃+σ₁) − (σ₂+σ₁)]/(σ₂+σ₁) = K·r².

Summing, K·r(1−r) + K·r² = K·r, and the remaining constant is 1. ∎

Substituting into 011's expression,

    a₁ = σ₁(σ₃+σ₁−c)/(c(σ₃−σ₁)) = σ₁·(rK)/(σ₃−σ₁) = r·Θ,

since Θ = σ₁K/(σ₃−σ₁). **a₁ = rΘ, derived.**

Verified: the identity holds on **15 571** random exact (σ, r) instances with
zero failures, and the resulting a₁ equals rΘ exactly at every attaining
instance tried.

## Outcome

`VERIFIED`. Combined with 010 (the field pattern is forced) and 011 (both
parameters follow from tangential continuity), **the attaining construction of
004 is now derived end to end** rather than fitted: the fields are forced by the
attainment condition plus interface continuity, a₅ = 1 − Θ and a₁ = rΘ follow,
and a₃ is a redundant encoding (005). Nothing in it rests on pattern-matching
any more.

**Not claimed.** That the construction is *optimal* or minimal-rank; 005 showed
it is rank 4 and 009's rank-3 sweep found nothing better at m₁₁, but neither is
a proof. Nothing new about the escalation of 006/009.

## Why it failed / what survived

Nothing failed. Worth recording as a small piece of method: the reduction that
had resisted two records collapsed once the terms were paired against the
pieces of f₃ = 1 − f₁ − f₂ rather than expanded independently. The grouping was
the whole difficulty.

## Leads generated

1. **Apply the same derivation to 007's coated-T² structure** in ℚ(√105), whose
   parameters were solved rather than derived. The interface argument of 011
   should transfer verbatim, and the extra step is the coating fraction.
2. **State the whole thing as one theorem** in a single record: given σ₁<σ₂<σ₃
   and r ∈ (0,1), the explicit rank-4 laminate attains the HS lower bound at
   f₁ = 2Θr(1−r), f₂ = r². All the pieces now exist across 004, 005, 010, 011
   and this record; assembling them is a writing task, not a research one.
3. **Carried from 010:** does (4.26)-active plus correct fractions *characterise*
   attainment? That would turn attainability into a local field criterion.

## References

- This repo: `attempts/002`, `attempts/004`, `attempts/005`, `attempts/010`,
  `attempts/011`
