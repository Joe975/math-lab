# 015 — Below m₁₁ the family's minimum tracks a fixed fraction of B2's claimed improvement

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (extends 006, 009, 013, 014)
- **Type:** exact computation
- **Tools:** the 006 family builder with exact isotropy bisection, plus a
  two-parameter minimisation refined to 1e−9; `explore/tp_cherkaev_bound.py`.
- **Sources:** this repo's 006, 009, 013, 014.

## Approach

Everything so far treated the below-m₁₁ structures qualitatively: they are
*below* B2. Nobody had asked *how far* below, or minimised properly. 013 and 014
classified thousands of structures without ever pushing the minimisation hard.
That is a cheap gap to close, and the answer turns out to be structured.

## What was done

Minimising the 006 family over its two free parameters (coarse grid then local
refinement to 1e−9, exact isotropy bisection at every evaluation, exact
arithmetic throughout), and reporting the position of the minimum between the
two bounds as

    ratio = (family min − HS_lo)/(B2 − HS_lo).

**The ratio is nearly constant in m₁**, at σ = (1,2,5):

| m₁ | m₁₁ − m₁ | ratio |
|---|---|---|
| 0.12490 | 0.00010 | 0.532374 |
| 0.12400 | 0.00100 | 0.532192 |
| 0.12300 | 0.00200 | 0.531984 |
| 0.12000 | 0.00500 | 0.531326 |
| 0.11500 | 0.01000 | 0.530111 |
| 0.11000 | 0.01500 | 0.528734 |
| 0.10000 | 0.02500 | 0.525381 |
| 0.09800 | 0.02700 | 0.524588 |

It rises smoothly to a limit as m₁ → m₁₁ from below. That limit depends on the
conductivities:

| σ | Θ | m₁₁ | limiting ratio c(σ) |
|---|---|---|---|
| (1,2,5) | 0.25000 | 0.125000 | 0.532369 |
| (2,5,11) | 0.19048 | 0.095238 | 0.532696 |
| (1,3,7) | 0.16667 | 0.083333 | 0.576304 |
| (1,4,9) | 0.12500 | 0.062500 | 0.598509 |
| (1,2,9) | 0.29167 | 0.145833 | 0.690125 |

So the family does not merely dip under B2 — it consistently realises only
about **53% to 69% of B2's claimed improvement over HS**, with a clean
σ-dependent limit strictly between 0 and 1. At σ = (1,3,7) the same limit was
already visible in 014's coarser data (ratios 0.5762 → 0.5685 falling away from
m₁₁), so this is not an artefact of the refinement.

Note the ratio is **not monotone in Θ**: (1,2,9) has the largest Θ and the
largest c, while (1,2,5) and (2,5,11) have quite different Θ and nearly
identical c to three decimals. Whatever governs c, it is not Θ alone.

## Outcome

`EVIDENCE`, scoped to the 006 family, m₂ = 1/4, at the (σ, m₁) listed, exact
arithmetic with the minimisation refined to 1e−9.

- The family's minimum below m₁₁ sits at a nearly m₁-independent fraction of the
  way from HS_lo to B2, approaching a σ-dependent limit c(σ) as m₁ → m₁₁⁻.
- Measured c ranges over 0.532 to 0.690 on five conductivity triples.

`SPECULATION`, labelled: that the sharp bound below m₁₁ is of the form
HS_lo + c(σ)·(B2 − HS_lo), i.e. that B2 overstates the true improvement by a
roughly constant factor of about 1.5 to 1.9 in this region. The evidence is one
family; a richer class could sit lower and change c.

**Not claimed.** That the family's minimum is the true minimum — it is one
topology at rank 5, and 013's exhaustive rank-4 sweep could not even reach B2,
so richer classes are not ruled out below. No closed form for c(σ); five data
points and no fit that survives them.

## Why it failed / what survived

Nothing failed. This converts the escalation from a qualitative oddity into a
quantitative one: the disagreement with B2 is not marginal or accidental but
tracks a stable curve, which is what one would expect if a genuine intermediate
bound sits between HS and B2 in this region and the family is近 attaining it.

Method note: three records classified 23 000 structures without once asking how
far below the bound the best one sits. **Classification is not measurement** —
when a search's purpose shifts from "does any exist" to "what is going on",
switch from enumerating to optimising.

## Leads generated

1. **Find a closed form for c(σ)**, or refute that one exists, with more triples.
   It is not Θ alone; try the invariants (σ₂σ₃−σ₁²)/((σ₂+σ₁)(σ₃−σ₁)) = 1−Θ, the
   ratios σ₂/σ₁ and σ₃/σ₁, and the Keller–Dykhne dual triple.
2. **Push a richer class below the family's minimum.** If a rank-6 or non-axis
   class drops materially below c(σ), the intermediate-bound reading weakens and
   the conflict with B2 grows. This is the single most informative next
   computation.
3. **Check c against Nesi's bound** if a transcribable form is ever found: an
   intermediate bound is exactly what Nesi's is supposed to be, and it would be
   a natural candidate for what the family is tracking.
4. Carried, unchanged: re-derive Remark 4.6, and the note to the author.

## References

- This repo: `attempts/006`, `attempts/009`, `attempts/013`, `attempts/014`
