# 016 — Breaking the C/B collapse makes things much worse: the family curve survives its first enrichment

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (executes 015 lead 2)
- **Type:** computational search — **negative result**
- **Tools:** the 006 family builder with the topology altered, exact isotropy
  bisection, two-parameter minimisation refined to 1e−8.
- **Sources:** this repo's 005, 013, 015.

## Approach

015 found that below m₁₁ the 006 family's minimum sits at a stable fraction
c(σ) ≈ 0.53–0.69 of the way from HS_lo to B2, and named as its most
informative next step whether a **richer class drops below that curve**. If one
does, the intermediate-bound reading weakens and the conflict with B2 grows.

The cheapest genuine enrichment comes from 005's structural finding. In the
family, nodes C and B laminate along the *same* normal e₂ in sequence and
therefore collapse, which is exactly why the parameter a₃ is redundant. Giving
node B the orthogonal normal e₁ breaks that collapse and turns a₃ into a
genuine degree of freedom — a strictly richer three-parameter class containing
no member of the original family.

## What was done

Same volume constraints, re-solved for the altered topology (a₄ from the
phase-2 constraint, a₁ from the phase-1 constraint, a₅ by exact isotropy
bisection), minimised over (a₀, a₃) by grid plus refinement, at σ = (1,2,5),
m₂ = 1/4:

| m₁ | HS_lo | B2 | original ratio | enriched min | enriched ratio |
|---|---|---|---|---|---|
| 0.1200 | 3.0268456376 | 3.0271504085 | 0.531326 | 3.0469217303 | **65.87** |
| 0.1100 | 3.0816326531 | 3.0845383760 | 0.528734 | 3.1011433028 | **6.71** |

The enriched class is **dramatically worse**, landing far *above* B2 rather than
below it — a ratio of 66 where the original family gives 0.53.

## Outcome

`REFUTED` as an improvement, for this enrichment: giving node B the orthogonal
normal does not beat the 015 curve, it loses to it by one to two orders of
magnitude in the gap. Scope: σ = (1,2,5), m₂ = 1/4, m₁ = 0.12 and 0.11, this
one topology change.

The 015 curve survives its first enrichment attempt.

**Not claimed.** That no richer class beats it — this is one enrichment out of
many, and 015 lead 2 stays open for non-axis normals, rank 6, and other
topologies.

## Why it failed / what survived

The interesting part is *why* it fails so badly. 005 showed the C/B collapse
makes a₃ redundant, and it was natural to read a redundant parameter as wasted
freedom worth recovering. The opposite is true: **the collapse is the right
structure, not a limitation.** Two consecutive same-normal laminations are what
lets the phase-2 material sit correctly relative to the p1|p3 core; forcing the
orthogonal normal destroys that and the value jumps by orders of magnitude.

That is worth carrying forward as a search heuristic: in this problem, an
apparent redundancy in a good structure signals that the structure has found the
right subspace, and "recovering" the freedom is a way of leaving it.

## Leads generated

1. **015 lead 2, still open**, but redirected: try enrichments that *preserve*
   the collapse — extra laminations elsewhere in the tree, or non-axis normals
   at the root — rather than ones that break it.
2. **Check whether the collapse is forced by the field conditions.** 010 derived
   the attainment fields from interface continuity; the same argument at the C/B
   interfaces may explain why same-normal stacking is required, which would turn
   this empirical dead end into a structural fact.

## References

- This repo: `attempts/005`, `attempts/013`, `attempts/015`
