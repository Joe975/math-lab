# 011 — The construction's parameters derived: a₅ = 1−Θ and a₁ = rΘ from tangential continuity

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (closes the open item of 005, 007 and 010 lead 1)
- **Type:** derivation
- **Tools:** `explore/tp_attain.py`, `explore/tp_fields.py`. Tests:
  `tests/test_three_phase_attain.py`.
- **Sources:** this repo's 002 (attainment field conditions), 004 (the
  construction), 005 (the skeptic pass that flagged the parameters as fitted),
  010 (the forced field pattern).

## Approach

004 found its parameters by pattern-matching solved numerics. 005 confirmed
they were correct and provably equal to an independently derived closed form,
but flagged that neither was *derived*; 007 and 010 carried the same item
forward. 010 supplied the missing ingredient by showing the attainment fields
are forced. With the field values known, the parameters follow from bookkeeping
at two interfaces — no fitting, no search.

## What was done

Write c = σ\* + σ₁, ρ_i = c/(σ_i+σ₁), and y = σ₃ρ₃/σ₁ (010: phase 1 carries the
mirror pair diag(ρ₃, y), diag(y, ρ₃); phases 2 and 3 are isotropic at ρ_i).
Recall the tree, in 005's collapsed rank-4 form:

    A    = laminate(p1 at a₁, p3),  normal e₂
    D    = laminate(p3 at a₅, p1),  normal e₁
    Y    = laminate(D  at w,  p2),  normal e₂
    root = laminate(A  at a₀, Y),   normal e₁

**a₅, from the interface inside Y.** In `Y`, laminated along e₂, the e₁
components are *tangential* and therefore equal across the interface. Phase 2
is isotropic at ρ₂, so `D`'s e₁ component must equal ρ₂. `D` is laminated along
e₁, so its e₁ component is the volume average a₅ρ₃ + (1−a₅)y. Hence

    a₅ = (y − ρ₂)/(y − ρ₃).

Reducing: y − ρ₂ = c(σ₂σ₃ − σ₁²)/(σ₁(σ₃+σ₁)(σ₂+σ₁)) and
y − ρ₃ = c(σ₃ − σ₁)/(σ₁(σ₃+σ₁)), so

    a₅ = (σ₂σ₃ − σ₁²)/((σ₂+σ₁)(σ₃−σ₁)) = 1 − Θ,

which is exactly 004's formula and 005's identity. **Derived.**

**a₁, from the root interface.** At the root, laminated along e₁, the e₂
components are tangential and equal; the applied field is Z₀ = I, so that
shared value is 1. `A` is laminated along e₂, so its e₂ component is
a₁y + (1−a₁)ρ₃, giving

    a₁ = (1 − ρ₃)/(y − ρ₃) = σ₁(σ₃+σ₁−c)/(c(σ₃−σ₁)),

which equals rΘ exactly at the attaining fractions. **Derived.**

Verified: both derived expressions equal the built parameters exactly, in ℚ, at
six (σ, r) instances spanning integer and rational conductivities; the a₅
reduction was additionally checked on 19 471 random exact (σ, c) instances with
zero failures.

## Outcome

`VERIFIED`. The construction of 004 is no longer a fitted formula: given 010's
forced field pattern, both non-trivial parameters follow from tangential
continuity at two interfaces, and a₅ reduces algebraically to 1 − Θ. Scope: the
topology of 004/005 as stated, at the instances listed.

**Not claimed.** A closed-form proof that a₁ = rΘ *as a function of r* — the
derivation gives a₁ = (1−ρ₃)/(y−ρ₃), and the reduction to rΘ uses σ\* = HS_lo at
the attaining fractions, verified instance-wise rather than reduced by hand.
Nothing new about the escalation of 006/009.

## Why it failed / what survived

Nothing failed. The item that three separate records carried as open turned out
to need only the field values plus the observation that a lamination's
tangential component is shared. The general lesson, already noted in 010 and
reinforced here: **when a construction's parameters look fitted, ask what the
interfaces force**. In this problem the interfaces determined everything except
one genuinely free parameter, and 005 had already shown that one to be a
redundant encoding.

## Leads generated

1. **Reduce a₁ = rΘ by hand**, using σ\* = HS_lo at f₁ = m₁₁, f₂ = r². The
   remaining step is showing (σ₃+σ₁)/c = 1 + r(σ₃−σ₂)/(σ₂+σ₁).
2. **Run the same derivation on 007's coated-T² structure** in ℚ(√105) — its
   parameters were also solved rather than derived, and the interface argument
   should apply verbatim.
3. **Use the forced pattern as a search filter** (carried from 010): reject any
   candidate whose phase-2 or phase-3 field is not isotropic at ρ_i before
   computing an effective tensor.

## References

- This repo: `attempts/002`, `attempts/004`, `attempts/005`, `attempts/007`,
  `attempts/010`
