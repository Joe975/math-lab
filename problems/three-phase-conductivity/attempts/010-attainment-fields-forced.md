# 010 — The attainment field pattern is forced, and it makes Cherkaev's (4.26) active by algebra

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (closes an open item of 005/007; explains 009's key observation)
- **Type:** derivation + exact verification
- **Tools:** `explore/tp_fields.py`, `explore/tp_attain.py`,
  `explore/tp_constrained_hunt.py`. Tests: `tests/test_three_phase_fields.py`.
- **Sources:** A. Cherkaev, preprint of Mech. Mater. 41 (2009) 411–433, §2.3
  and §4.4 (4.25)–(4.26) [T]; this repo's 002, 004, 005, 007, 009.

## Approach

004's construction was found by pattern-matching solved numerics, and 005 and
007 both left the same item open: derive its parameters from the field
conditions rather than fitting them. Separately, 009 observed that the
attaining structure makes Cherkaev's constraint (4.26) *exactly* active and
could not say whether that was a coincidence.

Both fall out of one calculation: work out what the attainment fields must be,
using interface continuity, and see how much freedom is left. The answer is
none, and (4.26) drops out as an identity.

## What was done

### 1. The attainment field pattern is fully forced

Work in Cherkaev's variables: Z = ∇u for a pair of potentials, applied
Z₀ = I. Write c = σ\* + σ₁ with σ\* = HS_lo, and ρ_i = c/(σ_i + σ₁).

- **Phases 2 and 3 are isotropic at ρ_i.** This is 002's attainment condition
  (uniform target field in every phase of nonvanishing polarization), so
  Z = ρ₂I and Z = ρ₃I there.
- **Phase 1's field is then forced.** Phase 1 borders phase 3 across a
  lamination normal to e₂. Tangential continuity makes its e₁ component equal
  phase 3's, namely ρ₃. Normal-flux continuity gives σ₁·y = σ₃·ρ₃, so its e₂
  component is

      y = σ₃ρ₃/σ₁.

  By symmetry the mirrored leaf carries diag(y, ρ₃). So phase 1 carries exactly
  the **mirror pair** diag(ρ₃, y), diag(y, ρ₃) — no freedom at all.
- **Consistency is automatic.** Phase 1's volume-mean must be ρ₁ = c/(2σ₁), so
  its smaller component must be 2ρ₁ − y. And

      2ρ₁ − y = c/σ₁ − σ₃c/(σ₁(σ₃+σ₁)) = (c/σ₁)·σ₁/(σ₃+σ₁) = ρ₃,

  an identity in σ and c. The pattern closes on itself for every σ\*.

Verified exactly against the built structures: at five (σ, r) instances the
observed leaf fields are precisely {(ρ₃, y), (y, ρ₃)} for phase 1 and {(ρ_i, ρ_i)}
for phases 2 and 3, in ℚ. Checked also on 19 623 random exact (σ, c) instances
for the two identities alone, zero failures.

### 2. (4.26) is active by algebra, not coincidence

With that pattern, phase 1 has

    S₁ = y + ρ₃,   D₁ = y − ρ₃,   and   ς_N = 2ρ₃ (phase 3's S).

Therefore

    S₁ − ς_N = (y + ρ₃) − 2ρ₃ = y − ρ₃ = D₁,

so **D₁² = (S₁ − ς_N)² exactly**. Cherkaev's constraint (4.26), D² ≤ (S − ς_N)²,
is *exactly active* at any attainment field pattern, for every σ and every σ\*.

That converts 009's slack-of-zero observation from an empirical coincidence
into an identity, and it explains why his framework describes the attaining
object so precisely: the equality case of his constraint and the equality case
of the bound are the same configuration.

### 3. A falsification run for 009

009's reading — that B2 separates on (4.26) — predicts no structure both
respects (4.26) and lies below B2. Extending the test beyond 009's single
family, an exhaustive axis-normal rank-3 sweep at σ = (1,2,5), m₂ = 1/4,
m₁ = 3/25 (`explore/tp_constrained_hunt.py`, 621 screened survivors):

| below B2 | respecting (4.26) | **both** |
|---|---|---|
| 0 | 432 | **0** |

No counterexample. The minimum among respecting structures is 3.0534786940
against B2 = 3.0271504085, a margin of +2.6e−2. A rank-4 sweep was launched and
its result is not in this record.

## Outcome

- `VERIFIED` (derivation plus exact checks): the attainment field pattern is
  forced by interface continuity and the phase-2/3 uniformity condition — phase
  1 carries the mirror pair with y = σ₃ρ₃/σ₁, and 2ρ₁ − y = ρ₃ identically.
  Matches the built structures exactly at five instances; identities hold on
  19 623 random exact instances.
- `VERIFIED` (two-line algebra): D₁² = (S₁ − ς_N)² at any attainment pattern, so
  (4.26) is exactly active. 009's observation is an identity.
- `EVIDENCE` (scope: axis-normal rank-3, one grid point, 621 survivors): no
  structure both respects (4.26) and beats B2, consistent with 009.

**Not claimed.** This derives the *field pattern*, not yet the lamination
*fractions*: a₁ = rΘ and a₅ = 1−Θ still follow from fraction bookkeeping that is
verified but not derived here, so 005's open item is narrowed rather than
closed. Nothing new about whether (4.26) is necessary at the optimum, which
remains the decisive question for 009's reading.

## Why it failed / what survived

Nothing failed. The pleasing part is that the two loose ends turned out to be
the same calculation: the reason 004's parameters looked like magic numbers is
that the fields leave no freedom, and the reason (4.26) sits exactly on its
boundary is that the same pattern forces it.

Method note: the derivation needed nothing beyond tangential and normal-flux
continuity at a single p1|p3 interface. When a construction's parameters look
fitted, ask what the *interfaces* force before trying to fit anything.

## Leads generated

1. **Finish the fraction bookkeeping** to derive a₁ = rΘ and a₅ = 1−Θ from the
   forced pattern, closing 005's open item. The field values are now known, so
   this is bookkeeping over the tree, not a search.
2. **Complete the falsification sweep** at rank 4 and rank 5, and at the other
   two m₁ values. 009's reading stands or falls on it.
3. **Does (4.26)-active characterise attainment?** §1–2 show attainment forces
   (4.26) active. Is the converse true within a suitable class — does an
   isotropic structure with (4.26) exactly active and the right fractions have
   to attain? If so, attainability gets a purely local field criterion.
4. **Use the forced pattern as a search filter.** Any structure whose phase-2 or
   phase-3 field is not isotropic at ρ_i cannot attain, so a screen can reject
   candidates before computing an effective tensor at all.

## References

- A. Cherkaev, preprint `math.utah.edu/~cherk/publ/newbounds7.pdf`, §2.3, §4.4 [T]
- This repo: `attempts/002`, `attempts/004`, `attempts/005`, `attempts/007`,
  `attempts/009`
