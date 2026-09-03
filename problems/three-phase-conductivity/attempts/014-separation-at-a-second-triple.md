# 014 — The (4.26) separation reproduces at a second conductivity triple

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (executes 013 lead 3)
- **Type:** falsification search
- **Tools:** `explore/tp_constrained_hunt.py` (the exact (4.26) classifier),
  dense sweep of the 006 family. Tests: `tests/test_three_phase_fields.py`.
- **Sources:** this repo's 006, 009, 013.

## Approach

Everything supporting 009's reading — that Cherkaev's B2 separates exactly on
his constraint (4.26) — had been computed at the single conductivity triple
σ = (1,2,5). That leaves one cheap way for the reading to be wrong: something
special about that triple. 013 flagged it as its lead 3. This runs the same
test at σ = (1,3,7).

## What was done

σ = (1,3,7), m₂ = 1/4 (so √m₂ = 1/2 and everything stays exact in ℚ), where
Θ = 1/6 and m₁₁ = 1/12 ≈ 0.083333. Dense sweep of the 006 family at three m₁
below m₁₁, classifying every member by the exact sign of (S − ς_N)² − D²:

| m₁ | HS_lo | B2 | tried | below B2 | respecting (4.26) | **both** | min respecting |
|---|---|---|---|---|---|---|---|
| 0.0820 | 4.3475935829 | 4.3476593044 | 785 | 14 | 395 | **0** | 4.3481960875 |
| 0.0780 | 4.3908355795 | 4.3919221648 | 793 | 58 | 328 | **0** | 4.3970367436 |
| 0.0700 | 4.4794520548 | 4.4867170363 | 825 | 176 | 222 | **0** | 4.5041758396 |

**Same pattern, new triple.** 248 structures below B2, every one violating
(4.26); 945 respecting structures, every one above B2; zero exceptions. The
margin of the respecting minimum above B2 again grows as m₁ falls
(+5.4e−4, +5.1e−3, +1.7e−2).

## Outcome

`EVIDENCE`, strengthening 009 and 013. Scope: σ = (1,3,7), m₂ = 1/4, at
m₁ = 0.082, 0.078, 0.070, within the 006 family; float screen then exact
classification.

The separation is **not an artefact of σ = (1,2,5)**. Combined with 013 the
total is over 23 000 structures across two conductivity triples, two exhaustive
axis-normal topology classes and six grid points, with no structure both
respecting (4.26) and beating B2.

**Not claimed.** Still a finite search, still silent on whether (4.26) is
necessary at the optimum — the one question no search can answer.

## Why it failed / what survived

The falsification failed again, at the last cheap place it could have
succeeded. 009's reading is now about as secure as computation can make it, and
the remaining disagreement is entirely about Cherkaev's structural-variation
argument in Remark 4.6.

## Leads generated

1. **Re-derive Remark 4.6.** Now the only thing left that can move this, and it
   has been the answer to that question since 008.
2. **Ask the author** (carried, for the human). Given how sharply the question
   is now posed — an explicit structure, a named condition it violates, and the
   observation that the bound separates exactly on that condition — this is a
   one-paragraph question with a definite answer.
3. Retire further searching on this line: 013 and 014 together give strongly
   diminishing returns.

## References

- This repo: `attempts/006`, `attempts/009`, `attempts/013`
