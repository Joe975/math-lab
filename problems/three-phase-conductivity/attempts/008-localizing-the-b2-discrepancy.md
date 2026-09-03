# 008 — The B2 discrepancy localized to one condition: constancy of the field in the most conducting phase

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (extends 006; reads the primary source)
- **Type:** obstruction analysis
- **Tools:** `explore/tp_fields.py` (per-leaf fields, exact),
  `explore/tp_below_m11.py`, `explore/tp_attain.py`. Data:
  `data/attained/below-m11-m1-*.json`, `data/attained/rank5-m11-1-2-5.json`,
  `data/cherkaev-primary-excerpts.md`.
- **Sources:** A. Cherkaev, preprint of Mech. Mater. 41 (2009) 411–433,
  `math.utah.edu/~cherk/publ/newbounds7.pdf`, §2.3, §4.1–4.4, §5.1 [T],
  text-extracted and read here.

## Approach

006 left an unresolved conflict: laminates below m₁₁ sitting under Cherkaev's
B2 bound, with transcription, scope, closed-form simplification and a bug in
our structures all ruled out. The adjudicator named the one unexamined place —
the field constraints the bound layers on top of plain translation — and then
ran out of session. This attempt goes there and asks the narrow question: **is
there a constraint in the derivation that our structures violate?**

Answering it does not by itself decide who is right, but it converts an
unexplained discrepancy into a named one, which is the deliverable.

## What was done

### 1. The constraints, from the primary text

Cherkaev's §2.3 writes the field matrix Z = ∇u for a *pair* of potentials in a
rotationally invariant basis:

    S = (Z₁₁+Z₂₂)/√2,  D\* = (Z₁₁−Z₂₂)/√2,  D\*\* = (Z₁₂+Z₂₁)/√2,  V = (Z₁₂−Z₂₁)/√2,
    D² = D\*² + D\*\*².

§4.4 then states "the minimal of all sets that satisfy conditions
(4.19)–(4.23)" [T]:

1. **(4.24)** V(x) = 0 everywhere.
2. **(4.25)** "Field in Ω_N is constant and isotropic, e_N = ς_N I/√2. Ω_N
   consists of one point: S_N = ς_N, D = 0." — Ω_N is the **most conducting**
   phase.
3. **(4.26)** For i < N: D² ≤ Σ_i(S) = (S − ς_N)².

§5.1 substitutes Σ₁ = (S₁−S_N)² to get the improved bound, so B2's entire
advantage over plain translation rests on these. The paper is explicit that
they are *optimality* conditions: §4.4 says "The toughest bound corresponds to
the smallest Σ_i ≥ 0" and Remark 4.6 says "Conditions of an optimal contact
coincide with the above conditions (4.25) and (4.26)", obtained by structural
variation.

### 2. Testing them on our two structures

Computed exactly, by applying both orthogonal fields and assembling Z per leaf:

| structure | phase | S | D² | |
|---|---|---|---|---|
| **at m₁₁** (004, attains) | p1 | 4.0 | 7.111 | |
| | p2 | 2.6667 | **0** | |
| | p3 = Ω_N | **1.3333 (single value)** | **0** | (4.25) **holds** |
| **below m₁₁** (006) | p1 | 4.0482, 4.0529 | 7.367, 7.425 | |
| | p2 | 2.6690 | 3.2e−05 | |
| | p3 = Ω_N | **1.3389 and 1.3461** | 5.4e−05, 2.5e−04 | (4.25) **VIOLATED** |

- **(4.24) holds for both.** V = 0 exactly in every leaf of both structures.
  Checked as exact rational equality, not numerically small.
- **(4.25) holds exactly for the attaining structure** — its most conducting
  phase carries a single, isotropic field value. That is a satisfying
  independent confirmation of Cherkaev's picture: the structure that attains the
  bound is precisely the one obeying his optimality conditions.
- **(4.25) fails for the below-m₁₁ structure.** Its phase-3 field takes two
  distinct values and is not isotropic. Phase 2 also carries a small nonzero D.

### 3. What that does and does not settle

**It names the discrepancy.** Our candidate counterexample violates exactly the
condition Cherkaev derives by structural variation and then imposes to lift the
bound above plain translation.

**It does not dissolve the conflict**, and it is important to say why. The
derivation's logic is: the *optimal* structure satisfies (4.24)–(4.26), so
minimising over the constrained field set still reaches the true minimum, hence
B2 lower-bounds every structure. Our structure need not be optimal — but its
conductivity is *achieved*, so the true minimum is at most our value, which is
below B2. If (4.25) really holds at the optimum, that is impossible. So the
tension is now squarely on (4.25) as an optimality condition, not on anything
else in the chain.

`SPECULATION`, labelled: the most economical explanation is that (4.25) is not
without loss of generality — that admitting a non-constant field in the most
conducting phase permits structures the constrained problem cannot see. Our
below-m₁₁ family would then be an explicit witness. The alternative, still open,
is an error on our side that six independent checks have not caught.

## Outcome

`MAP` — an obstruction analysis that localizes 006's escalation to a single
named condition. No claim is made that Cherkaev's Theorem 7.1 is wrong, and
none that it is right.

**Scope:** two structures, at σ = (1,2,5), m₂ = 1/4, exact. The constraint
statements are [T] from the primary text.

**Not claimed.** That violating (4.25) makes our structure invalid — (4.25) is
an optimality condition, not an admissibility condition, so a legitimate
composite may violate it. That (4.25) is wrong. Anything about attainability
below m₁₁, which remains as 003's coating lemma left it: unreachable by the
coated-T² route, since the required coating fraction exceeds 1 there.

## Why it failed / what survived

The hoped-for clean resolution — "our structure is outside the theorem's
hypotheses, no conflict" — did not arrive, because (4.25) constrains the
*derivation's* admissible fields rather than the *theorem's* stated hypotheses,
and the theorem is stated for every isotropic three-phase composite. What
survived is a much sharper object than 006 had: a specific equation, with a
specific published justification (structural variation, Remark 4.6), that a
specific verified structure violates.

The corroborating half is worth as much as the critical half: the *attaining*
structure of 004 satisfies (4.24) and (4.25) exactly, phase 3 carrying one
isotropic field value. Cherkaev's optimality conditions describe our attaining
object perfectly. That is strong evidence the framework is right where it is
claimed to bite, and isolates the question to whether the conditions are
necessary as well as descriptive.

## Leads generated

1. **Re-derive Remark 4.6's structural-variation argument** and check whether it
   establishes (4.25) as necessary at the optimum or only as a property of a
   particular optimal family. This is now the single decisive step; everything
   else has been eliminated.
2. **Search for a below-m₁₁ structure that satisfies (4.25)** and still beats
   B2. If one exists, the conflict no longer routes through (4.25) at all and
   becomes much more serious. If a constrained search provably cannot beat B2,
   that is strong evidence the bound is right and our structure is somehow
   inadmissible.
3. **Push the below-m₁₁ family lower** over the richer coated-T² class of 007,
   which the 006 family cannot express. A larger margin would make the question
   easier to adjudicate.
4. **Ask the author.** Recorded for the human, not an agent: a short note with
   the explicit structure would settle in one exchange what re-derivation may
   not.

## References

- A. Cherkaev, preprint `math.utah.edu/~cherk/publ/newbounds7.pdf`, §2.3,
  §4.1–4.4, §5.1, Remark 4.6 [T]; excerpts in
  `data/cherkaev-primary-excerpts.md`.
- This repo: `attempts/004`, `attempts/006`, `attempts/007`,
  `data/bound-adjudication.md`.
