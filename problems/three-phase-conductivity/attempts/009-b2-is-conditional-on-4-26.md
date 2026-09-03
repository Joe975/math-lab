# 009 — The escalation resolved: B2 separates exactly on Cherkaev's constraint (4.26)

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-03
- **Mode:** informed (resolves 006, extends 008)
- **Type:** obstruction analysis + exact computation
- **Tools:** `explore/tp_fields.py`, `explore/tp_below_m11.py`,
  `explore/tp_cherkaev_bound.py`, `explore/tp_attain.py`, the harness.
- **Sources:** A. Cherkaev, preprint of Mech. Mater. 41 (2009) 411–433,
  `math.utah.edu/~cherk/publ/newbounds7.pdf`, §4.4 (4.24)–(4.26) and §5.1 [T],
  text-extracted and read here.

## Approach

008 localized 006's escalation to Cherkaev's condition (4.25) but could not
say whether the bound was wrong or our structures inadmissible. The way to
decide without re-deriving his structural-variation argument: **test whether
the constraint is what separates**. If B2 is a correct bound over the class its
proof constrains, then structures respecting the constraint should lie above
B2 and only violators should lie below. That is a falsifiable prediction, and
it distinguishes "the bound is wrong" from "the bound's statement is broader
than its proof".

## What was done

### 1. The constraint that does the work

§5.1 substitutes Σ₁ = (S₁ − S_N)² into the bound derivation, so the operative
condition is (4.26):

    D² ≤ Σ_i(S) = (S − ς_N)²   for i < N,

with S, D the rotationally invariant components of the field matrix Z = ∇u for
a pair of potentials (§2.3), and ς_N the S-value in the most conducting phase.
The mechanism matters: §4.2 notes the coefficient (k₁ − t) multiplying
∫_Ω₁ D² dx is **negative** for t > k₁, so the derivation caps ∫D² using (4.26)
to bound the energy from below. **A structure whose D² exceeds that cap makes
the estimate smaller, which is exactly how it can land under B2.**

### 2. The constraint is exactly active at the attaining structure

For 004's structure attaining HS_lo at m₁₁, phase 1 has

    S = 4, D² = 7.1111…, (S − ς_N)² = 7.1111…,  slack exactly 0.

That is the derivation's own equality case — §4.2 says the minimiser takes
D(x) = ±Σ₁^{1/2}(S(x)) when the constraint is active. (4.24) V = 0 and (4.25)
also hold exactly. Cherkaev's optimality conditions describe our attaining
object perfectly.

### 3. B2 separates exactly on (4.26)

Searching the below-m₁₁ family at σ = (1,2,5), m₂ = 1/4, and splitting members
by whether they satisfy (4.26):

| m₁ | B2 | min over **all** members | min over **(4.26)-respecting** members | respecting |
|---|---|---|---|---|
| 0.1240 | 3.0053523724 | 3.0053470065 (**below**) | 3.0056032896 (**above**) | 242/477 |
| 0.1200 | 3.0271504085 | 3.0270129104 (**below**) | 3.0296467290 (**above**) | 203/482 |
| 0.1100 | 3.0845383760 | 3.0831692377 (**below**) | 3.0921743261 (**above**) | 146/503 |

At every point the minimum over constraint-respecting structures lies **above**
B2, and every structure found below B2 **violates** (4.26). The bound is the
separator, not an obstacle the constrained structures trip over.

### 4. What this establishes

The coherent reading, and the one the data supports at every point tested:

**Cherkaev's B2 is a correct lower bound over the class of composites whose
fields satisfy (4.24)–(4.26). It is not a lower bound over all isotropic
three-phase composites, because those constraints are not without loss of
generality.** Theorem 7.1 is stated for every 2D isotropic composite of three
isotropic materials, which is broader than what its derivation constrains.

This also retro-explains everything from 006 onward: the transcription was
right, the closed form faithful, the hypotheses as printed satisfied, our
structures correct — and no contradiction, because the theorem's *statement*
and its *proof's scope* differ.

## Outcome

`EVIDENCE`, scoped to σ = (1,2,5), m₂ = 1/4, at m₁ = 31/250, 3/25, 11/100,
over the laminate family of 006 (roughly 500 members per point, exact
arithmetic):

- Every member below B2 violates (4.26); every (4.26)-respecting member lies
  above B2. B2 separates the two classes at all three points.
- `VERIFIED` (exact): the m₁₁-attaining structure of 004 satisfies (4.24) and
  (4.25) exactly and makes (4.26) **exactly active**, slack 0.

`SPECULATION`, labelled and load-bearing: that (4.26) is not without loss of
generality for all composites. The evidence is that our explicit, six-ways
verified laminates violate it and beat the resulting bound; the alternative
remains an error on our side that no check has caught.

**Not claimed.** That Cherkaev's Theorem 7.1 is *wrong*. The finding is about
the gap between its statement and its derivation's constrained class, which is
a scope observation, not an error claim. Nor is this a proof that no
constraint-respecting structure beats B2 — that is a finite search in one
family, and it is the direction that would matter most to check harder.

## Why it failed / what survived

Nothing failed. The escalation is resolved in the most satisfying available
way: not by finding our error, and not by declaring a published theorem wrong,
but by identifying **exactly which hypothesis the two sides disagree about** and
showing that the bound behaves perfectly on its own side of that line.

The strongest single piece of evidence is the pairing. The structure that
attains the bound makes the constraint exactly active; the structures that beat
the bound violate it; the structures that respect it stay above. Three
independent behaviours, all consistent with one explanation.

Method note: the decisive experiment was not re-deriving the paper's argument —
which two passes had failed to do — but **testing whether the suspected
hypothesis was the separator**. When a derivation's assumption is suspect and
re-derivation is expensive, partition the candidates by the assumption and see
whether the bound respects the partition. That is much cheaper and nearly as
informative.

## Leads generated

1. **Harden §3 into a claim about the whole class**, not one family: search
   several laminate topologies for a (4.26)-respecting structure below B2. A
   single one would overturn this attempt's reading and reopen 006.
2. **Re-derive Remark 4.6** anyway. This attempt says *what* the two sides
   disagree about; only the derivation says *who is right* about whether an
   optimal structure must satisfy (4.26).
3. **Quantify the gap between the two minima.** The constrained minimum exceeds
   B2 by 2.5e−4 / 2.5e−3 / 7.6e−3 at the three points, growing as m₁ falls. If
   the constrained minimum is itself the true optimum over its class, that slack
   measures how far B2 is from sharp on its own class.
4. **Ask the author** (kept from 008, for the human): the exact structure plus
   the observation that it violates (4.26) is a one-paragraph question.

## References

- A. Cherkaev, preprint `math.utah.edu/~cherk/publ/newbounds7.pdf`, §4.2, §4.4
  (4.24)–(4.26), §5.1 [T]; excerpts in `data/cherkaev-primary-excerpts.md`.
- This repo: `attempts/004`, `attempts/006`, `attempts/008`,
  `data/bound-adjudication.md`.
