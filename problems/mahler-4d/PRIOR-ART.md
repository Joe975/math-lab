# Mahler in ℝ⁴ — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None.** This problem was onboarded and has not been worked. There is no prior
art to be informed by, so `blind` and `informed` mode are currently equivalent
here — which makes the first attempts worth running blind, since blind costs
nothing while the record is empty.

## Editorial view of the attack surface

The best exact-arithmetic fit of the three problems onboarded alongside it. For
a rational polytope everything in the question lives in ℚ: the polarity
conversion, the triangulation, both volumes, and the comparison against 32/3.
There is no certification problem to solve first, unlike the operator-theory
side of the lab — `fractions.Fraction` is sufficient, and a C kernel is
available if a census inner loop needs the speed.

What that buys is narrow, and worth being honest about up front. Every Hanner
polytope is a strict *local* minimizer, so no perturbative search near the
conjectured extremals can find anything; a census only has value if it reaches
bodies that are not near a Hanner polytope. Combined with the exponential gap
between the Bourgain–Milman/Kuperberg bounds and the conjecture, the realistic
deliverable is a mapped universe with counts, not an approach to the bound.

Concrete lines, if you want them:

- Self-test: reproduce vol·vol° for the cube, the cross-polytope, ball
  approximants, and a published ℝ³ value with our own tooling. Validates the
  polarity and triangulation code against known answers before anything else is
  claimed.
- Exact census over small symmetric polytopes in ℝ⁴ with vertices in {0, ±1}:
  is anything within ε of 32/3 that is not a Hanner polytope? `EVIDENCE` scoped
  to the generator universe, which must be stated in the record.
- Rational-perturbation check of local minimality at the Hanner polytopes in
  ℝ⁴ — independent confirmation of a published theorem, so its value is harness
  validation rather than progress. Count the Hanner types in ℝ⁴ before
  designing this; there are more of them than in ℝ³, and a plan that assumes
  two is wrong.

Kill/win condition, inline: any symmetric polytope with volume product below
32/3 refutes Mahler outright, so every search here is self-policing — and the
verification contract in `PROBLEM.md` deliberately makes such a claim expensive
to write down, which is the correct cost.

Cross-pollination note for the informed side: Viterbo's conjecture, the
symplectic statement that implies Mahler, was refuted in 2024, but by a
non-symmetric counterexample. Running `/ripple` on that refutation is a cheap
early exercise — the question worth asking is whether anything in the
refutation's mechanism transfers to the symmetric setting this problem lives
in.
