# Mahler's Volume-Product Conjecture in ℝ⁴

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** For a convex body K ⊂ ℝⁿ that is centrally symmetric (K = −K),
let K° = {y : ⟨x, y⟩ ≤ 1 for all x ∈ K} be its polar body, and let the volume
product be 𝒫(K) = vol(K)·vol(K°). Mahler (1939) conjectured

  𝒫(K) ≥ 4ⁿ/n!,

with equality exactly for the Hanner polytopes. The volume product is invariant
under invertible linear maps (vol(TK)·vol((TK)°) = vol(K)·vol(K°)), so the
question is about linear equivalence classes of bodies.

**Hanner polytopes** are the bodies obtained from segments by repeatedly taking
ℓ¹ and ℓ∞ sums; the cube and the cross-polytope are the two extreme cases, and
every Hanner polytope has volume product exactly 4ⁿ/n!. The conjectured minimum
is therefore attained by more than one body up to linear equivalence, which is
what makes the equality case delicate rather than incidental.

For n = 4 the conjectured bound is 4⁴/4! = 32/3.

## Published status

- **n = 2.** Mahler (1939).
- **n = 3.** H. Iriyeh, M. Shibata, *Symmetric Mahler's conjecture for the
  volume product in the three dimensional case*, arXiv:1706.01749, Duke Math.
  J. **169** (2020) — the symmetric case with its equality characterization. A
  streamlined proof emphasising equipartitions is due to
  Fradelizi–Hubard–Meyer–Roldán-Pensado–Zvavitch.
- **n ≥ 4. Open.** This is the frontier.
- **Local minimality.** Nazarov–Petrov–Ryabogin–Zvavitch (2010): the cube is a
  strict local minimizer of the volume product among symmetric convex bodies
  under the Banach–Mazur distance. J. Kim, *Minimal volume product near Hanner
  polytopes*, arXiv:1212.2544 (2014): every Hanner polytope is a strict local
  minimizer. So no counterexample can be found by perturbing a Hanner polytope.
- **Lower bounds.** Bourgain–Milman (1987): 𝒫(K) ≥ cⁿ·4ⁿ/n! for an absolute
  constant c > 0 (the "reverse Santaló inequality"). G. Kuperberg, *From the
  Mahler conjecture to Gauss linking integrals*, GAFA **18** (2008): an
  explicit constant, giving the conjectured bound up to a factor of order
  (π/4)ⁿ. The gap to the conjecture is therefore exponential in n, not
  constant.
- **Non-symmetric case.** The parallel conjecture is 𝒫(K) ≥ (n+1)^{n+1}/(n!)²
  with the simplex extremal; also open in general. Chen–Li–Xi–Xu, *The Mahler
  Conjecture in Three Dimensions*, arXiv:2605.09334 (May 2026), claim the full
  three-dimensional case including the non-symmetric bodies, via a "shadow
  flow" method, together with a new proof of the symmetric n = 3 case. Preprint,
  not yet peer reviewed; it does not address n ≥ 4.
- **Symplectic connection.** Artstein-Avidan–Karasev–Ostrover (2014) showed
  that Viterbo's symplectic capacity conjecture implies Mahler's. Viterbo's
  conjecture was refuted by Haim-Kislev–Ostrover, *A Counterexample to
  Viterbo's Conjecture*, arXiv:2405.16513 (to appear, Annals of Mathematics) —
  the counterexample is not centrally symmetric, and a symmetric variant that
  would still imply Mahler remains open.

## Verification contract

Any claim recorded against this problem must meet the bar in `CONTRIBUTING.md`.
Specifically here:

- Volume products of rational polytopes must be computed **exactly in ℚ**.
  Floating point may propose a candidate; it never supports a claim, because
  the entire question is a comparison against an exact rational threshold.
- A claim must state which representation was used (vertices or facets) and
  that the polarity conversion between them was exact. Polarity exchanges the
  two, so a body given by vertices with V of them yields a polar given by V
  facet inequalities, and the volume of each must be stated with the
  triangulation or decomposition that produced it.
- A **census** claim must state its generator universe precisely — the
  coordinate set, the vertex bound, the symmetry convention — and the counts
  it enumerated, including how linearly equivalent copies were identified
  rather than counted twice. A census is `EVIDENCE` about that universe and
  nothing wider.
- Any candidate with 𝒫(K) < 32/3 in ℝ⁴ would refute the conjecture. Such a
  claim requires, before it is written anywhere outside the attempt record:
  exact re-verification of vol(K) and vol(K°) by an **independent
  implementation** using a different triangulation and vertex ordering; exact
  verification that K = −K; and exact verification that the claimed vertex set
  is in convex position.
- Reproducing a known value (cube, cross-polytope, a published ℝ³ case) is
  `VERIFIED` with scope equal to the bodies actually recomputed — it says
  nothing about the conjecture.

## Harness (tier 0)

None yet. A contributor adding one should add it here.
