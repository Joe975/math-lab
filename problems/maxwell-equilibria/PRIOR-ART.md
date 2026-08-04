# Maxwell equilibria — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None.** This problem was onboarded and has not been worked. There is no
prior art to be informed by, so `blind` and `informed` mode are currently
equivalent here — which makes the first attempts worth running blind, since
blind costs nothing while the record is empty.

## Editorial view of the attack surface

This problem was onboarded days after the general conjecture was claimed
false (arXiv:2607.27197, July 29 2026), which is exactly why the budget is
high rather than why it should be low. What died is the (n−1)² count at
n = 5; what the refutation *opened* is a set of questions shaped like this
lab's tooling:

- **The refutation has no explicit witness.** The perturbation argument
  holds "for all sufficiently small ε" and certifies no concrete
  configuration. Producing a fully rational instance and certifying ≥ 24
  isolated equilibria in exact arithmetic would be the first independent
  verification of the counterexample — and would either confirm it or find a
  gap in a days-old result. Route note: the equilateral triangle embeds
  rationally in ℝ³ as e₁, e₂, e₃ (pairwise distances √2, plane x+y+z = 1,
  symmetry axis along (1,1,1)), so the whole five-charge configuration can
  be taken with rational coordinates and rational charges — the axial pair
  at (1/3, 1/3, 1/3) ± t·(1,1,1) with rational t, and a rational charge
  near the scaled (3/4)ε³ law. Everything the harness needs is then in ℚ.
  Expect the 21 bifurcated equilibria to be nearly degenerate for small t
  (they merge as t → 0): the certification cost grows as t shrinks and the
  equilibria may vanish for t large, so the work is finding the Goldilocks
  window. Kill condition: if no rational t admits certification at
  affordable subdivision depth, record the cost curve and the tightest
  bracketing achieved — that is a finding about the window, and `EVIDENCE`
  that the asymptotic regime is narrower than the preprint's framing
  suggests.
- **n = 3 is open between 4 and 6.** The community's intuition about this
  problem was just proven wrong at n = 5; the conjectured max 4 at n = 3
  deserves suspicion too. A census over three positive charges — two shape
  parameters and two charge ratios after similarity reduction — hunting for
  a configuration with ≥ 5 certified equilibria is self-policing in the
  same way the Mahler volume-product search is: a single certified witness
  settles a question people have worked on since 1873. The equal-magnitude
  case is closed (Tsai 2015, max 4), so weight the census toward unequal
  charges; near-degenerate strata (where equilibria merge and the count
  jumps) are both where new equilibria appear and where certification is
  most expensive.
- **The index-sum identity (Σ indices = 1 − n) is a completeness alarm.**
  Any complete census that fails it has missed an equilibrium or
  miscounted an index. It is implemented in the harness as a mandatory
  self-check; treat a violation as a bug until proven otherwise.

Concrete lines, if you want them:

- Self-test first (blind): collinear positive charges (exactly n−1
  equilibria, classical) and unit charges at an equilateral triangle
  (exactly 4, Tsai). Purpose is validating the certification machinery
  against known answers — a mismatch is a bug in our enclosures, not a
  finding.
- Then the explicit five-charge witness (informed; it consumes the
  arXiv:2607.27197 construction as a labelled input).
- Then the three-charge census.

Methodological warning carried over from the rest of the lab: the equilibria
of symmetric configurations sit at symmetric points, and subdivision grids
love to put cell faces through symmetric points. The harness offsets its
root box by non-dyadic rationals so certified zeros never land on cell
faces; keep that property if you modify the search design, or isolation
certificates will fail forever at the symmetric equilibria and the run will
report spurious UNRESOLVED regions.
