# Optimal three-phase conducting composites in 2D

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** Mix three isotropic conductors with conductivities
σ₁ < σ₂ < σ₃ in prescribed volume fractions (f₁, f₂, f₃) into a fine-scale
periodic microstructure in the plane. The homogenized medium has an effective
conductivity tensor σ*; call the composite isotropic when σ* is a multiple of
the identity. The open problem is the **exact attainable range of the
isotropic effective conductivity** — equivalently, sharp lower and upper
bounds attained by explicit microstructures for every (f, σ) — and, more
broadly, the G-closure of three isotropic phases at fixed volume fractions.
For **two** phases this is completely solved; for three phases it is open in
an intermediate parameter regime, which is the frontier here.

## Published status

- **Two-phase case: closed.** The Hashin–Shtrikman (HS) bounds (J. Appl.
  Phys. 33, 1962, with the 2D form due to Hashin) are attained for all volume
  fractions — by coated-disk assemblages and equally by second-rank
  laminates; the full G-closure of two isotropic phases is known
  (Lurie–Cherkaev; Tartar–Murat).
- **Three-phase HS bounds are not always attainable.** Milton (1981) showed
  attainability of the lower bound requires the smallest-phase fraction to be
  large enough; Nesi (Proc. Roy. Soc. Edinburgh 125A, 1995) proved bounds
  strictly tighter than HS in 2D when f₁ is small. So for three phases the
  optimal bound is genuinely different from HS in part of parameter space.
- **Improved bounds and structures.** Cherkaev (2009) derived bounds tighter
  than Nesi's for three-phase 2D conductivity by the translation method
  augmented with pointwise field constraints, and exhibited multi-rank
  laminate families attaining them in complementary sub-regimes;
  Albin–Cherkaev–Nesi (2007) constructed multiphase laminates of extremal
  conductivity; Cherkaev–Zhang (2011) and Cherkaev–Pruss produced optimal
  "wheel" assemblages for some three-material regimes. The current best
  bounds and best structures meet on part of the (f, σ) space and leave a gap
  elsewhere — the attainability frontier is not characterized.
- **Known exact identities.** Keller–Dykhne duality: in 2D, homogenization
  commutes with pointwise inversion composed with a 90° rotation, giving
  exact relations between the effective tensors of a structure and its dual;
  for two-phase self-dual geometries at f = 1/2 it pins σ* = √(σ₁σ₂)
  (checkerboard). Duality maps the three-phase lower-bound problem to the
  upper-bound problem with inverted conductivities.
- **Structure class of record.** Hierarchical (multi-rank) laminates:
  iterated layerings at separated length scales. Their effective tensors are
  exact rational functions of fractions, directions and phase conductivities,
  and in many optimal-design problems laminates are known to suffice; whether
  they suffice for the three-phase isotropic problem everywhere is itself
  unsettled.
- **Physical background.** This is the mathematics of designed composite
  media — metamaterial engineering of transport: which effective response
  can geometry alone produce from given constituents? The same G-closure
  machinery underlies bounds on polycrystals, multiphase elasticity, and the
  design problems of transformation optics.

## Verification contract

Any claim recorded against this problem must meet the bar in `CONTRIBUTING.md`.
This problem is an exact-arithmetic fit, and the contract exploits that:

- **Effective tensors are exact or they are nothing.** A laminate's tensor is
  computed over ℚ (rational fractions, rational conductivities, integer
  direction vectors), and any claimed value ships the tree that produces it.
  Floating-point optimization may search, never certify; every use is
  labelled a screen.
- **Attainment of a bound is an exact equality in ℚ**, or a certified
  rational enclosure with the residual gap stated. "Numerically on the
  bound" is a candidate, not a claim.
- **Bound values state their formula and source.** Wiener and HS bounds are
  classical and implemented in the harness; tighter bounds (Nesi, Cherkaev)
  must be transcribed from the papers with the transcription marked per
  CONTRIBUTING rule 6, and re-derived or cross-checked before anything is
  killed against them.
- **Census claims are scoped by the structure class** — lamination rank,
  direction set, parameter grid — all of which appear in the record. A
  bounded-rank search that stalls below a bound is `EVIDENCE` about that
  class and grid, never about all microstructures; only a proof separates a
  bound from attainability.
- **Every laminate record must pass the Keller–Dykhne identity** (the dual
  tree with inverted phases must produce the rotated inverse tensor,
  exactly), and re-verification by the independent implementation
  (`verify_laminate.py`) is required before any ledger mention.

## Harness (tier 0)

- `harness/three-phase-conductivity/laminate.py` — reference implementation.
  Exact effective tensors of hierarchical laminates over ℚ via the
  frame-free projection form of the lamination formula (scale-invariant in
  the integer direction vector, so no normalization square roots ever
  appear); exact per-phase volume fractions; Wiener and 2D multiphase HS
  bounds in comparison-medium form; an exact Keller–Dykhne duality check on
  every tree.
- `harness/three-phase-conductivity/verify_laminate.py` — independent
  re-computation. Each lamination is redone by exact rotation into the layer
  frame (rational because conjugation by the scaled rotation [[u, −v],[v, u]]
  divides out the norm) followed by the classical interface-condition
  averaging, and every field of a claimed record — tensor, fractions,
  bounds, duality — is recomputed from scratch and compared. Disagreement
  between the two routes is reported as requiring escalation rather than
  resolved by picking one.
