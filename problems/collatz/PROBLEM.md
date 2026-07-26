# Collatz Conjecture

> **Tier 0.** Published background only. Nothing below reflects what this lab
> has tried. See `AGENTS.md`.

**Statement.** Iterating n → n/2 (n even), n → 3n+1 (n odd) from any positive
integer eventually reaches 1.

## Published status

Open. Verified to ~2^68. Tao (2019): almost all orbits attain almost bounded
values (in logarithmic density). Known undecidability results for *generalized*
Collatz maps warn that fully general methods must fail — any approach has to
use something specific to the 3n+1 map.

Simons–de Weger established numeric exclusion bounds for nontrivial cycles.

## Verification contract

- This problem carries a deliberately small budget in this lab, and the bar
  for a claimed advance is correspondingly high: an approach must state up
  front which known obstruction it evades and why.
- Heuristic/statistical arguments are `EVIDENCE` at best, and the branching
  random walk model is a model, not a proof device.
- Any claimed cycle bound must be reproducible and stated against the
  Simons–de Weger frontier.

## Harness (tier 0)

None yet. A contributor adding one should add it here.
