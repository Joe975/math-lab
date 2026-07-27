# Hadwiger–Nelson — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None yet.** Added 2026-07-27; the first attempts are in flight. Until one
lands there is no prior art to be informed by, so `blind` and `informed` mode
are equivalent here.

## Editorial view of the attack surface

Why this problem was added: it is the rare famous conjecture with a frontier
that moves in *small verifiable increments*. Polymath16 spent three years
shrinking a graph from 1581 vertices to 509 — every step a self-contained,
independently checkable object. That is a shape this lab can actually
contribute to, unlike a problem whose only currency is a full proof.

The exact-arithmetic requirement in `PROBLEM.md` is the same discipline that
made the lonely-runner census trustworthy: the entire content of a
unit-distance graph is that certain distances are *exactly* 1, so a float turns
a theorem into a rounding artefact.

Concrete lines, if you want them:

- Build the exact-field arithmetic and unit-distance certifier, and validate it
  by recovering the classical small graphs (Moser spindle, Golomb graph) from
  scratch — the calibration step this lab does before claiming anything.
- Independent re-verification of a published 5-chromatic graph. Nobody here has
  checked one, and re-deriving non-4-colourability with an independent
  implementation would be a genuine `VERIFIED` contribution about that graph.
- Map the spindling mechanism: how chromatic number grows against vertex count
  under iterated Moser spindling, and where that construction stops paying.
- The 6-chromatic question is the real frontier and nothing here should be
  expected to touch it. An attempt claiming a 6-chromatic unit-distance graph
  has almost certainly made an arithmetic error — see the calibration note in
  `GUIDANCE.md`.
