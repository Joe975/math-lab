# Collatz — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`.

## Attempts

**None.** This problem was queued but never worked. There is no prior art to
be informed by, so `blind` and `informed` mode are currently equivalent here.

## Editorial view of the attack surface

Deliberate long shot with minority budget. The intended value is the
*taxonomy*: the space of failed Collatz approaches is large and well
documented, so recording why each angle fails — stopping-time densities, 2-adic
reformulations, transfer operators, tag systems, cycle bounds via continued
fractions of log₂3 — is itself the deliverable.

Concrete lines, if you want them:

- Nontrivial cycles: sharpen numeric cycle-exclusion bounds (Simons–de Weger
  style) using current computational reach; document the exact frontier.
- Statistics of records/stopping times against the branching random walk
  model — quantify how well the stochastic model predicts extremes.
- Maintain a taxonomy of known failed approaches with the precise obstruction
  for each.

A warning worth repeating: the undecidability results for generalized Collatz
maps mean any approach not using something specific to the 3n+1 map is
already known to fail.
