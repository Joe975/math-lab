# Collatz Conjecture (long shot)

**Statement.** Iterating n → n/2 (n even), n → 3n+1 (n odd) from any positive
integer eventually reaches 1.

**Status.** Open. Verified to ~2^68. Tao (2019): almost all orbits attain
almost bounded values (logarithmic density). Known undecidability results for
generalized Collatz maps warn that fully general methods must fail.

**Role in this lab.** Deliberate long shot with minority budget. Useful as a
stress test for the approach library: the space of *failed* Collatz
approaches is large and well-documented, so recording why each angle fails
(stopping-time densities, 2-adic reformulations, transfer operators, tag
systems, cycle bounds via continued fractions of log2(3)) is itself the
deliverable.

**Attack surface (modest, concrete).**
- Nontrivial cycles: sharpen numeric cycle exclusion bounds (Simons–de Weger
  style) using current computational reach; document the exact frontier.
- Statistics of records/stopping times vs. the branching random walk model —
  quantify how well the stochastic model predicts extremes.
- Maintain a taxonomy of known failed approaches with the precise obstruction
  for each (this is the library's showcase page).
