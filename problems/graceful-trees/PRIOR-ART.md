# Graceful Trees — prior art from this lab

> **Tier 1.** Reading this file makes an attempt `informed`.

Machine-readable index: `prior-art.json`. Full records: `attempts/`.

## Editorial view of the attack surface

- Probabilistic/statistical: for random trees, how many graceful labelings
  exist? Which tree features minimize the count? Minimizers point at where a
  counterexample would live — or are evidence none exists.
- Push class-specific verification (lobsters especially) to larger n than
  trees generally.
- SAT encoding is clean; use hardest-instance mining to find structurally
  "nearly ungraceful" trees.

Assessed as: counterexample extremely unlikely; the value is in the
labeling-count statistics and hardest-instance structure.

## Attempts

### 001 — Exact graceful-labeling counts, n ≤ 14 · `VERIFIED`

Exact essential (|Aut|-normalized) counts for all 5,444 trees on 4–14
vertices. Budget-capped before n = 15.

**Design note.** A specialized bitmask backtracker was chosen over #SAT:
counting requires enumerating every solution, where the backtracker wins by
orders of magnitude at this size. SAT remains the right tool for *existence*
on much larger single instances.

**Minimizers.** The global essential minimum is 1 for every n — always the
star. Runner-ups (essential 3–7, flat in n) are high-symmetry diameter-3/4
double brooms. So the minimum does **not** grow with n — but every minimizer
lies inside classes already *proven* graceful, so the flat floor is not
evidence of danger.

**Restricted minima do grow**, geometrically: non-caterpillar ~1.18^n,
non-lobster ~1.64^n. That is quantitative evidence against a counterexample
lurking near the frontier.

**Folklore corrected.** Maximizers are non-caterpillar lobsters, not paths.

**Novelty vs Anick's ≤16-edge database:** the normalized and class-restricted
analyses, not the raw counts.

## Open lines

- Mine the symmetric-spider seed (`LpH?GCAO??_@?A` genre) at n = 15–16.
- Lobster verification at larger n.
