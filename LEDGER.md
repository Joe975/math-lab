# Ledger

## TL;DR (updated 2026-07-25, cycle 3 = recovery cycle)

A second container restart (22:38Z) silently killed all three in-flight
agents (Erdős–Gyárfás n=22, Singmaster collisions, lonely-runner V=72);
all three resumed from checkpoints — no banked compute lost. No new lines
launched this cycle. Standing results: Erdős–Straus coverage map + Mordell
confirmation done (class-601 anomaly queued); union-closed idea B dead with
generalized no-go; lonely-runner validation recovered both Goddyn–Wong
tight instances at V≤40 (good tool signal). Next fresh launches: union-
closed idea C, Erdős–Straus 601-at-10^6, graceful-trees SAT.

## Problem status

| Problem | Status | Budget | Active line |
|---|---|---|---|
| Erdős–Gyárfás | active (run in flight) | high | girth≥5 exhaustion n=22 (→24); spectrum census next |
| Union-closed (Frankl) | active | high | ideas A/C/D remain (probability-charging constraints, dependent couplings); idea B closed |
| Erdős–Straus | active | medium | mine low-f(p) prime structure; explain class-601 anomaly |
| Singmaster | active | medium | multiplicity collision search (agent running) |
| Lonely runner | active | medium | k=8 near-tight scan (agent running) |
| Graceful trees | queued | low | SAT encoding; labeling-count statistics n ≤ 18 |
| Collatz | queued | low (long shot) | failed-approach taxonomy; cycle-bound frontier |

## Attempt queue (next cycles pull from the top)

1. [erdos-gyarfas] If attempt 002 record still missing: check/resume the n=22 search agent (slice checkpoints in attempts/erdos-gyarfas/data/ are reusable — never re-run completed slices).
2. [union-closed] Idea C (dependent/family-adaptive couplings): formalize the smallest nontrivial coupling class beyond Liu's conditionally-iid rung and test it against Sawin's family with the uc_weighted_kl.py machinery (adversarial test FIRST, per 002's protocol).
3. [erdos-straus] Class-601 anomaly: compute f(p) for primes to 10^6 in classes {1,49,73,97 mod 120} ∪ {601 mod 840}; is 601-enrichment real at scale or small-sample noise? (es_coverage.py has the C kernel.)
4. [graceful-trees] SAT encoding for graceful labeling; count labelings for all trees n ≤ 18.
5. [erdos-gyarfas] Cycle-spectrum realizability census from the n≤20 data (which length-sets occur?) — standalone interest.
6. [collatz] Failed-approach taxonomy page (library showcase; pure writing + citation verification).

## Verified results

- **[erdos-gyarfas] Cubic baseline n ≤ 20** (2026-07-25, attempt 001):
  every connected cubic graph on 4–20 vertices contains a cycle of length
  4, 8, or 16. 556,471 graphs; counts match OEIS A002851 at every n;
  spectra cross-checked against an independent implementation (networkx
  simple_cycles) on n=10,12 and generation cross-checked vs brute force
  (n ≤ 8). Computational evidence, cubic case only.
- **[erdos-gyarfas] Bridge bound** (2026-07-25, attempt 001): no near-cubic
  bridge-side block on ≤ 17 vertices avoids both C4 and C8 ⇒ a bridged
  cubic counterexample needs ≥ 38 vertices.
- **[union-closed] n ≤ 4 extremal check** (2026-07-25, attempt 001):
  exhaustive — minimum max element frequency over all union-closed families
  on ground sets of size ≤ 4 is exactly 1/2 (2/12/120/4958 families;
  optima triple-checked for closure).
- **[erdos-straus] Identity coverage + Mordell confirmation** (2026-07-25,
  attempt 001): 12 machine-verified polynomial identity families cover
  834/840 classes; the uncovered set mod 840 is exactly the 6 coprime
  quadratic residues {1,121,169,289,361,529}. All 9592 primes < 10^5
  solvable (0 failures), exact f(p) computed for each (C kernel,
  cross-checked vs independent Python counter for p < 3000).
- **[erdos-straus] Low-representation structure** (2026-07-25, attempt 001):
  the 24 lowest-f primes > 1000 are all ≡ 1 mod 24; bottom 50 concentrate
  in {1,49,73,97} mod 120; QR-mod-840 classes ~15× enriched but not
  characterizing (non-residue class 601 mod 840 holds 6 of bottom 50 —
  unexplained, queued).

## Insights / cross-problem notes

- Union-closed: the entropy method's ONLY use of closure is H(A∪B) ≤ log|F|
  for iid uniform A,B — an average-case fact, tight at (3−√5)/2 by
  Chase–Lovett's approximate family. Any advance must use worst-case
  closure of atypical/overlapping pairs, dependent couplings, or counting
  structure. See attempts/union-closed/001 §candidate-ideas (A–D).
- Erdős–Gyárfás: "{8}-only" near-misses share a bridged two-block anatomy;
  bridgeless assumption is now safe below 38 vertices. Heawood and
  Möbius–Kantor graphs are the canonical high-girth near-misses.
- Constrained graph generation (nauty geng, installed via apt) + SAT tooling
  is shared infrastructure for Erdős–Gyárfás and graceful trees.
- Ops: parallel subagents can die to 529 Overloaded during API load spikes;
  resume via SendMessage, and don't record a queue item as done until its
  files exist on disk.
- Ops (hardened after two container restarts): restarts kill agents, their
  compute, AND their waiters silently — an agent idling on a monitor dies
  without any notification. Orchestrator protocol every cycle: health-check
  in-flight agents via data-file mtimes + ps before launching new work;
  resume dead agents with a recovery message. Agent protocol: checkpoint
  every completed work unit to the repo immediately (per-slice files),
  never re-run completed units, and treat "the container can restart at any
  moment" as the design assumption.

## Dead ends

- **[union-closed] k-wise unions**: strictly worsen the entropy constant
  (0.382 → 0.318 → 0.276 for k = 2,3,4); recorded in attempt 001.
- **[union-closed] Gilmer's Conjecture 1** (strengthened entropy inequality):
  refuted by Sawin's construction — do not re-attempt as stated.
- **[union-closed] Weighted-KL ladder (idea B)** (2026-07-25, attempt 002):
  fully dead. Family-level version is vacuous for c ≤ 1, false for c > 1;
  distributional version killed for EVERY c ≥ 0 by Sawin's geometric-mixture
  family (c*(μ_n) → ∞ at marginals → ψ), with an exact finite-n certificate
  below the 0.38271 record. Root obstruction: KL charges escaping union-mass
  by log-likelihood (log(1/δ) for planted mass δ) vs Θ(n) entropy drop.
  Generalized no-go covers all smoothing-insensitive functionals
  Φ(law(U), μ) — see 002 before attempting ANY entropy-side strengthening.
