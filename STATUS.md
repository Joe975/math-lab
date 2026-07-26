# Ledger

## TL;DR (updated 2026-07-26 — LAB PAUSED by user; hourly cycle stopped)

The user stopped the loop after the in-flight work completed. All lines
are closed out; nothing is running; the cron trigger is deleted. To
resume: recreate an hourly routine with the standard cycle prompt (see
README "How the loop works") — all state needed is in this directory.

Final standing. HEADLINE: the union-closed dependent-couplings route
SURVIVED independent skeptic review — first live interface past the
entropy barrier in this lab (corrected ceiling 0.4315; three labeled
proof gaps are queue item 1). Complete: Singmaster census to 2.5×10^29
(3003 unique at mult 8); Erdős–Gyárfás verified for ALL cubic graphs to
n=22 plus all girth-≥5 cubics at n=24 (zero candidates); lonely-runner
k=8 census (Goddyn–Wong recovered; k=9 frontier); Erdős–Straus 601
refuted → real signal is QR-class identity-poverty (p ≈ 10^-110);
graceful census n ≤ 14. 12 verified results, 4 recorded dead ends,
7 problems, ~20 reusable tools.

## Problem status

| Problem | Status | Budget | Active line |
|---|---|---|---|
| Erdős–Gyárfás | n=24 done | high | next (on resume): 2-connected C16-free test (lead 2); n=26 is a ~16h run |
| Union-closed (Frankl) | ROUTE LIVE (verified) | high | close the three proof gaps in the couplings route (skeptic-corrected, ceiling 0.4315) |
| Erdős–Straus | 601 resolved | medium | next: prove the identity-poor mechanism (why QR classes force low f) |
| Singmaster | census done | medium | next: Diophantine curve table (search-deeper is now low value) |
| Lonely runner | k=8 done | medium | next: k=9 scan (k=8 likely settled by Rosenfeld preprint) |
| Graceful trees | census done (n≤14) | low | possible next: mine the symmetric-spider seed; lobster verification at larger n |
| Collatz | queued | low (long shot) | failed-approach taxonomy; cycle-bound frontier |

## Attempt queue (next cycles pull from the top)

1. [union-closed] Proof-gap attack, one gap per agent: (a) Plackett odds-ratio control for the tilt family; (b) the mutual-information tax bound; (c) the perturbative assembly at ρ≈1.03 around the proven c=0 argument. Use 004's corrected statements ONLY (ceiling 0.4315, restated mini-theorem, half-mixing coupling for the {0}∪[½,1) genre). Any claimed bound → skeptic review before ledger entry.
2. [erdos-straus] Prove the identity-poverty mechanism: why does QR-class membership mod 840 force fewer Type I covering congruences? Start from 001's obstruction analysis + 002's rate data; target a theorem "f(p) ≥ g(N_typeI(p))" or a disproof.
3. [lonely-runner] k=9 near-tight scan: reuse lonely_runner.py (threshold near 1/9, feasibility analysis first — 8-tuples grow fast; consider restricting to accelerations/near-APs of known structures plus a bounded full scan).
4. [singmaster] Diophantine curve table: which equations C(n,j)=C(m,k) (small j<k) are resolved vs open, per the census lead that all in-range coincidences come from known families.
5. [erdos-gyarfas] Cycle-spectrum realizability census from the n≤20/22/24 data (which length-sets occur?) — standalone interest.
6. [collatz] Failed-approach taxonomy page (library showcase; pure writing + citation verification).
7. [graceful-trees] Mine the symmetric-spider seed (LpH?GCAO??_@?A genre) at n = 15-16 targeted; lobster verification at larger n.

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
- **[singmaster] Complete multiplicity census to 2.5×10^29** (2026-07-25,
  attempt 001): 3003 is the unique multiplicity-8 value; exactly seven
  mult-6 values (120, 210, 1540, 7140, 11628, 24310, C(104,39)=C(103,40));
  mult 5 and 7 provably empty in range. Self-test + independent brute
  force to 10^7 + double re-verification of every hit + agreement with
  de Weger's coincidence list.
- **[erdos-gyarfas] n=22 exhausted** (2026-07-25, attempt 002): all 90,938
  girth-≥5 connected cubic graphs on 22 vertices contain C8 (count matches
  OEIS A014372 exactly; C kernel re-validated bit-for-bit vs the Python
  spectrum tool on 4,569 graphs). Conjecture verified for all cubic graphs
  through 22 vertices with our own reproducible tooling.
- **[erdos-gyarfas] n=24 girth-≥5 exhausted** (2026-07-26, attempt 002
  addendum): all 1,620,479 girth-≥5 connected cubic graphs on 24 vertices
  contain C8 — zero candidates. Total matches A014372; slice totals,
  girth histogram, and combo counts independently reconcile; the unique
  girth-7 graph (McGee) appears exactly once as predicted. Eight
  {8}-only graphs recorded (vs 1 at n=22) — bridgedness check folded
  into lead 2.
- **[union-closed] Dependent-couplings route VERIFIED LIVE (skeptic-
  corrected)** (2026-07-26, attempts 003+004): exact union-closure
  licenses every coupling of uniform marginals (lemma re-derived
  independently); the refuted Gilmer Conjecture 1 is exactly the
  diag⊕iid rung (two-way equivalence confirmed); overlap-tilt couplings
  separate Sawin and Chase–Lovett killer families (independent
  implementation matches to 3 decimals); the functional genuinely evades
  the 002 no-go. Skeptic corrections: single-λ model ceiling is 0.431496
  in closed form (NOT 0.445 — grid-floor artifact), extremal genre
  δ_∅⊕Bern(½+ε); mini-theorem restated (fails on {0}∪[½,1) mixtures,
  that genre now covered by a verified half-mixing C₃ coupling). Route
  status: LIVE — first interface to survive the full protocol; ceiling
  0.4315 ≫ current record 0.38271; three labeled proof gaps remain
  (Plackett odds-ratio control, mutual-information tax, perturbative
  assembly at ρ≈1.03).
- **[erdos-straus] Low-f characterization at 10^6; 601 refuted**
  (2026-07-26, attempt 002): exact f(p) for all 9732 primes ≡ 1 mod 24 to
  10^6 (kernel triple-validated). Class-601 anomaly was small-sample
  noise (0 of size-normalized bottom 49 vs 2.1 expected). Real signal:
  size-normalized bottom 2% is 98% inside Mordell's six QR classes mod
  840 (share 24.4%, p ≈ 10^-110). Mechanism: low-f primes are
  identity-poor (class mean f monotone in number of covering Type I
  families, r = 0.92; smooth p−1 further depresses f). Band-minima of f
  grow like (log p)^3.
- **[graceful-trees] Labeling census n ≤ 14** (2026-07-26, attempt 001):
  exact essential (|Aut|-normalized) graceful counts for all 5,444 trees;
  min = 1 always (the star). All global minimizers lie in proven-graceful
  rigid classes; restricted minima grow geometrically (non-caterpillar
  ~1.18^n, non-lobster ~1.64^n) — quantitative evidence against a nearby
  counterexample. Maximizers are non-caterpillar lobsters, not paths
  (folklore corrected). Novelty vs Anick's ≤16-edge database: the
  normalized and class-restricted analyses.
- **[lonely-runner] k=8 tightness census to V=72** (2026-07-25, attempt
  001): among all 1,473,109,704 speed 7-tuples with max speed ≤ 72,
  exactly 3 primitives have ML < 13/100, all with ML = 1/8 exactly — (1..7)
  plus the two Goddyn–Wong instances (recovered from scratch). Nothing
  below 1/8. The ML spectrum below 1/7 (V ≤ 40) is exactly {s/(7s+k),
  k∈{1,2}} — consistent with Fan–Sun's amended spectrum conjecture, new
  data point at n=7. Exact rational arithmetic, independently re-verified.

## Insights / cross-problem notes

- Union-closed: the entropy method's ONLY use of closure is H(A∪B) ≤ log|F|
  for iid uniform A,B — an average-case fact, tight at (3−√5)/2 by
  Chase–Lovett's approximate family. Any advance must use worst-case
  closure of atypical/overlapping pairs, dependent couplings, or counting
  structure. See problems/union-closed/attempts/001 §candidate-ideas (A–D).
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
