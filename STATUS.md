# Status

The live ledger: where every problem stands, what is queued next, and what has
already been ruled out. Read this before starting work.

## TL;DR (updated 2026-07-30)

**Open for contributions. No automated loop is currently running** — the hourly
cycle that produced these records was stopped once its in-flight work closed
out, and no cron trigger is active. Nothing is half-finished: every line is
written up, and the queue below is a list of starting points rather than
abandoned work.

To pick something up, take an item from the attempt queue and follow
`docs/CYCLE.md` — either by hand, one line at a time, or by restarting an
hourly routine with that file as the prompt. Everything needed is in this
directory. New attempts are welcome as pull requests; see `CONTRIBUTING.md`,
and `AGENTS.md` for whether to work blind or informed.

Current standing. HEADLINE: the union-closed dependent-couplings route
SURVIVED independent skeptic review — first live interface past the
entropy barrier in this lab (corrected ceiling 0.4315). Its gap (a) —
pointwise Plackett odds-ratio control — is now REFUTED as stated, in
both directions, with a sharp replacement range (005, confirmed by
006); the live restatement is *averaged* odds-ratio control, and the
perturbative lead survives untouched. The 2026-07-30 cycle also ran the
first two blind attempts of the lab: a mahler-4d census (minimum volume
product over all 1.77M centrally symmetric {0,±1}⁴-vertex bodies ≤ 22
vertices is EXACTLY 32/3, attained only by Hanner polytopes, spectral
gap ≈ 0.146 — no counterexample in that universe) and a billiards word
census (coverage staircase stalls at 135°; needed word length grows
hyperbolically, not super-exponentially, but certified area collapses —
length buys angle reach, not area; the ~112.5° constructive frontier
was rediscovered blind). Both skeptic-confirmed by from-scratch
re-implementations. Also complete: Singmaster census to 2.5×10^29;
Erdős–Gyárfás verified for ALL cubic graphs to n=22 plus all girth-≥5
cubics at n=24; lonely-runner k=8 census; Erdős–Straus 601 refuted →
QR-class identity-poverty; graceful census n ≤ 14. 15 verified
results, 5 recorded dead ends, 10 problems, ~35 reusable tools.
Crouzeix remains the one onboarded problem with no attempts (queue;
run blind).

## Problem status

| Problem | Status | Budget | Active line |
|---|---|---|---|
| Erdős–Gyárfás | n=24 done | high | next: 2-connected C16-free test (lead 2); n=26 is a ~16h run |
| Union-closed (Frankl) | ROUTE LIVE (verified) | high | gap (a) refuted as stated (005/006) → close the restated gaps: (a′) averaged odds-ratio control, (b) mutual-information tax, (c) perturbative assembly at ρ≈1.03 |
| Erdős–Straus | 601 resolved | medium | next: prove the identity-poor mechanism (why QR classes force low f) |
| Singmaster | census done | medium | next: Diophantine curve table (search-deeper is now low value) |
| Lonely runner | k=8 done | medium | next: k=9 scan (k=8 likely settled by Rosenfeld preprint) |
| Graceful trees | census done (n≤14) | low | possible next: mine the symmetric-spider seed; lobster verification at larger n |
| Collatz | queued | low (long shot) | failed-approach taxonomy; cycle-bound frontier |
| Triangular billiards | census done blind (skeptic-confirmed) | high | next: prove the W(a,b) death-angle laws (test W(10,9)→162.90°); pinch-gap scan ℓ=28–34; Taylor-model enclosure as a new tier-0 tool |
| Mahler in ℝ⁴ | census done blind (skeptic-confirmed) | medium | next: close k=12–20 (falsifiable: no proper mask with P<11); run the same pipeline on {0,±1}³ for the n=3 spectrum comparison |
| Crouzeix | onboarded, no attempts | medium | harness ready, certification risk retired; first attempt is the dim-3 landscape (run blind) |

## Attempt queue (next cycles pull from the top)

1. [union-closed] Proof-gap attack, one gap per agent: (a′) the AVERAGED odds-ratio control M_i ≥ λ from 005's restatement (the pointwise version is REFUTED — see 005/006; the averaged form survived its first falsifiable test on the crash family, but must first resolve the zero-mass bookkeeping for i ≥ 3 on small supports); (b) the mutual-information tax bound; (c) the perturbative assembly at ρ≈1.03 around the proven c=0 argument — 005's Prop 6 (first-order cancellation, λ + O(δ²) near product) is the natural starting point, noting its family-dependent asymptotic onset (006 part S7). Use 004+005+006's corrected statements ONLY. Any claimed bound → skeptic review before ledger entry.
2. [erdos-straus] Prove the identity-poverty mechanism: why does QR-class membership mod 840 force fewer Type I covering congruences? Start from 001's obstruction analysis + 002's rate data; target a theorem "f(p) ≥ g(N_typeI(p))" or a disproof.
3. [lonely-runner] k=9 near-tight scan: reuse lonely_runner.py (threshold near 1/9, feasibility analysis first — 8-tuples grow fast; consider restricting to accelerations/near-APs of known structures plus a bounded full scan).
4. [singmaster] Diophantine curve table: which equations C(n,j)=C(m,k) (small j<k) are resolved vs open, per the census lead that all in-range coincidences come from known families.
5. [erdos-gyarfas] Cycle-spectrum realizability census from the n≤20/22/24 data (which length-sets occur?) — standalone interest.
6. [collatz] Failed-approach taxonomy page (library showcase; pure writing + citation verification).
7. [graceful-trees] Mine the symmetric-spider seed (LpH?GCAO??_@?A genre) at n = 15-16 targeted; lobster verification at larger n.
8. [crouzeix] **Run blind.** Local-maxima census, dimension 3, polynomial degree ≤ 3: do all basins terminate at known extremal structure? `EVIDENCE` scoped by dimension, degree and search design, all of which must be in the record — and state how the search design differs from Greenbaum–Overton's before running it, since reproducing their finding is not a result. Kill condition (measured, not hypothetical): one certified ratio enclosure costs about 0.4s in dim 2, 0.75s in dim 3 and 1.4s in dim 4 at tolerance 1e-9 with 32 directions, so a wide census is affordable only in the low thousands of points; if the design needs more, narrow to structured families where the norm has a closed form and record the narrowing.
9. [billiards-triangles] Prove the W(a,b) death-angle laws from 001 (SPECULATION there; survived two out-of-sample rows in 002): the corridor endpoints are rational functions of the apex, so degeneration at γ = 180a/(a+1)° is a finite symbolic computation per a. Test predictions first: W(10,9) dies at 162.90°, W(12,12) at 166.154°. A proof turns the hyperbolic growth law from measured to established along the family curve.
10. [mahler-4d] Close the {0,±1}⁴ universe: k = 12–20 pairs (~30M orbits at k=12, improper fraction already 77% at k=9). Falsifiable: no proper mask with k ≥ 12 has P < 11. Needs the improper-detection shortcut or a streaming canonicalizer; see 001 lead 1. Cheap side quest, same pipeline: the {0,±1}³ census for the n=3 spectrum comparison (13 pairs, trivial) — does the non-Hanner gap grow or shrink with n?
11. [billiards-triangles] Pinch-gap scan: what covers [death(W(a,a)), birth(W(a+1,a))] (e.g. the measured [135.000°, 135.049°] where nothing ≤ 30 is alive)? Scan words of length 28–34 restricted to the family's letter statistics. If nothing bounded covers a neighbourhood of γ = 180a/(a+1)° on the arc-minimum curve, those angles are genuine accumulation points of the constructive problem.
12. [billiards-triangles] Coverage self-test: re-derive the acute and right-triangle cases as a scoped attempt record. Low value now that the harness self-test covers Fagnano and the orthic geometry and 001 mapped the obtuse side — take it only if something turns up that the certificate machinery cannot express.

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
- **[union-closed] Pointwise odds-ratio control REFUTED, both directions
  (skeptic-confirmed)** (2026-07-30, attempts 005+006): the tilt family's
  claimed pointwise control OR ≤ 2^λ fails upward on every diagonal
  history (OR ≥ 2^λ by Cauchy–Schwarz in the tilt kernel's inner product)
  and downward on an explicit 4-atom family with a positive-mass history
  where OR = 2^{λ(3−n)} → 0, with all marginals < 0.38271. Sharp universal
  range OR ∈ [2^{λ(1−m)}, 2^{λ(1+m)}], both ends attained. Near product
  measures log₂OR = λ + O(δ²) — the perturbative lead survives. 006
  re-derived every proposition by hand and re-computed with a from-scratch
  engine (sparse supports, off-grid λ); the Karlin–Rinott citation was
  checked against the paper and the binary-cube case re-proved via
  Ahlswede–Daykin. Live restatement: averaged odds-ratio control, which
  survived its first falsifiable test. Route stays LIVE.
- **[mahler-4d] {0,±1}⁴ census: minimum is exactly 32/3, Hanner-only
  (blind, skeptic-confirmed)** (2026-07-30, attempts 001+002): over all
  centrally symmetric conv(±S), S ≤ 11 antipodal pairs of nonzero {0,±1}⁴
  points — every centrally symmetric 4-polytope with ≤ 22 vertices in that
  lattice cube; 18,637,214 B4-orbits, Burnside-verified; 1,773,715 distinct
  bodies — the minimum volume product is exactly 32/3, attained by 1113
  orbits, each GL-certified Hanner; nearest non-attainer is a single
  10-vertex simplicial orbit at 32/3 + 4535/31104 ≈ +0.146. Skeptic
  re-derived the orbit counts, exactly recomputed all 1176 near-bound
  bodies from scratch, verified all 1113 certificates, and audited floats
  vs exact (max deviation 8.9e-15). No Mahler counterexample in this
  universe; k = 12–20 not covered.
- **[billiards-triangles] Obtuse word census to length 26 (blind,
  skeptic-confirmed)** (2026-07-30, attempts 001+002): 17,527 canonical
  translation words, 376 exact dyadic-box certificates (family orbits to
  length 66, γ ≤ 160°). Coverage staircase: ℓ=14 to 112.5°, 18 to 120°,
  22 to 130°, 26 to 135°, nothing ≤ 26 past 135°; the length-14 family
  dies at 112.4989° — a blind rediscovery of the published constructive
  frontier. Needed word length grows hyperbolically in the angle
  (ℓ ≈ 1440/(180−γ), death-law SPECULATION survived two out-of-sample
  prediction rows), but certified box width collapses exponentially in ℓ:
  length buys angle reach, not area. Negatives are sample-bounded; says
  nothing about unstable orbits.
- **[lonely-runner] k=8 tightness census to V=72** (2026-07-25, attempt
  001): among all 1,473,109,704 speed 7-tuples with max speed ≤ 72,
  exactly 3 primitives have ML < 13/100, all with ML = 1/8 exactly — (1..7)
  plus the two Goddyn–Wong instances (recovered from scratch). Nothing
  below 1/8. The ML spectrum below 1/7 (V ≤ 40) is exactly {s/(7s+k),
  k∈{1,2}} — consistent with Fan–Sun's amended spectrum conjecture, new
  data point at n=7. Exact rational arithmetic, independently re-verified.

## Insights / cross-problem notes

- Infrastructure (2026-07-27): cross-pollination layer added. `mechanisms.json`
  maps every approach tag to a field lens (`scripts/mechanisms.py` for
  gaps/matrix queries); `docs/IDEATE.md` = field-sweep ideation on one
  problem, `docs/RIPPLE.md` = propagate a new result across problems; both
  are hooked into `docs/CYCLE.md` and have session skills (`/ideate`,
  `/ripple`, `/mechanisms`, `/cycle`). No sweep or scan has been run yet.
- Onboarding (2026-07-28): three problems added, and two literature findings
  from the check that changed their framing. Triangular billiards — the
  constructive frontier is 112.3° (Garber–Marinov–Moore–Tokarsky 2018), not
  the 100° of Schwartz's theorem; and Forni's June 2026 preprint
  (arXiv:2606.10102) claims a periodic orbit in *every* polygon, which if it
  holds settles existence. It is non-constructive, so the certificate-covering
  frontier is untouched — but frame any attempt against the constructive
  question, not against existence. Mahler — the Viterbo counterexample that
  killed its physics motivation is non-symmetric, so the symmetric variant
  that implies Mahler still stands; `/ripple` on that refutation is a cheap
  early exercise. Crouzeix's certification risk is retired: stdlib-only
  certified enclosures are affordable (timings in queue item 10).
- Interval arithmetic is the wrong default for iterated geometry
  (2026-07-28, from the billiards harness). Unfolding recombines coordinates
  that earlier steps already widened, and an interval cannot see that those
  errors are the same error, so the enclosure grew about an order of magnitude
  per reflection — a six-step word was hopeless from a box of width 1e-6.
  First-order affine forms made the growth additive and the certificates
  possible. Two smaller traps in the same file: squaring must be its own
  operation (a generic product of x with itself loses non-negativity and a
  length comes out zero), and a quantity that is a geometric invariant must be
  computed once from the original data rather than recomputed from widened
  coordinates. Any future harness that iterates a map over a parameter box
  should start from affine forms rather than discover this again.
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
- Blind mode produced its first data points (2026-07-30): both blind
  attempts independently rediscovered published structure — the billiards
  census hit the ~112.5° constructive frontier to four decimal places, and
  the mahler census recovered the full Hanner equality case — without
  access to PRIOR-ART.md or the framing there. Early but real evidence
  that the harnesses alone are enough to orient an attempt, and a baseline
  for the anchor-vs-help question the mode field exists to answer.
- Ops (2026-07-30): background watcher loops can die silently even without
  a container restart — an agent that parks itself on "wake me when the
  slices finish" may never be woken although its compute completed fine.
  Orchestrator protocol unchanged but sharpened: on any agent-stopped
  notification, health-check the compute (data-file mtimes + ps) and, if
  the work is done but the agent idle, resume it with a state summary
  rather than waiting. Also: 529-Overloaded kills can hit the same agent
  repeatedly; resume works, but make restructuring/formatting jobs
  orchestrator-side rather than burning agent restarts on them.
- Ops (hardened after two container restarts): restarts kill agents, their
  compute, AND their waiters silently — an agent idling on a monitor dies
  without any notification. Orchestrator protocol every cycle: health-check
  in-flight agents via data-file mtimes + ps before launching new work;
  resume dead agents with a recovery message. Agent protocol: checkpoint
  every completed work unit to the repo immediately (per-slice files),
  never re-run completed units, and treat "the container can restart at any
  moment" as the design assumption.

## Dead ends

- **[union-closed] Pointwise Plackett odds-ratio control (gap (a) as
  stated)** (2026-07-30, attempts 005+006): dead in both directions — do
  not re-attempt any pointwise uniform-in-μ version; diagonal histories
  force the opposite inequality and structured 4-atom adversaries crash
  the ratio to 0 within the record-relevant marginal regime. The averaged
  restatement (a′) is the live replacement; see queue item 1.
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
