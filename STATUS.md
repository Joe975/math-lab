# Status

The live ledger: where every problem stands, what is queued next, and what has
already been ruled out. Read this before starting work.

## TL;DR (updated 2026-07-31)

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

Current standing. HEADLINE: the 2026-07-31 cycle attacked all three
gaps of the LIVE union-closed route in parallel; every verdict is
skeptic-confirmed. Gap (a′) — averaged odds-ratio control — is
REFUTED by an explicit 10-atom witness (M_5 = λ − 0.122, all
marginals ≤ 0.318 < 0.38271), certified float-free in exact rational
arithmetic and robust to the bookkeeping convention; the i-AGGREGATED
control survives certified (+1.84 on the witness itself) and is the
restated gap. Gap (c): the smoothness step 005's Prop 6 was missing is
now a proved theorem (P6′, quantitative IFT, skeptic-verified at fixed
n) — but the skeptic's n=32 census killed the naive n-uniform budget
(the flattening reverses past n ≈ 22) and showed the tax is δ-linear,
so the conditional assembly needs a corrected budget object. Gap (b)
got its first precise statement and survived everywhere tested (the
slice tilt is provably tax-free; a newly named "second tax"
~0.40·log₂n is the binding loss channel; λ-window law λ_max ≈
4.847/(n−3)). Route stays LIVE, ceiling 0.4315 vs record 0.38271.
Also this cycle: the billiards W(a,b) death-angle laws held both
pre-registered out-of-sample predictions, were unified into a general
law γ_d(a,b) = 180 − 90(a+b)/(a(b+1)), and PROVEN as a necessity
theorem (machine-certified for 15 members, skeptic-confirmed); and the
first cross-field ideation sweep ran on union-closed (6 lenses → 3
queue-worthy routes, 5 no-purchase verdicts, one proven product-weight
no-go). Standing library: Singmaster census to 2.5×10^29;
Erdős–Gyárfás for all cubics to n=22 + girth-≥5 at n=24; lonely-runner
k=8; Erdős–Straus identity-poverty; graceful n ≤ 14; mahler-4d and
billiards blind censuses. 19 verified results, 7 recorded dead ends,
10 problems, ~40 reusable tools. Crouzeix remains the one problem with
no attempts (queue; run blind).

## Problem status

| Problem | Status | Budget | Active line |
|---|---|---|---|
| Erdős–Gyárfás | n=24 done | high | next: 2-connected C16-free test (lead 2); n=26 is a ~16h run |
| Union-closed (Frankl) | ROUTE LIVE (verified) | high | gap (a′) refuted (007/013) → live gaps: aggregated OR control (probe n ≳ 20 BEFORE proof effort), δ-linear tax budget, corrected assembly budgets (012); plus three sweep leads (010) |
| Erdős–Straus | 601 resolved | medium | next: prove the identity-poor mechanism (why QR classes force low f) |
| Singmaster | census done | medium | next: Diophantine curve table (search-deeper is now low value) |
| Lonely runner | k=8 done | medium | next: k=9 scan (k=8 likely settled by Rosenfeld preprint) |
| Graceful trees | census done (n≤14) | low | possible next: mine the symmetric-spider seed; lobster verification at larger n |
| Collatz | queued | low (long shot) | failed-approach taxonomy; cycle-bound frontier |
| Triangular billiards | death law PROVEN as necessity (003/004, skeptic-confirmed) | high | next: sufficiency at γ_d (first-order corner argument); general-(a,b) proofs of I1–I3 / Lemma C; a > 2b+3 branch; pinch-gap scan |
| Mahler in ℝ⁴ | census done blind (skeptic-confirmed) | medium | next: close k=12–20 (falsifiable: no proper mask with P<11); run the same pipeline on {0,±1}³ for the n=3 spectrum comparison |
| Crouzeix | onboarded, no attempts | medium | harness ready, certification risk retired; first attempt is the dim-3 landscape (run blind) |

## Attempt queue (next cycles pull from the top)

1. [union-closed] The restated gap: i-AGGREGATED odds-ratio control (007's replacement after the per-i averaged form died — see 007/013). FIRST extend the certified probe to larger n: 012's orbit-symmetrized engine reaches n = 32, and its budget-growth reversal past n ≈ 22 is a standing warning against small-n optimism — the aggregate is certified positive only at n ≤ 7. If it survives n ≳ 20, attack the proof through 007's Gram/Frobenius machinery (PSD mass matrix; first-order-in-λ perfect square) and the margin-modulated variant. Any claimed bound → skeptic review before ledger entry.
2. [union-closed] Rebuild the assembly budgets per 012's corrections: the tax is δ-LINEAR (restate B2), the τ_half step needs s₀ ≤ 0.0843, corrected δ₀ ≈ 0.004. The open question is whether ANY n-uniform budget object exists — 008's conditional theorem survives structurally; its constants need the corrected inputs, and the budget census machinery (uc_pert.py + uc_pert_skeptic.py's orbit engine) is ready.
3. [union-closed] Sweep leads from 010, one per cycle: (i) pairwise-closure LP/degree-2 SOS certification on the n ≤ 4 census (all 4958 families) — kill: certified bound converges to ≈ 0.382; (ii) union-transfer-operator eigenvalue field — test log-supermodularity of λ_C on the n ≤ 5 census (pre-check the total-positivity records first); (iii) bipartite-MIS decomposition — geng census to n = 12, connectivity of extremals, compositionality across 1- and 2-cuts. Kill conditions recorded in 010.
4. [erdos-straus] Prove the identity-poverty mechanism: why does QR-class membership mod 840 force fewer Type I covering congruences? Start from 001's obstruction analysis + 002's rate data; target a theorem "f(p) ≥ g(N_typeI(p))" or a disproof.
5. [lonely-runner] k=9 near-tight scan: reuse lonely_runner.py (threshold near 1/9, feasibility analysis first — 8-tuples grow fast; consider restricting to accelerations/near-APs of known structures plus a bounded full scan).
6. [singmaster] Diophantine curve table: which equations C(n,j)=C(m,k) (small j<k) are resolved vs open, per the census lead that all in-range coincidences come from known families.
7. [erdos-gyarfas] Cycle-spectrum realizability census from the n≤20/22/24 data (which length-sets occur?) — standalone interest.
8. [collatz] Failed-approach taxonomy page (library showcase; pure writing + citation verification).
9. [graceful-trees] Mine the symmetric-spider seed (LpH?GCAO??_@?A genre) at n = 15-16 targeted; lobster verification at larger n.
10. [crouzeix] **Run blind.** Local-maxima census, dimension 3, polynomial degree ≤ 3: do all basins terminate at known extremal structure? `EVIDENCE` scoped by dimension, degree and search design, all of which must be in the record — and state how the search design differs from Greenbaum–Overton's before running it, since reproducing their finding is not a result. Kill condition (measured, not hypothetical): one certified ratio enclosure costs about 0.4s in dim 2, 0.75s in dim 3 and 1.4s in dim 4 at tolerance 1e-9 with 32 directions, so a wide census is affordable only in the low thousands of points; if the design needs more, narrow to structured families where the norm has a closed form and record the narrowing.
11. [billiards-triangles] Close the death-law sufficiency: certified alive points approaching γ_d (the first-order corner argument, 003 lead 2 — note 004's correction that a zero-width touching corridor AT γ_d is not excluded by the necessity theorem); prove I1–I3 and Lemma C for general (a,b) (per-member certified for 15 members; Lemma C numerically true for all b ≤ a ≤ 25); map the a > 2b+3 branch (W(6,3)-type members are unmeasured).
12. [mahler-4d] Close the {0,±1}⁴ universe: k = 12–20 pairs (~30M orbits at k=12, improper fraction already 77% at k=9). Falsifiable: no proper mask with k ≥ 12 has P < 11. Needs the improper-detection shortcut or a streaming canonicalizer; see 001 lead 1. Cheap side quest, same pipeline: the {0,±1}³ census for the n=3 spectrum comparison (13 pairs, trivial) — does the non-Hanner gap grow or shrink with n?
13. [billiards-triangles] Pinch-gap scan: what covers [death(W(a,a)), birth(W(a+1,a))] (e.g. the measured [135.000°, 135.049°] where nothing ≤ 30 is alive)? Scan words of length 28–34 restricted to the family's letter statistics. If nothing bounded covers a neighbourhood of γ = 180a/(a+1)° on the arc-minimum curve, those angles are genuine accumulation points of the constructive problem. The general law γ_d(a,b) from 003 now predicts where W(a+1,a) births/deaths sit — use it to target the scan.
14. [billiards-triangles] Coverage self-test: re-derive the acute and right-triangle cases as a scoped attempt record. Low value now that the harness self-test covers Fagnano and the orthic geometry and 001 mapped the obtuse side — take it only if something turns up that the certificate machinery cannot express.

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
- **[union-closed] Gap 2 (mutual-information tax) formalized; probe
  survives (skeptic-confirmed)** (2026-07-31, attempts 009+011): first
  precise candidate statement (TAX at p) via a chain-rule assembly value
  CR and a newly named second tax ST ≥ 0; smoothing-sensitive on 002's
  own certificate gadgets (Θ(n) positive — so the 002 no-go does not
  apply); the pure slice tilt is provably tax-free and the half-mixing
  coupling provably has ST = 0 (both proofs independently re-derived);
  survives all four adversary genres, and the skeptic closed 009's
  untested tilt-recipe/Sawin cell positively to n = 300. λ-window law:
  λ_max ≈ 4.847/(n−3), sup CR at λ = 0 for n ≥ 14 (crash genre pulls
  λ → 0, mixtures pull λ ≳ 1 — gap 3 in miniature). Second-tax scaling
  ~0.40·log₂n survives extension to n = 240 (slope is λ-dependent).
  Finite-n EVIDENCE; four reporting-level corrections in 011.
- **[union-closed] Gap (a′) — per-i averaged OR control — REFUTED
  float-free (skeptic-confirmed)** (2026-07-31, attempts 007+013):
  well-posedness first forced (proved degeneracy dichotomy; the
  normalization is pinned by the i=n and product-μ anchors), then the
  control PROVED for four subclasses (≤4-atom supports — upgrading 006's
  S8 to a corollary — potential-MTP₂ μ, i ∈ {1,n}, first order in λ with
  a perfect-square coefficient), then REFUTED in the record-relevant
  regime: an explicit 10-atom μ on 2^[7] (marginals ≤ 0.318 < 0.38271)
  has M_5 = λ − 0.122033, violating for every λ ≳ 0.03. Certified in
  exact rational arithmetic (at rational tilts t = 4, 16 the violation
  is a fully rational statement), perturbation-stable, and robust to the
  conjecture-friendliest alternative bookkeeping (still −0.067). The
  i-AGGREGATED control is certified positive (+1.844669) on the witness
  and survived seeded attacks — it is the restated gap.
- **[union-closed] Theorem P6′: the Prop-6 smoothness step closed
  (skeptic-confirmed at fixed n); naive n-uniform budgets killed**
  (2026-07-31, attempts 008+012): existence/uniqueness/differentiability
  of the symmetric Sinkhorn potential in μ proved via quantitative IFT;
  |log₂OR − λ| ≤ 924·min(p,1−p)^(−3n)·δ², re-derived line-by-line and
  stress-tested to the theorem boundary with zero violations.
  Calibration corrected: ρ*(0.383) = 1.0422 (the queue's 1.03 was
  wrong); centered kernel suppresses the crash OR to 2^{λ(3−n)/n}.
  Skeptic kills on the assembly half: the averaged-budget flattening
  REVERSES past n ≈ 22 (008's own pre-asymptotic worry was real — B1
  dead as stated), the tax is δ-linear not quadratic (B2 restated), the
  box-uniform τ_half step is analytically false at its own parameters
  (s₀ ≤ 0.0843), δ₀ → ≈ 0.004. The conditional theorem survives
  structurally with corrected constants.
- **[billiards-triangles] W(a,b) death-angle law PROVEN as necessity
  and generalized (skeptic-confirmed)** (2026-07-31, attempts 003+004):
  both pre-registered out-of-sample predictions held (W(10,9) → 162.90°,
  W(12,12) → 2160/13°, agreement ~1e-12); a division-free Laurent-ring
  unfolding collapses the corridor criterion to a glide-axis offset
  inside gate projections, with the binding functions factoring exactly;
  the general law γ_d(a,b) = 180 − 90(a+b)/(a(b+1)) unifies 001's two
  slice formulas and predicted W(5,3) → 144° and W(4,2) → 135° BEFORE
  measurement (mirror-half members 001/002 never scanned). Necessity —
  no positive-width corridor at any γ ≥ γ_d — machine-certified for 15
  members and fully re-verified by an independent ring, an independent
  interval stack, and a 1.94M-point case-tree search. Scope: zero-width
  touching at γ_d is not excluded; the general-(a,b) identities remain
  SPECULATION (per-member certified); a > 2b+3 uncovered.

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

- Verification standard (2026-07-31): when a claim is decidable over Q,
  exact rational arithmetic IS the skeptic pass — cross-engine float
  agreement is NOT independence. 007's first "witnesses" were
  shared-IEEE-underflow artifacts that two independent engines AND a
  round-trip check all agreed on; a physical parameter sweep exposed
  them, and 013's Fraction-plus-certified-log₂-enclosures settled the
  real witness beyond floats. Billiards 003/004 ran the same standard
  geometrically (exact-rational intervals, independent stacks). New
  mechanism tag: exact-rational-arithmetic. Guard adopted in engines:
  cell-floor diagnostics + dynamic-range clamps.
- Small-n flattening is not evidence of n-uniformity (2026-07-31): 012
  extended 008's budget census from n = 9 to n = 32 with an
  orbit-symmetrized engine and the "flattening" reversed past n ≈ 22 —
  the pre-asymptotic worry 008 itself recorded was real. Corollary
  applied to the queue: the restated aggregated-control gap must be
  probed at n ≳ 20 before proof effort is invested (its certified
  positives live at n ≤ 7).
- First ideation sweep ran (2026-07-31, union-closed 010): 6 lenses,
  3 queue-worthy routes, 5 recorded no-purchase verdicts. The filter
  step earned its keep — one lens agent's numeric "witness" did not
  reconcile with its own construction and was replaced by a one-line
  proof of the same no-go. Sweeps on other problems remain unrun.
- Infra (2026-07-31): tests/test_site.py checks one_line summaries with
  html.escape() (quote=True) while scripts/build_site.py renders them
  with quote=False, so an apostrophe in the first ~40 characters of a
  one_line fails the site test. Worked around by wording; a proper fix
  should align the two escapings. Related: leak_terms must name
  findings, not generic vocabulary — "certified enclosure" collided
  with tier-0 problem statements and tripped the leak test.

## Dead ends

- **[union-closed] Pointwise Plackett odds-ratio control (gap (a) as
  stated)** (2026-07-30, attempts 005+006): dead in both directions — do
  not re-attempt any pointwise uniform-in-μ version; diagonal histories
  force the opposite inequality and structured 4-atom adversaries crash
  the ratio to 0 within the record-relevant marginal regime.
- **[union-closed] Per-i averaged odds-ratio control M_i ≥ λ (gap (a′),
  005's restatement)** (2026-07-31, attempts 007+013): dead in the
  record-relevant regime — a 10-atom witness certified in exact rational
  arithmetic violates at every λ ≳ 0.03 and survives every bookkeeping
  convention. Do not re-attempt any per-i averaged form; the witness
  genre (near-proportional light slices, anti-aligned cross-ratios) is
  the reusable adversary. Live replacements: the i-aggregated and
  margin-modulated controls (queue item 1).
- **[union-closed] n-uniform δ²-budgets for the perturbative assembly
  (008's B1/B2 as stated)** (2026-07-31, attempt 012): the averaged
  downward budget's flattening reverses past n ≈ 22 (evidence to
  n = 32), and the tax side is δ-linear, not quadratic. Restate the
  budget object before any assembly re-attempt; 008's conditional
  theorem itself survives with corrected constants.
- **[union-closed] Product-weight functionals / the multiplicative-
  Dirichlet toolbox** (2026-07-31, attempt 010): product-weighted Frankl
  is trivially FALSE (the power set gives frequency x/(1+x) < 1/2 for
  any product weight x < 1), and every functional visible to the
  lcm/divisor encoding is a product-weight functional — the arithmetic
  toolbox aims at an invariant for which the conjecture fails. Pre-empts
  multiplicative repackagings the way the smoothing-insensitive no-go
  pre-empts functional ones.
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
