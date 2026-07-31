# 010 — Ideation sweep: six field lenses on union-closed

- **Problem:** union-closed, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed
- **Type:** survey (ideation sweep per `docs/IDEATE.md`; first sweep run in the lab)
- **Tools:** `python scripts/mechanisms.py gaps union-closed` (gap profile);
  six parallel lens agents, each given only `PROBLEM.md`, its field lens, the
  spent-tag list, and the dead-end one-liners (per IDEATE.md §3 — deliberately
  not the attempt prose). Orchestrator filter per IDEATE.md §4. Lens agents'
  quick numeric screens are marked unverified below; two load-bearing
  one-line claims were re-derived by hand during filtering.
- **Sources:** none consulted directly; routes naming literature results carry
  an explicit literature-check step. No PDFs read.

## Approach

The cross-pollination layer (2026-07-27 infrastructure) had never been
exercised; the union-closed route just survived its second adversarial cycle
but all its machinery lives in probability/information theory, and gap (a)
died as stated (005/006). That is exactly IDEATE.md's trigger: resist the
reflex "try a variant of the same thing" by sweeping the untried field lenses.

`mechanisms.py gaps union-closed` lists seven untried fields. Six were swept:
**algebra, analysis, dynamics, geometry-topology, graph-theory,
number-theory**. Operator-theory was deferred, not judged: IDEATE.md caps a
useful sweep at six lenses, and its natural object here (a transfer/functional
operator on the family) is exactly what the dynamics lens swept; revisit it if
Lead 2 below goes anywhere.

The alternative — another direct attack on gaps (a′)/(b)/(c) — is running as
attempts 007–009 in this same cycle; the sweep is the hedge against all three
sitting in a local optimum.

## What was done

Each lens agent proposed 1–3 routes with (object, falsifiable first step,
kill condition), or a reasoned no-purchase verdict. The filter discarded
routes that matched a spent mechanism tag, lacked a falsifiable first step,
or restated a recorded dead end in new vocabulary. Two lens claims were
checked by hand during filtering:

1. **Coatom-eigenvalue reformulation (dynamics) — CHECKED, holds.** For the
   union transfer operator `L f(A) = |F|⁻¹ Σ_{B∈F} f(A∪B)` and character
   `χ_C(A) = 1[A ⊆ C]`: since `A∪B ⊆ C ⟺ A ⊆ C and B ⊆ C`, each χ_C is an
   eigenvector with eigenvalue `λ_C = |F ∩ 2^C| / |F|` (one line). Frankl is
   then literally: some coatom eigenvalue `λ_{U∖{i}} ≤ 1/2`. A reframing,
   not a result — but it makes the whole down-set counting profile
   `C ↦ λ_C` the object, and closure constrains that profile worst-case
   (each `F ∩ 2^C` is itself union-closed).

2. **Product-weight no-go (number theory) — the lens agent's witness did not
   reconcile; replaced by a trivial proof of the same statement.** The agent
   reported weighted frequencies 0.326–0.345 for Dirichlet weight 1/n on the
   divisor-encoded power set of {101,103,107}; that number is inconsistent
   with its own construction (1/n weights give frequency 1/(p+1) ≈ 0.01).
   The intended no-go is however true and needs no computation: for ANY
   product weight `w(S) = Π_{i∈S} x_i` with `x_i < 1`, the full power set is
   union-closed and element i's weighted frequency is `x_i/(1+x_i) < 1/2`.
   So **product-weighted Frankl is false**, and every functional the
   multiplicative/Dirichlet-series toolbox can see on the lcm encoding is a
   product-weight functional — the arithmetic toolbox aims at an invariant
   for which the conjecture fails, except at the uniform point where the
   toolbox is empty. (Hand-verified one-line argument; recorded here as the
   structural reason the number-theory lens has no purchase.)

Unverified quick screens by lens agents (rerun before relying on them):
QR-translate adversarial families have max frequency *increasing* with p
(0.69→0.84 for p = 7…23) — closure explosion pushes away from 1/2, route
self-killed; univariate generating functions forget element identity and
collapse to averaged statistics (no-purchase verdict, analysis).

## Outcome

`MAP` — scope: a survey of candidate routes and no-purchase verdicts across
six field lenses; no claim of progress on the conjecture, no bound, nothing
here is a result. Three routes were judged queue-worthy (leads 1–3), six
recorded as secondary (leads 4–9), and five no-purchase verdicts recorded so
the next sweep does not re-run them.

**Not claimed:** none of the routes has been executed; the two hand-checks
above are the only verified statements in this record, and both are one-line
observations, not progress. Kill conditions are as proposed by lens agents
and may need sharpening on contact.

## Why it failed / what survived

Nothing was attempted, so nothing failed; what this record contributes is
triage. The recurring pattern across lenses: every viable route uses the
**exact, worst-case content of pairwise closure** (a constraint system, an
eigenvalue field, a matching structure, a cut structure) — precisely the
resource the 001 obstruction note says the entropy method discards, and none
of it is expressible as a functional `Φ(law(U), μ)`, so the 002 no-go does
not apply on its face. The routes to be most suspicious of are the two
flagged with mandatory pre-checks (leads 7 and 2: smoothing-no-go check and
total-positivity-overlap check respectively).

No-purchase verdicts (recorded to save future sweeps):

- **Girth restriction** (graph theory): closure manufactures C4s in the
  incidence graph (A and A∪B share A); girth ≥ 6 forces trivial families.
  The `girth-restriction` tag does not transfer.
- **Iteration/ergodic averaging** (dynamics): union is idempotent, so the
  union walk collapses to t-fold iid unions = the dead k-wise-union route;
  ergodic averages are the average-case entropy method again. Only the exact
  spectral side survives (lead 2).
- **Multiplicative/Dirichlet arithmetic** (number theory): dead by the
  product-weight no-go above. `residue-class-covering` and
  `mordell-qr-obstruction` have nothing to cover here.
- **Univariate generating functions** (analysis): forget element identity;
  averaged consequences only.
- **Compression dynamics** (dynamics): compression does not preserve
  union-closure (folklore); nothing falsifiable beyond re-deriving that.

## Leads generated

Queue-worthy:

1. **Pairwise-closure constraint certification** (algebra + analysis,
   merged). Closure is 2^n exact constraints "Σ_{x∨y=z} f(x)f(y) = 0 for
   z ∉ supp(f)", of which the entropy method keeps one scalar shadow. For
   every union-closed family on n ≤ 4 (all 4958) and a sample at n = 5,
   solve the LP/degree-2 SOS dual: what min-max-frequency bound do the
   pairwise constraints alone certify? Kill: certified bound converges to
   ≈ 0.382 (the constraints carry nothing beyond the entropy scalar);
   deliverable either way is the fractional witness.
2. **Union-transfer-operator eigenvalue field** (dynamics). With
   `λ_C = |F ∩ 2^C|/|F|` (checked above): test log-supermodularity
   `λ_C λ_{C'} ≤ λ_{C∪C'} λ_{C∩C'}` [SPECULATION] on the full n ≤ 5 census;
   if it holds, all-coatoms-> ½ forces `λ_C > 2^{-(n-|C|)}` for all C — then
   check that against the exact identity `Σ_C |F∩2^C| = Σ_{A∈F} 2^{n−|A|}`
   for a contradiction. Kill: one small log-supermodularity counterexample.
   Mandatory pre-check: read 003/005's total-positivity material first —
   correlation-inequality-adjacent, must not be a spent tag in disguise.
3. **Bipartite MIS decomposition** (graph theory). In the
   Bruhn–Charbit–Schaudt–Telle bipartite reformulation (literature check of
   the exact statement is step 0), define r(G) = min-vertex share of maximal
   stable sets per side. `geng -b -c` all connected bipartite graphs to
   n = 12, compute r(G), tabulate connectivity/treewidth of the extremal
   graphs, and directly test whether r composes across 1- and 2-cuts. Kill:
   extremal graphs already contain cut vertices / treewidth-2 examples, or an
   explicit 2-cut compositionality counterexample. Prize: "minimal
   counterexample is 3-connected"-type reductions — worst-case assembly with
   no mutual-information tax to pay.

Secondary (recorded, not queued this cycle):

4. **Topological Hall on union-injection graphs** (geometry-topology): for
   each element x, bipartite G_x on F_¬x ∪ F_x with S ~ S∪T, T ∈ F_x; a
   matching saturating F_¬x gives frequency ≥ ½ with a witness. Census
   n ≤ 5: does some x always satisfy Hall? Literature check first (injection
   variants of Frankl are studied; rediscovery is not a result).
5. **Union-escape graph expansion** (graph theory): edges = atypical-union
   pairs; test min-cut/algebraic connectivity vs max frequency on n ≤ 6
   census. Aimed directly at the 001 obstruction ("worst-case closure of
   atypical pairs").
6. **Crosscut/nerve truncation on Poonen's lattice form**
   (geometry-topology): inclusion–exclusion over atom filters is exact and
   internal (∩ of filters is a filter by closure); check what 2-term
   Bonferroni truncations force on small censuses and adversary families.
7. **Trilinear OR-inequality** (analysis): sharp constants in
   `E[f(x)g(y)h(x∨y)] ≥ ‖f‖_a‖g‖_b‖h‖_c` on the two-point space, then
   tensorize. Mandatory step 0: check the functional is smoothing-SENSITIVE
   and that two-point extremizers do not reduce it to Gilmer's inequality —
   both failure modes land in recorded dead ends.
8. **Transitive-symmetry families** (algebra): under a transitive
   automorphism group Frankl ⇔ average set size ≥ n/2 [check literature
   status of the transitive case first — SPECULATION that it is open].
   Cheap census of C_n-invariant families n ≤ 7 via
   `group-orbit-enumeration` transfer (mahler-4d/001).
9. **Lattice Möbius/congruence mining** (algebra + number theory, merged):
   tabulate μ(0̂,1̂), alternating sums, |F| mod m against max frequency on
   the n ≤ 5 census (`representation-count-mining` transfer,
   erdos-straus/002). Low prior; one afternoon; decisive kill (no separation
   beyond |F| and n).

## References

- Prior attempts: `problems/union-closed/attempts/001` (obstruction note:
  entropy method's only use of closure), `002` (smoothing-insensitive no-go),
  `003`–`006` (live route and its gaps — this sweep is the hedge against
  them being a local optimum).
- Procedure: `docs/IDEATE.md`; gap profile from `scripts/mechanisms.py`.
- Literature named by lens agents, to be verified at execution time, not
  consulted here: Bruhn–Charbit–Schaudt–Telle bipartite reformulation
  (lead 3); Poonen's lattice equivalence (lead 6); Aharoni–Berger–Meshulam
  topological Hall machinery (lead 4). Marked unread.
