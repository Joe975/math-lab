# 001 — Ideation sweep: eleven field lenses on Collatz, swarm-executed

- **Problem:** Collatz, `problems/collatz/PROBLEM.md`
- **Date:** 2026-08-15
- **Mode:** informed
- **Type:** survey (ideation sweep per `docs/IDEATE.md`; first attempt on this
  problem, and the first sweep executed on external swarm workers per
  `docs/SWARM.md`)
- **Tools:** `python scripts/mechanisms.py gaps collatz` (gap profile: every
  field untried); `scripts/swarm.py` fanning the brief in
  `explore/ideation-sweep-2026-08/` (template + values, committed) to
  **22 workers: gpt-5.6-luna, effort medium** — 11 field lenses × 2 stances
  (tool-transfer / obstruction-first), 12,019 tokens in / 68,193 out,
  $0.084, ~2.5 min wall clock. Director filter (session model) read all 22
  reports in full. Two load-bearing worker claims re-derived exactly by
  `explore/ideation-sweep-2026-08/barrier_checks.py` (stdlib, seconds,
  deterministic); one was **false as stated** and is corrected below.
- **Sources:** no PDFs consulted. Adjacent-work assignments come from worker
  guesses cross-checked against director memory only — every literature
  attribution below is [T]-grade and must be verified against the actual
  papers before anything is built on it.

## Approach

Collatz had no attempts and the thinnest genuine queue in the ledger (one
writing task). With no prior art, every lens is untried, so a sweep maps the
whole field cheaply before any route is spent — and the problem's own
verification contract (any approach must state which known obstruction it
evades) makes lens reports unusually easy to triage.

Deviation from IDEATE.md's 3–6 lens cap, deliberately: the cap prices lens
reports at one session-agent each. Swarm workers collapse that cost ~1000×,
so all 11 lenses were swept, doubled with an obstruction-first stance on the
theory that for a problem this heavily attacked, **barrier statements would
outvalue routes**. That theory held (see Outcome). The director filter — the
actual bottleneck per SWARM.md — was kept at full strength: every report
read, every load-bearing claim re-derived or discarded.

Worker briefs contained only the tier-0 problem statement and the lens
descriptions; but the brief was authored informed and the filter is
informed, so `mode: informed` (per SWARM.md mode discipline).

## What was done

Each worker proposed 1–3 routes with (core object, which known obstruction
it evades, sub-hour first calculation, kill condition, adjacent-work guess),
or a no-purchase verdict. Filter criteria per IDEATE.md §4: no falsifiable
first step → discard; known dead end in new vocabulary → discard (checked
against the literature rather than the empty index); duplicate of another
worker's route → merged.

**Convergence finding (the sweep's main empirical result): 17 of 22 reports
collapse onto two classical route families.**

**Family A — finite-state residue Lyapunov / ranking certificates**
(11 workers, across seven different lenses): potentials
`V(n) = log n + c_{n mod 2^M}` (or per-residue affine `V_r`), pointwise
descent demanded over the residue automaton of the accelerated map
`A(n) = (3n+1)/2^{v₂(3n+1)}`, searched by LP/SAT/difference constraints.
**KILLED at every finite memory, by re-derivation** (`barrier_checks.py`
check 1, exact to L = 200): `n = 2^L − 1` has its first `L−1` accelerated
valuations all equal to 1, staying ≡ −1 mod shrinking powers of 2, so for
**every** modulus `2^M` the residue automaton contains the self-loop at
`−1 mod 2^M` with log-multiplier `log(3/2) > 0`. Any difference-constraint
system demanding pointwise descent with bounded lookahead is therefore
infeasible at every `M` — no enlargement of modulus or block length escapes,
because the obstruction reproduces at all scales. This is the classical
unbounded-stopping-time fact (Terras 1976) in certificate-infeasibility
form: a **rediscovery, recorded as such**, not a new result. Its value here
is prophylactic — it forecloses the sweep's single most popular route
family, which would otherwise re-enter the queue in every future sweep.
Worker 022 (probability, obstruction-first) was the one report that stated
this kill precisely, naming the 2-adic fixed point −1 as the support.

**Family B — valuation-word cycle Diophantine sieve** (10 workers): the
cycle equation `(2^A − 3^k)·x₀ = C(a₁,…,a_k)`, enumerate valuation words,
sieve by congruence/QR obstructions. This is Böhm–Sontacchi / Steiner /
Simons–de Weger territory [T]; no worker produced a mechanism that evades
the recorded cycle-bound limitation (finite exclusion is not proof), and
several said so themselves. **No purchase as a proof route.** One salvage,
recorded as Lead 2: no worker or director knows of a published *decay-rate
map* for locally-admissible words — quantifying how much of word space
local sieving kills, as calibration data rather than as an attack.

**No-purchase verdicts with certified barriers (the obstruction-first
stance earning its seat):**

1. **Information theory (workers 015/016): parity statistics carry zero
   signal — CHECKED, holds after correction.** For the shortened map
   `T(n) = n/2 | (3n+1)/2`, residue → k-bit parity word is a bijection
   (Terras 1976 [T]; re-verified exhaustively for k ≤ 16, check 2), so the
   parity block of a uniform residue is *exactly* uniform: H = k bits,
   KL = 0. Entropy/KL functionals of parity data are dead on arrival — "KL
   controls mass, not individual trajectories" (worker 015, whose NO
   PURCHASE verdict was the cleanest report in the sweep). **Correction:**
   worker 016 asserted the bijection for the unshortened map, where it is
   FALSE (residues 1 and 3 collide at k = 2 — an odd step forces an even
   successor). The barrier only holds in the shortened form. One of two
   load-bearing worker claims re-derived; one was wrong. That error rate is
   itself sweep-calibration data: **worker arithmetic is not trustable
   without re-derivation, exactly as SWARM.md assumes.**
2. **Graph theory as unlabeled structure (worker 014): connectivity of the
   Collatz digraph is the conjecture itself**, and label-free invariants
   are shared by undecidable generalized maps, so any graph-theoretic
   purchase must retain the arithmetic edge labels. Sharp, and it
   constrains Lead 1 below to labeled/quantitative form.
3. **Algebra (worker 002): no nonconstant polynomial invariant exists** —
   `f(2x) = f(x)` on infinitely many points forces constancy (one line,
   director-checked). Finite-dimensional algebraic invariants are out.

**SPECULATION labels:** every surviving worker route leans on at least one
inline-labelled speculation; none was promoted by the filter.

## Outcome

`MAP`. Scope: a triage of 22 lens reports; two barrier statements re-derived
exactly in range (L ≤ 200, k ≤ 16 respectively — `EVIDENCE` about those
ranges; both are classical facts and carry [T] citations, not novelty
claims); three queue leads below.

**Not claimed:** no route is validated; no progress on the conjecture; no
worker computation was run except the two re-derivations; all
adjacent-work attributions are unverified [T]; the Family-A kill applies to
bounded-lookahead finite-memory certificates over residue states, not to
certificates carrying magnitude information or unbounded adaptive lookahead.

## Why it failed / what survived

Nothing "failed" — the sweep did what sweeps do. What it *revealed*: the
field-lens taxonomy has **weak purchase on Collatz**. Ten of eleven lenses
refracted into the same two classical families, because the problem's only
handles (parity words, the affine cocycle `3^o/2^e`, the cycle equation) are
few and every field grabs the same ones. Contrast union-closed, where the
same sweep design produced genuinely disjoint routes per lens. This
mono-handle behavior is worth knowing before anyone budgets a second sweep.

Survives, reusable:

- the **−1-tower certificate-infeasibility argument** (check 1) — cite it
  whenever a residue-Lyapunov proposal resurfaces;
- the **parity-uniformity barrier** in corrected shortened-map form
  (check 2) — same role for entropy proposals;
- `barrier_checks.py` — both checks exact, seconds, extendable;
- the swarm sweep protocol itself: 22 workers, $0.084, one real correction
  caught at the filter. The economics work; the filter load (~70k tokens
  read) is the true cost, as SWARM.md predicted.

## Leads generated

1. **Graph-structure census of truncated Collatz digraphs** (from workers
   013/014, the sweep's one nonstandard framing — no established
   treewidth/dominator-theoretic Collatz result known to worker or director
   [T], flagged for the literature check). For `B = 2^12 … 2^24`: build
   `G_B` with escape sink, compute (i) dominator tree from 1, (ii) minimum
   directed cuts separating `[1, 2^k]` from vertices whose orbit exceeds
   `2^{k+1}`, (iii) treewidth bounds via flow-cutter, tracking growth in
   `B`. Per worker 014's barrier, keep arithmetic labels: report cut
   composition by residue class. **Falsifiable:** bounded-vs-growing
   interface `|S_k|` is a definite dichotomy; growth kills the
   bounded-interface route (recorded no-go with rates — standalone
   interest), boundedness is a structure worth a follow-up attempt.
2. **Survivor-decay map for the cycle equation** (salvage from Family B):
   enumerate valuation words `k ≤ 20, A ≤ 120`, count words surviving local
   congruence tests by `(k, A)`, fit the decay slope. `EVIDENCE`-scoped
   calibration of what local sieving is worth; kill: no negative slope ⇒
   record that local obstructions do not thin the word space, closing the
   sieve family properly.
3. **Tilted-moment concentration map** (worker 021): exact
   `M_{L,K}(θ) = 2^{-K} Σ e^{θ D_L(n)}` over odd residues mod `2^K`,
   `K = 24, L = 20` — locates quantitatively where the pathwise barrier
   bites (which conditioning destroys the spectral gap). Adjacent to the
   Tao/Krasikov–Lagarias density machinery [T]; value is the map, not a
   proof route.

Not queued: everything in Families A and B as proof routes (killed/no
purchase above); operator-theory pseudospectral probes (019/020 — honest
but their own kill conditions predict their outcome and the information
gained duplicates check 1's).

## References

- Terras (1976), *A stopping time problem on the positive integers* — parity
  bijection, unbounded stopping times. [T]
- Böhm–Sontacchi (1978); Steiner (1977); Simons–de Weger (2005),
  *Theoretical and computational bounds for m-cycles of the 3n+1 problem* —
  Family B territory. [T]
- Lagarias, *The 3x+1 problem: an annotated bibliography* — general
  adjacency source. [T]
- Tao (2019), arXiv:1909.03562 — almost-all barrier context. [T]
- Yolcu–Aaronson–Heule (2021), arXiv:2105.14697 — automated termination
  proving for Collatz-like rewriting; nearest published relative of
  Family A. [T]
- Sweep materials: `explore/ideation-sweep-2026-08/` (template, values,
  `barrier_checks.py`); raw worker reports in `$MATHLAB_OUT/swarm/`
  (uncommitted working space; per-job meta with model + usage recorded).
