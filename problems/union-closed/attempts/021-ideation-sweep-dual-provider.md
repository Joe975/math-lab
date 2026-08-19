# 021 — Dual-provider ideation sweep: seven lenses × two model families

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** survey (ideation sweep per `docs/IDEATE.md`, swarm-executed per
  `docs/SWARM.md`)
- **Tools:** `scripts/swarm.py`; briefs committed as
  `explore/swarm020_ideation_brief.md` + `swarm020_ideation_values.txt`;
  7 lenses × 2 families = 14 workers (gpt-5.6-luna effort low,
  gemini-3.7-flash effort low), ≈ $0.06; outputs under
  `$MATHLAB_OUT/swarm/020-ideation-{gpt,gem}/` (prompt hashes in
  `.meta.json`). Director triage against `mechanisms.json` and the
  STATUS.md dead-end list. Two cheap census checks run inline (below).
- **Sources:** none consulted directly; worker-recalled citations are
  marked [T] where kept.

## Approach

Second ideation sweep on this problem (first: 010, 2026-07-31, six
lenses). This one covers the five fields the gap profile lists as never
tried here — dynamics, geometry-topology, graph-theory, number-theory,
operator-theory — plus two narrow-slice retries (analysis beyond
coupling inequalities; probability beyond dependence engineering), and
runs every lens on BOTH external model families. The family duplication
is deliberate: per `docs/SWARM.md`, cross-family agreement is a triage
signal (and per-family blind spots differ). Workers saw the tier-0
problem statement, the spent mechanism-tag list, and one-line dead-end
summaries — not the records.

## What was done

Each worker returned 1–3 routes (object / first calculation / kill
condition / escape argument) or an explicit no-purchase verdict. 14/14
returned cleanly. Director triage discarded routes matching spent tags
or recorded obstructions, then ranked the rest by cross-family
convergence and falsifiability. Two verdicts were checkable immediately
and are settled inline here.

### Cross-family convergences (the sweep's main signal)

1. **Union-transfer operator, sharpened** (dynamics + operator-theory,
   both families — four independent reports land on the same object).
   For μ on F define (T_μ f)(A) = E_{B~μ}[f(A ∪ B)] on ℓ²(F, μ). The
   sharpening, from the gemini operator-theory worker, is an exact and
   elementary fact the queued 010 lead lacked: **φ_i(A) = 1 − 1_{i∈A} is
   an eigenfunction of T_μ with eigenvalue 1 − p_i** (proof: E_B[(1 −
   1_{i∈A∪B})] = (1 − 1_{i∈A})(1 − p_i); verified by hand). So Frankl for
   the uniform measure is exactly: *some coordinate eigenvalue of T_μ
   lies in [0, 1/2]* — a spectral-gap statement about a non-normal
   operator whose coordinate eigenvalues are the frequencies themselves.
   This consolidates queue 3(ii) rather than replacing it: the value
   would come from spectral machinery (numerical range, compressions,
   resolvent bounds) controlling min_i λ_i, which no worker could supply.
   First calculation as queued in 010, now with the eigenfunction basis
   as the starting point.
2. **Möbius/zeta incidence structure** (three gpt reports + one gemini).
   LP search over low-order zeta/Möbius statistics on the n ≤ 5 census
   for any valid inequality implying the 1/2 bound. Concrete and cheap;
   kill condition well-stated (no nontrivial LP-feasible inequality on
   the exact census = no purchase). Folds into the queued 010 lead 3(i)
   LP machinery as an additional statistic family.
3. **Lorentzian / strong-Rayleigh generating polynomials** (probability,
   both families). The naive version is **settled no-purchase inline**:
   ultra-log-concavity of the slice sequence f_k fails trivially for
   union-closed families with internal zeros — e.g. F = 2^{{1,2}} ∪
   {[6]} is union-closed with f = (1,2,1,0,0,0,1). What survives is the
   feasibility question the gpt worker posed: does every union-closed
   family admit SOME strongly-Rayleigh measure with all marginals below
   1/2 (if yes, negative dependence alone cannot prove Frankl; if some
   family admits none, that family class is where negative dependence
   bites). Falsifiable by SDP/feasibility search on the n ≤ 5 census.
4. **Fourier/OR-convolution level-1 statements** (analysis, both
   families). Both propose census tests of level-1 spectral inequalities
   under the OR semigroup; both flag that Walsh analysis diagonalizes
   XOR, not OR, so the right harmonic analysis may not exist. The
   gemini majorization variant Φ_q = Σ p_i^q ≥ 1 + (n−1)(1/2)^q is
   **settled no-purchase inline**: it fails on the same F = 2^{{1,2}} ∪
   {[6]} example — n = 6, |F| = 5, frequency vector p = (3/5, 3/5, 1/5,
   1/5, 1/5, 1/5), so Φ₂ = 2·(9/25) + 4·(1/25) = 0.88 < 2.25 =
   1 + 5·(1/2)²; the family satisfies Frankl (p₁ = 3/5 ≥ 1/2), so the
   inequality is not even a valid necessary property.

### Non-convergent but retained

- **Discrete-Morse / fiber matching** (gemini geometry-topology): test on
  the n ≤ 5 census whether an acyclic Hasse-diagram matching always pairs
  F \ F_i* into F_i*. This is a strengthening in the folklore
  injection-conjecture family; cheap, falsifiable, likely to fail fast —
  worth one census pass at most.
- **Incidence-graph separators / treewidth census** (gpt graph-theory):
  compute separator/conductance invariants of near-extremal families; the
  stated kill condition (invariants generic ⇒ no purchase) makes it a
  bounded diagnostic, not a program.
- **Exchangeability at degrees 3–4** (gpt probability): Hoeffding
  components of iid samples beyond pairwise — the one probabilistic
  route not obviously inside the spent coupling family; speculative, but
  its kill condition (deficits vanish under replication) is concrete.

### No-purchase verdicts (recorded so the next sweep skips them)

- **Number theory, both families independently**: no purchase. Both
  argue the same way — union-closure has no multiplicative/congruence
  structure to grip; any counting-mod-p statistic is blind to the
  frequency vector. Cross-family agreement on a negative is the
  strongest no-purchase signal this protocol can produce.
- **Tropical convexity** (gemini geometry-topology): F is its own
  tropical vertex set, but no mechanism links tropical facets to
  frequencies; kill condition was not falsifiable as stated. Discarded.
- **Slice-sequence ULC / Lorentzian raw form**: killed inline (above).
- **Majorization Φ_q form**: killed inline (above).
- Worker-proposed compression/shifting routes duplicate the classical
  literature (shifts do not preserve union-closure) and were discarded
  at triage.

## Outcome

`MAP`. Scope: seven field lenses, two model families, one sweep. Three
queue-worthy consolidations (transfer-operator eigenfunction sharpening;
Möbius/zeta LP statistics; strong-Rayleigh feasibility question), two
inline no-purchase settlements with explicit counterexample family
F = 2^{{1,2}} ∪ {[6]}, one cross-family no-purchase (number theory), and
three bounded diagnostics retained. No claim of progress on the
conjecture; no worker text entered the repo unverified.

## Why it failed / what survived

The sweep did what sweeps are for: it did not produce a new attack
family, but it (a) upgraded a queued lead with an exact structural fact
(the eigenfunction identity — the kind of small true statement that
makes a vague lead concrete), (b) closed four idea-shapes cheaply, two
with an explicit counterexample, and (c) demonstrated that cross-family
convergence/divergence is a usable triage axis: every convergent object
was implementable, and both inline kills came from routes only one
family proposed. The recurring failure mode across lenses is the same
one 010 recorded: invariants of F that forget multiplicity (homology,
tropical hulls, downward closures) cannot see the frequency vector at
all.

## Leads generated

1. **[queue-merge] Transfer-operator lead (010 3(ii)), restated:** on the
   n ≤ 5 census, compute the full spectrum of T_μ alongside the
   coordinate eigenvalues 1 − p_i; test whether any spectral statistic
   of the non-coordinate part (numerical radius of the compression off
   span{φ_i}, Jordan structure) separates families with min eigenvalue
   ≤ 1/2 from hypothetical violators. Kill: coordinate and
   non-coordinate spectra vary independently.
2. **[queue-merge] Add zeta/Möbius statistics to the 010 3(i) LP pass.**
3. **[new, bounded] Strong-Rayleigh feasibility on the n ≤ 5 census:**
   for each family, does a measure with Rayleigh inequalities and all
   marginals ≤ 0.49 exist? Either outcome is informative (see above).
4. **[new, one census pass] Fiber-matching test** on n ≤ 5.
5. The two inline counterexample settlements should be cited by any
   future sweep that touches log-concavity or majorization shapes here.

## References

- This repo: `docs/IDEATE.md`, `docs/SWARM.md`, attempt 010 (first
  sweep), `mechanisms.json` gap profile for union-closed.
- Worker-recalled, unverified [T]: Plackett 1965; Karlin, *Total
  Positivity* 1968; Brändén–Huh Lorentzian polynomials; Forman /
  Babson–Hersh discrete Morse theory.
