# 004 — Skeptic review of 003 (dependent couplings): adversarial verification

- **Problem:** union-closed sets conjecture (Frankl), `problems/union-closed.md`
- **Date:** 2026-07-26
- **Type:** adversarial verification of `003-dependent-couplings.md` (default
  stance: refute). Every load-bearing claim re-derived and/or re-computed with
  an independent implementation (`tools/uc_skeptic.py`, no imports from
  `uc_couplings.py`, natural-log internals, different Sinkhorn
  parametrization; agreement is meaningful because Sinkhorn scalings of a
  fixed kernel are unique, so both codes must converge to the same coupling).
- **Outcome in one line:** the interface and the headline separations are
  REAL, but 003 contains one genuinely false proof claim (mini-theorem (e)
  bracket positivity), one wrong headline constant (recipe ceiling is
  ≈ 0.4315, not ≈ 0.445 — the part-E hunt's grid excluded the extremal
  adversary, and the corrected ceiling sits BELOW 0.44, contradicting 003's
  "does not affect the record-relevant regime p ≤ 0.44"), one
  non-reproducible data file, and several overstated sentences. LIVE status
  survives, with corrections.
- **Tools:** `tools/uc_skeptic.py` (parts L, I, E, S, U, M, K; deterministic;
  checkpoints `data/004_*.json`); full re-run of `tools/uc_couplings.py`
  (bit-identical checkpoints, log in scratch); WebFetch of arXiv:2306.08824
  (reference authenticity).

## Claims attacked

1. Licensing lemma ("exact closure legalizes every coupling; approximate
   closure legalizes none") and the implication (S-coup at p) ⟹ Frankl at p.
2. Derivative lemma (diag⊕iid rung ⟺ Gilmer Conjecture 1; ceiling exactly ψ).
3. Computational separations (Sawin n=300: −2.14 → +10.36; CL smoothed slice
   n=300: −4.62 → +10.85; part-D sign flips; union-closed controls ≤ 0;
   attempted construction of a union-closed family with positive tilt gain).
4. The "ceiling ≈ 0.445" for the single-λ recipe, and the mini-theorem (e)
   (block-adaptive couplings kill every low-mixing-entropy product mixture
   below 1/2).
5. The no-go evasion (Φ_C violates 002's `Φ ≤ H + O(D)` smoothing hypothesis
   genuinely, not just on small instances).

## Refutations found

### R1. Mini-theorem (e): the positivity claim is false at p = 0, and the
### stated consequence is not proved by the displayed bound
Location: 003 §4(e), lines "every bracket being > 0 whenever `p_{k,i} < 1/2`
(and = 0 otherwise)" and "Consequence: no mixture of products with O(1)-many
components … can be a counterexample to (S-coup) for the conditionally-iid
class at any threshold below 1/2."

- The bracket `h(min(2p−p², max(p,1/2))) − h(p)` is **0 at p = 0**, not
  positive (`004_partM.json`). The displayed Gain inequality itself is
  correct (re-derived: `H_π(U) ≥ H(U|k)` plus `H(μ) ≤ H(P) + Σ_k P_k H_k`,
  and the m-formula is the exact per-coordinate optimum of the
  conditionally-iid range m ∈ [p, 2p−p²]); only the positivity clause and
  the consequence fail.
- **Explicit unshielded family:** `μ_n = ½·δ_∅ + ½·Bern(0.6)^⊗n` — two
  product components, all marginals 0.3 < 1/2, H(μ) = Θ(n). The (e) bound is
  −1 bit for every n, and the (e) coupling's actual gain is **exactly 0**
  for every n (its U-law equals μ). So (e) proves nothing for precisely the
  `δ_∅ ⊕ Bern(>1/2)` genre that 003's own part (d) identifies as the hard
  instances — an internal inconsistency: part (e)'s consequence sentence
  claims a shield over a genre that part (e)'s coupling demonstrably does
  not touch.
- The kill is also **not uniform in the component parameters**: the
  n-threshold `H(P)/rate` diverges as components approach {0} ∪ [1/2, 1)
  (computed: n* ≈ 447 at p_low = 10⁻³, ≈ 26,000 at 10⁻⁵, ∞ at 0;
  `004_partM.json`). "Low-mixing-entropy" must be read as
  `H(P) < Σ_k P_k Σ_i bracket_{k,i}` — a joint condition on H(P), n, AND
  parameter non-degeneracy, not on the component count alone. (Ellis's n=2
  example, cited in the consequence's parenthetical, is likewise not killed
  by (e) at n=2 — though it is out of scope anyway, marginals = 1/2.)
- **Repair (this review), restoring the consequence for the exhibited
  genre:** the *half-mixing* coupling — draw `s ~ Bern(p_h)^⊗n`; with
  probability `2P₁` let (A,B) be iid uniform on {∅, s}; else A = B = s — is
  conditionally iid (components are point-mass products), has both marginals
  exactly μ (brute-force check at n=6: marginal error 2·10⁻¹⁶), and gives
  `law(U) = (P₁/2)δ_∅ + (1−P₁/2)Bern(p_h)^⊗n`, hence
  Gain = (P₁/2)·n·h(p_h) − [H₂(P₁) − H₂(P₁/2)] = Θ(n) > 0
  (verified exactly: +4.67 at n=20 up to +48.36 at n=200 for the μ_n above).
  So `δ_∅ ⊕ Bern(p_h)` mixtures are NOT counterexamples to (S-coup) even for
  C₃ — but by a different coupling than (e)'s, and a fully general
  O(1)-component shield remains unproven (see Residual risk).

**Corrected statement of (e):** brackets are > 0 exactly for
`p_{k,i} ∈ (0, 1/2)`; the shield covers mixtures in which coordinates with
`p_{k,i} ∈ (0,1/2)` carry total rate exceeding H(P); components with
coordinates in {0} ∪ [1/2, 1) require different couplings (e.g.
half-mixing), and are handled for the 2-block genre by the repair above.

### R2. The single-λ recipe ceiling is ≈ 0.4315, not ≈ 0.445 — and it IS in
### the "record-relevant regime p ≤ 0.44"
Location: 003 §4(d) ("genuine ceiling ≈ 0.445 against 2-block mixtures",
"none of this affects the record-relevant regime p ≤ 0.44", "the hardest
high-block parameter at cap→1/2 is p_h ≈ 0.62 ≈ φ") and Outcome ¶2
("first genuine obstruction at marginal cap ≈ 0.445 … engine-validated").

- Their model code is **correct**: my independent implementation (own
  Plackett solver checked against bisection, exact best-block-coupling via
  vertex enumeration of the transportation polytope) reproduces their J on
  all 8 reported hardest instances to 7 decimals (`004_partK.json`).
- The error is the **search domain**: `hunt()` floors `p_low` at 0.00025
  (grid starts at 0.001; refinement enforces `0 < q1`), excluding the legal
  adversary `p_low = 0` (δ_∅ is a bona-fide product component, and
  `μ = P_l δ_∅ + P_h Bern(p_h)^⊗n` satisfies every (S-coup) hypothesis).
  The reported "live" values at caps 0.42–0.44 (+0.008, +0.0026, +2·10⁻⁷)
  are all carried by the tilt channel's tail gain on the near-∅ block — a
  mirage that vanishes at `p_low = 0` (at p_low = 0 the tilt channel's sup
  is 0, approached only in the diagonal limit: never a strict gain).
- With `p_low = 0` exactly, 2-block model-killers exist at caps 0.432–0.44
  (all channels strictly negative; table in `004_partK.json`), and the
  boundary has a **closed form**: the binding channel is the
  anti-assortative/best-σ coupling, giving
  `cap*(p_h) = p_h·(2h(p_h) − h(q_h)) / (3h(p_h) − 2h(q_h))`,
  `q_h = 2p_h − p_h²`, minimized as p_h → 1/2⁺ at

      cap* = (2 − h(3/4)) / (2(3 − 2h(3/4))) = 0.431496…

- **Engine corroboration at finite n** (which 003's killers never had): at
  n = 300, μ = 0.1335·δ_∅ + 0.8665·Bern(0.502)^⊗n (cap 0.435), the actual
  Sinkhorn-tilt sweep gives gain ≤ −3.03 at every λ₂ ∈ [−0.4, 2.5] —
  killed, not marginally but decisively.
- A 3-block hunt (400 random instances per cap, channels generalized, exact
  best-σ) found **no** killers below the 2-block boundary at caps
  0.40–0.428, so `δ_∅ ⊕ Bern(1/2+ε)` appears to be the extremal genre and
  0.4315 the model boundary.
- Consequences for 003's prose: "≈ 0.445" is wrong (artifact of the grid
  floor); "none of this affects the record-relevant regime p ≤ 0.44" is
  wrong (0.4315 < 0.44); the "p_h ≈ 0.62 ≈ φ golden-ratio mechanics"
  observation is an artifact of the broken hunt (the corrected boundary
  curve is smooth with binding p_h → 1/2; nothing distinguished happens at
  φ); Outcome ¶2's "engine-validated" was never true of the killer
  instances themselves (only of the anti channel on a different instance).
  The *strategic* conclusion survives: 0.4315 is still far above the
  0.38271 record, and the killers are still separated at Θ(n) by two-scale
  members of the full class (indeed by the conditionally-iid half-mixing
  coupling of R1 — a stronger statement than 003's "trivially, but outside
  the one-parameter family").

### R3. Non-reproducible checkpoint
`data/003_partC_ext.json` (the [ext] +2.07-at-λ₂=6 row of the §4(a) table)
is not produced by any code path in `uc_couplings.py`. All other 003
checkpoints regenerate **bit-identically** from the checked-in tool. The
ext numbers themselves are correct — my independent engine reproduces all
five rows to 3 decimals (−2.144, −0.631, +1.965, +2.065, +1.853) — but the
generating sweep should have been committed.

### R4. Overstated sentence in the licensing discussion
"Approximate closure licenses **nothing** here" (003 §1; also the tool
docstring "legalizes NONE of this") is literally false as a universal: the
diagonal coupling satisfies `H_π(U) = H(μ) ≤ log|F|` for *every* family,
and any coupling supported on non-escaping pairs stays bounded. The correct
and load-bearing statement is existential: approximate closure fails to
legalize the bound *uniformly over couplings* — in particular the specific
gain-carrying tilted couplings put Θ(1) mass on escaping pairs (verified:
11% escape at the CL mini family; the pure-slice tilt sweeps the deleted
middle band). 003 only ever *applies* the lemma to exactly closed families,
and the barrier-evasion argument only needs the existential form, so no
downstream damage — but the sentence should be corrected before anyone
builds on the universal reading.

## Claims that survive (and what was done to break them)

### Licensing lemma (core) — re-derived from scratch
If both marginals of π equal μ = Unif(F), then supp(π) ⊆ supp(μ)² = F×F
(any atom of either coordinate has positive marginal mass), so exact
union-closure gives A∪B ∈ F a.s. and H_π(U) ≤ log|supp| ≤ log|F| = H(μ);
the equality log|F| = H(μ) uses uniformity — the quantifier 002's ladder
relaxed away. Contrapositive gives (S-coup at p) ⟹ Frankl at p, |F| = 1
edge case included. The diagonal-membership and Gain ≥ 0 remarks check out.
Airtight.

**Attempted highest-value refutation (positive tilt gain on an exactly
union-closed family): provably impossible, and empirically negative.**
The support argument above yields gain ≤ 0 for uniform μ on union-closed F
for *any* sub-probability kernel on F×F — independently of Sinkhorn
convergence or marginal fitting. (Corollary worth recording: 003's part-F
controls are *tautological* — they can only ever detect support-handling
bugs, never stress the functional; they are code checks, not evidence.)
Battery run anyway (`004_partU.json`): 12 exactly union-closed families —
closure of the CL slice+top (|F| = 848), a two-scale semilattice product
(up-set × chain on disjoint grounds, n = 12), {|S| ≥ 4} at n = 11
(closure of a full slice, |F| = 1816), a low-entropy 3-generator closure,
six random generator closures, powerset — probed with the tilt sweep AND
~200 random-kernel Sinkhorn couplings (far outside the tilt family). Max
gain over everything: −4.4·10⁻⁵ (all ≤ 0). CONFIRMED.

### Derivative lemma — full equivalence, both directions
Re-derived: f(δ) = H((1−δ)μ + δν) is concave with
f′(0) = CrossEnt(ν,μ) − H(μ) = H(ν) + D(ν‖μ) − H(μ); concavity makes
sup_δ f > f(0) ⟺ f′(0) > 0 an honest equivalence (f′(0) ≤ 0 ⟹ f
nonincreasing; D = ∞ and f′(0) = 0 edge cases checked). So the worry that
only one implication holds is dispelled — the rung genuinely IS Gilmer's
Conjecture 1 per-μ, not just implied by it. Verified numerically in 40
random trials (`004_partL.json`). Ceiling exactly ψ verified end-to-end:
lower bound from f_μ(1) ≥ f_μ(0) and the c = 0 theorem; upper bound by
independent recomputation of Sawin-gadget profiles — at
(n, ū, θ) = (60000, 0.3823, 0.001): max marginal 0.382536, D = 9.966,
f_μ(1) = −41.68, c* = 5.2, matching 002's table exactly (`004_partS.json`).
One nuance worth recording: at these parameters f_μ(1) turns negative only
for n ≳ 15000 (it is +5.67 at n = 5000), so the rung's kill — like all the
gadget kills — is asymptotic in n at fixed (ū, θ). CONFIRMED.

### Computational separations — reproduced twice over
- Full re-run of `uc_couplings.py`: all checkpoints bit-identical
  (determinism confirmed; ~2.5 min).
- Independent implementation (`uc_skeptic.py`): profile engine validated
  against a from-scratch atom-level engine (asymmetric alternating-scaling
  Sinkhorn) at n = 10 to 4·10⁻¹²; B-marginals of the fitted couplings
  explicitly verified (ln-error < 5·10⁻¹²; mass = 1 to 10⁻⁶).
- **Sawin n=300 ū=0.3823 θ=0.02:** iid gain −2.137 by closed form (no
  Sinkhorn: the iid U-law of a union-convolution-closed mixture is itself an
  explicit mixture) and −2.137 by engine; tilt +10.357 at λ₂ = 2.5.
  **Sign flip CONFIRMED.**
- **CL smoothed slice n=300 (w₀=120, t=2⁻¹⁰):** iid −4.623 by closed form
  (exact four-case decomposition slice²/slice·prod/prod²) and by engine;
  tilt +10.850 at λ₂ = 2.5. **Sign flip CONFIRMED.**
- Note: 003's "best tilt gain" values are sweep-edge values (λ₂ = 2.5 was
  the grid max); the sup is higher (+11.11 and +12.29 already at λ₂ = 3.5).
  Understatement only — harmless.
- Part-B pure-slice numbers (+1.25 iid / +407.98 tilted at n = 10⁴,
  p = 0.3823) and all four part-D flips reproduced exactly, including the
  large-n slice rows (+6695.99 at n = 2·10⁵, +79046.96 at n = 2·10⁶) via an
  independent mixture-entropy evaluator; legality of the fixed-intersection
  slice coupling re-derived (existence needs w₀ ≥ n/4 ✓; exchangeability
  gives U ~ Unif(slice n/2) exactly). CONFIRMED.

### No-go evasion — genuine, with closed-form growth
Independent n-sweep at (ū, θ) = (0.3823, 0.001), n = 5000…60000
(`004_partS.json`): the C₃ block-adaptive gain grows exactly linearly,
empirical slope 0.04031 bits/coordinate = the closed form
Σ_k P_k (h(m_k) − h(p_k)), while D(U_iid‖μ) is constant at 9.966 (Sawin's
O(1), reconfirmed). Hence Φ_C₃ − (H_iid + cD) = Θ(n) for every fixed c —
the violation of 002's smoothing hypothesis is structural, not a
small-instance artifact. The consistency argument (Liu's theorem would be
refuted by the gadgets if his functional obeyed the smoothing bound) is
logically sound. Reference authenticity: arXiv:2306.08824 fetched — Jingbo
Liu, "Improving the Lower Bound for the Union-closed Sets Conjecture via
Conditionally IID Coupling", constant ≈ 0.38271, matching 003's use.
CONFIRMED (with the [L] transcription caveat retained — see Residual risk).

### Per-coordinate calculus (§3)
ρ*(p) = p(3p−1)/(1−2p)² and p*(ρ) = (4ρ−1−√(4ρ+1))/(2(4ρ−3)) re-derived
analytically (discriminant 4ρ+1 checked) and confirmed by independent
bisection to 5·10⁻¹² (`004_partL.json`); the constant-odds-ratio claim for
product tilts is an identity (per-coordinate factorization, odds ratio
w₀₀w₁₁/w₀₁w₁₀ = 2^λ), confirmed by direct atom computation at n = 6.
The "unproven SPECULATION" labels on the odds-ratio stability lemma (Gap 1)
are accurate as placed. CONFIRMED.

## Verdict

| # | Claim | Verdict |
|---|-------|---------|
| 1a | Licensing lemma, exact side + (S-coup ⟹ Frankl) | **CONFIRMED** |
| 1b | "Approximate closure legalizes none" | **CORRECTED** (existential, not universal; no downstream damage) |
| 2 | Derivative lemma; diag⊕iid rung ⟺ Gilmer Conj. 1; ceiling exactly ψ | **CONFIRMED** (kill asymptotic in n — nuance recorded) |
| 3a | Sawin n=300 and smoothed-slice n=300 sign flips; part B/D figures | **CONFIRMED** (independent implementation, exact agreement) |
| 3b | Union-closed controls ≤ 0 | **CONFIRMED** (and shown tautological — support argument; controls are code checks, not evidence) |
| 3c | Reproducibility of checkpoints | **CORRECTED** (all bit-identical except `003_partC_ext.json`, which no code generates; its numbers independently verified) |
| 4a | Single-λ recipe ceiling ≈ 0.445; "does not affect p ≤ 0.44"; φ-numerology | **REFUTED** — corrected ceiling cap* = (2−h(3/4))/(2(3−2h(3/4))) = 0.431496 < 0.44, extremal genre δ_∅ ⊕ Bern(1/2+ε); engine-corroborated at n=300; strategic conclusion (ceiling ≫ record) survives |
| 4b | Mini-theorem (e): bracket positivity + "no O(1)-component product mixture below 1/2" consequence | **REFUTED as proved** (bracket(0) = 0; explicit unshielded μ_n); inequality itself correct; consequence repaired for the exhibited genre by the half-mixing C₃ coupling (this review) |
| 5 | No-go evasion (Φ ≰ H + O(D), Θ(n) violation) | **CONFIRMED** |

**Net assessment:** 003's LIVE verdict stands — the two mandated
adversarial families are genuinely separated, the licensing interface is
sound, and the no-go is genuinely evaded. But 003 ships one false proof
sentence, one wrong constant whose corrected value crosses its own
"record-relevant regime" line, and a numerological observation that
dissolves under a correct hunt. The part-E hunt's optimizer (grid floor +
local refinement) failing on its own extremal genre is the cautionary tale:
boundary claims from local search need the boundary adversary included in
the domain.

## Residual risk

- **Liu [L] details.** Only the abstract-level facts (record 0.38271,
  conditionally-iid class) were re-verified here. Clause-level details of
  Liu's construction in 003 §2 still rest on the machine transcription;
  the orchestrator spot-check flagged in 003 remains open.
- **The three proof gaps** (odds-ratio stability, mutual-information tax,
  recipe totality) are untouched by this review; nothing here upgrades the
  interface toward a theorem.
- **General O(1)-mixture shield.** After R1's repair, product mixtures with
  components on {0} ∪ [1/2, 1)ⁿ beyond the 2-block genre have no proven
  C₃ separation; the half-mixing trick plausibly generalizes (pair each
  degenerate-or-high component with samples of itself/others), but this is
  unproven. Until it is, "product mixtures are ruled out wholesale" in
  003's Outcome ¶3 should read "…except mixtures with degenerate/high
  components, ruled out for the 2-block genre".
- **The corrected 0.4315 boundary** is exact within the asymptotic model
  and corroborated by the engine at n = 300, but (like 003's screen) it is
  a necessary-condition boundary: unmodeled λ-schedules could in principle
  revive some killed instances (raising the true recipe ceiling above
  0.4315), and my 3-block hunt was randomized (400 instances/cap) rather
  than exhaustive, so a structured ≥3-block adversary below 0.4315 is not
  excluded — though the channel analysis (mid blocks in (0,1/2) always feed
  the assortative-Plackett channel) makes that unlikely.
- **Shared-reduction risk.** My engine and theirs share the mathematical
  reductions (exchangeability ⟹ profile bookkeeping; Sinkhorn uniqueness).
  Both were validated against brute-force atom enumeration only at n ≤ 10;
  an error in the shared *reduction* at large n would evade both codes.
  The closed-form iid checks (which bypass Sinkhorn entirely) and the exact
  B-marginal residuals make this remote.

## References

- 003-dependent-couplings.md; 002-weighted-kl-ladder.md;
  001-entropy-barrier-map.md (this repo).
- `tools/uc_skeptic.py`; checkpoints `data/004_part{L,I,E_engine,S,U,M,K}.json`.
- Liu, arXiv:2306.08824 (fetched 2026-07-26, abstract verified).
- Sawin arXiv:2211.11504; Chase–Lovett arXiv:2211.11689; Gilmer
  arXiv:2211.09055; Ellis arXiv:2211.12401 (as cited in 001/002).
