# 029 — Reviewer-independence pass on 025–028 (branch kills, certification, adaptive cell, recipe v2)

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** adversarial review (default stance: refute) of 025 + 026 +
  027 + 028 — the queue-(a′) reviewer pass those four records point at.
  **Independence, stated exactly:** the reviewer is a fresh-context
  subagent with **no conversation state shared with the author
  sessions**, working only from what is committed in the repo — but it
  was **spawned from the author's session with this review brief, and is
  the same model family**. Whether that meets the 024 fresh-session
  standard (024 was a separately started session) is for the human to
  judge; the implementation-level independence below (own code, own
  derivations) does not depend on that judgment.
- **Outcome in one line:** everything load-bearing survives — Lemmas
  EM/EM-coverage re-derive by hand in full, every re-implemented number
  agrees (38/38 checks), all six scripts re-run byte-identically and
  each certification kit *alone* certifies all of 026's statements — but
  026's dual-kit soundness sentence ("certified under either kit's
  audit") is wrong as an argument about intersections, 027's claimed
  n ≤ 8 refinement never ran (its code path is dead on the survival
  branch; true scanned range is n ≤ 6), and 025 misquotes its own
  skeptic log (agreement is ≤ 6.4e-13, not ≤ 4e-15).
- **Tools:** re-runs, unmodified, of `explore/uc_branch_attack.py`,
  `uc_branch_skeptic.py` (exit 0), `uc_branch_certify.py`,
  `uc_adaptive_cell_attack.py`, `uc_adaptive_cell_skeptic.py` (exit 0),
  `uc_recipe_v2.py` — every committed checkpoint reproduced
  byte-identically (clean `git status` after each);
  `explore/uc_reviewer029_reimpl.py` (new; stdlib only; shares **no
  code** with the modules under review — own coupling builders from the
  records' verbal definitions, own chain-rule evaluator computing
  z = P(U_i = 1 | prefix) with exact-Fraction coupling accumulation at
  the rational instances, own profile-sum entropies, own IPF Sinkhorn;
  output `data/uc_reviewer029_reimpl_out.txt`, exit 0);
  `explore/uc_reviewer029_certaudit.py` (new; imports the two log₂ kits
  exactly as 026 does, re-implements 026's interval algebra, and adds a
  from-scratch 200-bit exact-Fraction log₂ reference; output
  `data/uc_reviewer029_certaudit_out.txt`, exit 0). Python 3.11.15; no
  swarm, no external workers.
- **Sources:** none.

## Claims attacked

1. **Lemma EM** (025 §1): ST = 0 for the ε-mixing coupling family with
   arbitrary equal-marginal (ε_A, ε_B) joints and arbitrary s-laws;
   the concavity corollary and the surprisal criterion.
2. **Lemma EM-coverage** (025 §4): the G-criterion, monotonicity in P₁,
   the P₁ ≥ 1/2 case, the P₁ < 1/2 n ≥ 3 case, and every constant in it
   (h(2·0.38271), the log₂ ratio bound, max −x log₂ x, the n ≥ 3
   margin); the float-level n ∈ {1, 2} boundary sweeps.
3. **The half-mixing kill** (025 §2, 026 C1/C2): in-regime CR < 0 at the
   sliver instances; the repair; the P₁ window (0.61729, 2/3].
4. **The block kills** (025 §3, 026 C3/C4): the ψ component threshold,
   Θ(n) negativity, the comonotone CR = 0 identity, the ∅-mixing rescue.
5. **026's certification logic**: the interval helpers, the dual-kit
   intersection soundness argument, the scoping of the dependence on
   009 fact (i) and Lemma EM, the kits themselves.
6. **027's scan and probe claims**: the 510-instance survival, the
   tightest margins, the CR/H collapse, latent-mixing values and their
   EM limits, the stated scan range.
7. **028's battery**: the 19-instance table, the recipe-v2
   specialization claims, and whether the record's numbers ship code.
8. **Records-vs-scripts fidelity** for all four (overclaiming audit,
   status vocabulary per AGENTS.md).

## Refutations found

### R1. 026's dual-kit soundness sentence is wrong as an argument;
### the conclusion is rescued by single-kit certification (done here)

026 (Approach; echoed in "Why it failed") claims that intersecting the
two kits' enclosures means "the result is certified under either kit's
audit." **That inference is invalid for intersections.** If exactly one
kit is sound, the intersection need not contain the true value: a
tight-but-wrong enclosure from the unsound kit can pull the intersection
off the truth without ever becoming disjoint (e.g. truth 0.5, sound kit
[0.4, 0.6], unsound kit [0.55, 0.58] → intersection [0.55, 0.58]:
non-disjoint, no abort, wrong). Disjointness proves *at least one kit
unsound*; non-disjointness does **not** prove the intersection contains
the truth. An intersected-enclosure certificate is certified under the
**conjunction** of the two audits — strictly weaker than "either."

The conclusion 026 wanted is recoverable, and this review establishes
it: `uc_reviewer029_certaudit.py` re-runs the full 026 computation three
times — kit A (022 digit-extraction) **alone**, kit B (016 atanh-series)
**alone**, and intersected — with 026's interval algebra re-implemented
independently. **Each kit alone certifies all twelve statements** (the
nine named C1–C4 rows plus the three H(μ) lower bounds), and the
intersected enclosures reproduce 026's committed ones to the last digit
(9/9). So every 026 certificate is genuinely certified under either
kit's audit — by the single-kit computations, not by the intersection
argument. Kit widths: A ≤ 2.5e-32, B ≤ 2.6e-14; margins ~1e-3 dwarf
both. **Corrected pattern for future certifications:** run the kits
separately and require each to certify alone (the disjunction 026
claimed); use the intersection only as a mutual-consistency check.

The kits themselves were spot-tested against a from-scratch reference
(exact-Fraction square-and-shift bit extraction, 200 fractional bits,
floor-rounded only at 4096-bit denominators with the amplification
bounded): both bracket the reference on 40 random rationals, bracket
exact powers of two, and satisfy log₂(r) + log₂(1/r) ∋ 0. 026's own
interval helpers (`ivl_add`/`ivl_scale`/`ivl_sub`/`xlog2x_ivl`/
`H_profile`) are exact-Fraction endpoint operations with the sign split
in `ivl_scale` correct — audited by hand, and my re-implementation of
the same algebra reproduces all committed enclosures.

### R2. 027's "refinement at n ≤ 8" never ran; the scanned range is n ≤ 6

027's Outcome claims "refinement of the tightest instances at n ≤ 8 in
the checkpoint," and its index `range` says "refinement n <= 8." **No
such refinement exists.** In `uc_adaptive_cell_attack.py` the
refinement block (`OUT["refine_n"]`, n up to 8) sits inside
`if neg:` — it executes only when a *negative* instance is found. The
scan found none, so the code path is dead, `refine_n` is absent from
`data/adaptive_cell.json`, and nothing in 027 was ever evaluated at
n = 7 or 8. The correct scope of 027's survival claim is **n ≤ 6**
(510 instances) exactly as its own §"What was done" opens with.

Gap filled here rather than left hanging: my own builder/evaluator runs
the tightest scan instance (p_lo = 0.001, p_hi = 0.85, P_hi = 0.42) at
n = 3 and n = 8: CR = +0.001217 and +0.016046 — positive, so the
survival *extends* to n = 8 at that instance (one instance, not a
grid; the range correction stands).

### R3. 025 misquotes its own skeptic agreement: ≤ 6.4e-13, not ≤ 4e-15

025 §2 ("The skeptic reproduced the three spot instances … to 4e-15
agreement (S1)"), its Outcome ("two independent implementations
agreeing to ≤ 4e-15"), and its index one-line ("two independent
implementations agree to 4e-15") all misreport the committed skeptic
log: `data/branch_skeptic_out.txt` S1 deltas are 3.6e-16 (n = 2),
−3.35e-15 (n = 6), and **6.37e-13 (n = 16)**. The correct statement is
≤ 6.4e-13 — still ~10 orders below the 1.5e-2 kill margins, so nothing
downstream moves, but the quoted figure covers only two of the three
instances.

### R4. 028's cr_chain cross-checks are transcript-only (no committed code)

028 §3 reports inline independent checks of the genuinely new battery
numbers "via the 025 skeptic's `cr_chain` (transcript values in §3)" —
but `uc_recipe_v2.py` contains no such calls and no committed script
produces them, so those specific cross-check values fail
CONTRIBUTING's reproducibility bar as shipped. Discharged here with my
own (third) implementation: the 2block n = 4 closed-form value is
reproduced full-history at q\* with CR = +0.290240, ST = −1.3e-15, and
the 3block n = 8 / sliver n = 16 rows re-derive independently (R7
below) — 028's quoted values were right; they just shipped without
code.

### Minor misreports (recorded, no downstream effect)

- 025 Outcome: "slopes matching the per-coordinate prediction to 3
  decimal places at n = 40" — true for 2block only (−0.0993 vs
  −0.0993); 3block is −0.1271 vs prediction −0.1244 (2 dp), hi+lo
  −0.0526 vs −0.0521. The finite difference at n = 24…40 still carries
  the O(1) transient; Θ(n) negativity is unaffected.
- 025 §4: repair values "CR(q = 1/2) ≈ +0.067…+0.073 everywhere the
  stated branch is negative (A2)" — the A2 n = 1 row is +0.0599; the
  range should read +0.060…+0.073.
- 025 §2: "H(μ) ≈ 0.97" for the A2 family — the n = 1 row is 0.940.
- `uc_recipe_v2.py` docstring says "~90 s"; actual ~0.7 s (the record's
  own "~1 s" is right).

## Claims that survive (and what was done to break them)

### Lemma EM — re-derived by hand, in full; stress-tested on fresh instances

The proof is correct, and its crux survives scrutiny: with A = ε_A·s,
B = ε_B·s, (ε_A, ε_B) ⊥ s, the three nonzero-prefix cell types
(a, 0), (0, a), (a, a) each fix (ε_A, ε_B) and force s_{<i} = u, so by
independence the conditional law of U_i = s_i in each cell is
law(s_i | s_{<i} = u) — the same as conditioning on U_{<i} = u alone —
and the zero cell is the *identical σ-field atom* {A_{<i} = B_{<i} = 0}
= {(ε_A∨ε_B)·s_{<i} = 0} = {U_{<i} = 0}, so its ST summand vanishes
trivially. law(U) = qδ∅ ⊕ (1−q)ν since P(ε_A∨ε_B = 0) = q. Hence
ST = 0 and CR = Gain = F(q) − F(P₁). The Fréchet range
q ∈ [max(0, 2P₁−1), P₁] is exactly right for equal-marginal ε-joints.
The corollary also re-derives: F′(w) = −log₂ m_w(∅) + E_ν[log₂ m_w]
(the 1/ln 2 terms cancel because both measures are probability
measures), F strictly concave for ν ≠ δ∅, and with P₁ ∈ (0, 1) the
q-interval has nonempty interior below P₁, making
sup_q CR > 0 ⟺ F′(P₁) < 0 an honest iff — the surprisal criterion as
stated. Attack attempted: fresh ST = 0 instances 025/026 never ran
(non-product two-component s-law; a 3-atom arbitrary s-law with q at
the Fréchet lower edge; q strictly between the edges at n = 4), all
through my own builder and evaluator: max |ST| = 8.9e-16. No crack.

### Lemma EM-coverage — every constant re-derived; the two float legs re-swept

- The G-criterion re-derives: F′(P₁) = −G with
  G = (1−b)·log₂(m₀/(1−P₁)) + n·h(ph) + b·log₂ b exactly (checked
  symbolically against my own E_ν[log₂ μ] expansion).
- Monotonicity in P₁: m₀ increasing, 1−P₁ decreasing, so the bracket is
  strictly increasing; the other terms are P₁-free. Correct.
- **P₁ ≥ 1/2 case**: bracket ≥ 0 from m₀ ≥ P₁ ≥ 1−P₁. The sub-claim
  n·h(ph) + b·log₂ b > 0 for ph ∈ (0, 1) re-proves cleanly:
  −b·log₂ b = n(1−ph)ⁿ(−log₂(1−ph)) ≤ n(1−ph)(−log₂(1−ph)) < n·h(ph),
  the last step because h(ph) adds the positive −ph·log₂ ph term.
  Rigorous, as claimed.
- **P₁ < 1/2, n ≥ 3 case**: in-regime with ph ≥ 1/2 forces
  P₁ > 1 − 2(0.38271) = 0.23458 and (using P₁ < 1/2) ph < 0.76542.
  Constants, recomputed from scratch: h(2·0.38271) = h(0.76542) =
  **0.785910** ✓; bracket ≥ −log₂(0.76542/0.23458) = **−1.706172**
  (025's −1.7062 ✓; the (1−b) < 1 factor only shrinks the negative
  bracket, so the bound is valid); max_x(−x log₂ x) = log₂(e)/e =
  **0.530738** ✓; margin 3(0.785910) − 0.530738 − 1.706172 =
  **+0.120821** ✓ (025: +0.1208). Rigorous.
- **n ∈ {1, 2} legs**: my own 20 001-point sweep of my own G on the
  boundary P₁ = 1 − 0.38271/ph reproduces min G = +0.344848 (n = 1)
  and +0.865786 (n = 2), both at the ph = 0.5 edge, matching the
  engine's 200 001-point values to 1e-3 grid effects. These legs remain
  float-level, exactly as 025 labels them.
- One logical nicety checked and found sound: 025's "A3 grids found
  zero in-regime points with sup_q CR ≤ 0" — A3 evaluates the q-sweep
  only at half-mixing-failure points, but at every other grid point
  sup_q CR ≥ CR(q = P₁/2) > 0 trivially, so the claim follows.

**Consequence:** 025's own criterion for upgrading EM/EM-coverage from
candidate-VERIFIED was a reviewer-level re-derivation. This pass
provides one, subject to the independence caveat in the front matter:
both lemmas are re-derived by hand with no crack found, the n ∈ {1, 2}
legs staying computer-checked as stated.

### The half-mixing kill and repair — re-implemented exactly

Exact-Fraction couplings built from the ε-mixing definition at 026's
rational instances, evaluated by my own chain-rule accounting:
CR_hm(n = 2, ph = 1987/2000) = −0.002267927906 and CR_hm(n = 6,
ph = 1997/2000) = −0.007133120856, matching 025's A5 floats to 1e-11
and sitting on 026's enclosure midpoints to < 1e-11, with ST ≤ 2.2e-16,
marginal deviation exactly 0, marginals in-regime. Repairs at q = 1/2:
+0.073075033229 / +0.070953077915 — same agreement. The A2 family
reproduces at n = 1, 2, 6 full-history and at n = 16, 100 via my own
profile closed form (−0.029388 / −0.015386 / −0.015321 / −0.015321 /
−0.015321). The mechanism re-derives: fixed cost h(P₁) − h(P₁/2) > 0
for P₁ ∈ (0, 2/3) against a gain that vanishes as ph → 1, with the
in-regime window P₁ ∈ (0.61729, 2/3] exactly as stated (1 − 0.38271
and the Fréchet validity of q = P₁/2). A3's failure-region percentages
(3.91 / 0.15 / 0.03 at n = 1/4/16) and worst −0.067 verified against
the checkpoint.

### The block kills, ψ threshold, comonotone identity, rescue

- ψ is exact algebra: h(2p−p²) = h(p) with 2p−p² ≠ p forces
  2p − p² = 1 − p, i.e. p² − 3p + 1 = 0, root (3−√5)/2. Re-derived;
  no numerics needed.
- 2block Gain(n = 8) = −0.7940026768 by my own direct-enumeration
  entropies on the exact-Fraction mixture — matches 025's checkpoint
  and 026's enclosure to 1e-10. Full-history CR of my own iid-block
  build sits at Gain with ST = −4.4e-16 (for iid-given-k the U-law
  refinement costs nothing here), CR < 0 confirmed.
- Fact (i) (CR ≤ Gain) re-derived: σ(A_{<i}, B_{<i}) refines σ(U_{<i}),
  conditioning reduces entropy termwise, ST ≥ 0. 026's use of it in C3
  is correctly scoped (a proved lemma, cited, not recomputed).
- Comonotone CR = 0 re-derived as the chain-rule identity for U = A.
- The 3block ∅-mixing rescue: my exact-Fraction build at q = 1/5 gives
  CR = +0.692968847018, on 026's midpoint to 1e-11, ST = +4.4e-16.
  026's C4 scoping is honest (certified quantity is the entropy
  difference; Lemma EM supplies CR = Gain — a dependence this review's
  hand re-derivation now discharges).

### 027's survival, weakness finding, and latent-mixing — reproduced

Four tightest scan instances re-evaluated with my own adaptive-cell
builder (from 009 part D's rule as quoted): all CR > 0, agreement
1e-10. Scaling probe at p_lo ∈ {1e-5, 1e-3}, n = 6: adaptive CR
reproduces to 1e-9 with CR/H = 1.0e-5 at bounded H = 2.52 — the
CR/H → 0 collapse is real. Latent-mixing values reproduce from my own
build (marginal symmetry exact by construction, checked at 0) to 1e-9,
and the EM limits (+0.221002 / +1.189593 at n = 2/6) re-derive from my
profile form. 027's scope caveats (measured only; no ST = 0 lemma;
SPECULATION label on ST = O(p_lo·n)) are accurate as placed.

### 028's recipe v2 and battery — statement checks out, 6/19 rows re-derived

The three specialization claims re-derive from the definitions:
q = P_light makes the cross cells vanish leaving exactly the
adaptive-block coupling; a δ∅ light component turns the equal-label
light cell into the ∅ atom and the family into 025's ∅-mixing (Lemma
EM applies); and optimized q dominates half-mixing because q = P₁/2 is
an interior point of the swept range. Battery: the full 19/19 re-run
reproduces the committed table; independently re-derived here — sliver
n = 2 (+0.074173) and n = 16 (+0.067901) and 2block n = 4 (+0.290240)
and 3block n = 8 (+1.005836) via my profile closed form, the 2block
row also full-history at q\* (ST = −1.3e-15), and both G-gen floor rows
(floor:mmabskill +0.074203, floor:windowkill +0.223165) via **my own
IPF Sinkhorn** on the committed endpoint measures, agreeing to < 1e-6.
028's Outcome caveats (finite battery, measured-not-proved
latent-mixing, no totality) match what the script does; LIVE is the
right status word.

### Reproducibility — full

All six scripts re-run unmodified on Python 3.11.15: both skeptics exit
0 with zero refutations; `branch_attack.json`, `branch_skeptic_out.txt`,
`branch_certify.json`, `adaptive_cell.json`,
`adaptive_cell_skeptic_out.txt`, `recipe_v2_battery.json` all reproduce
byte-identically (the engines are genuinely deterministic — no RNG, no
hash-order dependence observed — in contrast to 024's finding on 023's
engine). 026: 98 dual-kit calls, zero disagreements, max width 2.47e-32,
H(μ) lower bounds 0.9832/0.9777/0.9687 (the record's "> 0.968" holds,
barely, at n = 16).

## Verdict

| # | Claim | Verdict |
|---|-------|---------|
| 1 | Lemma EM + corollary + surprisal criterion | **CONFIRMED** (hand re-derivation; fresh ST = 0 stress instances) |
| 2 | EM-coverage: criterion, monotonicity, both analytic cases, all constants | **CONFIRMED** (constants re-derived: 0.785910 / 1.706172 / 0.530738 / +0.120821); n ∈ {1, 2} legs float-level as labelled |
| 3 | Half-mixing kill + repair (025 §2, 026 C1/C2) | **CONFIRMED** (independent exact-Fraction re-implementation, ≤ 1e-11 agreement with both records) |
| 4 | Block kills, ψ threshold, comonotone identity, rescue (025 §3/B3, 026 C3/C4) | **CONFIRMED** (ψ by exact algebra; Gain/CR/rescue re-derived independently) |
| 5a | 026 interval helpers + kit usage + scoping of 009(i)/EM | **CONFIRMED** |
| 5b | 026 "certified under either kit's audit" via intersection | **REFUTED as an argument** — intersection is sound only under both audits; conclusion re-established here by single-kit certification (each kit alone passes 12/12) |
| 6 | 027 survival scan + CR/H collapse + latent-mixing | **CONFIRMED**, with the **range corrected to n ≤ 6** (claimed n ≤ 8 refinement never ran; my own n = 8 spot check stays positive) |
| 7 | 028 battery + recipe statement | **CONFIRMED** (6/19 rows independently re-derived incl. both G-gen floors; specializations re-derived); cr_chain cross-checks were transcript-only — now backed by committed code (this record) |
| 8 | Reproducibility of all six tools | **CONFIRMED** (byte-identical; skeptics exit 0) |

**Net assessment:** the 025–028 batch survives adversarial review with
reporting-level corrections only. The kills are real and certified, the
EM lemmas are correct, the recipe-v2 statement matches what its battery
tests, and no numerical claim moved. The pattern to keep: 026's
intersection idea is good *hardening* but was sold as the wrong
*soundness* argument — future dual-kit certifications should require
each kit to certify alone. The 025 upgrade condition ("until a
fresh-session pass re-derives EM/EM-coverage, treat them as
candidate-VERIFIED") is discharged at the level of independence stated
in the front matter; the human owner should confirm that level counts
before quoting the lemmas as fully VERIFIED under the 024 standard.

## Residual risk

- **Independence caveat, restated:** this reviewer had no shared
  conversation state with the author sessions but was spawned from the
  author's session with a task brief, same model family. A
  cross-family or separately-initiated pass would be strictly stronger;
  every artifact needed to run one is committed.
- **Shared-reduction risk:** my chain-rule evaluator shares the
  prefix-grouping *reduction* with both prior evaluators (as any
  chain-rule accounting must). Mitigation: the closed-form/profile
  paths and the entropy-difference identities (which bypass the
  chain rule) agree with the full-history values wherever both exist,
  and the certified instances pin the same numbers in exact arithmetic.
- 025's grids and non-certified families remain float-level EVIDENCE at
  the six-instance certified core, unchanged by this review.
- My kit-soundness spot tests are sanity checks, not audits; the kits'
  standing rests on their 022/024 (A) and 016/017 (B) audits, which
  this review did not repeat line-by-line.

## References

- This repo: 025, 026, 027, 028 (the records under review); 004 (review
  shape, half-mixing R1), 009/011 (CR frame, fact (i), part-D rule),
  022/024 (kit A + audits, fresh-session standard), 016/017 (kit B +
  audit), 023 (floor endpoints).
- Tools: `explore/uc_reviewer029_reimpl.py`,
  `explore/uc_reviewer029_certaudit.py`; outputs
  `data/uc_reviewer029_reimpl_out.txt`,
  `data/uc_reviewer029_certaudit_out.txt`.
- No external sources.
