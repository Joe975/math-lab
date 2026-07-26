# 001 — Map of the entropy method and its (3−√5)/2 barrier

- **Problem:** union-closed sets conjecture (Frankl), `problems/union-closed.md`
- **Date:** 2026-07-25
- **Type:** frontier survey / barrier analysis + computational baseline (no new theorem claimed)
- **Tools:** web survey of the primary literature (arXiv 2211.09055, 2211.11731,
  2211.11689, 2211.11504, 2211.13139, 2212.00658, 2212.12500, 2306.08824,
  2306.12351, 2301.09664, 2302.12276); `tools/uc_search.py` (written and run here)

Notation used throughout: `h(x) = −x log₂ x − (1−x) log₂(1−x)` (binary entropy),
`φ = (√5−1)/2 ≈ 0.6180340`, `ψ = (3−√5)/2 ≈ 0.3819660`. Key identities:
`ψ = 1 − φ`, `φ² = ψ` (equivalently `ψ² − 3ψ + 1 = 0`), and
`2ψ − ψ² = 1 − (1−ψ)² = 1 − ψ`. For a family `F` and element `i`, the
*frequency* of `i` is `Pr_{A~Unif(F)}[i ∈ A]`.

## Approach

Map the entropy-method frontier precisely enough that later cycles can (i)
re-derive the ψ bound from scratch, (ii) know exactly what the Chase–Lovett
approximate-counterexample barrier does and does not block, and (iii) start
from concrete, falsifiable ideas for injecting *exact* union-closure — the
property the current argument provably discards. In parallel, build the
reusable small-ground-set search tool requested in the problem file's attack
surface and establish the computational baseline.

## What was done

- Verified the technical content below against the primary sources (abstracts
  and full texts via ar5iv where available); statements of the key lemmas
  (Gilmer Thm 1 / Conjecture 1, AHS Claim 3, Chase–Lovett Cor 3.2, Claim 4.1,
  Thm 1.3, Example 1.4, Sawin's extension, Pebody's optimality, Yu/Cambie's
  0.38234 evaluation, Liu's 0.38271) are quoted from those sources.
- Independently re-derived the chain-rule argument (section a) and checked the
  two-variable lemma numerically at ~10 points including the equality point.
- Did two small original computations to sharpen the candidate-ideas section:
  the KL-interpolation threshold curve `p(c)` and the k-wise-union critical
  constants (both reproduced below; one-off scripts, results reproducible from
  the formulas given).
- Wrote and ran `tools/uc_search.py` (exhaustive n ≤ 4, heuristic n = 5, 6).

---

## (a) Gilmer's entropy argument (re-derivable form)

**Theorem (Gilmer 2022, sharpened by AHS / Chase–Lovett / Sawin / Pebody
2022).** If `F ⊆ 2^[n]` is union-closed and `F ≠ {∅}`, some `i ∈ [n]` has
frequency ≥ ψ = (3−√5)/2.

**Step 0 (setup).** WLOG `|F| ≥ 2` (if `F = {S}`, `S ≠ ∅`, any `i ∈ S` has
frequency 1). Let `A, B` be i.i.d. *uniform* on `F` and `U := A ∪ B`. Union
closure enters the argument at exactly one point:

    U ∈ F  always   ⟹   H(U) ≤ log|F| = H(A).                    (∗)

Goal: show that if every frequency is < ψ then `H(U) > H(A)`, contradicting (∗).

**Step 1 (chain rule + dropping conditioning).** Order coordinates `1..n`;
write `A_{<i}` for the first `i−1` indicator bits of `A`. Since `U_{<i}` is a
function of `(A_{<i}, B_{<i})`:

    H(U) = Σᵢ H(Uᵢ | U_{<i}) ≥ Σᵢ H(Uᵢ | A_{<i}, B_{<i}).

Conditioned on `A_{<i} = a, B_{<i} = b`, the bits `Aᵢ, Bᵢ` are independent
(A, B independent). Writing `x = Pr[Aᵢ = 0 | A_{<i} = a]`,
`y = Pr[Bᵢ = 0 | B_{<i} = b]`, we get `Pr[Uᵢ = 0 | a,b] = xy`, so

    H(Uᵢ | A_{<i} = a, B_{<i} = b) = h(xy).

**Step 2 (the key entropy inequality).** The sharp form (Chase–Lovett
Corollary 3.2; equivalent to AHS's Claim 3 after their reduction of the
measure-optimization to two-point measures):

    For all x, y ∈ [0,1]:   h(xy) ≥ (1/(2φ)) · (x·h(y) + y·h(x)),      (KEY)

with `1/(2φ) = (√5+1)/4 ≈ 0.809` and equality iff `x = y = φ` (or degenerate
endpoints). The one-variable diagonal form proved by AHS is: for `x ∈ [φ,1]`,
`φ·h(x²) ≥ x·h(x)`, equality only at `x ∈ {φ, 1}`. AHS verified it by
interval arithmetic on three subdivided regions; Boppana (arXiv:2301.09664)
later gave a computer-free proof by elementary calculus. The inequality is
*delicate*: on the diagonal it is within ~0.3% of equality across much of
`[0.4, 1]` (spot-checked numerically here), which is why crude estimates lose
so much (Gilmer's original constants) and why there is no slack to spend.

**Step 3 (take expectations, sum over coordinates).** Averaging (KEY) over
`(a,b)` with the product measure, and using
`E_a[x_i(a)] = Pr[Aᵢ = 0] = 1 − freq(i)` and
`E_a[h(x_i(a))] = H(Aᵢ | A_{<i})`:

    H(Uᵢ | A_{<i}, B_{<i}) ≥ (1/(2φ)) [ (1−freq(i))·H(Bᵢ|B_{<i}) + (1−freq(i))·H(Aᵢ|A_{<i}) ].

(This is Chase–Lovett Claim 4.1: if `Pr[Aᵢ=0], Pr[Bᵢ=0] ≥ p` for all `i`, then
`H(A∪B) ≥ (p/(2φ))(H(A) + H(B))`.) If every frequency ≤ q, summing over `i`:

    H(U) ≥ ((1−q)/(2φ)) (H(A) + H(B)) = ((1−q)/φ) H(A).

**Step 4 (contradiction).** `(1−q)/φ > 1 ⟺ q < 1 − φ = ψ`. So if all
frequencies were < ψ, then `H(U) > H(A)` (strictness propagates since
`H(A) = log|F| > 0` forces some coordinate to carry positive conditional
entropy), contradicting (∗). ∎

**Remarks for a re-deriver.**
- Gilmer's original paper proves the weaker Theorem 1: frequencies ≤ 0.01 ⟹
  `H(A∪B) ≥ 1.26·H(A)`, via case analysis (`h(p+p′−pp′) ≥ 1.4·(h(p)+h(p′))/2`
  for `p,p′ ≤ 0.1`; concavity bound `h(p+p′−pp′) ≥ (1−p)h(p′)`; Markov to
  control coordinates with large conditional bias). The follow-ups replace
  this with (KEY).
- The argument actually proves a purely *distributional* statement: for ANY
  distribution on `2^[n]` with `H(A) > 0` and all marginals `Pr[i ∈ A] < ψ`,
  i.i.d. copies satisfy `H(A∪B) > H(A)`. Union-closed families only enter
  via (∗) with the uniform distribution. This generality is exactly what the
  barrier (c) exploits.

## (b) The improvement chain and what each step added

| step | ref | constant | new ingredient |
|---|---|---|---|
| Knill 1994; Wójcik | — | Ω(1/log|F|) | pre-entropy era best |
| Gilmer, arXiv:2211.09055 (16 Nov 2022) | — | **0.01** | the whole framework: i.i.d. samples, `H(A∪B) ≤ log|F|`, coordinate chain rule, one-dim entropy estimate (crude). Identified the Bernoulli-product obstruction at ψ; posed the sharp inequality as a conjecture; posed "Conjecture 1" (see (c)) aiming at 1/2 |
| Alweiss–Huang–Sellke, 2211.11731 (≈1 wk later) | | **ψ = (3−√5)/2 ≈ 0.381966** | proved the sharp one-variable inequality `φh(x²) ≥ xh(x)` on `[φ,1]` (interval arithmetic) + reduction of the measure optimization to ≤2-point support |
| Chase–Lovett, 2211.11689 (same week) | | **ψ** | two-variable form (KEY); extension to *(1−ε)-approximately* union-closed families with loss `δ = 2ε(1 + log(1/ε)/log|F|)`; the tightness construction (the barrier, see (c)) |
| Sawin, 2211.11504 (same week) | | **ψ**, sketch of ψ+ε | sharp estimate independently; **disproved Gilmer's Conjecture 1**; observed that mixing the i.i.d. coupling with a "max-entropy" coupling (distributions supported on F but non-uniform / non-product structure) breaks the tightness of the Bernoulli(ψ) obstruction ⟹ strictly beyond ψ |
| Pebody, 2211.13139 (same week) | | **ψ**, proved optimal | computed the exact optimum of the base (single-measure, i.i.d.) Gilmer scheme: it is exactly ψ — the base method is closed out |
| Yu 2212.00658 (Entropy 2023); Cambie 2212.12500 | | **≈ 0.38234** | independently evaluated Sawin's mixture scheme explicitly; this is the *optimum* of the iid+max-entropy convex-combination class |
| Liu, 2306.08824 (ISIT 2024) | | **≈ 0.38271** | "conditionally i.i.d." couplings: A, B i.i.d. given an auxiliary variable; improvement partially rigorous, final constant rests on a numerically solved 9-dimensional optimization |
| status mid-2026 | | ~0.3827 | no rigorous movement of order beyond the third decimal found in this survey; 1/2 wide open. Survey: Cambie 2306.12351 |

Reading of the chain: everything after Gilmer is (i) sharpening his single
estimate to its true optimum (ψ, provably optimal per Pebody), then (ii)
enlarging the class of sampling distributions (Sawin → Yu/Cambie → Liu) for
diminishing returns: +0.0004, then +0.0004 again. Each enlarged class gets its
own computed ceiling far below 1/2.

## (c) The barrier: approximate union-closure is all the argument sees

**Definition (Chase–Lovett).** `F` is *c-approximately union-closed* if for at
least a `c` fraction of (ordered, uniform) pairs `(A,B) ∈ F²`, `A∪B ∈ F`.

**Positive part (CL Theorem 1.3).** Every `(1−ε)`-approximately union-closed
family has an element of frequency ≥ ψ − δ, `δ = 2ε(1 + log(1/ε)/log|F|)`.
(Proof: condition on the indicator `I = [A∪B ∈ F]`; `H(A∪B) ≤ H(I) +
Pr[I=1]·log|F| + Pr[I=0]·2log|F| ≤ log|F| + small`, then run (a).)

**Negative part — the blocking family (CL Example 1.4).** Let

    F₁ = { S ⊆ [n] : |S| = ψn + n^{2/3} }        (a single slice, slightly above ψn)
    F₂ = { S ⊆ [n] : |S| ≥ (1−ψ)n }              (the whole top of the cube)
    F  = F₁ ∪ F₂.

Mechanics (re-derived here, all verified against the identity `2ψ−ψ² = 1−ψ`):

1. `|F₂| / |F₁| ≈ φ^{−n^{2/3}} → 0`: each step from weight `ψn` toward the
   middle multiplies a binomial coefficient by ≈ `(1−ψ)/ψ = 1/φ ≈ 1.618`, so
   the `n^{2/3}` bump makes the slice exponentially-in-`n^{2/3}` larger than
   the entire top block. Hence `|F| = (1+o(1))|F₁|` and every element has
   frequency ≤ ψ + O(n^{−1/3}) = ψ + o(1) (the slice is symmetric).
2. *Approximately* union-closed: two uniform sets from the slice have
   per-element density `p = ψ + n^{−1/3}`; their union has density
   `1−(1−p)² = (1−ψ) + 2(1−ψ)n^{−1/3} + O(n^{−2/3})`, i.e. weight
   `(1−ψ)n + Θ(n^{2/3})` with fluctuation only `O(√n)`. So whp the union lands
   in `F₂`. Unions involving an `F₂` member are always in `F₂`. Total failure
   probability o(1). This is exactly why the bump `n^{2/3}` is chosen: large
   enough to push typical unions over the `(1−ψ)n` threshold, small enough to
   keep frequencies at ψ + o(1).
3. *Exactly* union-closed it is emphatically not: any two slice sets with
   large overlap have a union of weight in the deleted middle band
   `(ψn + n^{2/3}, (1−ψ)n)` — not in `F`. These are an o(1) fraction of
   pairs, so the entropy argument never sees them.

**Barrier statement (precise).** The Gilmer argument interfaces with
union-closure through exactly one fact: `Pr_{A,B iid Unif(F)}[A∪B ∈ F] = 1`
(used only to get `H(A∪B) ≲ log|F|`), and this fact degrades gracefully to
`1−ε`. Therefore *any* proof whose only use of closure is
measure-theoretic/average-case — i.e. any proof that remains valid for
(1−o(1))-approximately union-closed families, equivalently any proof
insensitive to destroying closure on an o(1) fraction of pairs — cannot
establish a constant above ψ, because CL's family is (1−o(1))-approximately
union-closed with max frequency ψ + o(1). ψ is not a weakness of the
estimates (they are sharp, per AHS/Pebody); it is the true answer to the
question the method asks.

At the distributional level the same obstruction is the i.i.d.
Bernoulli(ψ)-product measure: per coordinate, `h(2ψ−ψ²) = h(1−ψ) = h(ψ)`, so
`H(A∪B) = H(A)` exactly — the distributional theorem in (a) is tight at ψ, and
CL's family upgrades this from "the inequality can't be improved" to "no
average-case argument can be improved."

**Exactly which properties of exact closure the argument fails to use:**

1. **Worst-case (universal) pairwise closure.** The atypical pairs — highly
   overlapping slice sets, whose unions populate the deleted middle band —
   carry zero weight in the entropy bookkeeping. Exact closure would force
   the entire middle band's worth of unions into `F`, which (see Idea A) would
   blow up `|F|` and destroy the frequency count of the CL family.
2. **Closure under arbitrary re-weighting.** For exact `F`, `A∪B ∈ F` for
   *any* distribution supported on `F`, not just uniform — non-uniform or
   dependent samplers are legal. This is precisely the crack Sawin/Liu opened
   (uniform-pairs approximate closure does NOT license non-uniform samplers,
   which is why they could pass ψ at all), but so far only with
   family-oblivious samplers, for +0.0007 total.
3. **Iterated (k-wise) closure.** `A₁∪…∪A_k ∈ F` for all k. Caution — used
   naively this is *worse*: the k-wise Bernoulli obstruction sits at the root
   of `(1−p)^k = p`, computed here: p₂ = 0.381966 (=ψ), p₃ = 0.317672,
   p₄ = 0.275508, p₅ = 0.245122 — strictly decreasing in k (cf. the "almost
   k-union closed" line, arXiv:2302.12276). Any use of iteration must couple
   the constraints for all k jointly on the same F, not apply them separately.
4. **The exact identity `D(A∪B ‖ Unif(F)) = log|F| − H(A∪B)`.** For exact
   `F` and uniform `A`, this holds *identically* (A∪B is supported inside F).
   For the CL family it fails catastrophically: `A∪B` concentrates on `F₂`,
   whose relative size is `φ^{−n^{2/3}}`, so the divergence from uniform-on-F
   explodes. Gilmer's refuted Conjecture 1 lived exactly here: he conjectured
   that all marginals < 1/2 and `H(A) > 0` imply
   `H(A∪B) + D(A∪B‖A) > H(A)`; applied to uniform-on-F the left side *equals*
   `log|F| = H(A)` identically, so the strict inequality would have forbidden
   any union-closed family with all frequencies < 1/2 — the full conjecture.
   Sawin disproved the distributional statement (with a non-product
   distribution; note it is *exactly tight* for products at p = 1/2, since
   `d(q‖1/2) = 1 − h(q)` makes `h(q) + d(q‖1/2) ≡ 1`). The refutation kills
   the distributional route as stated, not the family-level identity.
5. **Semilattice structure.** `F` is a join-semilattice: it has a maximum
   element `∪F`; for each `S ∈ F` the map `T ↦ T∪S` is a retraction of `F`
   onto the up-set `up_F(S) = {T ∈ F : T ⊇ S}`. None of the order/lattice
   structure (Möbius counting, up-set sizes, chains) is used anywhere in the
   entropy chain.

## (d) Candidate ideas for injecting exact closure (all SPECULATION unless marked)

**Idea A — stability + slice endgame (dichotomy).** The equality analysis of
(KEY) is rigid: near-equality forces conditional marginals `x_i(a) ≈ φ`
throughout, i.e. `A` is close (in the chain-rule/KL sense) to the
Bernoulli(ψ) product — weight concentrated near ψn, `log|F| ≈ n·h(ψ)`,
"slice-like". Program: (i) prove a *stability version* of the ψ theorem —
either some frequency ≥ ψ + δ, or `F` is ε(δ)-close to the slice profile;
(ii) kill the slice profile *using exactness*: an exactly union-closed family
essentially concentrated on a weight band around ψn must contain unions of its
own overlapping pairs, which sweep out a positive-density middle band above
the slice, contradicting the entropy profile `log|F| ≈ n·h(ψ)` (the middle
band alone has more sets than that). First falsifiable sub-question: what is
the union-closure of the full ψn-slice, and does its max frequency return to
≥ 1/2? (Small-n versions computable with `tools/uc_search.py`-style code.)
SPECULATION at step (ii); step (i) stability statements are a known genre and
plausibly extractable from the AHS two-point-measure reduction.

**Idea B — weighted-KL interpolation (a quantitative ladder between ψ and
1/2).** For exact `F` and uniform `A`: `H(U) + c·D(U‖A) = log|F| − (1−c)D ≤
H(A)` for every `c ∈ [0,1]`. So any theorem of the form
"all marginals < p(c) and H(A) > 0 ⟹ H(U) + c·D(U‖A) > H(A)"
immediately gives the bound p(c) for exact families. Computed here, the
Bernoulli-product obstruction for coefficient `c` sits at the first root of
`h(2p−p²) + c·d(2p−p²‖p) = h(p)`:

    c    : 0       0.1      0.2      0.3      0.5      0.7      0.9      1.0
    p(c) : 0.38197 0.39264  0.40356  0.41473  0.43782  0.46192  0.48705  ≥ 0.5

Gilmer's c = 1 endpoint is refuted (Sawin), but the refuting example is
non-product and its *c-profile is unknown*: determine the largest c₀ such
that the c-weighted statement survives all known counterexamples, then try to
prove it for some c > 0. Even c = 0.1 would beat every known bound
(0.3926 > 0.38271). Key technical question: does the chain-rule machinery
tensorize for the functional `H + c·D`? (D against the uniform-on-F reference
is not coordinate-decomposable in general — that is where the difficulty
migrates.) SPECULATION, but sharply posed and each increment is checkable.

**Idea C — family-adaptive / dependent couplings on atypical pairs.** Exact
closure legalizes samplers that approximate closure does not: any *joint*
distribution of `(A,B)` with both marginals supported in `F` still has
`U = A∪B ∈ F` and `H(U) ≤ log|F|`. In the CL family the deleted middle band
is exactly the image of the *overlapping* pairs; a coupling biased toward
pairs with `|A∩B|` above typical (e.g. `B` re-sampled conditioned on
`|A∩B| ≥ θ|A|`) would detect its non-closure instantly, while for exact
families it yields new legal entropy constraints. The obstacle: Step 1 of (a)
uses independence of `A,B` to factor `Pr[Uᵢ=0|a,b] = xy`; for dependent
couplings one needs a conditional-independence surrogate — Liu's
conditionally-i.i.d. auxiliary-variable coupling (0.38271) is exactly the
first rung of this ladder, still family-oblivious. The proposal is to make
the auxiliary variable *depend on F* (e.g. on the weight distribution or on
`A` itself). SPECULATION; concrete first target: design any coupling class
whose Bernoulli-product obstruction exceeds 0.383.

**Idea D — up-set retraction constraints (local entropy inequalities).** For
every fixed `S ∈ F`, exactness gives `A∪S ∈ up_F(S)`, hence
`H(A∪S) ≤ log|up_F(S)|` — one constraint *per set*, and after averaging over
`S ~ B`: `H(A∪B | B) ≤ E_B[log|up_F(B)|]`. This is strictly finer than the
single global constraint (∗). Honest negative check (computed here): the CL
family satisfies these local constraints *with slack* — for a slice set `S`,
`up_F(S)` consists of the ≥(1−ψ)n-weight supersets of `S`, and since
`(1−2ψ)/(1−ψ) = ψ < 1/2`, `log|up_F(S)| ≈ (1−ψ)n` bits, comfortably above
`H(A∪S|S) ≈ (1−ψ)h(ψ)n ≈ 0.593n`. So up-set constraints alone do NOT kill
the barrier example; to bite they must be combined with a counting tradeoff,
e.g. `Σ_{S∈F} |up_F(S)| = Σ_{T∈F} |down_F(T)|` (which ties up-set sizes to
element frequencies) or with weight-conditioning. Recorded partly as a
warning: the naive version of the most obvious "use the lattice" idea is
already checked and insufficient. SPECULATION beyond the checked part.

Known dead ends to not re-tread: (1) sharpening (KEY) further — impossible,
it is tight at `x=y=φ` and Pebody closed out the whole base scheme at ψ;
(2) k-wise unions used separately — constants strictly worsen (table in
(c) item 3); (3) Gilmer's Conjecture 1 as stated — refuted by Sawin.

---

## Computational search (`tools/uc_search.py`)

Re-runnable tool, deterministic per seed. Representation: sets as bitmasks,
family as a set of masks. Objective: `max_i freq(i)/|F|` (exact `Fraction`),
minimized over union-closed `F ∉ {∅-family, {∅}}`.

- **n ≤ 4: exhaustive** over all `2^(2^n)` subfamilies, filtered for closure.
  Counts found: 2 / 12 / 120 / 4958 union-closed families for n = 1..4
  (n = 2 count verified by hand), of which 1 / 4 / 14 / 51 achieve exactly 1/2.
- **n = 5, 6: heuristic** — deterministic sweep of all ≤3-generator closures,
  structured seeds (power sets, `{∅,{1}}`, closure of the 2-slice), and
  hill-climbing local search (moves: delete a set not expressible as a union
  of two remaining members — checked before deletion; add a random set and
  re-close), 60–150 restarts.
- **Verification:** every reported optimum is independently triple-checked
  (pairwise closure re-check, closure-fixed-point re-check, per-element exact
  rational recount). Anything below 1/2 would be flagged as a potential
  counterexample and explicitly marked untrusted pending further independent
  verification. Nothing below 1/2 appeared.

**Result** (seeds 0 and 7, restarts 60/150, iters 400/800 — identical):

    n = 1..4 (exhaustive):  minimum = 1/2 exactly   (witness {∅,{1}}; power sets also achieve it)
    n = 5..6 (heuristic):   best found = 1/2        (no family below 1/2)

Consistent with the conjecture; the minimum 1/2 is provably exact for n ≤ 4
(exhaustive), evidence-only for n = 5, 6. Runtime ≈ 2 s at defaults.

## Outcome

Dead end *as an attempt to beat ψ* — deliberately: this cycle's goal was the
map, not a new constant. Deliverables: the barrier statement made precise
(any average-case use of closure stops at ψ; CL's slice+top family is the
universal blocker), a re-derivable writeup of the sharp argument, four
concrete injection ideas with two small original computations backing them
(the `p(c)` ladder; the k-wise decay table), one candidate idea (up-sets,
Idea D) already partially falsified honestly, and a verified computational
baseline plus reusable tool.

## Why it failed / what survived

- **Why the entropy method fails at ψ:** not loose estimates — the key
  inequality is tight and the base scheme's optimum is exactly ψ (Pebody).
  The method's only interface to closure is average-case (`H(A∪B) ≤ log|F|`
  under uniform i.i.d. sampling), and at that interface ψ is the *correct*
  answer (Chase–Lovett's approximately-union-closed slice+top family achieves
  ψ + o(1)). Passing ψ requires using worst-case/exact closure; the only
  moves that have done so (Sawin, Liu: non-uniform / conditionally-i.i.d.
  samplers, legal only under exact closure) bought ≈ 0.0007 total, each with
  its own computed ceiling.
- **Survived:** the (a)–(c) map (verified against sources); the p(c)
  interpolation table (product obstructions only — a ceiling, not a theorem);
  the k-wise decay computation; the negative check that naive up-set
  constraints don't kill the CL example; `tools/uc_search.py` + the exhaustive
  n ≤ 4 ground truth (min = 1/2, and the extremal counts 1/4/14/51).

## Leads generated

1. **(best)** Idea B: determine the c-profile of Sawin's counterexample to
   Gilmer's Conjecture 1; if some c₀ > 0 survives, attempt the c-weighted
   entropy argument. Even c = 0.1 ⟹ 0.3926, a new record.
2. Idea A first step: compute union-closures of slice families for small n
   (extend `uc_search.py`) — test "closure of a slice restores max frequency
   ≥ 1/2"; then look for a stability version of the AHS reduction.
3. Idea C first target: any family-adaptive coupling whose product
   obstruction beats 0.383.
4. Idea D refinement: combine up-set constraints with the counting identity
   `Σ|up_F(S)| = Σ|down_F(T)|` before discarding the lattice route.
5. Characterize the 51 extremal (exactly-1/2) families at n = 4 from the
   exhaustive run — structural data for the "equality case" of the
   conjecture itself (attack-surface item in the problem file).

## References

- Gilmer, *A constant lower bound for the union-closed sets conjecture*, arXiv:2211.09055
- Alweiss, Huang, Sellke, arXiv:2211.11731 — ψ via sharp one-variable inequality
- Chase, Lovett, *Approximate union closed conjecture*, arXiv:2211.11689 — ψ, approximate version, tightness (Example 1.4)
- Sawin, arXiv:2211.11504 — ψ, disproof of Gilmer's Conjecture 1, beyond-ψ sketch
- Pebody, *Extension of a Method of Gilmer*, arXiv:2211.13139 — optimality of ψ for the base scheme
- Yu, *Dimension-free bounds…*, arXiv:2212.00658 / Entropy 2023 — 0.38234 explicit
- Cambie, arXiv:2212.12500 — 0.38234 explicit; survey arXiv:2306.12351
- Liu, *…via conditionally IID coupling*, arXiv:2306.08824 / ISIT 2024 — ≈ 0.38271
- Boppana, arXiv:2301.09664 — computer-free proof of the key inequality
- *Almost k-union closed set systems*, arXiv:2302.12276
