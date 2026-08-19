# Status

The live ledger: where every problem stands, what is queued next, and what has
already been ruled out. Read this before starting work.

## TL;DR (updated 2026-08-19)

**2026-08-19 HEADLINE (union-closed 020/021): GAP 1 CLOSES — every
odds-ratio-control candidate is REFUTED, certified, and the whole
per-history-weighting family is ruled out.** A dual-family LLM swarm
(gpt-5.6-luna + gemini-3.7-flash, 26 briefs, ≈ $0.21, protocol
`docs/SWARM.md`) drafted attack directions; two cross-family convergences
supplied the designs, and every kill was implemented and certified by the
director in exact rational arithmetic. (i) The unsigned |σ|-weighted
control — the last margin-modulated reading — fails at n = 4 (7 atoms)
and n = 5 at healthy margins inside the λ-window: its first-order
coefficient is a sensitivity-kernel-weighted Gram form, the kernel is NOT
PSD in-regime (diagonal collapses at the z = 1/2 kink), and realizing the
violating margin pair with anti-aligned future biases gives certified
MM_abs < 0 with the plain aggregate certified positive on the same
tables. (ii) The λ-window-restricted aggregated control dies at n = 4:
prefix-signed cancellation (Σ_a ν(a)D_a = 0 at every coordinate, far from
products) admits measures with certified-negative λ² coefficient —
A < 0 certified at three rational tilts λ ≈ 0.07–0.11 ≪ λ_win(4),
marginals ≤ 0.2575, tensor-extending float-negative in-window at
n = 8, 12. (iii) The assembly's true aggregate requirement (the secant
condition, identified independently by both swarm families) is certified
violated inside the window, and R₊-pure certified witnesses rule out
EVERY per-history weight compatible with the chain rule. With
005/006 + 007/013 + 016/017 + 018/019, every stated form of Gap 1 is now
dead; the route stays LIVE with Gap 2 (mutual-information tax) as the
surviving bridge and a precise no-go boundary for what any new bridge
must avoid. **Skeptic pass 022 ran the same day (cross-family
drafted-verifier pattern): a gemini-3.7-flash-drafted exact verifier —
own digit-extraction log₂, own bisection, own census, zero shared code,
audited then run unmodified — reproduces all 12 certificates
digit-for-digit including the three R₊-purity legs; zero corrections.
The window kill is additionally certified at n = 8 (225-atom tensor,
t = 25/24, λ = 0.059 ≪ window 0.969, A certified in
[−2.5616e-7, −2.5616e-7]). And Gap 2's CR > 0 candidate SURVIVES: +0.40/
+0.80 on the kill geometries (022), then a seeded adversarial pass (023)
grinds it only to floors +0.074…+0.22 with no in-regime violation — the
sole "negative" is the degenerate H(μ) → 0 corner where CR → 0 correctly,
so any quantitative TAX statement must scale with H(μ). The signed/
λ-integrated per-history form is certified negative at two rational
tilts per witness (023), closing 020's last float-level piece. Gap 2 is
the route's sole surviving bridge and has earned H-scaled proof effort.
~~Residual for a fresh session: reviewer-level independence (022/023 were
same-session; scope caveats are in the records).~~ **DISCHARGED in 024
(same day, fresh session): all 12 certificates reproduce byte-for-byte
AND inside a third from-spec-only reimplementation; the verifier audit
and both hand re-derivations check out; Gap 2 survives three fresh
adversarial trajectories (floors +0.0697…+0.30, zero in-regime
violations). One reporting-level correction: 023's CR engine is NOT
deterministic (string-hash job seeding) — its checkpoint is not
regenerable and the "+0.0742 sharpest floor" is trajectory-specific.**
**Same-day follow-up (025): the recipe's two non-Sinkhorn branches — which
no adversary had ever exercised (023 lead 2) — both FAIL as stated,
in-regime: half-mixing (q = P₁/2) goes CR < 0 at every n on a p_h → 1
sliver (margins ~1.5e-2, two independent implementations to 4e-15), and
the iid-given-k block coupling is Θ(n)-negative on any mixture with a
component above ψ = (3−√5)/2 — the AHS constant resurfacing as a
component-wise coupling threshold (the adaptive variant instead flatlines
at CR = 0 exactly on all-above-half mixtures, failing strictness). One
repair covers both: the extended mixing lemma (ST = 0 for arbitrary
ε-joints AND arbitrary s-laws) plus concavity of w ↦ H(wδ∅ ⊕ (1−w)ν)
gives sup_q CR > 0 on the entire in-regime {0}∪[1/2,1) genre — proved
except two float-level 1-D boundary sweeps (n = 1, 2) — and ∅-mixing also
rescues every constructed block failure (+0.69 where Sinkhorn manages
+0.22). So Gap 2's candidate survives its third adversarial pass, but Gap
3's recipe draft is refuted-and-corrected; the branch-assignment rule
("never iid-couple above-ψ mass"; the surprisal criterion for ∅-mixing)
is now sharp. Floor pre-work for the proof effort: 023's in-regime
endpoints have CR/H ∈ [0.0387, 0.0834], and the anneal was grinding CR,
not H — so any CR ≥ c·H target has c ≤ 0.0387. Lemmas EM/EM-coverage are
same-session skeptic-passed; reviewer pass pending (024 standard).
Same-window follow-up (026): the kills and repairs are now CERTIFIED in
exact arithmetic at six rational instances (half-mixing CR < 0 at
n = 2/6/16 with H(μ) > 0.968 certified; block Gain < 0 at n = 8; both
rescues positive), every log₂ evaluated by BOTH audited kits (022
digit-extraction ∩ 016 atanh-series, 98 calls, zero disagreements,
widths ≤ 2.5e-32) — the dual-kit intersection is now the recommended
default for certifications.**
Also: the deep-anneal
float negatives at near-empty supports are certified ARTIFACTS (exact
recomputation flips −24.5 to +23.4) — 007's underflow lesson confirmed a
second time. A parallel dual-provider ideation sweep (021, MAP) sharpened
the queued transfer-operator lead with an exact eigenfunction identity
(coordinate frequencies ARE eigenvalues of T_μ), settled number theory
no-purchase by cross-family agreement, and killed two idea-shapes inline
with an explicit 5-member family.

## Previous TL;DR (2026-08-16)

**Open for contributions. No automated loop is currently running** — the hourly
cycle that produced these records was stopped once its in-flight work closed
out, and no cron trigger is active. Nothing is half-finished: every line is
written up, and the queue below is a list of starting points rather than
abandoned work.

**2026-08-16 HEADLINE (billiards 016): the seven S-lemmas of 015 §5 are
PROVEN — the window law holds for ALL (a,b), and with it the coverage
conjecture for the W family closes: every obtuse angle has an explicit
alive W(a,b) triangle.** The queued skeptic pass on 015 ran with the
mandated refute stance and did the queue's main task in the same record:
all seven fan-margin lemmas plus the τ row re-derived from scratch (the
swarm drafts are not in the repo and were never read — only the committed
briefs), each at k = 1/2 for all a ≥ b ≥ 1, a ≥ 2, each ≤ half a page:
D-AB/B-AB by progression-window bounds (sin L_s > sin at), D-BA/D-BC by
one cosine comparison each, D-AC by an exact product-to-sum pairing whose
only sign-varying factor φ₁ = c − A/2 + (b+3−a)t is nonnegative for every
member except W(2,2) (closed by a two-line domination), B-AC/B-BC jointly
via one new "birth core" inequality sin α·sin(H−(a−b)t) >
sin θ·cos(at)·sin((b+1)t) (chord + sinc floors, one b = a case split),
and τ via a new closed form R = −4·sinθ·sin(aα)·sin((b+1)β) — which
makes τ ≠ 0 AUTOMATIC at every alive point of every member, subsuming
the certifier's τ legs. T2 is confirmed with a one-line simplification
(n = a gives sin(z_a) = −sin((b+1)β) > 0, killing the aα ∈ (90,120)
branch instantly); N1 is confirmed by a new derivative-free proof.
Skeptic findings, neither load-bearing: 015 §1's s = 0 C-image identity
is −N1, not −N1−N2; and 015 §4 overclaims — wmargins_certify.py never
evaluates the death-segment BC0 (= N4) row (loop starts at j = 1, no N4
leg), so one margin per member rested on 005's I4 rather than that run —
repaired by an exact n4leg certificate, 104/104 members. Independent
validation: 1200/1200 sign agreements with 002's blind-era corridor
geometry, 62 fresh members ring-verified, all 32 proof-step inequalities
scanned to a = 3000 (zero violations), 453k-point T2 adversarial scan
clean. Remaining under the coverage chain: only the general-(a,b)
endpoint bookkeeping (015's narrow SPECULATION; formal ring pass queued)
and the standard skeptic pass on 016's own proofs — both queue 11's
replacement items.

**2026-08-16 (billiards 015): queue 11 broken open — the birth
law's mechanism is found and its necessity half is PROVEN for all (a,b),
and the full window law is EXACT for every member with a ≤ 14.** Every
gate margin of W(a,b) now has a closed form (one sinusoid around one of
two pivots; ring-verified bijectively). The birth edge is the N2/PB
double wall meeting at (α,β) = (90/a, 90/(b+1)) — 006's SPECULATION
birth law is the second corner of a two-wall system, windows-touch is a
one-line identity, and 013's Diophantine condition is literally
θ ∈ (θ_d, θ_birth). Birth-side necessity (alive ⟹ θ < θ_birth, hence
γ > γ_birth) is proven for ALL a ≥ b ≥ 1, a ≥ 2 by a five-step
sign-constant-progression pigeonhole — no case tree, no intervals
(swarm-drafted, director re-derived with one repair). Sufficiency: two
canonical segments realize every θ in the window, certified in exact
interval arithmetic for all 104 members a ≤ 14 (float-verified to
a = 60, 1,829 members; the death segment's alive range ends at the
sharp wall t* = 45/(a(b+1)), explaining 005's ρ=2 uniformity).
**Coverage of (90°, 180°) is now a theorem for every γ whose 013
witness has a ≤ 14, and for ALL obtuse γ it is conditional on exactly
seven elementary fan-margin lemmas** — each numerically bulletproof,
each with a swarm-drafted proof awaiting director/skeptic re-derivation
(five others, including the hardest two — parametric N1 and the T2
pigeonhole — are already director-verified). Verifying those seven
(≤ 2 pages of elementary trig each) closes the coverage conjecture; it
is the queue's new top item. Swarm: 13 drafts, ≈ $0.16, protocol per
`docs/SWARM.md`.

**2026-08-16, one line (billiards 014): 013's coverage Diophantine lemma
is skeptic-CONFIRMED — the proof, the constructive witness, and the
checker all survive an adversarial pass** (every step re-derived by hand
including the strictness edge cases and a cleaner minimal-a₀ closure of
the union-endpoint worry; a from-scratch implementation written from the
record's formulas alone reproduces every count and diffs
witness-for-witness with zero mismatches on all 27,397 Farey rationals;
novelty label holds after two fresh searches). Two prose errors
corrected, neither load-bearing: the γ = 144° spot-check names W(4,3),
which is not a witness there (the committed checker itself says W(11,2);
W(4,3) is *never* the constructive witness for any t), and Step 2's
parenthetical inclusion I(q−1,q) ⊂ (1/q, 2/q) is false (the interval is
disjoint from it, lower endpoint exactly 2/q). Queue 12's Diophantine
side is now closed for real; the sole blocker for conditional coverage
of (90°, 180°) remains the birth/window theorem (queue 11).

**2026-08-15, one line (collatz 002): 001's one nonstandard lead — the
graph-structure census of truncated Collatz digraphs, queue item 21 — is
CLOSED, REFUTED by a one-paragraph argument the census then confirmed.**
The Collatz map is a function, so `G_B` has out-degree one everywhere and
is a forest plus the `1→4→2→1` triangle; every unlabelled invariant the
lead proposed is therefore fixed before any arithmetic is consulted. The
dominator tree **is** the orbit tree (`idom(v) = T(v)`, zero mismatches to
B = 2^18 under a real Cooper–Harvey–Kennedy computation); the treewidth is
**exactly 2 for every B**, so the lead's "track growth in B" has nothing to
track; and the minimum cut has the closed form
`|S_k| = (2^{k−1} − (−1)^{k−1})/3 = J_{k−1} ~ 2^k/6`, decided by the single
inequality `3r + 1 > 2^{k+1}` — pure magnitude, no dynamics — with B
dropping out of the cut entirely. Keeping the arithmetic labels, as 001
instructed, shows nothing: the cut is an interval of odd integers, uniform
mod 3 by construction (233,017 each at k = 22). Verified three ways beyond
the census: Dinic max-flow (k ≤ 12), and **three independent
gemini-3.7-flash implementations** under different stances (k ≤ 14, all
agreeing) — the first use of the new cross-family skeptic rule, and one
worker returned a *cleaner proof* than the one the census was built on
(crossings above 2^k have in-degree zero outright). The pricing argument
generalizes to any iterated map and is recorded as the mechanism tag
`functional-graph-invariants`. One unexplained by-product survives as the
live lead: the escape fraction `#{n ≤ 2^k : orbit > 2^{k+1}}/2^k` is
non-monotone — 0.609 at k = 12, dropping to 0.379 at k = 13, then flat near
0.404 through k = 22 (brute-force confirmed).

**2026-08-15, one line (collatz 001): first attempt on Collatz — a
swarm-executed ideation sweep (22 external gpt-5.6-luna workers, 11 field
lenses × 2 stances, $0.084; protocol in `docs/SWARM.md`, new this cycle)
finds the lens taxonomy has weak purchase: 17/22 reports collapse onto two
classical families, both now closed** — finite-memory residue-Lyapunov
certificates are killed at every modulus by an exact re-derivation of the
−1-tower obstruction (`barrier_checks.py`, L ≤ 200, classical root Terras
1976), and cycle-equation sieving offers nothing beyond the recorded
cycle-bound limitation. The parity-uniformity entropy barrier is certified
(k ≤ 16) after correcting a false worker claim (bijection holds for the
shortened map only). Best surviving lead, queued: a dominator/treewidth
census of truncated Collatz digraphs — the sweep's one nonstandard framing.

**2026-08-06, one line (attempt 018, skeptic-confirmed 019): Gap 1's
margin-modulated candidate is REFUTED in both signed readings at n = 7 —
by 007's own 10-atom witness, the instance it was invented to survive.**
The h-sensitivity weighting of 007 lead 2 splits into three precise
readings; the assembly-exact secant form is certified negative
(−0.02225465820566 at the exactly-rational tilt t = 181/16, in-regime,
dichotomy exact — the skeptic's from-scratch engine, different census and
different enclosure kit, matches to width 5.9e-30) on the same tables
where the plain aggregate is certified +1.844669. Mechanism: the
sensitivity weight is signed, and in the record regime the surplus
histories' realized both-zero probabilities sit past the entropy peak
(z̃ > 1/2), so the modulation flips the aggregate's biggest positive
terms — it kills the surplus, not the deficit it was designed to
down-weight. What survives of Gap 1: the unsigned |σ|-weighted control
(positive on everything tried, including 016's three certified kill
instances; softest point +0.0101 at a free-support n = 8 endpoint) and
the λ-window-restricted control (positive across the θ-re-optimized
ladder genre to n = 320, with the first-order protection a₁(n)·n GROWING
4.7 → 18.1 and the joint (θ, λ) attack converging to the trivial λ → 0
boundary). Also settled: 016 lead 3 — the raw-weight ladder does NOT
cross through n = 320. Five corrections in 019, two substantive (raw
ladder MM signs flip on the n ≥ 96 dilution branch; a wrong golden-ratio
aside), neither load-bearing.

**2026-08-04, two cycles.** First (Astra follow-through): a
literature-novelty check joined the skeptic pass, the formalization lane
(`docs/FORMALIZE.md`, status `FORMALIZED`) went in above `VERIFIED`, the
ripple scan of the ten Astra results ran with **zero hits** (see
insights), and the pilot certificate landed (billiards L1, 009/010).
Second cycle HEADLINE: **the i-AGGREGATED odds-ratio control — Gap 1's
restatement — is REFUTED at large n** (016, skeptic-confirmed 017): the
10-atom witness replicates into a ladder MU(n,r) whose aggregate falls
like a − b·log n and crosses zero near n ≈ 90, certified in exact
rational arithmetic at n = 96/128/160 (λ = 2, marginals ≤ 0.309), the
skeptic's from-scratch engine matching to 18+ digits. The kill closes
the recorded ∀λ gap; the surviving Gap-1 candidates are the
margin-modulated control and the λ ≲ c/n window-restricted variant (the
violation lives at λ ∈ [2, 2.5], ~40× above the workable window at
n = 96). Probe-before-proof vindicated a second time: certified
positives at n ≤ 32 masked the asymptotic failure — and a **parallel
line reached exactly that ceiling independently** (014, skeptic-confirmed
015): the same control certified positive to n = 32 with the margin
growing, ~1900 instances, zero violations. The two agree everywhere they
overlap; 014's forward-looking call ("proof effort now justified") is the
part 016 overturns, and its equality-set structure (products are the
exact equality set for λ > 0) survives as the reason the ladder had to be
built to see the crossing. Also this cycle: the **entire 005 Laurent
block is now FORMALIZED** (011/012 — closed-form composition, I1–I4,
glide facts, specializations; 23 theorems, zero corrections, statement
review + independent rebuild complete), **crouzeix got its first
attempt** (blind, 001, skeptic-confirmed 002): a certified 286-start
local-maxima census at n = 3, deg ≤ 3 finds no unknown extremal
structure — every problem carried at the time now has an attempt (the two
onboarded below do not yet) — and the
**maxwell-equilibria escalation is DISCHARGED** (002 confirms 001's
certified 24-equilibria witness; see the standing note below).

2026-08-04: two problems onboarded from a conductance-physics scoping pass
(`docs/PLAN-conductance-problems.md`): **almost-mathieu** (Dry Ten Martini at
critical coupling; exact Chambers-polynomial harness over the real cyclotomic
field, smoke-tested through golden-mean convergent 21/34) and
**three-phase-conductivity** (attainability of optimal three-phase 2D
composite bounds; exact laminate algebra with dual-route verification).
Harnesses self-tested and cross-verified; **no attempts yet on either** — both
first attempts should run blind (queue 18–19).

To pick something up, take an item from the attempt queue and follow
`docs/CYCLE.md` — either by hand, one line at a time, or by restarting an
hourly routine with that file as the prompt. Everything needed is in this
directory. New attempts are welcome as pull requests; see `CONTRIBUTING.md`,
and `AGENTS.md` for whether to work blind or informed.

Current standing. HEADLINE (2026-08-04, second cycle): the
maxwell-equilibria escalation is DISCHARGED — skeptic review (002)
confirms 001's certified 24-equilibria witness: a shifted+enlarged-region
re-census whose tree shares no box with 001's again certifies exactly 24
(and the independent verifier PASSes it), an 80-bit re-run reproduces
every leaf verdict, and the fold brackets 001 predicted are now certified
counts — 12 at q = 4360/10⁶, 16 at q = 4400/10⁶. The centroid Hessian
closes exactly: signature decided by 6561t⁶ vs 128q², degeneracy at
q\* = 81√2·t³/16 ≈ 4.3968·10⁻³, certified index flip across it. Real
correction: verify_equilibria.py's tiling check (prefix-freeness + Kraft)
is UNSOUND for axis-labelled paths — coverage independence silently
rested on driver honesty; closed for 001's certificate by an
axis-consistency reconstruction (0 violations / 235,993 paths), hardening
queued. The explicit refutation witness is now settled library fact.
Union-closed: two independent lines ran the mandated large-n probe of the
i-AGGREGATED odds-ratio control in parallel, and together they settle it.
014 (skeptic-confirmed 015) certified the control POSITIVE in exact
rationals to n = 32 (15 certificates, softest +5.99e-7; independent
re-certification of 5, zero disagreement), float trends to n = 64,
minimum margin GROWING in n across ~1900 instances with zero in-regime
violations, and identified products as the exact equality set (λ > 0,
inclusion proven; converse open). 016 (skeptic-confirmed 017) pushed the
same object past that ceiling on a family 014 did not build — the 10-atom
witness replicated into the unit ladder MU(n,r) with re-optimized shared
weights — and certified it NEGATIVE at n = 96/128/160 (λ = 2, marginals
≤ 0.309). The two searches agree at n ≤ 32 (016 finds nothing negative
there either), so the merged reading is: the ∀λ aggregated control is
REFUTED, 014's positivity is a small-n fact rather than a licence for
proof effort, and what survives is the margin-modulated control and the
λ ≲ c/n window-restricted variant — 016's violation sits at λ ∈ [2, 2.5],
~40× above the workable window at n = 96, which 014's equality-set
structure says is exactly where the boundary question still lives.
Previous headlines: 2026-08-04 first cycle — maxwell-equilibria onboarded,
001 certified the explicit rational five-charge witness (24 > 16) for the
six-day-old arXiv:2607.27197 refutation; 2026-07-31 — the two-track
billiards plan, all four records skeptic-confirmed (six corrections
total, none load-bearing).
Conservative track (005+008): the W(a,b) death law is closed on both
sides — the half word composes in closed form to
R₀·Rot_A(2aα)·Rot_B(−2bβ), so I1–I4 and the glide facts are proven
for ALL (a,b) by one formal-ring check, Lemmas C and D fall to an
elementary monotonicity lemma, the case tree closes with the
a ≤ 2b+3 restriction REMOVED (necessity is now fully parametric for
a ≥ b ≥ 1, a ≥ 2), and certified alive segments into the death
corner give death(W(a,b)) = γ_d(a,b) EXACTLY (sup not attained) for
20 members, independently re-certified to γ_d − 9.3e-10. Exploratory
track (006+007): the 135° stall dissolves — 001's pinch gap
[135.000°, 135.049°] was a sampler artifact (W(4,3) certified alive
inside it; the mechanism of 001's error is pinned), births follow
γ_birth = 180 − 90(a+b+1)/(a(b+1)) (SPECULATION; survives every
out-of-sample test at the sampler floor), consecutive windows touch,
and W(5,2) plus a genuinely non-W four-block word are certified
alive at γ = EXACTLY 135° — the "135° is a constructive barrier"
kill condition is refuted. Union-closed (previous cycle, all
skeptic-confirmed): gap (a′) refuted float-free by a 10-atom witness;
the i-AGGREGATED control survives and is the restated gap; P6′
proved at fixed n but naive n-uniform budgets killed (flattening
reverses past n ≈ 22; tax is δ-linear); gap (b) precisely stated and
surviving; route LIVE, ceiling 0.4315 vs record 0.38271. Standing
library: Singmaster census to 2.5×10^29; Erdős–Gyárfás cubics to
n=22 + girth-≥5 at n=24; lonely-runner k=8; Erdős–Straus
identity-poverty; graceful n ≤ 14; mahler-4d and billiards blind
censuses. As of 2026-08-04: 28 verified results, 8 recorded dead
ends, 13 problems (11 with at least one attempt), ~84 reusable tools,
and a machine-checked certificate layer (24 Lean theorems) over the
billiards route. Almost-mathieu and three-phase-conductivity are the
problems with no attempts (queue 18–19; run blind).

## Problem status

| Problem | Status | Budget | Active line |
|---|---|---|---|
| Erdős–Gyárfás | n=24 done | high | next: 2-connected C16-free test (lead 2); n=26 is a ~16h run |
| Union-closed (Frankl) | ROUTE LIVE; **GAP 1 FULLY CLOSED (020, certified; 022 cross-family implementation-level skeptic: confirmed, zero corrections; 024 fresh-session reviewer pass: confirmed, one reporting-level correction to 023)** | high | every odds-ratio control dead: pointwise (005), per-i (007), aggregated (016), margin-modulated signed (018), unsigned \|σ\| (020), λ-window-restricted (020); per-history weightings ruled out by R₊-pure witnesses (020) → live: Gap 2 mutual-information tax (009/011 — first step: run 020's witnesses against CR/TAX/ST), non-per-history bridges (020 lead 3), corrected assembly budgets (012); sweep leads: transfer-operator eigenfunction sharpening + strong-Rayleigh feasibility + Möbius/zeta LP stats (021, merging into 010's queue items) |
| Erdős–Straus | 601 resolved | medium | next: prove the identity-poor mechanism (why QR classes force low f) |
| Singmaster | census done | medium | next: Diophantine curve table (search-deeper is now low value) |
| Lonely runner | k=8 done | medium | next: k=9 scan (k=8 likely settled by Rosenfeld preprint) |
| Graceful trees | census done (n≤14) | low | possible next: mine the symmetric-spider seed; lobster verification at larger n |
| Collatz | first sweep done (001, MAP) | low (long shot) | next: graph-structure census of truncated digraphs (001 lead 1, queue 21); residue-Lyapunov and cycle-sieve families closed by 001's certified barriers |
| Triangular billiards | WINDOW LAW PROVEN for ALL (a,b) (015+016): necessity by 005+T2, sufficiency by the seven S-lemmas — proven in 016's skeptic pass on 015 (confirmed with corrections) — plus N1/N4/τ; **coverage of (90°,180°) closes: every obtuse angle has an explicit alive W(a,b) triangle** (via 013/014). Last SPECULATION under the chain: general-(a,b) endpoint bookkeeping (67-member exact + geometry-validated) | high | next: **skeptic pass on 016's seven proofs** (attack surface in 016 lead 1) + **formal all-(a,b) endpoint-bookkeeping ring pass** (016 lead 2) — together they retire the chain's last caveats; formalize T2 + the lemmas (Lean lane: L1 + Laurent block already FORMALIZED, 009–012); sampler blind spot + minimal-witness staircase (queue 12) |
| Mahler in ℝ⁴ | census done blind (skeptic-confirmed) | medium | next: close k=12–20 (falsifiable: no proper mask with P<11); run the same pipeline on {0,±1}³ for the n=3 spectrum comparison |
| Crouzeix | dim-3 census done (blind, skeptic-confirmed) | medium | next: hunt the published intermediate-maxima basins (informed; seed at Overton's ≈1.185/≈1.433 configurations) — the census's recorded gap |
| Maxwell equilibria | 24-equilibria witness SETTLED: skeptic-confirmed, escalation discharged, fold brackets 12/16 certified, centroid degeneracy exact (001+002) | high | next: harden verifier tiling check (queue 16, tier-0 fix); blind 3-charge strata map (queue 15); n=3 census hunting 4-vs-6 (queue 17); certified window edges + q\* sliver (002 leads 3-4) |
| Almost Mathieu (critical) | onboarded, no attempts | low (long shot) | harness exact to q ≈ 34 in seconds; first attempt is the rational-flux gap census (run blind) |
| Three-phase conductivity | onboarded, no attempts | medium | dual-route laminate harness ready; first attempt is the two-phase ground-truth self-test (run blind) |

## Attempt queue (next cycles pull from the top)

1. [union-closed] ~~Attack both surviving Gap-1 candidates (unsigned |σ|-control; λ-window variant; assembly restatement)~~ **DONE in 020 — ALL REFUTED, certified** (see TL;DR; swarm-guided designs, director-certified; 018's queue items (a), (b), (c) all settled negatively: the |σ|-control dies via the non-PSD sensitivity kernel at n = 4–5, the window variant dies via a₁-cancellation measures with negative a₂ at n = 4, and the assembly's aggregate secant requirement is violated in-window, with R₊-pure witnesses ruling out every chain-rule-compatible per-history weighting). ~~Replacement (a) skeptic pass, (b) Gap-2 check, (c) n = 8 tensor certificate~~ **(a)–(c) DONE in 022 same-day** — cross-family drafted-verifier skeptic confirms all 12 certificates digit-for-digit with zero corrections (implementation independence: different model family, different log₂ algorithm, zero shared code; reviewer independence still open — 022 was same-session), the n = 8 tensor is certified in-window, and **Gap 2's CR > 0 candidate survives both kill geometries at margins +0.40/+0.80**. New replacement, in order — items (b') and (c) already executed same-day in 023 (the CR adversarial pass found NO violation: floors +0.074…+0.22 at every H-bounded in-regime endpoint seeded at the kill geometries, and the signed/integrand form is certified negative at two rational tilts per witness): ~~(a) **fresh-session reviewer pass on 020/022/023**~~ **DONE in 024** — everything load-bearing confirmed (byte-identical re-run + third-path reimplementation inside every enclosure; verifier audit and both hand re-derivations sound; three fresh CR trajectories, zero violations); one reporting-level correction: 023's engine is nondeterministic (string-hash seeding), so its floors are per-trajectory (observed +0.0697…+0.0776 across four runs) and future CR engines need stable seeds + recorded `PYTHONHASHSEED` alongside the H ≥ ε guard; ~~pre-step per 023 lead 2: exercise the recipe's block-adaptive and half-mixing branches adversarially first~~ **DONE in 025 — both branches REFUTED as stated, and repaired**: half-mixing dies in-regime at every n (p_h → 1 sliver), iid-block dies Θ(n) via the component-wise ψ threshold, adaptive-block flatlines at CR = 0 on all-above-half mixtures; the extended mixing lemma (EM: ST = 0 for arbitrary ε-joints and s-laws) + concavity covers the whole genre with optimized q (proved except two float 1-D sweeps). New sub-items, in order: (a') **fresh-session reviewer pass on EM/EM-coverage + the 025 kills** (024 standard; until then the lemmas are candidate-VERIFIED); (b') **restate the recipe (Gap 3) with the corrected assignments** — ∅-mixing (optimized q) on {0}∪[1/2,1) and all-above-half mixtures, never iid-couple above-ψ mass, surprisal criterion as the branch test — and re-pose (TAX at p) against it; (b) **proof effort on (TAX at p), H-scaled** (023 lead 1): target CR ≥ c(p)·H(μ) on the generic genre with c ≤ 0.0387 (025's ratio envelope; the anneal grinds CR, not H), test battery = the floor instances from 023/024 runs (sharpest observed: +0.0697 at n = 6, 4 atoms — per-trajectory, see 024) plus 025's sliver instances (best-coupling CR ≈ +0.068 at H ≈ 0.97); (c) exact certification path for CR (mechanical now: 022 log₂ kit + the tensor certificate's dyadic-accumulation pattern) — ~~extend to one 025 kill per family~~ **branch kills + repairs DONE in 026 (dual-kit, six instances)**; still open for CR of general Sinkhorn couplings (needs certified h(z) sums over a fitted coupling, not just closed forms); (d) map the a₁-cancellation manifold (020 lead 5: closed-form a₂ on the manifold would settle the every-n tensor claim without computation); ~~(e) the uncovered branch cell: mixtures with both below-half and above-ψ components (adaptive-block positive on the one family tested at n ≤ 10; prove or kill — 025 lead 4)~~ **DONE in 027 — no kill (510-instance in-regime scan, CR > 0 throughout: the adaptive coupling's degenerate limit is comonotone, approached from above, so the fixed-cost sliver mechanism cannot bite), but CR/H(μ) → 0 as p_lo → 0 at bounded H (measured to 1e-5) — the adaptive assignment cannot carry any H-scaled floor on this cell; the new latent-mixing family (labels coupled, free both-lo probability q — the continuous closure of ∅-mixing) holds CR ≥ +0.22…+1.19 there and converges to the EM limit. Recipe v2 (b') now has all its inputs; new lead: prove adaptive CR ≥ 0 on the cell / an ST = O(p_lo·n) bound for latent-mixing (027 leads 1-2).**
2. [union-closed] Rebuild the assembly budgets per 012's corrections: the tax is δ-LINEAR (restate B2), the τ_half step needs s₀ ≤ 0.0843, corrected δ₀ ≈ 0.004. The open question is whether ANY n-uniform budget object exists — 008's conditional theorem survives structurally; its constants need the corrected inputs, and the budget census machinery (uc_pert.py + uc_pert_skeptic.py's orbit engine) is ready.
3. [union-closed] Sweep leads from 010, one per cycle: (i) pairwise-closure LP/degree-2 SOS certification on the n ≤ 4 census (all 4958 families) — kill: certified bound converges to ≈ 0.382; (ii) union-transfer-operator eigenvalue field — test log-supermodularity of λ_C on the n ≤ 5 census (pre-check the total-positivity records first); (iii) bipartite-MIS decomposition — geng census to n = 12, connectivity of extremals, compositionality across 1- and 2-cuts. Kill conditions recorded in 010.
4. [erdos-straus] Prove the identity-poverty mechanism: why does QR-class membership mod 840 force fewer Type I covering congruences? Start from 001's obstruction analysis + 002's rate data; target a theorem "f(p) ≥ g(N_typeI(p))" or a disproof.
5. [lonely-runner] k=9 near-tight scan: reuse lonely_runner.py (threshold near 1/9, feasibility analysis first — 8-tuples grow fast; consider restricting to accelerations/near-APs of known structures plus a bounded full scan).
6. [singmaster] Diophantine curve table: which equations C(n,j)=C(m,k) (small j<k) are resolved vs open, per the census lead that all in-range coincidences come from known families.
7. [erdos-gyarfas] Cycle-spectrum realizability census from the n≤20/22/24 data (which length-sets occur?) — standalone interest.
8. [collatz] Failed-approach taxonomy page (library showcase; pure writing + citation verification). Seed material now exists: 001's Family A/B closures and certified barriers, plus its [T]-grade reference list to verify.
9. [graceful-trees] Mine the symmetric-spider seed (LpH?GCAO??_@?A genre) at n = 15-16 targeted; lobster verification at larger n.
10. [crouzeix] **Hunt the intermediate-maxima basins** (from 001/002's recorded gap; informed — the blind census is spent). Seed local maximization AT Overton's published intermediate configurations (ratios ≈ 1.185 and ≈ 1.433 at n = 3; re-derive the seeds from arXiv:2105.14176's descriptions, not the [L] transcriptions) and map their basins with the 001 pipeline + 002's equal-sample escape probe: are they genuine local maxima under this design's probe standard, and how do their basins sit relative to the 001 start families that never found them? Falsifiable either way, and either outcome sharpens the landscape SPECULATION ({1, 2}-only) recorded in 001.
11. [billiards-triangles] ~~Parametric sufficiency + the birth side~~ **BROKEN OPEN in 015; the seven S-lemmas PROVEN in 016** (2026-08-16, skeptic pass on 015: confirmed with corrections — 015 §1's s=0 C-image identity is −N1 not −N1−N2, and the a ≤ 14 run never certified the death BC0/N4 row, repaired by 016's exact n4leg). Window law now holds for ALL (a,b); coverage of (90°,180°) closes via 013/014; τ ≠ 0 shown automatic (R = −4·sinθ·sin(aα)·sin((b+1)β)); N1 re-proven derivative-free. **Replacement: (i) skeptic pass on 016's seven proofs** (attack surface pre-registered in 016 lead 1 — the D-AC φ₁ case split, the birth-core constant chain, the B-BC X_r lower bound, the F5 λ<1 case); **(ii) the formal all-(a,b) endpoint-bookkeeping ring pass** (016 lead 2 = 015's residual SPECULATION — the A-fan offset two-liner with s formal); (iii) Formalize T2 + the lemmas (015 lead 3 / 016 lead 3; N1's new proof is mathlib-friendly). (iv) Cheap side tasks unchanged: the a = 1 column; birth-edge touching strictness (mirror of 004 C2); simplify the certifier per 016 lead 4 (drop N1-derivative + τ legs, add n4leg).
12. [billiards-triangles] **The coverage conjecture, and the sampler blind spot** (from 006/007, 2026-07-31; absorbs the old pinch-gap item — its motivating gap [135.000°, 135.049°] is CLOSED, W(4,3) is certified alive inside it): ~~006 reduced "every obtuse angle has an alive W member" to an elementary Diophantine statement (unproven; float-checked at 157 + 25 arcs over 90.5°–165°, zero failures). Prove it, using the birth law as a labelled input where needed.~~ **The Diophantine sublemma is PROVEN in 013** (2026-08-15; elementary interval-chain proof, no continued fractions, constructive witness verified exactly for every rational with denominator ≤ 300 in two arithmetic layers; **skeptic-CONFIRMED in 014**, 2026-08-16 — proof, witness, and checker all survive; two non-load-bearing prose errors of record: the γ = 144° spot-check's W(4,3) is wrong (witness is W(11,2)), and the Step 2 I(q−1,q) parenthetical asserts a false inclusion). Coverage of (90°, 180°) is now conditional SOLELY on the birth/window law (queue 11); the minimal-witness/word-length staircase (continued-fraction structure) remains open as 013 lead 3. Note the certificates so far are POINTWISE (007's C2): window-interval continuity on sub-arcs is float + SPECULATION law only, and per-triangle coverage of a whole arc is a different (open) question — the windows are x-slivers at the corners. Separately falsifiable (007 lead): every sampler in use accumulates only at the 90/j window edges, so an interior-pinch alive window would hide from ALL current designs — build one targeted interior-accumulation test before trusting any negative screen again.
13. [mahler-4d] Close the {0,±1}⁴ universe: k = 12–20 pairs (~30M orbits at k=12, improper fraction already 77% at k=9). Falsifiable: no proper mask with k ≥ 12 has P < 11. Needs the improper-detection shortcut or a streaming canonicalizer; see 001 lead 1. Cheap side quest, same pipeline: the {0,±1}³ census for the n=3 spectrum comparison (13 pairs, trivial) — does the non-Hanner gap grow or shrink with n?
14. [billiards-triangles] Coverage self-test: re-derive the acute and right-triangle cases as a scoped attempt record. Low value now that the harness self-test covers Fagnano and the orthic geometry and 001 mapped the obtuse side — take it only if something turns up that the certificate machinery cannot express.
15. [maxwell-equilibria] **Run blind.** First attempt: certified counts for structured 3-charge families beyond the harness self-test knowns — collinear with unequal charges (does the count stay 2 or drop?), isoceles families, a coarse (shape × charge-ratio) sweep. Deliverable is the count strata map, `EVIDENCE` scoped by grid and region. Every complete count must pass the index-sum identity; treat a violation as a harness bug, not a finding.
16. [maxwell-equilibria] ~~Skeptic review of 001~~ **DONE in 002** (confirmed-with-corrections; escalation discharged; brackets 12/16 certified; centroid Hessian exact). Replacement — harden the verification chain per 002's findings: (i) fix `verify_equilibria.py`'s tiling check with the axis-consistency reconstruction (prefix-freeness + Kraft is provably insufficient for axis-labelled paths — 002 R1) and add the two tamper-demo certificates from `explore/skeptic_verifier_gap_demo.py` as regression tests. Tier-0 change: keep it generic (a checker fix, no findings in prose), tests must pass. (ii) Guard `krawczyk()`'s rounded midpoint at low precision / deep min-width (002's latent hazard: margin 2⁻³¹ at 001's params, untriggered). (iii) Optional certified add-ons: binary-search the window edges q ∈ (4.389, 4.400)·10⁻³ and the sliver census q ∈ (4.389, 4.3968)·10⁻³ separating the 24→16 edge from the centroid degeneracy q\* = 81√2·t³/16.
17. [maxwell-equilibria] Three-positive-charge census hunting the open 4-vs-6 gap: after 15's strata map, target the strata boundaries (where counts jump) with unequal charges. Win condition kept in view by every run: any configuration certified with ≥ 5 isolated equilibria refutes the conjectured max 4 and is a result people have sought since 1873 — it must survive verify_equilibria.py and be reported as requiring escalation. Equal-magnitude configurations are settled (Tsai 2015, max 4) — spend no effort there.
18. [almost-mathieu] **Run blind.** Rational-flux gap census, all p/q with q ≤ 30: certify every gap open except the even-q central touching (re-derives van Mouche / Choi–Elliott–Yui in range; expected `VERIFIED`, scope = the q range), then the golden-mean convergent table — exact minimal-gap widths and q·|σ| along Fibonacci p/q as far as tooling reaches, against the (unproven) Thouless constant 32C/π. Onboarding smoke runs: q·|σ| = 9.2509 / 9.3199 / 9.3608 at q = 13 / 21 / 34 vs 9.3299 conjectured. Every record re-verified with `verify_bands.py` before ledger entry. Kill condition: if exact arithmetic stalls before q ≈ 100 even with a C kernel to the same contract, record the wall — no asymptotic claims from small denominators.
19. [three-phase-conductivity] **Run blind.** Two-phase ground truth first: rank-2 laminates attaining the 2D HS bounds exactly in ℚ, duality checks, series/parallel forms (expected `VERIFIED`, harness validation). Then the three-phase attainability map: fixed rational (σ₁,σ₂,σ₃), rational grid on the fraction simplex, bounded-rank laminate optimization (float screen, exact certification), gap-to-HS charted per cell (`MAP`/`EVIDENCE`, scoped by rank + direction set + grid). Nesi/Cherkaev improved bounds enter as marked transcriptions cross-checked against the papers' examples before anything is killed against them. Kill condition: bounded-rank optima plateauing strictly inside bounds across the whole grid = one negative-map record, then cap the budget.
20. [union-closed] **Push the kill's frontier** (cheap, from 016 leads): direct θ-optimization of the MU ladder at target n to find the minimal violating n (currently bracketed (32, 96]). ~~Extend the raw-weight ladder past n = 128~~ and ~~re-run the ladder at the λ-window boundary~~ **DONE in 018** — the raw ladder does NOT cross through n = 320 (decaying along each dilution branch, min +0.104 at n = 256; the 016 kill is entirely the re-weighting), and the window boundary is positive at every n ≤ 320 tried, θ re-optimized there included.
21. [collatz] ~~Graph-structure census of truncated Collatz digraphs~~ **DONE in 002 — REFUTED, and the whole graph-decomposition family is priced out with it.** Out-degree one makes `G_B` a forest plus one triangle, so: dominator tree = orbit tree (checked to B = 2^18), treewidth = 2 for every B, and the min cut is `|S_k| = J_{k−1} ~ 2^k/6` in closed form, set by `3r+1 > 2^{k+1}` alone with B dropping out. Novelty pre-step (001 gap 2) discharged: no dominator/treewidth Collatz work found, but the inverted-graph-is-a-tree fact is classical (Wirsching; Ebert arXiv:1905.07575 [T]), which explains the empty literature better than novelty does. **Do not rebuild the flow-cutter plan** — no truncated Collatz digraph has treewidth other than 2. Replacements, both from 002: (a) **explain the escape-fraction anomaly** — map `#{n ≤ 2^k : max orbit > 2^{k+1}}/2^k` to k ≤ 26 and find the mechanism for the k = 12→13 drop (0.609 → 0.379) and the ~0.404 plateau; falsifiable either as a window-alignment artifact with an exact description or as a real density worth stating, and it is the only part of the census with dynamical content; (b) **close the graph family with a reason** — the one un-killed variant is the reverse map (out-degree > 1); compute the density of binary-branching vertices in the reverse tree on [1, 2^k], k ≤ 24, against the 1/3 that `n ≡ 4 (mod 6)` forces. Expected kill; if it is not 1/3 the shape carries arithmetic and the family reopens.

## Verified results

- **[union-closed] Gap 1 fully closed: both surviving candidates REFUTED
  and the per-history-weighting family ruled out (skeptic-confirmed:
  022 cross-family implementation pass + 024 fresh-session reviewer
  pass)** (2026-08-19, attempt 020): the unsigned |σ|-weighted
  control is certified negative in-regime at n = 4 (7 atoms, t = 6/5 and
  3/2) and n = 5 — deficit at healthy margins via the non-PSD
  sensitivity kernel, plain aggregate certified positive on the same
  tables; the λ-window-restricted aggregated control is certified
  negative at n = 4 at three rational tilts (λ ≈ 0.070–0.111, window
  4.85, marginals ≤ 0.2575) on an a₁-cancellation measure (cancellation
  residuals ≤ 1.7e-16, non-product score 0.91), float-extending
  in-window to n = 8, 12 by block tensoring; the aggregate secant
  condition (the assembly's true requirement, identified convergently by
  both swarm families) is certified violated in-window at n = 4, 5; and
  three R₊-pure witnesses (every Plackett target certified < 1/2) with
  certified-negative σ-weighted census close the whole compatibility
  class. Enclosure widths < 1e-14 throughout; dichotomy exact; deep-anneal
  near-empty float negatives certified to be artifacts (sign flips under
  exact recomputation). Attack directions swarm-drafted (gpt-5.6-luna +
  gemini-3.7-flash, 12 briefs); all mathematics director-verified.

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

- **[billiards-triangles] Death law closed on both sides: parametric
  necessity for all (a,b), death = γ_d exactly for 20 members
  (skeptic-confirmed)** (2026-07-31, attempts 005+008): the half word
  composes in closed form to R₀·Rot_A(2aα)·Rot_B(−2bβ), so adjoining
  e^{iaα}, e^{ibβ} as formal Laurent variables proves I1–I3, the glide
  facts, and a new identity I4 (gate 2a+2) for ALL integers a, b in one
  exact polynomial check, specialized by ring homomorphism (bridge to
  geometry cross-checked exactly on 40 members to a = 40 by the
  skeptic's own off-torus recomposition). Lemma C is proven for all
  a ≥ 2, b ≤ a and Lemma D proven — both via one elementary lemma
  (c·cot(cs) strictly decreasing), retiring 004's residual. The case
  tree closes for all a ≥ b ≥ 1, a ≥ 2: the a ≤ 2b+3 restriction is
  GONE (new branches re-derived by hand twice; ~1.56M + 447k
  adversarial sign-scan points incl. 10 members with a > 2b+3, zero
  hits). Sufficiency: along a universal segment into the death corner
  every corridor condition is certified strict, so death(W(a,b)) =
  γ_d(a,b) EXACTLY (sup not attained) for 20 members; the skeptic
  independently certified alive by full 2n-gate interval unfolding (no
  glide reduction) down to γ_d − 9.3e-10. Five new members measured at
  γ_d to ~5e-12: W(6,3)→146.25, W(6,1)→127.5, W(8,2)→142.5,
  W(8,1)→1035/8, W(7,1)→900/7 (windows on the mirror half). Correction
  (008-C1): four of the five new alive certificates reach only
  γ_d − 1e-4, not the 1e-6 one sentence claims — superseded by the
  segment certificates. Open: parametric sufficiency (generic fan-gate
  margins; the death corner is 3-fold degenerate, not 2-fold), the
  a = 1 column, the birth side.
- **[billiards-triangles] The 135° stall dissolves: pinch gap closed,
  birth law found, alive at exactly 135° (skeptic-confirmed)**
  (2026-07-31, attempts 006+007): inverting the death-law machinery
  into a structure-class design search (plus a small new fact: ANY odd
  half word u makes u² a translation word, so the glide reduction
  applies class-wide) shows 001's pinch gap [135.000°, 135.0486°] was
  a sampler artifact: W(4,3) is certified alive inside it, and the
  mechanism of 001's error is pinned — the alive window is a corner
  sliver thinner than its grid's terminal 6e-4 clearance (001's
  knowledge-scoped gap claim stays true as written; what falls is its
  birth measurements, off by 10× the stated noise). Births obey
  γ_birth(a,b) = 180 − 90(a+b+1)/(a(b+1)) (SPECULATION as a law;
  survives three genuinely out-of-sample members incl. a > 2b+3 at the
  ~3e-11 sampler floor), making consecutive windows touch. W(5,2)
  (length 30, never measured before) and a genuinely non-W four-block
  word N1 (length 46) are certified alive at γ = EXACTLY 135° via
  rational apexes on the arc (2x−1)² + (2y+1)² = 2 (on-arc identity
  hand-verified in Fractions; the W(5,2) triangle is irrational-angled
  by Niven). All five exact certificates re-established digit-for-digit
  by a third corridor implementation. 003's theorem extends to (5,2)
  and (6,2). The kill condition is REFUTED: 135° is not a constructive
  barrier — every sampled arc in 90.5°–165° (157 + 25 independent) has
  a W member float-alive. Negative screens (nothing of length ≤ 26
  alive past 135°; extra letter blocks collapse death angles; clean
  angle laws confined to the two-fan factorization) are sample-bounded.
  New structural fact (007-C4): W(a,b) and W(b+1,a−1) are the SAME
  canonical word — 003's "two distinct words dying on the same arc"
  was a relabeling. Corrections: length-30 universe is 566 words (not
  811); gap certificates are pointwise, not "alive throughout".

- **[billiards-triangles] Lemma L1 FORMALIZED — the lab's first
  machine-checked certificate (skeptic-confirmed)** (2026-08-04, attempts
  009+010): 005's L1 (for fixed s ∈ (0, π/2], c ↦ c·cot(cs) is strictly
  decreasing on (0, 1] — the lemma that retires Lemmas C and D) proved in
  Lean 4 (mathlib v4.32.2), zero sorries, axioms exactly {propext,
  Classical.choice, Quot.sound}. Independent statement review confirmed
  fidelity to 005 hypothesis-by-hypothesis with no narrowing (Ioc,
  StrictAntiOn, and Real.cot semantics read from mathlib source), the
  cheat scan is clean, the root module provably elaborates the theorem
  file, and the reviewer's own rebuild + axiom audit went green. Scope:
  exactly the formal statement — L1 only; the I1–I4 Laurent block
  followed the same day (011/012, next entry). Source in
  problems/billiards-triangles/formal/.
- **[billiards-triangles] The full 005 Laurent block FORMALIZED
  (skeptic-confirmed)** (2026-08-04, attempts 011+012): the closed-form
  half-word composition, pair collapses, u² = translation, glide action,
  identities I1–I4/D1/D2, and the ℤ-specialization corollaries — 23
  theorems in Lean 4 (mathlib v4.32.2), zero sorries, axioms exactly
  {propext, Classical.choice, Quot.sound}, with (a,b) carried as actual
  powers of ∀-quantified generators so one check covers all integers.
  Statement review (012): all 23 statements mapped to 005 verbatim
  (angle monomials, factor-of-2 conventions, signs; permutation/sign
  hunts empty), conjugation-as-substitution confirmed as 008's star
  involution with nothing smuggled, the dropped i² = −1 hypothesis
  independently re-proved safe (the block is homogeneous in i), 31/31
  identities re-derived exactly in a from-scratch stdlib Laurent engine,
  independent rebuild + axiom audit green. Method notes: no Laurent-ring
  mathlib API needed (arbitrary field + field_simp + ring_nf) — the
  anticipated obstruction never materialized. Scope: the ring-model
  statements only; the geometry bridge (unfolding = closed form) is
  permanently carried by 005/008's informal cross-checks.

- **[billiards-triangles] Window law exact for a ≤ 14; birth necessity
  proven for all (a,b); margins in closed form** (2026-08-16, attempt
  015, skeptic pass pending): every distinct gate margin of W(a,b) is
  pivot ± sidelength·sin(linear form) — 2a+2b+2 margins, ring-verified
  as an exact bijection for 25 members — turning the alive criterion
  into an explicit sinusoid system. The birth edge is the N2/PB double
  wall meeting at (90/a, 90/(b+1)): 006's birth law is a corner, like
  the death law, and 013's Diophantine condition is literally
  θ ∈ (θ_d, θ_birth). PROVEN for all a ≥ b ≥ 1, a ≥ 2 (elementary,
  five steps, no intervals): alive ⟹ aα < 90 and (b+1)β < 90, hence
  θ < θ_birth — the A-fan margins force 2a all-negative sines along an
  arithmetic progression (step < 90°), pinning it to one component and
  bounding its span; the B-fan margins interleave to give
  sin(aα+rβ) > 0 for every r ∈ [−(b−1), b], killing the (b+1)β ≥ 360
  branch (worker-drafted, director re-derived, Step-5 branch repaired).
  Sufficiency: the death segment (ρ=2, whose alive range ends at the
  sharp wall t* = 45/(a(b+1)) — the window midpoint, explaining 005's
  ρ=2 uniformity) and the birth diagonal segment jointly realize every
  θ in the window; certified in exact interval arithmetic for all 104
  members a ≤ 14 (N1 by derivative+bisection, τ ≠ 0 included), float-
  verified to a = 60 with zero failures plus a 468k-sample adversarial
  scan above θ_birth (zero alive). Result: window(W(a,b)) =
  (γ_birth, γ_d) EXACTLY for every a ≤ 14 — closing 006/007's float
  births into exact statements — and coverage of (90°, 180°) is a
  theorem for every γ whose 013-witness has a ≤ 14, conditional beyond
  that on exactly seven remaining fan-margin lemmas (worker-drafted,
  SPECULATION until re-derived; five of twelve, including parametric
  N1 and the T2 pigeonhole, are director-verified). Swarm: 13 drafts,
  gpt-5.6-luna + gemini-3.7-flash, ≈ $0.16.

- **[billiards-triangles] The seven S-lemmas PROVEN; window law for ALL
  (a,b); coverage of (90°,180°) closes (skeptic pass on 015: confirmed
  with corrections)** (2026-08-16, attempt 016): all seven fan-margin
  lemmas of 015 §5 plus the τ row re-derived from scratch — the swarm
  drafts are not in the repo and were never read — and proven for all
  a ≥ b ≥ 1, a ≥ 2 at k = 1/2: D-AB/B-AB by progression-window bounds,
  D-BA/D-BC by single cosine comparisons, D-AC by an exact
  product-to-sum pairing (sole sign-varying factor φ₁ = c − A/2 +
  (b+3−a)t ≥ 0 for every member but W(2,2), which closes by two-line
  domination), B-AC/B-BC jointly via one new birth-core inequality
  sin α·sin(H−(a−b)t) > sin θ·cos(at)·sin((b+1)t), and τ via the new
  closed form R = −4·sinθ·sin(aα)·sin((b+1)β), making τ ≠ 0 automatic
  at every alive point of every member. T2 confirmed (with a one-line
  simplification: n = a gives sin(z_a) = −sin((b+1)β) > 0, killing the
  aα ∈ (90°,120°) branch); N1 confirmed by a new derivative-free proof
  ((2b+1)·cos(at)·cos((2b+1)t) > ab·sin(90/a) via tan u ≥ u and
  sin(nx) ≤ n·sin x). With 005-necessity + T2 + 013/014: **window(W(a,b))
  = (γ_birth, γ_d) exactly for every member, and every obtuse γ has an
  explicit alive W(a,b) triangle with apex angle γ.** Corrections to 015,
  neither load-bearing: the s = 0 C-image margin is −N1 (not −N1−N2),
  and the a ≤ 14 certificate run never evaluated the death-segment
  BC0 (= N4) row — true via 005's I4, certified exactly by 016's n4leg
  (104/104 members), with 015 §5's reduction table missing that family.
  Validation: 1200/1200 sign agreements with 002's independent corridor
  geometry, 62 fresh members ring-verified (union 67), all 32 proof-step
  inequalities scanned over a ≤ 40 dense + 20k random members to
  a = 3000 with zero violations, 453k-point T2 beyond-wall scan clean.
  Remaining under the chain: the general-(a,b) endpoint bookkeeping
  (015's narrow SPECULATION; formal ring pass queued) and the standard
  skeptic pass on 016's own proofs.

- **[billiards-triangles] The coverage Diophantine lemma PROVEN
  (skeptic-confirmed)** (2026-08-15/16, attempts 013+014): for every
  t ∈ (0,1) there exist integers a, b ≥ 1 with a+b < t·a(b+1) < a+b+1 —
  proved by an elementary overlapping-interval-chain argument (fixing
  q = b+1, the intervals I(a,q) chain-overlap exactly for a ≥ q and
  union to (1/q, 2/q); no continued fractions needed), with the
  constructive witness q = ⌊1/t⌋+1, a = ⌊(q−1)/(qt−1)⌋+1 verified
  exactly for all 27,397 reduced rationals with denominator ≤ 300 in
  two independent arithmetic layers, 2,494 boundary adversaries, and
  denominators to 2·10⁵⁰. The skeptic (014) re-derived every step by
  hand (both integer edge cases in the strictness bookkeeping; a
  cleaner minimal-a₀ closure of the union-endpoint question),
  re-implemented the checker from the record's formulas alone before
  reading the committed code (all counts reproduced; witness diff on
  the full Farey sweep: zero mismatches), and re-searched the
  literature (nothing found; stays labelled elementary/likely
  folklore). Corrections of record (014): the γ = 144° spot-check's
  W(4,3) is not a witness there (constructive witness W(11,2), minimal
  (6,3)/(4,5); W(4,3) is never the constructive witness for any t),
  and Step 2's I(q−1,q) ⊂ (1/q, 2/q) parenthetical is false (the
  interval is disjoint, lower endpoint exactly 2/q). Billiards
  meaning stays conditional on 006's SPECULATION birth/window law
  (queue 11); the lemma itself is unconditional.

- **[union-closed] i-AGGREGATED odds-ratio control REFUTED at large n
  (skeptic-confirmed)** (2026-08-04, attempts 016+017): the control —
  certified positive at n ≤ 7 by 013 and surviving 58 seeded attacks at
  n ≤ 32 (and by the parallel line at 014/015, next entry) — fails
  asymptotically: replicating the 10-atom witness's
  crash core Θ(n) times with orbit-optimized shared unit weights gives
  a ladder MU(n,r) whose aggregate falls like a − b·log n, certified
  negative in exact rational arithmetic at n = 96 (−7.6e-4, enclosure
  width ≤ 1.3e-20), n = 128 (−0.0154), n = 160 (−0.0242) at t = 4
  (λ = 2 exactly), all marginals ≤ 0.309 < 0.38271, sign robust to the
  alternative degenerate-history bookkeeping. 015 re-derived the
  quantity from 007/013's prose, re-implemented everything from scratch
  (own log₂ enclosure via interval squaring), rebuilt the family
  atom-for-atom, and matched all three certificates to 18+ digits;
  positive controls confirmed, the surviving raw-weight ladder upgraded
  to certified +0.157 at n = 96. Corrections: n = 128 atom count (252);
  and the violation lives at λ ∈ [2, 2.5], NOT in the 009/011 workable
  window (≈ 0.05 at n = 96) — so the ∀λ gap closes but the
  window-restricted variant is untouched (queue 1). Recorded as 016/017
  rather than 014/015 because it ran in parallel with an independent
  probe of the same object (next entry) that reached the 014/015 slots
  first; nothing else about either record was changed on merge, and the
  two agree at n ≤ 32.

- **[crouzeix] Dim-3 local-maxima census: no unknown extremal structure
  (blind, skeptic-confirmed)** (2026-08-04, attempts 001+002): 286-start
  census at n = 3, deg p ≤ 3 from seven structured start families, with
  certified rational enclosures at every endpoint — scoped precisely as
  **0 endpoints certified above 2** (182 enclosure upper ends exceed 2,
  so this is a no-refutation statement, not a global bound). Every
  probe-surviving local maximum is known structure: 13 genuine ratio-1
  maxima in the ice-cream-cone configuration (eigenvalue on ∂W(A),
  |p|-peak there; probed hard by the skeptic — no ascent) and 94
  Jordan-type endpoints climbing toward ratio 2; all four intermediate
  candidates confirmed as optimizer stalls by the skeptic's independent
  equal-sample probe. Skeptic re-certified 30/286 with a structurally
  different certifier (Sylvester-minor bisection, Bernstein segment
  bounds), reproduced both calibration anchors to the last printed
  digit, and confirmed the denominator discretization errs in the sound
  direction. Corrections: near-2 class 94 not 97; the 2021 paper is
  Overton alone. Scope: this design only — Overton's published
  intermediate maxima (≈ 1.185, ≈ 1.433) never appeared and their
  basins are untested (queue 10).

- **[maxwell-equilibria] Explicit five-charge witness: exactly 24
  nondegenerate equilibria, certified (both engines)** (2026-08-04,
  attempt 001): unit charges at e₁, e₂, e₃ plus charges 4367/1000000 at
  (1/3,1/3,1/3) ± (17/200)(1,1,1) — a rational embedding of the
  arXiv:2607.27197 family at ε ≈ 0.18, charge within 2·10⁻⁸ of their law —
  has exactly 24 isolated equilibria, all nondegenerate with certified
  index signs (10/−14, sum −4 = 1−n as Poincaré–Hopf demands). Complete
  census: 472k boxes tile the localization ball, zero unresolved leaves,
  enclosure widths 5·10⁻⁷–2·10⁻⁴; independently re-verified leaf-by-leaf
  (different arithmetic, interval Newton vs Krawczyk), PASS in 924 s.
  24 > (5−1)² = 16: an explicit, checkable instance of the refutation of
  Maxwell's conjecture — the preprint's own result is asymptotic with no
  explicit ε. Float `EVIDENCE`, sampled grids: the 24-count window exists
  only for ε ≲ 0.18 and spans ~0.5% in q. Escalation open: skeptic review
  of the shared driver (queue 16) before this is settled library fact.

- **[maxwell-equilibria] Witness settled: 001 skeptic-confirmed, fold
  brackets certified, verifier gap found and closed for this certificate**
  (2026-08-04, attempt 002): the escalation is discharged. A
  shifted+enlarged-region re-census (472,775 boxes, tree sharing no box
  with 001's) certifies exactly 24 again and PASSes the independent
  verifier; an 80-bit re-run reproduces 001's tree verdict-by-verdict
  (splits are grid-independent — precision perturbs nothing, the region
  shift carries the independence weight). 001's float predictions are now
  certified counts: 12 at q = 4360/10⁶ (index +4/−8), 16 at q = 4400/10⁶
  (+6/−10; 2.2× tree, near-degenerate). Fresh-code structural checks on
  001's certificate all pass in exact ℚ (tiling proof, 24 enclosures
  pairwise disjoint with min gap 1.22e-3, index sum −4, region ⊇ ball).
  Correction of record (R1): the checker's prefix-freeness + Kraft tiling
  test is unsound for axis-labelled paths — a tampered certificate hiding
  a real equilibrium PASSes it; coverage independence rested on the
  (audited, correct) driver. Closed here by a strictly stronger
  axis-consistency tree reconstruction: 0 violations over all 235,993
  paths. Exact centroid Hessian (001 lead 5): signature decided by
  6561t⁶ vs 128q², degeneracy q\* = 81√2·t³/16 ≈ 4.3968·10⁻³, certified
  index flip +1 → −1 across it (matches C/D certs). Hardening leads in
  queue 16.

- **[union-closed] Aggregated odds-ratio control survives the mandated
  large-n probe (skeptic-confirmed)** (2026-08-04, attempts 014+015): the
  restated gap A(μ,λ) ≥ 0 (013's i-aggregated control) holds on every
  instance tried — ~1900 total: 840 orbit-census (14 families × λ grids
  incl. the 4.847/(n−3) window law, n ∈ {10,…,32}), 432 perturbative, 426
  sparse structured (the 007/013 witness genre re-laid-out per n,
  multi-gadget stacks, crash/mirror), 82 aggregate-objective hill-climbs,
  float trends to n = 64. Fifteen exact-rational certificates (Fraction
  end-to-end, directed log₂ enclosures, widths ≤ 5.9e-12), softest
  +5.988e-7 at n = 32. Minimum census margin GROWS monotonically with n
  (crashmix at small λ: +4.9e-4 at n = 10 → +1.2e-3 at n = 32) — the
  012-style reversal was hunted for and not found, including on 012's own
  genre. Structure: product measures are the exact equality set for
  λ > 0 (inclusion proven; converse open — and false at λ = 0), so
  perturbative departures rise ~quadratically and the tightest census
  points are boundary noise. Skeptic (015): independent certifier on a
  different algorithm re-certified 5/15 with zero disagreement (~200×
  tighter), orbit formulas re-derived and cross-checked at n = 8, 9, dip
  row (+1.4e-7 at n ≈ 48, recovers by 64) reproduced digit-for-digit and
  shown not to be a grid artifact, ~116 new adversarial instances in
  regions 014 could not see — all positive. Four reporting-level
  corrections (fitted-order range d^1.39–d^2.06, converse overstatement,
  docstring runtimes, instance count 1820). Finite-n `EVIDENCE`; proof
  program in queue 1. **Superseded where they differ by 016/017**: the
  parallel line built the one family this battery did not — the unit
  ladder MU(n,r) with re-optimized shared weights — and certified the
  aggregate negative at n = 96–160, so the ∀λ statement is refuted and
  the proof effort this record licensed is off. Everything certified
  here still stands: n ≤ 32 is positive on both searches, and the
  equality-set structure is what makes the λ-restricted variant the
  live question.

- **[union-closed] Margin-modulated OR control REFUTED in both signed
  readings at n = 7 (skeptic-confirmed)** (2026-08-06, attempts 018+019):
  007 lead 2's "h-sensitivity-weighted" candidate splits into three
  precise statements; the assembly-exact secant form
  E[h₂(z̃) − h₂(z_{2^λ}(x̃, ỹ))] and the derivative-at-target form are
  both negative on 007's own 10-atom witness — the instance the weighting
  was invented to survive — with the secant form certified float-free at
  t = 181/16 (MM_sec = −0.0222546582…, in-regime, dichotomy exact, the
  first certified enclosure of a nonlinear census functional in the
  library: dyadic-bisection Plackett roots + directed-log₂ binary
  entropy, margin-pair cached). 019 re-certified from scratch (own
  census from 007 §1's prose, exact-rational Newton with symbolically
  sign-checked brackets, interval-squaring log₂ — no shared code) to
  width 5.9e-30, matching; re-derived the secant identity, the
  product-measure zero, and the degenerate-history invariance by hand;
  confirmed no natural signed variant survives (realized-OR weighting
  dies too, −1.458); and newly certified the 2-digit tidy witness
  (−0.0222736748). Perturbation-stable (20/20 at 3%, plus 40/40 on the
  skeptic's fresh seeds), violating for every λ ∈ [1, 5] (crossing
  λ ≈ 0.83, 019 C4). Mechanism: the weight is signed and the surplus
  histories sit past the entropy peak. Survivors, EVIDENCE-scoped: the
  unsigned |σ|-weighted control (positive everywhere tried incl. 016's
  certified kill instances; θ re-optimized against it, transfers RISING
  with n; softest +0.0101 at a free n = 8 endpoint — only 3 of 8 free
  climbs ended in-regime, 019) and the λ-window-restricted control
  (positive at λ_win, λ_win/2, λ_win/4 across raw and θ\* ladders
  n ≤ 320, θ re-optimized at window λ, joint (θ, λ) climb converging to
  the trivial λ → 0 boundary; a₁(n)·n grows 4.7 → 18.1 θ\*, ≈ 66 raw).
  016 lead 3 settled: the raw-weight ladder never crosses through
  n = 320 (min +0.104 at n = 256). Corrections of record (019): raw
  ladder MM signs flip positive on the n ≥ 96 dilution branch (018's
  "negative at every n" is θ\*-only); the golden-ratio aside is false
  as stated (ρ = 1 crossing at marginal 1 − 2^{−1/2}, not (3−√5)/2);
  witness MM_sec crossing λ ≈ 0.83 not ≲ 0.7; §F certificates landed
  after review close (not covered by 019); runtime understated.

## Insights / cross-problem notes

- Prove the lattice, not the inequalities (2026-08-16, billiards 015):
  when a family of trig-margin conditions shares a lattice structure
  (fan index linear in the angles), the global statement is a
  pigeonhole about sign-constant arithmetic progressions — one span
  bound replaces 2a+2b interval estimates, and the "hard case tree"
  evaporates. New mechanism tag
  `sign-constant-progression-pigeonhole`; candidate reuse: rotation-
  number window arguments in almost-mathieu (015 lead 6). Same cycle,
  swarm calibration data point: 13 proof drafts for ≈ $0.16, the two
  hardest (parametric N1, the T2 case tree) both survived full
  re-derivation — with one repair the draft made findable — while a
  worker "REFUTED" verdict was a correct catch of a sign-convention
  gap in the brief. Drafting is now cheap; director re-derivation is
  the bottleneck, exactly as SWARM.md prices it.

- Signed sensitivity weights backfire (2026-08-06, union-closed 018/019):
  weighting a deviation control by the derivative of the gain it feeds
  does not neutralize an adversary that hides its deficit at degenerate
  margins — the derivative is SIGNED, and in the record regime the
  surplus histories sit on the decreasing side of the entropy curve, so
  the modulation flips the control's biggest positive terms instead of
  suppressing the deficit. Salvaging positivity needs the unsigned
  weight, which decouples the functional from the chain rule it was
  meant to serve — a modulated control can be plausible-looking and
  assembly-irrelevant at the same time. Same cycle, tooling note: a
  NONLINEAR census functional (algebraic root + entropy per history) is
  certifiable at ladder scale by enclosing the root with dyadic
  bisection and caching by the margin pair — unit exchangeability
  collapses thousands of rows to a handful of enclosures.

- Blind mode, data point #3 (2026-08-04, crouzeix 001/002): with prior
  art physically absent, the census blind-rediscovered both known
  extremal structures (the ratio-1 ice-cream-cone configuration and the
  ratio-2 Jordan limit). All three blind censuses to date (mahler,
  billiards, crouzeix) have now independently recovered published
  structure from the harness + problem statement alone. Caveat recorded
  by 002 and inherited by future blind assignments: the task *framing*
  (slice, cost-first sizing, design-difference requirement) came from
  the orchestrator's queue, so the blind label certifies prior-art
  absence, not task-selection independence.
- Escape probes need sample-count parity (2026-08-04, crouzeix 001/002):
  001's deep probe compared 128-sample escape values against 512-sample
  baselines — on a noisy nonsmooth objective that mismatch can
  manufacture ascent and turn a genuine local maximum into a "stall".
  002 re-ran every stall verdict at 512-vs-512 (all survived, so the
  census's conclusions stand), but the lesson generalizes: any
  stall-vs-maximum verdict must compare like against like, and a
  landscape census without an escape-probe pass at matched sampling
  reports optimizer artifacts as mathematics.
- Unit replication + re-weighting is a portable adversary genre
  (2026-08-04, union-closed 016/017): a family certified positive at
  small n can fail after replicating its adversarial core Θ(n) times
  and re-optimizing shared weights — the zero crossing near n ≈ 90 was
  invisible to every fixed-n search below it. Second vindication of
  probe-before-proof. Companion lesson from 017's correction C2: a kill
  must be checked against the workable parameter window before
  declaring a branch fully closed — 016's kill closes the recorded ∀λ
  statement; the window-restricted variant needed its own verdict and
  didn't get one from this family.
- Ripple scan of the Astra results ran (2026-08-04, informed):
  all ten results characterized (press-transcribed only — openai.com, the
  certificate repo, and erdosproblems.com were unreachable, so every scope
  statement is [L]-grade and must be re-checked against a manuscript
  before any transfer work) and scanned against all ten problems. **No
  purchase anywhere**: no recorded gap unblocked, no queued route
  foreclosed. Both queue-flagged bites are reasoned misses — Ehrhart vs
  mahler-4d fails on direction (upper vs lower bound), symmetry (Ehrhart's
  content is the non-symmetric case; the symmetric case is classical
  Minkowski), and hypotheses (Mahler has no lattice-point constraint);
  extremal/Ramsey vs erdos-gyarfas fails on regime (Turán/dense and
  complete-graph coloring vs bounded-degree cycle forcing — a cubic graph
  sits below every Turán threshold in play). One conditional seed,
  SPECULATION, deliberately not queued: if the Astra Ehrhart proof turns
  out complex-analytic (Berman–Berndtsson lineage — the toolbox
  historically shared with Bourgain–Milman), it becomes the seed for an
  analysis-lens /ideate sweep on mahler-4d; check the method (one hour)
  when a manuscript mirror is reachable. Two hazards recorded for future
  scans: the "Erdős–Gyárfás generalized Ramsey (p,q)-coloring" problem in
  Astra-adjacent press is NOT our cycle conjecture; and press sources
  conflict on which of Erdős #146/#180 is compactness vs degeneracy — do
  not cite that numbering from this scan.
- First FORMALIZE.md pass ran end-to-end (2026-08-04, billiards 009+010),
  ops notes: the toolchain is viable in-session (elan, lean4 v4.32.2,
  prebuilt mathlib via `lake exe cache get`, ~7.4 GB — never build mathlib
  from source; `formal/.lake/` must be git-ignored, and the docs lint now
  skips vendored deps). mathlib v4.32.2 has Real.cot but no derivative API
  for it — an inline quotient rule is three lines, not an obstruction. The
  statement-review step earned its keep procedurally: the semantic
  verdicts (interval endpoints, StrictAntiOn direction, cot convention)
  were read from mathlib source rather than assumed — exactly the
  informal→formal bridge risk the lane exists to cover.
- External event + process change (2026-08-04): OpenAI announced (2026-08-02)
  ten results on decade-plus-open problems (non-sofic groups, Connes rigidity
  counterexample, Ehrhart volume conjecture, three Erdős-catalog problems
  incl. multicolor Ramsey #183, sphere-packing and coding bounds, arithmetic
  circuits, quantum parallel repetition, closest-vector, extremal graphs),
  each with a model reasoning walkthrough and a Lean 4 certificate at zero
  sorries; unrefereed as of this note. Adopted here: (1) literature-novelty
  check added to the skeptic pass — the instructive contrast is OpenAI's own
  October 2025 episode, where correct derivations already in the literature
  were announced as new; (2) `FORMALIZED` status + `docs/FORMALIZE.md` lane
  for kernel-checked certificates above `VERIFIED`, with the statement-review
  step carrying the bridge risk we already hit in formal-parameter
  specialization; (3) ripple scan of the ten results queued as item 15,
  formalization pilot as item 16. No mathematical state changed by this note.
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

- Sampler blind spots are one-sided and shared (2026-07-31, billiards
  006/007): 001's birth angles were wrong by 0.0486° — 10× its stated
  noise — because the alive window near birth is a corner sliver
  thinner than the grid's terminal clearance, and its "deaths only
  under-estimated" hedge does not transfer to births. Cheap
  validations that caught it: floor-tracking (re-measure at several
  sampler depths and check whether the value tracks the sampler floor)
  and certifying an exact point inside the claimed-empty region.
  Standing residual: every billiards sampler in use accumulates only
  at the 90/j window edges, so an interior-pinch window would hide
  from all of them (queue item 12).
- Formal-parameter specialization proves infinite identity families in
  one check (2026-07-31, billiards 005/008): move the family
  parameters into the exponent lattice (adjoin e^{iaα}, e^{ibβ} as
  formal Laurent variables), verify the identity once by exact
  polynomial subtraction, specialize by ring homomorphism to every
  integer (a,b). The entire risk concentrates in the hand-derived
  bridge (geometry = closed form) — the formal check cannot see a
  wrong bridge — so cross-check the bridge exactly on many members;
  the specialization direction itself is safe (formal zero ⟹
  geometric zero). New mechanism tags: formal-parameter-specialization,
  structure-class-design-search.
- Canonicalize before counting evidence (2026-07-31, billiards 007):
  W(a,b) and W(b+1,a−1) are the same canonical word, which made a
  claimed "window pairing" vacuous (the fit spanned 9 canonical words,
  not 11) and turned "two distinct words dying on the same arc" into a
  relabeling. Any law fitted over family members must first quotient
  by the word-canonicalization symmetry.

- Onboarding (2026-08-04): maxwell-equilibria added per `docs/PLAN-em-problems.md`,
  and the Phase-0 literature check changed the framing mid-plan: Maxwell's
  general (n−1)² conjecture was refuted by a preprint SIX DAYS before
  onboarding (arXiv:2607.27197, ≥ 24 nondegenerate equilibria from 5 charges;
  companion 2607.28785 takes 3 positive charges from 12 to 6). The refutation
  is perturbative with no explicit certified witness — queue 16 targets
  exactly that gap, which our certificate machinery is unusually suited to.
  Frame all attempts against the surviving questions (n=3 max 4-vs-6,
  explicit witnesses, growth of the max), never against the dead statement.
  thomson-sphere (second problem in the plan) deferred to a follow-up
  onboarding session. New harness pattern worth reusing: certificate =
  subdivision tree with split paths + per-leaf certificates, so the
  independent checker can re-establish every leaf AND verify the tiling
  combinatorially (prefix-freeness + Kraft equality) — coverage claims stop
  being trust-me. Also: decimal's squareRoot ignores context rounding
  (always half-even, per spec) — directed sqrt bounds must be established
  by hand; caught by the verifier's own validation suite.

- Prefix-freeness + Kraft equality does NOT imply a tiling for
  axis-labelled split paths (2026-08-04, maxwell 002): {"01","11"} is
  prefix-free with Kraft sum 1 and tiles nothing — so a checker that
  verifies coverage combinatorially from split paths accepts non-covering
  trees, and a tampered certificate that hides a real equilibrium PASSes.
  The lesson generalizes: an "independent" verifier must re-establish
  coverage from leaf geometry (axis-consistency reconstruction), not from
  combinatorial summaries the driver emitted; otherwise coverage
  independence silently reduces to trusting the driver. Same cycle, same
  family: changing working precision does not perturb a width-driven
  subdivision tree at all (identical 471,985 boxes at 64 and 80 bits) —
  region offsets, not precision, are the real tree-independence test.
- Probe-before-proof paid for itself twice (2026-08-04, union-closed
  014/015): the mandated large-n probe reshaped the proof target — products are the exact equality set of the
  aggregated control (λ > 0), so the conjecture is boundary positivity of
  a functional vanishing identically on products, and the census minima
  are boundary noise rather than adversaries. It did NOT clear the proof
  gate, though this record read it as having done so: the parallel probe
  at 016/017 found the crossing past n ≈ 90 that this battery's families
  could not reach. Third lesson, from the pair: "survives every genre I
  built" is a statement about the genre list, and two independent probes
  of one object are worth more than one probe with more instances. Two engine-methodology
  notes that made the certificates cheap: integer-scaled rows make the
  aggregate denominator exact, and a from-scratch certifier on a
  different enclosure algorithm (binary digit-extraction log₂ vs atanh
  series) is a fast, decisive skeptic instrument (5 re-certifications,
  nested intervals, zero disagreement).

## Dead ends

- **[union-closed] The entire per-history odds-ratio-control family
  (Gap 1, all restatements)** (2026-08-19, attempt 020; skeptic pass
  pending): dead — pointwise (005), per-i averaged (007), i-aggregated
  ∀λ (016), margin-modulated signed (018), unsigned |σ|-weighted (020),
  and λ-window-restricted (020) all carry certified in-regime
  counterexamples, and R₊-pure certified witnesses rule out EVERY weight
  function that equals the sensitivity where it is positive (the
  chain-rule compatibility class). Two reusable obstructions: (i) any
  margin-dependent weight turns the first-order coefficient into a
  kernel-weighted Gram form, and the sensitivity kernel is NOT PSD
  in-regime (diagonal collapses at the z = 1/2 kink) — so weighting
  cannot be "shaped" out of trouble; (ii) first-order (Theorem-C)
  protection is voidable by prefix-signed cancellation on measures far
  from products, where the λ² coefficient has no sign protection. Killed
  as stated for per-history weightings of log₂OR deviations; NOT killed:
  coordinate-level budgets with cross-history structure, λ-integrated
  forms with weights outside the compatibility class, and Gap 2's
  chain-rule assembly value (009/011).

- **[collatz] Unlabelled graph invariants on the forward Collatz digraph**
  (2026-08-15, attempt 002): dead, and dead for *any* iterated map, not
  just this one. Because the map is a function the digraph has out-degree
  one, hence is a forest plus one cycle, hence the dominator tree is the
  orbit tree, the treewidth is the constant 2, and the minimum cut has a
  closed form set by a one-step size inequality. Do not re-attempt
  dominator, treewidth, cut-width or flow-cutter routes on this object, and
  do not expect arithmetic labels to rescue them — the cut is an interval of
  odd integers and every residue statistic on it is uniform by
  construction. The one-paragraph collapse argument is the reusable part
  (`functional-graph-invariants`): it prices any graph-decomposition route
  before machinery is built. NOT killed: labelled/weighted objects on the
  same digraph, where the arithmetic rides in the weights rather than the
  shape.
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
  margin-modulated controls (queue item 1). [Update 2026-08-04: the
  i-aggregated replacement is now also dead — see below.]
- **[union-closed] i-AGGREGATED odds-ratio control, unrestricted-λ form
  (013's restatement)** (2026-08-04, attempts 016+017): dead at large n —
  the replicated-witness ladder MU(n,r) with orbit-optimized shared unit
  weights is certified negative in exact rational arithmetic at
  n = 96/128/160 (λ = 2, marginals ≤ 0.309), independently re-derived to
  18+ digits by 017's from-scratch engine. Do not re-attempt any ∀λ
  aggregated form; the ladder + re-weighting genre is the reusable
  adversary (`explore/uc_agg_ctrl_probe2.py` builds it,
  `explore/uc_agg_ctrl_skeptic.py` certifies independently). NOT covered
  by the kill: the margin-modulated control and the λ ≲ c/n
  window-restricted variant (the violation sits at λ ∈ [2, 2.5], far
  above the workable window) — both live, queue item 1.
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
