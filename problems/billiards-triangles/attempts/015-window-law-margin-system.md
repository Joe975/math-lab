# 015 — The window law as a closed-form margin system: birth corner found, birth necessity proven, both window halves certified

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-08-16
- **Mode:** informed (follows 005 leads 1+3, 006/007, queue 11; uses the
  013/014 Diophantine lemma as the coverage interface)
- **Type:** structural closed forms (exact ring) + parametric proofs +
  certified interval computation + swarm-drafted proof candidates
- **Tools:** `explore/wmargins.py` (closed-form margin system, ring-verified;
  float alive test), `explore/wmargins_certify.py` (exact interval
  certification of both canonical segments; reuses 005's certified trig),
  `explore/wmargin_brief.md` + `explore/wmargin_tasks.txt` (swarm briefs).
  Data: `data/wmargins_certify_a14.json`.
- **Workers:** gpt-5.6-luna (effort medium, 10 briefs) and
  gemini-3.7-flash (effort medium, 4 briefs, 1 blocked by a spurious
  content filter), ~102k tokens total, ≈ $0.16. Division of labor per
  `docs/SWARM.md`: workers drafted proofs of named elementary lemmas;
  every statement entering this record as *proven* was re-derived
  line-by-line by the director (provenance marked per lemma below);
  un-re-derived drafts are recorded as SPECULATION with their claimed
  ranges.
- **Sources:** prior attempts 003–008, 013/014 (this repo). No external
  literature used.

## Approach

Queue 11 asked for (i) parametric sufficiency of the death law and (ii) a
birth-law theorem, "which gate pair binds at the birth edge". 005 had
reduced per-member aliveness to ~2a+2b+1 gate-margin sign conditions and
noted the fan margins are short trig polynomials via the prefix maps —
but certified them per member by bisection, leaving the parametric bound
open (its obstruction 1).

The move here: make every margin CLOSED FORM first, and only then prove
inequalities. The prefix maps of 005 Step 3 give, for each fan gate, far
endpoint = pivot + (side length) · (unit phase linear in the fan index),
so each margin is pivot margin ± sidelength·sin(integer combination of
α, β, aα, bβ). Why this rather than extending 005's per-member interval
campaign: a sinusoid system can be *reasoned* about globally (pigeonhole
over arithmetic progressions of angles), which is exactly what the
per-member route cannot do.

## What was done

### 1. The complete closed-form margin system (exact, ring-verified)

With δ = e^{i(−aα+(b−1)β)}, θ = α+β, and writing
N2 = cos(aα)sin((b+1)β)sinθ (identity I2), PB = cos((b+1)β)sin(aα)sinθ
(identity I3, = p(B) − m), the 2a+2b+2 distinct gate-endpoint margins
p(·) − m of the half word are exactly:

    A-fan pivot:              N2                                (> 0 req.)
    gate-1 far endpoint:      PB                                (> 0 req.)
    A-fan C-images s=0..a-1:  N2 + sinβ·sin((a-2s-1)α-(b+1)β)   (< 0 req.)
    A-fan B-images s=1..a:    N2 + sinθ·sin((a-2s)α-(b+1)β)     (< 0 req.)
    B-fan C-images j=0..b-1:  -PB + sinα·sin(aα+(b-2j)β)        (> 0 req.)
    B-fan A-images j=1..b:    -PB + sinθ·sin(aα+(b+1-2j)β)      (> 0 req.)

Structural facts established on the way: the B-fan pivot satisfies the
exact antisymmetry p(B₁) = 2m − p(B) (so its margin is −PB); the B-fan
pivot IS the A-fan B-image s=a (same unfolded point; equivalently
−I3 = I2 − sinθ·sin(aα+(b+1)β), a sine-addition identity); the s=0
C-image margin is −N1−N2 (I1) and the j=0 C-image margin is N4 (I4) — so
N1 and N4 are rows of the same system, not extras. Derivation of the fan
offsets (from 005's prefix maps R₀Rot_A(2sα), R₀Rot_A(2aα)Rot_B(∓2jβ)):
z − pivot factors as e₀e^{−2isα}·conj(X) (A-fan) resp.
−sinα·e^{i(−2aα+(2j−1)β)} and −sinθ·e^{i(−2aα+2(j−1)β)} (B-fan).

**Verification:** `wmargins.py selftest` — for 25 members up to W(15,7),
the closed forms match the symbolic engine's distinct differences as an
exact ring BIJECTION (zero tolerance), and the count is 2a+2b+2; the
float alive() verdict agrees with 002's independent corridor-width
implementation at 400/400 random near-corner triangles. The closed forms
for general (a,b) follow from 005's proven prefix-map structure plus the
endpoint bookkeeping above; the bijection is machine-checked per member,
25 members. [The general-(a,b) endpoint bookkeeping is short and
mechanical but is stated here from a hand derivation + 25-member exact
check, not a formal all-(a,b) proof; SPECULATION only in that narrow
sense, same status 005's Step 3 had before its formal pass.]

### 2. The birth corner: the window is a two-wall corner pair

Maximizing θ = α+β over the alive region (local search seeded on the
known-alive segment, 10 members) converges to θ_birth = 90(a+b+1)/(a(b+1))
with binding margins ALWAYS {N2, PB}. Mechanism, now obvious from the
closed forms: alive requires (in the principal branch) cos(aα) > 0 AND
cos((b+1)β) > 0, i.e. aα < 90 and (b+1)β < 90, so

    θ < 90/a + 90/(b+1) = θ_birth,   attained at the corner
    (α*, β*) = (90/a, 90/(b+1)).

The death corner (90/a, 90(a−1)/(a(b+1))) is where the walls N1/N2/N4
meet (005); the birth corner is where N2/PB meet. 006's SPECULATION birth
law γ_birth = 180 − θ_birth is exactly this corner. Windows-touch is now
one line: θ_birth(a+1,a) = 180/(a+1) = θ_d(a,a). And 013's Diophantine
condition a+b < t·a(b+1) < a+b+1 is LITERALLY θ ∈ (θ_d, θ_birth) for
θ = 90t — the interface needs no translation at all.

### 3. Birth-side necessity: PROVEN for all a ≥ b ≥ 1, a ≥ 2

**Theorem (T2).** If the corridor of W(a,b) has positive width at a
triangle with angles α, β > 0, α+β < 90, then aα < 90 and (b+1)β < 90;
hence θ < θ_birth, i.e. γ > γ_birth(a,b) = 180 − 90(a+b+1)/(a(b+1)).

Proof shape (drafted by gpt-5.6-luna on the [T2] brief; re-derived in
full by the director, with one repair): the A-fan margins force
sin(zₙ) < 0 for the 2a-term arithmetic progression zₙ = aα−(b+1)β−nα
(n = 1..2a), step α < 90, so all zₙ lie in ONE negative-sine component,
giving the span bound 2aα − α < 180 and hence aα < 120; the branch
aα ∈ [90, 120) then forces (via cos(aα) ≤ 0 and the pivot margins)
(b+1)β into (270, 360) mod 360, which splits z₁ and z_{2a} into
different negative-sine components — contradiction. So aα < 90. Then
sin((b+1)β) > 0, and the m ≥ 1 branch of (b+1)β ∈ (360m, 90+360m) dies
by the SAME pigeonhole applied to the B-fan progression: the B-fan
margins jointly give sin(aα + rβ) > 0 for EVERY integer r ∈ [−(b−1), b]
(the two families interleave parities), a 2b-term progression with step
β, so (2b−1)β < 180, while (b+1)β > 360 needs β > 360/(b+1) — and
360(2b−1) < 180(b+1) fails for every b ≥ 1 (b = 1 separately: 2β > 360
contradicts β < 90). [The m ≥ 1 repair is the director's; the worker's
Step 5 was flagged as imprecise and replaced.] ∎

With 005's death-side necessity this gives, for ALL a ≥ b ≥ 1, a ≥ 2:

    alive  ⟹  θ ∈ (θ_d, θ_birth)      (the window is at most the law).

### 4. Window sufficiency: two canonical segments, certified

The alive region reaches both corners. Two explicit segments:

    DEATH segment  (90/a − t, 90(a−1)/(a(b+1)) + 2t):  θ = θ_d + t
    BIRTH segment  (90/a − t, 90/(b+1) − t):           θ = θ_birth − 2t

both with t ∈ (0, 22.5/(a(b+1))]. Together they realize EVERY
θ ∈ (θ_d, θ_birth) (the death segment covers the lower half-window
openly, the birth segment covers the upper half INCLUDING the midpoint
θ_d + 45/(a(b+1))). Float sweeps: death segment alive on the whole range
for all 1,829 members a ≤ 60 — in fact all the way to the PB wall at
t = 45/(a(b+1)), where (b+1)β reaches 90 exactly (the sharp endpoint;
ρ=2 was 005's lucky accident: its segment exits the alive region
EXACTLY at the window midpoint); birth segment likewise 1,829/1,829;
window coverage checked at 3,200 (member, θ) pairs, zero failures;
adversarial search for alive triangles at θ ≥ θ_birth: 468,000 samples,
zero hits (consistent with T2's proof).

**Certified:** `wmargins_certify.py sweep --amax 14` — for ALL 104
members with a ≤ 14, every margin sign on BOTH segments is certified in
exact rational interval arithmetic (005's certified-trig layer; N1 via
the derivative-plus-bisection treatment; N2/PB/AB_a by rational range
checks — AB_a on the birth segment is cos(at)sin((b+1)t)sinθ > 0 by the
sine addition formula, exact), and Re(conj(δ)τ) ≠ 0 on both segments.
So: **for every member with a ≤ 14, the window of W(a,b) is EXACTLY
(γ_birth, γ_d) as a set of realized angles** — necessity by T2+005,
sufficiency by certificate. This includes every member 006/007 measured
and closes their float birth measurements into exact statements.

### 5. The parametric margin lemmas (towards ALL (a,b)): status table

The sufficiency for all (a,b) reduces to 9 elementary lemma families
(fixed k = 1/2, i.e. t ≤ 22.5/(a(b+1))). Status — **P** = proven (a
complete derivation the director has personally verified line-by-line;
worker provenance noted), **S** = SPECULATION (worker-drafted proof with
claimed range, not yet re-derived here; queued for the skeptic):

    N2, PB > 0 both segments      P  (rational ranges; director)
    AB_a birth (addition formula) P  (director)
    B-BA (B-fan A-images, birth)  P  (director; one-line comparison
                                      using b ≤ a)
    N1 > 0, death, k=1/2          P  (gemini-3.7-flash draft, director
                                      re-derived every step; also proves
                                      N1'(0) > 0 parametrically — the
                                      (2b+1)sinβ_c > a·sin(90/a)·sin(bβ_c)
                                      sub-lemma — via sin(x)/x
                                      monotonicity + a·sin(90/a) < π/2)
    T2 case tree                  P  (luna draft + director repair, §3)
    D-AC (death A-fan C), k=1/2   S  (gemini k=1/2; luna weaker range)
    D-AB (death A-fan B)          S  (luna, k=(b+1)/(b+2) ≥ 2/3)
    D-BC (death B-fan C)          S  (luna, claims k=1 — to the wall)
    D-BA (death B-fan A), k=1/2   S  (luna, exactly k=1/2)
    B-AC (birth A-fan C)          S  (luna + gemini, both full range)
    B-AB (birth A-fan B)          S  (luna, full range)
    B-BC (birth B-fan C)          S  (luna, full range)
    τ ≠ 0 parametric              S  (luna closed form + sign; flagged
                                      the δ vs −δ convention correctly)

Every S row is *numerically* true on the float sweeps above and
*certified* true for a ≤ 14; what is missing is only director/skeptic
verification of the drafted parametric arguments.

## Outcome

- **VERIFIED (exact ring, 25 members; float corridor 400/400):** the
  closed-form margin system of §1 — every distinct gate margin of
  W(a,b), bijectively.
- **VERIFIED (proof, all a ≥ b ≥ 1, a ≥ 2):** T2, birth-side necessity:
  alive ⟹ θ < θ_birth. With 005: the window is contained in
  (θ_d, θ_birth) for every member. 006's birth law is thereby proven AS
  AN UPPER BOUND on the window; it stops being SPECULATION on that side.
- **VERIFIED (certified interval arithmetic, all 104 members a ≤ 14):**
  both canonical segments alive on t ∈ (0, 22.5/(a(b+1))], realizing
  every θ in (θ_d, θ_birth): **window(W(a,b)) = (γ_birth, γ_d) EXACTLY
  for every a ≤ 14** — the full window law, both laws, exact, for the
  certified range (suprema not attained on either side).
- **VERIFIED (proof, all (a,b)):** the five P-rows of §5, including
  N1 > 0 on the k = 1/2 death segment — 005's obstruction-1 hardest
  case — and its parametric N1′(0) sub-lemma.
- **EVIDENCE (float):** segment aliveness a ≤ 60 (1,829 members × both
  segments, zero failures), t\* = 45/(a(b+1)) wall sharpness, window
  coverage 3,200 pairs, 468k-sample birth-side adversarial scan.
- **SPECULATION:** the seven S-rows of §5 (worker-drafted parametric
  proofs, ranges as claimed, not yet independently re-derived). The
  full ALL-(a,b) window law is proven MODULO exactly these; the
  coverage conjecture inherits this conditionality (below).
- **NOT claimed:** the coverage conjecture as a theorem (see the honest
  chain below); anything about arcs as intervals per-triangle (007's C2
  pointwise caveat stands — our windows are swept pointwise in θ, one
  triangle per θ, which IS the form 006/013 need); the a = 1 column
  (unmeasured, queue); unstable orbits; any claim at the window
  endpoints themselves.

**The coverage chain as it now stands.** For every obtuse γ ∈ (90°,180°),
013 (skeptic-confirmed by 014) gives an explicit W(a,b), a ≥ b+1 ≥ 2,
with θ = 180−γ strictly inside (θ_d, θ_birth). By T2+005 that window is
exactly the alive range of θ; by §4 every θ in it is realized by an
explicit alive triangle. Every link in that chain is now either PROVEN
for all (a,b) (013, T2, 005-necessity) or proven for a ≤ 14 and reduced
to the seven S-lemmas beyond — so **"every obtuse angle has an alive
W-member" is now a theorem for every γ whose 013-witness has a ≤ 14**
(an explicit, dense-in-(90,180) but not full set of angles), **and for
ALL obtuse γ it is conditional on exactly the S-rows of §5** — no other
SPECULATION remains anywhere in the chain. (For calibration: this
record deliberately does NOT declare the conjecture resolved. The
S-lemmas are elementary and numerically bulletproof, but this lab does
not count a proof drafted by a worker and unverified by anyone as a
result. That verification is the queued next step, and it is small.)

## Why it failed / what survived

Nothing pursued this cycle failed; the honest ledger:

1. **What made queue 11 tractable was refusing to fight the fans
   per-gate.** 005 saw "3–5-term trig polynomials with the fan index
   linear in the angles" and stopped before hand-proving ~2a+2b margin
   bounds. The closed forms reveal the margins are ONE sinusoid each
   around two pivots, and the fan conditions are statements about
   sign-constant arithmetic progressions of angles — which is a
   pigeonhole, not an estimate. T2's proof is five steps and uses no
   interval arithmetic at all. The lesson for the mechanisms index:
   when a family of inequalities shares a lattice structure, prove the
   lattice statement, not the inequalities.
2. **ρ = 2 was load-bearing all along and nobody knew.** 005 chose the
   death-segment slope for a derivative-sign convenience; it happens to
   be exactly the slope that exits the alive region at the window
   midpoint (the PB wall at t = 45/(a(b+1))). That the certified
   segments "all certified with the same ρ = 2, t₀ = 1/4" was 005's
   observed uniformity; the wall explains it and gives the sharp
   constant.
3. **The swarm's real value was drafts-at-parity, not search.** 13
   returned drafts for ≈ $0.16; two (N1, T2) survived full director
   re-derivation — including the batch's two hardest — one needed a
   repair the draft itself made findable (its Step 5 named the right
   objects, imprecisely), and one "REFUTED" verdict was a correct
   catch of a sign-convention gap in the brief. None of the 13 was
   trusted unread; the S-rows stay SPECULATION precisely because
   re-derivation, not drafting, is the bottleneck this protocol prices.
4. **Reusable:** the margin system (`wmargins.py` — any future W(a,b)
   question starts from closed forms now); the two-segment window sweep
   pattern; the certified-trig segment certifier generalized to
   arbitrary closed-form margins (`wmargins_certify.py`); the
   progression-pigeonhole (T2 Steps 2–3) as a template for branch
   exclusion in ANY vertex-pivoting word family (005 lead 5's families
   included); the B-fan interleaving fact sin(aα+rβ) > 0 ∀r ∈ [−(b−1),b].

## Leads generated

1. **Skeptic pass on this record + the seven S-lemmas** (default stance:
   refute). The S-lemma drafts are committed verbatim under
   `$MATHLAB_OUT` is gitignored, so: re-run the two swarm jobs from the
   committed briefs (`wmargin_brief.md` + `wmargin_tasks.txt`, prompt
   hashes in the meta files) or re-derive from scratch — the latter is
   the actual task; the drafts are scaffolding, and each is ≤ 2 pages of
   elementary trig. Attack surface for the proven parts: T2 Step-5
   repair (the r-interleaving and the b = 1 edge case); N1 proof's
   Ratio 2 (sin u ≥ u·cos u usage range); the ring-bijection scope
   (25 members vs all (a,b) — formalize the endpoint bookkeeping in the
   4-variable ring the way plaw_general did I1–I4); the certifier's
   birth-segment rational range checks (θ → 90 for W(2,1)).
   **Anyone completing the seven verifications closes the coverage
   conjecture for all obtuse angles.** That is the single
   highest-value small task in the lab right now.
2. **Extend the certified sweep** past a = 14 (pure compute, ~an hour
   to a = 25 at current speed) — widens the unconditional angle set
   while the S-lemmas await verification.
3. **Formalize T2** (`docs/FORMALIZE.md` lane): five elementary steps,
   no interval arithmetic, mathlib-friendly objects (sin sign windows,
   arithmetic progressions). A natural next Lean target after the
   Laurent block, and T2 is now load-bearing for everything.
4. **The a = 1 column, still unmeasured** (005 lead 4): T2's proof never
   uses a ≥ 2 until the span bound 2 − 1/a ≥ 3/2; check whether the
   argument closes for a = 1 (W(1,b) windows would be (θ_d(1,b), 90+...)
   — degenerate?) and settle the column.
5. **Sup non-attainment at the birth edge**: T2 gives strict inequality,
   and the certified segments approach but never reach θ_birth. Prove
   zero-width touching at θ_birth itself (mirror of 004's C2 at γ_d) to
   finish the boundary bookkeeping.
6. **Cross-problem ripple**: the progression-pigeonhole (all-negative
   sinusoid over an AP forces one component, hence a span bound) is a
   general mechanism for "rotation-number window" arguments; scan
   almost-mathieu's Chambers-polynomial gap conditions for the same
   shape before inventing new machinery there.

## References

- `problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`
  (prefix maps, I1–I4, certified-trig layer, obstruction 1, leads 1+3).
- `problems/billiards-triangles/attempts/003-death-angle-laws.md` and
  `004-skeptic-review-of-003.md` (glide reduction; C2 sup-non-attainment).
- `problems/billiards-triangles/attempts/006-design-family-past-135.md`
  and `007-skeptic-review-of-006.md` (birth-law SPECULATION, window
  claim, pointwise caveat).
- `problems/billiards-triangles/attempts/013-coverage-diophantine-lemma.md`
  and `014-skeptic-review-of-013.md` (the coverage interface).
- `docs/SWARM.md` (worker protocol); swarm briefs committed at
  `explore/wmargin_brief.md`, `explore/wmargin_tasks.txt`; job metas
  (model, tokens, prompt hashes) under `$MATHLAB_OUT/swarm/wmargin-*/`.
- Worker drafts: gpt-5.6-luna (10), gemini-3.7-flash (3; T2 brief
  blocked by provider content filter, luna's copy used).
