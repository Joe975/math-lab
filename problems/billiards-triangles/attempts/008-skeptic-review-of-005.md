# 008 — Skeptic review of 005 (complete death-law theorem): adversarial verification

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed (read `prior-art.json`, attempts 003, 004, 005 in full,
  001/002 skimmed for conventions, the two tools under review
  `explore/plaw_general.py` and `explore/plaw_suffice.py` line by line, the
  28 `data/plaw_*.json` files, and tier-0 `unfold.py` as cross-reference)
- **Type:** skeptic review of `005-complete-death-law-theorem.md` (default
  stance: REFUTE). Every load-bearing layer re-derived by hand and/or
  re-computed by code written from scratch for this review; the only prior
  code in any verdict path is 004's `deathlaw_skeptic.py` stack (skeptic
  property, reused for its interval trig and its death sampler).
- **Outcome in one line:** 005 survives. The closed-form bridge, the
  4-variable formal ring proof of I1–I4, Lemmas L1/C/D, the extended case
  tree (a ≤ 2b+3 removed, k ≥ 1 branch), the sufficiency segments, and the
  five new deaths are all independently confirmed; one numeric correction
  (four of the five new exact-alive certificates are at γ_d − 1e-4, not
  the claimed 1e-6), one wording caveat; no load-bearing claim refuted.
- **Tools:** new `explore/psk_review.py` (stdlib-only, deterministic, seed
  20260738; subcommands `bridge`, `ring`, `lemma`, `threshold`, `casetree`,
  `segment`, `trigaudit`; ~5 min total). Reused: 004's
  `deathlaw_skeptic.py` (its Machin-pi interval stack and its adaptive
  death sampler — written for skeptics, no code shared with 003/005), and
  002's `skeptic_family.width` as a float cross-reference only. Nothing
  imports `plaw_general`/`plaw_suffice` except `trigaudit`, where
  `plaw_suffice` is the object under test.
- **Sources:** repo only; no external papers.

My data lives in `data/psk_*.json`. No file of 005 (or any prior attempt)
was modified.

Reproduce everything (repo root):

```
python3 problems/billiards-triangles/explore/psk_review.py ring      --out .../data/psk_ring.json
python3 .../psk_review.py bridge    --out .../data/psk_bridge.json
python3 .../psk_review.py lemma     --out .../data/psk_lemma.json
python3 .../psk_review.py threshold --out .../data/psk_threshold.json
python3 .../psk_review.py casetree --a 6 --b 1 --out .../data/psk_casetree_W6_1.json   # + (8,1),(9,1),(10,1),(12,1),(20,1),(12,2),(15,2),(30,3),(11,2)
python3 .../psk_review.py segment   --out .../data/psk_segment.json      # ~40 s
python3 .../psk_review.py trigaudit --trials 150 --out .../data/psk_trigaudit.json
python3 .../deathlaw_skeptic.py death --a 6 --b 1 --full --uniform 800 --iters 44 --out .../data/psk_death_W61.json   # + (8,2),(7,1),(6,3),(8,1)
```

## Claims attacked

Ordered as 005's own skeptic list (a)–(e) plus the two headline layers and
the scope audit.

### 1. The closed-form bridge u = R₀ ∘ Rot_A(2aα) ∘ Rot_B(−2bβ) (Steps 1–3). **CONFIRMED (hand re-derivation + exact off-torus recomposition)**

This is the one layer the formal check cannot protect. Hand re-derivation,
from the base geometry (A = 0, B = sin(α+β) real, C = sin β e^{iα}):

- *Base reflections:* side AB is the real axis through 0, so R₂ = conj.
  Side CA passes through 0 with unit direction e^{iα}, so R₁ =
  e^{2iα} conj. C − B = sin β e^{iα} − sin(α+β) = sin α·e^{i(180−β)}
  (expand sin(α+β); the real parts cancel to −sin α cos β), so side BC has
  direction-squared e^{−2iβ} and R₀(z) = B + e^{−2iβ} conj(z − B), using
  conj(B) = B. All three match 005's Step 2.
- *Pair collapse:* R₁R₂ = e^{2iα} conj∘conj = rotation about A by +2α;
  R₀R₂(z) = B + e^{−2iβ}(z − B) = rotation about B by −2β. ✓
- *Composition:* u(z) = R₀(e^{2iaα}(B + e^{−2ibβ}(z−B))) expands to
  μ conj z + w with μ = e^{−2iaα+2i(b−1)β} and
  w = B[1 + e^{−2i(aα+β)} − e^{−2iβ} − μ] — 005's Step 3 exactly, and
  δ = e^{i(−aα+(b−1)β)} gives δ² = μ. The word-order convention
  (M_k = R_{s₁}∘…∘R_{s_k}, gate k = M_{k−1}(side s_k)) is the same
  composed-map unfolding 004 already tied to the harness corridor.
- *Gate 2a+2's C-endpoint:* letters 1..2a+1 are 0(12)^a, so
  M_{2a+1} = R₀(R₁R₂)^a = R₀ Rot_A(2aα) and the C-image endpoint of gate
  2a+2 (mirror side 0 = BC, endpoints B, C) is v₂ = R₀(Rot_A(2aα) C) —
  005's I4 anchor is the correct geometric object, and in
  `deathlaw_symbolic`'s gate list it is `gates[2a+1][1]` (SIDE_ENDS[0] =
  (1,2), v-endpoint = C-image), which is what both `plaw_general.py
  specialize` and `plaw_suffice.py` use. ✓
- *Glide action:* p(u(z)) = Im(conj(δ)(δ² conj z + w)) = −p(z) + p(w),
  p(τ) = 0, m = p(w)/2, hence p∘u = 2m − p. ✓

Independent computation (`bridge`, `data/psk_bridge.json`): my own
**off-torus pair algebra** (004's device, re-implemented from scratch:
track (value, star-value) at rational points x₀, y₀ off the unit torus;
star = the involution x→1/x, y→1/y, i→−i is a ring automorphism equal to
complex conjugation on the torus, so +, −, ×, ÷, conj all evaluate
exactly in Q(i)). The half word is composed by DIRECT nested application
of my base reflection maps — no closed form anywhere — and compared
exactly with 005's claimed μ, w, v₂ (specialized P→x₀^a, Q→y₀^b) plus
μ = δ², Im(conj(δ)τ) = 0 and the glide action p(u(z₁)) + p(z₁) − 2m = 0
at a generic point. **40 members — the 23 of 005's specialize list plus
(1,1), (40,1), (40,40), (37,13), (25,24) and 12 random pairs with
a ≤ 40 — at 2 random rational points each: every comparison is an exact
match.** (Each match is an exact rational identity between a
(2a+2b+1)-step reflection composition and the closed form; a wrong
bridge cannot survive this at random points.)

Float check: 005's own `floatcheck` (400 random (a,b,α,β,z), a ≤ 40)
reproduced, worst error 1.4e-14 — same number as the record.

### 2. The 4-variable formal ring proof and the specialization argument. **CONFIRMED (independent ring; the specialization direction is sound)**

- *Soundness of specialization:* the map P ↦ x^{2a}, Q ↦ y^{2b} (into the
  half-angle ring), followed by evaluation on the torus, is a composition
  of ring homomorphisms that also commutes with conj (exponent negation +
  coefficient conjugation on both sides) and hence with Im. A polynomial
  that is ZERO in Q(i)[x^±,y^±,P^±,Q^±] maps to zero under any such
  homomorphism — the dangerous direction ("formal identity fails to
  specialize") cannot occur: the formal ring has strictly fewer relations
  than the geometry, so formal zero ⟹ geometric zero, for every integer
  a, b and all angles. What COULD fail is the bridge (a wrong closed form
  being formally consistent) — that is claim 1, tested independently
  above. The half-integer subtlety is absent: all eight identities live in
  whole-angle monomials, and the specialize target ring's keys (m,n) are
  even multiples throughout.
- *My own formal ring* (`ring`, `data/psk_ring.json`): re-implemented
  from scratch, and — unlike 005, which hard-codes the closed forms — I
  RE-DERIVE u inside the ring by composing affine maps: R₀, R₁, R₂ as
  antilinear maps, the block products checked to be the two rotations
  (`rotA_block`, `rotB_block` true), the a-th/b-th powers introduced via
  the fixed-point identities (rot_B fixes B: checked as a ring identity;
  z ↦ λ(z−fix)+fix powers to λ^n(z−fix)+fix, the only lattice steps being
  (x²)^a = P², (y^{−2})^b = Q^{−2}). The derived μ, w match 005's closed
  forms as ring identities, and **I1, I2, I3, I4, D1, D2, μ = δ²,
  Im(conj(δ)τ) = 0 are all exact zero polynomials in my ring.** 14/14
  checks pass.
- *I4's geometric meaning:* covered exactly by the bridge's v₂ test
  (claim 1) plus `ring`'s I4 = p(v₂) − m identity; also note I1 = −(S+T),
  I4 = S − T with the same two products — verified structurally in my
  ring by building S, T once and using them in both.

### 3. Lemma L1, Lemma C, Lemma D. **CONFIRMED (full hand re-derivation; adversarial numerics clean)**

Hand re-derivation of the entire chain (radians):

- *L1:* d/dc[c cot(cs)] = [sin(2cs) − 2cs]/(2 sin²(cs)) < 0 for
  cs ∈ (0, π/2] since sin x < x for x > 0 and sin(cs) ≠ 0. Strict. ✓
- *Lemma C:* (log Φ)′ = 2/sin 2v − (1/a)cot((π/2−v)/a) − (b/a)cot(bv/a);
  the split 2/sin 2v = tan v + cot v is the identity
  1/(sin v cos v) = tan v + cot v ✓. First bracket
  [cot v − (b/a)cot((b/a)v)] ≤ 0 by L1 at s = v ∈ (0, v₀] ⊆ (0, π/4]
  (c = b/a ≤ 1 vs c = 1; equality iff a = b) ✓; second bracket
  [tan v − (1/a)cot((π/2−v)/a)] < 0 by L1 at s = π/2 − v ∈
  [π/2 − π/4, π/2) ⊂ (0, π/2), using cot(π/2−v) = tan v, strict since
  1/a ≤ 1/2 < 1 (this is where a ≥ 2 enters) ✓. So Φ strictly decreases
  on (0, v₀]. Endpoint: π/2 = (b+1)v₀ gives (π/2−v₀)/a = bv₀/a exactly,
  so Φ(v₀) = tan v₀ ✓. Factorization: cos v cos v₀ sin(bv/a)·Φ(v)
  = cos v₀ sin v sin((π/2−v)/a) and cos v cos v₀ sin(bv/a)·Φ(v₀)
  = cos v sin v₀ sin(bv/a), so H2 = cos v cos v₀ sin(bv/a)[Φ(v)−Φ(v₀)]
  > 0 on (0, v₀), with the prefactor positive because v < v₀ ≤ π/4
  (cos v > 0) and bv/a ≤ v < π/4 (b ≤ a enters here). H2(v₀) = 0 by
  substitution ✓. Domains: all cot arguments stay in (0, π/2] as
  required by L1; no degree/radian slip (005's subsection is consistently
  radian). The b = 1 path is not degenerate for Lemma C itself (v₀ = 45°);
  a = b = 1 fails only the second bracket's strictness — correctly
  excluded (H2 ≡ 0 there, as 005 says).
- *Lemma D:* (log D_b)′ = b cot(bx) − cot x < 0 on (0, π/(2b)] by L1 at
  s = bx (c = 1/b vs 1, strict for b ≥ 2) ✓. This retires 004's
  "Lemma D taken as classical" residual with a genuinely two-line proof.

Adversarial numerics (`lemma`, `data/psk_lemma.json`, my own radian
code): 14,183 checks — L1 monotonicity at 4000 random (c₁,c₂,s);
F(v) < 0 for a ∈ {2,3,4,7,20,100,1000,10⁶}, b ∈ {1,2,a/2,a−1,a}, v down
to 1e-9 and up to (π/2)(1−1e-12), incl. v near v₀ from both sides; H2 > 0
on (0,v₀) via the cancellation-free factored form plus agreement of the
factored and direct forms (identity check); Lemma D as tan(bx) > b tan x
for b up to 10⁵. **Zero violations.** (My first run reported 5 Lemma-D
"violations"; all were artifacts of MY test — x outside (0, π/(2b)] for
b = 10⁵, and float cancellation at x = 1e-9 — recorded in "Kill attempts
that failed" below.)

### 4. The extended case tree (a ≤ 2b+3 removed; the k ≥ 1 branch). **CONFIRMED (every branch re-derived by hand + 447k-point scan on a > 2b+3 members)**

Hand re-enumeration of all branches (Case I/II × residue), for general
a ≥ b ≥ 1, a ≥ 2, every k ≥ 0. Highlights, with every inequality
re-checked:

- *Exhaustiveness and boundary self-exclusion:* as in 004 (residues at
  multiples of 90 in aα or (b+1)β force N2 = 0 or N3 = 0; sin θ > 0
  always since θ < 180). ✓
- *Case II, r ∈ (0,90) ∪ (270,360)* — the A4′ fix: cos(aα) > 0 on BOTH
  intervals, so N2 < 0 forces (b+1)β > 180 and θ > β > 180/(b+1) ≥ θ_d
  ⟺ a ≥ b, never touching α. k-free. ✓ The other Case II branches:
  r ∈ (90,180) has slack exactly 90/(a(b+1)) ✓; r ∈ (180,270): the >360
  escape needs 3a ≥ b ✓ (always true), the (0,90)-residue path forces
  cos(α+(b+1)β) < 0 with α+(b+1)β ∈ (0,270), hence > 90, and α > 180/a
  beats θ_d by 90b/(a(b+1)) > 0 ✓; both k-free as claimed.
- *Case I, r ∈ (270,360)* — the A4 fix: cos(aα) > 0, sin(aα) < 0 force
  (b+1)β mod 360 ∈ (90,180) via N2, N3 > 0, so (b+1)β > 90; with
  α > 270/a: 270/a + 90/(b+1) − θ_d = (180b+270)/(a(b+1)) > 0 —
  unconditional, restriction gone. ✓ r ∈ (90,180) and (180,270) as in
  003, k-free ✓.
- *Case I, r ∈ (0,90) — the main case.* cos(aα) = sin v, sin(aα) = cos v
  for aα = 360k + 90 − v ✓; (b+1)β mod 360 ∈ (0,90) forced, > 360 escape
  needs 3a ≥ b ✓; G = sin v sin α sin(bβ) − cos v sin β sin w with
  w = α + (b+1)β − 90 ∈ (−90, 180) ✓ (α < 180, (b+1)β ∈ (0,90)).
  **Threshold algebra:** aw − bv = a(b+1)α + a(b+1)β − 90(a+b) − 360kb
  = a(b+1)(θ − θ_d) − 360kb — re-derived by hand and verified as an
  exact Fraction identity at 2000 random rational (α, β, a ≤ 60, b ≤ a,
  k ≤ 5) (`threshold`, `data/psk_threshold.json`). ✓
  **k ≥ 1:** θ ≤ θ_d gives aw ≤ bv − 360kb = b(v − 360k) < 0 since
  v < 90 < 360k, so w < 0 strictly (005's "w ∈ (−90,0], sin w ≤ 0" is if
  anything conservative), so G ≥ sin v sin α sin(bβ) > 0, contradicting
  Case I's G < 0. No lemmas, no bound on a. ✓
  **k = 0:** 003's argument verbatim; I re-checked the w ≤ 0 sub-case,
  Regime 2 (needs w < α < 90/a ≤ 45, i.e. a ≥ 2, and strict Lemma D for
  b ≥ 2 with β < 90/(b+1) strict), Regime 1 (needs w ≤ bv/a ≤ v, i.e.
  b ≤ a, then Lemma C strict + Lemma D; the b = 1 path closes via strict
  Lemma C with v₀ = 45, cot v₀ = 1), and the θ = θ_d equality case (still
  a strict contradiction). ✓ Preconditions a ≥ 2, b ≤ a are used exactly
  where 005 says and nowhere else; b > a is nowhere claimed.

Adversarial search (`casetree`, `data/psk_casetree_W*_*.json`): my own
design and seed (residue boundaries at ALL multiples 90k/a for every
reachable k, (b+1)β boundaries, death-corner accumulation, θ at θ_d
exactly and θ_d(1−10^{−j}) down to 10^{−12}, random fill), 10 members
chosen to stress a > 2b+3 and the k ≥ 1 region — (6,1), (8,1), (9,1),
(10,1), (12,1), (20,1), (12,2), (15,2), (30,3), (11,2); aα reaches 945°
for (20,1) and 742° for (30,3), so k = 1, 2 are genuinely exercised —
**447,000+ points, zero Case I/II hits, zero sub-margin near-hits, and
the positive control just above θ_d fires for all 10 members.**

### 5. Sufficiency: segments, certifier, converse reduction. **CONFIRMED (hand re-derivations + independent interval certification at γ_d − 9.3e-10)**

- *Segment algebra by hand:* θ(t) = 90/a + 90(a−1)/(a(b+1)) + t and the
  endpoint identity θ(0) = [90(b+1) + 90(a−1)]/(a(b+1)) = θ_d is exact;
  γ(t) = γ_d − t sweeps [γ_d − 1/4, γ_d) monotonically for t ∈ (0, 1/4]
  ✓ (with ρ = 2, the record's general (ρ−1)t reduces to t ✓). The corner
  identities aα_c = 90 and α_c + (b+1)β_c = 90 are exact ✓, giving
  N1(t) = cos(at) sin β(t) sin(ŵt) − sin(at) sin α(t) sin(bβ(t)),
  ŵ = (b+1)ρ − 1 ✓ (cos(90−at) = sin(at); cos(90+ŵt) = −sin ŵt), and
  N2, N4 all-factors-positive under the rational range checks, which I
  verified use exactly a·t₀ ≤ 45, (b+1)β(t) ∈ (90−90/a, 180),
  θ(t) < 90, bβ(t) < 180, ŵt₀ < 180 — each sound for the members run
  (90 − 90/a > 0 needs a ≥ 2, satisfied by all 20). ✓
- *N1′ re-differentiated by hand:* six terms; matches `n1p_iv` exactly
  including the dropped global π/180 (every term of dN1/dt carries exactly
  one such factor — sign-safe) and N1′(0) = ŵ sin β_c − a sin α_c
  sin(bβ_c) ✓. The N1 logic (N1(0) = 0 exactly; N1′ > 0 certified on
  [0, t₁] ⟹ N1 > 0 on (0, t₁]; bisection on [t₁, t₀]) is sound. ✓
- *Converse reduction re-derived:* if lo_k < m < hi_k strictly for all
  half gates, then L = max lo_k < m < min hi_k = H, and the full corridor
  [max(L, 2m−H), min(H, 2m−L)] has both upper candidates > m and both
  lower candidates < m — positive width. Uses gate n′+k = u(gate k)
  (functoriality — an identity of the composed-map construction),
  p∘u = 2m − p (claim 1), and ∩(2m − I_k) = 2m − ∩I_k ✓. τ ≠ 0 is needed
  for the projection direction to exist; the certifier certifies
  Re(conj(δ)τ) of constant nonzero sign along the whole segment, and
  Im(conj(δ)τ) = 0 in the ring, so τ ≠ 0 throughout ✓. Scaling to the
  harness's τ-normal projection is the nonzero-real-rescale argument 004
  already audited. ✓
- *Certifier audit* (`plaw_suffice.py` read line by line): the trig layer
  is sound — exact mod-360 midpoint reduction into [−180, 180) (Fraction
  floor-division), PI_LO/PI_HI radian bracketing with correct sign
  handling, dyadic pre-round with its error added to the slack (plus a
  conservative 1/DY), alternating-series remainders valid at 14 terms for
  |x| ≤ 4 > π (the term ratio x²/(30·31) < 1 makes the first-omitted-term
  bound legitimate), 1-Lipschitz widening in radians for the interval
  half-width (valid across reduction since |sin x − sin y| ≤ |x−y|
  globally). Every evaluation is a single non-iterated interval
  evaluation. `poly_terms`' lexicographic (m,n) > (0,0) selection counts
  each conjugate pair exactly once, and the (0,0) coefficient is real by
  conj-invariance (asserted). Gate straddling requires certified opposite
  signs at every gate's two endpoints; the three corner-vanishing margins
  get signs on the open segment only, which is all the t ∈ (0, t₀] claim
  needs. Failure paths are loud. One nit: `n1_iv` assigns an unused
  variable `t`.
- *Trig audit* (`trigaudit`, `data/psk_trigaudit.json`): plaw_suffice's
  sin_deg/cos_deg enclosures at 300 random rational degrees up to ±10⁷
  (raw argument to plaw, exactly-reduced argument to my reference — 004's
  independent Machin stack, width ~1e-84) — my tight enclosure sits
  inside plaw's in every case, certifying true-value containment; same
  for 148 random degree INTERVALS at 4 interior rational points each.
  **596 checks, zero violations; worst plaw point-enclosure width
  5e-16.** (My first run "FAILED" — entirely my audit's fault: I fed my
  reduction-free reference stack raw 3000°+ arguments, far outside its
  series' validity. See kill attempts.)
- *Independent alive certification* (`segment`, `data/psk_segment.json`)
  — the strongest new evidence: for (2,1), (6,3), (6,1), (8,2), (8,1),
  (7,1), (12,12) at t ∈ {1/4, 1/8, 1/97, 1/1024, 2^{−20}, 2^{−30}} I
  certify positive corridor width by MY OWN interval unfolding: the FULL
  2n-gate corridor (no glide reduction anywhere), gates composed directly
  from base reflections in interval complex arithmetic (004's certified
  trig, outward dyadic rounding at 2^{−320}), projection onto the normal
  of the full-word translation τ, with the linear part certified to
  enclose 1 and τ certified nonzero. **All 42 points certified ALIVE,
  including t = 2^{−30}, i.e. γ within 9.3e-10 of γ_d — with certified
  width lower bounds ~1e-10 clear of the interval slop by ~20 orders.**
  This independently confirms the segment construction where it matters
  most (the 3-fold-degenerate corner approach) and simultaneously
  exercises the glide reduction's converse (my full corridor is positive
  exactly where the certifier's half-gate criterion says it must be).
  Float cross-checks via 002's corridor agree in sign at all 42 points.

### 6. The five new float deaths. **CONFIRMED (independent re-measurement; windows located)**

004's sampler (`deathlaw_skeptic.py death --full`, my parameters, both
arc halves), `data/psk_death_W*.json`:

| word | predicted | 005 measured − pred | mine − pred | last-alive apex | window half near death |
|---|---|---|---|---|---|
| W(6,3) | 585/4 = 146.25   | −3.7e-12 | −3.7e-12 | (15, 18.75)      | **α < β only** (best α ≥ β width −6.1e-6) |
| W(6,1) | 255/2 = 127.5    | −5.3e-12 | −5.4e-12 | (15, 37.5)       | **α < β only** (−5.2e-6) |
| W(8,2) | 285/2 = 142.5    | −3.6e-12 | −3.8e-12 | (11.25, 26.25)   | **α < β only** (−7.4e-6) |
| W(8,1) | 1035/8 = 129.375 | −4.9e-12 | −4.9e-12 | (11.25, 39.375)  | **α < β only** (−4.9e-6) |
| W(7,1) | 900/7 = 128.5714 | −5.1e-12 | −5.1e-12 | (12.857, 38.571) | **α < β only** (−4.0e-8) |

All five die at γ_d at the predicted corner (90/a, 90(a−1)/(a(b+1))).
The near-death windows of all five b ≤ a−2 members live exclusively on
the mirror half α < β — the half 001/002 never scanned — extending 004
§5's finding to the new members (005's skeptic item (e)). The ~5e-12
agreements share the one-sided float-tolerance design bias 004
documented; the exact content is carried by the theorem + segments, as
005 itself says.

### 7. Scope honesty. **CONFIRMED with one correction (C1) and one caveat**

- "Sup not attained" is exactly justified: necessity gives no positive
  width at γ ≥ γ_d (004's C2 wording adopted verbatim — a zero-width
  touching corridor at γ_d stays possible and is explicitly not
  excluded); the segments give alive γ filling [γ_d − 1/4, γ_d); together
  sup = γ_d, not attained. ✓
- The NOT-claimed list is complete and correct: parametric sufficiency
  open (per-member certificates only), a = 1 column open, birth side
  open, other words/unstable orbits untouched, measured deaths float.
  The theorem's hypotheses (a ≥ b ≥ 1, a ≥ 2) match where they are used;
  b > a is nowhere claimed (identities hold for all a, b ≥ 1 — correct,
  the formal proof does not need b ≤ a). ✓
- Index entry (`prior-art.json` 005): one_line, range, leak_terms, gaps
  all match the record and the data files (casetree 13 members ×
  ~120k = 1.56M ✓; lemmac grid 914,500 ✓; specialize 23/23 ✓; suffice
  20/20 with ρ = 2, t₀ = 1/4 ✓; float cross-check 40/40 positive ✓;
  measure tables match the JSON to the last digit ✓). ✓
- **But the record's exact-certificate sentence is wrong — C1 below.**

## Refutations found

No load-bearing claim is refuted. One numeric correction, one caveat:

- **C1 (wrong constant, same genus as 004's C1).** 005, sufficiency
  section: "all five: exact width > 0, orbit verified, γ certified within
  **1e-6** of γ_d." False for four of the five: the certified lower
  bounds in the data files are γ_d − **1e-4** for W(6,3)
  (1462499/10000), W(8,2) (1424999/10000), W(8,1) (1293749/10000), and
  W(7,1) (8999993/70000); only W(6,1) (127499999/1000000) reaches
  γ_d − 1e-6. Cause (read from `deathlaw_exact.py cmd_alive`): the
  `--dg 1e-6` request sets the float target, but the certified bound is
  the largest 10^{−k} that CERTIFIES after apex snapping, and for these
  four the k = 6 certification fails, leaving k = 4. Nothing downstream
  is damaged — the segment certificates (independently confirmed here to
  γ_d − 9.3e-10) supersede these brackets entirely — but the sentence
  misstates its own data files by two orders of magnitude.
- **Caveat (not a refutation).** For the five NEW members the necessity
  side is the parametric theorem alone — no per-member ring obligations
  (003's `deathlaw_prove`) were run for them, unlike the fifteen 003
  members. That is by design (the general identities specialize, and
  were re-proven here), but a reader tallying per-member machine
  certificates should know the new five lean entirely on the general
  proofs.

## Claims that survive

| # | 005 claim | Verdict |
|---|-----------|---------|
| 1 | Closed-form bridge u = R₀ Rot_A(2aα) Rot_B(−2bβ), gate-(2a+2) anchor v₂ | **CONFIRMED** — hand re-derivation of Steps 1–3; exact off-torus recomposition, 40 members incl. a = 40, every comparison an exact match |
| 2 | I1–I4, D1, D2, glide facts formal in Q(i)[x,y,P,Q]; valid for ALL (a,b) by specialization | **CONFIRMED** — independent ring, u re-derived by symbolic map composition, 14/14 zero polynomials; specialization direction proven sound (ring hom; formal ring has fewer relations) |
| 3 | Lemma L1; Lemma C for a ≥ 2, b ≤ a; Lemma D for b ≥ 1 | **CONFIRMED** — full hand re-derivation incl. domains, strictness, endpoint identity, factorization; 14,183-point adversarial numerics clean to a = 10⁶ |
| 4 | Extended case tree: a ≤ 2b+3 removed; k ≥ 1 branch; threshold algebra | **CONFIRMED** — all branches re-derived by hand; threshold identity exact at 2000 rational points; 447k-point scan on 10 a > 2b+3 members (k up to 2 reachable), zero hits, controls fire |
| 5 | death(W(a,b)) = γ_d exactly (sup not attained), 20 members, via certified segments | **CONFIRMED** — segment algebra + N1′ re-derived by hand; certifier and trig layer audited sound; converse reduction re-proven; independent interval certification of alive at 42 segment points down to γ_d − 9.3e-10 |
| 6 | Five new deaths at γ_d (~5e-12), predicted corners | **CONFIRMED** — independent sampler agrees to ~1e-13 of 005's values; all five windows located on the mirror half α < β |
| 7 | Scope statements and index entry match certificates | **CONFIRMED with C1** (four of five new exact certificates are at γ_d − 1e-4, not 1e-6) |

**Kill attempts that failed,** for the record: (i) my Lemma-D stress
initially reported 5 violations — all artifacts of my own test
(out-of-domain x for b = 10⁵, float cancellation at x = 1e-9); the fixed
domain-respecting test is clean. (ii) My first trig audit reported
plaw_suffice enclosure violations — all artifacts of my own reference
stack, which has no mod-360 reduction and whose 30-term series is
invalid past ~660°; with exact pre-reduction (sin/cos are 360-periodic,
so the same value is tested) plaw_suffice contains the true value in all
596 checks. Two lessons about MY tools, zero about 005's. (iii)
Boundary- and k-region-targeted sign search for a dropped branch of the
extended tree: nothing, at 447k points with near-hit logging (zero even
below margin). (iv) Hunting for a failure of the segment construction at
extreme t: my independent full-corridor interval certification stays
positive down to t = 2^{−30}, with ~20 orders of margin over the
interval slop.

**Net assessment:** 005's headline — the death law completed on both
sides, necessity parametric for all a ≥ b ≥ 1, a ≥ 2, and
death(W(a,b)) = γ_d(a,b) exactly with the sup not attained for 20
members — stands as scoped, now with every load-bearing layer re-derived
or re-computed independently. The status upgrade to VERIFIED is
justified under the lab's definitions (formal-ring identities = proof;
Lemma proofs elementary and re-checked; per-member sufficiency
certified; parametric sufficiency correctly left open).

## Residual risk

- **Shared normalization.** My bridge/ring/segment work all live in the
  same circumdiameter normalization introduced in 003 (as do 004's
  checks). 004's apex-coordinate cross-checks (harness point boxes, exact
  width matches) bound the risk that the normalization itself is wrong,
  and my `segment` certification is cross-checked in sign against 002's
  apex-coordinate float corridor at all 42 points — but a conceptual
  error common to every angle-domain formulation would have to be caught
  by those apex-side ties; nothing here adds to them.
- **Sampling, not proof, on the extended tree's reachable set.** The
  hand proof covers all k; the scans (mine + 005's) only sample k ≤ 2.
  The k ≥ 3 region (roughly a ≳ 30b members) is proof-covered but has
  never been numerically exercised.
- **The five new members' necessity is parametric-only** (see caveat).
- **Shared-design bias in float death numbers** persists exactly as 004
  described; two same-design instruments agreeing at ~1e-12 is not two
  independent measurements of that digit.
- **Not covered here:** re-running the 20-member `plaw_suffice all`
  certification end to end (I reproduced selftest + one member, audited
  the code, and replaced the full sweep with my own 42-point independent
  certification on 7 members); the orbit realizations inside the five
  `plaw_exact_alive_*` files (002's simulator, validated there; I
  verified the alive criterion, not the trajectories); Lemma C's wider
  F < 0 claim on all of (0, 90) beyond v₀ (not load-bearing; spot checks
  only); the 3-fold (rather than higher) degeneracy count at the corner
  (not load-bearing for the certificates — the certifier does not assume
  it — only for 005's obstruction narrative).

## References

- `problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`
  (under review); `003-death-angle-laws.md` (the base theorem);
  `004-skeptic-review-of-003.md` (skeptic exemplar; its C1/C2, its
  composed-map unfolding and interval stack); `001-…`/`002-…` (skimmed,
  conventions).
- Code under review: `problems/billiards-triangles/explore/plaw_general.py`,
  `explore/plaw_suffice.py`; data `data/plaw_*.json` (28 files);
  supporting 003 tools `deathlaw_symbolic.py`, `deathlaw_exact.py` (read
  for C1's cause).
- Tier-0: `harness/billiards-triangles/unfold.py` (read as
  cross-reference only; `plaw_suffice` imports its Iv/pi — audited here).
- New code and data: `problems/billiards-triangles/explore/psk_review.py`;
  `problems/billiards-triangles/data/psk_{ring,bridge,lemma,threshold,segment,trigaudit}.json`,
  `data/psk_casetree_W{6_1,8_1,9_1,10_1,12_1,20_1,12_2,15_2,30_3,11_2}.json`,
  `data/psk_death_W{61,82,71,63,81}.json`.
- No external papers consulted.
