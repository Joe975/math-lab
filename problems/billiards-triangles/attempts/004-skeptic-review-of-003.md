# 004 — Skeptic review of 003 (death-angle laws): adversarial verification

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed (read `prior-art.json`, attempts 001, 002, 003 in full,
  the four `deathlaw_*.py` tools under review, the tier-0 harness, and the
  `STATUS.md` affine-forms insight)
- **Type:** skeptic review of `003-death-angle-laws.md` (default stance:
  REFUTE). Every load-bearing layer re-derived by hand and/or re-computed by
  code written from scratch for this review; nothing imported from the four
  tools under review.
- **Outcome in one line:** 003 survives. The glide reduction, the mod-360
  case tree, the identities I1-I3, Lemma C, and all four out-of-sample death
  angles are independently confirmed; two minor wording-level corrections
  (one wrong exponent in a certificate citation, one "empty corridor"
  overstatement), no load-bearing claim refuted.
- **Tools:** new `explore/deathlaw_skeptic.py` (stdlib-only, deterministic,
  seeded 20260731; ~10 min total, dominated by the Lemma C re-certification
  under my slow exact pi). Independence choices: unfolding by composed
  isometry maps applied to base-triangle sides (003 reflects vertex chains);
  my own Laurent-ring arithmetic; a second, representation-free check by
  exact rational evaluation OFF the unit torus in a formal-conjugate pair
  algebra; my own interval stack (Machin pi with alternating-series
  remainders, own Taylor enclosures) — `unfold.py` intervals are not used in
  any verdict path. 002's `skeptic_family.py` is imported only inside the
  selftest as a cross-reference.
- **Sources:** repo only; no external papers.

My data lives in `data/dlsk_*.json` (never overwriting 003's files). No file
of 003 was modified.

Reproduce everything (repo root):

```
python3 problems/billiards-triangles/explore/deathlaw_skeptic.py selftest
python3 .../deathlaw_skeptic.py identities --out .../data/dlsk_identities.json
python3 .../deathlaw_skeptic.py gridproof --a 2 --b 1 --out .../data/dlsk_gridproof_W21.json
python3 .../deathlaw_skeptic.py reduction --out .../data/dlsk_reduction.json
python3 .../deathlaw_skeptic.py casetree --out .../data/dlsk_casetree.json
python3 .../deathlaw_skeptic.py alivecheck --out .../data/dlsk_alivecheck.json
python3 .../deathlaw_skeptic.py death --a 10 --b 9  --out .../data/dlsk_death_W109.json   # + (12,12); (5,3),(4,2) with --full; (2,1),(3,3) calibration
python3 .../deathlaw_skeptic.py deadhunt --a 5 --b 3 --out .../data/dlsk_deadhunt_W53.json  # + (4,2),(10,9),(12,12)
python3 .../deathlaw_skeptic.py lemmac --a 12 --b 12 --out .../data/dlsk_lemmac_W1212.json  # + (2,1),(5,3),(10,9)
```

## Claims attacked

Ordered as 003's own skeptic list plus the scope audit.

### 1. The glide-to-corridor reduction. **CONFIRMED (re-derived by hand)**

The load-bearing bridge is: *corridor has positive width ⟹ the axis offset
m lies strictly inside the projection interval of every half-word gate*.
Hand re-derivation, from 001's corridor definition (intersection over all
2n' gates of the endpoint-projection intervals onto the normal of tau):

- *Gate functoriality:* the unfolding of u·u traverses u's gates and then
  the composed-isometry images u(gate k); so gate n'+k = u(gate k). I
  verified this structurally with my composed-map unfolding (gate k =
  M_{k-1}(base side), M_k = M_{k-1} ∘ R_{s_k}) — a different construction
  from 003's reflected-vertex chains that makes the functoriality an
  identity of the construction, and whose gates agree exactly with the
  harness corridor (selftest: exact corridor sign equals `unfold.certify`
  at point boxes; float corridor equals 002's independent evaluator on 120
  samples).
- *Glide action:* with u(z) = mu·conj(z) + w, mu = delta², |delta| = 1, and
  p(z) = Im(conj(delta) z), one line gives p(u(z)) = -p(z) + p(w), and
  tau ∥ delta gives m = p(w)/2 = Im(conj(delta)(w - tau/2))/2 — 003's
  formula. I verified p(u(P)) + p(P) - 2m = 0 as an exact ring identity on
  all base vertices for all 15 members (this is my own decomposition, not
  003's code path).
- *Strictness:* positive corridor width ⟹ the corridor, an interval of
  positive length, is contained in I_k ∩ (2m - I_k) for every half gate
  I_k = [lo_k, hi_k]; positive-length overlap means
  max(lo_k, 2m - hi_k) < min(hi_k, 2m - lo_k), whose two cross terms give
  lo_k < m < hi_k **strictly**. The boundary case is handled correctly: if
  m equals an endpoint, the overlap is at most a point, i.e. NOT positive
  width — so no alive triangle is wrongly excluded. The open-vs-closed
  bookkeeping is airtight.
- *Scaling:* the harness projects onto (-tau_y, tau_x); with tau a real
  multiple of delta this rescales all projections and m by the same nonzero
  real (possibly negative), and "min < m < max" is invariant under both.
  tau = 0 can only produce UNKNOWN/dead in both harness and float layers,
  never positive width, so the hypothesis excludes it.

Numerically: on 1818 triangles (300 random per member across six members
plus targeted near-death apexes), full-corridor width / |tau| equals the
reduced width min(hi, 2m-lo) - max(lo, 2m-hi) to 1e-8 relative, with **0
mismatches**, and every alive sample has m strictly inside every half-gate
interval (`dlsk_reduction.json`).

### 2. The mod-360 case tree. **CONFIRMED (re-enumerated by hand + 1.9M-point search)**

There is no code tree to diff against — the tree lives only in 003's prose;
the machine-certified inputs are I1-I3 and Lemma C. So I re-derived every
branch by hand (for general (a,b) in scope, which covers both the b = a and
b < a members) and then attacked the implication numerically.

Hand re-derivation highlights (full detail suppressed where it reproduces
003 verbatim; every inequality was re-checked):

- *Exhaustiveness of Case I/II:* m strictly inside gates 1 = [B,C] and
  2 = [C,A1] splits on the sign of m - p(C); m = p(C) is impossible since
  p(C) is an endpoint of gate 1's interval. So Case I (N1,N2,N3 > 0) and
  Case II (all < 0) are the only possibilities. ✓
- *Branch boundaries are self-excluding:* a·alpha ∈ {90, 180, 270, 360}
  forces N2 = 0 or N3 = 0, contradicting either case; likewise
  (b+1)·beta ∈ {90, 180, ...} — so the open-interval branches cover
  everything and no boundary case is silently dropped. ✓
- *The threshold algebra:* a·w - b·v = a(b+1)·theta - 90(a+b), so
  a·w < b·v ⟺ theta < theta_d, exactly as claimed. ✓
- *Each branch inequality:* A1' needs a ≥ b (180/(b+1) ≥ theta_d ⟺ a ≥ b);
  A2' gives theta > 90/a + 90/(b+1) > theta_d with strict slack
  90/(a(b+1)); A3' (the subtle one): residue of (b+1)beta in (0,90),
  the >360 escape needs 3a ≥ b, and the G-sign step gives
  theta > (b·alpha + 90)/(b+1) with alpha > 180/a, which beats theta_d by
  90b > 0 — re-derived including the cos < 0 window (the sum
  alpha + (b+1)beta lies in (0, 270), so cos < 0 really means > 90). ✓
  A4'/A4 use 270/a ≥ theta_d ⟺ **a ≤ 2b+3 — the scope precondition is
  load-bearing exactly here and nowhere else**. ✓ Case I A1: the >270
  escape needs 2a ≥ b; Regime 2 uses alpha < 90/a ≤ 45 (**this is where
  a ≥ 2 is used**), sin(w) < sin(alpha), and Lemma D's strict
  monotonicity (beta < 90/(b+1) strictly since (b+1)beta ∈ (0,90) open);
  Regime 1 uses w ≤ bv/a ≤ v (**b ≤ a used here**) and Lemma C; the b = 1
  degenerate D_1 ≡ 1 path closes via strict Lemma C as stated. Equality
  theta = theta_d also lands in a strict contradiction (via strict Lemma
  C / strict Lemma D), so the conclusion theta > theta_d is genuinely
  strict. ✓ Lemma D itself re-derived (tan addition induction) and
  numerically spot-checked.
- *Adversarial search* (`casetree`, `dlsk_casetree.json`): 1,938,477
  seeded points across all 15 members with theta ≤ theta_d — half spread,
  half concentrated at theta_d·(1 - 10^-k) and at the branch boundaries
  a·alpha ≈ 90/180/270/360 ± eps and the death corner alpha ≈ 90/a — checking
  the raw sign patterns of (N1, N2, N3) from the closed forms: **zero
  Case I or Case II hits** (float margin 1e-9). Positive control: just
  above theta_d at the predicted corner, Case I fires for every member —
  so the search is not vacuous.

### 3. The binding identities I1-I3. **CONFIRMED (independent software, two representations)**

Three independent layers, none sharing code with 003's Laurent ring:

- *My own ring* (`identities`, `dlsk_identities.json`): composed-map
  unfolding over my own Laurent arithmetic in Q(i)[X^±, Y^±],
  X = e^(i·alpha/2), Y = e^(i·beta/2). For all 15 members: I1, I2, I3
  subtract to the exact zero polynomial, as do the context identities D1,
  D2, the gate identifications (gate 1 = [B,C], gate 2 = [C,A1] — the
  second endpoint is genuinely R_0(A)), mu = delta² monomial,
  Im(conj(delta)·tau) = 0, and my own glide-action identities. My delta is
  derived independently (integer mirror-angle bookkeeping: mu =
  (-1)^r X^{2p} Y^{2q}, delta = X^p Y^q or i·X^p Y^q by parity of r) and
  reproduces 003's identities with the same signs — the sign convention is
  consistent, not fitted.
- *Off-torus grid proof for W(2,1)* (`gridproof`, `dlsk_gridproof_W21.json`):
  a representation-free certification that does not trust ANY polynomial
  code, mine included. Track pairs (value, star-value) at a fixed rational
  point (u0, v0) off the unit torus — star (the ring involution X→1/X,
  Y→1/Y, i→-i) is an automorphism, so the pair algebra evaluates every
  expression built from +, -, ×, star exactly in Q(i). Criterion: every
  residual is a Laurent polynomial whose exponents are bounded a priori by
  D = 6·n' + 8 (n' = half-word length; recursion bound |L_k| ≤ 4k,
  |t_k| ≤ 4(k-1)+6, base entities ≤ 6, delta ≤ 2n', so every checked
  residual is ≤ 6n'+6); a Laurent polynomial with exponents in [-D, D]²
  vanishing on a (2D+1)×(2D+1) grid of distinct rational abscissae is
  identically zero (clear denominators, interpolate one variable at a
  time). For W(2,1): D = 50, all 15 residuals vanish at all 101×101 =
  10,201 grid points. **This is a proof of the W(2,1) identities that
  shares no representation with either ring implementation.**
- *Spot extension to all 15 members:* the same pair algebra at 20 seeded
  random rational off-torus points per member — all residuals exactly
  zero (coincidence would require my ring bug and 003's ring bug to agree
  with 300 exact rational evaluations of the true geometry; not a proof,
  recorded as the cross-check it is; the per-member proof is the ring
  subtraction, which two implementations now give independently).

### 4. Lemma C certification, including the H2' leg. **CONFIRMED (independent interval stack)**

Re-certified from scratch for (2,1) [the b = 1 degenerate case], (5,3),
(10,9), (12,12) — `dlsk_lemmac_W*.json` — with my own intervals: Machin pi
(width < 1e-60, checked by rational bracket, not by float compare — the
double nearest pi actually sits below my PI.lo), my own alternating-series
sin/cos with explicit remainders. All three legs pass:

- *Endpoint algebra by hand:* H2(v0) = 0 because 90 - v0 = b·v0 (exact
  rational identity, asserted in code); H2(v) = -∫_v^{v0} H2' > 0 when
  H2' < 0 on [v0-d2, v0] — sound. I re-differentiated H2 by hand; my four
  terms match 003's `h2p_iv` exactly, and dropping the global pi/180
  factor is sign-safe since every term carries it.
- *The flagged H2' sign leg:* my certified upper bounds for H2' on
  [v0-d2, v0] are -0.055 (2,1), -0.034 (5,3), -0.0102 (10,9),
  -0.0069 (12,12) — comfortably negative, not marginal.
- *Near-0 leg re-derived:* sin(bv/a°) ≤ (b/a)·tan(v°) from
  sin x ≤ x ≤ tan x in radians; the endpoint-monotonicity of the bracket
  ((90-v)/a decreasing, 1/cos increasing) is used correctly.
- *Bisection leg:* my own adaptive bisection closes with 37-114 leaves.

003's `certify_lemma_c` code was also read line-by-line: the d2/d1
fallback loops fail loudly (`else` on the `while`), enclosures are single
non-iterated evaluations (the STATUS.md affine-forms lesson does not bite),
and `unfold.py`'s interval sin/cos midpoint-Lipschitz enclosures are sound
for these argument ranges.

### 5. Definitional identity of "death" between 001 and 003. **CONFIRMED (with one scoping note)**

001's `census.py cmd_family` bisects death as sup gamma over arcs where the
float corridor width (same tau-normal gate-projection intersection) is
positive, x scanned in [x_lo, 0.5]; 002's re-measurement is the same
criterion in independent code; the harness point criterion is the exact
version of the same corridor. My selftest ties all of these together: my
exact corridor sign equals `unfold.certify`'s verdict at point boxes and my
float corridor equals 002's to 1e-9. So 003's theorem (about corridor
width) bounds exactly the quantity 001 measured — the out-of-sample
"confirmations" compare like with like. One scoping note: 001's
death measurements scanned only x ≤ 1/2, which is fine for the members 001
measured (their alive windows lie there — my `death` runs confirm the
argmax at alpha ≥ beta for b ≥ a-1) but would have missed the b ≤ a-2
windows entirely; 003 says exactly this and used --full. Consistent.

### 6. The out-of-sample death angles and exact certificates. **CONFIRMED**

My own adaptive measurement (own corridor, own sampler with geometric
accumulation on BOTH sides of both window edges, 1200 uniform points, 46
bisection iters), `dlsk_death_W*.json`:

| word | predicted | 003 measured | my measured | mine - pred | window half |
|---|---|---|---|---|---|
| W(10,9)  | 162.9        | 162.8999999999984 | 162.89999999999833..36 | -1.7e-12 | alpha ≥ beta |
| W(12,12) | 2160/13      | 166.1538461538449 | 166.15384615384482..85 | -1.3e-12 | wide (a=b) |
| W(5,3)   | 144          | 143.9999999999959 | 143.99999999999594..96 | -4.0e-12 | **alpha < beta only** |
| W(4,2)   | 135          | 134.9999999999942 | 134.9999999999942..23  | -5.8e-12 | **alpha < beta only** |

Calibration members reproduce the law the same way (W(3,3): -5.8e-12,
W(2,1): -1.4e-11 — consistent with 003's "< 2e-11" calibration claim and
with the shared one-sided float-tolerance bias; see Residual risk). The
last-alive apexes sit at the predicted corners: (9, 8.1) for W(10,9),
(18, 18) for W(5,3), (22.5, 22.5) for W(4,2). **Mirror-half claim
verified:** near death, W(5,3) and W(4,2) have every alive sample at
alpha < beta (x > 1/2) and the best width on the alpha ≥ beta half is
negative — the half 001/002 scanned is genuinely dead there, and 001's
family list never included these members.

Exact certificates (`alivecheck`, `dlsk_alivecheck.json`): all 8
`deathlaw_exact_alive_W*.json` apexes re-checked — my exact corridor is
positive at every one, and the exact width **matches 003's recorded
rational digit-for-digit** in all 8 (hundreds of digits; my composed-map
path vs their reflected-chain path — same normalization, so this is a
strong two-implementation agreement, not a copy). The certified gamma
brackets (gamma > claimed bound, gamma < gamma_d) all re-certify under my
own pi/cos enclosures. The recorded orbit simulations (002's
`skeptic_orbit`) were not re-run — the corridor positivity, which is the
alive criterion, is what I re-established.

Kill attempts above the law (`deadhunt`, `dlsk_deadhunt_W*.json`): full-arc
scans with corner accumulation at gamma_d + {0, 1e-10, 1e-8, 1e-6, 1e-4,
1e-2, 0.5} for (5,3), (4,2), (10,9), (12,12): **no alive point** (best
float widths ≤ -2.6e-11 above gamma_d; at gamma_d exactly the best width
is float-zero ~ +2.7e-15 at the touching corner, see correction C2). An
exact probe at W(5,3)'s corner — rational apex with MY-certified
gamma ≥ 144 — has exact corridor width < 0.

### 7. Scope honesty. **CONFIRMED, two corrections below**

- The 15-member list is exactly `deathlaw_prove_all.json`'s content; all
  obligations there are `ok`, every member satisfies a ≥ b ≥ 1, a ≥ 2,
  a ≤ 2b+3, and the preconditions are used precisely where the record
  says (a ≤ 2b+3 only in A4/A4'; the a > 2b+3 hole is real — 270/a drops
  below theta_d at a = 2b+4, so branch A4' genuinely breaks).
- The certificate constants in the index/range and "Death brackets"
  sections match the data files (1e-6 everywhere except 1e-8 for (3,3)
  and (4,2); W(12,12) at 2160/13 - 13/13e6 = 2160/13 - 1e-6). One prose
  sentence contradicts this — correction C1.
- SPECULATION labels for general-(a,b) I1-I3 and Lemma C are present and
  placed at the right steps; the NOT-claimed list (no sufficiency limit,
  no unstable orbits, no other words, no birth side) matches what the
  certificates actually cover; float deadscan correctly labelled
  sample-bounded EVIDENCE and superseded by the theorem.

## Refutations found

No load-bearing claim is refuted. Two corrections:

- **C1 (wrong constant in one sentence).** "The general law" section:
  "both with exact alive certificates at gamma_d - 1e-8
  (`data/deathlaw_exact_alive_W53.json`, `..._W42.json`)" — the W(5,3)
  certificate is at gamma_d - **1e-6** (certified bound 143999999/1000000;
  re-certified here), not 1e-8; only W(4,2) (and W(3,3)) reach 1e-8. The
  record's own Outcome, range, and Death-brackets statements have it
  right; this sentence alone is wrong by two orders of magnitude.
- **C2 (wording overreach in the Outcome bullet).** "the W(a,b) corridor
  is **empty** at every triangle with gamma ≥ gamma_d" — the theorem
  proves the corridor has no **positive width** there (003's own theorem
  statement is phrased correctly). A zero-width corridor (a touching
  line, the limiting orbit) is not excluded at gamma = gamma_d — and
  numerically that is what the death corner looks like (float width
  ~ +3e-15 ≈ 0 at (18,18) on the gamma = 144 arc). Since "alive" is
  defined as positive width everywhere in 001-003, nothing downstream is
  damaged; the bullet should say "has no positive width" or "has empty
  interior".

Nits, no action needed: `deathlaw_exact.py cmd_alive` computes `theta`
unused; `deathlaw_symbolic.py cmd_structure` computes `lo_lbl` unused; the
`--dg 1e-8` in the T1 command line yields a certified bound of only 1e-6
for W(10,9) (the snap sits between 1e-8 and 1e-6 of the prediction), which
the record states correctly but a reader could misread as certifying 1e-8.

## Claims that survive

| # | 003 claim | Verdict |
|---|-----------|---------|
| 1 | Glide-to-corridor reduction (positive width ⟹ m strictly inside gates 1-2) | **CONFIRMED** — hand re-derivation incl. boundary strictness; ring-verified glide action; 1818-triangle numeric identity of full vs reduced corridor |
| 2 | Mod-360 case tree, all branches, scope preconditions | **CONFIRMED** — hand re-enumeration (boundaries self-excluding; threshold algebra exact); 1.94M-point adversarial sign search, 0 hits, positive control fires |
| 3 | I1-I3 (+ D1, D2, glide facts) exact for all 15 members | **CONFIRMED** — my own ring: 15/15; off-torus grid interpolation proof for W(2,1); 300 exact off-torus spot points |
| 4 | Lemma C certified; H2' leg sound | **CONFIRMED** — independent interval stack, 4 members incl. b = 1; H2' bounds comfortably negative; endpoint algebra re-derived |
| 5 | 001-death and corridor-death are the same criterion | **CONFIRMED** — census/skeptic/harness/my code tied together at point boxes; x ≤ 1/2 scoping consistent |
| 6 | W(10,9) → 162.90, W(12,12) → 2160/13, W(5,3) → 144, W(4,2) → 135 (~1e-12); mirror-half windows; 8 exact alive certificates | **CONFIRMED** — independent re-measurement agrees to ~1e-12; exact widths match digit-for-digit; my own gamma certification; deadhunt above gamma_d finds nothing |
| 7 | Scope statements match certificates | **CONFIRMED with C1, C2** |

**Kill attempts that failed,** for the record: (i) boundary-targeted
sign-pattern search for a dropped mod-360 branch (1.94M points, all 15
members) — nothing; (ii) full-arc corner-accumulated alive hunts at and
above gamma_d on four members, with exact confirmation armed for any float
positive — nothing (all float positives at gamma_d are ≤ 3e-15 ≈ 0, and
the exact corner probe is negative); (iii) attempting to break the
strictness of the reduction at window edges — the point-overlap case
provably cannot carry positive width; (iv) hunting for a sign-convention
mismatch between delta conventions — my independently derived delta
reproduces the identities with identical signs.

**Net assessment:** 003's headline — both 001 laws confirmed out-of-sample,
unified into gamma_d = 180 - 90(a+b)/(a(b+1)), and the necessity theorem
machine-certified for 15 members — stands as scoped. The status upgrade of
the death laws (upper side: theorem; lower side: exact alive at
gamma_d - 1e-6/1e-8) is justified. The general-(a,b) statements remain
SPECULATION exactly as labelled.

## Residual risk

- **Shared-design bias in the float death numbers.** My measurement and
  003's share the design "death = last arc with a positive-width sample"
  and an absolute alive tolerance (1e-12); with window widths pinching
  linearly, both under-estimate by ~1e-12..1e-11 deg in the same
  direction. The ~1e-12 agreements are therefore two instruments of the
  same type agreeing; the exact brackets [gamma_d - 1e-6, gamma_d] are
  what carries the result, and those are certified.
- **Lemma D** is taken as classical: re-derived by hand (tangent-addition
  induction) and numerically spot-checked here, but not machine-certified
  by 003 or by me.
- **Orbit realizations** (the "26/34/78/98 bounces confirmed" lines in
  003's alive files) rest on 002's simulator, which 002 validated; I
  re-verified the corridor criterion, not the trajectories.
- **The three ring implementations share the mathematical reduction**
  (circumdiameter normalization, division-free reflections). The off-torus
  grid proof (W(2,1)) and the apex-coordinate cross-checks (exact width
  match on all 8 certificates; harness point-box agreement) bound this
  risk, but a conceptual error in the shared normalization itself would
  have to show up in those apex-coordinate checks to be caught — it did
  not.
- **Nothing here touches the open gaps** 003 lists: sufficiency
  (death = gamma_d exactly), general-(a,b) I1-I3/Lemma C, a > 2b+3.

## References

- `problems/billiards-triangles/attempts/003-death-angle-laws.md` (under
  review); `001-word-census-coverage-map.md` (corridor/death definitions);
  `002-skeptic-review-of-001.md` (independent corridor reused only as a
  selftest cross-reference).
- Code under review: `problems/billiards-triangles/explore/deathlaw_{measure,exact,symbolic,prove}.py`;
  its data `data/deathlaw_*.json`.
- Tier-0: `harness/billiards-triangles/unfold.py` (point-box corridor
  cross-reference in the selftest only).
- New code and data: `problems/billiards-triangles/explore/deathlaw_skeptic.py`;
  `problems/billiards-triangles/data/dlsk_{identities,gridproof_W21,reduction,casetree,alivecheck}.json`,
  `data/dlsk_death_W{109,1212,53,42,33,21}.json`,
  `data/dlsk_deadhunt_W{53,42,109,1212}.json`,
  `data/dlsk_lemmac_W{21,53,109,1212}.json`.
- No external papers consulted.
