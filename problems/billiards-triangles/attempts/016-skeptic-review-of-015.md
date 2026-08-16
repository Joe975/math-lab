# 016 — Skeptic review of 015 + the seven S-lemmas proven: the window law closes for all (a,b)

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-08-16
- **Mode:** informed
- **Type:** adversarial verification of `015-window-law-margin-system.md`
  (default stance: refute), taking exactly the attack surface 015's lead 1
  pre-registered — the T2 Step-5 repair, the N1 proof, the ring-bijection
  scope, the certifier's range checks — PLUS the queued main task: re-derive
  or refute the seven S-lemmas of 015 §5 and the τ row from scratch.
- **Outcome in one line:** everything load-bearing SURVIVES and the seven
  S-lemmas (and τ) are now PROVEN for all a ≥ b ≥ 1, a ≥ 2 at k = 1/2 —
  so the window law `window(W(a,b)) = (γ_birth, γ_d)` and, via 013/014,
  coverage of (90°, 180°) hold for ALL members, closing queue 11 — but the
  review found two real record errors: 015 §1's identity for the s = 0
  C-image margin is wrong (it is −N1, not −N1−N2), and 015 §4's "every
  margin sign on BOTH segments is certified" is false as stated: the
  death-segment BC0 (= N4) row is never evaluated by
  `wmargins_certify.py` (its death B-fan C loop starts at j = 1 and no
  N4 leg exists). That row is true — and is certified exactly here for
  all a ≤ 14 — so no theorem falls, but the record's §5 reduction table
  is missing an N4-death family and its §4 sentence overclaims what the
  a ≤ 14 run certified.
- **Independence protocol:** the swarm drafts of the S-lemmas live only
  under the gitignored `$MATHLAB_OUT/swarm/wmargin-*/` and do not exist in
  this fresh checkout; they were never read. Only the committed lemma
  STATEMENTS (`explore/wmargin_brief.md`, `explore/wmargin_tasks.txt`)
  were used. Every proof below is a from-scratch derivation; every named
  step of every proof was then adversarially machine-scanned.
- **Tools:** `explore/wmsk015.py` (new, stdlib only, deterministic):
  `selftest` (~90 s) — my own margin implementation validated against the
  geometry (002's blind-era `skeptic_family.width`, 1200/1200 sign
  agreements incl. far-from-window triangles), identity checks, τ closed
  form vs the exact symbolic τ, fresh 62-member ring bijection;
  `lemmas` (~6 min) — 32 tracked inequalities (the 7 S-lemmas, τ, N4, and
  every intermediate proof step) scanned over all (a,b), a ≤ 40 × 96
  t-points plus 20,000 random members to a = 3000: 0 violations;
  `t2` (~2 min) — 453,605 adversarial beyond-the-wall candidates
  (including 360°-higher branches): 0 alive; `n4leg` — the missing
  certificate leg, exact (ring identity + rational ranges), 104/104
  members a ≤ 14. Committed tools re-run as published:
  `wmargins.py selftest` PASS, `wmargins_certify.py member --a 5 --b 3`
  PASS (its `death_fan_leaves = 13` is itself the finger-print of the
  missing BC0 row: 14 fan rows exist beyond N1/N2/PB/AB_a).
- **Sources:** prior attempts 003–008, 013/014, 015 (this repo). No
  external literature.

## Notation (degrees throughout)

A = 90/a, H = 90/(b+1), c = 90(a−1)/(a(b+1)), T = a(b+1), so
(b+1)c = 90 − A, θ_d = A + c = 90(a+b)/T, θ_birth = A + H, 90 − θ_d = bc.
Death segment: (α, β) = (A − t, c + 2t); birth segment: (α, β) =
(A − t, H − t); both with t ∈ (0, 22.5/T] (k = 1/2). On the death
segment aα = 90 − at, (b+1)β = 90 − A + 2(b+1)t, θ = θ_d + t; on the
birth segment aα = 90 − at, (b+1)β = 90 − (b+1)t, θ = θ_birth − 2t.

Small facts used repeatedly, all elementary (F1–F7): **F1**
at ≤ 22.5/(b+1) ≤ 11.25, (b+1)t ≤ 22.5/a = A/4, α ≥ (7/8)A,
2(b+1)t ≤ A/2. **F2** (chord) sin x ≥ (x/y)·sin y for 0 ≤ x ≤ y ≤ 180.
**F3** sin x ≥ 0.9·(π/180)·x for 0 ≤ x ≤ 45 (sinc is decreasing;
sinc(45°) = 0.9003). **F4** if X ∈ [L, 180−L] with 0 < L ≤ 90 then
sin X ≥ sin L. **F5** sin(λx) ≥ λ·sin(x)·cos(λx) for every λ > 0 with
0 < λx ≤ 90 (proof: tan u ≥ u ≥ λ·sin x in radians, u = λx·π/180).
**F6** sin(nx) ≤ n·sin x for n ∈ ℕ, x ∈ [0, 180]. **F7**
a·sin(90/a) < π/2 for all a ≥ 2 (increasing in a, limit π/2), so
a·b·sin A < (π/2)·b. On the birth segment θ < θ_birth ≤ 90 (equality
of θ_birth only at (2,1)); on the death segment θ < 90 and
sin θ = cos(bc − t).

## Claims attacked

1. **T2 (birth-side necessity), including the Step-5 repair** — the
   sign-constant-progression pigeonhole, the aα ∈ [90, 120) branch, the
   B-fan interleaving, and the b = 1 edge case.
2. **The closed-form margin system of 015 §1** — the 2a+2b+2 forms, the
   bijection claim, the stated identities (s = 0 C-image, j = 0 C-image,
   AB_a = −PB), and whether the closed forms really decide the geometric
   corridor (not just match the symbolic engine that produced them).
3. **The certified window theorem for a ≤ 14 (015 §4)** — does
   `wmargins_certify.py` actually certify every margin it says it does?
4. **The N1 P-row** (015's hardest director-verified lemma).
5. **The seven S-lemmas** [D-AC], [D-AB], [D-BC], [D-BA], [B-AC],
   [B-AB], [B-BC] **and the τ row** — re-derive from scratch or refute.

## Refutations found

None load-bearing. Three record errors / overclaims, all repaired here:

1. **015 §1: "the s=0 C-image margin is −N1−N2 (I1)" is FALSE.** The
   s = 0 A-fan C-image margin equals **−N1** exactly (N1 in the
   brief's/005's form: on the death segment N1 = cos(at)·sin β·
   sin((2b+1)t) − sin(at)·sin α·sin(bβ)). Machine check: max
   |AC0 + N1| = 3.2e-16 over 4000 random member/segment points, while
   −N1−N2 differs by N2 = O(t) ≫ float error. This matches
   `plaw_suffice.py`'s own I1 ("m − p(C) = N1"). Nothing downstream is
   damaged — N1 > 0 (with or without N2 > 0) still gives AC0 < 0 — the
   identity as printed is simply wrong.
2. **015 §4: "every margin sign on BOTH segments is certified in exact
   rational interval arithmetic ... for ALL 104 members" is OVERSTATED.**
   `wmargins_certify.py` never evaluates the death-segment BC0 row: its
   `death_margins` B-fan C loop is `range(1, b)`, and unlike N1/N2/PB
   (handled by separate legs, as its docstring says) there is no N4 leg
   anywhere in `certify_member`. Finger-print: for W(5,3) the death
   segment has 2a+2b+2 = 18 distinct margins; N2/PB are rational legs,
   AC0 is the N1 leg, AB5 = −PB, leaving 14 fan rows — the run reports
   `death_fan_leaves` covering only 13. The row itself is TRUE: BC0 = N4
   (identity I4; re-verified as an exact ring identity here), and on the
   death segment N4 = sin(at)·sin α·sin(bβ) + cos(at)·sin β·sin((2b+1)t),
   a sum of two products of factors that are positive by rational range
   checks — 005 §(iv)'s own argument, which is parametric in (a,b). So
   the a ≤ 14 window theorem stands, but one of its margins was proven by
   005's I4, not certified by that run. Repaired: `wmsk015.py n4leg`
   certifies BC0 = N4 (exact ring subtraction, zero tolerance) plus the
   positive-factor ranges (exact rational endpoints) for all 104 members
   a ≤ 14: ALL CERTIFIED. Consequently **015 §5's status table is missing
   a family**: the all-(a,b) reduction needs an "N4 > 0, death" P-row
   (via I4 + rational ranges, valid for every k < 2), and the record's
   "reduces to 9 elementary lemma families" / "seven S-rows" counts do
   not reconcile with its own table (8 S rows incl. τ; 13 rows total incl.
   the missing N4). Bookkeeping only — every family is (now) proven.
3. **Certifier docstring nit:** "excluding N1/N2/PB which are handled
   separately" should also name BC0/N4 (and AB_a on the death side,
   which is −PB). Same omission as (2), recorded so the tool's contract
   is honest.

Also checked and NOT errors: the τ brief's flagged δ-vs-−δ convention
(irrelevant — a sign flip of δ flips R's sign, not R ≠ 0, and the closed
form below settles the convention); luna's claimed D-AB range
k = (b+1)/(b+2) (my proof gives k = 1 there; a weaker claim, not a wrong
one); the windows-touch one-liner and the 013 interface identity
θ ∈ (θ_d, θ_birth) ⟺ a+b < tT < a+b+1 (both re-checked, exact).

## Claims that survive

1. **The margin system (015 §1).** My own implementation of the closed
   forms agrees in alive/dead sign with the *geometry* — 002's blind-era
   corridor width, an implementation predating and independent of the
   whole 003→015 stack — at 1200/1200 random triangles (a ≤ 12, near
   both corners and far from the window; the committed selftest's 400
   were near-corner only). The exact ring bijection (2a+2b+2 distinct
   diffs, zero tolerance) holds on a 62-member set (all 54 members with
   a ≤ 10 — 015 checked only 20 of them — plus (11,2), (12,7), (13,13),
   (14,1), (16,3), (17,16), (19,1), (20,7), all eight disjoint from
   015's list; a > 2b+3 shapes included). Union with 015's 25: 67
   members ring-verified.
   AB_a = −PB and BC0 = N4 verified to 4e-16; AC0 = −N1 as in
   Refutation 1. The general-(a,b) endpoint bookkeeping remains, as 015
   itself flagged, a hand derivation + per-member exact check
   (now 25+62 members) — see Leads.
2. **T2 and its Step-5 repair (015 §3): CONFIRMED,** re-derived in full
   with one simplification. My derivation: (i) the fan margins with
   N2 > 0, sin β > 0, sin θ > 0 force sin(z_n) < 0 for
   z_n = aα − (b+1)β − nα, n = 1..2a (AC rows give odd n, AB rows even
   n). (ii) A sign-constant progression with step α < 90 cannot cross a
   positive-sine arc of length 180, so all z_n lie in one arc
   (180+360m, 360+360m); the span gives (2a−1)α < 180, hence
   aα < 180a/(2a−1) ≤ 120. (iii) If aα ∈ (90, 120) (aα = 90 is killed
   by N2 ≠ 0): cos(aα) < 0, so N2 > 0 forces sin((b+1)β) < 0. But
   **n = a gives z_a = −(b+1)β, so sin(z_a) = −sin((b+1)β) > 0**,
   contradicting (i) directly — a one-line kill; the record's route
   (PB refines (b+1)β into (270, 360) mod 360 and z₁ lands outside
   z_{2a}'s arc) was also checked and is correct, just heavier. So
   aα < 90. (iv) Then N2, PB > 0 give (b+1)β ∈ (360m′, 360m′+90). (v)
   The B-fan rows give sin(aα + rβ) > 0 for every integer
   r ∈ [−(b−1), b] (BC supplies r ≡ b, BA supplies r ≡ b+1 mod 2 —
   the interleaving is exactly as repaired in 015); for b ≥ 2 the same
   pigeonhole gives (2b−1)β < 180, and m′ ≥ 1 needs β ≥ 360/(b+1):
   360(2b−1) < 180(b+1) ⟺ 3b < 3, false; b = 1: 2β ≥ 360 contradicts
   β < 90. So (b+1)β < 90. ∎ Numerically: 453,605 adversarial
   beyond-wall candidates (including (b+1)β ∈ 360+(90..162) branches),
   zero alive.
3. **The N1 P-row: CONFIRMED,** by a new derivative-free proof (015/005
   needed N1′(0) > 0 plus bisection; this needs neither). By F5 with
   x = at, λ = (2b+1)/a (note (2b+1)t < 22.5 ≤ 90):
   sin((2b+1)t) ≥ ((2b+1)/a)·sin(at)·cos((2b+1)t). By F6
   sin(bβ) ≤ b·sin β, and sin α < sin A, so
   N1 ≥ sin(at)·sin β·[((2b+1)/a)·cos(at)·cos((2b+1)t) − b·sin A], and
   the bracket is positive: by F7, a·b·sin A < (π/2)b, while
   (2b+1)·cos(at)·cos((2b+1)t) ≥ (2b+1)·cos(11.25°)·cos(22.5°) =
   (2b+1)·0.9062 > (π/2)·b since 0.9062·(2b+1) − 1.5708·b =
   0.2416·b + 0.9062 > 0. Strict for every t in (0, 22.5/T]. ∎
   (This also re-proves the N1′(0) sub-lemma statement en passant — it is
   the t → 0 limit — and would let the certifier drop its derivative leg.)
4. **The certified window theorem for a ≤ 14** — survives with
   Refutation 2's repair: their run + the n4leg here jointly certify
   every margin on both segments for all 104 members. The sharp wall
   t\* = 45/T and the two-segment sweep bookkeeping re-checked.
5. **The seven S-lemmas and τ: ALL PROVEN** — complete proofs below,
   every step machine-scanned (32 tracked inequalities, min margins all
   strictly positive over a ≤ 40 dense + 20k random members to
   a = 3000).

## The lemma proofs (all a ≥ b ≥ 1, a ≥ 2; t ∈ (0, 22.5/T])

**[D-AB]** Claim: sin(L_s) > sin(at)·sin((b+1)β) for s = 1..a−1,
L_s = (2s−1)A + (a+2b+2−2s)t. Since sin((b+1)β) ≤ 1 it is enough that
L_s ∈ (at, 180 − at), which gives sin L_s > sin(at) (at < 90).
Lower: L_s − at = (2s−1)A + 2(b+1−s)t; for s ≤ b+1 both terms are ≥ 0
(the first > 0); for s ≥ b+2, 2(s−b−1)t ≤ 45(s−b−1)/T =
(A/2)(s−b−1)/(b+1) ≤ (A/2)(s−2) < (2s−1)A. Upper: L_s + at =
(2s−1)A + (2a+2b+2−2s)t increases in s (step 2A − 2t > 0); at s = a−1
it is 180 − 3A + (2b+4)t and (2b+4)t ≤ (A/2)(b+2)/(b+1) ≤ (3/4)A < 3A. ∎
(The same numbers give k = 1: the claimed luna range was weaker.)

**[D-BA]** Claim: cos(at − rβ) > sin(A − 2(b+1)t)·cos(at) for
r = b+1−2j, j = 1..b (so |r| ≤ b−1). Since at + (b−1)β < 180,
cos(at − rβ) ≥ cos(at + (b−1)β). Now at + (b−1)β < 90 − A: using
(b−1)c = (90−A)(b−1)/(b+1), this is (a+2b−2)t < (90−A)·2/(b+1) =
180(a−1)/T, and 22.5(a+2b−2) < 180(a−1) ⟺ 2b + 6 < 7a, which holds for
all a ≥ 2, b ≤ a (at a = 2, b ≤ 2: 2b+6 ≤ 10 < 14). Hence
cos(at − rβ) > cos(90 − A) = sin A > sin(A − 2(b+1)t) ≥ RHS/cos(at),
and multiplying by cos(at) ≤ 1 finishes. ∎

**[D-BC]** (empty unless b ≥ 2) Claim: sin α·cos(at − rβ) > PB for
r = b−2j, j = 1..b−1 (so |r| ≤ b−2), PB = sin(A−2(b+1)t)·cos(at)·sin θ.
First sin α > sin(A − 2(b+1)t) (α = A − t > A − 2(b+1)t, both in
(0, 90)). Second, cos(at − rβ) ≥ cos(at + (b−2)β) > cos(bc − t) = sin θ:
the angle comparison at + (b−2)β < bc − t reduces, via
(b−2)β = (b−2)c + 2(b−2)t and bc = (b−2)c + 2c, to (a+2b−3)t < 2c =
180(a−1)/T, i.e. 22.5(a+2b−3) < 180(a−1) ⟺ 2b+5 < 7a, true for b ≤ a,
a ≥ 2; both angles lie in (0, 90). Multiply the two and use
cos(at) ≤ 1: sin α·cos(at−rβ) > sin(A−2(b+1)t)·sin θ ≥ PB. ∎
(Also holds at k = 1, consistent with the luna claim "to the wall".)

**[D-AC]** Claim: sin β·sin(M_s) > N2 for s = 1..a−1,
M_s = 2sA + (a+2b+1−2s)t, N2 = sin(at)·sin((b+1)β)·sin θ.
*Step 1 (window).* M_s increases in s (step 2α) and
M_{a−1} = 180 − m₀ exactly, where m₀ := 2A + (a−2b−3)t; also
M_1 ≥ m₀ and 0 < m₀ ≤ 90 ((2b+3)t ≤ (5/8)A gives m₀ ≥ (11/8)A; for
a = 2 the t-coefficient is negative, for a ≥ 3, 2A + at ≤ 71.25). By
F4, sin(M_s) ≥ sin(m₀). *Step 2 (core).* Dropping sin((b+1)β) ≤ 1, it
is enough that Δ := 2[sin β·sin(m₀) − sin(at)·sin θ] > 0. Writing
β = c + 2t, θ = A + c + t and converting products to sums, Δ pairs
EXACTLY (machine-checked identity) as
Δ = 2·sin((3/2)A − (b+2)t)·sin(φ₁) + 2·sin((3/2)A + c + (a−b)t)·
sin(A/2 − (b+1)t), with φ₁ = c − A/2 + (b+3−a)t. The second pair is
always positive: A/2 − (b+1)t ≥ A/4 by F1, and the other argument lies
in (0, 124°). The first factor of the first pair is positive
((b+2)t ≤ (3/8)A). And φ₁ = 45(2a−b−3)/T + (b+3−a)t ≥ 0 for every
(a,b) ≠ (2,2): if a ≥ 3 then 2a−b−3 ≥ a−3 ≥ 0 and if additionally
a > b+3, φ₁ ≥ 22.5(3a−b−3)/T > 0; at (2,1), φ₁ = 2t > 0. So Δ > 0
except possibly at (2,2). *Step 3, W(2,2)* (A = 45, c = 15, T = 6,
t ∈ (0, 3.75]): Δ/2 = sin(67.5−4t)·sin(3t−7.5) + sin(82.5)·sin(22.5−3t).
For t ≥ 2.5 both terms are ≥ 0 with the second > 0. For t < 2.5,
|first| ≤ sin 67.5°·sin 7.5° < 0.121 while
second ≥ sin 82.5°·sin 15° > 0.256. ∎

**[D-N1]** and the **N4-death row**: Claims-that-survive items 3 and
Refutation 2 respectively (N1 by the derivative-free proof; N4 =
sin(at)·sin α·sin(bβ) + cos(at)·sin β·sin((2b+1)t) with every factor
positive by rational ranges — bβ = bc + 2bt < 90 + 22.5, (2b+1)t < 22.5).
Together with N2, PB (rational ranges) and AB_a = −PB these close every
death-segment row.

**[B-AB]** Claim: sin(G_s) > sin(at)·cos((b+1)t) for s = 1..a−1,
G_s = 2sA + (a−b−1−2s)t. Enough that G_s ∈ (at, 180−at). Lower:
G_s − at = 2sA − (b+1+2s)t, minimized at s = 1: 2A − (b+3)t and
(b+3)t ≤ (A/4)(b+3)/(b+1) ≤ A/2 < 2A. Upper: G_s + at =
2sA + (2a−b−1−2s)t, at s = a−1: 180 − 2A + (1−b)t ≤ 180 − 2A. ∎

**[The birth core]** Claim: sin α·sin(H − (a−b)t) >
sin θ·cos(at)·sin((b+1)t). Proof: (a−b)t ≤ at ≤ 22.5/(b+1) = H/4, so
H − (a−b)t ≥ (3/4)H; also H − (a−b)t < θ (their difference is
A + (a−b−2)t ≥ A − 2t > 0) and θ ≤ 90, so by F2
sin(H − (a−b)t) ≥ [(H − (a−b)t)/θ]·sin θ. Meanwhile
sin((b+1)t) ≤ (π/180)(b+1)t ≤ (π/180)·A/4 and, by F1+F3,
sin α ≥ 0.9·(π/180)·(7/8)A. If b < a (so A ≤ H): LHS ≥
sin α·sin θ·(3H/4)/θ and it suffices that 0.9·(7/8)·(3/4)·H =
0.5906·H > θ/4, i.e. θ < 2.3625·H — true since θ < A + H ≤ 2H. If
b = a, (a−b)t = 0 and the same chain with H in place of (3/4)H needs
θ < 3.15·H — true since then A/H = (a+1)/a ≤ 3/2, so
θ < A + H ≤ (5/2)·H. ∎

**[B-AC]** Claim: sin β·sin(J_s) > N2 for s = 0..a−1,
J_s = (2s+1)A + (a−b−2−2s)t, N2 = sin(at)·cos((b+1)t)·sin θ.
*Step 1.* J_s increases in s; sin(J_{a−1}) = sin(A + (a+b)t) and both
mirror arguments A + (a−b−2)t = J_0-arg and A + (a+b)t lie in (0, 68°],
so by F4-type endpoint comparison sin(J_s) ≥ sin(J_0) for all s.
*Step 2 (exact identity, machine-checked).* sin β·sin(J_0) − N2 =
sin α·sin(H − (a−b)t) − sin θ·cos(at)·sin((b+1)t) — the two sides are
identical trig polynomials (product-to-sum; the (a−b−3)t terms cancel).
*Step 3.* Apply the birth core. ∎

**[B-BC]** Claim: sin α·cos(X_r) > PB for X_r = rH − (a+r)t,
r = b−2j, j = 0..b−1, PB = sin((b+1)t)·cos(at)·sin θ. X_r increases in
r (coefficient H − t > 0): X_r ≤ bH − (a+b)t < bH = 90 − H, and
X_r ≥ −(b−2)H − (a−b+2)t ≥ −(90 − 3H) − H/2 > −(90 − H) (since
(a−b+2)t ≤ at + 2t ≤ H/4 + H/(2a) ≤ H/2). So |X_r| < 90 − H, giving
cos(X_r) > sin H ≥ sin(H − (a−b)t), and the birth core finishes. ∎

**[B-BA]** (P-row, re-derived) Claim: cos(at − r′(H−t)) >
sin((b+1)t)·cos(at) for r′ = b+1−2j, j = 1..b (|r′| ≤ b−1). LHS ≥
cos(at + (b−1)(H−t)) = sin(2H − at + (b−1)t) (since 90 − (b−1)H = 2H);
the argument lies in (0, 90] and exceeds (b+1)t because
2H − at − 2t ≥ 2H − H/4 − H/(2a) > 0; so LHS > sin((b+1)t) ≥ RHS. ∎
Likewise **AB_a** (P-row): −AB_a = sin θ·[sin((a+b+1)t) −
sin(at)cos((b+1)t)] = sin θ·cos(at)·sin((b+1)t) > 0. ∎

**[TAU]** With μ = e^{2i(−aα+(b−1)β)}, δ² = μ,
w = sin θ·(1 + e^{−2i(aα+β)} − μ − e^{−2iβ}), τ = w + μ·conj(w):
conj(δ)·τ = conj(δ)w + conj(conj(δ)w) = 2·Re(conj(δ)w) is already real,
and multiplying out (four monomials; the ±cos(aα−(b−1)β) terms cancel)
gives the closed form

    R = Re(conj(δ)τ) = 2·sinθ·[cos(aα+(b+1)β) − cos(aα−(b+1)β)]
                     = −4·sinθ·sin(aα)·sin((b+1)β).

Machine-checked against the exact symbolic τ of `deathlaw_symbolic` on
41 members at random points (max error 8.9e-16). Consequence, stronger
than the S-row asked for: **R ≠ 0 whenever N2·PB ≠ 0** — N2 > 0 and
PB > 0 already force sin θ, sin(aα), sin((b+1)β) all nonzero — so
τ-nondegeneracy is automatic at EVERY alive point, for every (a,b),
with no range analysis; on both canonical segments R < 0 explicitly
(death: R = −4·sinθ·cos(at)·cos(A−2(b+1)t); birth:
R = −4·sinθ·cos(at)·cos((b+1)t)). The δ vs −δ convention only flips
R's sign. (The certifier's per-member τ bisection legs are hereby
subsumed.) ∎

## Consequences

With 015's T2 + 005's death-side necessity (containment) and the
segment sufficiency now parametric — every death-segment row by
N2/PB/N1/N4/AB_a + [D-AC]/[D-AB]/[D-BC]/[D-BA], every birth-segment row
by N2/PB/AB_a/B-BA + [B-AC]/[B-AB]/[B-BC], τ automatic — the two
canonical segments are alive for ALL a ≥ b ≥ 1, a ≥ 2 on
t ∈ (0, 22.5/T], and they realize every θ ∈ (θ_d, θ_birth). Hence:

- **Window law, all members:** window(W(a,b)) = (γ_birth(a,b), γ_d(a,b))
  exactly, as a set of realized angles, for every a ≥ b ≥ 1, a ≥ 2
  (suprema not attained on either side).
- **Coverage:** combined with the 013 lemma (skeptic-confirmed in 014):
  for EVERY obtuse γ ∈ (90°, 180°) there is an explicit W(a,b) and an
  explicit alive triangle with apex angle exactly γ — the coverage
  conjecture for the W family, no longer conditional on any S-row.

Scope and honesty: (i) the closed-forms-are-the-corridor bridge for
general (a,b) still rests on 005's proven prefix-map structure plus the
mechanical endpoint bookkeeping that 015 flagged — hand-derived,
exact-ring-verified per member (now 67 members), geometry-validated at
1200 random points, but not yet a formal all-(a,b) ring pass
(SPECULATION only in that narrow sense; lead 1). (ii) Per the lab's
standard, the proofs in THIS record are themselves unreviewed until a
skeptic pass re-derives them; every step is machine-scanned here, but
that is evidence, not review. (iii) The 007-C2 pointwise caveat stands:
one triangle per θ, which is the form the conjecture needs. (iv) The
a = 1 column and behavior AT the window endpoints remain open.

## Leads generated

1. **Skeptic pass on this record** (default stance: refute). Attack
   surface, pre-registered: the D-AC pairing identity and its φ₁ case
   split (the W(2,2) paragraph especially); the F5 lemma's λ < 1 case;
   the birth-core constant chain (0.9·7/8·3/4 vs θ/H ≤ 2 — small
   numeric floors, check them in interval arithmetic); the B-BC lower
   bound on X_r (the −(b−2)H − (a−b+2)t step); the T2 n = a one-liner
   (is z_a really in the forced index range for every parity?  n = a ∈
   [1, 2a] always — check the row-to-index map); every claimed "both
   angles in (0,90)" domain check.
2. **Formal all-(a,b) endpoint bookkeeping** (unchanged from 015 lead
   1): do for the six margin families what plaw_general did for I1–I4 —
   the A-fan offset is the two-liner R₀Rot_A(2sα)X − R₀Rot_A(2sα)A =
   e₀e^{−2isα}·conj(X) with s formal; this would retire caveat (i).
3. **Formalize (Lean lane): T2 + the seven lemmas.** All are sine-sign
   windows, arithmetic progressions and the F1–F7 facts — mathlib-
   friendly; T2 first (015 lead 3), then D-N1 (now derivative-free,
   much easier than 005's treatment), then the rest.
4. **Simplify the certifier** using this record: drop the N1 derivative
   leg (item 3's proof), drop the τ legs (closed form), add the N4 leg
   (n4leg), and the whole a ≤ 14 certificate re-runs leaner; or retire
   per-member certification entirely now the law is parametric.
5. **Ripple (015 lead 6 stands, sharpened):** the τ computation
   pattern — glide-projected translation collapsing to a three-sine
   product — and the progression pigeonhole are both candidates for the
   almost-mathieu Chambers-polynomial gap conditions; check before
   building new machinery there.

## References

- `problems/billiards-triangles/attempts/015-window-law-margin-system.md`
  (the record under review) and its committed briefs
  `explore/wmargin_brief.md`, `explore/wmargin_tasks.txt`.
- `problems/billiards-triangles/attempts/005-complete-death-law-theorem.md`
  and `008-skeptic-review-of-005.md` (prefix maps, I1–I4, N4 treatment,
  certified-trig layer); `003/004` (glide reduction);
  `013/014` (coverage Diophantine lemma).
- `explore/wmsk015.py` (this review's engine); committed tools re-run:
  `explore/wmargins.py`, `explore/wmargins_certify.py`,
  `explore/plaw_suffice.py` (N4 leg, read not re-run),
  `explore/skeptic_family.py` (002; the independent geometry).
- No external literature used; the trig facts F1–F7 are classical.
