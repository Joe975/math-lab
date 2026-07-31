# 003 — Death-angle laws of the W(a,b) family: predictions confirmed, closed forms found, necessity proven

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed (read `prior-art.json`, attempts 001 and 002 in full, and
  the 002 reviewer's tools under `explore/`; queue item 9 of `STATUS.md`)
- **Type:** computation + formalization — out-of-sample prediction tests,
  exact symbolic unfolding, and a certified proof of the death-angle upper
  bound
- **Tools:** new `explore/deathlaw_measure.py` (adaptive float death
  bisection), `explore/deathlaw_exact.py` (exact rational alive/dead
  brackets with certified angle enclosures), `explore/deathlaw_symbolic.py`
  (division-free exact unfolding over the ring Q(i)[e^{i alpha/2},
  e^{i beta/2}]), `explore/deathlaw_prove.py` (per-member certified proof
  obligations).  Reused: `skeptic_family.py` / `skeptic_orbit.py` (002's
  independent implementations — deliberately not census.py's code path),
  tier-0 `unfold.py` (interval sin/cos with rational pi enclosure only; the
  certificate path of the harness is not used here).  All stdlib-only,
  deterministic; total compute ≈ 10 min.
- **Sources:** none external.  Repo: attempts 001, 002; `STATUS.md` insight
  2026-07-28 on affine forms vs intervals (heeded: every certified interval
  evaluation here is a single non-iterated evaluation, where plain intervals
  are sound; nothing iterates a map over a box).

**Family and conventions (from 001).** W(a,b) = (0 (12)^a (02)^b)^2, length
4a+4b+2; side s is opposite vertex s (side 0 = BC, 1 = CA, 2 = AB); alpha,
beta are the angles at A, B; gamma = 180 - alpha - beta the angle at C; a
constant-gamma "arc" is the locus of apexes with that gamma.  "Alive at a
triangle" means the unfolding corridor of the word has positive width there
(the harness's TRUE criterion at a point); "death angle" = sup of gamma over
alive triangles.  001 measured (SPECULATION there):
death(W(a,a)) = 180a/(a+1) deg and death(W(a,a-1)) = 180(1-(2a-1)/(2a^2))
deg; 002 confirmed two out-of-sample rows, measured only.

## Approach

Queue item 9, in the ordered discipline it demands: (T1) test the two
stated out-of-sample predictions — W(10,9) dies at 162.90 deg, W(12,12) at
166.154 deg — before touching a proof, with exact arithmetic where a claim
is made; (T2) turn the corridor-endpoint geometry into exact symbolic
closed forms per a; (T3) attempt the parametric theorem.

Why symbolic-in-angles rather than the obvious alternative (the harness's
apex-coordinate rational functions): in apex coordinates every reflection
divides by a squared side length, so corridor endpoints are rational
functions whose structure is opaque and whose degeneration locus is buried
in resultants.  Normalizing the triangle to circumdiameter 1 instead
(A = 0, B = sin(alpha+beta), C = sin(beta) e^{i alpha} in the complex
plane) makes every mirror direction a unit complex number e^{i(...)} with
(...) an integer combination of alpha, beta — so reflection z -> u +
d^2 conj(z - u) is DIVISION-FREE in the Laurent ring
Q(i)[e^{i alpha/2}, e^{i beta/2}], and every unfolded vertex, the
translation, and every corridor projection is an exact trigonometric
polynomial.  That is what makes the death law a finite computation.

A second structural simplification found on the way: the half word
u = 0(12)^a(02)^b composes to a GLIDE REFLECTION (odd number of
reflections; W = u^2 a translation), so the full 2n-gate corridor
criterion collapses to "the glide-axis offset m lies strictly inside the
projection interval of every half-word gate" — and near death only gates
1 = [B, C] and 2 = [C, A_1] bind (A_1 = the reflection of A across BC).

## What was done

### T1 — out-of-sample predictions first

`deathlaw_measure.py death` bisects the death angle with an adaptive alive
test (uniform alpha-grid plus geometric grids accumulating at the two
window edges alpha = 90/a and beta = 90/(b+1); positives are exact-checkable
apexes, DEAD verdicts are bounded by this sampling).  This removes the
fixed-grid bias that made 001/002's measured deaths sit 0.001-0.006 deg low
(their alive windows pinch linearly, ~0.075 deg of x-window per deg of
gamma, below any fixed grid).  Calibration on the six previously measured
members reproduces the law to < 2e-11 deg.  The two predictions:

| word | len | predicted | measured death | meas - pred |
|---|---|---|---|---|
| W(10,9)  | 78 | 162.900000 = 1629/10  | 162.8999999999984 | -1.6e-12 |
| W(12,12) | 98 | 166.153846 = 2160/13 | 166.1538461538449 | -1.2e-12 |

```
python3 problems/billiards-triangles/explore/deathlaw_measure.py death --a 10 --b 9  --out problems/billiards-triangles/data/deathlaw_measure_W109.json
python3 problems/billiards-triangles/explore/deathlaw_measure.py death --a 12 --b 12 --out problems/billiards-triangles/data/deathlaw_measure_W1212.json
```

**Both predictions hold to ~1e-12 deg (float measurement).**  Exact
brackets (`deathlaw_exact.py`):

- `alive --a 10 --b 9 --dg 1e-8`: a rational apex (denominator 2^40) whose
  exact Fraction corridor is positive AND whose 78-bounce periodic orbit is
  re-derived by 002's independent exact simulator, with obtuse angle
  CERTIFIED > 162.899999 deg (certified < 162.9 too).  So
  death(W(10,9)) >= 162.899999 exactly.  Same for W(12,12):
  death >= 2160/13 - 1e-6 exactly (files
  `data/deathlaw_exact_alive_W109.json`, `..._W1212.json`).
  The angle certification compares the exact rational
  cos^2(gamma) = (CA.CB)^2/(|CA|^2|CB|^2) against an interval enclosure of
  cos^2 of the rational-degree bound (tier-0 `unfold.py` interval cosine,
  rational pi enclosure) — single evaluations, no iteration.
- `deadscan --dg 1e-4`: 160 (resp. 174) rational apexes with certified
  gamma in [D+5e-5, D+2e-4], concentrated at the last-life window edge:
  every exact corridor empty (max width -5e-5).  Sample-bounded EVIDENCE
  upward — superseded by the theorem below.

### T2 — exact closed forms

`deathlaw_symbolic.py` implements the division-free unfolding
(selftest: ring identities, triangle side lengths, and corridor
alive/dead verdicts against `skeptic_family`'s independent float corridor
on 4 family members x 60 random obtuse triangles — all agree).  For the
half word it extracts the glide axis direction delta (a monomial), the
axis offset m, and all gate-endpoint projections p(z) = Im(conj(delta) z)
as exact trig polynomials.  Structure found (then proven exact in the
ring, `deathlaw_prove.py`):

With p(B), p(C) the projections of gate 1 = [B, C] and p(A_1) of the fan
pivot in gate 2 = [C, A_1], the three binding functions factor:

- **I1**: m - p(C) = -[cos(a alpha) sin(alpha) sin(b beta)
  + sin(a alpha) sin(beta) cos(alpha + (b+1) beta)]
- **I2**: p(A_1) - m = cos(a alpha) sin((b+1) beta) sin(alpha+beta)
- **I3**: p(B) - m = cos((b+1) beta) sin(a alpha) sin(alpha+beta)

(also, not load-bearing: p(A_1) - p(C) = -sin(beta)
sin((a-1)alpha - (b+1)beta) and p(B) - p(C) = sin(alpha)
sin(a alpha - b beta), exact for all 15 members via the same prover —
these are the "collapsing gate" loci that the float layer sees as the
pinch.)

These identities were verified EXACTLY (polynomial subtraction equals the
zero polynomial in the ring) for all 15 members
(2,1),(2,2),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5),(6,6),(7,7),
(8,8),(10,9),(12,12), together with the glide facts (half-word isometry is
orientation-reversing with monomial linear part delta^2, and translation
tau parallel to the axis: Im(conj(delta) tau) = 0 in the ring):

```
python3 problems/billiards-triangles/explore/deathlaw_prove.py all --out problems/billiards-triangles/data/deathlaw_prove_all.json   # ~4 min, "ALL OBLIGATIONS CERTIFIED"
```

**SPECULATION (labelled):** I1-I3 hold for ALL integers a >= b >= 1.  They
are exact for the 15 members above; the general statement is a finite
geometric-series computation over the two mirror fans that I did not carry
out symbolically in (a,b).

### T3 — the theorem

**Setup facts (proven).**  For a translation word w = u^2 with u an
orientation-reversing isometry whose square is a translation, u is a glide
reflection with axis parallel to tau, and gate n'+k of the unfolding is
the u-image of gate k (unfolding functoriality).  A glide reflection acts
on axis-normal projections by p -> 2m - p (elementary; uses the glide
facts verified in the ring).  Hence if the corridor (intersection of all
2n gate projection intervals, in the direction of tau — the harness's
criterion; the rescaling between the harness's normalization and mine is
nonzero real, and every condition used below is invariant under scaling
and sign flip of the projection) has positive width, then for every
half-word gate g:  min_g < m < max_g, strictly — positive width forces
gate g's interval and its glide image [2m - max_g, 2m - min_g] to overlap
in more than a point.  Applied to gates 1 and 2:

either (Case I)  p(C) < m < min(p(A_1), p(B)),
or     (Case II) max(p(A_1), p(B)) < m < p(C).

By I1-I3 these are sign conditions on three explicit trig products.  Write
theta = alpha + beta (= 180 - gamma) and

  theta_d(a,b) := 90 (a+b) / (a (b+1)) deg,
  gamma_d(a,b) := 180 - theta_d(a,b).

**THEOREM (necessity).**  Let a >= b >= 1 be integers, a >= 2,
a <= 2b + 3, and let (a,b) be a member for which I1-I3 and Lemma C below
hold (machine-certified for the 15 members listed; conjectured for all).
If the corridor of W(a,b) has positive width at a triangle with angles
(alpha, beta), then theta > theta_d(a,b), i.e. gamma < gamma_d(a,b).

For b = a this is gamma < 180a/(a+1); for b = a-1 it is
gamma < 180(1 - (2a-1)/(2a^2)) — exactly 001's two SPECULATION laws.

**Proof.**  All angles in degrees; alpha, beta > 0, alpha + beta < 180;
sin(alpha), sin(beta), sin(alpha+beta) > 0.  Note I1 rewrites (exact
product-to-sum, verified in the ring along with I1) as

  m - p(C) = -G,  G := cos(a alpha) sin(alpha) sin(b beta)
                       + sin(a alpha) sin(beta) cos(alpha + (b+1) beta).

*Case II.*  N2 := p(A_1) - m < 0 and N3 := p(B) - m < 0 and m < p(C).
Branch on a alpha mod 360:
 (A1') a alpha in (0,90): cos, sin of a alpha > 0, so N2 < 0 forces
   sin((b+1)beta) < 0, i.e. (b+1)beta > 180, so theta > 180/(b+1)
   >= theta_d  (equivalent to a >= b).  Done.
 (A2') a alpha in (90,180): N3 < 0 forces cos((b+1)beta) < 0 so
   (b+1)beta > 90; with alpha > 90/a: theta > 90/a + 90/(b+1) > theta_d.
 (A3') a alpha in (180,270): N2 < 0 forces sin((b+1)beta) > 0, N3 < 0
   forces cos((b+1)beta) > 0, so (b+1)beta < 90 mod 360; if
   (b+1)beta > 360 then theta > 360/(b+1) >= theta_d (4a >= a+b); else
   (b+1)beta < 90, so b beta < 90 and sin(b beta) > 0.  Then both terms
   of G have cos(a alpha) < 0, sin(a alpha) < 0, so m < p(C) (G > 0)
   forces cos(alpha + (b+1)beta) < 0, i.e. alpha + (b+1)beta > 90, so
   theta > (b alpha + 90)/(b+1) > (b(180/a) + 90)/(b+1) > theta_d
   (the last step is 180b + 90a > 90(a+b)).
 (A4') a alpha in (270,360): theta > alpha > 270/a >= theta_d, which is
   a <= 2b+3.  (A5') a alpha >= 360: a fortiori.

*Case I.*  N1 := m - p(C) > 0, N2 > 0, N3 > 0.  Branch on a alpha:
 (A2) a alpha in (90,180): cos(a alpha) < 0, so N2 > 0 forces
   sin((b+1)beta) < 0, (b+1)beta > 180, theta > 180/(b+1) >= theta_d.
 (A3) a alpha in (180,270): sin(a alpha) < 0, so N3 > 0 forces
   cos((b+1)beta) < 0; and cos(a alpha) < 0, so N2 > 0 forces
   sin((b+1)beta) < 0: together (b+1)beta > 180, done as in (A2).
 (A4/A5) a alpha > 270: theta > 270/a >= theta_d (a <= 2b+3).
 (A1) MAIN CASE: a alpha in (0,90).  sin(a alpha) > 0, so N3 > 0 gives
   cos((b+1)beta) > 0; if (b+1)beta > 270 then
   theta > 270/(b+1) >= theta_d (2a >= b), so assume (b+1)beta in (0,90).
   Then sin((b+1)beta) > 0 (consistent with N2 > 0) and
   sin(b beta) > 0.  Let
     v := 90 - a alpha in (0, 90),   w := alpha + (b+1) beta - 90,
   so that G = sin(v) sin(alpha) sin(b beta) - cos(v) sin(beta) sin(w),
   and N1 > 0 means G < 0.  Elementary algebra:
     a w < b v  <=>  theta < theta_d  (and = iff =).
   Suppose toward contradiction theta <= theta_d, i.e. w <= b v / a.
   - If w <= 0: G >= sin(v) sin(alpha) sin(b beta) > 0, contradiction.
   - Regime 2, v >= 90/(b+1): from (b+1)beta < 90, w < alpha, and
     alpha < 90/a <= 45, so sin(w) < sin(alpha) and
     G > sin(alpha)[sin(v) sin(b beta) - cos(v) sin(beta)]
       = sin(alpha) sin(beta) cos(v) [tan(v) D_b(beta) - 1] >= 0,
     because tan(v) >= tan(90/(b+1)) and D_b(beta) := sin(b beta)/
     sin(beta) >= D_b(90/(b+1)) = cot(90/(b+1)) (Lemma D below), so the
     bracket is >= tan(90/(b+1)) cot(90/(b+1)) - 1 = 0.  Contradiction.
   - Regime 1, v < 90/(b+1): 0 < w <= bv/a <= v < 90, so
     sin(w) <= sin(bv/a) and, with alpha = (90-v)/a,
     G >= sin(beta) cos(v) [tan(v) sin((90-v)/a) D_b(beta)
                            - sin(bv/a)].
     By Lemma D (strict for b >= 2, since beta < 90/(b+1) strictly),
     D_b(beta) > cot(90/(b+1)) for b >= 2, and by LEMMA C below
     tan(v) sin((90-v)/a) >= tan(90/(b+1)) sin(bv/a), so the bracket is
     > 0 (for b = 1, D_1 = 1 and Lemma C is strict on the open interval,
     giving the same strict conclusion).  Contradiction.
   Hence theta > theta_d.  QED (modulo Lemmas C, D).

**LEMMA D** (classical).  D_b(x) = sin(bx)/sin(x) is strictly decreasing
on (0, 90/b); equivalently tan(bx) > b tan(x) for 0 < bx < 90 (induction
on b via the tangent addition formula, using tan A tan B < 1 for
A + B < 90); and D_b(90/(b+1)) = cos(90/(b+1))/sin(90/(b+1)).

**LEMMA C** (the one nontrivial input; per-member certified).  With
v0 = 90/(b+1), for all v in (0, v0):
  H2(v) := sin(v) sin((90-v)/a) cos(v0) - sin(v0) sin(bv/a) cos(v) > 0,
and H2(v0) = 0 identically (because (90-v0)/a = b v0/a — this equality is
exactly what places the death angle at theta_d).  Certified in
`deathlaw_prove.py` by: (i) on (0, d1], the bound sin(bv/a) <=
(b/a) sin(v)/cos(v) (from sin x <= x <= tan x) reduces positivity to one
certified-positive constant bracket; (ii) adaptive interval bisection on
[d1, v0-d2] (43-130 leaves per member); (iii) on [v0-d2, v0], a certified
sign H2' < 0, which with H2(v0) = 0 forces H2 > 0.  All enclosures are
exact-rational interval sin/cos.  **SPECULATION:** Lemma C for all
integers a >= b >= 1 with (a,b) != (1,1) (numerically true for all
1 <= b <= a <= 25 at 4000 samples each; h == 0 identically at a=b=1).

### The general law, and two more out-of-sample confirmations

theta_d(a,b) = 90(a+b)/(a(b+1)) is defined for ALL a >= b, not just
b in {a, a-1}: the theorem predicts death(W(5,3)) = 144 deg exactly and
death(W(4,2)) = 135 deg exactly — angles 001 never measured (it did not
know these members were alive at all).  Measured (`deathlaw_measure.py
death --full`, since these members' alive windows sit on the mirror half
alpha < beta that 001/002 never scanned):

| word | len | predicted | measured | meas - pred | last-alive apex |
|---|---|---|---|---|---|
| W(5,3) | 34 | 144 | 143.9999999999959 | -4.1e-12 | (18.000000, 18.000000) |
| W(4,2) | 26 | 135 | 134.9999999999942 | -5.8e-12 | (22.500000, 22.500000) |

both at the corner (alpha, beta) = (90/a, 90(a-1)/(a(b+1))) the theorem
predicts, and both with exact alive certificates at gamma_d - 1e-8
(`data/deathlaw_exact_alive_W53.json`, `..._W42.json`).  Note
death(W(4,2)) = death(W(3,3)) = 135: two distinct words dying on the same
arc — relevant to 001's pinch-gap lead.

### Death brackets now on record (exact side)

Per-member, combining the theorem (upper bound, exact) with the exact
alive certificates (lower bound, exact rational apex + certified angle +
independently simulated orbit):

  death(W(a,b)) in [gamma_d - 1e-6, gamma_d]   for
  (2,1), (2,2), (3,2), (3,3), (5,3), (4,2), (10,9)   (1e-8 for (3,3),(4,2));
  death(W(12,12)) in [2160/13 - 1e-6, 2160/13].

(Lower-bound constants are the certified per-file values; see
`data/deathlaw_exact_alive_W*.json`.)

## Outcome

- **VERIFIED (T1 — the queue's question):** both out-of-sample predictions
  of 001's SPECULATION laws hold: W(10,9) dies at 162.90 deg and W(12,12)
  at 2160/13 = 166.1538 deg, to ~1.5e-12 deg in the float measurement,
  with exact-rational alive certificates (orbit re-derived by independent
  simulation) at gamma >= gamma_d - 1e-6 and, per the theorem, no corridor
  at gamma >= gamma_d.  The laws are NOT refuted; they are now derived.
- **VERIFIED (per member):** for the 15 members
  (2,1),(2,2),(3,2),(3,3),(4,2),(4,3),(4,4),(5,3),(5,4),(5,5),(6,6),(7,7),
  (8,8),(10,9),(12,12): the identities I1-I3, the glide facts, and Lemma C
  are machine-certified (exact ring arithmetic; exact-rational interval
  enclosures), and hence the necessity THEOREM holds for each: **the
  W(a,b) corridor is empty at every triangle with gamma >= gamma_d(a,b)
  = 180 - 90(a+b)/(a(b+1))** — for every obtuse angle at or beyond the
  law's value, not merely on sampled arcs.  This upgrades 001's death-angle
  laws from measured to established (upper side), and pins the family's
  length-vs-angle growth law l(gamma) = 1440/(180-gamma) - 6 along
  W(a,a): its ordinates are now theorems on the upper side and certified
  to 1e-6 on the lower.
- **VERIFIED (new law):** the two 001 laws are the b = a and b = a-1
  slices of one formula gamma_d(a,b) = 180 - 90(a+b)/(a(b+1)), confirmed
  out-of-sample on W(5,3) -> 144 and W(4,2) -> 135 (predicted before
  measurement, ~5e-12 agreement, death corner exactly as predicted).
- **EVIDENCE:** the float death measurements (adaptive sampling as stated;
  DEAD verdicts sample-bounded); the exact deadscans above gamma_d
  (160/174 apexes, superseded by the theorem).
- **SPECULATION (labelled inline):** (i) I1-I3 for general (a,b) beyond
  the 15 verified members; (ii) Lemma C for general (a,b) (checked
  numerically to a <= 25); with (i)+(ii) the theorem is fully parametric
  for a >= b >= 1, a >= 2, a <= 2b+3.
- **NOT claimed:** that death = gamma_d exactly (the remaining gap is
  one-sided: alive points APPROACHING gamma_d; certified alive is at
  gamma_d - 1e-6, and the limit argument — first-order margins at the
  death corner — is Lead 2, not executed); anything about words outside
  the W(a,b) family, about unstable orbits, or about existence/absence of
  periodic orbits at gamma >= gamma_d (death of THIS corridor construction
  only — a triangle past gamma_d may well have orbits from other words);
  anything about the birth side of the alive windows; Lemma C or I1-I3
  for (a,b) outside the stated scope (in particular a > 2b+3, where the
  branch analysis has a genuine hole and W(6,3)-type members are
  unmeasured).

## Why it failed / what survived

Nothing failed in the headline sense — the predictions survived and the
proof went through — but the honest obstruction ledger:

1. **The sufficiency limit is not closed.**  The theorem is one-sided.
   Proving death = gamma_d exactly needs alive triangles at every
   gamma < gamma_d near gamma_d, i.e. a curve into the alive corner with
   all ~n gate margins positive; each margin vanishes linearly at the
   corner, so this is a first-order computation with an explicit
   second-order remainder — routine but n-gates-wide, and I stopped at
   certified alive points at gamma_d - 1e-6 rather than do it loosely.
2. **The parametric theorem rests on two per-member-verified inputs.**
   I1-I3: the two-fan structure makes the general identity a geometric
   sum any patient hand can do; I verified it exactly member-by-member
   instead (15 members, including both queue targets), because the ring
   check is unconditional while a hand derivation of the general case
   would itself need a skeptic.  Lemma C: certified per member; the
   general v-monotonicity route (log-derivative of
   tan(v) sin((90-v)/a)/sin(bv/a)) looked provable but I did not finish
   it.  Both are cleanly stated, finite, and attackable.
3. **A method lesson worth the ledger:** the death mechanism was invisible
   in apex coordinates (001 called the endpoints "rational functions",
   true but useless) and became three one-line products in the
   circumdiameter normalization where reflections are division-free.
   Choosing the exact domain to kill denominators did more than any
   amount of computing.
4. **What survived for reuse:** the Laurent-ring unfolding engine
   (`deathlaw_symbolic.py`, selftested against 002's independent
   corridor); the glide-reflection reduction (any u^2 translation word:
   corridor = half-word interval I cap (2m - I)) — this halves every
   future corridor computation for doubled words; the certified 1-D
   inequality pattern in `deathlaw_prove.py` (endpoint-algebra + interval
   bisection + derivative-sign closure); the adaptive death-measurement
   tool; and the binding-gate identities I1-I3, D1, D2 themselves.

The skeptic should attack, in order: (a) the glide-to-corridor reduction
step (is "positive width => m strictly inside gates 1 and 2's projection
intervals" airtight against the harness's exact criterion? — it rests on
gate n'+k = u(gate k) plus the ring-verified glide facts); (b) the case
tree of the theorem (a dropped mod-360 branch would be silent); (c) the
interval certification in `deathlaw_prove.py` (a sign error in H2' would
void the near-v0 leg); (d) re-derive I1-I3 for one member with different
software; (e) the claim that measured "death" and corridor-death coincide
(001's family bisection used the same corridor criterion, so this is
definitional, but say so).

## Leads generated

1. **Close the parametric gaps.**  (i) Prove I1-I3 for all (a,b) by the
   two-fan geometric sum (finite hand computation; the fan projections
   are arithmetic progressions of sines — the Dirichlet-kernel forms are
   already visible in the per-member output).  (ii) Prove Lemma C for all
   a >= b >= 1: show Phi(v) = tan(v) sin((90-v)/a)/sin(bv/a) is
   decreasing on (0, 90/(b+1)]; its log-derivative is
   2/sin(2v) - (1/a)cot((90-v)/a) - (b/a)cot(bv/a), and the claim is that
   it is negative there.  Either yields the fully parametric theorem.
2. **Close the sufficiency limit.**  At the death corner the three binding
   margins have explicit nonzero gradients (from I1-I3); all other gate
   margins are bounded below by exact cyclotomic values.  A first-order
   expansion along (alpha, beta) = corner + t*(direction into the wedge)
   with a certified quadratic remainder proves alive for all
   gamma in (gamma_d - eps, gamma_d), giving death = gamma_d exactly.
   Concrete and finite per member.
3. **W(a,b) for a > 2b+3.**  The branch analysis genuinely breaks at
   a = 2b+4 (the a alpha in (270,360) branch can undercut theta_d).  Is
   W(6,3) alive anywhere, and does it die at 146.25 = gamma_d(6,3) or at
   something new?  A --full measurement run answers the first question in
   minutes.
4. **Pinch gaps with the general family.**  death(W(4,2)) = 135 =
   death(W(3,3)) and the W(4,2) tile lives on the mirror half: re-examine
   001's pinch gap [135.000, 135.049] with BOTH halves scanned — queue
   item 11's scan should include mirror-half words before concluding a
   gap is empty.
5. **Birth angles.**  The same machinery (which gate pair binds at the
   birth edge?) should give closed-form birth angles; 001 measured e.g.
   birth(W(4,3)) = 135.0486.  If births also have closed forms, the
   family's alive windows are fully algebraic and the staircase geometry
   of 001 becomes exact.
6. **Kill the corridor enclosure loss instead (001 lead 4).**  I2/I3 show
   the binding functions are smooth products with explicit factors; a
   Taylor-model certifier for JUST these two functions would certify
   family tiles at cell scale rather than 2^-57, directly testable
   against 001's anchor data.

## References

- `problems/billiards-triangles/attempts/001-word-census-coverage-map.md`
  (family definition, measured laws; its "leads generated" item 1 is this
  attempt) and `002-skeptic-review-of-001.md` (out-of-sample rows W(5,4),
  W(7,7); independent tools reused here).  Both were read — mode is
  informed.
- Tier-0: `harness/billiards-triangles/unfold.py` (corridor criterion;
  interval sin/cos used for certifications), `verify_cover.py` (not used;
  orbit re-derivation used 002's `skeptic_orbit.py` instead).
- New code and data: `explore/deathlaw_{measure,exact,symbolic,prove}.py`;
  `data/deathlaw_measure_W*.json` (11 members),
  `data/deathlaw_exact_alive_W*.json` (8 members),
  `data/deathlaw_exact_dead_W{109,1212}.json`,
  `data/deathlaw_prove_all.json` (15 members, all obligations certified).
- No external papers consulted.
