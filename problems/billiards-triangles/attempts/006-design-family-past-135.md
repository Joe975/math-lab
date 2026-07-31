# 006 — Design hunt past 135°: the stall dissolves — birth law, W(5,2) alive at exactly 135°, a new four-block family

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed (read `prior-art.json`; attempts 003 and 004 in full;
  001 in full; skimmed 002's tools; queue item 12 of `STATUS.md` — the
  exploratory track of the two-track plan)
- **Type:** computational search (structure-class enumeration + float
  screening) + design inversion of 003's symbolic machinery + exact
  certificates
- **Tools:** new `explore/design_neighbours.py` (enumeration of the glide
  structure class, numeric glide-reduced corridor, adaptive window
  measurement), `explore/design_factor.py` (exact elementary-factor
  extraction over 003's division-free Laurent ring, for arbitrary doubled
  words), `explore/design_certify.py` (exact alive certificate with
  certified rational gamma bracket, any word), `explore/design_exact135.py`
  (exact certificate ON the gamma = 135 arc via its rational-point
  parametrization).  Reused: `deathlaw_symbolic.py` (the substrate; selftest
  run first), `deathlaw_prove.py` (`member` mode, new members),
  `skeptic_orbit.py` / `skeptic_family.py` (002's independent exact/float
  corridor and simulator — every float hit and every certificate is
  cross-checked against them), `deathlaw_exact.py` (gamma certification
  helpers; tier-0 `unfold.py` interval cosine underneath).  All stdlib-only,
  deterministic (fixed seeds where any randomness exists); total compute
  ≈ 50 min.  New data: `data/design_*.json` (16 files).
- **Sources:** none external (Niven's theorem cited as classical).

**Conventions** (001/003): W(a,b) = (0 (12)^a (02)^b)^2, length 4(a+b)+2;
side s opposite vertex s; alpha, beta the angles at A, B; gamma at C;
"alive at a triangle" = the word's unfolding corridor has positive width
there; window = the set of gamma with an alive apex on the constant-gamma
arc.  003 (skeptic-confirmed by 004) proved per-member:
death gamma_d(a,b) = 180 − 90(a+b)/(a(b+1)).

## Approach

Queue item 12 asks for a hunt through half-word/gate structures for
families alive past 135°, with the death-law machinery run in reverse.
Before hunting, the item demands the frontier semantics be pinned down,
because the W family itself has death angles far above 135° (certified
alive members to gamma ≈ 159.25° in 001).  Reading 001's coverage
staircase precisely, the "census stall at 135°" is three separate facts:

- **(F1) length stall:** no word of length ≤ 26 was known alive at any
  sampled point with gamma ≥ 135° (sample-bounded);
- **(F2) pinch gap:** the measured gap [135.000°, 135.049°] between
  death(W(3,3)) = death(W(4,2)) = 135 and measured birth(W(4,3)) =
  135.0486, in which nothing of any length was known alive;
- **(F3) window coverage:** above 135° the known alive set was the union
  of measured W-windows [birth, death] — with the birth values known only
  from 001's fixed 400-point arc scans.

So "a family past 135°" must mean: alive in the pinch gap, or alive at
gamma ≥ 135 with length ≤ 26, or covering gamma-arcs the W family does
not.  That framing directed the work, and it is where the surprise came
from: **(F2) and most of (F3) turn out to be artifacts of 001's sampler.**

Structure class searched (and why it is the natural one): 003's glide
reduction applies to w = u^2 with u an odd word.  A small structural
observation makes this class self-contained: the linear part of any
orientation-reversing isometry composed of an odd number of reflections
is a unit mu with mu·conj(mu) = 1, so u^2 ALWAYS has trivial rotation
part — every doubled odd word is automatically a translation word (no
integer condition to check; verified numerically in the selftest).  The
class searched is: all doubled odd half words |u| ≤ 15 (word length
≤ 30), the "2-alternating" subclass u = y_0 y_1 2 y_2 2 ... y_k 2
(y_i ∈ {0,1}; W(a,b) is y = 0 1^a 0^b) up to k = 12 (length ≤ 50), and
the three-block family u = 0(12)^a(02)^b(12)^c exhaustively to
a,b,c ≤ 3.  Why this rather than the full even-word universe at length
28–34 (queue 14's job): the glide reduction halves the corridor cost,
membership is parametrized (so verdicts become statements about named
families, not word lists), and the symbolic factor machinery applies.

## What was done

All commands from the repo root.  Selftests first:
`deathlaw_symbolic.py selftest` (PASSED), then

```
python3 problems/billiards-triangles/explore/design_neighbours.py selftest
python3 problems/billiards-triangles/explore/design_factor.py selftest
```

The first checks: doubled odd words have translation unfoldings; the
numeric glide-reduced width times |tau| equals `skeptic_family`'s
independent full-corridor width (rel. err < 5e-14 on 600 samples, verdict
agreement 600/600); W(3,3) and W(4,2) windows reproduce 001/003 values.
The second re-derives 003's identity I2 for W(3,3) by blind factor
extraction: p(A_1) − m = sin(4b) cos(3a) sin(a+b) × const, exactly (the
factor list is found numerically but each division is EXACT in the ring,
and the reconstruction product is re-verified by exact multiplication).

### 1. The pinch gap (F2) is not a gap: W(4,3)'s window touches 135

`design_neighbours.py window` measures birth/death with a sampler that
adds geometric accumulation at every candidate window edge alpha = 90/j,
beta = 90/j (001/002/003 used uniform arc grids, plus — in 003 — the two
known W-window edges).  Result: birth(W(4,3)) measures 135.0000076, and
the residual 7.6e-6 is exactly the sampler's deepest edge offset
(0.5·2^-16 deg): the alive window pinches onto the CORNER
(alpha, beta) = (22.5, 22.5) and 001's 400-point scan could not see it.
Deep corner accumulation finds W(4,3) float-alive at gamma = 135.04,
135.02, ..., down to 135.000001 (width ~ 0.34·(gamma−135), confirmed at
every point by `skeptic_family`'s independent corridor).  Exact:

```
python3 .../design_certify.py --word 012121212020202012121212020202 \
    --gamma 135.02  --corner-alpha 22.5 --glo 135001/1000 --ghi 135048/1000 \
    --out .../data/design_cert_W43_pinchgap.json
python3 .../design_certify.py --word 012121212020202012121212020202 \
    --gamma 135.0001 --corner-alpha 22.5 --glo 1350000001/10000000 \
    --ghi 135001/1000 --bits 64 --out .../data/design_cert_W43_gapdeep.json
```

Both PASS: exact Fraction corridor positive at a rational apex, the
30-bounce periodic orbit re-derived by 002's independent exact simulator,
and gamma certified in (135.001, 135.048) resp. (135.0000001, 135.001) —
strictly inside 001's recorded gap.  **The pinch-gap lead (001 lead 2,
STATUS queue 14's premise) is closed: length 30 is alive throughout the
measured gap.**

### 2. Closed-form birth law (003 Lead 5, executed)

`design_factor.py bind` on W(4,3) near the birth corner shows the binding
difference there is exactly 003's I2 with a = 4: cos(a alpha) ·
sin((b+1) beta) · sin(alpha+beta) — the SAME factor cos(a alpha) that
sets the death edge alpha = 90/a also pinches the birth, now jointly with
I3's cos((b+1) beta) at beta = 90/(b+1).  The birth corner is
(90/a, 90/(b+1)), giving

  theta_b(a,b) = 90/a + 90/(b+1),
  **gamma_birth(a,b) = 180 − 90 (a+b+1) / (a (b+1))**,

so the window is [gamma_birth, gamma_d], of gamma-length 90/(a(b+1)).
Float test at 11 members (a,b) = (2,1),(2,2),(3,2),(3,3),(4,2),(4,3),
(4,4),(5,3),(5,4),(5,5),(6,6): alive at gamma_birth + 1e-6 (widths
3e-8..1e-7), dead (sampled) at gamma_birth − 1e-3, every member.
Adaptive bisection on W(5,2) hits the law to 1.5e-11 deg.  Consequences,
each checked:

- **Touching identity:** gamma_birth(a+1,a) = 180a/(a+1) =
  gamma_d(a,a) — successive staircase windows TOUCH (with zero width) at
  the nice angles rather than leaving gaps; 001's measured birth values
  (135.0486, 127.5064, ...) were all sampler-biased upward.  Verified in
  float at the 120, 135, 144, 150 touch points, and exactly (certified
  bracket (144.0000001, 144.001)) for W(5,4) at the 144 touch:
  `design_cert_W54_touch144.json`.
- **Window pairing:** a+b and a(b+1) are invariant under
  (a,b) → (b+1, a−1), so members pair up with IDENTICAL windows:
  W(3,3)/W(4,2) share [127.5, 135], W(4,4)/W(5,3) share [139.5, 144] —
  exactly the double-death coincidences 003 noticed.

**SPECULATION (labelled):** the birth law for general (a,b), and that the
window is exactly the open interval (gamma_birth, gamma_d).  What is
exact here: the two certificates above, plus (necessity side) 003's Case
I branch A1 already contains theta < 90/a + 90/(b+1) as its standing
hypothesis-set — but the other branches are not excluded here, so no
birth-necessity theorem is claimed.

### 3. The screens: what the structure class contains

```
python3 .../design_neighbours.py screen --mode all --max-half 13 \
    --arcs 135.0,135.02,136.5,138.0,140.0,140.7,142.5,144.0,145.0,148.0,151.0,155.0,159.0 \
    --out .../data/design_screen_all13.json
python3 .../design_neighbours.py screen --mode all --max-half 15 \
    --arcs 135.0,135.02,136.5,138.0,140.0,140.7,142.5,144.0 \
    --out .../data/design_screen_all15.json
python3 .../design_neighbours.py screen --mode s1 --klo 8 --khi 12 \
    --arcs 135.0,141.0,144.0,145.5,147.0,148.5,150.0,151.5,154.0,156.0,158.0 \
    --out .../data/design_screen_s1_k8_12.json
python3 .../design_neighbours.py triple --amax 3 --out .../data/design_triple_table.json
```

Every arc-hit is cross-checked against `skeptic_family`'s independent
corridor (all agree).  Results:

- **|u| ≤ 13 (length ≤ 26): 245 canonical words, ZERO alive at all 13
  arcs ≥ 135.**  The length-26 stall (F1) survives this class, now with
  edge-accumulating sampling (sample-bounded EVIDENCE, design as stated).
- **|u| = 15 (length 30): 811 canonical words, exactly TWO alive on
  135–144:** W(4,3) — and **W(5,2)**, a member NOBODY had measured
  (001 scanned b ≥ a−1; 003 added b = a−2; W(5,2) is b = a−3), alive AT
  gamma = 135.0 with large width 0.19.  Its window measures
  [132.0000000000, 137.99999979] — the law's prediction [132, 138]
  straddling 135, dissolving the "exceptional touch arc" at 135.
- **2-alternating class, k = 8..12 (length 34–50): 2008 canonical words,
  16 alive**, of which 14 are W(a,b) members (including a > 2b+3 members
  like W(8,2), W(9,2) — outside 003's theorem scope, alive as predicted
  by the laws) and **TWO are genuinely new four-block structures** (not
  W(a,b) under rotation/reversal/0↔1-swap, checked against all (a,b)
  with matching length):
    N1: u = 0(12)^3(02)^3(12)^2(02)^3, length 46, window
        [134.25, 136.8265] — alive AT 135.0, width 0.082;
    N2: u = 0(12)^3(02)^3(12)^2(02)^4, length 50, window
        [140.264, 141.396].
- **Three-block family 0(12)^a(02)^b(12)^c, all a,b,c ≤ 3: 25/27 members
  dead everywhere obtuse (sampled)**; the two alive are (3,1,1) with
  window [90, 105.0000] and (3,2,1) with [109.25, 111.42] — the extra
  (12)^c block collapses the death angle far below the c = 0 family.

### 4. Exact certificates AT gamma = 135 exactly

The apexes with gamma = 135 form the arc of (2x−1)^2 + (2y+1)^2 = 2,
y > 0 — a circle with rational point (1,0), hence DENSE rational points,
parametrized by rational chord slope.  At such an apex, gamma = 135
exactly is a pure Fraction identity (2·(CA·CB)^2 = |CA|^2|CB|^2 with
CA·CB < 0) — no interval arithmetic at all.

```
python3 .../design_exact135.py --word 012121212120202012121212120202 \
    --out .../data/design_exact135_W52.json
python3 .../design_exact135.py \
    --word 0121212020202121202020201212120202021212020202 \
    --out .../data/design_exact135_N1.json
```

Both PASS: **W(5,2) (length 30, exact corridor width ≈ 0.76) and N1
(length 46, width ≈ 0.60) have exact-arithmetic stable periodic orbits,
independently re-simulated (30 resp. 46 bounces), at rational-apex
triangles whose obtuse angle is EXACTLY 135°.**  At the W(5,2) apex,
tan(alpha) = 300824218724645284407543/990593882979969795299081 is
rational and not in {0, 1}, so by Niven's theorem alpha (hence beta =
45° − alpha) is an irrational multiple of pi: this is an irrational-
angled obtuse triangle at exactly 135°.  Note the trick needs cos^2
(gamma) rational, so it works at touch angles 120, 135, 150 but NOT at
144 (cos^2 144° = (3+√5)/8 — no rational apexes on that arc; hence the
bracket-style certificate for the 144 touch in §2).

Necessity side for the new members: 003's prover run on (5,2) and (6,2)

```
python3 .../deathlaw_prove.py member --a 5 --b 2 --out .../data/design_prove_W52.json
python3 .../deathlaw_prove.py member --a 6 --b 2 --out .../data/design_prove_W62.json
```

ALL OBLIGATIONS CERTIFIED for both — so 003's theorem now covers them: no
positive-width corridor at gamma ≥ 138 (W(5,2)) resp. 140 (W(6,2)).
Death brackets: death(W(5,2)) ∈ [135, 138] with exact alive at 135 and
float death 138 − 2e-6; W(6,2) window [135.0000076, 140.0] (float)
touches 135 from above like W(4,3).

### 5. Binding factorization of the new family, and why its angles are dirty

`design_factor.py bind` on N1 near its death (gamma = 136.8, argmax
alpha = 26.47 — an INTERIOR point, not a 90/j corner):

The binding differences do NOT factor into elementary sin/cos terms: the
leading one is sin(beta)·sin(alpha)·[an 11-term cosine sum]
(`design_bind_N1_death.json`).  So the four-block family's death angle is
a root of a genuine multi-term trigonometric polynomial — consistent with
the measured dirty values 136.8265, 141.3960 — while W's death/birth
angles are clean rational multiples of 90° precisely because its two-fan
gate structure makes every binding difference a THREE-FACTOR product of
elementary terms.  This is the design-space answer to "invert the
machinery": the inversion is clean exactly on the class where the
factorization stays elementary, and that class, in everything searched
here, is the W family itself.

### 6. Coverage: the W family covers every sampled obtuse angle

With windows [gamma_birth, gamma_d] in closed form, member selection
becomes arithmetic: gamma is strictly inside W(a,b)'s window iff
a+b < t·a(b+1) < a+b+1 for t = (180−gamma)/90.  Check: for 157 arcs
(gamma = 90.5 to 165 step 0.5, plus 112.5, 120, 130, 135, 144, 150,
154.3), pick the shortest such member and test aliveness in float
(`design_coverage_check.json`): **all 157 alive, zero failures**, longest
member needed: length 94.  Every "accumulation angle" of 001 — including
112.5, 120, 135, 144, 150 — is strictly inside the window of some W(a,b)
with b ≤ a−2 (mirror-half members 001 never scanned).  **SPECULATION
(labelled):** the union of open windows covers all of (90, 180); by the
window arithmetic this reduces to an elementary Diophantine statement
(lead 2 below).

## Outcome

- **VERIFIED (exact arithmetic, independently re-simulated orbits):**
  five certificates, each an exact-Fraction positive corridor at an
  explicit rational apex plus an independent exact orbit simulation:
  (i) W(4,3) at certified gamma ∈ (135.001, 135.048);
  (ii) W(4,3) at certified gamma ∈ (135.0000001, 135.001) — together
  closing 001's recorded pinch gap;
  (iii) W(5,4) at certified gamma ∈ (144.0000001, 144.001) — the 144
  touch;
  (iv) W(5,2), length 30, at gamma = EXACTLY 135° (exact on-arc
  identity; irrational-angled by Niven);
  (v) N1 = (0(12)^3(02)^3(12)^2(02)^3)^2, length 46, a word outside the
  W family, at gamma = EXACTLY 135°.
- **VERIFIED (003's machine-certified obligations, two new members):**
  I1–I3 + glide facts + Lemma C for (5,2) and (6,2), hence no positive
  width at gamma ≥ 138 resp. 140.
- **EVIDENCE (float, sample-bounded by the stated samplers/arc lists):**
  the screens (245 words ≤ 26: nothing alive at 13 arcs ≥ 135; 811
  length-30 words: only W(4,3), W(5,2) on 135–144; 2008 words of the
  2-alternating class 34–50: only W members plus N1, N2); the three-block
  table (25/27 dead everywhere); the birth law at 11 members; the window
  measurements of W(5,2), W(6,2), N1, N2; the 157-arc coverage check.
- **SPECULATION (labelled inline):** the general birth law
  gamma_birth(a,b) = 180 − 90(a+b+1)/(a(b+1)) and exact-window claim;
  full coverage of (90, 180) by the union of windows; non-cyclotomicity
  of the four-block death angles.
- **NOT claimed:** any birth-side necessity theorem (only 003's death
  side is proven, now for 17 members); that no length ≤ 26 word is alive
  past 135 (float screen, finite arcs, doubled words only — non-doubled
  even words of length ≤ 26 were screened by 001, not here); area/box
  coverage (all claims are about arcs and points; 001's area-collapse
  obstruction stands); anything about unstable orbits; aliveness of
  every W(a,b) in its formula window (tested only at the 157 chosen
  arcs + measured members).

## Why it failed / what survived

The queue's kill condition ("every family in the class has gamma_d ≤
135°, evidence that 135° is a real barrier") is REFUTED rather than
confirmed — but not by a new family: by the discovery that the frontier
semantics themselves were wrong.  Specifically:

1. **What was wrong:** 001's pinch gap and birth values were artifacts of
   uniform arc sampling; windows pinch onto corners (90/a, 90/(b+1)) and
   need geometric edge accumulation to see.  The staircase windows touch;
   b ≤ a−2 members interleave and cover the touch points; nothing about
   135° obstructs the W family, let alone the class.  The 135° stall
   reduces to a pure LENGTH statement: nothing ≤ 26 found alive past 135
   (that part survives everything thrown at it here, and is sharp:
   length 30 is alive at 135 with corridor width 0.76 — not marginal).
2. **The obstruction that survived, precisely:** in the searched
   neighbourhood of W's structure (three-block: 27 members; four-block
   and beyond within the 2-alternating class to length 50: 2008 words;
   all doubled odd words to length 30: 1056 words), aliveness in the
   obtuse region is RARE — everything alive past 135 is either W(a,b) or
   the two four-block words N1, N2, whose binding differences do not
   factor into elementary terms and whose windows are short with dirty
   algebraic endpoints.  Clean angle laws — and hence design-by-formula —
   appear confined to the two-fan W geometry.  The design tool works
   (it re-derives I2 blind, and diagnoses N1's non-factorability), but
   the design SPACE is nearly empty near W.
3. **What survived for reuse:** the closed-form window law (birth AND
   death) with the (a,b) ↔ (b+1, a−1) pairing; the coverage-member
   selector (arithmetic, no search); the edge-accumulating sampler
   (any future window measurement must use it — uniform arc grids
   under-estimate every pinching window); the rational-points-on-arc
   certification trick for any gamma with rational cos^2; the arbitrary-
   word certifier and the exact factor-extraction tool; the observation
   that doubled odd words are automatically translation words (kills the
   integer translation check in this class, and makes "length ≡ 2 mod 4"
   the natural census axis for stable structures).

A skeptic should attack, in order: (a) the five certificates (re-run
them; check the on-arc identity algebra and that `verify_orbit`'s
positivity + simulation really are independent of the float layer);
(b) the birth-law fit — is 135.0000076 really the sampler floor and not
a genuine offset (re-measure W(4,3) birth with a deeper sampler);
(c) the screens' negative verdicts (rerun with different arc lists and
denser/other samplers — a thin window between my arcs would be missed);
(d) the canonical-form check that N1, N2 are not W members under some
missed equivalence; (e) the claim that 001's samplers could not have
seen the corner windows (read `census.py cmd_family`'s grid).

## Leads generated

1. **Prove the birth law.**  Necessity candidate: alive + Case I branch
   A1 forces theta < 90/a + 90/(b+1) already in 003's case tree; what is
   missing is excluding the other branches below theta_b, plus a
   sufficiency corner expansion at (90/a, 90/(b+1)) exactly like 003's
   Lead 2 at the death corner.  Same machinery, finite per member.
   Falsifiable: certify alive at gamma_birth + 1e-8 and dead-scan below
   gamma_birth − epsilon for several members.
2. **Coverage theorem.**  Prove: for every t ∈ (0,1) there exist
   integers a ≥ 1, b ≥ 1 with a+b < t·a(b+1) < a+b+1 (then, assuming the
   window law, every obtuse gamma has an alive W member; with the b=1
   and b=2 rows alone the union already covers t ∈ (1/3, 1) up to
   boundary points).  Elementary number theory; the needed member length
   at gamma is governed by the continued-fraction structure of t — an
   explicit, provable refinement of 001's hyperbolic length law.
3. **Certify the irrational-cos^2 touch arcs.**  At gamma = 144 there
   are no rational apexes; certify a tier-0 harness box (TRUE verdict)
   whose interior provably straddles the arc — needs box half-width
   ~2^-40 at length 38 and corner gamma brackets tighter than the box's
   gamma variation.  Concrete and finite.
4. **The length frontier at 135.**  Is any word of length ≤ 26 alive at
   some gamma ≥ 135?  The doubled-odd class is now screened (negative);
   the remaining universe is non-doubled even words at 24, 26 — queue
   14's letter-statistics scan, now pointed at a sharp question with a
   known answer at length 30 (both W(4,3) and W(5,2)).
5. **Four-block systematics.**  N1's corridor at 135 is ~8x wider than
   W(4,3)'s anywhere in the old gap; short dirty-angle windows with FAT
   corridors are exactly what area-coverage (001's real obstruction)
   needs.  Enumerate 0(12)^a(02)^b(12)^c(02)^d with a+b+c+d ≤ 14,
   measure window x corridor-area, and test whether four-block tiles
   fill 001's mid-region holes at 92–135° that W tiles leave.
6. **Update the coverage staircase.**  With the window law and W(5,2)-
   type members, recompute 001's minimal-length-per-arc staircase from
   the formula (min 4(a+b)+2 subject to strict window containment) and
   compare against the census data — any arc where the census beats the
   formula signals a non-W short word worth naming.

## References

- `problems/billiards-triangles/attempts/003-death-angle-laws.md` (the
  machinery inverted here; its Leads 4, 5 are executed by this attempt)
  and `004-skeptic-review-of-003.md` (scope of the glide reduction; the
  C2 zero-width caveat is why "touch" points need the strict-interior
  members).
- `problems/billiards-triangles/attempts/001-word-census-coverage-map.md`
  (frontier semantics; the pinch-gap lead closed here) and 002's tools
  (`skeptic_orbit.py`, `skeptic_family.py` — the independent layer under
  every certificate).
- Tier-0: `harness/billiards-triangles/unfold.py` (interval cosine used
  in the gamma-bracket certifications, via `deathlaw_exact.py`).
- New code: `explore/design_{neighbours,factor,certify,exact135}.py`.
  New data: `data/design_screen_{all13,all15,s1_k8_12}.json`,
  `data/design_triple_table.json`, `data/design_windows_new.json`,
  `data/design_window_W43.json`,
  `data/design_cert_W43_{pinchgap,gapdeep}.json`,
  `data/design_cert_W54_touch144.json`,
  `data/design_exact135_{W52,N1}.json`,
  `data/design_bind_{W43_birth,N1_death}.json`,
  `data/design_prove_W{52,62}.json`, `data/design_coverage_check.json`.
- Niven's theorem (rational multiples of pi with rational tangent are
  exactly those with tan ∈ {0, ±1}) — classical; used only for the
  "irrational-angled" remark, not load-bearing.
