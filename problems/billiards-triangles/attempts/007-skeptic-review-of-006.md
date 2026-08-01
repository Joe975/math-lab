# 007 — Skeptic review of 006 (design hunt past 135°): adversarial verification

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed (read `prior-art.json`, attempts 001, 003, 004 in full,
  002 skimmed via its tools, the four `design_*.py` tools and all
  `design_*.json` data under review, `census.py cmd_family` (001's sampler),
  `deathlaw_skeptic.py` (004's stack, reused here), STATUS queue items 12/14,
  tier-0 `unfold.py` for cross-reference only)
- **Type:** skeptic review of `006-design-family-past-135.md` (default
  stance: REFUTE). Every load-bearing layer re-established by code written
  for this review or by 004's previously-adversarial stack; nothing from
  `design_*.py` is imported in any verdict path.
- **Outcome in one line:** 006 survives where it is exact — all five
  certificates re-established from scratch (third corridor implementation,
  own simulator, widths match digit-for-digit), the pinch-gap mechanism of
  001's error is reproduced and pinned, the birth law holds on three
  genuinely out-of-sample members to ~3e-11, and the (5,2)/(6,2) extension
  re-proves — but four corrections: a universe-count mislabel (811 is
  |u| ≤ 15, the true length-30 count is 566; "1056" double-counts), an
  "alive throughout the gap" wording overreach, the refutation semantics
  against 001 (001's gap claim was knowledge-scoped and stays true; what
  falls is 001's birth *measurements*), and the discovery that 006's
  "window pairing" (a,b) ↔ (b+1, a−1) is a word identity — the paired
  members are the SAME canonical word — so that consistency check carries
  zero evidential weight for the birth law (and quietly halves the member
  count of several lists).
- **Tools:** new `explore/dsk_design_skeptic.py` (stdlib-only,
  deterministic, seed 20260731; ~25 min total compute). Independence
  choices: exact unfolding by composed affine maps with 2x2 rational
  matrices — a third implementation path (006's exact layer is 002's
  `skeptic_orbit.py`, reflected-vertex chains; 004's is complex composed
  maps) with the same tau-normal projection, so exact widths must and do
  agree digit-for-digit across all three; my own exact rational billiard
  simulator (own ray/segment solve); certified gamma brackets and ring
  identities via `deathlaw_skeptic.py` (004's adversarial stack: own
  Machin pi, own Taylor enclosures, own Laurent ring — never 006's code,
  never tier-0 intervals); own word enumeration and canonicalization; own
  seeded samplers; a replica of `census.py cmd_family`'s exact grid design
  for the 001-mechanism question.
- **Sources:** repo only; Niven's theorem cited as classical (as in 006).

My data lives in `data/dsk_*.json` (distinct from 004's `dlsk_*` and the
parallel agent's `plaw_*`; nothing of 006 was modified).

Reproduce everything (repo root):

```
python3 problems/billiards-triangles/explore/dsk_design_skeptic.py selftest
python3 .../dsk_design_skeptic.py certs      --out .../data/dsk_certs.json
python3 .../dsk_design_skeptic.py birthdepth --out .../data/dsk_birthdepth_W43.json
python3 .../dsk_design_skeptic.py censusgap  --out .../data/dsk_censusgap.json
python3 .../dsk_design_skeptic.py birthlaw   --out .../data/dsk_birthlaw.json
python3 .../dsk_design_skeptic.py counts     --out .../data/dsk_counts.json
python3 .../dsk_design_skeptic.py deadsample --n 60 --out .../data/dsk_deadsample.json
python3 .../dsk_design_skeptic.py coverage   --out .../data/dsk_coverage.json
python3 .../dsk_design_skeptic.py newmembers --out .../data/dsk_newmembers.json
python3 .../deathlaw_skeptic.py lemmac --a 5 --b 2 --out .../data/dsk_lemmac_W52.json
python3 .../deathlaw_skeptic.py lemmac --a 6 --b 2 --out .../data/dsk_lemmac_W62.json
# W(6,3)/W(8,2) out-of-sample births: inline driver recorded in data/dsk_birthlaw2.json
```

Selftest: Fagnano closes under my simulator, obtuse Fagnano rejected,
3-way exact-corridor agreement (mine vs `skeptic_orbit` vs
`deathlaw_skeptic`) 100/100 at seeded rational apexes, W(3,3) alive at
134.99 / dead at 135.01 under my float corridor.

## Claims attacked

### 1. The five exact certificates. **CONFIRMED (re-established from scratch)**

`certs` (`dsk_certs.json`) re-verifies each of 006's five certificate
files with: (a) my matrix-form exact corridor positive at the recorded
rational apex — and the exact rational width equal DIGIT-FOR-DIGIT to the
recorded string in all five (same tau-normal normalization, third
implementation); (b) my own exact simulation closing the periodic orbit
(30, 30, 38, 30, 46 bounces; own intersection solve, strict-interior
bounce check, exact closure of position and direction, struck-sequence
check); (c) the word string equal to the claimed family member (W(4,3),
W(4,3), W(5,4), W(5,2), and N1's four-block structure — rebuilt from the
(a,b,c,d) exponents, not copied); (d) for the two bracket pairs, my
gamma brackets re-certified with 004's interval stack (exact rational
cos^2 against Machin-pi Taylor enclosures): gamma in (135.001, 135.048),
(135.0000001, 135.001), (144.0000001, 144.001) all hold — the first two
strictly inside 001's recorded [135.000, 135.0486] gap.

The exactly-135 layer was re-derived by hand before running anything:

- Parametrization: with s = 2(1−t)/(1+t²), x = 1−s/2, y = ts/2 one gets
  (2x−1)² + (2y+1)² = (1−s)² + (ts+1)² = 2 + s(s(1+t²) − 2 + 2t) = 2,
  identically. ✓
- On the circle, 4(x²−x+y²) + 4y = 0, i.e. dot := CA·CB = x²−x+y² = −y,
  automatically negative for y > 0; then |CA|² = x−y, |CB|² = 1−x−y and
  (x−y)(1−x−y) = (x−y) − (x²−y²) = 2y² = 2·dot², so
  cos γ = −y/√(2y²) = −1/√2 exactly, i.e. γ = 135° (γ obtuse from
  dot < 0; the γ = 45° branch of cos² = 1/2 is excluded). So 006's pure
  Fraction criterion "2(CA·CB)² = |CA|²|CB|² and CA·CB < 0" is exactly
  γ = 135°, no interval arithmetic needed. Both recorded apexes satisfy
  it exactly (checked in Fractions, plus the circle identity directly).
- Niven: at the W(5,2) apex tan α = y/x = 512242667139/1686780588413
  (006 quotes the unreduced ratio of numerators — same rational, a nit),
  rational and not in {0, ±1}, so α (hence β = 45° − α, hence both
  non-gamma angles) is an irrational multiple of π. The
  irrational-angled claim stands.

Verdict: **all five certificates are real.** An irrational-angled
triangle with obtuse angle exactly 135° carrying a length-30 stable
periodic orbit is now double-verified.

### 2. The pinch-gap finding vs 001. **CONFIRMED — mechanism identified; semantics corrected (C3)**

Two independent attacks, both in `dsk_birthdepth_W43.json` /
`dsk_censusgap.json`:

- *Is 135.0000076 a sampler floor or a genuine offset?* My own birth
  bisection of W(4,3) with corner accumulation at (22.5, 22.5), at four
  sampler depths: birth − 135 = 7.63e-6 (floor 7.63e-6), 2.98e-8 (floor
  2.98e-8), 1.18e-10 (floor 1.16e-10), 3.6e-12 (floor 4.5e-13, float
  noise regime). The measured "birth" tracks the probe floor at every
  depth — it is the sampler, not the window. Combined with the exact
  gapdeep certificate (an alive point with certified γ < 135.001 and
  > 135.0000001), birth(W(4,3)) ≤ 135.0000001 exactly.
- *Could 001's design have seen it?* I read `census.py cmd_family`: its
  alive test is 400 uniform x-samples on [0.02, 0.5] (largest sample
  x = 0.4994, terminal clearance 6.0e-4), bisected in γ. My replica of
  exactly that design, with MY corridor, reproduces 001's number: birth
  between 135.04857 and 135.04858 (001 recorded 135.0486). The measured
  alive x-window of W(4,3) is a sliver [0.5 − w(γ), 0.5) hanging on the
  corner x = 1/2 (α = β = 22.5), with w = 6.2e-5 at γ = 135.005 growing
  to 6.0e-4 at γ = 135.0486 — exactly the census grid's terminal
  clearance, where the window first swallows a grid point
  (`contains_grid_point` flips true precisely there). **Mechanism
  pinned: 001's uniform grid stops one half-cell short of the corner
  apex the window pinches onto; its "birth" is the γ at which window
  width equals grid clearance.**

Semantics (correction C3): 001 scoped the gap claim as knowledge —
"measured gap ... in which nothing ≤ 30 is *known* alive", births
tabulated as measurements. That claim stays literally true as written;
what 006's certificates falsify is 001's birth *values* (and their
implied "±~0.005° sampling noise", which for births is off by 10x —
the true error at W(4,3) is 0.0486°; 001's own "deaths only ever
under-estimated" hedge silently does not transfer to births, where
finite sampling errs the other way). So "the pinch-gap **refutation**
of 001" (the framing this review was handed) and the record's
"artifact" language are right about the measurements and wrong to the
extent they suggest 001 claimed the gap was empty: 001 claimed only
that nothing was known alive there, and 006 is the attempt that made
something known. "Closed/corrected", not "refuted". 006's own record
mostly says exactly this ("the pinch-gap lead ... is closed"); the
index entry's "001 pinch gap was a sampler artifact" is fair if read
about the measured gap object.

### 3. The birth law γ_birth = 180 − 90(a+b+1)/(a(b+1)). **CONFIRMED out-of-sample; one evidential deflation (C4)**

First the deflation. 006's "window pairing" consequence — (a,b) and
(b+1, a−1) have identical windows because a+b and a(b+1) are invariant —
is a **word identity**: canonical(W(a,b)) = canonical(W(b+1, a−1)) under
cyclic rotation, reversal and the 0↔1 swap. Verified exactly for
(3,3)~(4,2), (5,2)~(3,4), (4,4)~(5,3), (5,5)~(6,4), (7,3)~(4,6),
(2,2)~(3,1); (2,1), (4,3), (5,4) are self-paired. The swap is the
x → 1−x mirror, which fixes γ, so the two labels denote one canonical
word and their γ-windows coincide *for any window quantity whatsoever*.
Consequences: (i) the pairing "check" in 006 §2 is vacuous as evidence
for the law (it is a necessary consistency property of ANY correct
formula, and would equally "check" for many wrong ones); (ii) 003's
"two distinct words dying on the same arc" (W(3,3)/W(4,2)) are not
distinct canonical words; (iii) 006's 11-member float test spans only 9
distinct canonical words, and my first out-of-sample pick W(6,4)
collapsed into fit-member W(5,5) — replaced below.

The law itself, attacked on THREE genuinely out-of-sample canonical
words (none in the fit classes), own sampler, both arc halves, corner
accumulation at (90/a, 90/(b+1)) (`dsk_birthlaw.json`,
`dsk_birthlaw2.json`):

| word | len | predicted birth | measured − predicted | dead at birth − 1e-3 |
|---|---|---|---|---|
| W(7,3) | 42 | 1012.5/7 = 144.642857... | +3.0e-11 | yes (best −4.3e-5) |
| W(6,3) | 38 | 142.5 | +3.0e-11 | yes (best −6.1e-12) |
| W(8,2) | 42 | 138.75 | +3.0e-11 | yes (best −2.2e-11) |

Every +3.0e-11 equals my depth-34 sampler floor (0.5·2⁻³⁴ = 2.9e-11) —
the same floor-tracking signature as W(4,3) at 135. Note W(8,2) has
a > 2b+3, *outside* 003's theorem scope: the birth law holds there
anyway, as 006's s1-screen hits suggested. The W(5,2) window edges also
re-measure to the law: alive at 132.001, dead-sampled at 131.999, alive
at 135 with width 0.69. The touching identity at 144 is carried by the
re-verified W(5,4) certificate. The law remains SPECULATION as a general
statement — 006 labels it so, correctly — but it now survives
out-of-sample tests it was not fitted to, including outside the death
theorem's scope.

### 4. Coverage (157 arcs). **CONFIRMED at my own 25 arcs**

`dsk_coverage.json`: my own arc list (25 rational γ including the
historical trouble spots 112.5, 120, 135, 135.02, 135.0486, 144, 150,
1080/7 = 154.2857..., and 157.1, 160.9, 162.5, 164.8), my own
strict-window member selector (Fraction arithmetic on
a+b < t·a(b+1) < a+b+1), my own float corridor and sampler: **25/25
alive**, longest member length 90. My selector independently lands on
the mirror-labelled members (e.g. (3,4) ≡ W(5,2) at γ = 135), matching
006's b ≤ a−2 story. The coverage conjecture for all of (90, 180)
remains SPECULATION, labelled as such in 006.

### 5. The negative screens. **CONFIRMED with counting corrections (C1)**

Independent recount (`dsk_counts.json`), own generator and
canonicalization: |u| ≤ 13 by half-length 1/2/6/16/52/168, total **245**
— matches; 2-alternating k = 8..12: **2008** — matches; but |u| = 15 is
**566**, not 811: 811 is the count for |u| ≤ 15 (245 + 566; 006's
screen file `design_screen_all15.json` indeed enumerates `--max-half 15`
= everything from 3 to 15, so the *screen's coverage* is fine and its
own n_words = 811 is correct for what it ran). Corrections: the record's
"**|u| = 15 (length 30): 811 canonical words**" mislabels the universe;
"all doubled odd words to length 30: 1056 words" double-counts (245 +
811, but the 245 are inside the 811) — the true number is 811; the index
one_line's "3000+ canonical doubled words to length 50" is likewise
inflated by the same double count (811 + 2008 + 27 three-block = 2846).

Screen verdicts spot-attacked: (a) the s1 hit list reclassified with my
canonicalizer — exactly 14 W members (including W(8,2), W(9,2) with
a > 2b+3, as claimed) plus N1, N2; (b) N1 and N2 match no W(a,b)
canonical form at their lengths (all a+b = 11 resp. 12 checked — and
this check now correctly includes the 0↔1-swap equivalence, the missed-
equivalence worry 006 itself raised); (c) `deadsample`: 60 seeded-random
words of the 245 at 006's 13 arcs plus 6 arcs of my own choosing, own
sampler with random probes added: **zero alive, best width 5.6e-17**
(float zero); the 19 recorded near-misses attacked harder (denser
sampler, ±0.35° arc shifts): the only positive is the W(3,3)≡W(4,2)
canonical word at γ = 134.65 — *below* 135, inside its known window,
exactly as it should be. The length-≤26 stall at γ ≥ 135 survives this
attack. Not re-verified: the full 566-word length-30 sweep, the full
2008-word s1 sweep (only its positives were reclassified), and the
three-block table — all remain 006's sample-bounded EVIDENCE, correctly
labelled.

### 6. The (5,2)/(6,2) extension of 003's theorem. **CONFIRMED (own ring, own intervals)**

`dsk_newmembers.json`, `dsk_lemmac_W{52,62}.json`: scope preconditions
hold (a ≥ b ≥ 1, a ≥ 2, a ≤ 2b+3: 5 ≤ 7, 6 ≤ 7). In 004's independent
Laurent ring (my driver, not 006's `deathlaw_prove` path): all 15
residuals — I1, I2, I3, D1, D2, mu monomial, tau parallel, gate
identifications, glide action on all four base points — are exactly zero
for both members, plus 10 seeded off-torus rational spot points each in
the pair algebra. Case-tree adversarial sign search at θ ≤ θ_d (273k /
277k seeded points, corner- and boundary-accumulating): zero Case I/II
sign-pattern hits. Lemma C re-certified for both members with 004's
interval stack (H2' upper bounds on the endpoint zone: −0.037, −0.032 —
comfortably negative; endpoint identity 90 − v0 = b·v0 exact). With
004's hand-verified general-(a,b) case tree, the necessity theorem
genuinely extends to (5,2) and (6,2): no positive-width corridor at
γ ≥ 138 resp. 140. The death bracket death(W(5,2)) ∈ [135, 138] is
then exact on both sides (alive-at-135 certificate + theorem).

### 7. Scope honesty. **CONFIRMED with corrections C1-C4**

SPECULATION labels are present and correctly placed (general birth law
and exact-window claim in §2; coverage conjecture in §6;
non-cyclotomicity in Outcome); the negative screens are labelled
sample-bounded EVIDENCE with their designs stated; the NOT-claimed list
is accurate and pre-empts the right overclaims (no birth-necessity
theorem, no area coverage, no unstable orbits, non-doubled words ≤ 26
not rescreened). Two wording-level overreaches found: **(C2)** §1's
"length 30 is alive **throughout** the measured gap" — the certificates
prove alive points AT γ = 135 exactly and at one point in each of
(135.0000001, 135.001) and (135.001, 135.048); "throughout" (every γ in
the gap) is float evidence (W(5,2)'s window [132, 138] straddling the
gap) plus the SPECULATION window law, and should be so worded — the
Outcome bullet's "together closing 001's recorded pinch gap" inherits
the same point-vs-continuum slack. And **(C1/C4)** as above. The index
entry's status VERIFIED describes the certificate layer accurately;
its range field states the float/certificate split honestly.

## Refutations found

No load-bearing claim is refuted. Four corrections:

- **C1 (counting).** "|u| = 15 (length 30): 811 canonical words" — the
  length-30 universe is 566; 811 is |u| ≤ 15. "All doubled odd words to
  length 30: 1056 words" double-counts; the true count is 811. The
  index one_line's "3000+" is 2846 counted correctly. (The screens
  themselves covered what they claim to cover; only the counts are
  wrong.)
- **C2 (wording overreach).** "Length 30 is alive throughout the
  measured gap": certificates establish alive points at γ = 135 exactly
  and at one γ in each recorded bracket — pointwise, not throughout; the
  continuum statement is float + SPECULATION law.
- **C3 (refutation semantics).** What is overturned in 001 is its birth
  *measurements* (off by 10x its stated ±0.005° noise; the "deaths only
  under-estimated" hedge does not apply to births) and the emptiness
  *reading* of the gap; 001's actual recorded claim ("nothing ≤ 30
  known alive") was knowledge-scoped and remains true as written. The
  relationship is correction-of-measurement + lead-closure, not
  refutation of a recorded claim.
- **C4 (vacuous consistency check + member double-listing).** The
  window pairing (a,b) ↔ (b+1, a−1) is a canonical-word identity (0↔1
  swap + reversal/rotation; verified exactly on six pairs), so: it is
  not evidence for the birth law; W(3,3)/W(4,2) and W(4,4)/W(5,3) are
  single canonical words, making 006's 11-member fit really 9 canonical
  members and 003's "two distinct words dying on the same arc" a
  relabeling; and "W(5,2), a member NOBODY had measured" is, as a
  canonical word, the mirror W(3,4) — new to the measured record all
  the same.

Nits, no action needed: the W(5,2) tan α is quoted in unreduced form;
§3's width 0.19 (glide-reduced normalization) vs §4's 0.76 (full
corridor) for W(5,2) at 135 are different normalizations of the same
aliveness, unflagged; `design_certify.best_apex` imports
`skeptic_family.width` lazily for non-doubled words but is only ever
run on doubled ones here.

## Claims that survive

| # | 006 claim | Verdict |
|---|-----------|---------|
| 1 | Five exact certificates (W(4,3) x2 in the gap, W(5,4) at 144-touch, W(5,2) and N1 at exactly 135°) | **CONFIRMED** — third-implementation exact corridor, digit-for-digit widths; own exact simulation; own gamma brackets; on-arc identity hand-derived; Niven applies |
| 2 | 001's pinch gap [135.000, 135.0486] not empty; birth values sampler artifacts | **CONFIRMED** (as measurement-correction, C3) — floor-tracking at 4 depths; 001's grid design reproduced and its 135.0486 regenerated; corner-sliver mechanism pinned |
| 3 | Birth law γ_b = 180 − 90(a+b+1)/(a(b+1)) (SPECULATION, float) | **SURVIVES stronger tests than 006 ran** — 3 out-of-sample canonical members incl. a > 2b+3, all +3e-11 = floor; pairing check deflated (C4) |
| 4 | W-family member alive at every sampled obtuse arc | **CONFIRMED** at 25 arcs of my own choosing incl. all trouble spots |
| 5 | Screens: 245 words ≤ 26 dead ≥ 135; only W(4,3)/W(5,2) alive at len 30 on 135-144; s1 hits = W + N1/N2; N1/N2 non-W | **CONFIRMED** modulo counts (C1) — 245 and 2008 recounted; hits reclassified; N1/N2 non-W incl. swap; 60-word random re-test + near-miss attack, zero kills |
| 6 | 003's obligations machine-certified for (5,2), (6,2) | **CONFIRMED** — own ring 15/15 residuals, spots, 550k-point case-tree search, own Lemma C intervals |
| 7 | Scope honesty | **CONFIRMED with C1-C4** |

**Kill attempts that failed,** for the record: (i) hunting a
normalization or transcription error in the certificates by demanding
digit-for-digit width equality from a third implementation — all five
match; (ii) trying to expose 135.0000076 as a genuine birth offset by
quadrupling sampler depth — it tracked my floor at every depth; (iii)
attempting to falsify the birth law on members outside its fit set,
including outside 003's theorem scope — held to the sampler floor all
three times; (iv) sign-pattern search for a case-tree breakdown on the
two new members (550k adversarial points) — nothing; (v) random-
subsample and near-miss re-testing of the dead screens with a different
seeded sampler and extra arcs — the only positive was a known-alive W
word below 135; (vi) coverage at adversarial arcs (touch angles, 001's
accumulation points, γ > 160) — 25/25 alive.

## Residual risk

- **Shared sampler design on the negative side.** My float samplers,
  like 006's, accumulate at 90/j edges (plus uniform and random
  probes). A window pinching onto an interior non-90/j point thinner
  than the uniform/random resolution could hide from both. N1/N2 (whose
  argmaxes ARE interior points, found by the same class of sampler)
  bound this risk but do not eliminate it; all negative verdicts remain
  sample-bounded exactly as labelled.
- **Not re-swept:** the 566-word length-30 universe and 2008-word s1
  universe full negatives (positives were reclassified, negatives only
  subsampled at |u| ≤ 13), the three-block 25/27-dead table, the
  157-arc coverage run (25 arcs re-done), and `design_factor.py`'s
  binding-factorization claims (§5 of 006: the I2 blind re-derivation
  and N1's non-elementary factorization) — the latter feed no VERIFIED
  claim, only the design-space narrative.
- **Common mathematical core.** All three exact corridors implement the
  same corridor definition in the same normalization; agreement rules
  out implementation error, not a shared conceptual error. That risk is
  bounded by 004's earlier tie of this corridor to tier-0
  `unfold.certify` at point boxes and by the simulations (mine and
  002's) re-deriving actual billiard orbits with no unfolding at all.
- **Float birth measurements share the one-sided floor bias** (mine and
  006's both report birth ≈ truth + floor); the exact statements are
  the certificates, which are one-sided (alive at points) — no exact
  birth bracket from below exists for any member except via death(W)
  of the touching predecessor.

## Leads generated

1. **Record the pairing as structure.** Prove
   canonical(W(a,b)) = canonical(W(b+1, a−1)) for all a ≥ 2, b ≥ 1
   (finite word computation per member; the swap+reversal bijection is
   visible in the block structure and six cases are verified here) and
   fold it into the census tooling: quotienting by it halves every
   W-member list and deduplicates future fits. Falsifiable by a single
   canonical-form comparison.
2. **Certify one interior point of the gap continuum.** C2's slack
   closes with one more certificate: W(5,2) at a rational apex with
   certified γ ∈ (135.02, 135.03) (its window is wide there — width
   ~0.69 in float — so `design_certify`-pattern certification is easy).
   Would upgrade "alive throughout the gap" from float to three-point
   + wide-window evidence.
3. **Birth-side exact brackets.** The floor-tracking signature suggests
   certifying alive at γ = γ_b + 1e-8 for one out-of-sample member
   (e.g. W(6,3) at 142.5 + 1e-8, corner (15, 22.5)) — same pattern as
   006's gapdeep certificate; would give the birth law its first exact
   out-of-sample bracket.
4. **The interior-pinch blind spot.** Both 006's and my samplers are
   strongest at 90/j corners. Build one screen pass whose accumulation
   points are the measured interior argmaxes of near-miss words (the 19
   recorded near-misses are the natural test set); if any near-miss
   flips alive, the negative screens' shared-design risk is real and
   every dead verdict needs a wider sampler family.

## References

- `problems/billiards-triangles/attempts/006-design-family-past-135.md`
  (under review) and its tools/data
  `explore/design_{neighbours,factor,certify,exact135}.py`,
  `data/design_*.json`.
- `001-word-census-coverage-map.md` (the corrected birth measurements;
  `explore/census.py cmd_family` read for the mechanism),
  `003-death-angle-laws.md` (laws, identities, Lemma C),
  `004-skeptic-review-of-003.md` (the adversarial stack
  `explore/deathlaw_skeptic.py` reused here for ring/interval layers),
  002's `skeptic_orbit.py`/`skeptic_family.py` (cross-referenced in
  selftest only).
- Tier-0: `harness/billiards-triangles/unfold.py` read as
  cross-reference; not used in any verdict path.
- New code: `explore/dsk_design_skeptic.py`. New data:
  `data/dsk_{certs,birthdepth_W43,censusgap,birthlaw,birthlaw2,counts,deadsample,coverage,newmembers}.json`,
  `data/dsk_lemmac_W{52,62}.json`.
- Niven's theorem — classical, as cited by 006.
