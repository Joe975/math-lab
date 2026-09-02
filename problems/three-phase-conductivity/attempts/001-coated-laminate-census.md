# 001 — Coated-laminate ground truth, Milton's attainability rule, and a bounded-rank laminate census of the gap region

- **Problem:** Optimal three-phase conducting composites in 2D, `problems/three-phase-conductivity/PROBLEM.md`
- **Date:** 2026-09-02
- **Mode:** informed
  (The problem had zero prior attempts, so there was no route-specific prior
  art to read. The director session read `STATUS.md` (the queue line for this
  problem: "first attempt is the two-phase ground-truth self-test, run
  blind") and `docs/CYCLE.md` before starting, which are tier 1, so the
  honest label is `informed`; treat it as effectively blind for the
  blind-vs-informed dataset. The four delegated workers were briefed only on
  tier-0 material plus this attempt's own scripts and never read
  `attempts/`.)
- **Type:** computational search + exact construction + literature skeptic
- **Tools:** all standard-library Python, run with the repo venv.
  `explore/tp_coated.py` (exact coated-laminate constructions, `--selftest`,
  `--scan`), `explore/tp_search.py` (seeded float screen over bounded-rank
  laminate trees), `explore/tp_certify.py` (exact rationalisation of a
  screened tree + eigenvalue enclosures + `verify_laminate.py` pass),
  `explore/tp_shapes.py` (closed forms of the two winning shapes, `--selftest`,
  `--report`). Data: `data/screen/` (float screens + `summary.csv`),
  `data/certified/` (27 exact records, every one passes
  `harness/three-phase-conductivity/verify_laminate.py`),
  `data/shapes-notes.md`, `data/literature-check.md`.
  Every screen is deterministic (seed 1); rank-3 screens take ~150 s per
  grid point, rank-4 (axis normals only) ~15 min. Test:
  `tests/test_three_phase_certify.py`.
- **Sources:** Milton 1981 (Appl. Phys. A 26, 125–130) — *not reached*, cited
  through Albin–Cherkaev–Nesi 2007 (J. Mech. Phys. Solids 55, 1513–1553) and
  the Cherkaev 2009 preprint (arXiv:1009.3060 / Mech. Mater. 41, 411–433), both
  [T]; Nesi 1995 (Proc. Roy. Soc. Edinburgh A 125, 1219–1239) — *not reached*.
  All quotations and the transcribed Cherkaev bound are in
  `data/literature-check.md`, marked [T].

## Approach

Three steps, each falsifiable and each exact where the contract asks for it.

1. **Ground truth.** Rank-2 coated laminates must attain the two-phase
   Hashin–Shtrikman (HS) bounds as equalities in ℚ. If the harness or my
   coating algebra is wrong, this is where it shows.
2. **Extend the coating to three phases and find where it breaks.** The
   obvious extension (one rank-1 core of phases 2|3 coated by phase 1)
   never reaches HS; the Milton-type assemblage (each better phase coated
   separately, core fractions chosen so every coated composite has the same
   conductivity, then mixed) does, and its feasibility condition is a clean
   inequality. Derive it, test it exactly on grids, and compare it with the
   literature.
3. **Where the assemblage is infeasible, ask what bounded-rank laminates do.**
   An exhaustive seeded float screen over every rank-3 laminate tree (leaf
   labels × normals from {e₁, e₂, (1,1), (1,−1)} × continuous fractions) and
   over rank-4 trees with axis normals, followed by exact certification of the
   winners and closed-form analysis of their shapes.

Why this rather than transcribing the improved bounds first: the contract
says bounds must be re-derived or cross-checked before anything is killed
against them, and the harness's HS bounds are the only ones already trusted.
Building structures first gives certified *upper* estimates of the attainable
minimum that any later bound must respect; the literature pass then tells us
which published bound they should be compared with.

## What was done

### 1. Two-phase ground truth (`tp_coated.py --selftest`) — exact

Rank-1 coating identity, checked exactly against `laminate.laminate()` on
every use: for a core tensor T (fraction c) coated by the isotropic matrix sI
in integer normal n,

    (σ' − sI)⁻¹ = (1/c)(T − sI)⁻¹ + ((1−c)/(c s)) n nᵀ/|n|².

Iterating in e₁ then e₂ (core fractions c₁, c₂, total F = c₁c₂) gives
(σ* − sI)⁻¹ = (1/F)(T − sI)⁻¹ + ((1−F)/(F s)) diag(ρ, 1−ρ) with
ρ = (1−c₁)/(1−F). Isotropy is *linear* in ρ, so around any diagonal core the
isotropic rank-2 coated laminate is a rational object (`iso_coat`).

Statement tested: for σ ∈ {(1,3), (2/7, 11), (1, 1000)} and core fraction
f ∈ {1/10, 1/3, 1/2, 7/9, 99/100}, the rank-2 coated laminate with the worse
phase as matrix equals the HS lower bound exactly, and with the better phase
as matrix equals the HS upper bound exactly (`t[0] == lo`, tensor isotropic,
Keller–Dykhne identity holds). All 30 equalities hold in ℚ.

### 2. Three-phase: the naive coat fails, the Milton assemblage works

`tp_coated.py --scan-naive`: one rank-1 core of phases 2|3 (normal e₁),
coated by phase 1 in e₁ then e₂ with ρ chosen for isotropy. At every point of
the k/8 grid with σ = (1,2,5) it is either infeasible (ρ < 0) or strictly
above HS_lo (e.g. 3.4113 vs 3.3636 at f = (1,1,6)/8). Same on the upper
side with phase 3 as coat.

`tp_coated.py --scan` (`milton()`): coat phases 2 and 3 *separately* by phase
1 with core fractions φ₂, φ₃ chosen so both coated laminates have the same
scalar conductivity v, then laminate the two (any normal; identical tensors).
Summing the two-phase HS identities with weights w_i = f_i/φ_i shows v is
exactly the three-phase HS lower bound, and feasibility is φ_i ≤ 1, i.e.
v ≤ σ_i for both cores, i.e.

    HS_lo(f, σ) ≤ σ₂   (lower bound attainable),
    HS_hi(f, σ) ≥ σ₂   (upper bound attainable, phase 3 as coat).

Tested exactly: at every interior grid point for (σ, N) ∈ {((1,2,5), 8),
((1,3,20), 40), ((1/7, 5/3, 101/2), 30)} the construction is feasible iff the
inequality holds (assertion in `scan_milton`), and when feasible the exact
tensor equals HS to the last digit. Boundary case φ₂ = 1 (phase 2 uncoated)
occurs when HS_lo = σ₂, e.g. f = (3,2,3)/8, σ = (1,2,5), where the rank-3
screen below also finds the exact value 2 (`data/certified/r3_3-8_2-8_3-8_lower`).

**Literature (`data/literature-check.md`).** ACN 2007 and Cherkaev 2009 both
attribute to Milton 1981 the condition [T] m₁ ≥ 2Θ(1 − m₂),
Θ = k₁(k₃−k₂)/((k₂+k₁)(k₃−k₁)), for attainability of the HS lower bound.
That is the same condition: substituting f₃ = 1 − f₁ − f₂ into
HS_lo ≤ σ₂ ⇔ f₁/(2σ₁) + f₂/(σ₁+σ₂) + f₃/(σ₁+σ₃) ≥ 1/(σ₁+σ₂) makes it linear
in f₁ and rearranges to f₁ ≥ (1−f₂)·2σ₁(σ₃−σ₂)/((σ₁+σ₂)(σ₃−σ₁)) = 2Θ(1−f₂).
Cross-checked on 19 924 random exact instances with zero disagreements. So
step 2 is a **rediscovery of Milton 1981** (via secondary sources; the
primary paper was not reached).

### 3. Bounded-rank laminate census in the gap region

`tp_search.py`: for fixed target fractions and side, enumerate every binary
tree with R internal nodes, every leaf labelling by {p1,p2,p3}, every
assignment of normals from {(1,0),(0,1),(1,1),(1,−1)} (rank 3: 11 520
topologies; rank 4 with axis normals: 33 600), optimise node fractions by
Nelder–Mead (2 restarts, seed 1) on ±trace/2 + 200·(anisotropy + fraction
error), keep results with both residuals < 1e−6.

`tp_certify.py`: per-phase fractions are multilinear in node fractions;
fix rank−2 of them at `limit_denominator(10⁶)` rationalisations and solve
the remaining two exactly through an (ancestor, cherry) pair (method proved
in its docstring). Then exact tensor, eigenvalue enclosures by bisection to
width < 2⁻⁴⁰, exact HS bounds at the exact fractions, and a
`verify_laminate.py` pass. All 27 screened trees solved exactly; all 27
records verify. Isotropy residuals are ≤ 2⁻³⁶ except the rank-4 record
(≤ 1.3e−8); gaps below are the certified enclosure of λ − HS.

Headline numbers, σ = (1,2,5), lower side, f = (1,1,6)/8 (HS_lo = 37/11):

| structure | rank | σ* | σ* − HS_lo |
|---|---|---|---|
| naive single-core coat (exact) | 3 | 3.4113 | 4.8e−2 |
| nested coat: σ₂ coats σ₃, σ₁ coats that (exact closed form) | 4 | 3.4074 | 4.4e−2 |
| best rank-3 tree, all normals (certified) | 3 | 3.38480017 | 2.116e−2 |
| best rank-4 tree, axis normals (certified) | 4 | 3.36422 | 5.8e−4 |
| same shape, exact 1-parameter optimisation (`tp_shapes.py`) | 4 | 3.3641819563 | 5.456e−4 |
| Milton assemblage | 5 | infeasible (HS_lo > σ₂) | — |

**Rank-3 winner** (`tp_shapes.py`, corrected from the screen's printout,
whose `normals` tuple is in DFS-preorder node order): L12 = p1|p2 rank-1 in
e₁; L13 = p1|p3 rank-1 in e₂; result = L12 | L13 laminated in e₁. With the
fraction constraints m_L = 1 − f₂/m_top, m_R = 1 − f₃/(1−m_top), isotropy in
the sole free parameter m_top is a **quadratic**, 2917x² − 3374x + 520 = 0 at
this point, discriminant 2²·3⁴·61·269 with 61·269 = 16409 squarefree, so the
optimal m_top = (1687 − 9√16409)/2917 ≈ 0.18310642 is irrational and the
certified record is a rational enclosure, not an equality (as the contract
allows). The field computation in `data/shapes-notes.md` §2 shows phases 2
and 3 each see one uniform field while phase 1 is split across two branches
with two different fields — the assemblage's requirement that phase 1 coat
both others uniformly is exactly what is infeasible here, and this shape
never imposes it.

Along f₂:f₃ = 1:6 this rank-3 shape gets *worse* relative to HS as f₁
shrinks (gap 0.021 → 0.038 → 0.051 at f₁ = 1/8, 1/16, 1/32) and has **no
isotropic member at all** below f₁* ∈ (0.03007757, 0.03007758).

**Rank-4 winner**: X = p3|p1 rank-1 in e₂; Z = p1|p3 rank-1 in e₁ (a
different mixture, not the same core coated twice); Y = Z | p2 in e₂;
result = X | Y in e₁. Two free parameters after the fraction constraints,
one after isotropy; golden-section over the exact closed form with exact
bisection for the isotropic point gives 3.3641819563, gap 5.456e−4. The gap
is real for this shape, not an optimiser artefact.

**Wide screen** (`data/screen/summary.csv`, `summary.md`; 33 screened
points). Controls behave: at f = (3,2,3)/8, where HS_lo = σ₂ = 2 is Milton's
boundary, the rank-3 screen returns exactly 2 and certification confirms an
exact hit. Milton-attainable does not imply rank-3-attainable: at
f = (6,1,1)/8 lower the assemblage attains HS_lo but rank 3 stalls 8.1e−3
above it, and on the upper side at (1,1,6)/8 rank 3 stalls 2.0e−2 below HS_hi.
Along every ray tried the gap shrinks monotonically as f₁ grows toward
Milton's threshold, on both sides. Larger contrast means larger gap: at
σ = (1,3,20) the rank-3 gaps are about an order of magnitude larger than at
(1,2,5) for the same fractions. Rank 4 cuts the gap 40-60× at f₁ = 1/8
(gaps 5.8e−4, 4.6e−4, 3.3e−4 at f = (1,1,6), (1,2,5), (1,3,4)/8) but
**returned no improvement at f₁ = 2/8**, reproducing the rank-3 value to every
printed digit at (2,1,5)/8 and (2,2,4)/8; since the rank-3 winner uses axis
normals the rank-4 axis space contains it, so that is a failure to improve
rather than a contradiction. SPECULATION: two-restart non-convergence on a
harder landscape rather than a feature of the family; falsifiable by a rank-4
rerun at (2,2,4)/8 with more restarts and a different seed.

### 4. Which bound is the right comparison

Cherkaev 2009, Theorem 7.1 (transcribed [T] in `data/literature-check.md`,
evaluated by an independent Python transcription and cross-checked by
maximising his one-parameter family on a 200 001-point grid): the bound
improves on HS only for m₁ < m₁₁ = 2√m₂(1−√m₂)·Θ; above m₁₁ it "degenerates
into the Hashin–Shtrikman bound" [T]. At σ = (1,2,5): m₁₁ = 0.1143 for
m₂ = 1/8 and 0.1250 for m₂ = 1/4, so both f = (1,1,6)/8 and (2,2,4)/8 lie in
the HS region of that bound, whereas Milton's threshold there is
2Θ(1−m₂) = 0.4375 and 0.375. **The interval m₁₁ ≤ m₁ < 2Θ(1−m₂) is where the
best published lower bound we could evaluate is plain HS and the classical
assemblage cannot reach it.** Nesi 1995's bound could not be transcribed
(paywalled); Cherkaev states it is "generally not achievable by a
structure" [T]. Whether Cherkaev 2009 or ACN 2007 exhibit a structure
attaining HS on that interval was not settled by the literature pass (the
ACN six-parameter family was not unpacked; Cherkaev–Zhang 2011 PDF failed
to extract).

## Outcome

- `VERIFIED` (scope: the listed instances, exact in ℚ, two harness routes):
  rank-2 coated laminates attain both two-phase HS bounds; the Milton
  assemblage attains the three-phase HS lower/upper bound exactly whenever
  HS_lo ≤ σ₂ / HS_hi ≥ σ₂, on the three grids above; 27 certified laminate
  records with stated eigenvalue enclosures.
- `VERIFIED` (algebra, plus 19 924 exact random instances): the condition
  HS_lo ≤ σ₂ is identical to Milton's m₁ ≥ 2Θ(1−m₂). Rediscovery.
- `EVIDENCE` (scope: rank ≤ 3 with the four-normal set, rank 4 with axis
  normals, σ = (1,2,5) and (1,3,20), grid k/8, seed 1, 2 restarts): the best
  isotropic laminate values per grid point in `summary.csv`; at f = (1,1,6)/8
  the rank-4 value 3.36418 sits 5.5e−4 above HS_lo, forty times closer than
  rank 3.
- `MAP`: the interval m₁₁ ≤ m₁ < 2Θ(1−m₂) as the place where bound (HS) and
  known construction (Milton) do not meet, per the sources reached.

**Not claimed.** Nothing about all microstructures: a screen that stalls
above HS is evidence about that rank, normal set and optimiser only. No
claim that HS is or is not attainable on the interval. No claim about Nesi's
bound. The rank-4 screen used axis normals only and two restarts, so its
value is an upper estimate of the rank-4 optimum. The Milton attribution
rests on two secondary sources.

## Why it failed / what survived

Nothing in the exact part failed. The census part stops at an obstruction
worth naming: **bounded-rank laminates with axis normals close in on HS from
above in the gap region but each shape has a finite feasibility floor in
f₁** (the rank-3 winner loses isotropy below f₁* ≈ 0.0301 at f₂:f₃ = 1:6),
so "does some finite rank attain HS on the interval" is not answerable by
this screen. SPECULATION (analyst, labelled): the 40× drop from rank 3 to
rank 4 is consistent with either convergence to HS or to a strictly larger
floor; two data points cannot distinguish them.

Reusable: the rank-1 coating identity and `iso_coat` (isotropic coating of
any diagonal core is rational); `milton()` as an exact HS-attaining
construction generator; `tp_search.py` + `tp_certify.py` as a
screen-then-certify pipeline for any laminate family; the closed forms in
`tp_shapes.py`; Cherkaev's Theorem 7.1 transcription with its region
thresholds m₁₁, m₁₂.

## Leads generated

1. **Rank 5 and 6 on the "split the worst phase" family** at f = (1,1,6)/8,
   σ = (1,2,5): if the gap ratio per rank stays near 40× the family converges
   to HS (then hunt the exact structure: a value equal to 37/11 in ℚ or a
   proof it is a limit); if the gap stalls or the shape loses feasibility,
   record the floor. Definite either way; ~1 h of screen at rank 5 with axis
   normals.
2. **Does a published structure attain HS on m₁₁ ≤ m₁ < 2Θ(1−m₂)?** Unpack
   ACN 2007's six-parameter family (their Fig. 6) and Cherkaev 2009 §8 and
   evaluate at (1,1,6)/8. If yes, our rank-4 value must be ≥ it and the lead
   1 target is known; if no, lead 1 is genuinely open.
3. **Transcribe Nesi 1995's bound** (needs the PDF) and evaluate at the
   census points; compare with the certified laminate values. A certified
   laminate value below a transcribed bound is a transcription error, and a
   large margin says which bound is tighter.
4. **Skeptic pass on this record**: re-implement `milton()` from the
   coated-sphere picture instead of the coating identity; re-derive the
   quadratic for the rank-3 shape; re-run one rank-3 screen with a different
   optimiser (e.g. coordinate descent) and seed, and check the same topology
   wins.
5. **Feasibility floor map**: for the rank-3 winner, compute f₁*(f₂:f₃) as a
   curve and check whether it coincides with Cherkaev's m₁₂ or m₁₁ threshold
   (falsifiable by evaluation).

## References

- G. W. Milton, Appl. Phys. A 26 (1981) 125–130. Not reached; condition
  quoted via [T] ACN 2007 and Cherkaev 2009.
- P. Albin, A. Cherkaev, V. Nesi, J. Mech. Phys. Solids 55 (2007) 1513–1553. [T]
- A. Cherkaev, Mech. Mater. 41 (2009) 411–433; preprint arXiv:1009.3060. [T]
- V. Nesi, Proc. Roy. Soc. Edinburgh A 125 (1995) 1219–1239. Not reached.
- Harness: `harness/three-phase-conductivity/laminate.py`, `verify_laminate.py`.
