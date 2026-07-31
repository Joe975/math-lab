# 002 — Skeptic review of 001 ({0,±1}^4 census)

- **Problem:** mahler-4d
- **Date:** 2026-07-30
- **Mode:** blind
- **Type:** skeptic review (adversarial re-implementation, default stance REFUTE)
- **Tools:** Python 3 stdlib only, three from-scratch skeptic implementations
  (`problems/mahler-4d/explore/skeptic_burnside.py`, `skeptic_exact.py`,
  `skeptic_hanner.py` — sharing no code with `census.c`, the harness, or the
  attempt's Python tools); `awk`/`cmp` for data-file integrity; the tier-0
  harness `verify_product.py` used once as a third opinion, never as my
  evidence; the attempt's own `hanner_match.py` used only as the *object under
  adversarial test* in claim 5.
- **Sources:** ATTEMPT.md, tier-0 `problems/mahler-4d/PROBLEM.md`, tier-0
  harness, the attempt's `explore/` code (read, not trusted) and `data/`
  (1.2 GB of eval slices, level files, candidates, certificates, claims).
- **Stance:** adversarial; every headline number re-derived independently.

My differences from the record's implementations: Burnside counting and
canonical forms via my own B4-on-pairs action with byte-table mask images;
volume products with normals by 4D/3D generalized cross products (not rref
nullspace), integer-only arithmetic scaled from Q, 3D volumes by
cone-from-centroid (not origin-based divergence), 2D base case by centroid
fan over supporting lines (not 1D max−min recursion), and polar facets
computed two independent ways (brute-force hull of the polar vertex set vs
vertex–facet duality) that were cross-checked against each other on 40 seeded
bodies (`skeptic_exact.py crosscheck`) with exact agreement.

## Claims attacked

**Claim 1 — level completeness (18,637,214 orbits, per-level counts).
Verdict: CONFIRMED.**
My own Burnside count over the 384 signed permutations acting on the 40
antipodal pairs (`skeptic_burnside.py counts 11`) reproduces all eleven
claimed level counts exactly (4, 21, 123, 756, 4310, 22567, 103649, 415920,
1456466, 4478554, 12154844), and the kernel of the action on pairs is exactly
{±I} (order 2). `wc -l data/level*.masks` matches every count. Canonicality
audit with my own canonical form: **every** mask in levels 1–9 (2,003,816
masks) is canonical and the files are duplicate-free — since the number of
distinct canonical forms equals my independently computed orbit count, levels
1–9 are *proven* complete orbit transversals. Levels 10–11 were sampled
(50,000 masks each, seed 20260730): 0 non-canonical, 0 duplicates, counts
match Burnside — strong evidence, not proof, at those two levels
(`data/skeptic/audit_small.txt`, `audit_large.txt`). Separately, the eval
slices were checked to cover each level file exactly, line-for-line in order
(`cmp` of concatenated first fields vs level files, all 8 levels), the
per-level ok/improper/degen counts match the ATTEMPT table exactly
(ok = 1,773,715; degen = 1,988; improper = 16,861,363), and zero `st=punt`
lines exist anywhere.

**Claim 2 — exact near-bound values (1113 at exactly 32/3, 63 above, 0
below, gap 336311/31104). Verdict: CONFIRMED.**
My own exact-rational implementation recomputed **all 1176** bodies in
`data/candidates_all_4_9.txt`: 1113 have P = 32/3 exactly, 63 lie strictly
above, none below; the minimum non-attainer is P = 336311/31104 =
32/3 + 4535/31104 at mask `81101001` (vol = 73/12, vol° = 4607/2592, 10
vertices, 22 facets), matching the record digit-for-digit
(`data/skeptic/exact_all_1176.txt`). Full non-attainer spectrum reproduced:
1 × 336311/31104, 2 × 29165/2688, 2 × 6251/576, 4 × 56399/5184 (all k=5),
54 × 98/9 (all k=6). Attainer counts by vertex count: 517/414/163/19 at
8/10/12/16 vertices. I also independently verified that the candidate set is
the right set: an awk scan of all 61+ eval slices finds exactly 1176 `st=ok`
lines with float P < 10.9, all at k = 4, 5, 6, 8 and none at k = 7, 9, 10,
11, and the mask set is identical to `candidates_all_4_9.txt`. The exact
level minima claimed for the non-Hanner levels were re-derived from my own
argmin scan and recomputed exactly: k=7: 11 (`16b40`), k=9: 799/72
(`28ae50`), k=10: 133/12 (`32c028b000`), k=11: 4775/432 (`72c028b000`)
(`data/skeptic/level_minima.txt`). All 7 stored claim JSONs have vertex sets
consistent with their masks and values matching my recomputation; the gap
body was additionally confirmed by the harness Fubini verifier (third
implementation). My sanity anchors: cube (16, 2/3), cross-polytope (2/3, 16),
24-cell (P = 16) all exact with my code.

**Claim 3 — every attainer is Hanner (10 B4-direct + 1103 GL certificates).
Verdict: CONFIRMED.**
I checked **all 1113 certificates, not a sample** (`skeptic_hanner.py`,
`data/skeptic/hanner_cert_check.txt`): all 10 B4 lines verified with my own
canonical form; all 1103 GL lines have det(T) ≠ 0 (my own Laplace expansion
over Fractions) and T mapping the mask's vertex set *bijectively* onto the
target Hanner vertex set (4 distinct GL targets). Targets were validated
against **my own** Hanner generator built from the definition (segments under
ℓ¹/ℓ∞ sums): it yields exactly 22 vertex sets in R^4 in 10 B4-orbits with
vertex counts {8, 10, 12, 16}, matching the harness family the record used
(note: that is 22 vertex sets / 10 orbits / 4 combinatorial types — the
record correctly never claims "8 combinatorial types"). Every one of my
Hanner masks has P = 32/3 by my own exact code. The set of 1113 attainer
masks in the certificate file is *identical* to the set my claim-2
recomputation found at exactly 32/3, and attainers occur only at Hanner
vertex counts.

**Claim 4 — float-filter soundness. Verdict: CONFIRMED.**
Stratified seeded sample of **240 bodies** (20 per slice from 12 slices
covering every level k = 4..11, seed 20260730), each recomputed exactly with
my own code and compared against the kernel's printed float P: maximum
|float − exact| = **8.9e-15** (`data/skeptic/float_vs_exact.txt`), with
nv and nf also matching on every sampled line. No deviation approached 1e-6,
let alone the 0.23 margin between the 10.9 cutoff and 32/3 ≈ 10.667; no body
with float P ≥ 10.9 and exact P < 10.9 was found. The record's stated error
bounds (~1e-13 relative kernel accumulation, ~1e-11 in the Outcome) are
conservative relative to my measurements. I additionally attacked the
bookkeeping that the filter rests on: 40 sampled `st=improper` masks are all
genuinely improper by my own vertex test (tight-normal rank), 8 sampled
`st=degen` masks all genuinely have rank < 4, and the record's two improper
spot-check examples reproduce exactly (`207f`→P=32/3 = its extreme submask
`205a`; `20bf`→P=221/18 = `20ba`), so skipping impropers hides nothing.

**Claim 5 — the fixed hanner_match cannot be fooled the way the disclosed
mid-run bug was. Verdict: CONFIRMED.**
Code inspection: the fixed `gl_equivalent` has two independent guards — an
explicit rank-4 check on T and a set-*equality* (not subset) test
`img == Wset` — either of which kills the singular-map false positive.
Adversarial probes running *their* code (`data/skeptic/adversarial_probe.txt`):
the gap body `81101001` (10 vertices — a Hanner-compatible count), the 98/9
body `1680180` (12 vertices), and a 16-vertex non-attainer `21b7` (P = 40/3)
are all classified NOT-HANNER by the exhaustive basis-image search, while a
genuine attainer (`11001010`) still produces a valid certificate. The search
logic is exhaustive as claimed: any linear equivalence maps 4 independent
basis pairs of V to 4 distinct pairs of W up to signs, and ±T normalization
loses nothing.

**Claim 6 — scope honesty. Verdict: CONFIRMED (one minor wording flag).**
The k ≤ 11 / ≤ 22-vertex truncation, the {0,±1} coordinate restriction, and
the two-tier float-then-exact design are all stated explicitly and early; the
"NOT claimed" paragraph is accurate and appropriately broad; `VERIFIED` is
used only for the 1176 individually recomputed bodies "in the harness's
sense" while the census-scope statement is labeled EVIDENCE, per the tier-0
verification contract; SPECULATION is labeled inline where used (Leads §3).
The universe-equivalence sentence ("every centrally symmetric polytope in R^4
with at most 22 vertices, all in {0,±1}^4") is mathematically correct: the
vertex set of a symmetric body is symmetric, 0 cannot be a vertex, and B4
reduction preserves P. Internal arithmetic (18,637,214 = Σ levels;
18,637,066 = k ≥ 4; 1,773,715 + 1,988 + 16,861,363 partition) all checks out.

## Refutations found

None that touch any claim. One wording inconsistency worth a corrective note
in the index rather than an edit (per the no-edit rule):

- Outcome §4 says the level minima beyond k = 8 sit "without a monotone
  trend (799/72 ≈ 11.097 > 133/12 ≈ 11.083 > 4775/432 ≈ 11.053)" — but the
  parenthetical itself exhibits a strictly monotone *decrease* over
  k = 9, 10, 11, and the "Why it failed" section correctly calls the same
  numbers "slightly DECREASING toward 11". The sentence is only true if it
  is read as including k = 7 (minimum exactly 11, below all three). The
  values themselves are correct (I recomputed all four exactly); only the
  prose is self-contradictory. Severity: cosmetic.

Also noted, not a refutation: level-10/11 canonicality is sampled evidence
(50,000 masks each), not the proof-grade full audit I achieved for k ≤ 9;
the record's own Burnside cross-check has the same character. Anyone wanting
proof-grade transversal completeness at k = 10, 11 should run the full audit
(~12 min of compute with `skeptic_burnside.py audit`).

## Claims that survive

All six attacked claims survive, and therefore the headline survives:

> Over all conv(±S), S = k antipodal pairs of nonzero {0,±1}^4 points,
> 1 ≤ k ≤ 11, up to B4 (18,637,214 orbits; 1,773,715 proper bodies), the
> minimum volume product is exactly 32/3, attained by exactly 1113 orbits,
> every one GL-equivalent to a Hanner polytope; the nearest non-attainer is
> 336311/31104 = 32/3 + 4535/31104; nothing lies below the bound.

| # | Claim | Verdict | Independent basis |
|---|-------|---------|-------------------|
| 1 | Level completeness, 11 level counts | CONFIRMED | my Burnside counts match all 11; full canonicality audit k≤9, sampled k=10,11; slice coverage exact |
| 2 | 1113 = 32/3, 63 above, 0 below, gap 336311/31104 | CONFIRMED | all 1176 recomputed with my own exact code; spectrum identical |
| 3 | All attainers Hanner via certificates | CONFIRMED | all 1113 certificates re-verified; my own Hanner generator (22 sets / 10 orbits) |
| 4 | Float filter cannot hide a sub-bound body | CONFIRMED | 240 stratified samples, max dev 8.9e-15 vs 0.23 margin; improper/degen flags verified |
| 5 | Fixed matcher rejects non-Hanner bodies | CONFIRMED | adversarial probes at 10/12/16 vertices all NOT-HANNER; dual guards in code |
| 6 | Scope stated honestly, EVIDENCE-language | CONFIRMED | read-through vs tier-0 contract; one cosmetic wording flag above |

Scope of this review: everything sampled is listed with its seed
(20260730 unless stated; crosschecks 424242, 31337); everything else named
above was checked exhaustively. My review compute stayed within ~25 minutes.
Skeptic artifacts: code in `problems/mahler-4d/explore/skeptic_*.py`, outputs
in `problems/mahler-4d/data/skeptic/`. This working copy is not a git
repository, so nothing was committed here; no existing attempt file or data
file was modified.

## References

- `ATTEMPT.md` (the record under review, root of this working copy).
- Tier-0 `problems/mahler-4d/PROBLEM.md` (statement, verification contract).
- Tier-0 harness `harness/mahler-4d/polytope.py`,
  `harness/mahler-4d/verify_product.py` (selftests re-run; Fubini verifier
  used once as third opinion on the gap body).
- Attempt tooling under review: `problems/mahler-4d/explore/census.c`,
  `aggregate.py`, `exact_check.py`, `hanner_match.py`, `burnside_check.py`;
  data under `problems/mahler-4d/data/`.
- Skeptic re-implementations (this review):
  `problems/mahler-4d/explore/skeptic_burnside.py`,
  `problems/mahler-4d/explore/skeptic_exact.py`,
  `problems/mahler-4d/explore/skeptic_hanner.py`; outputs in
  `problems/mahler-4d/data/skeptic/` (burnside_counts.txt, audit_small.txt,
  audit_large.txt, exact_all_1176.txt, hanner_cert_check.txt,
  float_vs_exact.txt, improper_degen_check.txt, level_minima.txt,
  adversarial_probe.txt, near_masks_myscan.txt).
