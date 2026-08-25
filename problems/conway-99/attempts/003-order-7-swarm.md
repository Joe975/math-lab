# 003 — The order-7 case has no orbit matrix, and a local-model swarm as a third worker family

- **Problem:** conway-99, `problems/conway-99/PROBLEM.md`
- **Date:** 2026-08-25
- **Mode:** informed. Read `prior-art.json`, `PRIOR-ART.md`, 002 in full and
  001's lemma section (L1–L4). No `blind.sh` checkout was used.
- **Type:** computational elimination (orbit matrices for an assumed
  automorphism of order 7), plus the harness work that made it run, plus a
  calibration of a new bulk-worker family.
- **Target shape:** object — an orbit matrix, then a 99-vertex graph, either
  of which a checker settles instantly. **Push rounds: 4.** The first search
  did not reach row 5 in 20 million nodes; each round below added one
  structurally different prune, and the fourth (a repaired symmetry break)
  is what finished it.
- **Tools:** all standard-library Python 3.13, deterministic. New, under
  `problems/conway-99/explore/`:
  - `orbit7.py` — orbit-matrix enumerator with row-inner-product pruning,
    lex-leader symmetry breaking, an optional trace filter and an optional
    eigenvalue-projector PSD/rank prune. CLI documented in the docstring.
  - `char_blocks.py` — exact arithmetic in Q(ζ_p); verifies the character-
    block multiplicity lemma below on four real graphs.
  - `swarm-003/make_briefs.py` and the briefs it writes — the worker prompts
    (provenance for every swarm return quoted here).
  - `tests/test_conway99_orbit7.py` — 19 tests, ~6 s (`MATHLAB_SLOW=1` runs
    the headline computation to exhaustion, ~30 s).
  - Harness change (tier 0, methodological): `scripts/swarm.py` gained a
    `local` provider for an OpenAI-compatible llama-server on this machine.
- **Workers:** `Qwen3.8-27B-UD-Q4_K_XL` (unsloth GGUF) on llama-server, port
  1234, two slots, thinking **off** (`--effort minimal`), 12 verification
  briefs at 12k tokens and 18 ideation briefs at 3k. Every return is a
  candidate; what survived the director's filter is stated per item below.
- **Sources:** `problems/conway-99/PROBLEM.md` only, machine transcribed
  `[T]`. No papers read. In particular Behbahani–Lam (2011) and
  Cesarz–Woldar (2023), which PROBLEM.md summarises as leaving order 7 open,
  were **not** consulted; see "What is not claimed".

## Approach

002 left five leads. Its lead 3 — *order 7 is the smallest search this
problem admits and remains untouched* — is the one whose whole case split is
already pinned by proved-or-nearly-proved facts: 001's L3 gives exactly one
fixed vertex, so the orbit profile is (1, 7¹⁴) and two full rows of the
orbit matrix are forced before any search. Its lead 1 — *find the t ≥ 4
successor of the forced-autocorrelation trick* — turned out to have a
partial answer that feeds straight into lead 3: decomposing the adjacency
matrix over the characters of Z₇ ties the eigenvalue multiplicities of the
quotient to those of the graph, which fixes the trace of the orbit matrix.

Why orbit matrices rather than the alternatives: a SAT encoding of the
lifted graph is what the published attack does and reports as hard; an
annealer has no positive control at this size (002's ceiling table); and
the orbit-matrix step is the one part of the assumed-symmetry route that is
*complete in principle* — 002 closed Z₃₃ with it. The question was only
whether the 15-orbit profile could be made to terminate. It could, but not
by faster search: by three exact constraints and one correct symmetry break.

The swarm was used the way `docs/SWARM.md` prescribes — drafts of
independent implementations, skeptic passes on the lemmas the case rests
on, an ideation sweep — with one change of provider: a local 27B model,
free and private, to see whether it can carry the breadth tier at all.

## What was done

### 1. The order-7 profile, and what is forced before searching

Let τ have order 7. By 001's L3, Fix(τ) = {v₀}; the other 98 vertices form
14 orbits O₁..O₁₄ of size 7. N(v₀) is τ-invariant of size 14 with no fixed
vertex, so it is two whole orbits, O₁ ∪ O₂. A matching edge inside O₁ would
have 7 distinct images inside O₁, i.e. 14 endpoint incidences on 7 vertices
of matching-degree 1; so every matching edge joins O₁ to O₂. Hence

    r₀₁ = r₀₂ = 7,  r₁₀ = r₂₀ = 1,  r₀ⱼ = rⱼ₀ = 0 (j ≥ 3),
    r₁₁ = r₂₂ = 0,  r₁₂ = r₂₁ = 1,  r₀₀ = 0,

and among the 7-orbits R is symmetric with even diagonal (S_ii = −S_ii ⊆
Z₇ \ {0}). The (S), (R), (Q) conditions are as in 002. This is the seed
`orbit7.z7_profile()`; it is the same seed `orbit_general.py` used.

### 2. Rows 1 and 2 are forced up to relabelling (by hand and by machine)

Row 1's free entries (columns 3..14) have sum 12 and, from (Q) at (1,1),
sum of squares 18; so the multiset is {2³,1⁶,0³} or {3,1⁹,0²}. Row 2 must
satisfy Σⱼ r₁ⱼr₂ⱼ = 6 (from (Q) at (1,2): 14 − r₁₂ − r₁₀r₀₂ = 6) with the
same sum and sum of squares. For {3,1⁹,0²} the three cases r₂₃ ∈ {0,1,2}
each force the two columns where r₁ = 0 to carry sum ≥ 6 with squares > 18
— no row 2 exists. For {2³,1⁶,0³}, writing a, b, c for row 2's sums over
the columns where row 1 is 2, 1, 0: 2a + b = 6, a + b + c = 12 and the
square budget kill everything except a = 0, b = 6, c = 6 with all entries
in the b-block equal to 1 and the c-block equal to (2,2,2). So, relabelling
orbits 3..14:

    row 1 = (1, 0, 1 | 2,2,2, 1,1,1,1,1,1, 0,0,0)
    row 2 = (1, 1, 0 | 0,0,0, 1,1,1,1,1,1, 2,2,2)

The enumerator reproduces this: with `--dump-rows 3` exactly one (row 1,
row 2) pair reaches row 3 (`tests/test_conway99_orbit7.py`, slow branch).
The 12 remaining orbits fall into classes A (3 orbits, r₁ = 2, r₂ = 0), B
(6 orbits, r₁ = r₂ = 1) and C (3 orbits, r₁ = 0, r₂ = 2); every row j ≥ 3
then has Σ_{l≥3} r_jl = 12, and (Q) at (1,j) and (2,j) give, for a vertex
in any of the 12 orbits, x neighbours in the A-orbits, 12 − 2x in the
B-orbits and x in the C-orbits.

### 3. The enumerator, and why 002's did not terminate

`orbit7.py` fills rows in order. Because (S) makes the column below the
diagonal a function of the rows above it, *once rows a and i are both
complete the (Q) equation at (a,i) is fully determined* — so while row i
is being filled, each earlier row contributes a linear constraint with
known non-negative coefficients on row i's free entries, and the diagonal
(Q) equation is a fixed weighted sum of squares. Every one of these is
bound-checked after every assignment (suffix maxima). `orbit_general.py`
checked (Q) only when a row was complete, which is why it never reached a
first solution at 15 orbits.

Round 1 (this pruning alone) still did not reach row 5 in 20 M nodes, for a
reason that only became visible later: the ~19 000 arrangements of row 1
over columns 3..14 are pure relabelling, and nothing was rejecting them.

### 4. Round 2 — the character-block multiplicity lemma, and the trace

**Lemma (character blocks).** Let τ be an automorphism of prime order p of
a graph with adjacency A, with t orbits of size p and f fixed points. Then
for every eigenvalue θ,

    mult_A(θ) = a_θ + (p − 1)·r_θ,

where a_θ is the multiplicity of θ in the quotient (orbit) matrix R and r_θ
is the multiplicity of θ in the t×t matrix M_ij = Σ_{m : o_i ~ τ^m o_j} ζ^m
over Q(ζ_p), the same for every non-trivial character.

*Proof.* C^n splits into the isotypic components of ⟨τ⟩. A commutes with τ
so preserves each. On the trivial component (dimension t + f, spanned by
orbit indicator vectors) A acts as the quotient, whose spectrum is that of
R. On the χ = ζ^u component (dimension t, spanned by f_i = Σ_d ζ^{−ud}
e_{τ^d o_i}) a direct computation gives A f_i = Σ_j M^{(u)}_{ji} f_j, and
M^{(u)} is the entrywise Galois conjugate of M^{(1)}, so has the same
rational-eigenvalue multiplicities. ∎

Checked exactly in `char_blocks.py` (rank over Q(ζ₃) by Gaussian
elimination) on Paley(9)/Z₃, rook(3)/Z₃ (semiregular) and rook(4)/Z₃,
Petersen/Z₃ (one fixed point each): all twelve (θ, graph) identities hold
with r_θ equal across the two non-trivial characters.

For the 99-graph with p = 7: 54 = a + 6r where a is the multiplicity of 3
in the 15×15 quotient, so **a ∈ {0, 6, 12}** and tr R = 14 + 3a − 4(14 − a)
= 7a − 42 ∈ {−42, 0, 42}; with tr R ≥ 0, **tr R ∈ {0, 42}**. (Q) alone only
gives tr R ≡ 0 (mod 7) and the parity rule makes it even, i.e. tr R ∈
{0, 14, 28, 42, 56, 70}; the lemma kills four of the six.

### 5. Round 3 — the eigenvalue projectors

Write B = D^{1/2} R D^{−1/2} (D = diag of orbit sizes), symmetric, with
B = 14ww^T + 3P₃ − 4P₋₄ for w_i = √(n_i/99) and P₃, P₋₄ orthogonal
projectors of ranks a and 14 − a. Then P₋₄ = (3I + 11ww^T − B)/7 has
diagonal entry (34/9 − r_ii)/7 on a 7-orbit, which is ≥ 0 only if
r_ii ≤ 3, i.e. r_ii ∈ {0, 2}. In the a = 12 case tr R = 42 needs Σ r_ii =
42 over twelve orbits with r_ii ≤ 2 — impossible. Hence **a = 6, tr R = 0,
every r_ii = 0**, and the 2×2 minors of P₃ (diagonal 30/77, off-diagonal
(r_ij − 14/11)/7) give **r_ij ≤ 4** for every pair of 7-orbits.

`orbit7.py --spectral a` checks, on every completed leading block, that the
integer-scaled congruent forms of both projectors are PSD with rank ≤ a and
≤ t − 1 − a (exact LDLᵀ over Fractions). Calibration: with `a` set from the
real quotient's rank, the real orbit matrices of rook(4)/Z₃, Petersen/Z₃
and Paley(9)/Z₃ survive, and with the wrong `a` they are cut.

### 6. Round 4 — the symmetry break, and the bug it hid

Orbits with the same size and the same seed row are interchangeable.
`orbit7.py` requires R ≥ P R Pᵀ in row-major order for every adjacent
transposition P of two such orbits, comparing only positions already
decided and stopping undecided otherwise — the lex-maximal representative
of every class satisfies this, so no class is lost, and `canonical` (an
exact canonical form by refinement and cell permutation) dedupes the rest.

The first version of this looked right and changed nothing: the Z₇ run's
row-completion profile was identical with and without it. The `--diag-zero`
seed had entered each orbit's *own column index* into its seed signature,
so no two orbits were ever classed as interchangeable. Solution counts (zero
either way) could not have shown this; the per-row completion counts did.
With the signature fixed the search finished.

### 7. The runs at 99

All on the (1, 7¹⁴) profile with the §1 seed, symmetry breaking on unless
stated, single-threaded Python:

| prunes | solutions | nodes | rows reached | completions per row 3,4,5,6,7,8,9,10 | time |
|---|---|---|---|---|---|
| none beyond (S)(R)(Q)+parity (run D) | **0** | 1 494 160 | 11 | 33, 414, 1135, 551, 57, 12, 7, 2 | 26.5 s |
| trace ∈ {0,42} (run C) | 0 | 1 493 998 | 9 | 33, 414, 1135, 551, 57, 6 | 22.4 s |
| all r_ii = 0 (run B) | 0 | 561 596 | 9 | 20, 194, 124, 76, 7, 1 | 13.7 s |
| r_ii = 0, r_ij ≤ 4, spectral a=6 | 0 | 492 085 | 9 | 20, 194, 124, 76, 4, 1 | 14.0 s |
| spectral a=12, trace 42 (run A) | 0 | 370 188 | 3 | — (33 of 37 row-3 completions cut by rank) | 17.5 s |
| trace ∈ {0,42}, symmetry from rows 0–2 only (run E) | see §8 | | | | |

Run D is the one that matters: **with nothing assumed beyond the orbit-
matrix conditions, the parity rule and the forced N(v₀) structure, the
order-7 profile admits no orbit matrix.** Rounds 2–3 were needed to *find*
the search that terminates, but the final statement does not lean on them.

Reproduce:

```bash
python problems/conway-99/explore/orbit7.py --profile z7                      # run D
python problems/conway-99/explore/orbit7.py --profile z7 --trace 0,42         # run C
python problems/conway-99/explore/orbit7.py --profile z7 --diag-zero          # run B
python problems/conway-99/explore/orbit7.py --profile z7 --spectral 12 --trace 42   # run A
python problems/conway-99/explore/orbit7.py --profile z7 --dump-rows 3 --out x.json  # rows 1-2
python problems/conway-99/explore/char_blocks.py
MATHLAB_SLOW=1 python -m pytest tests/test_conway99_orbit7.py -q
```

### 8. Verification performed

- **Contains real orbit matrices.** Paley(9)/Z₃, Paley(9)/Z₂ (one fixed
  point), rook(4)/Z₃ (one fixed point, 5 orbits — the same shape as the 99
  case), Petersen/Z₃ (one fixed point): the real orbit matrix is in the
  enumeration with symmetry breaking on.
- **Agrees with 002's independent enumerator** (`orbit_general.py`, symmetry
  off) on 11 profiles, including [2⁸] for SRG(16,6,2,2) with 39 620
  solutions in 12 classes, and every Z₃₃ profile (1 / 0 / 6, as 002
  reported).
- **Symmetry breaking loses no class** on four profiles (canonical class
  sets equal with it on and off).
- **Cross-family re-implementation.** Of three Qwen enumerator drafts,
  `enum-3` (plain DFS with forward checking, unpatched) reproduces every
  calibration count including the mixed-size ones (2, 1, 6, 3, 80). It was
  then run on the Z₇ profile seeded with rows 0–2; without symmetry breaking
  that search is far larger (see run E), and its status at write-up is in
  §10.
- **Run E** — the same enumerator with the lex prune restricted to rows 0–2
  (i.e. exactly the trivially valid "sort rows 1 and 2" and nothing else) —
  was still running after an hour at write-up. So the rows-3..14 symmetry
  break has been validated only on small profiles, not on this instance.

### 9. The swarm

Throughput first, because it set every other decision: the 27B dense model
at Q4 fills the 24 GB card, decodes at ~6 tokens/s per slot with two slots,
and in thinking mode a 150-word question ran 738 s to the 8192-token cap
without answering. Everything below is **thinking off**. Verification
briefs took 210–550 s each; the first version of the code briefs, capped at
6000 tokens, was cut off writing prose and had to be re-issued as
"code only, under 250 lines".

| brief | return | director verdict |
|---|---|---|
| enum-1 (row-by-row, inner products) | 370 lines | crashes (`NameError`, closure scoping) — **discarded** |
| enum-2 (row multisets) | 208 lines | `continue` outside a loop ×3; patched to `return` (documented); then correct on equal-size profiles, **0 on every mixed-size profile** (expected 6, 3, 80) — **fails calibration** |
| enum-3 (plain DFS, forward checking) | 311 lines | **passes all five calibrations unpatched**; used as the cross-family check |
| lift-1 (Z_m connection-set lifter) | 517 lines | returns 0 graphs, 0 nodes, `exhausted: false` on its own Paley(9) positive control — **fails control** |
| skeptic L3-a | confirms, step by step | agrees with the director's pass (below) |
| skeptic L3-b | opens "critical gap … FLAWED", re-derives the step, retracts, confirms | agrees; the false alarm is recorded as the model's characteristic failure shape |
| skeptic L5-a/b/c, skeptic r₁₂, z33-a/b | pending at write-up | see §10 |
| ideation, 18 lenses | pending at write-up | see §10 |

Director's own skeptic pass on what the result rests on: **L3** — the fixed
set is independent (a fixed vertex's neighbourhood carries a 7-cycle of
matching edges, so no fixed neighbours), closed under common neighbours
(odd order, set of size ≤ 2), and two fixed non-adjacent vertices would
force two fixed common neighbours adjacent to them: contradiction with
independence; |Fix| ≡ 1 (mod 7) gives exactly one. Confirmed. The
**N(v₀) structure** (§1): confirmed. The **parity rule**: confirmed. The
one director error found on this side: the enumerator briefs told workers
that [3,3,3] with k=4, λ=1, μ=2 has "exactly 3" solutions; it has 2. A
draft that trusted the brief over its own output would have looked wrong
for the wrong reason; enum-3 reported 2.

## Outcome

`LIVE` — **no SRG(99,14,1,2) admits an automorphism of order 7.** The chain
is: L3 (001, proof; director skeptic pass here, two worker passes concur)
→ the forced N(v₀) structure (§1, proof) → parity of r_ii (proof) → run D
(exhaustive enumeration, zero orbit matrices). The computation has been
calibrated against real graphs of the same shape and against two
independent implementations on smaller profiles, and repeated under five
different prune sets with the same answer; **it has not yet been reproduced
on this instance by a different implementation**, which is the bar for
`VERIFIED`. What would upgrade it: enum-3 (or any independent enumerator)
completing the rows-0–2-seeded profile, or run E completing.

`VERIFIED` (ranges): the calibration facts of §8; the character-block
lemma's twelve identities on four real graphs; rows 1–2 forced (hand
derivation §2 and enumeration agree).

`LIVE` (proved here, no skeptic pass): the character-block lemma (§4) and
its consequences a ∈ {0,6,12}, tr R ∈ {0,42}; the projector argument (§5)
giving a = 6, all r_ii = 0, r_ij ≤ 4. None of these is load-bearing for
the headline — run D does without them — but they are what makes the
enumeration finish in seconds rather than never, and they transfer.

**What is not claimed.** Nothing about orders 2, 3, 6, 9, 11 or the trivial
group. Nothing about the literature: PROBLEM.md's `[T]` summary says
Cesarz–Woldar (2023) leave G ≅ Z₇ as the only possibility when 7 divides
|G|, and Behbahani–Lam (2011) worked with exactly these orbit matrices; if
order 7 is already excluded in print, this record is an independent
re-derivation and should be marked as such by whoever reads the paper.
The swarm's ideation returns are candidates and enter no claim.

## Why it failed / what survived

Nothing failed at the target; what follows names the obstruction to the
*next* cases and what generalises.

- **Where the same machine stops.** Orders 3, 9 and 11 have no fixed point
  (L2, L5) and so no seed: 33, 11 and 9 interchangeable orbits, all rows
  free. For Z₁₁ the lemma gives 54 = a + 10r, so a = 4, r = 5, tr R = 10 —
  a single value, the best case; for Z₉ the two Galois classes of
  characters (six of order 9, two of order 3) give a ≡ 0 (mod 2) and
  nothing more. Whether the 9-orbit case terminates is in §10.
- **The multiplicity lemma is the successor 002 asked for**, at least for
  the trace: 002 noted that at t = 3 the diagonal equations form a square
  system and the trick evaporates at t = 9. The lemma is a constraint of a
  different kind — on the *spectrum* of R rather than its entries — and it
  bites at any t. The projector PSD/rank prune is its entrywise form.
- **The methodological lesson** is the same as 001's and 002's, in a new
  shape: a broken symmetry break is invisible in solution counts (zero
  either way). The row-completion profile caught it. Record such profiles.
- **Reusable:** `orbit7.py` (general sizes, seeds, all prunes optional),
  `char_blocks.py` (exact Q(ζ_p) arithmetic, orbit matrices from real
  automorphisms), the four-graph calibration set, and the `local` swarm
  provider with its measured throughput.
- **What the local worker family is good for, measured:** one correct
  enumerator in three drafts, zero correct lifters in one, two of two
  correct skeptic passes (one with a self-retracted false alarm). At 6
  tokens/s and thinking off, it is a slow floor tier for programs under 300
  lines and for paragraph returns; it is not a place to send anything that
  needs multi-step reasoning inside the answer. A faster serving
  configuration (the 35B-A3B MoE at ~100 tokens/s on this card, per the
  local-llm-server notes) would change that arithmetic.

## Leads generated

1. **Reproduce run D independently.** Cheapest route: seed `enum-3` (or a
   new implementation) with rows 0–2 *and* with each of the 33 row-3
   completions from run D, one job each; every job must return zero. That
   only checks rows 4–14, so add: enumerate row 3's completions with an
   independent implementation and compare with the 33. Definite outcome;
   upgrades the headline to `VERIFIED`.
2. **Order 11 with the same machine**: 9 orbits of 11, tr R = 10 forced,
   spectral a = 4. `orbit7.py --sizes 11,…,11 --trace 10 [--spectral 4]`.
   Either it terminates (then lift or eliminate), or its row profile says
   precisely where the 9-orbit case is wider than the 15-orbit one.
3. **Order 3**: profiles (3³³) and (1³, 3³²) by L5. Try the seeded one
   first: three fixed vertices forming a triangle (L5's f = 3 case) give
   three known rows.
4. **Skeptic pass on the character-block lemma and the projector argument**
   as written in §4–5 — small, self-contained, and the transferable part.
5. **Generalise the lemma to composite order** (Z₉, Z₃₃): characters fall
   into Galois classes with a shared multiplicity per class; write the
   constraint set and check whether it pins tr R for Z₉.
6. **Read Behbahani–Lam (2011) and Cesarz–Woldar (arXiv:2308.02978)** and
   settle whether order 7 was open. Pure reading; changes the label on this
   record either way.

## References

- `attempts/001-forced-pair-model-and-partner-regular-split.md` — L1–L4,
  the pair model. L3 is load-bearing here.
- `attempts/002-assumed-symmetry-construction.md` — orbit-matrix
  conditions, the Z₃₃ elimination, `orbit_general.py`, the leads followed.
- `problems/conway-99/PROBLEM.md` — tier-0 statement `[T]`; its summary of
  Makhnev–Minakova (2004), Behbahani–Lam (2011), Cesarz–Woldar
  (arXiv:2308.02978), Keramatipour (arXiv:2604.23037). None read.
- `docs/SWARM.md` — the director/worker protocol; the `local` provider
  added in this cycle.
- Data: `out/conway99-003/` (gitignored run logs and swarm drafts; the
  briefs are committed under `explore/swarm-003/`), swarm metas under
  `out/swarm/conway99-003-*/`.

## Addendum — results that landed after the write-up above

### 10. Late results

**Run E (symmetry break from rows 0–2 only, trace ∈ {0,42}):** did not
finish — 400 M nodes in 1689 s, no solution, row completions 712, 99 715,
595 197, 337 695, 35 739, 4041 for rows 3–8. Compared with run D's 33, 414,
1135, 551, 57, 12 this measures what the rows-3..14 lex prune buys (a
factor of 20 at row 3 rising to ~500 by row 6). Inconclusive as a check of
that prune on this instance; the validation of it remains the four small
profiles in §8.

**enum-3 on the rows-0–2-seeded profile:** 30-minute cap reached without
finishing or reporting. Inconclusive.

**Positive control at the same shape (`bvls_control.py`).** BvLS(243) is
the coset graph of the cyclic ternary Golay code, so the cyclic coordinate
shift is an automorphism of order 11 = k/2 with exactly one fixed vertex,
22 orbits of 11, and N(v₀) = two orbits joined by the matching — the exact
analogue of the order-7 situation. Its real orbit matrix satisfies
(S)(R)(Q) and **matches the forced seed entry for entry** (r₀₁ = r₀₂ = 11,
r₁₂ = 1, r₁₁ = r₂₂ = 0), which checks the §1 derivation on a graph that
exists. Its trace is 20, which is what the multiplicity lemma predicts
(132 = a + 10r forces a = 12 and tr R = 9a − 88 = 20) — the lemma checked
at p = 11 with a fixed point. Whether `orbit7.py` *finds* this matrix under
the same prunes as at 99 is reported in the line below, written after the
run:

- BVLS CONTROL RESULT. The *full* enumeration of the 23-orbit profile did
  not finish in this session: a 30 M-node probe reached row 11 of 23 with
  hundreds of completions per row (this profile is far wider than the
  15-orbit one — row 1 alone has 20 free columns), so the plain control is
  a multi-hour run and is left as lead 8. Two cheaper controls of the same
  machinery were run instead, both in `bvls_deep_seed.py`: seeding the
  enumerator with rows 0..11 of the real orbit matrix (symmetry off) it
  finds **exactly one completion, the real matrix, exhaustively** (12 104
  nodes); seeded with rows 0..7 only, likewise **exactly one completion,
  the real matrix, exhaustively** (373 265 nodes, fifteen rows searched).
  Seeded with rows 0..5 it did not finish in 280 s. So the row-by-row
  search — (S) transposition, inner-product constraints, sum of squares —
  loses nothing on a 23-orbit fixed-point profile from depth 8 onward; the
  symmetry break is validated only on the small profiles of §8.

**Order 11, same machine.** By 001's L2 an order-11 automorphism is
fixed-point-free, so its profile is nine orbits of 11 with no seed. The
lemma gives 54 = a + 10r, a ≤ 8, hence a = 4, r = 5, tr R = 10 — one
value. `orbit7.py --sizes 11×9` finishes **exhaustively in 0.2 s with zero
orbit matrices** (120 465 nodes; row completions 7, 41, 47, 11, 2), with
or without `--trace 10`, with or without `--spectral 4`. Without symmetry
breaking it does not finish in 50 M nodes. So, at the same status as the
order-7 statement and resting on L2 instead of L3: **no SRG(99,14,1,2)
admits an automorphism of order 11.** 002 reported that its enumerator
"does not reach a first solution in minutes" on this profile; the
difference is entirely the row-inner-product pruning and the symmetry
break. Combined with the published restrictions as PROBLEM.md transcribes
them `[T]` — |G| divides 2·3³·7·11, and |G| even ⟹ |G| divides 6 — these
two eliminations would leave |G| dividing 2·3³ with the even case being
Z₂, Z₆ or S₃. That combination is a remark, not a claim of this record.

**Order 9, same machine.** By 002's L5 an order-9 automorphism is
semiregular: eleven orbits of 9. The characters of Z₉ fall into two Galois
classes (six of order 9, two of order 3), so the lemma gives only
54 = a + 6r₁ + 2r₂, i.e. a even, and tr R = 7a − 26 ∈ {2, 16, 30, 44}.
`orbit7.py --sizes 9×11 --trace 2,16,30,44` finishes exhaustively in 355 s
(906 251 nodes; row completions 8, 62, 307, 369, 111, 30, 10, 4, 3, 1) with
**exactly one orbit matrix class**, saved in `out/conway99-003/z9.json`
and reproduced by the command in §7's style. This is the first surviving
orbit matrix at 99 outside Z₃₃, and it is small: lifting it is a search
over connection sets in Z₉ for 11 orbits — the next concrete lead.

**Swarm returns that landed late.**

| brief | return | director verdict |
|---|---|---|
| skeptic L5-a | REFUTED: claims the mixed-degree case is "based on a false model" (degree-8 vertices adjacent only to degree-8 vertices) | **wrong.** The counting identity Σ_{x∈N_F(u)} deg_F(x) = 2f − 2 (re-derived here: closure puts every common neighbour of two F-vertices in F, so the sum is deg_F(u)(1+λ−μ) + μ(f−1)) with degrees in {2,8} forces f ≤ 9 from a degree-2 vertex and f ≥ 9 from a degree-8 one, hence f = 9 and the bipartite (8,2) structure L5 states; the worker's contrary structure violates the identity. |
| skeptic L5-b | CONFIRMED, with a correct minimal dependency list | agrees with the director's pass |
| skeptic L5-c | REFUTED: "asserts f = 9 without derivation and fails to rule out f = 15, 21, …" | **wrong** for the same reason: f = 9 is derived in two lines from the identity. Verified the SRG(33,8,1,2) exclusion itself correctly. |
| skeptic r₁₂ | CONFIRMED, and tightens the argument (a matching on 7 vertices has at most 3 edges, so an orbit of 7 matching edges cannot lie inside O₁) | agrees; the tighter form is adopted in §1's reading |
| z33-a, z33-b | both programs, run here, print count = 1 and [[2,6,6],[6,2,6],[6,6,2]] | **cross-family re-derivation of 002's unique orbit matrix** (Qwen; 002's two were the same agent). The index for 002 should be read with this in mind. |

Tally for the local family on lemma-skeptic briefs: 4 confirmations, 2
false refutations, 0 correct refutations, on lemmas the director judges
correct. A false refutation costs director time to adjudicate; it is the
same failure shape as L3-b's retracted alarm, one step further.

**Ideation sweep (18 lenses at 3000 tokens; 8 were cut off at that cap
and re-issued at 4000; the 12 that landed before commit are triaged here,
the rest sit in `out/swarm/conway99-003-ideate/` unread).** Three lenses —
group rings, spectral, linear algebra — independently propose the
character decomposition / spectral pruning of R that §4–5 carried out;
the channel converging on the route that worked is weak evidence it points
somewhere real, and no evidence of anything else. Association schemes:
refine the orbit partition to a coherent configuration and enumerate
connection sets with intersection-number constraints — the natural
propagation for **lifting the Z₉ orbit matrix** (lead 7). Coding theory:
A² = 12I + A + 2(J − I − A) gives A² ≡ A (mod 2), so A is idempotent over
F₂ and its 2-rank is constrained; a p-rank obstruction to a Z₉-invariant
lift is a falsifiable first step nobody here has taken (lead 9). Fixed-
point theory (the only T3 return): for an involution a fixed vertex's seven
matching edges are each fixed pointwise, swapped, or moved, and the
(fixed, moved) neighbour split is constrained by regularity — a concrete
enumeration for the Z₂ case (lead 10). SDP/Lovász-theta: needs a solver
and is moot for Z₇ now. Finite geometry: correctly finds no GQ/partial-
geometry reading, then speculates about Hadamard matrices — no purchase.
Integer-programming cuts on r_ii ∈ {0,2,4,6}: subsumed (the diagonal is
zero). Local search on a Fourier-domain energy: dead for Z₇, possibly a
fallback for the Z₉ lift. Extremal triangle counts by orbit type: subsumed
by (Q). Structural (distance-2 orbits): restates §2, no purchase.
Net: two new leads (9, 10) and one propagation idea for lead 7 from
eighteen briefs, which is about the yield `docs/SWARM.md` predicts.

7. **Lift the unique Z₉ orbit matrix** (from §10): eleven orbits of 9,
   connection sets in Z₉ with the group-ring equations, coherent-
   configuration refinement as the propagation. Definite outcome: a
   99-vertex graph with an order-9 automorphism, or the elimination of
   order 9 (which with L5 would also settle Z₂₇ redundantly and leave the
   3-part of |Aut| at most 9).
8. **Finish the BvLS(243) positive control** (`bvls_control.py`, multi-hour
   on one core): the full 23-orbit enumeration with symmetry breaking must
   contain the real order-11 orbit matrix. This is the one control that
   exercises the lex prune at the shape of the 99 runs; until it is done the
   symmetry break's validation rests on the small profiles of §8.
9. **2-rank obstruction to the Z₉ lift** (from the coding-theory lens):
   A² ≡ A (mod 2) makes A idempotent over F₂; compute what the 2-rank of a
   Z₉-invariant lift of R = 3I + J must be (block-circulant structure over
   F₂) and whether the SRG parameters admit it. Small exact computation.
10. **Involutions by fixed-neighbour counting** (from the fixed-point lens):
   for a fixed vertex v of an involution, classify its seven matching edges
   as pointwise fixed / swapped / moved and derive the allowed (fixed,
   moved) neighbour splits from regularity and μ = 2; then run the same
   orbit-matrix machinery on each profile (1ᶠ, 2ᵐ). The published `[T]`
   claim that an even |G| divides 6 makes this the whole even story.
