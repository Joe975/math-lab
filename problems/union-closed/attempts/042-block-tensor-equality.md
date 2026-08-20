# 042 — The block-tensor equality family holds and resists attack at its own boundary; crash8 and the 0.497 saturation close out

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** equality-family verification + adversarial attack +
  close-outs (041 lead 1; 039 correction lead; 040 lead 3).
- **Tools:** `explore/uc_hu_blocks.py` (new; block-tensor builder,
  equality verification under arbitrary orders, anneal+polish attacks
  seeded at equality points, the two close-outs; deterministic, seeds
  951000+; checkpoint `data/hu_blocks.json`);
  `explore/uc_hu_blocks_skeptic.py` (new; own tensor builder, every
  value through the independent 037-skeptic stack; exit 0). Reproduce:
  run the two in that order.
- **Sources:** none.

## Approach

041 noted that CR_HU is additive over independent blocks. The short
proof, recorded here since the equality family leans on it: if
μ = μ₁ ⊗ μ₂ on disjoint coordinate blocks, then for any revelation
order each coordinate's cell conditionals x, y depend only on the
revealed prefix WITHIN its own block (independence), so each step's
transition matrix is the one the block coupling alone would use, and
the joint cell process is the independent product of the two block
processes. Hence Σ w·h(z) adds per step over blocks, H adds, and
CR_HU(μ₁ ⊗ μ₂) = CR_HU(μ₁) + CR_HU(μ₂) for every order. Since
diag_m(p) contributes h(max(1/2, 1−2p)) − h(p) regardless of m (040's
DIAG identity) and Bern(p) is the m = 1 case, **every tensor of
diagonal blocks and Bernoulli factors at a common p has
CR/H = c*(p) exactly**: the equality set of (HU-TAX) is at least this
whole combinatorial family — one member per partition of [n] into
blocks — not two isolated points.

That makes the family's neighborhood the sharpest falsification
surface known, so this attempt (i) verifies the equality claim
numerically under arbitrary orders including the p < 1/4 clamp branch,
(ii) attacks rollout-HU with 040's anneal seeded AT equality points,
and (iii) clears the two cheap queued close-outs.

## What was done

**A. Equality verification.** Eight block-tensor cases (partitions of
n = 4/6 into diagonal blocks and Bernoulli factors; p from 0.20 —
the clamp branch of c* — to 0.49), each under identity, rollout, and
three random interleaved orders: max |CR/H − c*(p)| = **4.4e-16** over
all 40 case-order pairs. The additivity proof above is exercised
across non-contiguous interleavings, which is where it could have
failed.

**B. Attack at the equality surface.** Anneal (4000 steps) + descent
polish against rollout, seeded at block-tensors at p = cap − 0.003:

    cap 0.49  n=4 d2⊗d2       → +0.000488   cap 0.497 n=4 d2⊗d2       → +0.000104
    cap 0.49  n=6 d3⊗d3       → +0.000488   cap 0.497 n=6 d2⊗d2⊗d2    → +0.000061
    cap 0.49  n=6 d2⊗d2⊗d2    → +0.000364          [extremals +0.000289 / +0.000026]

  Zero violations. The attack does exactly one thing: it slides ALONG
  the equality family toward the cap (the polished endpoints are
  single diagonals at p closer to the cap — e.g. +0.000364 ≈ c\*(0.4887)),
  never below it. Seeded at the family's own points, the adversary
  finds no direction off the family that loses value.

**C. Close-outs.** (i) **crash8** — the one standing start never
descended under rollout (039's correction): floors +0.000952 (cap
0.49) and +0.000310 (0.497), positive, grinding toward the extremal
like everything else. (ii) **0.497 saturation** (040 lead 3): the
diagonal-seeded anneal at cap 0.497 ends at +0.00010389 =
c\*(0.494) to 1e-12 — that cap's floor is now family-saturated, same
as 0.49 and 0.499 were in 040.

## Outcome

- **VERIFIED (elementary proof + machine-precision check): the
  equality set of (HU-TAX) contains the full block-tensor family**
  {⊗ⱼ diag_{mⱼ}(p) : partitions of [n]} at every p < 1/2, every
  order. (Additivity proof same-session; joins the reviewer batch.)
- **EVIDENCE: rollout-HU survives attack seeded at the equality
  surface itself** — five anneal+polish campaigns at caps 0.49/0.497,
  zero violations, every endpoint ON the family (single diagonals
  nearer the cap).
- **CLOSED: 039's crash8 gap** (positive under rollout descent at
  both caps) **and 040 lead 3** (0.497 family-saturated). Every
  standing start of record has now been descended under rollout, and
  every attacked cap's floor sits on the diagonal family.
- **Not claimed:** that the block-tensor family is the WHOLE equality
  set (equality at p < 1/2 requires every cell at the clamp boundary;
  characterizing that set exactly is lead 1); anything n > 8.

## Why it failed / what survived

The falsification surface held. The picture after 040/041/042 is
coherent: the extremal landscape of (HU-TAX) at cap p̄ is the
block-tensor family, every adversary of record converges onto it from
every seeding tried (random, hostile, embedded kills, the family
itself), and the margin grows in every off-family direction measured.
For the proof effort (040 lead 2), the equality characterization now
has a concrete conjecture shape: CR/H = c*(p) iff μ is a block-tensor
of diagonals at marginal p — "iff every history cell sits at its clamp
boundary and the posterior collapses after one live reveal per block".

## Leads generated

1. **Characterize the equality set**: prove CR_HU/H ≥ c*(p) with
   equality iff block-tensor (the "iff" is the new content; the "≥"
   is (HU-TAX) itself). The cell-boundary framing above is the
   suggested route; a Lean-able lemma if it lands.
2. The reviewer batch for 041/042 should re-derive: the count-pair DP
   (041), the additivity proof and its interleaving claim (042 A),
   and re-run the attack endpoints (042 B).
3. Exhausted here: no more cheap close-outs stand — the next
   substantive moves are proof-shaped (this lead 1, 038's sandwich,
   031's b2) or scale-shaped (n > 8 needs either the DP genre or new
   engine work).

## References

- This repo: 041 (the additivity observation and the lead), 040 (DIAG
  identity, anneal machinery), 039 (the crash8 correction), 037/038
  (rollout), 034 (c*), 031/030 (HU). `data/hu_blocks.json`.
- No external sources.
