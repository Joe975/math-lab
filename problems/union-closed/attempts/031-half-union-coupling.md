# 031 — The half-union coupling: a total, closed-form coupling that beats the λ-sweep on every G-gen hard instance of record

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** construction + adversarial sweep (queue 1(b), opened by 030's
  HU-mix theorem). Persistence window (round 2), same author as 025–028/030; 029's reviewer scope ends at 028 — this record joins the next reviewer batch.
- **Tools:** `explore/uc_hu_probe.py` (new; the HU builder, floor/cap/
  slice/crash sweeps; stdlib; deterministic — the only RNG is
  fixed-seeded random-support generation for the scans, seeds recorded
  in-file; checkpoint `data/hu_probe.json`);
  `explore/uc_hu_skeptic.py` (new; independent recursive-descent
  builder, independent state-major slice DP, natural-log evaluator via
  the 025 skeptic's `cr_chain`; exit 0). Reproduce: run both.
- **Sources:** Sawin arXiv:2211.11504 (the max-entropy coupling
  ingredient; as transcribed and verified in 003/004 — no new fetch of
  the PDF; the framing below also checked against the public abstract
  listing this cycle).

## Approach

### The coupling

For any μ and the fixed revelation order, define the **half-union (HU)
coupling** per history cell: with x = μ(A_i = 0 | a), y = μ(B_i = 0 | b),
set the conditional both-zero probability

    z(a, b) = clip(1/2; [max(0, x + y − 1), min(x, y)])
             = min( max(1/2, x + y − 1), x, y ),

i.e. push each coordinate's conditional union marginal to 1/2 whenever
the cell's Fréchet bounds allow, else to the nearest feasible point.
Because the cell margins are exactly μ's conditionals for every partner
history, both processes are μ-distributed: the construction is a genuine
coupling, **total** (defined for every μ — no genre detection), and
**closed-form** (no Sinkhorn iteration, no λ parameter).

**Provenance, per the novelty gate:** the per-coordinate move "maximize
union-marginal entropy" is Sawin's max-entropy-coupling ingredient
(2211.11504) and is already in this repo as 009 part D's adaptive
m_k = 1/2 rule. What is new here is the clamped per-history totalization
to arbitrary μ, and the measurements/theorem below. This record claims a
construction-level synthesis, not a new idea at the per-coordinate level.

On the [1/4, 1/2]-component mixture sub-genre the clamp never fires and
030's HU-mix theorem applies verbatim: CR = n − H(μ) ≥
(1−h(p̄))/h(p̄)·H(μ), tight exactly at products Bern(p̄)^⊗n (where
CR = n(1−h(p̄)) and H = n·h(p̄) — verified as an identity to 1e-9 by the
skeptic at p ∈ {0.3, 0.38271, 0.45}).

## What was done

Every number below is reproduced by the independent builder/evaluator
(`uc_hu_skeptic.py`, exit 0).

**1. The recorded hard instances — HU dominates everything.**

    instance          n   best prior branch     CR_HU     CR_HU/H
    windowkill floor  4   sinkhorn +0.2232     +0.3090    0.1155
    mmabskill floor   5   sinkhorn +0.0742     +0.0822    0.0429
    mmabskill_n6      6   sinkhorn +0.1013     +0.3007    0.1289
    random0/1 floors  4/5 sinkhorn +0.126/.197 +0.156/.299 0.065/.106
    crash family      5–17 λ-sweep (window→0)  +0.19919   0.1130 (flat in n)
    slices w0/n=.42   100 tilt (Θ(n) gain)     +3.442     0.0364

  The committed sharpest floor rises 0.0742 → 0.0822 (ratio 0.0387 →
  0.0429). The crash family — which forces every λ-based branch into the
  λ → 0 corner (009's window law) — reads as a flat, n-independent
  +0.199 because HU has no λ to crash; it also beats the iid point
  (+0.111, 009). Marginal exactness of the build: ≤ 1e-16 everywhere.

**2. Adversarial scans, in-regime (cap 0.38271):** 300 fixed-seed random
supports (n ≤ 6, ≤ 14 atoms, δ∅-projected into regime): zero violations,
worst CR/H = +0.0496. Deterministic pattern-search descent on CR/H from
the floor endpoints stalls at +0.0429/+0.0430 — right at the
(1−h(marg))/h(marg) scale — and cannot push below.

**3. Above ψ.** Same scans + descents at caps 0.40/0.43/0.46/0.49: zero
violations; worst attacked CR/H = 0.0331/0.0166/0.0059/0.00065, each
sitting just above the product-extremal value (1−h(cap))/h(cap) =
0.0299/0.0144/0.0046/0.00029. Slices by DP (exchangeable states
(w_a, w_b), two independent implementations agreeing to 2e-3 at
n = 300): CR_HU > 0 at every w0/n ∈ {0.40, 0.42, 0.45, 0.48} up to
n = 300 (0.40-series: CR/H → 0.0386 at n = 300, per-coordinate slope
→ ≈ 0.032 ≥ the product slope 1 − h(0.40) = 0.029).

**4. The conjecture this data writes down.**

  **(HU-TAX at p̄), CONJECTURE — labeled SPECULATION at every use:**
  for every μ with all element marginals ≤ p̄ < 1/2,
  CR_HU(μ) ≥ ((1−h(p̄))/h(p̄)) · H(μ), with equality exactly at
  Bern(p̄)^⊗n. Status: proved on the [1/4,1/2]-component mixture
  sub-genre (030); tight at products by identity; consistent with every
  attack above (the attacked floor tracks the extremal constant from
  above at all five caps); wide open in general.

## Outcome

**EVIDENCE**, scope stated precisely: the HU coupling has CR > 0 on
every instance attacked — the full 023/024/025/027/028 hard-instance
set, 300 + 4×200 fixed-seed random supports at caps up to 0.49 (n ≤ 6),
descent floors from every start, slices to n = 300 at fractions to 0.48,
and the crash family to n = 17 — and its measured floor ratio tracks the
product-extremal constant (1−h(p̄))/h(p̄) at every cap tested. Nothing
here is a bound for unattacked μ; the conjecture is SPECULATION.

**Not claimed, emphatically:** any proof of (HU-TAX) or even of
CR_HU ≥ 0 beyond the sub-genre theorem; any statement for n beyond the
ranges above; progress on Frankl itself. Calibration note: empirical
positivity of one coupling past ψ is *consistent with the known
landscape* (coupling gains above ψ are exactly why the dependent-
couplings route exists — Sawin's improvement lives there, and 003/004
verified Θ(n) slice gains at p ≈ 0.42 to n = 2·10⁶); the conjecture, if
proved for all p̄ < 1/2, would be Frankl via licensing — which is
precisely why it will not fall easily and why every use of it must stay
labeled. The route's difficulty has not shrunk; it has *relocated* into
one concrete inequality about one concrete coupling.

**5. Scope checks run before any dominance claim.** (i) Head-to-head
against 030's optimized K-mixing on the 030 battery (n = 4): K-mixing
wins 3 of 5 (lo+mid+hi +0.833 vs HU +0.595; spread +0.582 vs +0.546;
near0+2hi +0.403 vs +0.269); HU wins all-low+1mid and ties near-ψ. So HU
does NOT subsume the mixture branch — on structured mixtures the label
information buys real CR that the history-only HU rule cannot see. The
dominance claim is strictly: over the λ-swept Sinkhorn branch on the
G-gen hard instances of record, and over every branch on the crash
family and slices tested. (ii) HU depends on the coordinate revelation
order: over all orders at the floor instances, CR ranges e.g.
[+0.0822, +0.3382] at mmabskill (identity order is the min there) and
[+0.1180, +0.3399] at random0 — **min over orders stayed positive at
every instance checked** (full enumeration n ≤ 5, 120-order sample at
n = 6). Order is a free lever the λ-sweep never had; (HU-TAX) below is
posed for every order, and all sweeps in this record used the identity
order (the conservative reading at the instances where that is the
minimum, e.g. mmabskill).

## Why it failed / what survived

The attack failed to find any crack, which is itself the finding: the
G-gen branch (λ-swept Sinkhorn) is dominated on every recorded G-gen
hard instance by a coupling with no free parameter at all, while the
structured-mixture side keeps 030's K-mixing as its stronger branch. Two structural lemmas fell out of the proof attempt (both proved here,
pending reviewer pass):

  **Lemma HU-notax.** For the HU coupling, T_A = T_B = 0 identically:
  each cell joint's A-margin is exactly μ(A_i = 0 | a) for every partner
  history b, so x̃ᵢ is A_{<i}-measurable and Σᵢ E[h(x̃ᵢ)] =
  Σᵢ H(A_i | A_{<i}) = H(μ). (The λ-tilts never had this; HU's only
  loss channel is the second tax.) Hence, exactly:

      CR_HU = Σᵢ E[ h(z_i) − (h(x_i) + h(y_i))/2 ] ,

  a pure per-cell ledger. **Per-cell sign analysis:** in cells with both
  x, y ≥ 1/2 the term is ≥ 0 always (z = 1/2 gives h = 1 ≥ avg; z =
  x + y − 1 ∈ [1/2, min(x, y)] gives h(z) ≥ max(h(x), h(y))); the only
  deficit cells are mixed (x < 1/2 ≤ y) and both-below cells, where
  z = min(x, y) and h(min) can undershoot the average. In-regime the
  unconditional zero-probabilities are ≥ 1 − p̄ ≥ 0.617, so deficit
  cells exist only where history drags a conditional below 1/2 — mass-
  limited by exactly the entropy structure the bound is about. Where the
  proof difficulty now sits: bound Σ(deficit-cell mass × deficit) by
  Σ(surplus-cell mass × surplus) − c·H(μ). The same averaging
  difficulty as L6, but with a fixed coupling, zero first tax, a
  conjectured sharp constant, and known extremals — a materially
  better-posed target than 023's open-ended "CR ≥ c(p)·H on the generic
  genre".

Recipe consequence (v2 → v3 candidate): replace the G-gen branch's
λ-swept Sinkhorn with HU (strictly better on every recorded G-gen hard
instance, total, and parameter-free); keep 030's K-mixing on declared
mixtures (it beats HU there, §5). Two branches instead of recipe v2's
two — but the generic one is now closed-form and total, and both floors
rose.

## Leads generated

1. **Adversarial campaign against HU as the primary target** (the 023
   pattern: anneal with H ≥ ε guard, stable seeds, PYTHONHASHSEED
   recorded — plus the coordinate-order degree of freedom: HU depends on
   the revelation order; scan orders, since an order-dependent kill
   would refute totality-as-stated).
2. **Prove CR_HU ≥ 0 in-regime** before the H-scaled version: the
   clamp structure gives Σ E h(z_i) ≥ Σ E h(max(x̃, 1/2-feasible…))
   — a monotone-rearrangement argument may close ≥ 0 even if the sharp
   constant stays open.
3. **Order-optimized HU**: a canonical order rule (greedy per-step
   entropy?) and a joint (instance, order) adversarial descent — §5's
   4× spread at mmabskill says the identity order leaves real CR on the
   table, and an order-adversarial kill would sharpen (HU-TAX)'s
   quantifier.
4. **Exact certification** of CR_HU at the floor instances (the cell
   tree is finite and the z-values are rational functions of μ's
   conditionals — the 026 dual-kit machinery applies).
5. If (HU-TAX)'s averaging inequality resists: swarm brief for
   cross-family drafts of averaging schemes (the 020 pattern), since
   the statement is now compact enough to brief.

## References

- This repo: 030 (HU-mix theorem, near-ψ pinning), 009 (identity (ii),
  adaptive rule, crash family), 023/024 (floor instances), 025–028 (the
  branch landscape this replaces), 004/003 (Sawin/Liu transcriptions),
  `data/hu_probe.json`.
- Sawin arXiv:2211.11504 (max-entropy coupling ingredient; per §The
  coupling, used at idea-provenance level only).
