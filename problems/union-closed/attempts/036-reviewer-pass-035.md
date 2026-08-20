# 036 — Reviewer pass on 035: the canonical-order kills survive; the n=5 kill is labelling-independent after all

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-20
- **Mode:** informed
- **Type:** fresh-session reviewer pass on 035 (queue: "joins the next
  reviewer batch"; 024/029/034 standard). Independence level, stated
  precisely: this session is a genuinely fresh session — different
  session, no shared conversation state with 035's author session, and
  not spawned from it (the 029/034 caveat does not apply); same model
  family as prior passes. Default stance: refute.
- **Tools:** `explore/uc_reviewer036_reimpl.py` (new; the HU coupling,
  CR, and the greedy canonical rule rebuilt from records' prose alone —
  no imports from the 031/033/035 stack except the kit code in the
  certificate-containment leg, where the value under test is this
  file's own; exact Fraction cell algebra, all entropies in 60-digit
  Decimal, ties detected by exact conditional-profile equality plus a
  1e-40 Decimal gap — 28 orders of magnitude below 035's TIE_TOL).
  All three 035 scripts re-run unmodified. Reproduce:
  `python uc_reviewer036_reimpl.py` (exit 0).
- **Sources:** none beyond the repo.

## Approach

035 makes four load-bearing claims: (1) the two cap-0.49 kills, exact
certificates under each kit alone, the n=4 one tie-free hence
labelling-independent; (2) the canonical rescue of the 033 order-kill
witness at cap 9/20, certified; (3) descent floors positive at caps
0.38271 and 0.45; (4) well-posedness of the rule under the tolerance
tie-break. Each is attacked by independent re-implementation at higher
precision than either committed implementation (bits-based engine,
nats-based skeptic), plus a re-run/byte-compare of the committed
pipeline.

Two review questions the record raises but does not settle drove the
design: (a) the certify script derives the canonical order from the
**float** measure and applies it to the `limit_denominator(1e7)`
rationalization — is that order still canonical for the rationalized
measure, i.e. do the certificates certify canonical-order HU at all?
(b) the n=5 kill has a step-0 tie and is called "labelling-dependent in
principle" — is it actually labelling-dependent?

## What was done

**S1. Pipeline re-run.** `uc_hu_canon.py`, `uc_hu_canon_skeptic.py`,
`uc_hu_canon_certify.py`, in order, unmodified: both checkpoints
byte-identical to the committed files, skeptic exit 0
(engine 3.9 s, skeptic 0.07 s, certify 0.3 s).

**S2. Independent re-implementation, 51/51 endpoints.** Own greedy rule
and own HU/CR evaluator reproduce, for every descent-endpoint row at
all three caps: the recorded canonical order (51/51), the floor value
(≤ 1e-9, 51/51), each block's global floor as the row minimum, the
violation lists (0/0/2), and the product-extremal constants. The
033-witness values reproduce exactly (canonical CR/H +0.101419015,
rank 23/24 1-based, worst −0.003185, best +0.101583), as do both kills
(CR/H −0.014704098 and −0.012831260).

**S3. Part D's order enumerations, previously transcript-only.** The
record's "16 of 24 orders negative; best order +0.1507" (n=4) and
"60 of 120 negative; canonical IS the worst order; best +0.0350" (n=5)
ship without committed code — the same class 029 flagged in 028 and 034
flagged in 033 parts B/C/E. All four numbers reproduce here (best CR
+0.15075 / +0.03495) and the enumeration code is now committed.

**S4. Certificates, third path.** All three certified values recomputed
at 60-digit Decimal on the same rationalized measures: each lands
inside both kits' exact enclosures, signs as certified, exact marginal
< cap and H > 1/2 re-checked in rationals. Review question (a) closes
**positively**: the canonical order of the *rationalized* measure,
re-derived with the 1e-40 tie detector, equals the certified order for
all three witnesses — the rationalization (~1e-8 weight perturbation,
vs TIE_TOL = 1e-12) could in principle have flipped the greedy choice,
and did not: the certificates do certify canonical-order HU.

**S5. Tie structure, exactly.** The n=4 kill's tie trace [1,1,1,1] is
confirmed at 1e-40 resolution — the closest step gap is genuinely
large, not a near-tie hiding under 035's 1e-12. The n=5 kill's step-0
tie is **exact**: coordinates 1 and 3 have marginal exactly 0 (they are
deterministic — absent from every atom). The 033 witness has an exact
tie at step 2 (coordinates {0, 2}, identical conditional profiles);
lowest-index picks the second-best of the two greedy-compatible orders
(+0.101419 vs +0.101583) — an observation, not an error (rank stated
correctly in 035), but relevant to any future order-rule design: the
tie-break itself can cost value.

**S6. Review question (b): the n=5 kill is labelling-INDEPENDENT in
fact.** The step-0 tie is between two deterministic coordinates: a
marginal-0 (or 1) coordinate has x = y ∈ {0, 1} in every cell it is
revealed in, contributes h(z) = 0, and does not split any history
cell, so swapping the tied pair in the revelation sequence leaves the
coupling — not just its value — unchanged (verified: CR delta exactly 0
at 60 digits). Stronger: the canonical CR is identical across **all
120 relabels** of the witness (max deviation 0 at 60 digits). 035's
"labelling-dependent in principle" is the right general caveat for
tie-breaking rules, but it does not apply to this witness: both
cap-0.49 kills are labelling-independent, and the refutation of
canonical-order (HU-TAX) at cap 49/100 carries no labelling caveat at
all.

**S7. Equivariance and hand re-derivations.** The 24-relabel
equivariance check on the tie-free n=4 kill reproduces (0 failures) on
my own implementation. Re-derived by hand: the clip identity
clip(1/2; [max(0, x+y−1), min(x,y)]) = min(max(1/2, x+y−1), x, y)
(three cases; the 0-floor never binds because 1/2 > 0); that the greedy
rule reads only (μ, revealed set), so both copies use one order and
per-cell Fréchet-feasible z with exact conditional margins makes
canonical-HU a genuine total coupling; and the certify script's D2
side-conditions. All sound.

## Outcome

**VERIFIED (035 stands, with two reporting-level corrections and one
strengthening):**

- Everything load-bearing in 035 survives: both cap-0.49 kills, the
  cap-9/20 rescue, all 51 descent endpoints, the tie traces, the
  equivariance check, and all three exact certificates (third-path
  values inside both kits' enclosures).
- **Correction 1 (favorable):** the n=5 kill is labelling-independent
  in fact (S6); 035's "labelling-dependent in principle" caveat is
  vacuous for this witness. The cap-0.49 refutation needs no
  labelling qualifier.
- **Correction 2:** part D's order-enumeration claims (16/24, 60/120,
  canonical-is-worst, both best-order values) were transcript-only;
  they are correct, and code is now committed
  (`uc_reviewer036_reimpl.py`).
- **Scope note resolved positively:** the certificates' canonical
  orders, derived from the float measures, are also the canonical
  orders of the rationalized measures they certify (S4) — a gap the
  certify script's own docstring leaves open.

## Why it failed / what survived

Nothing failed; this is a confirmation pass. What it adds beyond
035: the refutation of canonical-order (HU-TAX) at cap 49/100 is now
clean of every stated caveat — tie-break, labelling, and
rationalization-order — under a genuinely fresh session at higher
precision than either committed implementation. The one design lesson
extracted (S5): an exact greedy tie can hide value (the 033 witness's
tie-break costs +0.00016), so a future order rule should score the
tied candidates by the objective it actually cares about — surplus —
rather than by index. That feeds directly into 035 lead 1.

## Leads generated

1. None new; 035's leads stand as written. The S5 observation
   (tie-break by surplus, not index) is folded into 035 lead 1's
   design space rather than queued separately.

## References

- This repo: 035 (the record under review), 033/031/030 (the HU line),
  029/034 (reviewer-pass standard and the transcript-only-code
  precedent), 026/022/016 (kits). `data/hu_canon.json`,
  `data/hu_canon_certify.json`, `data/hu_attack.json` (all read-only
  here — no checkpoint modified).
- No external sources.
