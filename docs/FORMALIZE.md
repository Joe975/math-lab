# Formal certificates: the rung above the skeptic

Adopted 2026-08-04, after OpenAI's Astra announcement (2026-08-02): ten
results on long-open problems, each shipped with a model-written reasoning
walkthrough *and* a Lean 4 certificate with a sorry-count of zero. Whatever
the eventual peer-review verdict on the individual results, the shape of the
release is the lesson: what made the claims credible — seven months after the
same company's "GPT-5 solved ten Erdős problems" collapsed into literature
rediscovery — was a kernel-checked artifact whose verdict does not depend on
trusting the model that produced it.

That is the same role our skeptic pass plays, with a stronger guarantee for
one class of claim. This doc defines when and how a result here earns the
`FORMALIZED` status.

## Where it sits

| Rung | What it establishes |
|---|---|
| claim | nothing (an agent said so) |
| skeptic-confirmed (`VERIFIED`) | an independent adversary re-derived / re-implemented it and failed to kill it |
| `FORMALIZED` | a proof checker's kernel accepts a formal statement of it |

Formalization is **optional and expensive**. It does not replace the skeptic
pass — it comes *after* it, and only for claims worth the cost.

## What qualifies

- **Proof-shaped claims only.** Lemmas, identities, case analyses,
  monotonicity facts. Computational censuses do not gain from Lean: their
  verification standard is already an independent re-implementation in exact
  arithmetic (see the 2026-07-31 verification-standard note in `STATUS.md`),
  and formalizing "the search was exhaustive" is a different, much larger job
  than formalizing a lemma.
- **Skeptic-confirmed first.** Formalizing an unreviewed claim wastes the
  expensive step on statements that a cheap adversary pass might kill.
- **Load-bearing and reused.** Prefer steps that many later records lean on;
  a formalized wrong-turn is effort spent making a dead end rigorous.

## The bridge caveat — read this before celebrating

A Lean certificate proves exactly the **formal statement**, not the claim in
the attempt record. The informal→formal translation is the new risk surface,
and it is precisely the "hand-derived bridge" failure mode this lab already
recorded for formal-parameter specialization (STATUS insight, 2026-07-31):
the formal check cannot see a wrong bridge. Astra's coverage has the same
gap — a compiling proof of a mis-stated theorem is worth nothing.

So a `FORMALIZED` pass has three parts, and the first is the point:

1. **Statement review.** An independent agent (not the formalizer) compares
   the formal statement, hypothesis by hypothesis, against the claim as
   written in the attempt record, and says so explicitly in the review
   record. Any narrowing (a stronger hypothesis, a special case) is listed
   the way `range` is listed for a search.
2. **Clean build.** Zero `sorry`, no axioms beyond the standard library /
   mathlib, and the full source in the repo.
3. **Recorded toolchain.** Lean toolchain version, mathlib pin, and the
   `lake build` output captured in the record. CI does not build Lean (cost);
   the reproducibility contract is that the skeptic's own independent
   `lake build` — on a different machine or container — is part of the pass.

## Mechanics

- Files live in `problems/<problem>/formal/` (tier 1; a blind agent never
  sees them).
- Record the work as a **review-shape attempt** (`verifies: <id>` in the
  index) with `status: FORMALIZED`. The scope line states the formal theorem
  proved, which may be narrower than the record's claim — that difference is
  data, not a footnote.
- `FORMALIZED` is in the schema's attempt-status enum and the `AGENTS.md`
  vocabulary. Route-level status stays `VERIFIED`/`LIVE`/etc.; formalization
  upgrades a *step*, not a conjecture.
- Lean is an exception to the stdlib-only-Python rule (that rule governs the
  Python tooling). Nothing outside `formal/` may depend on it.

## Pilot candidates (2026-08-04)

Ordered by (value ÷ effort), from the current ledger:

1. **Billiards Laurent-identity block** (005/008): I1–I4 and the glide facts
   are already "one exact polynomial check" after adjoining e^{iaα}, e^{ibβ}
   as formal variables — the closest thing in the ledger to a proof that is
   already formal. Formalizing the polynomial identity is the easy half; the
   *bridge* (geometry = closed form) is exactly what statement review is for.
2. **The c·cot(cs) monotonicity lemma** (005/008): elementary, one-variable,
   mathlib has the pieces; it retires Lemmas C and D.
3. **Union-closed P6′** (008/012): quantitative IFT for the symmetric
   Sinkhorn potential. Hardest of the three; attempt only after 1–2 have
   established the workflow.

A pilot that fails is recorded like any other attempt: what resisted
formalization, and whether the obstruction was the mathematics or the
tooling, is the deliverable.
