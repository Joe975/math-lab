# Entry point for agents

You are here to attack an open conjecture. Pick a problem, pick a **mode**,
work, then submit an attempt record. This file tells you exactly what to read.

## Pick a mode first

The repo is cut into two tiers. Which one you read determines what your result
means, so choose deliberately.

### `blind` — attack the problem independently

```
scripts/blind.sh <problem> ../work-<problem>
```

That materializes a working copy with the problem statement, the shared
harness, and the rules — and physically none of this lab's prior art. Work
there. You cannot re-tread our dead ends because you cannot see them, and you
cannot be anchored by our framing.

Choose this when you want your result to be **independent evidence**, or when
you suspect the existing routes are a local optimum.

### `informed` — build on what is here

Work in a normal clone. Read, for your problem:

1. `problems/<problem>/PROBLEM.md` — the statement and published background.
2. `problems/<problem>/prior-art.json` — **read this before the prose.** One
   line per attempt, with `mechanism` tags, `status`, and `gaps`. Use it to
   check whether your idea's mechanism family has already been tried.
3. `problems/<problem>/PRIOR-ART.md` — the narrative version, if the index
   says something relevant.
4. Only then the specific `attempts/NNN-*.md` records you actually need. They
   are long; do not read them all.
5. `STATUS.md` — live state across all problems, and the queue.

Choose this when you want to push an existing route further, or to attack a
labeled gap.

The informed side also has a cross-pollination layer: `mechanisms.json` at
the root maps every approach-family tag to a mathematical field, and
`scripts/mechanisms.py` queries where each family has been tried and which
fields are untried on a problem. `docs/IDEATE.md` (sweep field lenses across
one problem) and `docs/RIPPLE.md` (scan one new result across all problems)
are the procedures built on it. All tier 1 — none of it exists in a blind
copy, by design.

**Do not read tier-1 files and then report `mode: blind`.** The mode field is
data (see below); a wrong one corrupts it.

**Path note.** Attempt records written before the July 2026 restructure refer
to the old layout. `tools/<x>.py` is now either `harness/<problem>/<x>.py` or
`problems/<problem>/explore/<x>.py`, and `attempts/<problem>/data/` is now
`problems/<problem>/data/`. The records are kept as written rather than
rewritten, so translate as you read.

## The tiers

| Tier | Files | Blind agents see it |
|---|---|---|
| 0 | `problems/*/PROBLEM.md`, `harness/**`, `README.md`, `AGENTS.md`, `CONTRIBUTING.md` | yes |
| 1 | `problems/*/PRIOR-ART.md`, `prior-art.json`, `attempts/**`, `explore/**`, `data/**`, `STATUS.md`, `GUIDANCE.md` | no |

The cut is machine-readable in `tiers.json` and enforced by
`tests/test_tier_isolation.py`, which fails the build if a tier-1 finding
appears in tier-0 prose or in a blind checkout.

**Tier 0 is published background only.** Anything in it that isn't in the
literature is a bug — please report it.

## Why the mode matters beyond your own attempt

Every attempt records `mode`. Over time that turns the library into a dataset
on a question nobody has clean data for: **does prior knowledge help or
anchor?** Do blind agents rediscover known dead ends (prior art is
load-bearing), or do they find things informed agents were steered away from
(prior art is a local optimum)? Recording your mode honestly is what makes
that measurable, which is why a mislabelled mode is worse than no attempt.

## Ground rules, in short

Full version in `CONTRIBUTING.md`. The non-negotiable parts:

- **Nothing is a "result" until it survives adversarial verification.** Claimed
  proofs get an independent skeptic pass; computational claims get re-run.
- **Dead ends are first-class output.** An attempt that failed and says
  precisely *why* is a successful contribution. Most attempts here are those.
- **Computational evidence is `EVIDENCE`, never proof.** Say the range checked.
- **Label speculation `SPECULATION`.** Every unproven step, explicitly.
- Famous conjectures are long shots. The deliverables are the approach
  library, the tooling, and progress on the smaller open problems.

## Status vocabulary

Use these exact words; they are what the indexes are filtered on.

| Term | Means |
|---|---|
| `VERIFIED` | Independently re-derived or re-computed. State the range or scope. |
| `EVIDENCE` | Computational support, not proof. State the range. |
| `LIVE` | A route that has survived every adversary tried, with no bound claimed yet. |
| `SPECULATION` | Plausible, unproven, and load-bearing. Must be labelled inline. |
| `REFUTED` | Killed, with the obstruction recorded. Do not re-attempt as stated. |
| `MAP` | A survey/barrier analysis; no claim of progress. |
