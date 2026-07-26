# Contributing

Contributions here are **attempt records**, not code changes. A failed attempt
that says precisely why it failed is a full contribution — most of this repo is
those, deliberately.

## The bar

A claim enters the library only if it clears all of these. They exist because
an AI-generated mathematics repo that skips them is worse than nothing.

1. **Adversarial verification.** Anything claimed as a result must survive an
   independent attempt to refute it. For a proof step: a skeptic pass whose
   default stance is *refute*. For a computation: an independent
   re-implementation, not a re-run of the same code. Say who or what verified
   it and how.
2. **Reproducibility.** Every computational claim ships the command that
   produces it and the tool that produced it. A number nobody can regenerate is
   not evidence.
3. **Labelled speculation.** Every unproven step that the argument leans on is
   marked `SPECULATION` inline. Not in a caveats section at the end — at the
   step.
4. **Evidence is not proof.** A finite search is `EVIDENCE`, and the record
   states the range checked. `VERIFIED` describes the range, never the
   conjecture.
5. **Honest mode.** `mode: blind` or `mode: informed`, matching what you
   actually read. See `AGENTS.md`.
6. **Transcription is marked.** If a source detail came from a machine
   transcription of a paper rather than the paper, mark it and say so, as
   attempt 003 does with `[L]`.

Use the status vocabulary in `AGENTS.md` exactly: `VERIFIED`, `EVIDENCE`,
`LIVE`, `SPECULATION`, `REFUTED`, `MAP`.

## Submitting an attempt

1. Number it next in sequence for that problem:
   `problems/<problem>/attempts/NNN-<slug>.md`.
2. Use the required sections:
   - **front-matter** — problem, date, `mode`, type, tools used, sources
   - **Approach** — what you tried and why *this* rather than something else
   - **What was done** — enough that someone can follow and check it
   - **Outcome** — with a status term
   - **Why it failed / what survived** — the most valuable section. Be
     specific about the obstruction; "it didn't work" helps nobody
   - **Leads generated** — concrete, falsifiable next steps
   - **References**
3. Add an entry to `problems/<problem>/prior-art.json`, including:
   - `mechanism` — tags for the approach family. These are what lets a future
     agent check "has this family been tried?" without reading the prose.
     Reuse existing tags where they fit.
   - `leak_terms` — distinctive strings naming your *findings* (not tool
     vocabulary, not published results). CI uses these to keep tier 0 clean.
   - `gaps` — what blocks a proof, if the route is `LIVE`.
4. Put route-specific tooling in `problems/<problem>/explore/`. Only add to
   `harness/` if it verifies or enumerates the objects themselves — harness is
   tier 0, so a blind agent will see it.
5. Run `python -m pytest tests/ -q`. It fails if your record leaks into tier 0
   or if the index and `attempts/` disagree.

## Correcting an existing record

**Do not edit the original.** Records are kept as written, including their
errors — the corrections are part of the history and are evidence about how
the loop behaves.

Add a new attempt that reviews it, set `verifies` or `refutes` in the index,
and list the corrections there. `004-skeptic-review-of-003.md` is the worked
example: it found a false proof claim and a wrong constant in 003, and both
records stand with the relationship recorded between them.

## Adding a problem

`problems/<slug>/` with `PROBLEM.md` (published background only — no steer),
`PRIOR-ART.md`, and a `prior-art.json` with an empty `attempts` list. Add any
verifier to `harness/<slug>/`, and state the verification contract in
`PROBLEM.md`: what a claim about this problem has to do to be checkable.
