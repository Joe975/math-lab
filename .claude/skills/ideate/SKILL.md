---
name: ideate
description: Run a cross-field ideation sweep on one math-lab problem — apply untried mathematical field lenses (one agent per lens) to generate new attack routes. Use when the user wants fresh approaches, cross-disciplinary ideas, or new queue leads for a specific problem, e.g. "/ideate union-closed" or "sweep other fields for ideas on collatz".
---

# Ideation sweep

You are running the field-sweep flow of math-lab's cross-pollination layer:
**many fields, one problem**. The authoritative procedure is `docs/IDEATE.md`
— read it now, in full, before doing anything. This file only orients you.

## Argument

The skill takes a problem slug (`ls problems/` for the list). If none was
given, ask which problem to sweep — do not pick one silently.

## Shape of the work

1. `python scripts/mechanisms.py gaps <problem>` and `matrix` — the untried
   lenses. Then the problem's `prior-art.json` and the `STATUS.md` dead-end
   section. Not the full attempt records.
2. Pick 3–6 lenses per `docs/IDEATE.md` §2.
3. One subagent per lens, in parallel, with the exact brief in
   `docs/IDEATE.md` §3. Give each agent only the inputs the brief lists —
   the narrow context is deliberate.
4. Filter (§4: spent-tag check against the index, falsifiable first step,
   not-a-dead-end-renamed), then record ONE attempt via
   `python scripts/new_attempt.py <problem> ideation-sweep-<slug>`, status
   `MAP`, `mode: informed`. Surviving leads go into the `STATUS.md` queue.
5. `python -m pytest tests/ -q`, commit, push.

## Hard guardrails

- Informed mode only. Never run this in a blind working copy, and the
  resulting attempt is always `mode: informed`.
- The sweep produces *candidate leads*, never results. No status stronger
  than `MAP`/`SPECULATION` may come out of it.
- "This field has no purchase here, because X" is a success — record it.
- Do not edit existing attempt records, and register new mechanism tags in
  `mechanisms.json` only when a real attempt first uses them.
