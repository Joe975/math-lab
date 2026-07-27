---
name: mechanisms
description: Query math-lab's cross-problem approach taxonomy — which approach families exist, where each has been tried, and which field lenses are untried on a problem. Use for questions like "has an entropy method been tried on X?", "what approaches haven't we tried on Y?", or "/mechanisms gaps union-closed".
---

# Mechanism taxonomy queries

`mechanisms.json` (repo root) defines every approach-family tag used in the
`problems/*/prior-art.json` indexes, each under a mathematical field lens.
`scripts/mechanisms.py` queries it joined with live usage:

```bash
python scripts/mechanisms.py list [--field FIELD]  # all tags by field, with usage
python scripts/mechanisms.py where <tag>           # attempts that used a tag, with status
python scripts/mechanisms.py gaps <problem>        # untried fields/tags for one problem
python scripts/mechanisms.py matrix                # problems x fields overview
```

Answer the user's question from these outputs plus, where needed, the
`one_line`/`gaps` entries in the relevant `prior-art.json` — open full
attempt records only if the index genuinely does not answer it.

Interpretation rules:

- A tag's absence means the *family* is untried there, not that it would
  work. Check the problem's dead ends before presenting a gap as a lead —
  some obstructions foreclose a whole field.
- `logic-methodology` tags are process machinery, not attack lenses; never
  offer them as approach ideas.
- Follow-ups route to the sibling skills: generating new routes for a
  problem is `/ideate`; propagating a new result across problems is
  `/ripple`.

Maintenance: every tag used in an index must be defined in
`mechanisms.json` with a field and description (`tests/test_mechanisms.py`
enforces this); new tags are registered when the first real attempt uses
them. This is all tier 1 — never surface it into `PROBLEM.md`, `harness/`,
or a blind working copy.
