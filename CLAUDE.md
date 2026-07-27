# Working in this repo

This is a library of **attempted approaches** to open mathematical conjectures.
The failed attempts are the point, not a by-product. Read `AGENTS.md` before
doing anything else — it defines the two working modes and which files you may
read in each.

## Before you touch a problem, pick a mode

- **`blind`** — `scripts/blind.sh <problem> ../work` gives you the statement and
  harness with all prior art physically removed. Your attempt is then
  independent evidence.
- **`informed`** — read `problems/<problem>/prior-art.json` *first* (it is short),
  then only the specific attempt records it points you to. The records run to
  hundreds of lines; do not read them all.

Record which mode you used. A mislabelled mode corrupts the one dataset this
repo is uniquely able to produce.

## Hard rules

- **Never edit an existing attempt record**, including to fix an error in it.
  Add a new attempt that reviews it and set `verifies`/`refutes` in the index.
  The corrections are evidence about how the loop behaves.
- **Nothing is a result until something independent tried to refute it.** Proof
  steps get a skeptic pass; computations get re-run by a *different*
  implementation, not the same code again.
- **`VERIFIED` describes a range, never a conjecture.** A search to 10^6 is
  `EVIDENCE` about 10^6.
- **Label `SPECULATION` inline**, at the step that needs it.
- Tier 0 (`problems/*/PROBLEM.md`, `harness/**`) is published background only.
  Putting one of our findings there breaks blind mode and fails CI.

## Doing the work

- Running a research cycle: `docs/CYCLE.md` — the operational prompt.
- Starting an attempt: `python scripts/new_attempt.py <problem> <slug>`
  scaffolds the record and its index entry with the next number.
- Recording it: `docs/attempt-template.md`, schema in
  `docs/prior-art.schema.json`, full rules in `CONTRIBUTING.md`.
- Current state and the queue: `STATUS.md`.
- Cross-pollination (informed side only): `python scripts/mechanisms.py
  gaps <problem>` for untried field lenses; `docs/IDEATE.md` to sweep fields
  across one problem; `docs/RIPPLE.md` to scan a new result across problems.
  Session skills exist for each: `/ideate`, `/ripple`, `/mechanisms`,
  `/cycle`.

## Always

```bash
python -m pytest tests/ -q          # must pass before you commit
```

Python tools are standard-library only. C kernels build with `cc -O2`; see the
README. Harness runs write to `$MATHLAB_OUT` (default `./out`).

Commit each logical change with a clear message. Do not leave a large
uncommitted pile for the next agent.
