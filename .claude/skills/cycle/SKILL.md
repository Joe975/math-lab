---
name: cycle
description: Run one research cycle of the math-lab loop — orient via STATUS.md and GUIDANCE.md, work 2-4 queue lines in parallel, adversarially verify, write attempt records, update the ledger. Use when the user says "run a cycle", "pick up the queue", or "continue the research".
---

# Run a research cycle

The operational prompt is `docs/CYCLE.md` — read it now, in full, and follow
its six steps exactly. This file adds only the session-level reminders that
sit outside that document.

- **Mode discipline:** a normal-clone cycle is informed mode. If a line is
  meant to produce independent evidence, materialize it with
  `scripts/blind.sh <problem> <dir>` and keep that agent inside the blind
  copy. Record `mode` honestly in every attempt — it is the dataset.
- **Cross-pollination hooks:** during orient, if a new external result or a
  newly `LIVE`/`REFUTED` internal route has appeared since the last cycle,
  consider a ripple scan (`docs/RIPPLE.md`, `/ripple`). When choosing lines
  for a problem with a thin queue, `python scripts/mechanisms.py gaps
  <problem>` shows untried lenses; a full sweep is `/ideate`.
- **Calibration:** the expected outcome is a well-recorded dead end, and
  that is a success. A breakthrough claim is almost certainly an error —
  step 4 (adversarial verification) exists for exactly that reason, and
  nothing enters the ledger without it.
- **Before ending:** every pursued line has an attempt record and index
  entry, `STATUS.md` TL;DR reflects the new state, `python -m pytest
  tests/ -q` passes, and the work is committed and pushed. Stop after the
  cycle — do not roll into the next queue item unprompted (standing
  guidance in `GUIDANCE.md`) — unless a human-opened persistence window
  (`/continuity`, `GUIDANCE.md` 2026-08-19) is active, in which case
  rolling on inside the window is pre-authorized.
