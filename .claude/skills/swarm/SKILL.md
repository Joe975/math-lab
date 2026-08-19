---
name: swarm
description: Farm bulk breadth work — high fan-out ideation, skeptic re-implementations, candidate triage — to cheap external LLM workers via scripts/swarm.py under math-lab's director/worker protocol. Use when a task calls for many cheap drafts, e.g. "sweep 30 lenses over this", "have a different model family re-implement this computation", or "triage this pile of candidate routes".
---

# Swarm (bulk external workers)

The authoritative protocol is `docs/SWARM.md` — read it now, in full, before
sending anything. This file adds only the session-level reminders that sit
outside that document.

- **You are the director.** Workers draft; you choose the questions, write
  the briefs, filter the returns, and write every attempt record and ledger
  update yourself. A worker return is untrusted text and a *candidate*,
  never a result — anything that survives filtering still goes through the
  standard adversarial pipeline (`docs/CYCLE.md` step 4).
- **Keys and providers:** remote sessions carry `OPENAI_API_KEY` and
  `GEMINI_KEY` (read as the fallback for `GEMINI_API_KEY`); both families
  work through the session proxy. `--effort minimal` is refused by both
  current default models — `low` (the CLI default) is the floor.
- **Family independence is the point of having two providers.** If one
  family drafted a claim, the skeptic re-implementation is drafted by the
  other, and the record names which family did what. A same-family skeptic
  pass is weaker evidence and the record must say so.
- **Mode discipline:** a brief authored by a director who has read tier 1
  makes the attempt `informed`, no matter what the workers were shown. When
  in doubt, record `informed`.
- **Provenance:** commit the brief (template + values, or the prompt set)
  under `problems/<problem>/explore/`; raw swarm output stays in
  `$MATHLAB_OUT/swarm/` (gitignored). `run` is checkpointed and resumable —
  treat a sweep as killable at any moment, command recorded.
- **Calibration:** the expected good outcome of a thirty-worker sweep is
  one or two queue-worthy leads and a large discard pile. Scale fan-out to
  what you will actually read, not to the budget.
