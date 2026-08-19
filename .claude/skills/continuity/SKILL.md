---
name: continuity
description: Open a bounded autonomous research window — work the STATUS.md queue continuously, with scheduled self-encouragement pings and no per-step human approval, for a human-set duration. Use when the user says e.g. "/continuity 3h", "keep at it for the afternoon", or "run a persistence window on the queue". Only ever on explicit human request — never self-invoked.
---

# Persistence window

Authorized by `GUIDANCE.md` (2026-08-19, "Persistence windows"). The
invocation is the human's request; inside the window, rolling from one
queue line to the next without asking is pre-authorized. The calibration
bar is unchanged — verification (`docs/CYCLE.md` step 4) still gates every
claim.

## Arguments

Duration, default **3h**; ping cadence, default **15m**. The human may
override either ("/continuity 2h 30m-cadence"). Never open a window without
an explicit invocation, and never extend one past its deadline.

## Opening the window

1. Note the UTC start and deadline in your working notes and to the user.
2. Arm the pings for the whole window up front, one message per cadence
   step until the deadline:
   - **Remote (claude.ai/code) sessions:** `send_later` (claude-code-remote
     MCP), one call per ping. Cron triggers won't work here — their minimum
     interval is hourly — and one-shot `send_later` messages have minute
     granularity, so arm N = duration/cadence of them. Cache note: this
     environment's prompt-cache TTL is 1 hour, so any cadence ≤ 60m resumes
     on a warm cache; 15m is comfortably inside.
   - **Local sessions:** `/loop 15m` with the ping text, stopping itself at
     the deadline.
3. Ping message template (fill in the deadline and index):

   > Window ping k/N — active until <deadline UTC>. Good momentum: keep
   > going. Continue the current line, or if it just closed, write it up,
   > commit, and pull the next queue item. Standards unchanged; checkpoint
   > to disk as you go. If genuinely blocked on human input, say so and
   > hold — otherwise do not wait for approval.

## During the window

- Work the queue per `docs/CYCLE.md`, but rolling: line → write-up →
  commit → next line, without stopping to ask.
- Keep units small and committed — assume the process can die at any
  moment. A ping that arrives mid-line is license to continue, not a
  status-report request.
- Verification is not deferred to the end of the window: skeptic passes
  ride along with the claims they check, as usual.

## Closing the window

At the deadline (or when the human closes it early):

1. Finish the current write-up — do not start a new line.
2. Update `STATUS.md`, run `python -m pytest tests/ -q`, commit, push.
3. Delete any still-pending pings (`list_triggers` → `delete_trigger`).
4. Report what the window produced, then stop.
