# Standing Guidance

Human-set direction for the lab. `docs/CYCLE.md` has agents read this at the
start of every cycle, and **it overrides the queue** where the two disagree —
so it has to describe how things actually are, not how they once were.

## Current guidance (2026-07-26)

**Operating mode.** No continuous loop is running. Work happens on demand: a
cycle is run deliberately, by this lab or by an outside contributor's agent,
rather than on a timer. There is no orchestrator waiting to collect results, so
an agent that finishes a line should write it up, update the ledger and stop —
not roll on into the next queue item unprompted.

**The repository is public.** Everything written here is published
immediately. That does not change the standard, since these records were always
meant to be read by strangers, but it does mean overclaiming is now a public
act. The calibration rules below are load-bearing, not stylistic.

**Direction.**

- Rebalance effort toward whatever generates real traction. The budgets in the
  problem-status table are guidance, not quota.
- Famous long shots (Collatz) keep a minority share, deliberately.
- Keep per-cycle cost modest: many small verified steps beat rare huge
  fan-outs. Depth here comes from accumulation across cycles.
- The `STATUS.md` TL;DR must always reflect current state, so someone dropping
  in cold is oriented by one paragraph.

**Calibration — unchanged and non-negotiable.**

- Nothing becomes a result without surviving an independent attempt to refute
  it.
- Dead ends are the expected output and a full contribution.
- `VERIFIED` describes a range. `EVIDENCE` is not proof. `SPECULATION` is
  labelled inline.
- An agent reporting a breakthrough has almost certainly made an error.

## Changing this file

Update it when the direction actually changes, and date the change. Stale
guidance is worse than none, because an agent will follow it in preference to
the queue.
