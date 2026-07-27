---
name: ripple
description: Scan a new mathematical result (external paper/preprint, or an internal LIVE/VERIFIED/no-go finding) across all math-lab problems to find where it unblocks a recorded gap, kills a queued route, or opens an untried lens. Use when the user mentions a new paper, breakthrough, or lab result and wants its consequences propagated, e.g. "/ripple <result>" or "does this preprint help any of our other problems?".
---

# Ripple scan

You are running the propagation flow of math-lab's cross-pollination layer:
**one result, many problems**. The authoritative procedure is
`docs/RIPPLE.md` — read it now, in full, before doing anything. This file
only orients you.

## Argument

The skill takes a result: a paper/preprint reference, or a pointer to an
internal attempt record. If none was given, ask what result to propagate.

## Shape of the work

1. Characterize the result per `docs/RIPPLE.md` §1: what it establishes
   (actual scope — verify external sources really say it), its mechanism as
   `mechanisms.json` tags, and the *hypothesis list* — the structural
   features a problem must have for the argument to apply.
2. Scan EVERY other problem (§2): `python scripts/mechanisms.py where <tag>`,
   then each `prior-art.json` — gaps first, queue second, untried lenses
   third. One line per miss, a paragraph per hit.
3. Record (§3): misses → dated note in `STATUS.md` insights; hits → queue
   items citing the result (gap-unblocking hits go top of that problem's
   queue); a scan that changed the picture → a `MAP` attempt on the problem
   most affected.
4. `python -m pytest tests/ -q`, commit, push.

## Hard guardrails

- Informed mode only; the external result itself may enter tier-1 prose
  only. Never write anything of ours into `PROBLEM.md` or `harness/`.
- A hit must name the specific recorded gap or queue item it meets AND why
  the result's hypotheses hold there. "Morally applies" is not a hit —
  analogy inflation is this flow's known failure mode.
- Transfer hypotheses get `SPECULATION` labels inline like everything else.
- Killing a queue item on the strength of a no-go is in scope; killing or
  editing an attempt *record* is never in scope.
