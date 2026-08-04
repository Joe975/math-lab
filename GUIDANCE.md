# Standing Guidance

Human-set direction for the lab. `docs/CYCLE.md` has agents read this at the
start of every cycle, and **it overrides the queue** where the two disagree —
so it has to describe how things actually are, not how they once were.

## Current guidance (2026-08-04)

Everything in the 2026-07-26 section below still stands — operating mode,
public-repo posture, and the calibration rules. This section adds direction
drawn from OpenAI's Astra announcement (2026-08-02): ten results on
decade-plus-open problems in mathematics and TCS, each published with a
model-written reasoning walkthrough and a Lean 4 certificate at zero sorries,
for a reported ~$2,000 of tokens. None are peer-reviewed yet; treat the
individual results as reported, not settled. The *shape* of the release is
what matters here, and three things follow for this lab:

1. **Run the ripple scan first.** Before the next cycle pulls anything else
   from the queue, run `docs/RIPPLE.md` on the Astra results (queue item 15).
   The Ehrhart-volume and extremal-graph/Ramsey entries are the plausible
   bites on mahler-4d and erdos-gyarfas; the scan decides, not this note.
2. **The novelty gate is now part of the bar.** The October 2025 "GPT-5
   solved ten Erdős problems" collapse was a labeling failure, not a
   derivation failure: correct arguments, already in the literature, claimed
   as new. Our blind mode *wants* rediscoveries — but recorded as
   rediscoveries. The skeptic pass now includes a literature check on any
   claimed-new result (`docs/CYCLE.md` step 4, `CONTRIBUTING.md` bar item 7).
3. **Formalization is the new top verification rung.** For skeptic-confirmed,
   load-bearing, proof-shaped steps, a kernel-checked Lean certificate is the
   one thing stronger than a second independent derivation — it is what made
   Astra's claims land after the October episode. `docs/FORMALIZE.md`
   defines the lane and the `FORMALIZED` status; queue item 16 is the pilot.
   Expensive; use deliberately.

**Portfolio note.** Astra's target list — decade-open but
specialist-tractable problems whose deliverables are checkable objects
(constructions, counterexamples, quantitative bound improvements) — is a
validation of this lab's existing tilt toward the smaller open problems, not
a reason to chase the famous ones harder. When onboarding new problems,
prefer that shape: a problem whose answer, if found, is an object a verifier
can check (the Erdős catalog is a good source). Famous long shots keep their
deliberate minority share.

## Previous guidance (2026-07-26)

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
