# Swarm: bulk external workers under a director

The lab's agent fan-out has always been priced in director-class attention:
`docs/IDEATE.md` caps a sweep at three to six lenses because each lens costs
a full agent. A cheap external reasoning model changes that arithmetic — a
thirty-worker sweep costs pennies. What it does *not* change is who is
accountable for the output. This file is the protocol for using bulk workers
without corrupting the dataset the repo exists to produce.

The tool is `scripts/swarm.py`; run `python scripts/swarm.py --help`.

## Providers

Two API families are wired in. The provider is inferred from the model name,
so in practice `--model` is the only flag that matters; `--provider` exists
for OpenAI-compatible proxies serving a Gemini model name.

| provider | models       | key                              | cheap tier as of 2026-08          |
|----------|--------------|----------------------------------|-----------------------------------|
| `openai` | `gpt-*`      | `OPENAI_API_KEY`                 | `gpt-5.6-luna` — $0.20/$1.20 per 1M |
| `gemini` | `gemini-*`, `gemma-*` | `GEMINI_API_KEY`, else `GEMINI_KEY` | `gemini-3.1-flash-lite` — $0.25/$1.50 per 1M |

Step-up tiers when a brief needs more than the floor: `gpt-5.6-terra`
($2/$12) and `gemini-3.7-flash` ($0.75/$3.75). `--effort` maps onto Gemini's
`thinkingLevel` under the same four names; not every model supports every
level (`minimal` is refused by `gemini-3.7-flash` with a 400 naming it).
Two operational notes from the Gemini side: `503 UNAVAILABLE / high demand`
is routine and is absorbed by the retry ladder, and thinking tokens bill at
the output rate, so `status` folds them into the output count.

**The two families are not interchangeable, and that is the point.** Which
family drafted a return is a fact about how independent that return is —
see the skeptic rule below.

## Division of labor

The **director** is the session model running the cycle. It:

- chooses the questions and writes the briefs (worker prompts);
- filters the returns — this is where the value is created;
- writes every attempt record and ledger update itself;
- routes anything that survives filtering into the standard adversarial
  pipeline (`docs/CYCLE.md`, step 4).

The **swarm** drafts. Workers see only their brief, have no repo access, no
state, and no authority. A worker's return is untrusted text — treat it
exactly like a lens report in an ideation sweep: a *candidate*, never a
result, and nothing in it is exempt from verification.

Never delegate to the swarm: line selection, attempt records, ledger
writes, or any claim that would enter the repo unverified.

## What to farm out

- **Ideation at high fan-out.** The `plan` subcommand exists for this:
  template = the lens brief from `docs/IDEATE.md` step 3, values = one line
  per lens (or per sub-question within a lens). The filter step of IDEATE.md
  is unchanged and is still the bottleneck — thirty lens reports are only
  worth generating if the director actually triages them against the
  mechanism index and the recorded obstructions.
- **Skeptic re-implementations.** The verification rule demands that
  computations be re-run by a *different implementation* (`docs/CYCLE.md`,
  step 4). Having a different model family write that re-implementation is
  strictly stronger independence than a second agent of the same family.
  Record which model wrote it in the verification note. With two providers
  available this is now a choice the director must make deliberately: if the
  original computation was drafted by a `gpt-*` worker, the skeptic
  re-implementation should be drafted by a `gemini-*` one, and vice versa.
  A same-family skeptic pass is weaker evidence and the record must say so.
- **Cross-family agreement as a triage signal, not a result.** Running the
  same brief on both families and keeping only what both return is a cheap
  filter on the discard pile. It is still a filter on *candidates*: two
  families agreeing on a false claim is a documented failure mode (001's
  false bijection claim was unanimous within one family), and agreement
  never substitutes for step 4.
- **Bulk triage.** Scoring a pile of candidate routes for "does this have a
  falsifiable first step" before the director reads them; classifying
  search-output anomalies; drafting variations on a construction.

Poor fits: anything that needs repo context beyond what fits in a brief,
anything where the answer will be trusted rather than checked, and proofs —
a swarm-drafted proof step is `SPECULATION` until the skeptic pass, like
everything else.

## Mode discipline

`mode` is data, and swarm work can silently corrupt it. The rule:

- An attempt is **`blind`** only if every word the workers saw traces to
  tier-0 files — e.g. `PROBLEM.md` verbatim plus generic instructions,
  assembled from a `scripts/blind.sh` checkout. A brief *authored* by a
  director who has read tier 1 carries the lab's framing into the prompt,
  and the attempt is `informed` no matter what the workers were shown.
- Sweeps built on the mechanism index or dead-end summaries are
  informed-side by construction, exactly as `docs/IDEATE.md` says.

When in doubt, record `informed`. A mislabelled blind is worse than no
attempt.

## Provenance

- The attempt record names the external model and effort next to `mode`
  (e.g. "workers: gpt-5.6-luna, effort low, 28 briefs"). When a record leans
  on more than one family, name each and say which did what — the split is
  the independence claim, so it has to be legible without opening the metas.
- If a record leans on swarm output, commit the brief — the template and
  values file, or the prompt set — under `problems/<problem>/explore/`.
  Swarm output itself lives in `$MATHLAB_OUT/swarm/` (gitignored working
  space); each job's `.meta.json` carries model, token usage, and the
  prompt hash, so a committed brief plus the meta is enough to re-run.
- Re-runs are cheap and checkpointed: `run` skips completed jobs and
  retries failures, so treat a swarm like any harness computation — killable
  at any moment, resumable, command recorded.

## Calibration

High fan-out mostly buys noise, deliberately: the expected good outcome of
a thirty-worker sweep is one or two queue-worthy leads, several recorded
"no purchase" verdicts, and a large discard pile. A swarm that appears to
have solved the problem has made an error — that is what step 4 is for. At
current prices generation is never the constraint; director filtering is.
Scale the fan-out to what the director will actually read, not to what the
budget allows.

The provider choice has its own calibration. The floor tiers of the two
families are within 25% of each other on price, so cost is not the reason to
prefer one; independence is. Spend the cheap tier on breadth (ideation,
triage, drafting variations) and reserve the step-up tiers for briefs where
a wrong answer would cost director attention to detect — and keep the family
that drafted a claim off that claim's verification.
