# Swarm: bulk external workers under a director

The lab's agent fan-out has always been priced in director-class attention:
`docs/IDEATE.md` caps a sweep at three to six lenses because each lens costs
a full agent. A cheap external reasoning model (currently `gpt-5.6-luna`,
$0.20/$1.20 per 1M tokens) changes that arithmetic — a thirty-worker sweep
costs pennies. What it does *not* change is who is accountable for the
output. This file is the protocol for using bulk workers without corrupting
the dataset the repo exists to produce.

The tool is `scripts/swarm.py`; run `python scripts/swarm.py --help`.

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
  Record which model wrote it in the verification note.
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
  (e.g. "workers: gpt-5.6-luna, effort low, 28 briefs").
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
