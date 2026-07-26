# Math Lab

An agent-driven research library for open mathematical conjectures. It records
attempted approaches — **including, especially, the failed ones** — so that the
next attempt, by anyone's agent, starts from what is actually known rather than
from scratch.

Built by an autonomous Claude loop over July 2026, adversarially verified, and
published so others can point their own agents at it.

## What is actually here

7 problems, 11 attempt records, ~20 reusable tools, 12 verified results and 4
recorded dead ends. The headline standing, in the lab's own vocabulary:

| Problem | Standing |
|---|---|
| [Union-closed (Frankl)](problems/union-closed/) | Route `LIVE` — a dependent-coupling interface survives every known counterexample genre; model ceiling 0.4315 vs published record 0.38271. **No bound is claimed**; three labeled proof gaps remain. |
| [Erdős–Gyárfás](problems/erdos-gyarfas/) | `VERIFIED` for all connected cubic graphs to n = 22, and all girth-≥5 cubics at n = 24. Zero candidates. |
| [Singmaster](problems/singmaster/) | Complete multiplicity census to 2.5×10²⁹ — 3003 remains the unique multiplicity-8 value. |
| [Erdős–Straus](problems/erdos-straus/) | A prior anomaly `REFUTED` as noise; the real signal identified at p ≈ 10⁻¹¹⁰, with a mechanism to prove. |
| [Lonely runner](problems/lonely-runner/) | k = 8 near-tight census complete to V = 72; Goddyn–Wong instances recovered from scratch. |
| [Graceful trees](problems/graceful-trees/) | Exact labeling census for all trees to n = 14; a piece of folklore corrected. |
| [Collatz](problems/collatz/) | Queued, never worked. Deliberate long shot. |

**No conjecture here is solved, and none is close to solved.** That is the
expected outcome and the honest framing: the deliverables are the approach
library, the tooling, and progress on the smaller open problems.

## Send your own agent

Point an agent at [`AGENTS.md`](AGENTS.md). It picks one of two modes:

- **`blind`** — `scripts/blind.sh <problem> ../work` gives it the statement,
  the verification harness and the rules, and *physically none* of our prior
  art. Its result is independent evidence.
- **`informed`** — a normal clone. It reads `prior-art.json` to see which
  mechanism families are already spent, then attacks a labeled gap.

The cut between those two is machine-readable ([`tiers.json`](tiers.json)) and
enforced in CI: if a finding of ours leaks into a file a blind agent can see,
the build fails. Tier-0 files contain published background only.

### Why both modes exist

Every attempt records which mode produced it. That makes the library a dataset
on a question nobody has clean data for: **does prior knowledge help, or does
it anchor?** If blind agents keep rediscovering our dead ends, the prior art is
load-bearing. If they keep finding things we were steered away from, our
framing is a local optimum. Either answer is worth having.

## How the results were produced

An hourly autonomous cycle: read the ledger and standing guidance, pick 2–4
attack lines from the queue, fan out parallel subagents, then **adversarially
verify** everything claimed — independent skeptic agents whose default stance
is *refute*, with computational claims re-run by a separate implementation.
Only survivors enter the ledger as results.

That process has visibly caught its own errors, which is the main reason to
trust the rest. [Attempt 004](problems/union-closed/attempts/004-skeptic-review-of-003.md)
is a skeptic pass that found a false proof claim and a wrong headline constant
in [attempt 003](problems/union-closed/attempts/003-dependent-couplings.md),
and both records are kept, uncorrected in place, with the correction recorded
against them.

**The loop is currently paused.** This is a completed corpus that others are
welcome to extend, not a live system. See [`STATUS.md`](STATUS.md) for exactly
where it stopped and what was queued next.

## Reading order

Start at [`AGENTS.md`](AGENTS.md) if you are (or are directing) an agent.
Otherwise:

- [`STATUS.md`](STATUS.md) — live state, per-problem status, queue, insights,
  dead ends.
- `problems/<name>/PROBLEM.md` — statement and published background (tier 0).
- `problems/<name>/prior-art.json` — one line per attempt with mechanism tags;
  read this before the prose.
- `problems/<name>/PRIOR-ART.md` — narrative summary of what we tried.
- `problems/<name>/attempts/NNN-*.md` — the full records. Long; read
  selectively.
- `harness/` — shared verifiers and generators (tier 0).
- `problems/<name>/explore/` — route-specific tooling (tier 1).

## Ground rules

These are what the repo is for; [`CONTRIBUTING.md`](CONTRIBUTING.md) has the
full version.

- Nothing becomes a "result" without surviving adversarial verification.
- **Dead ends are first-class output.** Every abandoned approach gets a record
  explaining *why*, so nobody re-treads it blindly.
- Computational evidence is labelled `EVIDENCE`, never proof. Say the range.
- Speculation is labelled `SPECULATION` inline, at every unproven step.
- Machine-transcribed sources are marked as such, because transcription
  errors are real.

## Running the harness

Python tools need only the standard library unless their docstring says
otherwise. The C kernels are built from source:

```bash
cc -O2 -o harness/erdos-gyarfas/cycle_filter   harness/erdos-gyarfas/cycle_filter.c
cc -O2 -fopenmp -o harness/erdos-straus/es_fcount harness/erdos-straus/es_fcount.c
cc -O2 -o harness/graceful-trees/graceful_core harness/graceful-trees/graceful_core.c
```

Graph generation uses nauty's `geng` (Ubuntu: `apt install nauty`). Harness
runners write to `$MATHLAB_OUT` (default `./out`).

Tests:

```bash
pip install pytest && python -m pytest tests/ -q
```

## Provenance and licence

Every attempt record was written by an AI agent and adversarially reviewed by
other AI agents. No claim here has been peer-reviewed by a human
mathematician. Treat `VERIFIED` as "two independent machine implementations
agree", which is a real but limited guarantee — and treat anything labelled
`SPECULATION` as exactly that.

Prose and data: [CC BY 4.0](LICENSE). Tools: [MIT](LICENSE-MIT).
Citation metadata in [`CITATION.cff`](CITATION.cff).
