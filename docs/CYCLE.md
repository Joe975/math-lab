# The research cycle

The operational prompt for running a cycle of this lab. `STATUS.md` promises
that everything needed to resume is in this directory; this file is the part
that promise depends on.

Give an agent this file plus "run a cycle". Or run it yourself.

## One cycle

### 1. Orient

Read, in this order:

- `STATUS.md` — current standing, the attempt queue, insights, dead ends.
- `GUIDANCE.md` — standing human direction. It overrides the queue.
- `git log --oneline -15` — what the last cycles actually did.

If guidance and queue disagree, guidance wins; note the divergence in the
ledger.

If something *new* has landed since the last cycle — an external paper
touching a problem or a mechanism family, or an internal route newly `LIVE`,
`VERIFIED`, or killed by a general no-go — consider a ripple scan
(`docs/RIPPLE.md`) before choosing lines: it is cheap, and a result that
unblocks a recorded gap or forecloses a queued route changes what is worth
running this cycle.

### 2. Choose two to four lines

Pull from the top of the `STATUS.md` queue, adjusted for guidance. Keep cycles
small — depth in this lab comes from accumulation across many cycles, not from
one large fan-out. Prefer:

- a line with a **falsifiable first step** over one that needs a whole theory;
- a line whose failure would be *informative* over one whose failure teaches
  nothing;
- balance across problems, so a single hard problem does not eat the budget.

Give long shots (Collatz here) a minority share, deliberately.

When a problem's queue is thin, or its routes all live in one or two fields,
`python scripts/mechanisms.py gaps <problem>` shows which field lenses are
untried there; a full ideation sweep (`docs/IDEATE.md`) is itself a valid
line for a cycle, and produces a `MAP` attempt.

### 3. Work the lines in parallel

One agent per line. Each agent must:

- state its approach and **why this rather than the obvious alternative**;
- keep every computation re-runnable, with the command recorded;
- label `SPECULATION` inline at each unproven step it leans on;
- stop and write up when it hits an obstruction, rather than thrashing. An
  obstruction found and described precisely is the deliverable.

Checkpoint completed work units to disk immediately. Treat "this process can
die at any moment" as a design assumption — a long run that loses everything on
restart is worse than a short one that does not.

### 3b. Push, before you allow anyone to write up

Step 3 tells an agent to stop at the obstruction. That is right for a proof
step and wrong for an object, and conflating the two costs results.

Before accepting a partial result, ask: **is the target something a checker
could settle?** A counterexample, a graph, an integer, a point set — anything
that can be handed to an exhaustive verifier in exact arithmetic. If so:

- the agent builds the checker first, then keeps going;
- a conditional or nearly-complete construction is **not** a stopping point.
  It is the most dangerous kind of result, because it reads like progress and
  ends the search;
- push again for the complete object, and keep pushing while a *structurally
  different* idea remains. A retry of the same search is not a push;
- stop when the distinct ideas run out, not when patience does.

The failure mode this prevents is real and recent: in July 2026 a ~30-year-old
conjecture in flow theory fell to a seven-node counterexample that a model
produced only on the fourth push. Its first three answers were exactly the
"precisely-described obstruction" this document otherwise rewards.

None of this touches step 4. Ambition is affordable *because* verification is
cheap for these targets, so the bar does not move to accommodate it.

Record `push_rounds` and `target_shape` in the index entry. Whether results
arrive on the first try or only under persistence is a second dataset this
library can produce, alongside blind-vs-informed.

### 4. Adversarially verify — the step that makes this worth doing

Nothing enters the ledger as a result because an agent claimed it.

For every load-bearing claim, spawn a **skeptic whose default stance is
refute**. Not "check this" — *try to kill it*. The skeptic must:

- re-derive proof steps independently rather than reading and agreeing;
- **re-implement** computations from scratch, not re-run the same code. Two runs
  of one program agreeing tells you nothing about whether the program is right;
- check that headline constants come from the actual extremum, not from where a
  search grid happened to stop;
- verify that cited sources say what they are claimed to say;
- **check claimed novelty against the literature.** Search for the statement
  and its neighborhood before any record calls a finding *new*. Rediscovering
  a published result is a fine outcome — blind-mode rediscovery is calibration
  data this lab wants — but it is recorded as a rediscovery, with the
  citation. The failure mode is the October 2025 "GPT-5 solved ten Erdős
  problems" episode: real derivations, already in the literature, announced
  as new.

If the skeptic finds something, that finding is recorded, and the original
record is **left as written**. See `problems/union-closed/attempts/004-*` for
what a good skeptic pass looks like: it confirmed the main interface, refuted a
mini-theorem, corrected a headline constant, and flagged a data file that could
not be reproduced.

For a proof-shaped claim that survives the skeptic and is load-bearing, there
is one rung above this: a machine-checked certificate (`docs/FORMALIZE.md`,
status `FORMALIZED`). Optional and expensive — reserve it for steps many later
records will lean on, and read the bridge caveat there first.

### 5. Write it up

Every line pursued gets a record, whether it worked or not:

```bash
python scripts/new_attempt.py <problem> <slug>
```

Fill in the scaffold (`docs/attempt-template.md` explains each section) and
complete the `prior-art.json` entry — `mechanism` tags, `status`, `gaps`, and
`leak_terms` naming your findings so CI can keep them out of tier 0. Reuse
`mechanism` tags from `mechanisms.json` where they fit; a genuinely new tag
gets an entry there (field + method-level description), and
`tests/test_mechanisms.py` fails until it does.

The **"Why it failed / what survived"** section is the most valuable thing you
will write. Be specific about the obstruction. "It didn't work" helps nobody;
"KL charges escaping mass by log-likelihood while the entropy drop is Θ(n)"
stops the next agent dead before they waste a cycle.

### 6. Update the ledger and commit

Update `STATUS.md`: TL;DR, per-problem status, queue (remove what was done, add
leads generated), insights, dead ends. Then:

```bash
python -m pytest tests/ -q
git add -A && git commit && git push
```

The TL;DR must always reflect current state, so that someone dropping in cold
is oriented in one paragraph.

## Choosing what to queue next

Good leads are concrete and falsifiable. "Investigate couplings" is not a lead.
"Expand E[h(z_ρ)] − (1/2φ)(x h(y) + y h(x)) to O(ε²) at the AHS equality point
and check the sign" is a lead — a specific calculation with a definite outcome
either way.

When a route dies, ask what the *obstruction* rules out, not just what failed.
A no-go that covers a whole family of functionals is worth more than a single
refutation.

## Calibration

These are famous open problems. The expected outcome of any cycle is a recorded
dead end, and that is a success. The deliverables are the approach library, the
tooling, and progress on the smaller problems.

An agent that reports a breakthrough **it cannot hand to a checker** has almost
certainly made an error, which is exactly why step 4 exists. An agent that
hands you an object a checker accepts is a different situation entirely, and
the correct response is to verify it rather than to disbelieve it: run an
independent implementation, and if it holds, it holds. Calibration is about
where to place doubt, not about capping what an attempt is allowed to aim at.

Two failure modes, not one. Overclaiming is the obvious one and step 4 catches
it. Under-reaching is the quiet one — stopping at a partial construction that
a few more pushes would have completed — and nothing catches it, because a
tidy dead-end record looks identical whether or not the answer was two ideas
away. Step 3b exists for that reason.
