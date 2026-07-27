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

### 4. Adversarially verify — the step that makes this worth doing

Nothing enters the ledger as a result because an agent claimed it.

For every load-bearing claim, spawn a **skeptic whose default stance is
refute**. Not "check this" — *try to kill it*. The skeptic must:

- re-derive proof steps independently rather than reading and agreeing;
- **re-implement** computations from scratch, not re-run the same code. Two runs
  of one program agreeing tells you nothing about whether the program is right;
- check that headline constants come from the actual extremum, not from where a
  search grid happened to stop;
- verify that cited sources say what they are claimed to say.

If the skeptic finds something, that finding is recorded, and the original
record is **left as written**. See `problems/union-closed/attempts/004-*` for
what a good skeptic pass looks like: it confirmed the main interface, refuted a
mini-theorem, corrected a headline constant, and flagged a data file that could
not be reproduced.

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
tooling, and progress on the smaller problems. An agent that reports a
breakthrough has almost certainly made an error, which is exactly why step 4
exists.
