# Ideation sweep: many fields, one problem

The first cross-pollination flow: pick **one problem** and sweep the
mathematical field lenses across it, one agent per lens, to generate attack
routes the queue does not already contain. The output is a `MAP` attempt — a
survey with candidate leads — not a claim of progress.

Run it when a problem's queue is thin, when its live routes all sit in one or
two fields, or when a route just died and the reflex "try a variant of the
same thing" needs to be resisted.

This is an **informed-side** procedure by construction: its whole value is
reading what has been tried. Never run it in a blind working copy, and any
attempt that came out of a sweep is `mode: informed`.

## The procedure

### 1. Compute the gap profile

```bash
python scripts/mechanisms.py gaps <problem>
python scripts/mechanisms.py matrix
```

That lists the field lenses never tried on this problem and the tags proven
useful elsewhere. Then read, in order: `problems/<problem>/prior-art.json`
(statuses, `gaps`, dead ends), the *Dead ends* section of `STATUS.md`, and
`problems/<problem>/PROBLEM.md`. Do not read full attempt records yet — the
sweep should generate ideas before anchoring on the prose.

### 2. Choose the lenses

Take every untried field from the gap output, minus any the dead-end record
already forecloses. Include one or two *tried* fields if the tags used there
are a narrow slice of the field. Three to six lenses is the useful range —
below that it is not a sweep, above it the write-ups get thin.

### 3. One agent per lens, in parallel

Each agent gets the same brief:

- **Inputs:** `PROBLEM.md`, the field lens (name + description from
  `mechanisms.json`), the list of mechanism tags already spent on this
  problem, and the one-line dead-end summaries. Not the full records.
- **Task:** propose 1–3 attack routes on this problem *from inside this
  field*. For each route: the core object or quantity the field would look
  at, the first concrete calculation or search, and what result would kill
  the route ("no falsifiable first step" disqualifies an idea).
- **Standard:** honest triage over enthusiasm. "This field has no purchase
  here, and here is why" is a valid and useful lens report.
- Label every unproven assumption `SPECULATION` inline, as always.

Known-transfer check, per agent: does an existing tag from another problem
(`mechanisms.py where <tag>`) instantiate the route directly? Reusing a tool
that exists beats inventing one.

### 4. Filter and record

Collect the lens reports. Discard routes that:

- match a `mechanism` tag already spent here (check the index, not memory);
- have no falsifiable first step;
- are a known dead end wearing new vocabulary — compare against the
  *obstruction*, not the name. This is the failure mode of the whole
  exercise, so check it deliberately.

Record the sweep as one attempt (`python scripts/new_attempt.py <problem>
ideation-sweep-<something>`), status `MAP`, with one subsection per lens —
including the "no purchase" verdicts, which save the next sweep from
re-running that lens. Surviving routes go into the `STATUS.md` queue in the
standard lead format: concrete, falsifiable, with a definite outcome either
way.

### 5. Register vocabulary

Any genuinely new approach family the sweep produced gets a tag in
`mechanisms.json` (field + description) when its first real attempt is
recorded — not before. The sweep itself is tagged with the lenses it swept
plus `barrier-analysis` where a lens verdict was "no purchase, and here is
the obstruction".

## Calibration

A sweep that produces six exciting routes has probably produced six shallow
ones. The expected good outcome is one or two queue-worthy leads and several
recorded "no purchase" verdicts. Lens agents see less context than a full
attempt would, deliberately — their routes are *candidates for the queue*,
never results, and nothing they produce is exempt from the adversarial
pipeline that everything else goes through.
