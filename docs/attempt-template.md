# NNN — <short title of the approach>

<!--
Scaffold this with:  python scripts/new_attempt.py <problem> <slug>

This is the ATTEMPT shape. If you are reviewing someone else's record rather
than making progress yourself, use the REVIEW shape instead -- see
CONTRIBUTING.md, and problems/union-closed/attempts/004-* as the worked
example. tests/test_records.py enforces whichever one you pick.

Every section below is required. Delete these comments as you fill it in.
-->

- **Problem:** <name>, `problems/<slug>/PROBLEM.md`
- **Date:** YYYY-MM-DD
- **Mode:** blind | informed
  <!-- blind = you worked from a scripts/blind.sh copy and did not read prior
       art. informed = you read any of it. Be honest; this field is data. -->
- **Type:** <survey | formalization | computational search | skeptic review | ...>
- **Tools:** <what you wrote or ran, with paths; note determinism and runtime>
- **Sources:** <papers consulted. Mark [T] anything taken from a machine
  transcription rather than the PDF itself. -->

## Approach

What you tried, and **why this rather than the obvious alternative**. If you
chose a specialized method over a general one, say what the general one would
have cost.

If informed: which prior attempt or gap this follows from.

## What was done

Enough that a reader can follow and check it. Include:

- the actual statements under test, precisely enough to be falsifiable;
- commands that reproduce every number quoted;
- cross-checks performed, and against what independent thing;
- `SPECULATION` labels inline at every unproven step the argument leans on.

## Outcome

Lead with a status term: `VERIFIED` / `EVIDENCE` / `LIVE` / `REFUTED` / `MAP`
(reviews that produce a machine-checked certificate use `FORMALIZED`).

State the **range or scope** — `VERIFIED` describes what was checked, never the
conjecture itself.

Then state explicitly **what is not claimed**. This section is where
overclaiming happens, so pre-empt it.

## Why it failed / what survived

The most valuable section in the file. Be specific about the obstruction: name
the step that breaks and the quantity that goes the wrong way. "It didn't work"
helps nobody.

If nothing failed, say what remains unproven and label each gap.

List what is reusable: tooling, closed forms, lemmas, hard-instance families.

## Leads generated

Numbered, concrete, falsifiable. Each should have a definite outcome either way.
"Investigate X" is not a lead; "compute Y and check the sign of Z" is.

## References

Papers with arXiv ids. Prior attempts in this repo by path. Mark machine
transcriptions.
