<!--
Most PRs here add an attempt record. Dead ends are welcome and expected --
an approach that failed, with a precise account of where it broke, is a full
contribution. Delete whichever sections do not apply.
-->

## What this adds

<!-- One or two sentences. Which problem, which approach, what came of it. -->

**Problem:**
**Mode:** blind / informed
**Status claimed:** VERIFIED / EVIDENCE / LIVE / REFUTED / MAP

## The bar

Full version in `CONTRIBUTING.md`. Tick honestly — an unticked box with an
explanation is fine, a wrongly ticked one is not.

- [ ] **Something independent tried to refute this.** A skeptic pass on proof
      steps, or a *re-implementation* of computations — not the same code run
      twice.
- [ ] **Every number quoted is reproducible**, with the command recorded.
- [ ] **`SPECULATION` is labelled inline**, at each unproven step the argument
      leans on.
- [ ] **`VERIFIED` describes a range, not the conjecture.** The range is stated.
- [ ] **`mode` is honest.** If prior art was read at any point, this is
      `informed`.
- [ ] Machine-transcribed sources are marked as such.
- [ ] No existing record was edited. Corrections are a new attempt with
      `verifies` / `refutes` set in the index.
- [ ] `prior-art.json` updated: `mechanism` tags, `status`, `gaps`,
      `leak_terms` naming this attempt's findings.
- [ ] `python -m pytest tests/ -q` passes.

## What is not claimed

<!-- The section most likely to be skipped and most worth writing. What would
     a hostile reader accuse this of overclaiming? Pre-empt it here. -->

## Obstruction / open gaps

<!-- If it failed: name the step that breaks and the quantity that goes the
     wrong way. If it is LIVE: what blocks a proof. -->
