# Math Lab

An agent-driven research loop that attacks open mathematical conjectures and
builds a library of attempted approaches — including (especially) the failed
ones. Run autonomously by Claude via an hourly scheduled cycle, steered by
high-level human guidance.

## How the loop works

Every cycle (hourly), the orchestrating session:

1. Pulls the latest branch state and reads `LEDGER.md` (current state) and
   `GUIDANCE.md` (standing human guidance), plus any new chat messages.
2. Picks a small number of attack lines (typically 2–4) from the queue.
3. Fans out parallel subagents — e.g. computational counterexample searches,
   proof-sketch attempts from a named angle, structure surveys.
4. **Adversarially verifies** anything an agent claims: independent skeptic
   agents try to refute any claimed proof step or counterexample before it
   may enter the library as a result. Computational claims are re-run
   independently.
5. Writes a structured attempt record for every line pursued (success or
   dead end), updates `LEDGER.md`, commits, and pushes.

Cycles are deliberately modest (a few agents each) so the loop can run for a
month. Depth comes from accumulation and from the ledger steering later
cycles away from known dead ends.

## Layout

- `LEDGER.md` — live state: TL;DR, per-problem status, attempt queue, insights.
- `GUIDANCE.md` — standing human guidance; updated when the human steers.
- `problems/` — one file per conjecture: statement, known results, attack surface.
- `attempts/<problem>/NNN-<slug>.md` — one record per attempt. Required
  sections: Approach, What was done, Outcome, Why it failed / what survived,
  Leads generated.
- `tools/` — reusable computational scripts (searches, verifiers). Scripts are
  kept re-runnable so results can always be independently reproduced.

## Ground rules

- Nothing enters the ledger as a "result" without surviving adversarial
  verification.
- Dead ends are first-class output: every abandoned approach gets a record
  explaining *why*, so no future cycle re-treads it blindly.
- Computational evidence is labelled as evidence, never as proof.
- Honest calibration: famous conjectures are long shots; the deliverables are
  the approach library, the tooling, and possible progress on the smaller
  open problems.
