# 002 — Graph-structure census of truncated Collatz digraphs

- **Problem:** Collatz, `problems/collatz/PROBLEM.md`
- **Date:** 2026-08-15
- **Mode:** informed (follows 001 lead 1 = `STATUS.md` queue item 21; 001 was
  read in full, and the swarm brief below was authored from it)
- **Type:** computational census + structural no-go
- **Tools:** `problems/collatz/explore/collatz_graph_census.py` (stdlib,
  deterministic; `degeneracy` 0.4 s to B = 2^18, `interface` 7 s to k = 22,
  `verify-mincut` seconds to k = 12). Cross-family skeptic pass via
  `scripts/swarm.py` — **3 workers, gemini-3.7-flash, effort high**, 1,334
  tokens in / 49,744 out, $0.188, brief committed under
  `explore/graph-census-skeptic-2026-08/`. Provider chosen against the rule
  in `docs/SWARM.md`: the census was written by the session model, so its
  verification went to the other model family.
- **Sources:** literature pre-step done by web search only, no PDFs — every
  attribution here is [T]. Wirsching's monograph and Ebert (arXiv:1905.07575)
  are the standard references for the inverted Collatz graph being a rooted
  tree; no dominator-theoretic or treewidth-theoretic Collatz result was
  found.

## Approach

Queue item 21 asked for three unlabelled invariants on the truncated digraph
`G_B` — dominator tree from 1, minimum directed cuts separating `[1, 2^k]`
from orbits exceeding `2^{k+1}`, and treewidth via flow-cutter — computed for
`B = 2^12 … 2^24` and tracked for growth in `B`, with the falsifiable
dichotomy **interface size `|S_k|` bounded vs growing**.

The pre-step 001 required (its gap 2) was a novelty check on the
dominator/treewidth framing. It found no such work, but it also surfaced the
classical fact that the inverted Collatz graph is a rooted **tree**
(Wirsching; Ebert [T]). That inverts the reading of the novelty flag: the
likelier explanation for an empty literature is that these invariants are
trivial on this graph, not that nobody thought of them.

So the census was run smallest-first, with each invariant's degeneracy stated
as a claim and checked by a real algorithm, before any flow-cutter machinery
was built. The alternative — building the `B = 2^24` pipeline first, as the
lead specified — would have cost a cycle to reach the same place. It turned
out that scale is not merely unnecessary here but *meaningless*: see C3.1.

## What was done

**The object.** `G_B` has vertices `{1, …, B}` plus an escape sink, and one
edge per vertex: `v → T(v)` when `T(v) ≤ B`, else `v → sink`, with `T` the
tier-0 map. **Every vertex has out-degree exactly one.** That single fact
drives everything below.

Reproduce with:

    python problems/collatz/explore/collatz_graph_census.py degeneracy --b-exp 10 12 14 16 18
    python problems/collatz/explore/collatz_graph_census.py interface --k 2 22
    python problems/collatz/explore/collatz_graph_census.py verify-mincut --k-max 12

**C1 — the dominator tree is the orbit tree.** Forward orbits are unique, so
in the reverse graph rooted at 1 every vertex has in-degree 1 and there is
exactly one root→v path; `idom(v)` must be `T(v)`. Checked rather than
assumed, with a Cooper–Harvey–Kennedy iterative dominator computation:
`idom(v) = T(v)` for **every** vertex, zero mismatches, at
`B = 2^10, 2^12, 2^14, 2^16, 2^18` (up to 104,575 vertices reaching 1).
Data: `data/graph_census_degeneracy.json`.

**C2 — the treewidth is 2, for every B.** `|V| = B+1`, `|E| = B`. A vertex
reaches either 1 or the sink but never both, so those are separate
components: the sink component is a tree, and the 1-component is unicyclic
(the `1→4→2→1` triangle). A forest plus one cycle has treewidth ≤ 2, and the
triangle forces ≥ 2. Checked by degree-≤2 elimination: the graph eliminates
completely with maximum elimination degree 2 at every `B` above. Constant in
`B`, so "growth in B" cannot be measured — there is nothing growing.

**C3 — the minimum cut.** `|S_k|` = the fewest edges whose deletion leaves no
path from `[1, 2^k]` to anything above `2^{k+1}`.

- **C3.1 — `B` drops out entirely.** An orbit from `[1, 2^k]` stays at or
  below `2^{k+1}` until the step that crosses, so every cut edge lies inside
  `[1, 2^{k+1}]` and `B` never enters. The lead's "track growth in `B`" is
  vacuous, not merely uninformative.
- **C3.2 — crossings are odd.** From even `v`, `T(v) = v/2 ≤ 2^k`, never a
  crossing. So every candidate is an odd `r` with `3r + 1 > 2^{k+1}`.
- **C3.3 — crossings reachable from `[1, 2^k]` lie in `[1, 2^k]`.** Let `c` be
  odd with `c > 2^k`. A predecessor `u` of `c` is either odd, giving
  `T(u) = 3u+1` even ≠ `c`, or even, giving `u = 2c > 2^{k+1}` — outside the
  window. So `c` has **in-degree zero**: it is unreachable from anywhere.
  *(This is the skeptic's argument, cleaner than the inductive one the census
  was built on; see the verification note.)*
- Hence `S_k = {odd r ≤ 2^k : 3r + 1 > 2^{k+1}}`, those edges are pairwise
  edge-disjoint one-step paths so no smaller cut exists, and

      |S_k| = #{odd r : (2^{k+1} − 1)/3 < r ≤ 2^k} = (2^{k−1} − (−1)^{k−1}) / 3

  the Jacobsthal number `J_{k−1}`. Asymptotically `|S_k| ~ 2^k / 6`.

**Cross-checks, all independent of each other:**

1. A linear-time memoized first-crossing census computes `|S_k|` for
   `k = 2 … 22` (`data/graph_census_interface.json`). The closed form is
   **exact at every k**, and the computed cut *set* equals the predicted
   one-step set exactly at every `k` (`predicted_set_matches`).
2. A unit-capacity Dinic max-flow on the explicit graph — a different
   algorithm with no shared code path — agrees exactly for `k = 2 … 12`
   (`data/graph_census_mincut_check.log`).
3. **Cross-family skeptic pass**, reproduced with:

       python scripts/swarm.py plan explore/graph-census-skeptic-2026-08/jobs \
         --template problems/collatz/explore/graph-census-skeptic-2026-08/template.md \
         --values   problems/collatz/explore/graph-census-skeptic-2026-08/values.txt
       python scripts/swarm.py run explore/graph-census-skeptic-2026-08/jobs \
         --provider gemini --effort high --max-output-tokens 40000 --workers 3

   Three `gemini-3.7-flash` workers were given
   only the definitions (no algorithm, no numbers, no closed form) under three
   different methodological stances, and asked for self-contained programs.
   All three programs were reviewed, then run: **all three agree with the
   census on every `k = 2 … 14`.** Worker 001 wrote its own Dinic and used the
   max-flow min-cut theorem with no graph-specific shortcut. Worker 002
   independently derived the out-degree-one argument (distinct exits give
   edge-disjoint paths, so the exit count *is* the minimum). Worker 003
   independently derived the same closed form `J_{k−1}`, with the in-degree
   argument of C3.3, which is **better than what the census was built on** and
   is adopted above.
4. The escape fraction quoted under "what survived" was re-derived by
   unmemoized brute force at `k = 6, 8, 10, 12, 13, 14`, matching the census
   exactly (including the non-monotone drop at `k = 13`).

**Residue composition** (001's "keep the arithmetic labels" instruction). At
`k = 22`: the cut is 100% odd, and splits 233,017 / 233,017 / 233,017 across
residues mod 3 — exactly uniform, as any interval of odd integers must be.
The labels are carried and they show nothing.

## Outcome

**REFUTED** — the graph-structure route of 001 lead 1 / queue item 21 is dead,
and the dichotomy it was built on resolves in the degenerate direction.

Scope, stated precisely:

- C1 and C2 are proved in general and machine-checked for `B ≤ 2^18`.
- C3's closed form is **proved for all `k ≥ 2`** (C3.1–C3.3 is a complete
  argument, and it survived a cross-family skeptic pass that reproduced it
  independently), and machine-checked at `k = 2 … 22` by the census, `k ≤ 12`
  by max-flow, and `k ≤ 14` by three external implementations.
- `|S_k|` is not merely growing: it is `Θ(2^k)` with an exact closed form, a
  constant `1/6` fraction of the window.

**Not claimed.** Nothing whatsoever about the Collatz conjecture. These are
facts about a finite digraph; they constrain what a *method* can see, not what
is true. In particular C2 does not say Collatz is "simple" — it says treewidth
is measuring the wrong thing. No claim is made about labelled or weighted
objects built on the same digraph, which this attempt does not touch.

## Why it failed / what survived

**The obstruction, named exactly.** The Collatz map is a function, so `G_B`
has out-degree one everywhere, so it is a forest plus one triangle. Dominator
trees, treewidth and minimum cuts all measure *branching complexity*, and a
forest has none. These invariants therefore cannot depend on the `3n+1`
arithmetic at all — they are determined by the shape, and the shape is fixed
before any arithmetic is consulted.

This is 001's "unlabelled-invariant barrier" (worker 014), but sharper, and
the sharpening is the deliverable. 001 recorded that label-free structure has
"weak purchase". The truth is stronger and provable: two of the three
invariants are **constant by construction**, and the third is decided by the
single inequality `3r + 1 > 2^{k+1}` — pure magnitude, no dynamics. The
quantity the lead proposed to measure the dynamics with is a counting problem
about an interval of odd integers.

Keeping the arithmetic labels, as 001 explicitly instructed, does not rescue
it: the cut is a contiguous block of odds, so every residue statistic is
uniform *by construction*. There is no label distribution to look at.

**Reusable:**

- `collatz_graph_census.py` — census, max-flow cross-check, dominator and
  elimination checks. The Dinic and the elimination routine are generic.
- The closed form `|S_k| = J_{k−1}` and the C3.1–C3.3 argument.
- **The pricing argument generalizes**, and this is the part worth carrying
  to other problems: for *any* iterated map, the digraph is functional, so
  this entire family of unlabelled invariants collapses the same way. That is
  a one-paragraph check that prices a graph-decomposition route before any
  machinery is built. Recorded as the new mechanism tag
  `functional-graph-invariants`.
- One genuinely dynamical quantity fell out of the census as a by-product:
  the **escape fraction** `#{n ≤ 2^k : orbit exceeds 2^{k+1}} / 2^k`. It is
  non-monotone and does not settle where a naive orbit-max heuristic predicts
  — 0.359 (k=6), 0.461, 0.570, 0.609 (k=12), then a **drop to 0.379** at
  k=13, after which it sits in 0.401–0.406 through k=22. Brute-force
  confirmed. This is not explained here.

**What this does not kill:** labelled or weighted objects on the same digraph,
where the arithmetic rides in the weights rather than the shape — transfer
operators, and 001 lead 3's tilted moments. The collapse is a statement about
unlabelled invariants only, and it is not evidence about those.

## Leads generated

1. **Close the graph family properly, cheaply.** The only version of a graph
   route not killed above is one on a digraph with out-degree > 1 — i.e. the
   *reverse* map, where `n` has preimage `2n` always and `(n−1)/3` when
   `n ≡ 4 (mod 6)`. Falsifiable first step: does the branching set have any
   density beyond what that single congruence forces? Compute the density of
   binary-branching vertices in the reverse tree restricted to `[1, 2^k]` for
   `k ≤ 24` and compare to `1/3`. **Expected kill** (the congruence is the
   whole story, so it should be `1/3` to noise); if it is *not* `1/3`, the
   shape carries arithmetic after all and the family reopens. Either way the
   graph family is then closed with a reason rather than a shrug.
2. **Explain the escape-fraction anomaly** (from the by-product above). Map
   `#{n ≤ 2^k : max orbit > 2^{k+1}} / 2^k` for `k ≤ 26` and find the
   mechanism for the `k = 12 → 13` drop and the plateau near 0.404.
   Falsifiable: either the drop is a window-alignment artifact with an exact
   description (most likely — `2^{k+1}` crossing a specific orbit-max cluster),
   or the plateau is a real density and the limit, if it exists, is a
   statement about stopping-time distributions worth stating. This is
   `EVIDENCE`-scoped calibration in the same family as 001 leads 2–3, and it
   is the only part of this census with dynamical content.
3. **Retire the flow-cutter tooling plan from the queue.** No treewidth
   computation on any truncated Collatz digraph can return anything but 2,
   so the tool never needs to be built. Recorded so a later cycle does not
   re-derive the lead from the same 001 report.

## References

- Wirsching, *The Dynamical System Generated by the 3n+1 Function* — standard
  development of the inverted Collatz graph as a rooted tree. [T]
- Ebert (2021), *A Graph Theoretical Approach to the Collatz Problem*,
  arXiv:1905.07575 — inverted-perspective rooted tree. [T]
- Lengauer–Tarjan (1979); Cooper–Harvey–Kennedy (2001) — dominator tree
  algorithms; the iterative form is what `dominator_check` implements. [T]
- Prior attempt: `problems/collatz/attempts/001-ideation-sweep-fields.md`
  (lead 1, and gap 2's novelty check, both discharged here).
- Protocol: `docs/SWARM.md` (cross-family skeptic rule); brief and worker
  stances committed at `explore/graph-census-skeptic-2026-08/`; raw worker
  returns in `$MATHLAB_OUT/swarm/collatz-graph-census-skeptic/`
  (gitignored working space, per-job meta carries model and usage).
