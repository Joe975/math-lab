"""Generate the swarm briefs for attempt 003 (order-7 case, local Qwen workers).

    python problems/conway-99/explore/swarm-003/make_briefs.py

Writes jobs-verify/ (enumerator and lifter drafts, skeptic passes, Z_33
re-derivations) and the ideation template + lens values. The prompts are the
provenance for the attempt record; the swarm output lives in $MATHLAB_OUT.
Tier 1: the briefs quote this lab's lemmas, so anything built on them is
`informed`.
"""

from pathlib import Path

base = Path(__file__).resolve().parent
jv = base / "jobs-verify"
jv.mkdir(parents=True, exist_ok=True)

STATEMENT = (
    "PROBLEM (Conway's 99-graph). Does there exist a graph on 99 vertices in "
    "which every two adjacent vertices have exactly one common neighbour and "
    "every two non-adjacent vertices have exactly two common neighbours? "
    "Equivalently a strongly regular graph SRG(99,14,1,2): 14-regular, every "
    "edge in a unique triangle, every non-edge the diagonal of a unique 4-cycle. "
    "Open. Eigenvalues 3 and -4 with multiplicities 54 and 44. Published: the "
    "automorphism group G has order dividing 2*3^3*7*11; if 7 divides |G| then "
    "G is cyclic of order 7; if |G| is even then |G| divides 6. The "
    "neighbourhood of every vertex induces a perfect matching 7K2."
)

ORDER7 = (
    "THE ORDER-7 SUB-CASE (facts established earlier, take them as given). Let "
    "tau be an automorphism of order 7 of such a graph. Then tau has exactly one "
    "fixed vertex v0 and fourteen orbits O1..O14 of size 7. N(v0) is the union "
    "of two orbits, say O1 and O2, and the perfect matching on N(v0) pairs each "
    "vertex of O1 with one of O2 (so no matching edge lies inside O1 or inside "
    "O2). Write r_ij for the number of neighbours a vertex of O_i has in O_j "
    "(constant on O_i). With O0 = {v0}: r_00 = 0, r_01 = r_02 = 7, r_0j = 0 for "
    "j >= 3, r_10 = r_20 = 1, r_i0 = 0 for i >= 3, r_11 = r_22 = 0, "
    "r_12 = r_21 = 1. The 15x15 matrix R = (r_ij) is the ORBIT MATRIX. It "
    "satisfies, with n_i = |O_i|, k=14, lam=1, mu=2:\n"
    "  (S)  r_ij * n_i = r_ji * n_j\n"
    "  (R)  sum_j r_ij = 14\n"
    "  (Q)  sum_l r_il * r_lj = 14*[i=j] + lam*r_ij + mu*(n_j - [i=j] - r_ij)\n"
    "Among the 7-orbits R is symmetric and every diagonal entry r_ii is even "
    "(the connection set of a vertex within its own orbit under a cyclic group "
    "of odd order is closed under negation and excludes 0). Every entry is "
    "between 0 and 7. Each orbit matrix must then be LIFTED: r_ij becomes a "
    "connection set S_ij, a subset of Z_7 of size r_ij, with S_ji = -S_ij and "
    "S_ii = -S_ii, 0 not in S_ii; vertices are (i, x) with x in Z_7, and "
    "(i,x) ~ (j,y) iff y - x is in S_ij."
)

RULES = (
    "RULES FOR YOUR ANSWER. Be concrete and honest. Label every unproven step "
    "SPECULATION. Do not claim to have solved the problem. If you have no "
    "purchase, say so and say why -- that is a useful answer. Plain text, no "
    "more than about 700 words."
)

# --- enumerator drafts -----------------------------------------------------

ENUM_SPEC = STATEMENT + "\n\n" + ORDER7 + """

TASK. Write ONE standalone Python 3 program (standard library only, no numpy)
that enumerates ALL orbit matrices for a given orbit-size profile, with the
(S), (R), (Q) conditions above, the parity rule for odd-size orbits under a
cyclic group, and an optional dictionary of forced entries (a "seed"). It must
be a complete, runnable file.

Interface, exactly:
  enumerate_orbit_matrices(sizes, k, lam, mu, seed=None, odd_diag_even=True,
                           limit=None, node_budget=None) -> (matrices, stats)
    sizes: list of orbit sizes (ints). seed: dict {(i,j): value}. limit: stop
    after this many solutions. node_budget: stop after this many search nodes
    (return what was found so far and set stats["exhausted"]=False).
    Returns a list of t x t integer matrices (lists of lists) and a dict stats
    with at least {"nodes": int, "exhausted": bool}.
  check_orbit_matrix(R, sizes, k, lam, mu) -> (bool, reason): an INDEPENDENT
    direct check of (S), (R), (Q) by plain loops, not reusing the search code.
  A CLI: python prog.py --sizes 3,3,3 --k 4 --lam 1 --mu 2
         [--seed "0,0=0;0,1=3"] [--limit N] [--node-budget N]
    printing the count, the node count, whether the search was exhaustive, and
    each matrix as JSON on its own line.

Correctness requirements: the enumeration must be EXHAUSTIVE (no valid matrix
may be missed) -- prune only with conditions implied by (S),(R),(Q) and the
parity rule. Do NOT apply any symmetry breaking or isomorph rejection: output
every solution, including relabellings. Every returned matrix must pass
check_orbit_matrix (assert this before returning).

Known answers you can test against: sizes [3,3,3] with k=4,lam=1,mu=2 (the
Paley graph on 9 vertices under a Z_3 action) has exactly 3 solutions, one of
them [[0,2,2],[2,0,2],[2,2,0]]. sizes [33,33,33] with k=14,lam=1,mu=2 has
exactly 1 solution [[2,6,6],[6,2,6],[6,6,2]].

DESIGN EMPHASIS FOR THIS DRAFT: {emphasis}

Reply with the complete program in a single ```python code block and NOTHING
else: no prose before or after it. Put any explanation in code comments. Keep
the program under 250 lines. The reply is cut off at 12000 tokens, so be
economical.
"""

EMPHASES = [
    "Fill the matrix row by row. Because the matrix restricted to equal-size "
    "orbits is symmetric, once rows a and b are both complete the (Q) equation "
    "for the pair (a,b) is fully determined, and while filling row i the inner "
    "products with every earlier row are fixed constants -- prune each partial "
    "row against all of those linear constraints (using min/max bounds on the "
    "remaining entries) and against the row's fixed sum of squares from (Q) at "
    "(i,i).",
    "Precompute, for each row, the set of feasible row multisets from (R) and "
    "the diagonal (Q) equation (sum and sum of squares fixed), and enumerate "
    "rows from those multisets; check pairwise (Q) equations as soon as both "
    "rows are complete.",
    "Plain recursive depth-first search over the upper triangle in row-major "
    "order with forward checking: after each assignment, verify every (Q) "
    "equation all of whose terms are already known, and bound-check partially "
    "known ones.",
    "Treat it as a constraint-satisfaction problem: keep for every unassigned "
    "entry a domain (a set of allowed values), propagate (S), (R) and the "
    "bounds of every (Q) equation to shrink domains after each assignment, and "
    "branch on the entry with the smallest domain.",
    "Column-oriented: assign one full column of the upper triangle at a time "
    "(entries r_0j..r_jj), which completes column j and row j simultaneously "
    "by symmetry, then check every (Q) equation between row j and earlier rows "
    "exactly.",
]

# --- lifter drafts ---------------------------------------------------------

LIFT_SPEC = STATEMENT + "\n\n" + ORDER7 + """

TASK. Write ONE standalone Python 3 program (standard library only) that tries
to LIFT an orbit matrix to an actual graph for a cyclic group Z_m acting with
orbits of equal size m plus optionally one fixed vertex. Input: m, the orbit
matrix R (t x t; orbit 0 may be the fixed point with size 1, all others size
m), and k, lam, mu. Search over connection sets S_ij (subset of Z_m of size
r_ij, S_ji = -S_ij, S_ii symmetric without 0; the fixed point is adjacent to
all of an orbit or none of it) and output every graph found as a list of
adjacency bitmasks, checking each candidate with an INDEPENDENT strongly
regular checker (count common neighbours for every pair of vertices
directly). The search must be exhaustive up to an optional node budget, and
must prune using the group-ring form of the SRG equation: for each pair of
orbits (i,j) and each difference d in Z_m, the number of common neighbours of
(i,0) and (j,d) must equal lam if d in S_ij else mu (and k when i=j and d=0),
computable as soon as the sets involved are decided.

Interface: lift(m, R, k, lam, mu, node_budget=None) -> (graphs, stats) with
stats containing "nodes" and "exhausted".
CLI: python prog.py --m 3 --R "[[0,2,2],[2,0,2],[2,2,0]]" --k 4 --lam 1 --mu 2
Known answer: that call must find graphs isomorphic to the Paley graph on 9
vertices (SRG(9,4,1,2)); at least one graph must be found.

DESIGN EMPHASIS FOR THIS DRAFT: {emphasis}

Reply with the complete program in a single ```python code block and NOTHING
else: no prose before or after it. Put any explanation in code comments. Keep
the program under 250 lines. The reply is cut off at 12000 tokens, so be
economical.
"""

LIFT_EMPHASES = [
    "Decide the connection sets pair by pair in an order that makes as many "
    "group-ring equations checkable as early as possible, and check them "
    "incrementally.",
    "Use bitmasks over Z_m for the sets and precomputed difference tables; the "
    "core must be exact integer counting.",
]

# --- skeptic briefs --------------------------------------------------------

L5 = """CLAIM L5 (to be attacked). Let Gamma be an SRG(99,14,1,2) and tau an
automorphism of order 3. Then tau has exactly 0 or 3 fixed points.
PROOF OFFERED. Let F = Fix(tau), f = |F|, and suppose f > 0. For v in F, tau
permutes the 7 matching edges of N(v) in orbits of size 1 or 3, so the number
of fixed edges is congruent to 7 mod 3, i.e. 1, 4 or 7. A tau-fixed edge has
both endpoints fixed (tau has odd order, so it cannot swap the two endpoints).
Seven fixed edges would fix N(v) pointwise, and an automorphism fixing a
vertex and its whole neighbourhood pointwise is the identity (every other
vertex u is determined by its unique pair of neighbours in N(v): u is
non-adjacent to v so has exactly mu=2 common neighbours with v). Hence
|F intersect N(v)| is 2 or 8 for every v in F. F is closed under taking common
neighbours of two of its vertices: if x,y in F then tau permutes the set of
common neighbours of x and y, which has size at most 2, and an odd-order
permutation of a set of size at most 2 is trivial. Now for any induced
subgraph F closed in that sense, count paths u-x-w with u in F fixed, x in
N_F(u), w in F, w != u: sum over x in N_F(u) of deg_F(x) =
deg_F(u)*(1+lam-mu) + mu*(f-1) = 2f-2 for (lam,mu)=(1,2). If deg_F(u)=2 for
all u then 2*2 = 2f-2 gives f=3 (a triangle). If deg_F(u)=8 for all u then
f=33 and F is an SRG(33,8,1,2), which fails the integrality of eigenvalue
multiplicities. Mixed degrees force f=9 with every degree-8 vertex adjacent
only to degree-2 vertices and vice versa, which lam=1 kills (the third vertex
of the triangle on an (8,2) edge would need both degrees). Hence f in {0,3}."""

L3 = """CLAIM L3 (to be attacked). Let Gamma be an SRG(99,14,1,2) and sigma an
automorphism of order 7. Then sigma has exactly one fixed vertex.
PROOF OFFERED. |Fix(sigma)| is congruent to 99 mod 7, i.e. to 1, so Fix is
non-empty. For v in Fix(sigma), sigma permutes the 7 edges of the perfect
matching on N(v); orbits have size 1 or 7. If all 7 edges were fixed, sigma
would fix N(v) pointwise (odd order cannot swap an edge's endpoints), and an
automorphism fixing v and N(v) pointwise is the identity (each vertex u not
adjacent to v is determined by its unique pair of common neighbours with v,
and each vertex of N(v) is fixed already). So sigma acts on the 7 edges as a
single 7-cycle and fixes no vertex of N(v). Hence Fix(sigma) is an independent
set. Fix is closed under common neighbours (sigma permutes the
common-neighbour set of two fixed vertices, size at most 2, and an odd-order
permutation of at most 2 points is trivial). Two distinct fixed vertices are
non-adjacent, so have exactly 2 common neighbours, both fixed -- but that
contradicts independence of Fix. So |Fix| <= 1, hence |Fix| = 1."""

R12 = """CLAIM (to be attacked). Let sigma be an automorphism of order 7 of an
SRG(99,14,1,2) with a unique fixed vertex v0 (take that as given). Then N(v0)
is a union of exactly two sigma-orbits O1, O2 of size 7, and the perfect
matching on N(v0) has every edge going between O1 and O2 (none inside O1 or
O2). Consequently in the orbit matrix r_11 = r_22 = 0 and r_12 = r_21 = 1,
and r_01 = r_02 = 7.
PROOF OFFERED. N(v0) is sigma-invariant of size 14 and contains no fixed
vertex (sigma fixes only v0), so it is a union of orbits of size 7: exactly
two. If a matching edge lay inside O1, its 7 images under sigma would be 7
distinct matching edges inside O1, giving 14 endpoint incidences on the 7
vertices of O1, but each vertex is in exactly one matching edge, so at most 7
incidences: contradiction. So every matching edge joins O1 to O2, each vertex
of O1 has exactly one neighbour in O2 (its partner) and none in O1 (N(v0)
induces exactly the matching)."""

SKEPTIC_TAIL = """
YOUR JOB: try to REFUTE this proof. Attack each step separately: is it
correct as written, does it use anything not established, is any case
missing? If you find a genuine gap, state it precisely and if possible give a
concrete counter-scenario. If every step holds, say CONFIRMED and restate, in
your own words and in order, the minimal list of facts the proof actually
depends on. Default to skepticism: a proof that "looks fine" without you
checking each step is not confirmed. Do not use outside results beyond the
definition of a strongly regular graph and the parameters. Plain text, at
most 700 words."""

FOCUS = {
    "a": "Focus especially on the counting identity and the case analysis of "
         "the degrees.",
    "b": "Focus especially on the closure claim (fixed set closed under common "
         "neighbours) and on the rigidity claim (fixing a vertex and its "
         "neighbourhood pointwise forces the identity).",
    "c": "Focus especially on whether the argument silently assumes that a "
         "fixed vertex v has its fixed neighbours forming matching edges, and "
         "on the SRG(33,8,1,2) exclusion: verify the multiplicity computation "
         "yourself.",
    "r12": "Focus on whether the orbit structure of N(v0) and the "
           "endpoint-counting really force r_12 = 1 and r_11 = r_22 = 0.",
}

# --- Z_33 re-derivations ---------------------------------------------------

Z33 = STATEMENT + """

CLAIM TO RE-DERIVE INDEPENDENTLY. Suppose an SRG(99,14,1,2) has an
automorphism of order 33 acting semiregularly, i.e. with three orbits of size
33. Its orbit matrix R (r_ij = number of neighbours in O_j of a vertex of O_i)
is a symmetric 3x3 matrix of integers 0..33 with even diagonal, rows summing
to 14, satisfying sum_l r_il r_lj = 14[i=j] + r_ij + 2(33 - [i=j] - r_ij) for
all i,j. Claim: there is EXACTLY ONE such matrix, namely
[[2,6,6],[6,2,6],[6,6,2]].

TASK. {task} Reply with the complete program in a single ```python code
block, followed by the exact output you expect it to print.
"""

Z33_TASKS = {
    "a": "Write a short standalone Python 3 program that brute-forces every "
         "symmetric 3x3 integer matrix with entries 0..33, even diagonal and "
         "row sums 14, tests the equations by direct evaluation, and prints "
         "every solution and the total count.",
    "b": "First derive by hand what the equations force (use the row sum to "
         "eliminate variables and reduce to a small system), then write a "
         "short standalone Python 3 program that enumerates ONLY over the "
         "diagonal entries, solves for the rest, and prints every solution and "
         "the total count. Also state whether your hand derivation agrees with "
         "the claim.",
}

# --- ideation template -----------------------------------------------------

IDEATE = STATEMENT + "\n\n" + ORDER7 + """

WHAT THIS LAB HAS ALREADY TRIED (one line each; do not repeat these):
- Exact feasibility arithmetic, a forced "pair model" of the graph around one
  vertex, and a partner-regular case split: reduced but not closed.
- Constraint-propagation extension search from one vertex: fails its own
  positive control (cannot rebuild the known SRG(243,22,1,2)), so says
  nothing about 99.
- Orbit-matrix enumeration: complete for a semiregular order-33 automorphism
  (exactly one orbit matrix, which does not lift -- the three diagonal
  equations force all off-diagonal autocorrelations of the connection sets,
  and no 6-subset of Z_33 has the required autocorrelation), but the
  backtracking does not reach a first solution for 9 or more orbits (Z_11
  with 9 orbits, Z_9 with 11, the order-7 profile with 15, Z_3 with 33).
- Lemmas: an order-3 automorphism has 0 or 3 fixed points; no automorphism of
  order 27; order-9 automorphisms are semiregular; an order-7 automorphism has
  exactly one fixed point; an order-5 automorphism would have fixed set the
  3x3 rook graph with every moved vertex having exactly one fixed neighbour
  (no contradiction found there).
- Simulated annealing over regular graphs, over the pair model, over
  Z_m-invariant graphs (connection sets), and over unions of pair-orbits of
  an assumed group: every engine tops out below 45 vertices on known SRGs, so
  their failure at 99 is uninformative.

THREE OPEN TARGETS. (T1) The order-7 case: enumerate the 15x15 orbit matrices
and lift them, or prove no lift exists, or prove no order-7 automorphism
exists. (T2) Find the successor of the forced-autocorrelation trick for more
than 3 orbits: with t orbits there are t diagonal equations but t(t-1)/2
unknown off-diagonal autocorrelation functions; the t(t-1)/2 off-diagonal
group-ring equations are so far unused. (T3) Involutions: even order forces
|Aut| to divide 6; find a fixed-point lemma for an involution (the odd-order
trick "a permutation of a set of size <= 2 is trivial" fails).

YOUR LENS: {{value}}

TASK. From inside that lens only, propose 1 to 3 attack routes on T1, T2 or
T3. For each route give: (a) the core object or quantity your field would
compute, (b) the first concrete calculation or search, small enough to run in
an afternoon in Python, (c) the outcome that would KILL the route. A route
without a falsifiable first step does not count. If your field has no
purchase here, say so and say precisely why -- that is a valid answer.
""" + RULES

LENSES = [
    "Group rings and character theory: view the connection sets S_ij as "
    "elements of Z[Z_7] and the SRG equation as a matrix equation over the "
    "group ring; use the 7 characters of Z_7 (the Fourier transform over Z_7) "
    "to diagonalise.",
    "Linear algebra and rank: set up the full system of diagonal and "
    "off-diagonal equations in the unknown autocorrelation functions and "
    "determine its rank; what is forced and what is free.",
    "Spectral graph theory and eigenvalue interlacing: the quotient (orbit) "
    "matrix's eigenvalues interlace those of the graph (3, -4, 14); use this "
    "to prune orbit matrices.",
    "Coding theory: linear codes spanned by the adjacency matrix mod 2 or "
    "mod 7; dimensions and weight constraints.",
    "Design theory and partial geometries: the unique triangle per edge and "
    "unique 4-cycle per non-edge as a partial linear space; the 33 lines "
    "(triangles) as a configuration.",
    "SAT and symmetry breaking: encode the order-7 lifted graph (14 orbits "
    "over Z_7 plus a fixed point) as a SAT instance with the symmetry already "
    "quotiented out; estimate the variable count and what symmetry breaking "
    "remains.",
    "Difference sets and cyclotomy in Z_7: which subsets of Z_7 have which "
    "autocorrelation vectors; the multiplier group of Z_7 (units 1..6) acting "
    "on connection sets.",
    "Association schemes and coherent configurations: refine the orbit "
    "partition under Z_7 into a coherent configuration and apply the "
    "intersection-number constraints beyond the orbit-matrix level.",
    "Extremal and probabilistic combinatorics: counting arguments on the "
    "number of triangles and 4-cycles through orbits; second-moment estimates "
    "for whether a random lift could exist.",
    "Voltage graphs and covering-graph theory: the graph is a Z_7-voltage "
    "graph over the 15-vertex quotient; constrain lifts.",
    "Semidefinite programming and the Lovasz theta function: SDP relaxations "
    "of the existence question restricted to Z_7-invariant matrices.",
    "Finite geometry: interpret 99 = 1 + 14 + 84 and the 7K2 neighbourhoods "
    "as points and lines of a known geometry (generalized quadrangles, "
    "partial geometries, polar spaces) and ask what the order-7 symmetry "
    "would have to be there.",
    "Integer programming: formulate the orbit-matrix enumeration as an "
    "integer quadratic program and identify valid inequalities (cuts) that "
    "make it small.",
    "Fixed-point theory for group actions on graphs (Burnside-type counting): "
    "count fixed edges, fixed triangles and fixed 4-cycles of an involution or "
    "an order-3 element to derive fixed-point lemmas (target T3).",
    "Number theory: divisibility and congruence conditions (mod 7, mod 3, "
    "mod 2) on entries of the orbit matrix and on the counts of orbits of "
    "triangles and 4-cycles under the group.",
    "Local search and heuristics literature for strongly regular graphs "
    "(tabu, annealing, genetic): what move sets and energy functions have "
    "succeeded for SRGs near 100 vertices, and how to exploit Z_7 invariance.",
    "Polynomial methods: encode the connection sets as 0/1 polynomials in "
    "Z[x]/(x^7-1) and use Groebner bases or resultants over the small system.",
    "Structural graph theory: forced local structure two steps from the fixed "
    "vertex v0 (the 84 vertices at distance 2 and their 7-orbits) and what the "
    "unique-4-cycle condition forces about how orbits O3..O14 attach to O1 "
    "and O2.",
]


def main():
    # Drafts 1-3 and lifter 1 run; the rest are kept as briefs for a later
    # sweep -- at ~6 tokens/s per slot the local model affords about a dozen
    # 12k-token jobs in a session.
    extra = base / "jobs-verify-extra"
    extra.mkdir(exist_ok=True)
    for n, e in enumerate(EMPHASES, 1):
        d = jv if n <= 3 else extra
        (d / f"enum-{n}.md").write_text(ENUM_SPEC.replace("{emphasis}", e), encoding="utf-8")
    for n, e in enumerate(LIFT_EMPHASES, 1):
        d = jv if n <= 1 else extra
        (d / f"lift-{n}.md").write_text(LIFT_SPEC.replace("{emphasis}", e), encoding="utf-8")
    for name, claim in [("skeptic-L5-a", L5), ("skeptic-L5-b", L5), ("skeptic-L5-c", L5),
                        ("skeptic-L3-a", L3), ("skeptic-L3-b", L3), ("skeptic-r12", R12)]:
        focus = FOCUS[name.split("-")[-1]]
        (jv / f"{name}.md").write_text(
            STATEMENT + "\n\n" + claim + "\n\n" + focus + SKEPTIC_TAIL, encoding="utf-8")
    for key, task in Z33_TASKS.items():
        (jv / f"z33-{key}.md").write_text(Z33.replace("{task}", task), encoding="utf-8")
    (base / "ideate-template.md").write_text(IDEATE, encoding="utf-8")
    (base / "ideate-lenses.txt").write_text("\n".join(LENSES) + "\n", encoding="utf-8")
    print(f"wrote {len(list(jv.iterdir()))} verify briefs, {len(LENSES)} lenses")


if __name__ == "__main__":
    main()
