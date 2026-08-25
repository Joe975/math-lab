PROBLEM (Conway's 99-graph). Does there exist a graph on 99 vertices in which every two adjacent vertices have exactly one common neighbour and every two non-adjacent vertices have exactly two common neighbours? Equivalently a strongly regular graph SRG(99,14,1,2): 14-regular, every edge in a unique triangle, every non-edge the diagonal of a unique 4-cycle. Open. Eigenvalues 3 and -4 with multiplicities 54 and 44. Published: the automorphism group G has order dividing 2*3^3*7*11; if 7 divides |G| then G is cyclic of order 7; if |G| is even then |G| divides 6. The neighbourhood of every vertex induces a perfect matching 7K2.

THE ORDER-7 SUB-CASE (facts established earlier, take them as given). Let tau be an automorphism of order 7 of such a graph. Then tau has exactly one fixed vertex v0 and fourteen orbits O1..O14 of size 7. N(v0) is the union of two orbits, say O1 and O2, and the perfect matching on N(v0) pairs each vertex of O1 with one of O2 (so no matching edge lies inside O1 or inside O2). Write r_ij for the number of neighbours a vertex of O_i has in O_j (constant on O_i). With O0 = {v0}: r_00 = 0, r_01 = r_02 = 7, r_0j = 0 for j >= 3, r_10 = r_20 = 1, r_i0 = 0 for i >= 3, r_11 = r_22 = 0, r_12 = r_21 = 1. The 15x15 matrix R = (r_ij) is the ORBIT MATRIX. It satisfies, with n_i = |O_i|, k=14, lam=1, mu=2:
  (S)  r_ij * n_i = r_ji * n_j
  (R)  sum_j r_ij = 14
  (Q)  sum_l r_il * r_lj = 14*[i=j] + lam*r_ij + mu*(n_j - [i=j] - r_ij)
Among the 7-orbits R is symmetric and every diagonal entry r_ii is even (the connection set of a vertex within its own orbit under a cyclic group of odd order is closed under negation and excludes 0). Every entry is between 0 and 7. Each orbit matrix must then be LIFTED: r_ij becomes a connection set S_ij, a subset of Z_7 of size r_ij, with S_ji = -S_ij and S_ii = -S_ii, 0 not in S_ii; vertices are (i, x) with x in Z_7, and (i,x) ~ (j,y) iff y - x is in S_ij.

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

DESIGN EMPHASIS FOR THIS DRAFT: Treat it as a constraint-satisfaction problem: keep for every unassigned entry a domain (a set of allowed values), propagate (S), (R) and the bounds of every (Q) equation to shrink domains after each assignment, and branch on the entry with the smallest domain.

Reply with the complete program in a single ```python code block and NOTHING
else: no prose before or after it. Put any explanation in code comments. Keep
the program under 250 lines. The reply is cut off at 12000 tokens, so be
economical.
