PROBLEM (Conway's 99-graph). Does there exist a graph on 99 vertices in which every two adjacent vertices have exactly one common neighbour and every two non-adjacent vertices have exactly two common neighbours? Equivalently a strongly regular graph SRG(99,14,1,2): 14-regular, every edge in a unique triangle, every non-edge the diagonal of a unique 4-cycle. Open. Eigenvalues 3 and -4 with multiplicities 54 and 44. Published: the automorphism group G has order dividing 2*3^3*7*11; if 7 divides |G| then G is cyclic of order 7; if |G| is even then |G| divides 6. The neighbourhood of every vertex induces a perfect matching 7K2.

THE ORDER-7 SUB-CASE (facts established earlier, take them as given). Let tau be an automorphism of order 7 of such a graph. Then tau has exactly one fixed vertex v0 and fourteen orbits O1..O14 of size 7. N(v0) is the union of two orbits, say O1 and O2, and the perfect matching on N(v0) pairs each vertex of O1 with one of O2 (so no matching edge lies inside O1 or inside O2). Write r_ij for the number of neighbours a vertex of O_i has in O_j (constant on O_i). With O0 = {v0}: r_00 = 0, r_01 = r_02 = 7, r_0j = 0 for j >= 3, r_10 = r_20 = 1, r_i0 = 0 for i >= 3, r_11 = r_22 = 0, r_12 = r_21 = 1. The 15x15 matrix R = (r_ij) is the ORBIT MATRIX. It satisfies, with n_i = |O_i|, k=14, lam=1, mu=2:
  (S)  r_ij * n_i = r_ji * n_j
  (R)  sum_j r_ij = 14
  (Q)  sum_l r_il * r_lj = 14*[i=j] + lam*r_ij + mu*(n_j - [i=j] - r_ij)
Among the 7-orbits R is symmetric and every diagonal entry r_ii is even (the connection set of a vertex within its own orbit under a cyclic group of odd order is closed under negation and excludes 0). Every entry is between 0 and 7. Each orbit matrix must then be LIFTED: r_ij becomes a connection set S_ij, a subset of Z_7 of size r_ij, with S_ji = -S_ij and S_ii = -S_ii, 0 not in S_ii; vertices are (i, x) with x in Z_7, and (i,x) ~ (j,y) iff y - x is in S_ij.

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

DESIGN EMPHASIS FOR THIS DRAFT: Decide the connection sets pair by pair in an order that makes as many group-ring equations checkable as early as possible, and check them incrementally.

Reply with the complete program in a single ```python code block and NOTHING
else: no prose before or after it. Put any explanation in code comments. Keep
the program under 250 lines. The reply is cut off at 12000 tokens, so be
economical.
