PROBLEM (Conway's 99-graph). Does there exist a graph on 99 vertices in which every two adjacent vertices have exactly one common neighbour and every two non-adjacent vertices have exactly two common neighbours? Equivalently a strongly regular graph SRG(99,14,1,2): 14-regular, every edge in a unique triangle, every non-edge the diagonal of a unique 4-cycle. Open. Eigenvalues 3 and -4 with multiplicities 54 and 44. Published: the automorphism group G has order dividing 2*3^3*7*11; if 7 divides |G| then G is cyclic of order 7; if |G| is even then |G| divides 6. The neighbourhood of every vertex induces a perfect matching 7K2.

THE ORDER-7 SUB-CASE (facts established earlier, take them as given). Let tau be an automorphism of order 7 of such a graph. Then tau has exactly one fixed vertex v0 and fourteen orbits O1..O14 of size 7. N(v0) is the union of two orbits, say O1 and O2, and the perfect matching on N(v0) pairs each vertex of O1 with one of O2 (so no matching edge lies inside O1 or inside O2). Write r_ij for the number of neighbours a vertex of O_i has in O_j (constant on O_i). With O0 = {v0}: r_00 = 0, r_01 = r_02 = 7, r_0j = 0 for j >= 3, r_10 = r_20 = 1, r_i0 = 0 for i >= 3, r_11 = r_22 = 0, r_12 = r_21 = 1. The 15x15 matrix R = (r_ij) is the ORBIT MATRIX. It satisfies, with n_i = |O_i|, k=14, lam=1, mu=2:
  (S)  r_ij * n_i = r_ji * n_j
  (R)  sum_j r_ij = 14
  (Q)  sum_l r_il * r_lj = 14*[i=j] + lam*r_ij + mu*(n_j - [i=j] - r_ij)
Among the 7-orbits R is symmetric and every diagonal entry r_ii is even (the connection set of a vertex within its own orbit under a cyclic group of odd order is closed under negation and excludes 0). Every entry is between 0 and 7. Each orbit matrix must then be LIFTED: r_ij becomes a connection set S_ij, a subset of Z_7 of size r_ij, with S_ji = -S_ij and S_ii = -S_ii, 0 not in S_ii; vertices are (i, x) with x in Z_7, and (i,x) ~ (j,y) iff y - x is in S_ij.

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

YOUR LENS: Finite geometry: interpret 99 = 1 + 14 + 84 and the 7K2 neighbourhoods as points and lines of a known geometry (generalized quadrangles, partial geometries, polar spaces) and ask what the order-7 symmetry would have to be there.

TASK. From inside that lens only, propose 1 to 3 attack routes on T1, T2 or
T3. For each route give: (a) the core object or quantity your field would
compute, (b) the first concrete calculation or search, small enough to run in
an afternoon in Python, (c) the outcome that would KILL the route. A route
without a falsifiable first step does not count. If your field has no
purchase here, say so and say precisely why -- that is a valid answer.
RULES FOR YOUR ANSWER. Be concrete and honest. Label every unproven step SPECULATION. Do not claim to have solved the problem. If you have no purchase, say so and say why -- that is a useful answer. Plain text, no more than about 700 words.