PROBLEM (Conway's 99-graph). Does there exist a graph on 99 vertices in which every two adjacent vertices have exactly one common neighbour and every two non-adjacent vertices have exactly two common neighbours? Equivalently a strongly regular graph SRG(99,14,1,2): 14-regular, every edge in a unique triangle, every non-edge the diagonal of a unique 4-cycle. Open. Eigenvalues 3 and -4 with multiplicities 54 and 44. Published: the automorphism group G has order dividing 2*3^3*7*11; if 7 divides |G| then G is cyclic of order 7; if |G| is even then |G| divides 6. The neighbourhood of every vertex induces a perfect matching 7K2.

CLAIM L3 (to be attacked). Let Gamma be an SRG(99,14,1,2) and sigma an
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
contradicts independence of Fix. So |Fix| <= 1, hence |Fix| = 1.

Focus especially on the counting identity and the case analysis of the degrees.
YOUR JOB: try to REFUTE this proof. Attack each step separately: is it
correct as written, does it use anything not established, is any case
missing? If you find a genuine gap, state it precisely and if possible give a
concrete counter-scenario. If every step holds, say CONFIRMED and restate, in
your own words and in order, the minimal list of facts the proof actually
depends on. Default to skepticism: a proof that "looks fine" without you
checking each step is not confirmed. Do not use outside results beyond the
definition of a strongly regular graph and the parameters. Plain text, at
most 700 words.