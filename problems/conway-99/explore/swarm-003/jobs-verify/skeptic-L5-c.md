PROBLEM (Conway's 99-graph). Does there exist a graph on 99 vertices in which every two adjacent vertices have exactly one common neighbour and every two non-adjacent vertices have exactly two common neighbours? Equivalently a strongly regular graph SRG(99,14,1,2): 14-regular, every edge in a unique triangle, every non-edge the diagonal of a unique 4-cycle. Open. Eigenvalues 3 and -4 with multiplicities 54 and 44. Published: the automorphism group G has order dividing 2*3^3*7*11; if 7 divides |G| then G is cyclic of order 7; if |G| is even then |G| divides 6. The neighbourhood of every vertex induces a perfect matching 7K2.

CLAIM L5 (to be attacked). Let Gamma be an SRG(99,14,1,2) and tau an
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
of the triangle on an (8,2) edge would need both degrees). Hence f in {0,3}.

Focus especially on whether the argument silently assumes that a fixed vertex v has its fixed neighbours forming matching edges, and on the SRG(33,8,1,2) exclusion: verify the multiplicity computation yourself.
YOUR JOB: try to REFUTE this proof. Attack each step separately: is it
correct as written, does it use anything not established, is any case
missing? If you find a genuine gap, state it precisely and if possible give a
concrete counter-scenario. If every step holds, say CONFIRMED and restate, in
your own words and in order, the minimal list of facts the proof actually
depends on. Default to skepticism: a proof that "looks fine" without you
checking each step is not confirmed. Do not use outside results beyond the
definition of a strongly regular graph and the parameters. Plain text, at
most 700 words.