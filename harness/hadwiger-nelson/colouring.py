"""Exact chromatic-number tools, with three mutually independent methods.

The three are deliberately different algorithms, not three call sites of one
routine, so that agreement between them is evidence:

1. ``k_colouring``      -- DSATUR-ordered backtracking with symmetry breaking.
                           Fast; returns a colouring or None.
2. ``brute_force_k``    -- plain enumeration of k^n assignments, no pruning
                           beyond an early conflict check.  Only for small n,
                           but it shares no logic with (1).
3. ``chromatic_polynomial`` -- deletion/contraction.  Computes the whole
                           polynomial P(G,k); the chromatic number is the least
                           k with P(G,k) > 0.  Structurally unrelated to both.

``verify_colouring`` is the linear checker: it is what should be run on a
colouring produced elsewhere.
"""

from __future__ import annotations

from itertools import product


def adjacency(n, edges):
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    return adj


def verify_colouring(n, edges, colouring):
    """Proper-colouring check.  Deliberately the dumbest possible loop."""
    if len(colouring) != n:
        return False
    for a, b in edges:
        if colouring[a] == colouring[b]:
            return False
    return True


# --- method 1: DSATUR backtracking ------------------------------------------


def k_colouring(n, edges, k):
    """Return a proper k-colouring as a list, or None if none exists.

    Complete search: a None return is a proof of non-k-colourability modulo the
    correctness of this routine, which is why method 2 exists.
    """
    if n == 0:
        return []
    adj = adjacency(n, edges)
    colour = [-1] * n
    # Symmetry breaking: colour a maximal greedy clique with distinct colours.
    clique = _greedy_clique(n, adj)
    if len(clique) > k:
        return None
    for idx, v in enumerate(clique):
        colour[v] = idx
    forbidden = [set() for _ in range(n)]
    for v in range(n):
        if colour[v] >= 0:
            for w in adj[v]:
                forbidden[w].add(colour[v])
    for v in clique:
        if colour[v] in forbidden[v]:
            return None

    def pick():
        best, best_key = -1, None
        for v in range(n):
            if colour[v] >= 0:
                continue
            key = (-len(forbidden[v]), -len(adj[v]))
            if best_key is None or key < best_key:
                best, best_key = v, key
        return best

    def rec(assigned):
        if assigned == n:
            return True
        v = pick()
        avail = [c for c in range(k) if c not in forbidden[v]]
        if not avail:
            return False
        # Only ever open one fresh colour: colours above max-used+1 are symmetric.
        used = max((c for c in colour if c >= 0), default=-1)
        avail = [c for c in avail if c <= used + 1]
        for c in avail:
            colour[v] = c
            touched = []
            for w in adj[v]:
                if colour[w] < 0 and c not in forbidden[w]:
                    forbidden[w].add(c)
                    touched.append(w)
            if rec(assigned + 1):
                return True
            for w in touched:
                forbidden[w].discard(c)
            colour[v] = -1
        return False

    if rec(sum(1 for c in colour if c >= 0)):
        return list(colour)
    return None


def _greedy_clique(n, adj):
    order = sorted(range(n), key=lambda v: -len(adj[v]))
    clique = []
    for v in order:
        if all(v in adj[u] for u in clique):
            clique.append(v)
    return clique


def chromatic_number(n, edges, lo=1, hi=None):
    if n == 0:
        return 0
    if hi is None:
        hi = n
    for k in range(lo, hi + 1):
        if k_colouring(n, edges, k) is not None:
            return k
    return None


# --- method 2: brute force ---------------------------------------------------


def brute_force_k(n, edges, k):
    """Enumerate every assignment in {0..k-1}^n.  Returns a colouring or None.

    Fixes vertex 0 to colour 0 (harmless: colours are interchangeable).  No
    other pruning, no ordering heuristic, no shared code with ``k_colouring``.
    Feasible to about n = 16 for k = 3.
    """
    if n == 0:
        return []
    if k <= 0:
        return None
    for tail in product(range(k), repeat=n - 1):
        colouring = (0,) + tail
        for a, b in edges:
            if colouring[a] == colouring[b]:
                break
        else:
            return list(colouring)
    return None


# --- method 3: chromatic polynomial via deletion/contraction ------------------


def chromatic_polynomial(n, edges):
    """Coefficients of P(G, k), lowest degree first.  Exponential; small graphs only."""
    edges = sorted({(min(a, b), max(a, b)) for a, b in edges})
    return _chrom_poly(n, tuple(edges))


_CP_CACHE = {}


def _chrom_poly(n, edges):
    key = (n, edges)
    if key in _CP_CACHE:
        return _CP_CACHE[key]
    if not edges:
        poly = [0] * n + [1]  # k^n
    else:
        e = edges[0]
        rest = edges[1:]
        deleted = _chrom_poly(n, rest)
        # contract e = (a, b): merge b into a, relabel
        a, b = e
        mapping = {}
        nxt = 0
        for v in range(n):
            if v == b:
                continue
            mapping[v] = nxt
            nxt += 1
        mapping[b] = mapping[a]
        con_edges = set()
        for u, w in rest:
            uu, ww = mapping[u], mapping[w]
            if uu != ww:
                con_edges.add((min(uu, ww), max(uu, ww)))
        contracted = _chrom_poly(n - 1, tuple(sorted(con_edges)))
        poly = _sub(deleted, contracted)
    _CP_CACHE[key] = poly
    return poly


def _sub(p, q):
    m = max(len(p), len(q))
    p = p + [0] * (m - len(p))
    q = q + [0] * (m - len(q))
    return [x - y for x, y in zip(p, q)]


def poly_eval(poly, k):
    total = 0
    for c in reversed(poly):
        total = total * k + c
    return total


def chromatic_number_via_polynomial(n, edges):
    poly = chromatic_polynomial(n, edges)
    for k in range(0, n + 1):
        if poly_eval(poly, k) > 0:
            return k
    return None


# --- method 4: inclusion-exclusion over independent-set counts ---------------


def independent_set_counts(n, edges):
    """i[S] = number of independent subsets of G[S], for every S in 2^V.

    DP on subsets: with v the lowest vertex of S, an independent subset of G[S]
    either omits v (i[S \\ v]) or contains it (i[S \\ (v + N(v))]).
    """
    nbr = [0] * n
    for a, b in edges:
        nbr[a] |= 1 << b
        nbr[b] |= 1 << a
    i = [0] * (1 << n)
    i[0] = 1
    for S in range(1, 1 << n):
        v = (S & -S).bit_length() - 1
        i[S] = i[S & ~(1 << v)] + i[S & ~((1 << v) | nbr[v])]
    return i


def count_k_covers(n, edges, k, counts=None):
    """Number of ordered k-tuples of independent sets whose union is V.

    Inclusion-exclusion on the uncovered part:

        cov_k(G) = sum_{S subset V} (-1)^{|V| - |S|} i(S)^k.

    This is *not* the number of proper colourings -- covers may overlap -- but
    it is zero exactly when G is not k-colourable: a partition into independent
    sets is a cover, and from any cover one recovers a partition by assigning
    each vertex to the first set containing it (a subset of an independent set
    is independent).  So it is a colourability oracle that shares no logic with
    the backtracking search: it never constructs a colouring.

    Time and memory are Theta(2^n), so it is usable to roughly n = 22.
    """
    if counts is None:
        counts = independent_set_counts(n, edges)
    total = 0
    for S in range(1 << n):
        term = counts[S] ** k
        if (n - bin(S).count("1")) & 1:
            total -= term
        else:
            total += term
    return total


def is_k_colourable_ie(n, edges, k, counts=None):
    return count_k_covers(n, edges, k, counts) > 0


def chromatic_number_ie(n, edges):
    counts = independent_set_counts(n, edges)
    for k in range(0, n + 1):
        if count_k_covers(n, edges, k, counts) > 0:
            return k
    return None


# --- critical subgraph extraction --------------------------------------------


def vertex_critical_subgraph(n, edges, k):
    """Shrink to a vertex-minimal subgraph that is still not k-colourable.

    Greedy: try deleting each vertex; keep the deletion when the remainder is
    still non-k-colourable.  The result is (k+1)-vertex-critical in the sense
    that removing any one of its vertices makes it k-colourable.  Returns
    (kept_indices, relabelled_edges) or None if the input *is* k-colourable.
    """
    alive = list(range(n))
    if k_colouring(n, edges, k) is not None:
        return None
    edge_set = {(min(a, b), max(a, b)) for a, b in edges}
    changed = True
    while changed:
        changed = False
        for v in list(alive):
            trial = [u for u in alive if u != v]
            idx = {u: i for i, u in enumerate(trial)}
            sub = [(idx[a], idx[b]) for a, b in edge_set if a in idx and b in idx]
            if k_colouring(len(trial), sub, k) is None:
                alive = trial
                changed = True
    idx = {u: i for i, u in enumerate(alive)}
    sub = sorted((idx[a], idx[b]) for a, b in edge_set if a in idx and b in idx)
    return alive, sub
