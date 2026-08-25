"""A complete k-colouring solver that scales past the exhaustive methods.

`colouring.py` holds four small, obviously-correct algorithms whose job is to
be cross-checkable. This module is the one that has to run on graphs with
thousands of vertices, so it is written for speed and is therefore the one that
must never be trusted alone: anything it decides on a graph small enough for
`colouring.py` should be checked there too, and anything it decides on a large
graph should be reduced to a core small enough to check.

Three reductions, all colourability-preserving, applied before search:

* **k-core.** A vertex of degree < k can always be coloured last, so G is
  k-colourable iff the graph left after repeatedly deleting such vertices is.
  On sparse geometric graphs this is the difference between hopeless and easy.
* **Components.** Solved independently.
* **Unit propagation.** A vertex whose domain drops to one colour is assigned
  immediately rather than branched on.

Search is DSATUR order (smallest remaining domain, then largest degree) with
colour-symmetry breaking: a fresh colour may only ever be the lowest unused
one, which removes the k! redundancy.
"""

from __future__ import annotations

import sys


def build_adjacency_masks(n, edges):
    nbr = [0] * n
    for a, b in edges:
        if a == b:
            raise ValueError(f"self-loop at {a}")
        nbr[a] |= 1 << b
        nbr[b] |= 1 << a
    return nbr


def k_core(n, edges, k):
    """Vertices surviving repeated deletion of degree-<k vertices.

    Returns (kept_sorted, relabelled_edges, index_map, shed_order). G is
    k-colourable iff the returned graph is: a shed vertex had fewer than k
    neighbours among {core} u {shed later}, so colouring in reverse shed order
    always finds it a free colour.
    """
    adj = [set() for _ in range(n)]
    for a, b in edges:
        adj[a].add(b)
        adj[b].add(a)
    alive = [True] * n
    shed_order = []
    stack = [v for v in range(n) if len(adj[v]) < k]
    while stack:
        v = stack.pop()
        if not alive[v]:
            continue
        alive[v] = False
        shed_order.append(v)
        for w in adj[v]:
            if alive[w]:
                adj[w].discard(v)
                if len(adj[w]) < k:
                    stack.append(w)
    kept = [v for v in range(n) if alive[v]]
    idx = {v: i for i, v in enumerate(kept)}
    sub = sorted({(min(idx[a], idx[b]), max(idx[a], idx[b]))
                  for a, b in edges if alive[a] and alive[b]})
    return kept, sub, idx, shed_order


def components(n, edges):
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    seen = [False] * n
    out = []
    for s in range(n):
        if seen[s]:
            continue
        stack, comp = [s], []
        seen[s] = True
        while stack:
            v = stack.pop()
            comp.append(v)
            for w in adj[v]:
                if not seen[w]:
                    seen[w] = True
                    stack.append(w)
        out.append(sorted(comp))
    return out


def greedy_clique(n, nbr):
    order = sorted(range(n), key=lambda v: -bin(nbr[v]).count("1"))
    clique, mask = [], 0
    for v in order:
        if all(nbr[v] >> u & 1 for u in clique):
            clique.append(v)
            mask |= 1 << v
    return clique


def _search(n, nbr, k, node_limit):
    """Complete DSATUR search. Returns a colouring, None, or raises Timeout."""
    full = (1 << k) - 1
    domain = [full] * n
    colour = [-1] * n
    nodes = 0

    clique = greedy_clique(n, nbr)
    if len(clique) > k:
        return None
    trail = []

    def assign(v, c):
        """Assign, propagating; returns list of (w, bit) removals or None on wipeout."""
        removed = []
        colour[v] = c
        bit = 1 << c
        m = nbr[v]
        w = 0
        while m:
            low = m & -m
            w = low.bit_length() - 1
            m ^= low
            if colour[w] < 0 and (domain[w] & bit):
                domain[w] &= ~bit
                removed.append((w, bit))
                if domain[w] == 0:
                    return removed, False
        return removed, True

    def undo(v, removed):
        for w, bit in removed:
            domain[w] |= bit
        colour[v] = -1

    # seed the clique with distinct colours (safe: colours are interchangeable)
    seeded = []
    for i, v in enumerate(clique[:k]):
        if not (domain[v] >> i & 1):
            for v2, r2 in reversed(seeded):
                undo(v2, r2)
            return None
        removed, ok = assign(v, i)
        seeded.append((v, removed))
        if not ok:
            for v2, r2 in reversed(seeded):
                undo(v2, r2)
            return None

    def rec(assigned, used):
        nonlocal nodes
        nodes += 1
        if node_limit and nodes > node_limit:
            raise TimeoutError("node limit exceeded")
        if assigned == n:
            return True
        # unit propagation first, then DSATUR choice
        best, best_key = -1, None
        for v in range(n):
            if colour[v] >= 0:
                continue
            d = domain[v]
            pc = bin(d).count("1")
            if pc == 0:
                return False
            key = (pc, -bin(nbr[v]).count("1"))
            if best_key is None or key < best_key:
                best, best_key = v, key
                if pc == 1:
                    break
        v = best
        d = domain[v]
        limit = min(used + 1, k - 1)
        c = 0
        while d:
            low = d & -d
            c = low.bit_length() - 1
            d ^= low
            if c > limit:
                break
            removed, ok = assign(v, c)
            if ok and rec(assigned + 1, max(used, c)):
                return True
            undo(v, removed)
        return False

    start_assigned = sum(1 for c in colour if c >= 0)
    start_used = max((c for c in colour if c >= 0), default=-1)
    if rec(start_assigned, start_used):
        return list(colour)
    return None


def k_colourable(n, edges, k, node_limit=0, reduce=True):
    """Complete decision. Returns a colouring of the ORIGINAL vertex set, or None.

    The colouring returned is proper on the reduced graph and extended greedily
    back over the deleted low-degree vertices, so it is proper on all of G.
    Raises TimeoutError if node_limit is set and exceeded.
    """
    if n == 0:
        return []
    if k <= 0:
        return None
    edges = sorted({(min(a, b), max(a, b)) for a, b in edges})

    shed_order = []
    cur_n, cur_edges, mapping = n, edges, list(range(n))
    if reduce:
        kept, sub, idx, shed_order = k_core(n, edges, k)
        if len(kept) < n:
            cur_n, cur_edges, mapping = len(kept), sub, kept
        else:
            shed_order = []
    if cur_n == 0:
        core_colouring = []
    else:
        nbr = build_adjacency_masks(cur_n, cur_edges)
        parts = components(cur_n, cur_edges)
        core_colouring = [-1] * cur_n
        for comp in parts:
            ci = {v: i for i, v in enumerate(comp)}
            ce = [(ci[a], ci[b]) for a, b in cur_edges
                  if a in ci and b in ci]
            cn = build_adjacency_masks(len(comp), ce)
            sol = _search(len(comp), cn, k, node_limit)
            if sol is None:
                return None
            for v, c in zip(comp, sol):
                core_colouring[v] = c

    colouring = [-1] * n
    for i, v in enumerate(mapping):
        colouring[v] = core_colouring[i]
    if shed_order:
        adj = [[] for _ in range(n)]
        for a, b in edges:
            adj[a].append(b)
            adj[b].append(a)
        # Reverse shedding order: when v was shed it had < k neighbours among
        # the core and the vertices shed after it, and those are exactly the
        # ones already coloured here, so a free colour always exists.
        for v in reversed(shed_order):
            used = {colouring[w] for w in adj[v] if colouring[w] >= 0}
            free = [c for c in range(k) if c not in used]
            if not free:
                raise AssertionError(
                    "k-core extension failed; the shedding invariant is broken")
            colouring[v] = free[0]
    return colouring


# --- heuristic search for the YES direction ----------------------------------
#
# A colouring is self-certifying: it is checked in linear time by `verify`, so
# a heuristic that finds one has *proved* k-colourability.  Only the NO
# direction needs the complete search.  Splitting the two is what makes graphs
# of this size tractable at all -- the complete search is then only ever run on
# instances a heuristic has failed to colour, which are the interesting ones.


def greedy_random(n, edges, k, rng, nbr=None):
    """Randomised DSATUR: most-saturated first, ties broken at random."""
    if nbr is None:
        nbr = build_adjacency_masks(n, edges)
    colour = [-1] * n
    used_mask = [0] * n
    degree = [bin(x).count("1") for x in nbr]
    remaining = set(range(n))
    while remaining:
        best, best_key = None, None
        for v in remaining:
            key = (-bin(used_mask[v]).count("1"), -degree[v], rng.random())
            if best_key is None or key < best_key:
                best, best_key = v, key
        v = best
        remaining.discard(v)
        free = [c for c in range(k) if not (used_mask[v] >> c & 1)]
        if not free:
            colour[v] = rng.randrange(k)  # leave a conflict for local search
        else:
            colour[v] = rng.choice(free)
        bit = 1 << colour[v]
        m = nbr[v]
        while m:
            low = m & -m
            w = low.bit_length() - 1
            m ^= low
            used_mask[w] |= bit
    return colour


def local_search(n, edges, k, colour, rng, max_steps=200000):
    """Min-conflicts with a random-walk escape.  Returns (colour, conflicts)."""
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    # conflict count per vertex
    def count_at(v, c):
        return sum(1 for w in adj[v] if colour[w] == c)

    # Incremental conflict set: a vertex is "bad" iff it shares its colour with
    # a neighbour. Rebuilding this from scratch each step is what makes a naive
    # min-conflicts loop quadratic, so membership is maintained in place.
    bad_flag = [False] * n
    conflicted = []
    for v in range(n):
        if any(colour[w] == colour[v] for w in adj[v]):
            bad_flag[v] = True
            conflicted.append(v)

    def refresh(u):
        now = any(colour[w] == colour[u] for w in adj[u])
        if now and not bad_flag[u]:
            bad_flag[u] = True
            conflicted.append(u)
        elif not now:
            bad_flag[u] = False

    steps = 0
    while conflicted and steps < max_steps:
        steps += 1
        i = rng.randrange(len(conflicted))
        v = conflicted[i]
        if not bad_flag[v]:
            conflicted[i] = conflicted[-1]
            conflicted.pop()
            continue
        cur = count_at(v, colour[v])
        if cur == 0:
            bad_flag[v] = False
            conflicted[i] = conflicted[-1]
            conflicted.pop()
            continue
        order = list(range(k))
        rng.shuffle(order)
        best_c, best_n = colour[v], cur
        for c in order:
            if c == colour[v]:
                continue
            nc = count_at(v, c)
            if nc < best_n:
                best_c, best_n = c, nc
        if best_n >= cur and rng.random() < 0.3:
            best_c = order[0]
        colour[v] = best_c
        refresh(v)
        for u in adj[v]:
            refresh(u)
    bad_edges = sum(1 for a, b in edges if colour[a] == colour[b])
    return colour, bad_edges


def tabucol(n, edges, k, rng, max_iters=200000, tenure=10, tenure_frac=0.6):
    """Tabu search over improper k-colourings, minimising conflicting edges.

    The standard formulation, and far stronger than min-conflicts on graph
    colouring: it keeps gamma[v][c] = how many neighbours of v have colour c,
    so the change in conflicts from recolouring v to c is a table lookup, and
    it forbids returning (v, old colour) for a few iterations to escape the
    plateaus that trap a greedy descent.

    Returns (colour, conflicts).  Conflicts zero means a proper colouring, and
    that is a *proof* of k-colourability once `verify` confirms it; a nonzero
    result proves nothing at all.
    """
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    colour = [rng.randrange(k) for _ in range(n)]
    gamma = [[0] * k for _ in range(n)]
    for v in range(n):
        for w in adj[v]:
            gamma[v][colour[w]] += 1
    conflicts = sum(1 for a, b in edges if colour[a] == colour[b])
    best = conflicts
    tabu = {}
    for it in range(max_iters):
        if conflicts == 0:
            break
        # candidate moves: conflicting vertices only
        best_move, best_delta = None, None
        cand = [v for v in range(n) if gamma[v][colour[v]] > 0]
        if not cand:
            break
        # sample a bounded number of candidates so an iteration stays cheap
        if len(cand) > 60:
            cand = rng.sample(cand, 60)
        for v in cand:
            cur = gamma[v][colour[v]]
            for c in range(k):
                if c == colour[v]:
                    continue
                delta = gamma[v][c] - cur
                is_tabu = tabu.get((v, c), 0) > it
                if is_tabu and not (conflicts + delta < best):
                    continue
                if best_delta is None or delta < best_delta:
                    best_move, best_delta = (v, c), delta
        if best_move is None:
            v = cand[rng.randrange(len(cand))]
            c = rng.randrange(k)
            best_move, best_delta = (v, c), gamma[v][c] - gamma[v][colour[v]]
        v, c = best_move
        old = colour[v]
        colour[v] = c
        conflicts += best_delta
        for w in adj[v]:
            gamma[w][old] -= 1
            gamma[w][c] += 1
        tabu[(v, old)] = it + tenure + int(tenure_frac * len(cand))
        if conflicts < best:
            best = conflicts
    bad = sum(1 for a, b in edges if colour[a] == colour[b])
    return colour, bad


def find_colouring(n, edges, k, tries=8, steps=200000, seed=0):
    """Heuristic search for a proper k-colouring.  None means 'not found', NOT
    'does not exist' -- only the complete search can say that."""
    import random

    if n == 0:
        return []
    nbr = build_adjacency_masks(n, edges)
    rng = random.Random(seed)
    for t in range(tries):
        col = greedy_random(n, edges, k, rng, nbr)
        if verify(n, edges, col, k):
            return col
        col2, bad = local_search(n, edges, k, list(col), rng, max_steps=steps)
        if bad == 0 and verify(n, edges, col2, k):
            return col2
        col3, bad3 = tabucol(n, edges, k, rng, max_iters=steps)
        if bad3 == 0 and verify(n, edges, col3, k):
            return col3
    return None


def decide_k_colourable(n, edges, k, tries=8, steps=200000, seed=0,
                        node_limit=0, quick_nodes=30000):
    """(verdict, colouring, how) for k-colourability.

    Order matters for speed, and the order is: cheap-and-complete, then
    heuristic, then expensive-and-complete.

    1. The complete DSATUR search under a small node budget. With propagation
       it is itself a strong heuristic, so easy YES instances fall out at once,
       and small NO instances are *settled* here rather than merely suspected.
    2. If that budget is exhausted, randomised greedy plus local search. A
       colouring it finds is verified, so a YES from here is as good as a YES
       from anywhere.
    3. Only if the heuristics also fail does the complete search run under the
       full budget, which is the only way to reach a NO on a hard instance.

    ``verdict`` is None when nothing settled it; that is reported, never guessed.
    """
    edges = sorted({(min(a, b), max(a, b)) for a, b in edges})
    try:
        col = k_colourable(n, edges, k, node_limit=quick_nodes)
        return (col is not None), col, "complete-quick"
    except TimeoutError:
        pass

    kept, sub, idx, shed = k_core(n, edges, k)
    if kept:
        core_col = find_colouring(len(kept), sub, k, tries=tries,
                                  steps=steps, seed=seed)
        if core_col is not None:
            full = _extend_core_colouring(n, edges, k, kept, core_col, shed)
            if full is not None and verify(n, edges, full, k):
                return True, full, "heuristic"

    try:
        full = k_colourable(n, edges, k, node_limit=node_limit)
    except TimeoutError:
        return None, None, "timeout"
    return (full is not None), full, "complete"


def _extend_core_colouring(n, edges, k, kept, core_col, shed_order):
    """Push a colouring of the k-core back out over the shed vertices."""
    colouring = [-1] * n
    for i, v in enumerate(kept):
        colouring[v] = core_col[i]
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    for v in reversed(shed_order):
        used = {colouring[w] for w in adj[v] if colouring[w] >= 0}
        free = [c for c in range(k) if c not in used]
        if not free:
            return None
        colouring[v] = free[0]
    return colouring


def chromatic_number(n, edges, lo=1, hi=8, node_limit=0):
    for k in range(lo, hi + 1):
        if k_colourable(n, edges, k, node_limit=node_limit) is not None:
            return k
    return None


def verify(n, edges, colouring, k=None):
    if colouring is None or len(colouring) != n:
        return False
    if any(c < 0 for c in colouring):
        return False
    if k is not None and any(c >= k for c in colouring):
        return False
    return all(colouring[a] != colouring[b] for a, b in edges)


def contract(n, edges, u, v):
    """Merge v into u. Returns (n-1, edges) with vertices relabelled."""
    if u == v:
        raise ValueError("cannot contract a vertex with itself")
    u, v = min(u, v), max(u, v)
    idx = {}
    nxt = 0
    for w in range(n):
        if w == v:
            continue
        idx[w] = nxt
        nxt += 1
    idx[v] = idx[u]
    out = set()
    for a, b in edges:
        aa, bb = idx[a], idx[b]
        if aa != bb:
            out.add((min(aa, bb), max(aa, bb)))
    return n - 1, sorted(out)


sys.setrecursionlimit(100000)
