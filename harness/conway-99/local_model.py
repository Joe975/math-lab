"""The local decomposition of a strongly regular graph with lam = 1, mu = 2.

Around any fixed vertex v0 such a graph splits into three layers, and the two
outer layers carry forced combinatorial structure.  Writing k = 2m:

  * N(v0) is 1-regular, so it is a perfect matching m*K2 on k vertices.  Write
    x' for the matching partner of x.
  * A vertex d at distance 2 has exactly mu = 2 neighbours in N(v0), and those
    two are never partners (a partner pair together with v0 would give an
    adjacent pair two common neighbours).
  * Conversely a non-partner pair {x, y} in N(v0) is a non-adjacent pair whose
    common neighbours are exactly v0 and one further vertex, at distance 2.

  So d -> {its two neighbours in N(v0)} is a **bijection** from the distance-2
  layer onto the non-partner pairs of N(v0).  The counts agree identically:
  C(k,2) - k/2 = n - k - 1 whenever k(k-2) = 2(n-k-1).

  * Consequently, for each x in N(v0) the set A_x of distance-2 vertices whose
    pair contains x has size k - 2, and the subgraph it induces is a perfect
    matching M_x, which we read as a perfect matching on N(v0) \\ {x, x'}.

This module builds that decomposition from an explicit graph and checks every
one of those statements, so the reduction can be validated on graphs whose
existence is published before it is used on a parameter set where no graph is
in hand.
"""


def partner_map(adj, v0):
    """The matching on N(v0).  Returns {x: x'} or raises if it is not 1-regular."""
    nb = [u for u in range(len(adj)) if adj[v0] >> u & 1]
    nbset = set(nb)
    partner = {}
    for x in nb:
        inside = [u for u in nb if adj[x] >> u & 1]
        if len(inside) != 1:
            raise ValueError(f"N({v0}) is not a perfect matching: "
                             f"vertex {x} has {len(inside)} neighbours inside")
        partner[x] = inside[0]
    assert all(partner[partner[x]] == x for x in nb)
    assert set(partner) == nbset
    return partner


def decompose(adj, v0):
    """Full local decomposition around v0.

    Returns a dict with:
      nb        -- list of N(v0)
      partner   -- {x: x'}
      far       -- list of distance-2 vertices
      pair_of   -- {d: frozenset({x, y})}, the two N(v0)-neighbours of d
      vert_of   -- inverse map {frozenset({x,y}): d}
      match     -- {x: {y: z}}, the perfect matching M_x on N(v0)\\{x,x'},
                   given as an involution on labels
    Raises ValueError if any forced property fails.
    """
    n = len(adj)
    partner = partner_map(adj, v0)
    nb = sorted(partner)
    k = len(nb)
    far = [u for u in range(n) if u != v0 and not (adj[v0] >> u & 1)]

    pair_of, vert_of = {}, {}
    for d in far:
        common = [x for x in nb if adj[d] >> x & 1]
        if len(common) != 2:
            raise ValueError(f"vertex {d} has {len(common)} neighbours in N({v0}), expected 2")
        x, y = common
        if partner[x] == y:
            raise ValueError(f"vertex {d} sees the partner pair {{{x},{y}}}")
        key = frozenset((x, y))
        if key in vert_of:
            raise ValueError(f"pair {set(key)} is realised twice: {vert_of[key]} and {d}")
        pair_of[d] = key
        vert_of[key] = d

    expected = k * (k - 1) // 2 - k // 2
    if len(far) != expected:
        raise ValueError(f"distance-2 layer has {len(far)} vertices, "
                         f"non-partner pairs number {expected}")
    if len(vert_of) != expected:
        raise ValueError("pair map is not a bijection")

    # the perfect matchings M_x
    match = {}
    for x in nb:
        A = [d for d in far if x in pair_of[d]]
        if len(A) != k - 2:
            raise ValueError(f"A_{x} has size {len(A)}, expected {k-2}")
        m = {}
        for d in A:
            inside = [e for e in A if adj[d] >> e & 1]
            if len(inside) != 1:
                raise ValueError(f"A_{x} is not a perfect matching: "
                                 f"{d} has {len(inside)} neighbours inside")
            e = inside[0]
            (y,) = pair_of[d] - {x}
            (z,) = pair_of[e] - {x}
            m[y] = z
        assert all(m[m[y]] == y and m[y] != y for y in m)
        assert set(m) == set(nb) - {x, partner[x]}
        match[x] = m

    return {"v0": v0, "nb": nb, "partner": partner, "far": far,
            "pair_of": pair_of, "vert_of": vert_of, "match": match, "k": k}


def cycle_type(partner, m, x):
    """Cycle type of M_x against the partner matching, as a partition.

    The union of the two perfect matchings M_x and the partner matching, both
    on N(v0) \\ {x, x'}, is a disjoint union of even cycles; halving each cycle
    length gives a partition of (k-2)/2.  A part equal to 1 means M_x and the
    partner matching share an edge.
    """
    dom = set(m)
    seen, parts = set(), []
    for y in sorted(dom):
        if y in seen:
            continue
        c, cur = 0, y
        while cur not in seen:
            seen.add(cur)
            cur = m[cur]
            seen.add(cur)
            c += 1
            cur = partner[cur]
        parts.append(c)
    return tuple(sorted(parts, reverse=True))


def verify_decomposition(adj, v0=0, verbose=False):
    """Run `decompose` and report; returns (ok, info-or-error)."""
    try:
        dec = decompose(adj, v0)
    except ValueError as e:
        return False, str(e)
    types = {}
    for x in dec["nb"]:
        t = cycle_type(dec["partner"], dec["match"][x], x)
        types[t] = types.get(t, 0) + 1
    dec["cycle_types"] = types
    return True, dec
