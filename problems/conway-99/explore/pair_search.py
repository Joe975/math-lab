"""Constraint search for SRG(n, k, 1, 2) in the forced pair model.

The seed is completely forced up to isomorphism (see harness/conway-99/local_model.py):
v0, its k neighbours P carrying a fixed perfect matching (k/2 disjoint edges),
and n-k-1 further vertices D in bijection with the non-partner pairs of P.
Everything except the adjacency *within* D is fixed, so the search space is
exactly the C(|D|,2) binary variables on D x D.

The propagator is exact-count based, not clausal.  For each pair (u,v) of
D-vertices it maintains

    lo(u,v) = |yes(u) & yes(v)|      certain common neighbours inside D
    hi(u,v) = |poss(u) & poss(v)|    possible common neighbours inside D

against the target, which is forced by the strongly-regular conditions plus the
common neighbours the pair already has inside P:

    u,v share a P-label  -> exactly 1 common neighbour in P
    u,v disjoint         -> 0 common neighbours in P
    target = (1 if adjacent else 2) - (1 if sharing a label else 0)

k is a parameter so that the identical machinery can be run on the parameter
sets in this family whose answer is published (k=4 -> n=9, k=22 -> n=243)
before it is run on k=14 -> n=99.

NOTHING here claims exhaustion.  A "budget" result means the search was cut
off; an "exhausted" result from this code alone is not treated as a proof of
nonexistence, per the verification contract.
"""

import argparse
import sys
import time
from itertools import combinations


class Model:
    """The forced pair model for SRG((k*k+2)/2, k, 1, 2)."""

    def __init__(self, k):
        assert k % 2 == 0
        self.k = k
        self.n = (k * k + 2) // 2
        self.npos = k // 2
        self.pairs = [frozenset((x, y)) for x, y in combinations(range(k), 2)
                      if y != (x ^ 1)]
        self.nd = len(self.pairs)
        assert self.nd == self.n - k - 1, (self.nd, self.n - k - 1)
        self.index = {p: i for i, p in enumerate(self.pairs)}
        self.full = (1 << self.nd) - 1
        self.deg_d = k - 2
        # shared-label table
        self.share = [[None] * self.nd for _ in range(self.nd)]
        for i, j in combinations(range(self.nd), 2):
            s = self.pairs[i] & self.pairs[j]
            v = next(iter(s)) if len(s) == 1 else None
            self.share[i][j] = self.share[j][i] = v
        # A_x
        self.A = [0] * k
        for i, p in enumerate(self.pairs):
            for x in p:
                self.A[x] |= 1 << i

    def target(self, i, j, adjacent):
        return (1 if adjacent else 2) - (1 if self.share[i][j] is not None else 0)


def popcount(x):
    return bin(x).count("1")


class Contradiction(Exception):
    pass


class State:
    __slots__ = ("M", "yes", "no", "trail", "props")

    def __init__(self, M):
        self.M = M
        self.yes = [0] * M.nd
        self.no = [1 << i for i in range(M.nd)]
        self.trail = []
        self.props = 0

    def poss(self, i):
        return self.M.full & ~self.no[i]

    def unknown(self, i):
        return self.M.full & ~self.yes[i] & ~self.no[i]

    def status(self, i, j):
        if self.yes[i] >> j & 1:
            return True
        if self.no[i] >> j & 1:
            return False
        return None

    def mark(self):
        return len(self.trail)

    def undo(self, m):
        while len(self.trail) > m:
            i, j, val = self.trail.pop()
            if val:
                self.yes[i] &= ~(1 << j)
                self.yes[j] &= ~(1 << i)
            else:
                self.no[i] &= ~(1 << j)
                self.no[j] &= ~(1 << i)

    def assign(self, i, j, val, queue):
        cur = self.status(i, j)
        if cur is not None:
            if cur != val:
                raise Contradiction()
            return
        if val:
            self.yes[i] |= 1 << j
            self.yes[j] |= 1 << i
        else:
            self.no[i] |= 1 << j
            self.no[j] |= 1 << i
        self.trail.append((i, j, val))
        self.props += 1
        queue.append((i, j))


def propagate(st, queue):
    M = st.M
    deg = M.deg_d
    while queue:
        i, j = queue.pop()
        for v in (i, j):
            y = popcount(st.yes[v])
            if y > deg:
                raise Contradiction()
            p = popcount(st.poss(v))
            if p < deg:
                raise Contradiction()
            if y == deg:
                unk = st.unknown(v)
                while unk:
                    b = unk & -unk
                    unk ^= b
                    st.assign(v, b.bit_length() - 1, False, queue)
            elif p == deg:
                unk = st.unknown(v)
                while unk:
                    b = unk & -unk
                    unk ^= b
                    st.assign(v, b.bit_length() - 1, True, queue)
        for v in (i, j):
            _layer_rule(st, v, queue)
            for w in range(M.nd):
                if w != v:
                    _pair_rule(st, v, w, queue)


def _layer_rule(st, d, queue):
    """Enforce the mu-condition between the D-vertex d and each P-vertex z.

    For z in P not on d's pair, d is non-adjacent to z, so d and z have exactly
    2 common neighbours.  N(z) = {v0} u {z'} u A_z, and v0 is never adjacent to
    d, while z' is adjacent to d exactly when z' lies on d's pair.  Hence

        |N(d) & A_z| = 2 - [z' on d's pair].

    (For z on d's pair the same argument with lambda = 1 gives |N(d) & A_z| = 1,
    which the matching M_z already encodes.)

    These are cardinality constraints across a whole layer, and they are NOT
    implied by the D x D pairwise rules -- without them the propagator misses
    forced non-adjacencies.
    """
    M = st.M
    pair = M.pairs[d]
    yes, poss = st.yes[d], st.poss(d)
    for z in range(M.k):
        if z in pair:
            continue
        az = M.A[z]
        t = 2 - (1 if (z ^ 1) in pair else 0)
        lo = popcount(yes & az)
        hi = popcount(poss & az)
        if lo > t or hi < t:
            raise Contradiction()
        if lo == t and hi > t:
            s = (poss & az) & ~yes
            while s:
                b = s & -s
                s ^= b
                st.assign(d, b.bit_length() - 1, False, queue)
            yes, poss = st.yes[d], st.poss(d)
        elif hi == t and lo < t:
            s = (poss & az) & ~yes
            while s:
                b = s & -s
                s ^= b
                st.assign(d, b.bit_length() - 1, True, queue)
            yes, poss = st.yes[d], st.poss(d)


def _pair_rule(st, u, v, queue):
    M = st.M
    stt = st.status(u, v)
    yu, yv = st.yes[u], st.yes[v]
    pu, pv = st.poss(u), st.poss(v)
    both_y = yu & yv
    both_p = pu & pv
    lo = popcount(both_y)
    hi = popcount(both_p)
    if stt is None:
        t1 = M.target(u, v, True)
        t0 = M.target(u, v, False)
        ok1 = lo <= t1 <= hi
        ok0 = lo <= t0 <= hi
        if not ok1 and not ok0:
            raise Contradiction()
        if not ok1:
            st.assign(u, v, False, queue)
        elif not ok0:
            st.assign(u, v, True, queue)
        return
    t = M.target(u, v, stt)
    if lo > t or hi < t:
        raise Contradiction()
    if lo == t and hi > t:
        s = both_p & ~both_y
        while s:
            b = s & -s
            s ^= b
            w = b.bit_length() - 1
            if (yu >> w) & 1:
                st.assign(v, w, False, queue)
            elif (yv >> w) & 1:
                st.assign(u, w, False, queue)
    elif hi == t and lo < t:
        s = both_p & ~both_y
        while s:
            b = s & -s
            s ^= b
            w = b.bit_length() - 1
            st.assign(u, w, True, queue)
            st.assign(v, w, True, queue)


# ---------------------------------------------------------------------------
# seeding: M_{x0} up to the seed stabiliser
# ---------------------------------------------------------------------------

def perfect_matchings(S):
    if not S:
        yield {}
        return
    a = S[0]
    for i in range(1, len(S)):
        b = S[i]
        rest = S[1:i] + S[i + 1:]
        for m in perfect_matchings(rest):
            m = dict(m)
            m[a] = b
            m[b] = a
            yield m


def cycle_type(m):
    """Cycle type of m against the partner matching: a partition of |m|/2."""
    seen, parts = set(), []
    for y in sorted(m):
        if y in seen:
            continue
        c, cur = 0, y
        while cur not in seen:
            seen.add(cur)
            cur = m[cur]
            seen.add(cur)
            c += 1
            cur = cur ^ 1
        parts.append(c)
    return tuple(sorted(parts, reverse=True))


def matching_seeds(M):
    """One canonical M_{x0} per orbit of the stabiliser of x0 = 0.

    The stabiliser of the seed (v0, P, its matching) is the hyperoctahedral
    group B_m, m = k/2: permute the m matching edges and flip signs inside
    each, order m! * 2^m.  Its stabiliser of the vertex x0 = 0 is B_{m-1} on
    the remaining m-1 edges.  Two perfect matchings of P \\ {0,1} are in the
    same orbit exactly when their union with the partner matching has the same
    cycle type, so the partitions of m-1 index the orbits exhaustively.
    """
    S = [p for p in range(M.k) if p not in (0, 1)]
    seeds = {}
    for m in perfect_matchings(S):
        t = cycle_type(m)
        seeds.setdefault(t, m)
    return seeds


def apply_matching(st, x, m, queue):
    M = st.M
    ax = [i for i in range(M.nd) if (M.A[x] >> i) & 1]
    for i, j in combinations(ax, 2):
        (yi,) = M.pairs[i] - {x}
        (yj,) = M.pairs[j] - {x}
        st.assign(i, j, m.get(yi) == yj, queue)


# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------

def choose_branch(st, order="vertex"):
    """Pick the next undecided pair to branch on.

    order="vertex"  -- canonical extension: complete the neighbourhood of
                       D-vertex 0, then 1, then 2, ...  This makes "how many
                       D-vertices are settled" a meaningful depth measure, so
                       the search can be reported per extension level.
    order="mrv"     -- most-constrained-vertex first; faster at finding a
                       contradiction but gives no level structure.
    """
    M = st.M
    if order == "vertex":
        for v in range(M.nd):
            unk = st.unknown(v)
            if unk:
                return v, (unk & -unk).bit_length() - 1
        return None
    best, bestv = None, None
    for v in range(M.nd):
        unk = st.unknown(v)
        if not unk:
            continue
        need = M.deg_d - popcount(st.yes[v])
        key = (popcount(unk) - need, popcount(unk))
        if best is None or key < best:
            best, bestv = key, v
    if bestv is None:
        return None
    unk = st.unknown(bestv)
    return bestv, (unk & -unk).bit_length() - 1


def settled_prefix(st):
    """Number of leading D-vertices whose neighbourhood is fully decided."""
    M = st.M
    c = 0
    while c < M.nd and not st.unknown(c):
        c += 1
    return c


def search(st, budget, stats, depth=0, order="vertex", values=(False, True)):
    if stats["nodes"] >= budget:
        return "budget"
    stats["nodes"] += 1
    if depth > stats["max_depth"]:
        stats["max_depth"] = depth
    lvl = settled_prefix(st)
    if lvl > stats["max_settled"]:
        stats["max_settled"] = lvl
    # nodes spent while the settled prefix had this length
    stats["level_nodes"][lvl] = stats["level_nodes"].get(lvl, 0) + 1

    br = choose_branch(st, order)
    if br is None:
        stats["solutions"] += 1
        stats["witness"] = [st.yes[v] for v in range(st.M.nd)]
        return "solved"
    v, w = br
    for val in values:
        mk = st.mark()
        try:
            q = []
            st.assign(v, w, val, q)
            propagate(st, q)
        except Contradiction:
            st.undo(mk)
            stats["fails"] += 1
            stats["level_fails"][lvl] = stats["level_fails"].get(lvl, 0) + 1
            continue
        r = search(st, budget, stats, depth + 1, order, values)
        st.undo(mk)
        if r in ("solved", "budget"):
            return r
    return "exhausted"


def run_seed(M, ct, m, budget, order="vertex"):
    st = State(M)
    stats = {"nodes": 0, "fails": 0, "max_depth": 0, "max_settled": 0,
             "solutions": 0, "level_nodes": {}, "level_fails": {},
             "witness": None}
    t0 = time.time()
    try:
        q = []
        apply_matching(st, 0, m, q)
        propagate(st, q)
    except Contradiction:
        return {"cycle_type": ct, "result": "contradiction-at-seed",
                "time": time.time() - t0, **stats}
    decided = sum(popcount(st.yes[v]) + popcount(st.no[v]) - 1
                  for v in range(M.nd)) // 2
    stats["decided_pairs_after_seed"] = decided
    stats["settled_after_seed"] = sum(1 for v in range(M.nd) if not st.unknown(v))
    r = search(st, budget, stats, 0, order)
    return {"cycle_type": ct, "result": r, "time": time.time() - t0, **stats}


def reconstruct(M, yes):
    """Turn a completed D-adjacency into the full n-vertex graph.

    Vertex 0 is v0, 1..k are P (partner of 1+t is 1+(t^1)), and k+1..n-1 are D
    in the order of M.pairs.  Returned as a list of int bitmasks, ready for the
    independent checker in harness/conway-99/srg.py.
    """
    n = M.n
    adj = [0] * n

    def link(a, b):
        adj[a] |= 1 << b
        adj[b] |= 1 << a

    for x in range(M.k):
        link(0, 1 + x)
    for x in range(0, M.k, 2):
        link(1 + x, 1 + (x ^ 1))
    for i, p in enumerate(M.pairs):
        for x in p:
            link(1 + x, 1 + M.k + i)
    for i in range(M.nd):
        b = yes[i]
        while b:
            lb = b & -b
            b ^= lb
            j = lb.bit_length() - 1
            if j > i:
                link(1 + M.k + i, 1 + M.k + j)
    return adj


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-k", type=int, default=14)
    ap.add_argument("--budget", type=int, default=50000)
    ap.add_argument("--only", type=str, default=None)
    ap.add_argument("--order", type=str, default="vertex", choices=["vertex", "mrv"])
    ap.add_argument("--levels", action="store_true", help="print per-level node counts")
    args = ap.parse_args()

    M = Model(args.k)
    print(f"# k={M.k}  n={M.n}  |D|={M.nd}  pair variables={M.nd*(M.nd-1)//2}")
    seeds = matching_seeds(M)
    print(f"# canonical M_x0 seeds (orbits of the seed stabiliser): {len(seeds)}")
    for t in sorted(seeds, reverse=True):
        if args.only and ",".join(map(str, t)) != args.only:
            continue
        res = run_seed(M, t, seeds[t], args.budget, args.order)
        print(f"  seed {str(t):24s} result={res['result']:12s} "
              f"nodes={res['nodes']:8d} fails={res['fails']:8d} "
              f"maxdepth={res['max_depth']:5d} "
              f"settled={res['max_settled']:4d}/{M.nd} "
              f"sols={res['solutions']} t={res['time']:.1f}s")
        if args.levels:
            ln = res["level_nodes"]
            print("      nodes by settled-prefix length: " +
                  ", ".join(f"{k}:{ln[k]}" for k in sorted(ln)))
        sys.stdout.flush()


if __name__ == "__main__":
    main()
