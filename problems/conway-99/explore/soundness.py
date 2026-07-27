"""Soundness control for the pair-model propagator.

The propagator is only informative if it never rejects a graph that really
exists.  This script takes the realised members of the lam=1, mu=2 family,
transports each into the pair model's own labelling, asserts the whole
D-adjacency into a fresh search state, and checks that

  (a) propagation reports no contradiction,
  (b) the state comes out complete and consistent,
  (c) reconstructing from that state returns the original graph up to the
      relabelling used,

so that any pruning the search does on k=14 is pruning the propagator is
entitled to do.
"""

import sys
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg
import constructions as C
import local_model as L
import pair_search as PS


def transport(adj, v0=0):
    """Relabel an SRG(n,k,1,2) into the pair model's labelling.

    Returns (Model, label_of_P_vertex, dindex_of_far_vertex).
    P is numbered so that partner(x) == x ^ 1, matching the model.
    """
    dec = L.decompose(adj, v0)
    partner = dec["partner"]
    k = dec["k"]
    M = PS.Model(k)
    # order the matching edges arbitrarily but consistently
    edges = sorted({tuple(sorted((x, partner[x]))) for x in dec["nb"]})
    label = {}
    for t, (a, b) in enumerate(edges):
        label[a] = 2 * t
        label[b] = 2 * t + 1
    assert all(label[partner[x]] == label[x] ^ 1 for x in dec["nb"])
    dindex = {}
    for d in dec["far"]:
        key = frozenset(label[x] for x in dec["pair_of"][d])
        dindex[d] = M.index[key]
    assert len(set(dindex.values())) == M.nd
    return M, label, dindex, dec


def plant(adj, name, v0=0):
    M, label, dindex, dec = transport(adj, v0)
    st = PS.State(M)
    q = []
    far = dec["far"]
    try:
        for a in range(len(far)):
            for b in range(a + 1, len(far)):
                u, v = far[a], far[b]
                val = bool(adj[u] >> v & 1)
                st.assign(dindex[u], dindex[v], val, q)
        PS.propagate(st, q)
    except PS.Contradiction:
        return False, f"{name}: PROPAGATOR REJECTED A REAL GRAPH"
    # completeness
    incomplete = [v for v in range(M.nd) if st.unknown(v)]
    if incomplete:
        return False, f"{name}: state incomplete at {len(incomplete)} vertices"
    # round trip
    rebuilt = PS.reconstruct(M, [st.yes[v] for v in range(M.nd)])
    ok, got = srg.check_srg(rebuilt)
    if not ok:
        return False, f"{name}: reconstruction is not an SRG: {got}"
    return True, f"{name}: accepted, reconstructed as SRG{got}"


def main():
    print("=== propagator soundness: real graphs must survive ===")
    cases = [("paley9", C.paley9()), ("rook3", C.rook(3)),
             ("bvls243", C.bvls243()[0])]
    allok = True
    for name, adj in cases:
        # test around several base vertices, not just one
        for v0 in (0, 1, 5):
            ok, msg = plant(adj, f"{name} @v0={v0}", v0)
            print(("  OK   " if ok else "  FAIL ") + msg)
            allok &= ok
    print()
    print("=== does propagation ALONE (no branching) finish a real graph "
          "from only its 22/4 matchings? ===")
    for name, adj in cases:
        M, label, dindex, dec = transport(adj, 0)
        st = PS.State(M)
        q = []
        try:
            # assert only the label-sharing adjacency, i.e. every M_x
            far = dec["far"]
            for a in range(len(far)):
                for b in range(a + 1, len(far)):
                    u, v = far[a], far[b]
                    i, j = dindex[u], dindex[v]
                    if M.share[i][j] is not None:
                        st.assign(i, j, bool(adj[u] >> v & 1), q)
            PS.propagate(st, q)
        except PS.Contradiction:
            print(f"  {name}: contradiction (should not happen)")
            allok = False
            continue
        unk = sum(PS.popcount(st.unknown(v)) for v in range(M.nd)) // 2
        tot = M.nd * (M.nd - 1) // 2
        print(f"  {name}: all {M.k} matchings given -> {tot-unk}/{tot} pairs "
              f"decided, {unk} still open")
    print()
    print("ALL SOUNDNESS CHECKS PASSED" if allok else "SOUNDNESS FAILURE")
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
