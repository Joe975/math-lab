"""The A/B control that makes the k=14 search measurement mean anything.

Run the *identical* search under the *identical* hypothesis (partner-regularity
at v0, i.e. M_x = the partner matching for every x in P) on:

    k = 22  ->  n = 243   -- SATISFIABLE: the Berlekamp-van Lint-Seidel graph
                             is partner-regular, so a solution exists and the
                             search is allowed to find one.
    k = 14  ->  n =  99   -- the open case.

If the search cannot make progress on k=22 either, then any statement about
where it dies on k=14 is a statement about this search, not about 99.  That is
the difference between an obstruction and a tooling limit, and the verification
contract requires it be settled explicitly.
"""

import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg
import pair_search as PS
from partner_regular import impose_partner_regular, stats


def run(k, budget, order="vertex"):
    M = PS.Model(k)
    t0 = time.time()
    try:
        st = impose_partner_regular(M)
    except PS.Contradiction:
        return {"k": k, "n": M.n, "result": "contradiction-at-seed",
                "time": time.time() - t0}
    s = stats(M, st)
    stt = {"nodes": 0, "fails": 0, "max_depth": 0, "max_settled": 0,
           "solutions": 0, "level_nodes": {}, "level_fails": {}, "witness": None}
    r = PS.search(st, budget, stt, 0, order)
    out = {"k": k, "n": M.n, "nd": M.nd, "result": r, "time": time.time() - t0,
           "open_after_seed": s["open"], **stt}
    if stt["witness"] is not None:
        adj = PS.reconstruct(M, stt["witness"])
        ok, got = srg.check_srg(adj)
        out["verified"] = (ok, got)
    return out


def main():
    budget = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    order = sys.argv[2] if len(sys.argv) > 2 else "vertex"
    print(f"# hypothesis: partner-regular at v0 (M_x = partner matching, all x)")
    print(f"# budget = {budget} nodes per instance, branch order = {order}")
    for k in (4, 22, 14):
        r = run(k, budget, order)
        tag = {4: "SATISFIABLE (Paley 9)", 22: "SATISFIABLE (BvLS 243)",
               14: "OPEN (the target)"}[k]
        print(f"\nk={k:3d}  n={r['n']:4d}  {tag}")
        print(f"  result={r['result']}  nodes={r.get('nodes')}  "
              f"fails={r.get('fails')}  maxdepth={r.get('max_depth')}  "
              f"settled_prefix={r.get('max_settled')}/{r.get('nd')}  "
              f"t={r['time']:.1f}s")
        if "verified" in r:
            print(f"  witness independently checked: {r['verified']}")
        ln = r.get("level_nodes", {})
        if ln:
            print("  nodes by settled-prefix length: " +
                  ", ".join(f"{a}:{b}" for a, b in sorted(ln.items())[:12]))


if __name__ == "__main__":
    main()
