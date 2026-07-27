"""Regenerate problems/conway-99/data/pair-model-evidence.json.

Deterministic: no randomness, no wall-clock in the payload.  Run with

    python problems/conway-99/explore/make_evidence.py

Takes a few minutes, dominated by the k=22 control.
"""

import json
import os
import sys
from itertools import combinations, product

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg
import constructions as C
import local_model as L
import pair_search as PS
from partner_regular import impose_partner_regular, stats

BUDGET = 20000
out = {}

# --- feasibility of the family -------------------------------------------
fam = []
for (n, k, lam, mu, rep) in srg.enumerate_feasible(1, 2, 2_000_000):
    s = rep["spectrum"]
    fam.append({"n": n, "k": k, "r": s["r"], "f": s["f"], "s": s["s"],
                "g": s["g"], "conference": s.get("conference", False)})
out["lam1_mu2_family_n_le_2e6"] = fam
ok, rep = srg.feasibility(99, 14, 1, 2)
out["target_feasibility"] = {
    "params": [99, 14, 1, 2], "feasible": ok,
    "conditions": [{"name": c["name"], "ok": c["ok"]} for c in rep["conditions"]],
    "spectrum": {k2: v for k2, v in rep["spectrum"].items() if k2 != "reason"},
    "note": ("multiplicities are 54 and 44; they satisfy f+g=n-1, tr(A)=0 and "
             "tr(A^2)=nk, each checked independently"),
}

# --- calibration: constructions ------------------------------------------
cal = {}
for name, (build, params) in C.CATALOGUE.items():
    adj = build()
    good, got = srg.check_srg(adj, *params)
    cal[name] = {"claimed": list(params), "verified": bool(good),
                 "measured": list(got) if good else str(got)}
out["construction_calibration"] = cal

# --- the local reduction, validated on realised graphs --------------------
red = {}
for name, adj in [("paley9", C.paley9()), ("rook3", C.rook(3)),
                  ("bvls243", C.bvls243()[0])]:
    allok, types = True, {}
    for v0 in range(len(adj)):
        o, dec = L.verify_decomposition(adj, v0)
        if not o:
            allok = False
            break
        for t, c in dec["cycle_types"].items():
            types[str(list(t))] = types.get(str(list(t)), 0) + c
    red[name] = {"n": len(adj), "decomposition_holds_at_every_vertex": allok,
                 "M_x_cycle_types_vs_partner_matching": types}
out["local_reduction_validation"] = red

# --- partner-regular closure ---------------------------------------------
pr = {}
for k in (4, 14, 22):
    M = PS.Model(k)
    try:
        st = impose_partner_regular(M)
    except PS.Contradiction:
        pr[str(k)] = {"contradiction": True}
        continue
    s = stats(M, st)
    need = sorted({M.deg_d - PS.popcount(st.yes[v]) for v in range(M.nd)})
    cand = sorted({PS.popcount(st.unknown(v)) for v in range(M.nd)})
    pr[str(k)] = {"n": M.n, "D": M.nd, **s,
                  "neighbours_still_needed_per_vertex": need,
                  "candidates_per_vertex": cand}
out["partner_regular_closure"] = pr

# --- the n_B profile case split ------------------------------------------
R = list(range(5))
blocks = list(combinations(R, 2))
sols = []
for cur in product((0, 1, 2), repeat=10):
    if sum(cur) != 10:
        continue
    degs = [0] * 5
    for m, (a, b) in zip(cur, blocks):
        degs[a] += m
        degs[b] += m
    if all(d == 4 for d in degs):
        sols.append(cur)


def canon(sol):
    from itertools import permutations
    best = None
    for p in permutations(range(5)):
        mp = {}
        for m, (a, b) in zip(sol, blocks):
            x, y = sorted((p[a], p[b]))
            mp[(x, y)] = m
        t = tuple(mp[b] for b in blocks)
        if best is None or t < best:
            best = t
    return best


from collections import Counter
orb = Counter(canon(s) for s in sols)
out["nB_profile_case_split"] = {
    "blocks": [list(b) for b in blocks],
    "labelled_solutions": len(sols),
    "orbits_under_relabelling_the_5_positions": [
        {"profile": list(o), "orbit_size": c,
         "multiplicity_2_on": [list(blocks[i]) for i, m in enumerate(o) if m == 2],
         "multiplicity_1_on": [list(blocks[i]) for i, m in enumerate(o) if m == 1]}
        for o, c in sorted(orb.items(), key=lambda kv: -kv[1])],
    "exhaustive": "brute force over all 3^10 multiplicity functions",
}

# --- the A/B search control ----------------------------------------------
ab = {}
for k in (4, 22, 14):
    M = PS.Model(k)
    try:
        st = impose_partner_regular(M)
    except PS.Contradiction:
        ab[str(k)] = {"contradiction_at_seed": True}
        continue
    stt = {"nodes": 0, "fails": 0, "max_depth": 0, "max_settled": 0,
           "solutions": 0, "level_nodes": {}, "level_fails": {}, "witness": None}
    r = PS.search(st, BUDGET, stt, 0, "vertex")
    entry = {"n": M.n, "D": M.nd, "result": r, "nodes": stt["nodes"],
             "fails": stt["fails"], "max_depth": stt["max_depth"],
             "settled_prefix": stt["max_settled"],
             "level_nodes": {str(a): b for a, b in sorted(stt["level_nodes"].items())},
             "satisfiable": {4: "yes (Paley 9)", 22: "yes (BvLS 243)",
                             14: "unknown -- the open case"}[k]}
    if stt["witness"] is not None:
        adjw = PS.reconstruct(M, stt["witness"])
        good, got = srg.check_srg(adjw)
        entry["witness_verified"] = [bool(good), list(got) if good else str(got)]
    ab[str(k)] = entry
out["search_control_ab"] = {"budget_nodes": BUDGET, "branch_order": "vertex",
                            "value_order": "False first", "runs": ab}

dest = os.path.join(ROOT, "problems", "conway-99", "data",
                    "pair-model-evidence.json")
os.makedirs(os.path.dirname(dest), exist_ok=True)
with open(dest, "w") as fh:
    json.dump(out, fh, indent=1, sort_keys=True)
print(f"wrote {dest} ({os.path.getsize(dest)} bytes)")
for k in ("4", "22", "14"):
    e = out["search_control_ab"]["runs"][k]
    print(f"  k={k:>2}: n={e['n']:4d} settled {e['settled_prefix']:3d}/{e['D']:3d} "
          f"nodes={e['nodes']} result={e['result']}  [{e['satisfiable']}]")
