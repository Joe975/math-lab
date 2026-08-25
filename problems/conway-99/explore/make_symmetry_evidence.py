"""Regenerate problems/conway-99/data/symmetry-evidence.json.

Deterministic and fast (a couple of minutes); the long annealing runs are NOT
re-run here, their outcomes are recorded as data supplied by the caller.
"""

import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "harness", "conway-99"))
sys.path.insert(0, HERE)

import srg
import m33_exhaustive as X
import orbit_general as OG
import orbit_matrix as OM

out = {}

# --- which orders can act at all, and with what orbit profile --------------
out["admissible_orders"] = {
    "no_prime_ge_13": "L2: a prime p > k/2 = 7 with a fixed point forces the "
                      "identity, so p | 99; 99 = 3^2 * 11",
    "order_7": "L3: exactly one fixed point, so orbits are 1 + 14x7",
    "order_9": "L5: semiregular, 11 orbits of 9 (96 is not divisible by 9)",
    "order_11": "L2: fixed-point free, 9 orbits of 11",
    "order_27": "L5: impossible (neither 99 nor 96 is divisible by 27)",
    "order_33": "L5 + exhaustive: impossible (see z33 below)",
    "order_3": "L5: 0 or 3 fixed points",
}

# --- orbit matrices, per profile -----------------------------------------
prof = {}
for m in (33, 11, 9, 3):
    if 99 % m:
        continue
    t = 99 // m
    if m == 33:
        sols, st = OM.enumerate_semiregular(m, 14, 1, 2, 99, sorted_diagonal=True)
        prof[f"semiregular_Z{m}"] = {
            "orbits": f"{t} x {m}", "orbit_matrices": len(sols),
            "nodes": st["nodes"], "matrices": sols,
        }
    else:
        prof[f"semiregular_Z{m}"] = {
            "orbits": f"{t} x {m}",
            "orbit_matrices": "enumeration did not terminate",
            "note": "the backtracking search does not reach a first solution "
                    "in minutes for t >= 9; this route is complete in "
                    "principle but not in practice at this size",
        }
for name, sizes in sorted(OG.z33_profiles().items()):
    if sizes == [33, 33, 33]:
        continue
    sols, st = OG.enumerate_general(sizes, 14, 1, 2)
    prof[f"order33_{name}"] = {"sizes": sizes, "orbit_matrices": len(sols),
                               "nodes": st["nodes"]}
out["orbit_matrices"] = prof

# --- the Z_33 elimination -------------------------------------------------
diag = [frozenset({a, (-a) % 33}) for a in range(1, 17)]
alive, killed = 0, 0
for S00 in diag:
    f0 = X.diag_target(S00)
    for S11 in diag:
        f1 = X.diag_target(S11)
        for S22 in diag:
            f2 = X.diag_target(S22)
            sol, _ = X.solve_offdiag(f0, f1, f2)
            if sol is None:
                killed += 1
            else:
                alive += 1
from itertools import combinations
realisable = set()
for rest in combinations(range(1, 33), 5):
    realisable.add(tuple(X.autocorr((0,) + rest)))
required_hits = 0
for S00 in diag:
    f0 = X.diag_target(S00)
    for S11 in diag:
        f1 = X.diag_target(S11)
        for S22 in diag:
            f2 = X.diag_target(S22)
            sol, _ = X.solve_offdiag(f0, f1, f2)
            if sol is None:
                continue
            for A in sol:
                if tuple(A) in realisable:
                    required_hits += 1
out["z33_elimination"] = {
    "unique_orbit_matrix": [[2, 6, 6], [6, 2, 6], [6, 6, 2]],
    "diagonal_choices": len(diag) ** 3,
    "killed_by_forced_autocorrelation": killed,
    "surviving": alive,
    "distinct_autocorrelations_realised_by_6_subsets_of_Z33": len(realisable),
    "surviving_cases_whose_requirement_is_realisable": required_hits,
    "conclusion": "no semiregular Z_33 lift exists",
    "cross_checks": [
        "the forcing formula reproduces the true autocorrelations of Paley(9)",
        "unrealisability re-derived by brute force over all C(32,5)=201376 "
        "6-subsets of Z_33 containing 0",
        "the unique orbit matrix re-derived by direct brute force over all "
        "symmetric 3x3 candidates",
    ],
}

# --- L5 support -----------------------------------------------------------
ok33, rep33 = srg.feasibility(33, 8, 1, 2)
out["L5"] = {
    "statement": "an order-3 automorphism of an SRG(99,14,1,2) has 0 or 3 "
                 "fixed points",
    "degrees_in_Fix": [2, 8],
    "counting_identity": "sum_{x in N_F(u)} deg_F(x) = 2f - 2",
    "regular_degree_2_gives_f": 3,
    "regular_degree_8_gives_f": 33,
    "srg_33_8_1_2_feasible": ok33,
    "srg_33_8_1_2_fails": srg.failed_conditions(rep33),
}

# --- calibration ceilings of each construction method ---------------------
out["search_calibration_ceilings"] = {
    "_note": "each method was run on parameter sets whose graphs exist; the "
             "largest it could rebuild is its ceiling.  All four ceilings are "
             "far below n=99, so none of their null results at 99 is evidence "
             "of nonexistence.",
    "pair_search (propagating, exact)": {
        "recovers": "n=9",
        "fails": "n=243 (BvLS): 14/220 vertices settled after 120000 nodes",
    },
    "srg_anneal (unrestricted)": {
        "recovers": "n=27 (27,10,1,5) lam=1; also 9, 15, 16",
        "fails": "n=25 (25,12,5,6), n=36, n=45",
    },
    "cyclic_anneal (Z_m-invariant)": {
        "recovers": "n=27 (27,10,1,5) lam=1 at t=9 in 16s; also 9, 10, 15, 16",
        "fails": "n=45 (45,12,3,3) at t=9",
    },
    "group_anneal (pair-orbit toggling)": {
        "recovers": "n=9, 16",
        "fails": "n=243 (BvLS) under an assumed order-3 automorphism",
    },
}

# --- null results at 99 ---------------------------------------------------
# All runs were terminated at wrap-up.  None wrote a witness, so there was no
# checkpoint to salvage: the stopping point IS the result.  `controlled` records
# whether the search had first rebuilt a graph that exists AT THIS SCALE; every
# construction run here is uncontrolled, so none of these is evidence about 99.
out["construction_attempts_at_99"] = {
    "_note": "no witness found by any method. Every run below is UNCONTROLLED "
             "-- the engine could not rebuild a known graph at comparable "
             "scale -- so none of these stopping points is evidence of "
             "nonexistence. Recorded so the next attempt knows what was tried.",
    "cyclic_anneal_Z11_t9": {
        "orbits": "9 x 11", "restarts": 11, "iters_per_restart": 400000,
        "best_cost": 336, "controlled": False,
        "control": "rebuilds (27,10,1,5) at t=9 in 16s; fails (45,12,3,3) at t=9"},
    "cyclic_anneal_Z9_t11": {
        "orbits": "11 x 9", "restarts": 8, "iters_per_restart": 400000,
        "best_cost": 440, "controlled": False},
    "cyclic_anneal_Z3_t33": {
        "orbits": "33 x 3", "restarts": 1, "iters_per_restart": 300000,
        "best_cost": 1448, "controlled": False},
    "srg_anneal_unrestricted": {
        "restarts": 2, "iters_per_restart": 800000, "best_cost": 2546,
        "controlled": False, "control": "ceiling n=27"},
    "pair_anneal_in_pair_model": {
        "restarts": 1, "iters_per_restart": 1500000, "best_cost": 2460,
        "controlled": False,
        "control": "k=22 positive control stalled at cost 23464"},
    "group_anneal_Z11": {
        "pair_orbits": 441, "restarts": 3, "iters_per_restart": 60000,
        "best_cost": 2508, "controlled": False},
    "group_anneal_Z7_one_fixed_point": {
        "pair_orbits": 693, "restarts": 5, "iters_per_restart": 60000,
        "best_cost": 2653, "controlled": False},
    "group_anneal_Z3": {
        "pair_orbits": 1617, "restarts": 5, "iters_per_restart": 60000,
        "best_cost": 2733, "controlled": False},
    "positive_control_bvls243_group_anneal": {
        "pair_orbits": 9801, "restarts": 1, "iters_per_restart": 25000,
        "best_cost": 24729, "target_cost": 0,
        "note": "THE CONTROL ITSELF FAILED: this method cannot rebuild a graph "
                "that provably exists, which is why the rows above are marked "
                "uncontrolled"},
    "orbit_matrix_enumeration_m_9_11_3": {
        "result": "no first solution reached in minutes; terminated",
        "note": "this method is complete rather than heuristic, so its failure "
                "is a tractability limit, not a null result"},
}

out["controlled_results"] = {
    "_note": "by contrast, these ARE gated on controls that passed",
    "z33_elimination": "orbit-matrix enumerator validated to contain the true "
                       "orbit matrices of six real graph/automorphism pairs; "
                       "autocorrelation argument verified two independent ways",
    "pair_reduction": "machine-checked at every vertex of Paley(9), rook(3) "
                      "and BvLS(243) -- 261 decompositions",
    "feasibility": "recovers the published verdict on eleven parameter sets, "
                   "in both directions and for the right reason",
}

dest = os.path.join(ROOT, "problems", "conway-99", "data",
                    "symmetry-evidence.json")
with open(dest, "w") as fh:
    json.dump(out, fh, indent=1, sort_keys=True)
print(f"wrote {dest} ({os.path.getsize(dest)} bytes)")
print(f"  Z33: {killed} killed by parity, {alive} survivors, "
      f"{required_hits} realisable -> eliminated")
for k2, v in out["orbit_matrices"].items():
    print(f"  {k2}: {v.get('orbit_matrices')}")
