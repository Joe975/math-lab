"""Cross-check the three formulations of Giuga's congruence against each other.

Usage:
    python calibrate_faces.py [--direct-limit N] [--agoh-limit M]

Checks, for every n in range:
  A. power_sum_mod(n)  vs  power_sum_mod_structural(n)   (exact residue, squarefree n)
  B. satisfies_congruence(n) vs satisfies_congruence_structural(n) vs is_prime(n)
  C. agoh_holds(n) vs is_prime(n)                        (n <= agoh limit)
  D. is_counterexample(n) == (composite and satisfies_congruence(n))
"""

import argparse
import sys
import os

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))
from giuga.conditions import (  # noqa: E402
    agoh_holds, bernoulli, factorize, is_counterexample, is_prime,
    power_sum_mod, power_sum_mod_structural, satisfies_congruence,
    satisfies_congruence_structural,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--direct-limit", type=int, default=2000)
    ap.add_argument("--agoh-limit", type=int, default=150)
    a = ap.parse_args()

    bad = 0
    residue_checks = 0
    for n in range(2, a.direct_limit + 1):
        s_direct = power_sum_mod(n)
        s_struct = power_sum_mod_structural(n)
        if s_struct is not None:
            residue_checks += 1
            if s_direct != s_struct:
                print(f"MISMATCH residue n={n}: direct={s_direct} structural={s_struct}")
                bad += 1
        c_direct = satisfies_congruence(n)
        c_struct = satisfies_congruence_structural(n)
        if c_direct != c_struct:
            print(f"MISMATCH boolean n={n}: direct={c_direct} structural={c_struct}")
            bad += 1
        if c_direct != is_prime(n):
            print(f"CONGRUENCE HOLDS FOR COMPOSITE (or fails for prime) n={n}")
            bad += 1
        composite = not is_prime(n)
        if is_counterexample(n) != (composite and c_direct):
            print(f"MISMATCH counterexample-criterion n={n}")
            bad += 1
    print(f"direct vs structural: n = 2..{a.direct_limit}, "
          f"{residue_checks} exact residues compared, all agree" if bad == 0
          else f"{bad} mismatches")

    b = bernoulli(a.agoh_limit)
    agoh_bad = 0
    for n in range(2, a.agoh_limit + 2):
        if n - 1 > a.agoh_limit:
            break
        if agoh_holds(n, b) != is_prime(n):
            print(f"AGOH MISMATCH n={n}: agoh={agoh_holds(n, b)} prime={is_prime(n)}")
            agoh_bad += 1
    print(f"agoh vs primality: n = 2..{a.agoh_limit + 1}, "
          + ("all agree" if agoh_bad == 0 else f"{agoh_bad} mismatches"))

    # spot-check the "no even solution" shortcut the Bernoulli side gives
    evens = [n for n in range(4, a.agoh_limit, 2) if agoh_holds(n, b)]
    print(f"even n in 4..{a.agoh_limit} passing Agoh: {evens}")

    return 1 if (bad or agoh_bad) else 0


if __name__ == "__main__":
    sys.exit(main())
