import sys, os, time
sys.path.insert(0, "/home/user/work-crouzeix/harness/crouzeix")
import ratio
from ratio import Cx, Q

def mat(entries, den):
    return [[Cx(Q(a, den), Q(b, den)) for (a, b) in row] for row in entries]

# dense 3x3, denominator 32 (what we'd round float endpoints to)
A32 = mat([[(7,-3),(19,5),(-11,2)],[(3,8),(-5,-13),(17,1)],[(2,-9),(6,4),(-1,11)]], 32)
p32 = [Cx(Q(5,32),Q(-3,32)), Cx(Q(-7,32),Q(2,32)), Cx(Q(11,32),Q(1,32)), Cx(Q(9,32),Q(-4,32))]

for label, den, A, p in [("den32", 32, A32, p32)]:
    for tol_exp, db in [(6, 2), (9, 3)]:
        t0 = time.time()
        r = ratio.crouzeix_ratio(A, p, dir_bound=db, tol=Q(1, 10**tol_exp))
        dt = time.time() - t0
        print(f"{label} tol=1e-{tol_exp} dirs={db}: {dt:.2f}s ratio=[{float(r['ratio_lo']):.6f},{float(r['ratio_hi']):.6f}] refutes={r['refutes']}")
