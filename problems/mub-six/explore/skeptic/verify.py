import json, sys, itertools
B = json.load(open(sys.argv[1]))
B = [[[complex(re, im) for re, im in col] for col in basis] for basis in B]
k, d = len(B), len(B[0])
ip = lambda u, v: sum(x.conjugate()*y for x, y in zip(u, v))
od = max(abs(ip(b[i], b[j]) - (i == j)) for b in B for i in range(d) for j in range(d))
L = 0.0; per = {}
for a, b in itertools.combinations(range(k), 2):
    s = sum((abs(ip(B[a][i], B[b][j]))**2 - 1/d)**2 for i in range(d) for j in range(d))
    per[(a, b)] = s; L += s
print("k=%d d=%d  max orthonormality deviation %.3e" % (k, d, od))
print("L recomputed = %.15e  ASD = %.15f" % (L, 1 - L/(k*(k-1)/2*(d-1))))
for p, s in per.items(): print("pair", p, "%.6e" % s)
