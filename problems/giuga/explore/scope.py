"""Quick sizing: how fast does the admissible reciprocal sum grow?"""
from fractions import Fraction
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "..", "harness"))

def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0:2] = b"\x00\x00"
    i = 2
    while i * i <= n:
        if sieve[i]:
            sieve[i*i::i] = bytearray(len(sieve[i*i::i]))
        i += 1
    return [i for i in range(n + 1) if sieve[i]]

LIM = 5_000_000
P = primes_upto(LIM)
ODD = [p for p in P if p > 2]
print(f"pool: odd primes < {LIM}, {len(ODD)} of them")

# baseline: odd primes, no pairwise rule
s = Fraction(0); t = 0
for p in ODD:
    s += Fraction(1, p); t += 1
    if s > 1:
        break
print(f"odd primes, no pairwise rule: {t} primes, largest {p}")

def greedy(skip=()):
    S, s = [], Fraction(0)
    marks = []
    for p in ODD:
        if p in skip:
            continue
        if any((p - 1) % q == 0 for q in S):
            continue
        S.append(p); s += Fraction(1, p)
        if len(S) in (10, 20, 40, 80, 160, 320, 640, 1280):
            marks.append((len(S), p, float(s)))
        if s > 1:
            return S, s, marks
    return None, s, marks

S, s, marks = greedy()
print(f"greedy from 3: reached sum {float(s):.6f} with {len(marks) and ''}"
      f"{'SUCCESS' if S else 'FAILED'}")
for m in marks:
    print(f"    |S|={m[0]:>5}  largest={m[1]:>9}  sum={m[2]:.6f}")
