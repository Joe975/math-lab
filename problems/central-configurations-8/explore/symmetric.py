"""Exact scalar certificate for two aligned concentric squares.

Uniqueness is the analytic monotonicity argument in the attempt record.
This module certifies a bracket and independently tests the reduction against
the unreduced eight-body force computation. No general D4 classification.
"""
from fractions import Fraction as Q
import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'harness/central-configurations-8'))
from central import Interval


def h(t):
    t = Q(t)
    if t <= 1:
        raise ValueError('outer to inner radius must exceed 1')
    return 2*t*t*(3*t*t+1)/((t**3-1)*(t*t-1)**2)


def c4():
    return (Interval(1)+2*Interval(2).sqrt())/4


def bracket(steps=70):
    lo, hi = Q(2), Q(3)
    c = c4()
    if not h(lo) > c.hi or not h(hi) < c.lo:
        raise ValueError('initial bracket fails')
    for _ in range(steps):
        mid = (lo+hi)/2
        v = h(mid)
        if v > c.hi:
            lo = mid
        elif v < c.lo:
            hi = mid
        else:
            raise ValueError('increase square-root precision')
    return {'lo': str(lo), 'hi': str(hi), 'c4': [str(c.lo), str(c.hi)]}


def verify_bracket(cert):
    try:
        lo, hi = Q(cert['lo']), Q(cert['hi'])
        c = c4()
        return 1 < lo < hi and h(lo) > c.hi and h(hi) < c.lo
    except (ValueError, KeyError, ZeroDivisionError):
        return False


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = bracket()
    assert verify_bracket(result)
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result))
