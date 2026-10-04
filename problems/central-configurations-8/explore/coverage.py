"""Conservative exact coverage pilot on a published-bound gauge domain.

With unit masses and lambda=1, published bounds give max|q_i|<=14 and
min r_ij>1/1568. Rotate/relabel a farthest body onto x0>=1, y0=0.
This pilot only excludes boxes by necessary geometric conditions; every
remaining box is UNRESOLVED. It must never report a complete classification.
The binary tree records actual axes and rational cuts, not just path lengths.
"""
from fractions import Fraction as Q
import argparse
from collections import Counter
import json
from pathlib import Path


def square(lo, hi):
    return (Q(0) if lo <= 0 <= hi else min(lo*lo, hi*hi), max(lo*lo, hi*hi))


def points(box):
    return [[box[0], (Q(0), Q(0))]]+[[box[i], box[i+1]] for i in range(1, 15, 2)]


def reason(box):
    p = points(box)
    for axis in range(2):
        lo, hi = sum(q[axis][0] for q in p), sum(q[axis][1] for q in p)
        if lo > 0 or hi < 0:
            return 'center-of-mass'
    for q in p[1:]:
        lower = sum(square(*v)[0] for v in q)
        if lower > box[0][1]**2:
            return 'farthest-anchor'
    for i in range(8):
        for j in range(i):
            upper = sum(square(p[i][k][0]-p[j][k][1], p[i][k][1]-p[j][k][0])[1] for k in range(2))
            if upper <= Q(1, 1568)**2:
                return 'minimum-separation'
    return None


def build(depth):
    root = [(Q(1), Q(14))]+[(Q(-14), Q(14))]*14
    counts = Counter()

    def split(box, remaining):
        why = reason(box)
        if why or remaining == 0:
            counts[why or 'UNRESOLVED'] += 1
            return {'leaf': why or 'UNRESOLVED'}
        axis = max(range(15), key=lambda k: box[k][1]-box[k][0])
        cut = sum(box[axis])/2
        left, right = list(box), list(box)
        left[axis] = (box[axis][0], cut)
        right[axis] = (cut, box[axis][1])
        return {'axis': axis, 'cut': str(cut), 'left': split(left, remaining-1), 'right': split(right, remaining-1)}

    tree = split(root, depth)
    return {'version': 1, 'root': [[str(a), str(b)] for a, b in root],
            'max_depth': depth, 'tree': tree, 'leaf_counts': dict(counts),
            'complete_classification': False}


def verify(cert):
    expected_root = [['1', '14']]+[['-14', '14']]*14
    if cert.get('root') != expected_root or cert.get('complete_classification') is not False:
        return False
    counts = Counter()

    def visit(node, box):
        if 'leaf' in node:
            if set(node) != {'leaf'}:
                return False
            why = node['leaf']
            if why != 'UNRESOLVED' and reason(box) != why:
                return False
            counts[why] += 1
            return True
        if set(node) != {'axis', 'cut', 'left', 'right'}:
            return False
        axis = node['axis']
        if type(axis) is not int or not 0 <= axis < 15:
            return False
        cut = Q(node['cut'])
        if not box[axis][0] < cut < box[axis][1]:
            return False
        left, right = list(box), list(box)
        left[axis] = (box[axis][0], cut)
        right[axis] = (cut, box[axis][1])
        return visit(node['left'], left) and visit(node['right'], right)

    try:
        return visit(cert['tree'], [(Q(a), Q(b)) for a, b in cert['root']]) and dict(counts) == cert['leaf_counts']
    except (ValueError, TypeError, KeyError, ZeroDivisionError):
        return False


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--depth', type=int, default=12)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    if not 0 <= a.depth <= 20:
        p.error('depth must be between 0 and 20')
    cert = build(a.depth)
    assert verify(cert)
    a.output.write_text(json.dumps(cert, indent=2)+'\n')
    print(json.dumps(cert['leaf_counts']))
