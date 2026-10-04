"""Independent certificate rerun and rational distance separation audit."""
import importlib.util
import itertools
import json
from fractions import Fraction as F
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / 'problems/central-configurations-8/data'
spec = importlib.util.spec_from_file_location(
    'polynomial', ROOT / 'harness/central-configurations-8/verify_polynomial.py')
polynomial = importlib.util.module_from_spec(spec)
spec.loader.exec_module(polynomial)


def distance_orders(cert):
    """Independent endpoint formula, without importing the direct checker."""
    x = [F(t) for t in cert['center'][:15]]
    radius = F(cert['radius'])
    x.insert(1, F(0))
    bounds = []
    for i, j in itertools.combinations(range(8), 2):
        low = high = F(0)
        for axis in range(2):
            # Gauge y0 is exactly zero; every other coordinate has this radius.
            uncertainty = sum(F(0) if (k, axis) == (0, 1) else radius
                              for k in (i, j))
            diff = abs(x[2*i+axis]-x[2*j+axis])
            low += max(F(0), diff-uncertainty)**2
            high += (diff+uncertainty)**2
        bounds.append((low, high))
    return sorted(b[0] for b in bounds), sorted(b[1] for b in bounds)


def main():
    checked = []
    catalogue = []
    # Freeze this audit to the two batches in attempt 001; later attempts may
    # add other certificates without changing the scope of this record.
    paths = sorted(path for batch in ('run-001', 'catalogue-001')
                   for path in (DATA/batch).glob('*-polynomial.json'))
    for path in paths:
        cert = json.loads(path.read_text())
        valid = cert['n'] == 8 and polynomial.verify(cert)
        checked.append({'path': path.relative_to(ROOT).as_posix(), 'n': cert['n'], 'valid': valid})
        if path.parent.name == 'catalogue-001':
            catalogue.append((path.name, distance_orders(cert)))
    separations = []
    for (name_a, (al, ah)), (name_b, (bl, bh)) in itertools.combinations(catalogue, 2):
        gaps = [(k, max(bl[k]-ah[k], al[k]-bh[k])) for k in range(28)]
        index, gap = max(gaps, key=lambda item: item[1])
        separations.append({'a': name_a, 'b': name_b, 'order_index': index,
                            'positive_gap': str(gap), 'separated': gap > 0})
    result = {'polynomial_count': len(checked), 'all_valid_n8': all(c['valid'] for c in checked),
              'catalogue_count': len(catalogue), 'pair_count': len(separations),
              'all_pairs_separated': all(p['separated'] for p in separations),
              'scope': 'Local roots and inequivalence only; not completeness or novelty.',
              'certificates': checked, 'pairwise_distance_witnesses': separations}
    (DATA / 'verification-summary.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in {'certificates', 'pairwise_distance_witnesses'}}))
    assert len(checked) == 39 and len(catalogue) == 20 and len(separations) == 190
    assert result['all_valid_n8'] and result['all_pairs_separated']


if __name__ == '__main__':
    main()
