"""Compare discovery to the published numerical catalogue, then certify it.

The comparison is numerical isometry matching, not an exact identity theorem.
The catalogue supplies candidate coordinates, never trusted certificate flags.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path

import search as s


def read_published(path):
    section = path.read_text().split('*** COORDINATES ***')[1].split('*** SOLUTION TESTS ***')[0]
    rows = []
    for line in section.splitlines():
        fields = line.split()
        if len(fields) == 19 and fields[0].isdigit():
            rows.append((int(fields[0]), s.np.array([float(x) for x in fields[3:]]).reshape(8, 2)))
    if not rows:
        raise ValueError('no coordinates parsed')
    return rows


def run(source, discovery, output):
    output.mkdir(parents=True, exist_ok=True)
    found = json.loads((discovery/'configurations.json').read_text())
    comparison = []
    direct = []
    for ident, raw in read_published(source):
        q = s.solve(s.gauge(raw))
        if q is None:
            raise ValueError(f'published candidate {ident} did not converge')
        match = [row['id'] for row in found if s.equivalent(q, s.np.array(row['points']))]
        c = s.central.certify(q)
        p = s.verify_polynomial.certify(q)
        direct.append(c)
        for label, cert in [('direct', c), ('polynomial', p)]:
            (output/f'published-{ident:02d}-{label}.json').write_text(json.dumps(cert, indent=2)+'\n')
        comparison.append({'published_id': ident, 'discovered_matches': match,
                           'points': q.tolist(), 'direct_valid': s.central.verify(c),
                           'polynomial_valid': s.verify_polynomial.verify(p)})
        print(f'published {ident}: match {match}', flush=True)
    result = {'source_url': 'https://raw.githubusercontent.com/AlexandruDoicu/Balanced-and-Central-Configurations/main/CentralConfig/CCNMAS8.dat',
              'sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'catalogue_count': len(comparison), 'complete': False,
              'all_pairs_certifiably_distinct': all(s.central.distinct_by_distances(a, b) for a, b in itertools.combinations(direct, 2)),
              'comparison': comparison}
    (output/'comparison.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'comparison'}))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--discovery', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    run(a.source, a.discovery, a.output)
