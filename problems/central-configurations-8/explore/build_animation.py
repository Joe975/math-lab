"""Build the inline animation from the saved exact-certificate centers.

Rigid rotation, not a numerical integration or a stability demonstration.
Displayed coordinates are rounded to ten decimal places. All configurations
use the same lambda=1 spatial scale and the same angular speed.
"""
from fractions import Fraction
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent/'data'
NAMES = {3: 'Heptagon and center', 4: 'Twisted squares', 5: 'Collinear',
         7: 'Aligned squares', 8: 'Regular octagon'}


def configuration_data():
    rows = []
    for i in range(1, 21):
        cert = json.loads((DATA/'catalogue-001'/f'published-{i:02d}-direct.json').read_text())
        if cert['n'] != 8 or cert['normalization'] != 'unit-masses-lambda-1-y0-0':
            raise ValueError('unexpected certificate convention')
        flat = [round(float(Fraction(v)), 10) for v in cert['center']]
        flat.insert(1, 0.)
        rows.append({'id': i, 'name': NAMES.get(i, f'Configuration {i}'),
                     'points': [flat[k:k+2] for k in range(0, 16, 2)]})
    return rows


def build():
    template = (HERE/'rotation-template.html').read_text(encoding='utf-8')
    html = template.replace('__CONFIGURATION_DATA__', json.dumps(configuration_data()))
    html = html.replace('__ROTATION_MATH__', (HERE/'rotation.js').read_text(encoding='utf-8'))
    target = DATA/'rotating-configurations.html'
    target.write_text(html, encoding='utf-8')
    return target


if __name__ == '__main__':
    print(build())
