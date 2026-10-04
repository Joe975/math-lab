# Rotating configurations

`rotating-configurations.html` is an inline interactive visualization of the
20 saved catalogue configurations. Build it with:

```powershell
& C:/Repos/math-lab/.venv/Scripts/python.exe C:/Repos/math-lab/problems/central-configurations-8/explore/build_animation.py
```

Coordinates come directly from the direct-certificate centers in catalogue-001,
rounded to ten decimal places for display. All shapes retain the same spatial
scale (unit masses, lambda=1). Motion is the prescribed relative equilibrium
q(t)=R(t)q(0), with playback mapping ten wall-clock seconds to one revolution
at 1x. It is not a numerical integration or a stability experiment. A stationary
central body, when present, stays at the origin. Catalogue numbering matches the
published source; body numbers identify displayed points only.

Validation, 2026-09-06:

- Two deterministic tests execute the actual animation JavaScript, check
  pause/speed time mapping, preservation of all pairwise distances, and the
  acceleration relation a=-q for all 20 displayed configurations (residual <1e-8).
- The generated fragment matches the source template and embedded certificate data.
- Browser inspection: eight bodies, 20 choices, moving phase, exact pause,
  configuration selection, speed display, and reduced-motion initial pause.
- Light desktop (736px) and dark mobile (360px) inspected; no browser errors or
  mobile horizontal overflow. Same absolute spatial scale across selections.

```powershell
& C:/Repos/math-lab/.venv/Scripts/python.exe -m pytest C:/Repos/math-lab/tests/test_cc8_animation.py -q -p no:cacheprovider
```
