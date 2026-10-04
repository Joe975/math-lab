/* Rigid relative equilibrium: q(t) = R(t)q(0), lambda = omega^2 = 1.
   Playback rate only maps wall-clock seconds to dimensionless physical time.
   This prescribes the relative equilibrium; it does not integrate perturbations. */
function rotateConfiguration(points, angle) {
  const c = Math.cos(angle), s = Math.sin(angle);
  return points.map(([x, y]) => [c * x - s * y, s * x + c * y]);
}
function advanceRotation(angle, elapsedSeconds, speed, playing) {
  if (!playing) return angle;
  return (angle + Math.max(0, elapsedSeconds) * speed * (2 * Math.PI / 10)) % (2 * Math.PI);
}
if (typeof module !== 'undefined' && module.exports) {
  module.exports = { rotateConfiguration, advanceRotation };
}
