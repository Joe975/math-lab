# 024 — Fresh-session reviewer pass on 020/022/023 (the Gap-1 closure batch)

- **Problem:** Union-closed sets conjecture, `problems/union-closed/PROBLEM.md`
- **Date:** 2026-08-19
- **Mode:** informed
- **Type:** adversarial review (default stance: refute) of 020 + 022 + 023,
  run in a **fresh session** with no shared state with the 020/022/023
  director — the reviewer-level-independence pass those records queued
  (022's stated residual; queue item 1a). Scope as queued: re-run the
  frozen verifier, audit it independently, re-derive the two proof-shaped
  reductions, re-run the CR adversary. Not a rebuild.
- **Outcome in one line:** everything load-bearing survives — all 12
  certificates reproduce byte-for-byte AND land inside a third,
  from-spec-only reimplementation's values; the verifier audit and both
  hand re-derivations check out — but 023's "deterministic" claim for
  `uc_cr_attack.py` is FALSE (string-hash seeding), its committed
  checkpoint is not regenerable, and its floor +0.0742 is
  trajectory-specific (reporting-level; the qualitative finding is
  confirmed and strengthened by three fresh trajectories).
- **Tools:** `explore/uc_gap1_skeptic_gemini.py` re-run unmodified on the
  frozen `data/gap1deep_witnesses.json` (2.2 s, byte-identical to
  `data/gap1deep_skeptic_gemini_out.txt`);
  `explore/uc_reviewer024_thirdpath.py` (new, written for this review from
  `explore/swarm022_skeptic_brief.md` ALONE — exact Fractions for
  coupling/census/marginals/dichotomy, exact 70-step bisection for the
  Plackett root, float log2 only in the final weight/dev products; output
  `data/reviewer024_thirdpath_out.txt`); three re-runs of
  `explore/uc_cr_attack.py` (~45 s each, log
  `data/reviewer024_cr_reruns.log`). Python 3.11.15, no swarm, no
  external workers.
- **Sources:** none.

## Claims attacked

1. **The 12 certificates of 020 as re-established by 022**: six MM_abs
   kills (n = 4, 5 at t ∈ {6/5, 3/2}), three window kills (n = 4 at
   t ∈ {21/20, 16/15, 27/25}), three R₊-purity legs (t = 7/5) — signs,
   in-regime flags, dichotomy counts, and the R₊ z < 1/2 legs.
2. **022's verifier itself** (`uc_gap1_skeptic_gemini.py`): the audit that
   cleared it was done by the same-session director who produced 020, so
   the program's soundness points were re-audited from scratch here.
3. **The two proof-shaped steps**: the weight identity
   dh₂(z_ρ)/dλ = h₂′(z)·(dz/dρ)·ρ·ln 2 with 020's dz/dρ formula, and the
   R₊ reduction underlying the per-history-weighting no-go.
4. **023's Gap-2 adversarial pass**: the claim of a deterministic engine,
   the recorded floors (+0.074…+0.22), and the headline "no in-regime
   violation".

## Refutations found

### R1. 023's determinism claim is false; its checkpoint is not regenerable

`uc_cr_attack.py` (line 89) seeds each seeded-job anneal with
`random.Random(808000 + hash(name) % 9999)` where `name` is a **string**
— and Python randomizes string hashing per interpreter run unless
`PYTHONHASHSEED` is fixed. Consequences, all verified here:

- A default re-run produces a *different* trajectory from the committed
  `data/cr_attack_run.log` (e.g. windowkill n=4 floor +0.1222 vs
  committed +0.2232; global best +0.0697 vs +0.0742).
- Two runs under `PYTHONHASHSEED=0` are value-identical to each other
  (timestamps aside) — confirming the hash seeding is the *only*
  nondeterminism.
- The committed `data/cr_attack.json` was produced under an unrecorded
  hash seed and cannot be regenerated; 023's "deterministic; ~40 s"
  header line and its specific "sharpest floor +0.0742 (4 atoms, n = 5)"
  are therefore trajectory-specific, not reproducible facts.

Severity: **reporting-level, not load-bearing.** 023 is `EVIDENCE` from a
stochastic search; what the status rests on — no in-regime violation
found — is exactly what every fresh trajectory reproduces (see S4). But
the record's test battery ("023's floor instances") should be read as
"floor instances from *some* run of the engine", and queue 1(b)'s
"sharpest known test +0.0742" is superseded: this review's runs found an
in-regime endpoint at **+0.0697** (windowkill_n6, n = 6, 4 atoms, H = 1.87)
— notably on the very seed the committed run wasted in the degenerate
H → 0 corner. Fix for future engines (beyond 023's own H ≥ ε lead):
seed with a stable function of the job name, and record
`PYTHONHASHSEED` in any log that claims determinism.

No other refutation was found. In particular the certificate layer
(020/022) survives every check attempted.

## Claims that survive

### S1. All 12 certificates — now at reviewer-level independence

- **Byte-identical re-run**: `uc_gap1_skeptic_gemini.py` on the frozen
  witnesses reproduces `data/gap1deep_skeptic_gemini_out.txt` exactly
  (self-tests included). The verifier is genuinely deterministic (exact
  rational arithmetic throughout; no hash-order or float dependence).
- **Third-path agreement**: `uc_reviewer024_thirdpath.py`, written for
  this review from the committed spec alone (shares no code with either
  the 013-kit certificates or the gemini verifier; different root
  enclosure, float log2 instead of certified log2), lands inside every
  certified enclosure (tolerance 1e-11; the enclosures have width
  < 1e-14) and matches every metadata field: in-regime flags (max
  marginals 0.257466–0.382669, all exactly < 0.38271), nondegenerate
  table counts (14/15/21/37), `mixed_zero_tables = 0` everywhere (the
  degeneracy dichotomy holds row-by-row in a third census), and
  `all_rplus = True` on all three no-go witnesses. All 12 numerators
  certified negative stand: the kill signs are real, not an artifact of
  either implementation.

### S2. The verifier audit — re-done fresh, no defects

Soundness points checked independently against the program text:
(i) the two-track log2: the invariant v ≥ 1 holds *exactly* on both
tracks (floor(v·2^200)/2^200 ≥ 1 iff v ≥ 1; squaring preserves ≥ 1;
halving only fires at v ≥ 2), so frac_lo ≤ log2(r′) ≤ frac_hi + 2^−D
with no boundary leak — the monotonicity argument 022 sketched is
airtight, not just plausible; (ii) bisection maintains F(lo) ≤ 0 < F(hi)
with a strict initial bracket check, so the enclosure always contains
the root; (iii) the weight enclosure rounds outward at every step —
hp monotone decreasing used in the right direction, dz numerator
decreasing on z < min(x,y), denominator positivity checked loudly, box
product with straddle-safe absolute value; (iv) the R₊ leg uses the
conservative end (z_hi < 1/2); (v) in-regime and dichotomy checks are
exact and strict. The program matches its spec clause-for-clause; run
unmodified.

### S3. Both hand re-derivations

- **Weight identity**: implicit differentiation of the Plackett quadratic
  F(z; ρ) = z(1−x−y+z) − ρ(x−z)(y−z) gives ∂F/∂z = (1−x−y+2z) + ρ(x+y−2z)
  and ∂F/∂ρ = −(x−z)(y−z), hence dz/dρ = (x−z)(y−z) / ((1−x−y+2z) +
  ρ(x+y−2z)) — 020's formula exactly — and the chain rule with ρ = 2^λ
  supplies the ρ·ln 2 factor. The dropped global constant t·ln 2 is one
  positive number per witness multiplying every term, so it cannot move
  any sign (and cancels in the MM_abs ratio). Sound.
- **R₊ reduction**: in-bracket, dz/dρ > 0 (numerator positive on
  z < min(x,y); denominator certified positive), so sign(σ_λ) =
  sign(log₂((1−z)/z)), i.e. σ_λ > 0 ⟺ z < 1/2. Any w ≥ 0 equal to σ_λ
  on R₊ therefore satisfies E[w·dev] = E[σ_λ·dev] = E[|σ_λ|·dev] on any
  measure all of whose nondegenerate histories are certified in R₊ — so
  the three certified-negative all-R₊ witnesses rule out the entire
  compatibility class at once. Sound, with 020/022's boundary intact:
  "compatible" means equals σ on R₊; nothing is claimed about weights
  that differ from σ somewhere on R₊.

### S4. 023's qualitative finding — confirmed and strengthened

Three fresh, independent trajectories of the CR adversary (one default
hash seed, two at `PYTHONHASHSEED=0`), on top of the committed fourth:
**every in-regime, H-bounded endpoint across all four runs has
sup_λ CR > 0** — observed floors +0.0697, +0.0733, +0.0742 (committed),
with per-job floors up to +0.30. This run's trajectories also found
in-regime endpoints on jobs the committed run lost (windowkill_n6
in-regime at +0.0697 instead of degenerating; random2/random3 in-regime
at +0.16/+0.29 instead of "no in-regime endpoint") — more adversarial
coverage, zero violations. Gap 2's (TAX at p) candidate survives its
reviewer pass. The H(μ) → 0 structural note (no μ-uniform floor; any
quantitative TAX must scale with H(μ)) stands unchanged.

### Scope, and what is not claimed

The certificates' validity is environment-independent (exact rational
arithmetic), but the byte-identical re-run was same-repo, Python 3.11;
the third-path check is float-level (it confirms the enclosures bracket
the true values to ~1e-11, it is not itself a certificate). Unchanged
from 020/022/023 and not addressed here: the n = 12 tensor extension and
the "every n" tensor claim (float + SPECULATION), the measure-zero
completion of the λ-integrated bound (informal), Gap 2 beyond the
Sinkhorn branch (block-adaptive and half-mixing branches still
unexercised — 023 lead 2), and any exact certification of CR. With this
pass, 020's headline — every stated form of Gap 1 dead, with certified
witnesses — carries implementation-level AND reviewer-level
independence; the residual 022 named is discharged.

## References

- This repo: attempts 020, 022, 023 (under review), 013 (exact-census
  standard the witnesses derive from), 009/011 (Gap 2 candidate),
  `explore/swarm022_skeptic_brief.md` (the spec this review's third-path
  checker was written from), `docs/SWARM.md` (independence rules).
- No external sources.
