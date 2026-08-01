# 005 — The death law completed: parametric necessity for all (a,b), exact death per member

- **Problem:** billiards-triangles, `problems/billiards-triangles/PROBLEM.md`
- **Date:** 2026-07-31
- **Mode:** informed (read `prior-art.json`, attempts 003 and 004 in full,
  the tier-0 statement, and 001/002 skimmed for conventions; queue item 11
  of `STATUS.md`, the conservative track of the two-track plan)
- **Type:** formalization + computation — parametric proofs of 003's
  per-member inputs (I1–I3 for all integers (a,b); Lemma C and Lemma D for
  all a ≥ 2, b ≤ a), removal of the a ≤ 2b+3 scope restriction from the
  case tree, and certified sufficiency (death = γ_d exactly) for 20 members
- **Tools:** new `explore/plaw_general.py` (formal 4-variable Laurent-ring
  proofs, per-member exact specialization cross-checks, extended-tree
  adversarial scans, Lemma-C sanity grid), new `explore/plaw_suffice.py`
  (certified alive segments into the death corner: exact mod-360-reduced
  interval trig, per-gate 1-D sign certification).  Reused as substrate:
  003's `deathlaw_symbolic.py` (glide data; its selftest ties it to 002's
  independent corridor and 004 confirmed it digit-for-digit),
  `deathlaw_measure.py`, `deathlaw_exact.py` (which itself re-derives
  orbits with 002's `skeptic_orbit.py`), tier-0 `unfold.py` (rational pi
  enclosure; interval sin/cos as cross-reference in the selftest).  004's
  `deathlaw_skeptic.py` deliberately NOT used — it stays reserved for
  skeptics.  All stdlib-only, deterministic; total compute ≈ 15 min.
- **Sources:** none external.  Repo: attempts 003 (the record extended
  here), 004 (its corrections C1, C2 heeded: "alive" = positive corridor
  width throughout, and a zero-width touching corridor at γ_d is not
  excluded — the sup below is proven not attained).

**Conventions** (001/003): W(a,b) = (0 (12)^a (02)^b)^2, angles α at A, β
at B, γ = 180 − α − β obtuse at C, θ = α + β; θ_d(a,b) = 90(a+b)/(a(b+1))
deg, γ_d = 180 − θ_d; "alive at a triangle" = the word's unfolding
corridor has positive width there; death(W) = sup of γ over alive
triangles.  Circumdiameter-1 normalization: A = 0, B = sin(α+β),
C = sin(β)e^{iα}; p(z) = Im(conj(δ)z) is the projection onto the normal of
the glide axis direction δ of the half word u = 0(12)^a(02)^b; m is the
axis offset.  All angle arithmetic below is in degrees unless marked.

## Approach

Queue item 11 in its stated order: (i) sufficiency, (ii) parametric
inputs, (iii) the a > 2b+3 branch.  I did (ii) first because it changes
what (i) and (iii) even mean: with I1–I3 and Lemma C proven for general
(a,b), the necessity theorem stops being a per-member statement, and both
the sufficiency segments and the new-branch members inherit it for free.

The key methodological choice, and why it beats the obvious alternative:
003's Lead 1 proposed proving I1–I3 by carrying out the two-fan geometric
sum by hand (Dirichlet-kernel forms).  That proof would be pages of
sign-sensitive trigonometric algebra, exactly the kind a skeptic distrusts.
Instead I changed the domain so that the general statement becomes a
FINITE exact computation: since the half word composes to

    u = R_0 ∘ (R_1 R_2)^a ∘ (R_0 R_2)^b = R_0 ∘ Rot_A(2aα) ∘ Rot_B(−2bβ),

the parameters a, b enter every relevant quantity only through e^{iaα} and
e^{ibβ}.  Adjoining these as INDEPENDENT formal variables P, Q turns each
identity into one Laurent-polynomial identity in Q(i)[x^±, y^±, P^±, Q^±]
(x = e^{iα}, y = e^{iβ}), checkable once by exact arithmetic and then true
for every integer pair (a,b) and all angles by ring-homomorphism
specialization P → x^a, Q → y^b.  The only hand steps left are the short
derivation of the closed form of u (three lines, written out below) and
the specialization argument; both are cross-checked exactly, member by
member, against 003's independent gate-by-gate symbolic unfolding.

For Lemma C, 003's Lead 1 route (log-derivative of Φ) turned out to close
with one elementary monotonicity lemma; the same lemma also proves Lemma D,
removing 004's "taken as classical" residual.  For sufficiency, instead of
Lead 2's first-order expansion with a quadratic remainder (which needs
all-gate Taylor models), I certify an explicit rational SEGMENT of
triangles ending at the death corner, reducing "alive on the whole
segment" to ~n one-dimensional certified sign problems — the pattern of
`deathlaw_prove.py`, single non-iterated interval evaluations only.

## What was done

### (ii)-A. I1–I3 for ALL integers a, b ≥ 1 — formal proof

**Step 1 (composed-map form of the unfolding; standard, one line).**  If
M_k is the isometry placing the (k+1)-st unfolded copy, and gate k is the
image of the mirror side s_k, then reflecting across the current image
side g_k = M_{k−1}(side s_k) is M_{k−1} R_{s_k} M_{k−1}^{−1}, where R_s is
the reflection across side s of the BASE triangle; composing on the left
with M_{k−1} gives

    M_k = M_{k−1} ∘ R_{s_k},   so   M_k = R_{s_1} ∘ R_{s_2} ∘ … ∘ R_{s_k},

and gate k = M_{k−1}(side s_k).  (004 already re-derived the unfolding in
exactly this form and tied it to the harness corridor; here it is also
verified exactly per member, Step 4.)  Consequently the composed map of
the half word is the word-ordered product u = R_0 (R_1 R_2)^a (R_0 R_2)^b,
and since W = u·u, gate n′+k = u(gate k) — 003's functoriality — is an
identity of the construction.

**Step 2 (pair collapse).**  In the circumdiameter normalization the base
reflections are R_2: z ↦ conj z (side AB = real axis through A = 0), R_1:
z ↦ e^{2iα} conj z (side CA through A = 0 with direction e^{iα}), R_0:
z ↦ B + e^{−2iβ} conj(z − B) (side BC through B = sin θ with direction
e^{i(180−β)}).  Hence R_1∘R_2 = rotation about A by +2α and R_0∘R_2 =
rotation about B by −2β (compose the two formulas; the conjugations
cancel).

**Step 3 (closed form of u).**  Composing u = R_0 ∘ [z ↦ e^{2iaα} z] ∘
[z ↦ B + e^{−2ibβ}(z − B)] and using conj(B) = B:

    u(z) = μ conj(z) + w,
    μ = e^{−2iaα + 2i(b−1)β} = P^{−2} Q^2 y^{−2},
    w = B · (1 + P^{−2} y^{−2} − y^{−2} − μ),      B = sin(α+β),

so δ = P^{−1} Q y^{−1} = e^{i(−aα + (b−1)β)} satisfies δ² = μ, and
τ = w + μ conj(w), m = Im(conj(δ)(w − τ/2))/2.  The word length 2a+2b+1
of u is odd, so u is orientation-reversing; u² is automatically a
translation (u²(z) = |μ|² z + μ conj(w) + w = z + τ).

**Step 4 (formal verification — this is the proof of the identities).**
`plaw_general.py formal` builds these closed forms in the formal ring
Q(i)[x^±, y^±, P^±, Q^±] and checks by exact polynomial subtraction:

  - glide facts: μ = δ² and Im(conj(δ)τ) = 0 (τ parallel to the axis);
  - **I1**: m − p(C) = −[cos(aα) sin α sin(bβ) + sin(aα) sin β cos(α+(b+1)β)]
  - **I2**: p(A₁) − m = cos(aα) sin((b+1)β) sin(α+β)
  - **I3**: p(B) − m = cos((b+1)β) sin(aα) sin(α+β)
  - **I4** (new; needed for sufficiency): with v₂ = R_0(Rot_A(2aα) C) the
    C-image endpoint of gate 2a+2,
    p(v₂) − m = cos(aα) sin α sin(bβ) − sin(aα) sin β cos(α+(b+1)β)
    (note I1 = −(S+T) and I4 = S−T for the same two products S, T);
  - context: D1, D2 (collapsing-gate loci) as in 003.

All eight are the zero polynomial.  Since P → x^a, Q → y^b followed by
x → e^{iα}, y → e^{iβ} is a ring homomorphism, **I1–I4, D1, D2 and the
glide facts hold for every integer a, b ≥ 1 and all real angles.**  This
discharges 003's SPECULATION (i) — and the glide-to-corridor reduction,
whose only member-specific inputs were the glide facts, is now itself
fully parametric.

```
python3 problems/billiards-triangles/explore/plaw_general.py formal --out problems/billiards-triangles/data/plaw_formal.json
```

**Cross-checks of the hand steps 1–3** (the formal check cannot catch a
wrong closed form, so the bridge to the geometry is tested separately):

  - exact, per member: specializing the formal μ, w, δ, τ, m, p(B), p(C),
    p(A₁), p(v₂) at (a,b) reproduces `deathlaw_symbolic`'s independently
    computed ring elements (gate-by-gate reflected-vertex chains — a
    different construction, 004-audited) EXACTLY, for all 15 members of
    003 plus (6,3),(6,1),(8,2),(8,1),(7,1),(11,2),(14,4),(17,9):
    `plaw_general.py specialize` → 23/23 MATCH
    (`data/plaw_specialize.json`);
  - float: closed-form u(z) vs directly composed reflections at 400
    random ((a,b), α, β, z), a ≤ 40 — worst error 1.4e-14
    (`plaw_general.py floatcheck`).

### (ii)-B. Lemma C and Lemma D for general (a,b) — elementary proof

The whole of both lemmas reduces to (radians throughout this subsection):

**Lemma L1.**  For fixed s ∈ (0, π/2], the map c ↦ c·cot(cs) is strictly
decreasing on (0, 1].
*Proof.*  d/dc [c cot(cs)] = cot(cs) − cs/sin²(cs) =
[sin(2cs) − 2cs] / (2 sin²(cs)) < 0, since sin x < x for x > 0 and
cs ∈ (0, π/2] keeps sin(cs) > 0.  ∎

**Lemma C (general).**  For integers a ≥ 2, 1 ≤ b ≤ a, with
v₀ = π/(2(b+1)): H2(v) := sin v sin((π/2−v)/a) cos v₀ −
sin v₀ sin(bv/a) cos v > 0 for all v ∈ (0, v₀), and H2(v₀) = 0.
*Proof.*  Let Φ(v) = tan v · sin((π/2−v)/a) / sin(bv/a) (all three factors
positive on (0, v₀]: v ≤ v₀ ≤ π/4).  Then, using tan v + cot v =
1/(sin v cos v) = 2/sin 2v,

    (log Φ)′(v) = 2/sin(2v) − (1/a) cot((π/2−v)/a) − (b/a) cot(bv/a)
                = [cot v − (b/a) cot((b/a)·v)]
                + [tan v − (1/a) cot((π/2−v)/a)].

The first bracket is ≤ 0 by L1 with s = v, comparing c = b/a ≤ 1 against
c = 1 (equality iff a = b).  The second is < 0 by L1 with s = π/2 − v ∈
(0, π/2), comparing c = 1/a ≤ 1/2 against c = 1 (strict since a ≥ 2),
because cot(π/2 − v) = tan v.  So Φ is strictly decreasing on (0, v₀].
At the endpoint, π/2 − v₀ = b v₀ exactly, so Φ(v₀) = tan v₀; and
H2(v) = cos v cos v₀ sin(bv/a) · [Φ(v) − Φ(v₀)] > 0 on (0, v₀).  ∎

This is precisely the log-derivative route sketched in 003's Lead 1, and
it discharges SPECULATION (ii) — with a strictly larger scope than
conjectured (F < 0 holds on all of (0, π/2), not just (0, v₀]).

**Lemma D (upgrade from "classical, hand-derived").**  D_b(x) =
sin(bx)/sin(x) strictly decreasing on (0, π/(2b)] for b ≥ 2:
(log D_b)′ = b cot(bx) − cot x < 0 by L1 with s = bx, c = 1/b vs c = 1.
D_1 ≡ 1.  ∎  (004 listed Lemma D as a residual risk; it is now proven by
the same two-line mechanism.)

Sanity net (not load-bearing — the proofs above are elementary):
`plaw_general.py lemmac` checks F < 0 and H2 > 0 at 914,500 grid points,
a ≤ 60: zero violations (`data/plaw_lemmac_grid.json`).

### (iii). The a > 2b+3 branch: the restriction is REMOVABLE

**Measurements first** (the queue's order).  `deathlaw_measure.py death
--full` on the unmeasured members, including W(6,3) (b = a−3, inside old
scope but never measured) and the genuinely out-of-scope W(6,1), W(8,2)
(a = 2b+4), W(7,1) (a = 3b+4), W(8,1) (a = 3b+5):

| word | len | γ_d(a,b) | measured death | meas − pred | last-alive apex |
|---|---|---|---|---|---|
| W(6,3) | 38 | 146.25 = 585/4    | 146.2499999999963 | −3.7e-12 | (15, 18.75) |
| W(6,1) | 30 | 127.5 = 255/2     | 127.4999999999947 | −5.3e-12 | (15, 37.5) |
| W(8,2) | 42 | 142.5 = 285/2     | 142.4999999999964 | −3.6e-12 | (11.25, 26.25) |
| W(8,1) | 38 | 129.375 = 1035/8  | 129.3749999999951 | −4.9e-12 | (11.25, 39.375) |
| W(7,1) | 34 | 900/7 = 128.5714… | 128.5714285714234 | −5.1e-12 | (12.857, 38.571) |

```
python3 problems/billiards-triangles/explore/deathlaw_measure.py death --a 6 --b 3 --full --out problems/billiards-triangles/data/plaw_measure_W63.json   # + (6,1),(8,2),(8,1),(7,1)
```

Every member dies at γ_d, at the predicted corner (90/a, 90(a−1)/(a(b+1)))
— the law does NOT change past a = 2b+3.  That prompted a re-derivation of
the case tree, which closes it:

**Extended case tree (hand; replaces 003's A4/A4′/A5/A5′, keeps
everything else verbatim).**  Standing notation of 003's proof: positive
width forces (Case I) N1, N2, N3 > 0 or (Case II) N1, N2, N3 < 0, where
N1 = m − p(C) = −G, N2 = p(A₁) − m, N3 = p(B) − m, via I1–I3 all three
are the stated products, and sin θ > 0 always.  Branch on the residue
r = aα mod 360 and k = ⌊aα/360⌋ ≥ 0; residues at multiples of 90 force
N2 = 0 or N3 = 0 (excluded), likewise (b+1)β at multiples of 90.

*Case II.*
  - r ∈ (0,90) ∪ (270,360): cos(aα) > 0, so N2 < 0 forces
    sin((b+1)β) < 0, so (b+1)β > 180, so θ > β > 180/(b+1) ≥ θ_d
    (⟺ a ≥ b).  **The (270,360) half is the fix of A4′: 003 used
    θ > α > 270/a ≥ θ_d there, which is what required a ≤ 2b+3; the sign
    of cos(aα), already positive on this branch, makes the A1′ argument
    apply verbatim and never touches α.**
  - r ∈ (90,180): N3 < 0 forces cos((b+1)β) < 0, so (b+1)β > 90; and
    aα > 90 gives α > 90/a; θ > 90/a + 90/(b+1) > θ_d (slack
    90/(a(b+1))).  Works for every k (α only enters via α > 90/a).
  - r ∈ (180,270): N2 < 0 and N3 < 0 force sin((b+1)β) > 0,
    cos((b+1)β) > 0, so (b+1)β mod 360 ∈ (0,90).  If (b+1)β > 360:
    θ > 360/(b+1) ≥ θ_d (⟺ 3a ≥ b).  Else (b+1)β ∈ (0,90), so
    sin(bβ) > 0; both cos(aα) < 0 and sin(aα) < 0, so G > 0 (Case II)
    forces cos(α+(b+1)β) < 0, and α + (b+1)β ∈ (0, 270) (α < 180,
    (b+1)β < 90) makes that α + (b+1)β > 90; with α > 180/a:
    θ > (bα+90)/(b+1) > (180b/a + 90)/(b+1) > θ_d (⟺ 90b > 0).
    Works for every k.  (This is 003's A3′ unchanged; restated to show no
    hidden k-dependence.)

*Case I.*
  - r ∈ (90,180): cos(aα) < 0, N2 > 0 forces sin((b+1)β) < 0,
    (b+1)β > 180, θ > 180/(b+1) ≥ θ_d.  (003's A2; k-free.)
  - r ∈ (180,270): as 003's A3: sin((b+1)β) < 0 and cos((b+1)β) < 0
    together give (b+1)β > 180, done.  (k-free.)
  - r ∈ (270,360): cos(aα) > 0 and sin(aα) < 0.  N3 > 0 forces
    cos((b+1)β) < 0 and N2 > 0 forces sin((b+1)β) > 0, so
    (b+1)β mod 360 ∈ (90,180), hence (b+1)β > 90; with α > 270/a:
    θ > 270/a + 90/(b+1) ≥ θ_d ⟺ 3(b+1) + a ≥ a + b ⟺ 2b + 3 ≥ 0.
    **This is the fix of A4: the forced (b+1)β > 90, which 003 did not
    use on this branch, converts the deficit 270/a < θ_d (a > 2b+3) into
    an unconditional surplus.**
  - r ∈ (0,90), aα = 360k + 90 − v with v ∈ (0,90): cos(aα) = sin v > 0,
    sin(aα) = cos v > 0.  N2, N3 > 0 force (b+1)β mod 360 ∈ (0,90); if
    (b+1)β > 360 then θ > 360/(b+1) ≥ θ_d (⟺ 3a ≥ b); so assume
    (b+1)β ∈ (0,90), hence sin β, sin(bβ) > 0.  As in 003,
    G = sin v sin α sin(bβ) − cos v sin β sin w with w = α + (b+1)β − 90 ∈
    (−90, 180), and the threshold algebra generalizes to

        a·w − b·v = a(b+1)(θ − θ_d) − 360·k·b.

    Suppose θ ≤ θ_d.
      - **k ≥ 1 (the genuinely new branch):** then aw ≤ bv − 360kb < 0
        (v < 90 < 360), so w < 0, so sin w ≤ 0 (w ∈ (−90, 0]) and
        G ≥ sin v sin α sin(bβ) > 0 — contradicting N1 > 0 (G < 0).  No
        Lemma C, no Lemma D, no bound on a.
      - k = 0: 003's main-case argument verbatim (w ≤ 0 sub-case;
        Regime 2 via Lemma D, using α < 90/a ≤ 45, i.e. a ≥ 2; Regime 1
        via Lemmas C + D, using w ≤ bv/a ≤ v, i.e. b ≤ a) — both lemmas
        now proven in general above.
    Hence θ > θ_d, strictly (each escape is strict; the k = 0 core is
    strict by strict Lemmas C/D as in 003/004).

**THEOREM (parametric necessity).**  For all integers a ≥ b ≥ 1 with
a ≥ 2: if the corridor of W(a,b) has positive width at a triangle with
angles (α, β), then θ > θ_d(a,b), i.e. γ < γ_d(a,b).  — The a ≤ 2b+3
hypothesis of 003 is gone; no member-specific input remains.

Adversarial support for the new branches (the k ≥ 1 and r ∈ (270,360)
regions are reachable): seeded sign-pattern scans at θ ≤ θ_d, boundary-
and corner-accumulated, 13 members with a up to 3b+9 — (6,1), (7,1),
(8,1), (9,1), (12,1), (8,2), (11,2), (14,2), (10,3), (17,3), (6,3),
(16,4), (12,12) — ~120k points each, ~1.56M total: **zero Case I or
Case II hits**; the positive control just above θ_d fires for every
member (`plaw_general.py casetree --a A --b B`,
`data/plaw_casetree_W*_*.json`).

### (i). Sufficiency: death(W(a,b)) = γ_d(a,b) EXACTLY, certified for 20 members

**Reduction (converse direction, new — 003 only needed ⟹).**  With the
glide facts (now parametric), the full-word corridor is the intersection
[L, H] ∩ [2m − H, 2m − L] over the n′ half-word gate projection intervals
I_k = [lo_k, hi_k], L = max lo_k, H = min hi_k (gate n′+k = u(gate k) and
p ∘ u = 2m − p).  If lo_k < m < hi_k STRICTLY for every k, then L < m < H,
so min(H, 2m−L) > m > max(L, 2m−H): positive width.  Combined with 003's
forward direction: positive width ⟺ m strictly inside every half-gate
interval (τ ≠ 0 assumed; it is certified along the segments below).  So
"alive on a segment" is exactly 2n′ one-dimensional sign conditions.

**The segment.**  Corner (α_c, β_c) = (90/a, 90(a−1)/(a(b+1))) — where
every measurement in 003 and above found the last-alive apexes.  Take

    α(t) = α_c − t,  β(t) = β_c + ρt,  t ∈ (0, t₀],  ρ = 2,

so θ(t) = θ_d + (ρ−1)t and γ(t) = γ_d − (ρ−1)t ↑ γ_d as t ↓ 0.  At the
corner exactly three of the distinct difference functions
D = p(endpoint) − m vanish (all members inspected; everything else is
bounded away):

  - p(C) − m (gates 1, 2) = −N1; by I1, on the segment
    N1(t) = cos(at) sin β(t) sin(ŵt) − sin(at) sin α(t) sin(bβ(t)),
    ŵ = (b+1)ρ − 1 (uses aα_c = 90 and α_c + (b+1)β_c = 90 exactly);
  - p(A₁) − m = N2 (the A-fan pivot, one endpoint of gates 2…2a+2); by
    I2, N2(t) = sin(at) sin((b+1)β(t)) sin θ(t) — every factor positive
    on (0, t₀] by RATIONAL range checks alone (at ≤ 45,
    (b+1)β(t) ∈ (90−90/a, 180), θ(t) ∈ (θ_d, 90));
  - p(v₂) − m (the C-image endpoint of gate 2a+2) = N4; by the new I4,
    N4(t) = sin(at) sin α(t) sin(bβ(t)) + cos(at) sin β(t) sin(ŵt) —
    BOTH products positive, again rational range checks only.

N1 is the one genuinely two-sided margin: N1(0) = 0 exactly and
`plaw_suffice.py` certifies N1′ > 0 on [0, t₁] (interval enclosure of the
explicit derivative; ρ = 2 is chosen so the certified
N1′(0) = ŵ sin β_c − a sin α_c sin(bβ_c) > 0) and N1 > 0 on [t₁, t₀] by
adaptive bisection.  Every other distinct difference (5–102 bisection
leaves per member) plus the nondegeneracy function Re(conj(δ)τ) is
certified to have constant sign on [0, t₀] by adaptive interval bisection,
and each gate's two endpoint differences are checked to straddle m.  The
three identities I1, I2, I4 are re-verified per member by exact ring
subtraction before use.  Certified trig: exact mod-360 degree reduction,
dyadic pre-rounding, alternating-series remainders, the harness's rational
pi enclosure; every evaluation is a single non-iterated interval
evaluation (the 2026-07-28 affine-forms lesson does not bite);
cross-checked against `unfold.sin_iv/cos_iv` and floats in `--selftest`.

**Result.**  For all 20 members — 003's fifteen plus (6,3), (6,1), (8,2),
(8,1), (7,1) — the certification succeeds with ρ = 2, t₀ = 1/4:

    for every t ∈ (0, 1/4]:  W(a,b) is alive at (α_c − t, β_c + 2t),
    with γ(t) = γ_d − t sweeping [γ_d − 1/4, γ_d).

Hence death(W(a,b)) ≥ γ_d; with the necessity theorem, and since alive
γ < γ_d always, **death(W(a,b)) = γ_d(a,b) exactly, and the sup is not
attained** (004's C2: at γ_d itself the corridor closes to at most a
touching line).  Per member this replaces 003's brackets
[γ_d − 1e-6, γ_d] by the exact value.

```
python3 problems/billiards-triangles/explore/plaw_suffice.py selftest
python3 problems/billiards-triangles/explore/plaw_suffice.py all --out problems/billiards-triangles/data/plaw_suffice_all.json   # ~4 min, ALL MEMBERS CERTIFIED
```

Float cross-check at t₀/2 and t₀/97 via 002's independent corridor:
positive width at all 40 points.  Independent exact-orbit confirmations
for the five new members (rational apex, exact positive corridor, orbit
re-derived by 002's exact simulator, certified γ enclosure):

```
python3 problems/billiards-triangles/explore/deathlaw_exact.py alive --a 6 --b 3 --dg 1e-6 --out problems/billiards-triangles/data/plaw_exact_alive_W63.json   # + (6,1),(8,2),(8,1),(7,1)
```

all five: exact width > 0, orbit verified, γ certified within 1e-6 of γ_d.

## Outcome

- **VERIFIED (formal, all integers a, b ≥ 1):** the identities I1, I2,
  I3, I4, D1, D2 and the glide facts (μ = δ², τ ∥ axis) hold for every
  member of the family and all angles — proven by exact zero-polynomial
  checks in Q(i)[x^±, y^±, P^±, Q^±] plus the three-line closed-form
  derivation of u above; the derivation is exactly matched against the
  independent symbolic unfolding for 23 members and float-matched at 400
  random (a,b) up to a = 40.  003's SPECULATION (i) is discharged.
- **VERIFIED (proof, all integers a ≥ 2, 1 ≤ b ≤ a):** Lemma C, and
  Lemma D for all b ≥ 1, both via Lemma L1 (c·cot(cs) decreasing);
  003's SPECULATION (ii) is discharged and 004's Lemma-D residual risk
  is retired.  (a = b = 1 stays excluded: H2 ≡ 0 there.)
- **VERIFIED (theorem, parametric):** for all integers a ≥ b ≥ 1, a ≥ 2,
  the corridor of W(a,b) has no positive width at any triangle with
  γ ≥ γ_d(a,b) = 180 − 90(a+b)/(a(b+1)).  The former scope condition
  a ≤ 2b+3 is removed by sharpening two branches (the forced sign of
  (b+1)β) and closing the aα > 360 residue-(0,90) branch (θ ≤ θ_d forces
  w < 0 there).  The k = 0 core of the main case is 003's argument
  verbatim, now resting on the general lemmas.
- **VERIFIED (exact death, 20 members):** death(W(a,b)) = γ_d(a,b)
  EXACTLY — sup not attained — for (2,1), (2,2), (3,2), (3,3), (4,2),
  (4,3), (4,4), (5,3), (5,4), (5,5), (6,6), (7,7), (8,8), (10,9),
  (12,12), (6,3), (6,1), (8,2), (8,1), (7,1): certified alive segments
  (α, β) = (90/a − t, 90(a−1)/(a(b+1)) + 2t), t ∈ (0, 1/4], all corridor
  conditions certified strict throughout, γ(t) filling [γ_d − 1/4, γ_d).
- **EVIDENCE:** the five new float death measurements (all equal to γ_d
  to ~5e-12, argmax at the predicted corner; DEAD verdicts sample-bounded
  as always); the 1.56M-point adversarial sign scans supporting the
  extended tree (zero hits, positive controls fire).
- **NOT claimed:** parametric sufficiency (death = γ_d for ALL (a,b) at
  once) — the corner-segment construction is certified per member only;
  see the obstruction below.  Nothing about words outside W(a,b), about
  unstable orbits, about the birth edge, or about triangles at γ = γ_d
  itself beyond "no positive width" (a zero-width touching corridor at
  γ_d remains possible, per 004's C2 — that is what "sup not attained"
  means here).  The measured deaths remain float statements; the exact
  content is carried by the theorem + segments.  a = b = 1 excluded
  throughout.

## Why it failed / what survived

Nothing on the queue item's list failed; the honest ledger of what is
still open and why:

1. **Parametric sufficiency is structured but unfinished.**  The
   certified segments use, per member, the signs of ~n′ = 2a+2b+1 generic
   gate margins along [0, 1/4].  The uniformity is striking — every
   member certifies with the SAME ρ = 2 and t₀ = 1/4 — and the three
   corner-vanishing margins are now closed-form for all (a,b) (I1, I2,
   I4), so exactly one step is missing: a parametric positive lower bound
   on the OTHER gate margins along the segment.  These are fan
   projections: by Step 3 the prefix maps are R_0 Rot_A(2kα) and
   R_0 Rot_A(2aα) Rot_B(−2jβ), so each margin is a 3–5-term trig
   polynomial with fan index k (or j) entering linearly in the angles — an
   explicit two-parameter family of 1-D inequalities (Dirichlet-kernel
   style), not an unstructured n-gate mess.  That is the precise residue
   of Lead 2; per-member certification sidestepped it, and I stopped
   there rather than start a second hand-proof campaign in one cycle.
2. **The corner is 3-fold degenerate, not 2-fold.**  003's "only gates 1
   and 2 bind near death" is correct for the necessity direction, but at
   the death corner itself a THIRD margin vanishes: the gate-(2a+2)
   C-image endpoint (identity I4; N1 = −(S+T), N4 = S−T).  Anyone
   attempting Lead 2's Taylor route without noticing this would have
   certified a spurious first-order system.  Found because the segment
   certifier failed loudly on it at every t₀.
3. **A method lesson worth the ledger:** "prove it for all (a,b)" became
   a finite computation the moment the parameters were moved into the
   exponent lattice (P = e^{iaα}, Q = e^{ibβ} as formal variables).  The
   two-fan geometric sum that 003 expected to do by hand was never
   needed; the three-factor composition u = R_0 Rot_A(2aα) Rot_B(−2bβ)
   is the entire structural content of the family.  The same trick should
   apply to ANY one-parameter word family whose repeated blocks pivot on
   a common vertex — e.g. the design-track family hunt (queue item 12).
4. **What survived for reuse:** the formal-ring engine and closed forms
   (`plaw_general.py`); the certified-segment machinery with exact
   mod-360 trig (`plaw_suffice.py` — the trig layer is a drop-in upgrade
   for any future 1-D certification, ~40× tighter than needed here); the
   L1 lemma (one line, kills both Lemma C and Lemma D and plausibly any
   future cot-comparison); the extended case tree; identity I4; the
   converse reduction (m strictly inside all half gates ⟺ alive), which
   halves every future alive-certification for doubled words.

The skeptic should attack, in order: (a) the closed-form derivation
Steps 1–3 (the formal check is only as good as the bridge to geometry;
the 23-member exact specialization and 400-point float check are the
defense — re-derive u independently); (b) the extended case tree,
especially the k ≥ 1 threshold algebra aw − bv = a(b+1)(θ−θ_d) − 360kb
and the claim w ∈ (−90, 0] ⟹ sin w ≤ 0 given the branch's constraints;
(c) the segment certifier's trig layer (mod-360 reduction, dyadic
rounding, Lipschitz widening) and its N1′ derivative formula —
re-differentiate by hand; (d) the converse reduction's τ ≠ 0 bookkeeping;
(e) whether "alive" as certified here (positive corridor width) matches
001's death-measurement criterion on the mirror half for the b ≤ a−2
members (004 §5 says yes; re-confirm for the five new members).

## Leads generated

1. **Parametric sufficiency (finish the death law for all (a,b)).**
   Compute the generic gate margins along the universal segment
   (90/a − t, 90(a−1)/(a(b+1)) + 2t) in closed form via the prefix maps
   R_0 Rot_A(2kα), k ≤ a, and R_0 Rot_A(2aα) Rot_B(−2jβ), j ≤ b (each a
   short trig polynomial with k, j linear in the angles), and prove them
   ≥ c(a,b) > 0 on [0, 1/4].  Definite outcome: either death(W(a,b)) =
   γ_d for ALL a ≥ b ≥ 1, a ≥ 2, or a member family where some fan margin
   pinches — both publishable-grade facts for this repo.
2. **The 003+005 theorem is now design-grade for queue item 12.**  The
   death corner, the exact three binding margins (I1/I2/I4), and the
   uniform ρ = 2, t₀ = 1/4 segment geometry quantify the alive wedge at
   death for every member; the design track can select (a,b) with γ_d
   above any target and KNOWS its last-life window.  Cross-check against
   the design track's independent findings before merging conclusions.
3. **Birth angles by the same formal method.**  003's Lead 5, now with a
   concrete route: the birth edge should be a different pair (or triple —
   cf. obstruction 2) of binding gates; compute which margins vanish at
   the measured birth apexes (e.g. birth(W(4,3)) ≈ 135.0486), extract
   closed forms in the 4-variable ring, and derive the birth law.  If it,
   too, is algebraic, the family's alive windows become fully exact.
4. **W(1,1) and the a = 1 column.**  The parametric theorem needs a ≥ 2;
   W(1,b) for b = 1 has H2 ≡ 0 and θ_d(1,1) = 90.  Measure W(1,1)'s
   corridor on the obtuse side and settle the degenerate column (likely
   dead everywhere obtuse — then the theorem's a ≥ 2 is not a gap but a
   fact; either way one measurement run decides).
5. **Port the P,Q-formalization to other word families.**  Any family
   (w_0 B_1^a B_2^b)^2 with B_1, B_2 vertex-pivoting blocks has
   u = (isometry) ∘ Rot^a ∘ Rot^b and hence formal-ring identities in the
   same four variables.  Concretely: take the design track's best
   candidate family and attempt the same I1–I4 extraction; success gives
   its death law with the same proof skeleton.

## References

- `problems/billiards-triangles/attempts/003-death-angle-laws.md` (the
  record extended here: theorem, I1–I3, Lemmas C/D, case tree, Leads 1–3)
  and `004-skeptic-review-of-003.md` (corrections C1/C2; the composed-map
  view of the unfolding used in Step 1; both read in full — mode is
  informed).  Conventions from `001-word-census-coverage-map.md` and
  `002-skeptic-review-of-001.md` (skimmed).
- Tier-0: `harness/billiards-triangles/unfold.py` (rational pi enclosure;
  interval sin/cos as selftest cross-reference only).
- New code: `explore/plaw_general.py`, `explore/plaw_suffice.py`.
  New data: `data/plaw_formal.json`, `data/plaw_specialize.json`,
  `data/plaw_lemmac_grid.json`, `data/plaw_casetree_W*_*.json` (13),
  `data/plaw_measure_W{63,61,82,81,71}.json`,
  `data/plaw_suffice_all.json`, `data/plaw_exact_alive_W{63,61,82,81,71}.json`.
- Reused (003/002): `explore/deathlaw_symbolic.py`,
  `explore/deathlaw_measure.py`, `explore/deathlaw_exact.py`,
  `explore/skeptic_family.py`, `explore/skeptic_orbit.py`.
- No external papers consulted.
