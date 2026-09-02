# Literature check: attainability rule (Claim A) and laminate structures (Claim B)

Method: full-text PDF extraction (pypdf) of the freely-available preprints below,
grepped and read directly — not memory. Two papers were fully accessible as PDF;
the primary Milton (1981) paper itself was **not** reached (paywalled, no free
preprint found) — its content below is reconstructed only via how later papers
quote/restate it, which is explicitly flagged.

Sources actually read in full text:
- N. Albin, A. Cherkaev, V. Nesi, "Multiphase laminates of extremal effective
  conductivity in two dimensions," J. Mech. Phys. Solids 55 (2007) 1513-1553.
  Free PDF: https://www.math.k-state.edu/~albin/pubs/albin_cherkaev_nesi_2007_mle.pdf
  (cited below as **ACN2007**)
- A. Cherkaev, "Bounds for effective properties of multimaterial two-dimensional
  conducting composites, and fields in optimal composites" (preprint of the
  Mechanics of Materials 41 (2009) 411-433 paper), free PDF:
  http://www.math.utah.edu/~cherk/publ/newbounds7.pdf (cited below as
  **Cherkaev2009-preprint**; I could not confirm this preprint is byte-identical
  to the published journal version, only that title/author/abstract match and it
  cites itself as forthcoming in Mechanics of Materials — flagged as uncertain
  where it matters)

Sources I could **not** access (paywalled, no free copy found by search):
- Milton (1981), Appl. Phys. A 26, 125-130 — original paper
- Milton (2002), *The Theory of Composites* (book) — could not reach the
  attainability section
- Nesi (1995), Proc. Roy. Soc. Edinburgh A 125, 1219-1239 — original paper
  (only abstract/metadata pages found, e.g. Cambridge Core, Semantic Scholar)
- Cherkaev-Zhang (2011), Int. J. Solids Struct. 48(20), 2800-2813 — the WebFetch
  tool returned "cannot decode PDF" for the arXiv mirror (arXiv:1009.3060); a
  local pypdf re-extraction was not attempted for this one given time budget, so
  treat Claim B question (2) re: this paper as **not checked**, not "checked and
  absent."

---

## Claim A (attainability rule)

**VERDICT: partially confirmed, in a different form than stated.**

I could not reach Milton (1981) itself, so I cannot quote it directly — everything
here is via secondary restatement, explicitly labelled.

ACN2007 (p.10-11, "Milton's structure") describes the construction exactly as the
brief describes it — split phase 1 into two parts, coat phase 2 with one part and
phase 3 with the other so the two coated-circle composites have the same
conductivity, then mix — and gives the resulting attainability condition **as a
volume-fraction threshold on m1**, not as "bound value vs sigma2":

> [T, ACN2007 p.11] "This construction requires that (16) m1 ≥ 2Θ(1 − m2) where Θ
> is a constant defined below in (17)... Similar structures are optimal for the
> opposite bound (13) with K3 taking the role of the 'coating', and K1 and K2 the
> inclusions."
>
> [T, ACN2007 p.11, eq. 17] "Θ = k1(k3 − k2) / [(k2 + k1)(k3 − k1)] ≤ 1/2"

Cherkaev2009-preprint independently confirms this is attributed to Milton 1981,
in the same threshold-on-m1 form, not a bound-vs-sigma2 form:

> [T, Cherkaev2009-preprint p.4] "The first optimal three-material structure was
> found by Milton [29] who considered two kinds of Hashin-Shtrikman coated
> circles [18], mixed together. The structures realize the Hashin-Shtrikman
> bound (a.k.a the isotropic translation bound) in a region of parameters where
> the volume fraction m1 of the best material κ1 is larger than a threshold
> value."

(Ref [29] in Cherkaev2009-preprint's bibliography is listed as "G. W. Milton.
Concerning bounds on the transport and optical properties of a two-component
composite material" — note "two-component," not "multicomponent," and "optical"
not "mechanical." This looks like a bibliographic slip in Cherkaev's reference
list, conflating Milton's 1981 multicomponent paper (Appl. Phys. A 26) with his
different 1981 two-component paper in J. Appl. Phys. 52. ACN2007's bibliography
correctly cites "Concerning bounds on the transport and mechanical properties of
multicomponent composite materials, Appl. Phys., A 26:125-130" — matching your
brief exactly — so I'm confident the *substance* both papers attribute to Milton
is about the multicomponent paper, but flag the citation mismatch in case it
matters for exact sourcing.)

**So: neither secondary source states Milton's condition as "sigma_HS_lower <=
sigma2."** Both state it as a volume-fraction threshold m1 ≥ 2Θ(1−m2) (lower
bound) / the mirror condition with K3 as coating (upper bound). Whether that
threshold condition is algebraically equivalent to "sigma_HS_lower <= sigma2" is
a derivation I did not attempt — it's plausible (both are conditions carving out
the same "phase-1-is-scarce" region) but I have not verified the equivalence, and
neither paper phrases it that way, so treat the equivalence as **unconfirmed,
not disproved**. This needs an independent algebraic check before being recorded
as literature support for the exact form in Claim A, and any attempt record
should cite the m1-threshold form as what's actually published, with the
bound-vs-sigma2 form marked as a to-be-verified reformulation.

I could not check Milton's book (2002) directly (no accessible text found).

---

## Claim B (structures in the non-attainable region)

### (1) Best published bound at the two numeric points

**VERDICT: at both requested points, the "best known bound" — per Cherkaev2009's
own formula — turns out to equal the plain Hashin-Shtrikman lower bound, not a
strictly tighter value.** This is a numeric finding worth flagging back, since it
changes what "close to the best known bound" means for point (ii) of the claim.

Transcribing Cherkaev2009-preprint's explicit three-material theorem:

> [T, Cherkaev2009-preprint p.36, Theorem 7.1] "The effective conductivity k∗ of
> a two-dimensional isotropic composite of three isotropic materials with
> conductivities k1 < k2 < k3 taken in the fractions m1, m2 and m3, m1+m2+m3=1,
> is bounded from below by the bound kL = B(m1,m2): k∗ ≥ B(m1,m2) where
> B(m1,m2) = B1 if m11 ≤ m1 ≤ 1; B2 if m12 ≤ m1 ≤ m11; B3 if 0 ≤ m1 ≤ m12."

with (p.35-36, eqs. 7.5-7.18, [T]):

```
m11 = 2*sqrt(m2)*(1-sqrt(m2))*k1*(k3-k2) / [(k3-k1)*(k1+k2)]
m12 = (1-sqrt(m2))/(4*k2*(k3-k1)) * Z0
Z0  = 2*k2*(k3-k1) + sqrt(m2)*(k1+k2)*(2*k3-k1-k2) - sqrt(Z2)
Z2  = 4*k2^2*(k3-k1)^2 + 4*sqrt(m2)*Z3 + m2*Z4
Z3  = k2*(k1-k2)*(k1-k3)*(k1-k2+2*k3)
Z4  = (k1-k2)^2*(k1^2+6*k1*k2-4*k1*k3-4*k2*k3+4*k3^2+k2^2)

B1 = -k1 + [ m1/(2*k1) + m2/(k1+k2) + m3/(k1+k3) ]^-1     (= Hashin-Shtrikman lower bound, verbatim)
B2 = k2 + (1-sqrt(m2))^2 * Z5/Z6
     Z5 = m1*k1^2 - m1*k2^2 + 2*m3*k1*(k3-k2)
     Z6 = [(1-sqrt(m2))^2 + (1-m1-sqrt(m2))^2]*k1 + m1*(1-sqrt(m2))^2*k2 + m1*m3*k3
B3 = -k2 + [ m2/(2*k2) + Z7 ]^-1
     Z7 = [(k1-k2)*m1^2 + (2*k1-k2+k3)*m1*m3 + 2*k1*m3^2] / [(k1^2-k2^2)*m1 + 2*k1*(k2+k3)*m3]
```

I evaluated this in Python at both requested points (k1,k2,k3)=(1,2,5), and cross-checked
by directly maximizing the underlying one-parameter family B(t) = -t + 1/H1(t)
(eq. 7.1-7.2) over t on a 200,001-point grid, rather than trusting only the
closed-form region split — the two methods agree to 1e-9:

| f=(m1,m2,m3) | region (by m1 vs m11,m12) | Cherkaev bound B | HS lower bound | difference |
|---|---|---|---|---|
| (1/8, 1/8, 3/4) | m11=0.114277, m12=0.087290 → **m1=0.125 is in B1** | 3.363636 | 3.363636 | **0** |
| (1/4, 1/4, 1/2) | m11=0.125000, m12=0.096045 → **m1=0.25 is in B1** | 2.428571 | 2.428571 | **0** |

Both points land in region B1 — where, per the theorem, "B1 ... degenerates into
Hashin-Shtrikman bound... This happens when m1 ≥ m11" (p.36). So **Cherkaev's
2009 bound does not improve on Hashin-Shtrikman at either of your two requested
points** — it coincides with HS exactly. If your rank-3 and rank-4 laminate
results were compared against Cherkaev's bound expecting a value strictly below
HS, that comparison target is wrong: at these (f, sigma), Cherkaev's own theorem
says the best known bound from that paper is just HS itself. Your rank-4 result
(5.8e-4 short of HS lower at f=(1/8,1/8,3/4)) is then evidence about approaching
the *unimproved* HS lower bound in a regime the paper's own machinery does not
tighten — worth double-checking against the harness's HS computation
independently, since this literature-side computation used a from-scratch
Python transcription, not the repo's `laminate.py`.

I could not obtain Nesi's (1995) bound formula in transcribable form (paywalled,
no free text found), so I cannot evaluate it numerically at either point. However
Cherkaev2009-preprint states plainly, and I read this directly:

> [T, Cherkaev2009-preprint, Remark 4.4] "Nesi bound is generally not achievable
> by a structure."

and

> [T, Cherkaev2009-preprint p.4] "Later, the structures have been found in [9]
> that attain Nesi's bound in an asymptotic case when one material has infinite
> conductivity. Simultaneously, evidences were provided that the bound is not
> exact in the general case."

So Nesi's bound is a valid lower bound but is described by Cherkaev (2009) as
generally *not itself attained by any known structure* — it functions as a
comparison bound, not a construction target, at generic (f, sigma). This means
"best known lower bound in this regime" is ambiguous: Nesi's bound may
numerically be tighter than HS at these points even though not attainable, while
Cherkaev's improved, structurally-motivated bound coincides with HS at your two
points specifically. I was not able to resolve which is numerically tighter at
your two points because I could not get Nesi's explicit 3-phase formula from an
accessible source — this is a real gap, not a "not found" verdict on the
existence of a tighter bound.

### (2) Do published optimal structures include the rank-3 or rank-4 families in (i)/(ii)?

**VERDICT: not found as described, in ACN2007 — closest analogues are
different constructions.** (Cherkaev-Zhang 2011 not checked — see access note
above; Cherkaev2009-preprint's Section 8 "Optimal three-material structures" was
not read closely enough in the time available to rule out a match there either —
flag as **not fully checked**, not "absent.")

ACN2007 defines a hierarchy of named laminate families in its Section 4 ("New
optimal structures"), and none of the ones I found match your rank-3
construction (K1+K2 laminated with normal e1, K1+K3 laminated with the *same*
normal e1, then those two combined with normal e2) exactly:

- **T-structure** (ACN2007 p.15-16, Fig. 7a): "It is assembled as a sequence of
  laminates which depends upon two parameters. First, K1 and K3 are laminated
  with normal in the x1-direction. Then, the resulting structure is laminated
  with K2 with the normal in the x2-direction." [T] — this only ever laminates
  K1 with K3 once; K2 is added neat, not as a second K1-coated composite. Not a
  match: only one of your two rank-1 sub-laminates (K1+K3) appears, and K2 never
  gets coated with K1 at all in this structure.

- **Milton-Kohn matrix laminate** (ACN2007 p.11-12, Fig. 4c, Theorem 1): does
  build *two* separate K1-coated composites — K* = L(K1, L(K1,K2,n1,c1), n2,c2)
  and K*′ = L(K1, L(K1,K3,n3,c3), n4,c4) — closer in spirit to your (i), but (a)
  the two composites use *different, mutually orthogonal* normal pairs (n1⊥n2
  and n3⊥n4), not the shared single normal e1 your (i) specifies for both
  sub-laminates, and (b) the two composites are then mixed as a convex
  combination in the bound's linear parameter space ("if K* = K*′... the
  linearity of the bounds... allows us to mix the two constructions together in
  any way we wish"), not laminated together with a third normal e2. So this is
  structurally adjacent but not the same construction as your (i).

- **T²-structure** (ACN2007 p.20-21, Fig. 7c, eq. 26, Theorem 5): laminates K1
  and K3 with normal n1 to get K13, laminates K2 with K13 using normal n2 to get
  KT (this is the T-structure), separately laminates K1 and K3 again with
  normal n2 to get K13′, then laminates KT with K13′ using normal n1. Four
  design parameters, two separate K1-K3 sub-laminates with *swapped* normals
  relative to each other, not the K1-K2 / K1-K3-with-shared-normal pattern in
  your (i).

None of ACN2007's named families is literally "coat K2 with K1 along e1, coat K3
with K1 along e1, then laminate those two composites together along e2" — if
this is really what your rank-3 screen found, it reads as a genuine variant not
named in ACN2007's own catalogue (though it may still be implicit in their
general 6-parameter Figure 6 structure, which I did not fully unpack — see
below).

ACN2007's Section 4 opens by defining a fully general "orthogonal laminate of
high rank with six design parameters" (Fig. 6, p.14-15) of which T, coated-T, and
T² are named special cases; I did not have time to check whether setting its six
parameters a particular way reproduces your exact rank-3 construction as an
un-named special case. That would need a direct parameter-matching exercise
against eq. (19)/(26) in ACN2007, which I'm flagging as unfinished rather than
guessing at.

I did not reach Cherkaev-Zhang (2011) or check the "wheel assemblage" papers
(Cherkaev-Pruss) referenced in your problem statement's own background section —
out of scope for the time spent here, but worth a follow-up pass if the rank-3
structure turns out to be genuinely novel, since those are the most likely places
a different explicit construction could already exist.

---

## Summary for the ledger

- **Claim A**: Milton's 1981 condition is real and correctly attributed, but is
  published as a **volume-fraction threshold** m1 ≥ 2Θ(1−m2), Θ =
  k1(k3−k2)/[(k2+k1)(k3−k1)] (secondary-sourced via ACN2007 and
  Cherkaev2009-preprint, not the primary paper) — **not** in the "bound value
  vs sigma2" form your candidate states. Equivalence of the two forms is
  unverified algebra, not yet literature-confirmed. Record any use of the
  bound-vs-sigma2 phrasing as your own reformulation pending that check.
- **Claim B(1)**: at both your numeric points, Cherkaev's 2009 explicit bound
  (Theorem 7.1, transcribed and numerically verified above) **equals** the HS
  lower bound exactly — it does not improve on HS there. Nesi's bound could not
  be evaluated (no accessible formula) and is explicitly described by Cherkaev
  as "generally not achievable by a structure," so it's a real gap which bound
  is tightest at your points, not a closed question.
- **Claim B(2)**: your rank-3 and rank-4 constructions do not match any of the
  T-structure / Milton-Kohn matrix laminate / T²-structure families named in
  ACN2007 (2007). Cherkaev-Zhang (2011) and Cherkaev2009's own Section 8 were
  not checked closely enough to rule out a match there — genuine open item, not
  a clean "not found."
