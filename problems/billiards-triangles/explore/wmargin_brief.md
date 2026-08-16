You are a careful research mathematician. Work alone from this brief; you
have no other context. All angles are in DEGREES. Parameters: integers
a >= b >= 1, a >= 2 (the case a = b = 1 is excluded). Write out a complete,
elementary, self-contained proof — every inequality justified — or, if the
claim is FALSE, give an explicit numeric counterexample (a, b, s or j, t).
Do not hand-wave: a wrong "proof" is worse than an honest partial result.
If you can only prove the claim for a sub-range (e.g. smaller t, or a
restricted s), say exactly which sub-range and prove that cleanly.

Shared notation.
  A     = 90/a                (degrees)
  Bh    = 90/(b+1)
  bc    = 90(a-1)/(a(b+1))    (so A + (b+1)*bc = 90 exactly)
  theta_d = 90(a+b)/(a(b+1)) = A + bc  ("death" angle sum)
  DEATH SEGMENT:  alpha(t) = A - t,  beta(t) = bc + 2t,
                  theta(t) = alpha + beta = theta_d + t,
                  (b+1)*beta(t) = 90 - A + 2(b+1)t,
                  range  0 < t <= 45 k /(a(b+1)) with k in (0,1) a constant
                  you may choose (state the largest k your proof gives;
                  k = 1/2 is a good target; the claim FAILS at k = 1).
  BIRTH SEGMENT:  alpha(t) = A - t,  beta(t) = Bh - t,
                  theta(t) = A + Bh - 2t,
                  range  0 < t <= 22.5/(a(b+1)).
  On both segments 0 < alpha, 0 < beta, theta < 90.
  N2(t) [death] = sin(at) sin((b+1)beta(t)) sin(theta(t))
  N2(t) [birth] = sin(at) cos((b+1)t)      sin(theta(t))
  PB(t) [death] = sin(A - 2(b+1)t) cos(at) sin(theta(t))
                  [uses cos((b+1)beta) = sin(A - 2(b+1)t)]
  PB(t) [birth] = sin((b+1)t) cos(at) sin(theta(t))
Useful elementary facts you may use without proof: sin is concave on
[0,180]; sin(x) >= x/90 for x in [0,90] (x in degrees); sin(x) <= pi*x/180;
angle addition; product-to-sum identities.

YOUR TASK: {{value}}

Output format: (1) VERDICT: PROVED / PROVED-SUBRANGE / REFUTED / STUCK.
(2) The proof or counterexample, complete. (3) The exact constant k or
sub-range your argument covers. (4) Any step you are less than certain of,
flagged explicitly. Plain text/markdown, no code.
