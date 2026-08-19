# Swarm brief — attempt 020 attack pack (union-closed live candidates), template

You are a careful research mathematician. Work alone from this brief; you
have no other context and no web access. Give complete derivations for
anything you claim; flag every uncertain step SPECULATION. A wrong
confident claim is worse than an honest partial answer. Plain text/math,
no code.

## Setup (all of this you may take as given)

Context: entropy attacks on the union-closed sets conjecture. mu is a
probability measure on {0,1}^n (bit strings = sets; bit value 1 means
"element present"; every element-marginal P(bit_i = 1) <= 0.38271). Two
copies (A, B) of mu are coupled coordinate-by-coordinate: at coordinate
i, conditionally on the histories (prefixes) a = A_{<i}, b = B_{<i}, the
pair (A_i, B_i) has a 2x2 distribution ("table") with

  zero-margins  x = P(A_i = 0 | a),  y = P(B_i = 0 | b)   (both in (0,1)
                at a NONDEGENERATE history; degenerate histories are
                excluded from all averages below),
  both-zero probability z~ (Frechet range: max(0, x+y-1) < z~ < min(x,y)),
  realized odds ratio  OR = z~ (1-x-y+z~) / ((x-z~)(y-z~))  in (0, inf).

The coupling carries a target parameter lambda > 0. The Plackett table at
odds ratio rho > 0 with margins (x, y) has both-zero probability
z_rho(x,y) = the unique root of z(1-x-y+z) = rho (x-z)(y-z) in the
Frechet range. Define h2(z) = binary entropy in bits of {z, x-z, y-z,
1-x-y+z}'s first cell only, i.e. h2(z) = -z log2 z - (1-z) log2(1-z).
The per-history "gain sensitivity" is

  sigma_lambda(x, y) = d/dlambda [ h2( z_{2^lambda}(x, y) ) ]
                     = h2'(z) * dz/drho * rho * ln 2   at rho = 2^lambda.

Facts you may use freely (all established in this lab or standard):
(F1) sigma_lambda > 0 iff z_{2^lambda}(x,y) < 1/2; the sign flips at
     z = 1/2. For equal margins x = y = m, z_1(m,m) = m^2.
(F2) For product measures mu = Bern(p)^n, every history has OR = 2^lambda
     exactly (the coupling realizes its target), so any average of
     f(OR)-type deviations vanishes identically at products.
(F3) Mass-weighted averages below are over all nondegenerate histories
     (a,b) of all coordinates i < n, weighted by the coupling's history
     mass m_ab.
(F4) [Perfect-square first order] For fixed mu, the per-coordinate
     average M_i(lambda) - lambda = lambda * || sum_{a in N_i} nu(a) D_a ||^2
     + O(lambda^2), where nu is normalized prefix mass and
     D_a = E_mu[future coords | prefix a, bit 1] - E_mu[future | a, bit 0]
     (a vector in R^{n-i}). So the first-order-in-lambda coefficient a1 of
     the aggregated deviation A(mu, lambda) = E[log2 OR - lambda] is a
     mass-weighted sum over i of squared NORMS of SIGNED VECTOR SUMS —
     it can vanish by cancellation across prefixes a, not only at
     product measures.
(F5) [Known kill] Replicating a fixed 10-atom two-block witness r times
     (shared "block marker" coordinate; re-optimized shared weights)
     drives A(mu, lambda) < 0 for n >= ~90 at lambda in [2, 2.5], with
     all marginals <= 0.309. The deficit channel is a "light slice":
     histories with NEAR-DEGENERATE margins (x or y close to 0 or 1)
     carrying a Theta(n) deficit of log2 OR below lambda.
(F6) [Sensitivity-weighted kills] The signed weightings
     E[sigma * (log2 OR - lambda)] (secant and at-target readings) are
     both NEGATIVE on the witness itself: the surplus histories sit at
     z~ > 1/2 where sigma < 0, so the weighting flips the biggest
     positive terms. Only the UNSIGNED weighting
       MM_abs = E[ |sigma_lambda| (log2 OR - lambda) ] / E[ |sigma_lambda| ]
     survives everything tried: it is >= +0.97 on the replicated ladders,
     ground to +0.246 by adversarial re-weighting at n = 48, and to
     +0.0101 by a free-support search at n = 8. It resists the ladder
     because the light-slice deficit has near-degenerate margins, where
     |sigma| -> 0.
(F7) [Window law] The chain-rule assembly can only use
     lambda <= lambda_win(n) ~ 4.847/(n-3). Restricted to that window,
     the aggregate A on all known families behaves as
     A ~ a1(n) lambda + a2(n) lambda^2 with a1 >= 0 (by F4), a2 < 0,
     and on the killing ladder genre a1(n)*n GROWS with n (4.7 at n=48
     to 18.1 at n=320), so the window variant survives there. For the
     window variant to die at lambda = c/n one needs a family with
     a1(n) = O(1/n) (e.g. via F4's cancellation) while a2 stays
     negative and bounded away from 0.

## YOUR TASK

{{value}}

Output format:
VERDICT: one line (e.g. CONSTRUCTION PROPOSED / PROVED / PARTIAL / STUCK / NO-GO ARGUED).
BODY: the construction with explicit numbers, or the proof, complete.
FLAGS: every step you are less than certain of.
