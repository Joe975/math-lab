#!/usr/bin/env python3
# Attempt 022 skeptic verifier -- DRAFTED BY gemini-3.7-flash (effort medium)
# from the spec in swarm022_skeptic_brief.md alone (no repo code seen),
# then audited line-by-line by the director; run unmodified.
# Cross-family rule per docs/SWARM.md: the 020 certificates it re-checks
# were produced by director-written code on the 013 atanh-series log2 kit;
# this program uses digit-extraction log2 with directed dyadic rounding,
# its own bisection, and its own census -- nothing shared.
import sys
import json
from fractions import Fraction
import math

def certified_log2(r: Fraction, digits: int = 105, mantissa_bits: int = 200) -> tuple[Fraction, Fraction]:
    """
    Certified enclosure [lo, hi] of log2(r) for r > 0 using exact rational arithmetic
    with directed rounding mantissa truncation to avoid numerator blowup.
    """
    assert r > 0, "log2 only defined for positive rationals"
    
    # 1. Shift by integer power of 2 so r' = r / 2^k is in [1, 2)
    k = 0
    r_prime = r
    while r_prime >= 2:
        r_prime /= 2
        k += 1
    while r_prime < 1:
        r_prime *= 2
        k -= 1
        
    M = 1 << mantissa_bits
    
    # Lower track: floor rounding
    v_lo = r_prime
    lo_val = Fraction(0, 1)
    for d in range(1, digits + 1):
        v_lo = v_lo * v_lo
        # floor(v_lo * M) / M
        v_lo = Fraction((v_lo.numerator * M) // v_lo.denominator, M)
        if v_lo >= 2:
            lo_val += Fraction(1, 1 << d)
            v_lo = v_lo / 2

    # Upper track: ceil rounding
    v_hi = r_prime
    hi_val = Fraction(0, 1)
    for d in range(1, digits + 1):
        v_hi = v_hi * v_hi
        # ceil(v_hi * M) / M
        num = (v_hi.numerator * M + v_hi.denominator - 1) // v_hi.denominator
        v_hi = Fraction(num, M)
        if v_hi >= 2:
            hi_val += Fraction(1, 1 << d)
            v_hi = v_hi / 2

    rem = Fraction(1, 1 << digits)
    lo_bound = k + lo_val
    hi_bound = k + hi_val + rem
    return lo_bound, hi_bound

def plackett_bracket_and_root(x: Fraction, y: Fraction, t: Fraction, target_width: Fraction = Fraction(1, 1 << 90)) -> tuple[Fraction, Fraction]:
    """
    Enclose the unique root z_t(x, y) of F(z) = z*(1-x-y+z) - t*(x-z)*(y-z)
    in (x*y, min(x, y)) using rational bisection to width <= target_width.
    """
    def F(z: Fraction) -> Fraction:
        return z * (1 - x - y + z) - t * (x - z) * (y - z)
    
    lo = x * y
    hi = min(x, y)
    
    f_lo = F(lo)
    f_hi = F(hi)
    
    if not (f_lo < 0 and f_hi > 0):
        raise ValueError(f"Plackett bracket invalid: F(lo)={f_lo}, F(hi)={f_hi}")
        
    while (hi - lo) > target_width:
        mid = (lo + hi) / 2
        f_mid = F(mid)
        if f_mid <= 0:
            lo = mid
        else:
            hi = mid
            
    return lo, hi

def compute_weight_enclosure(x: Fraction, y: Fraction, t: Fraction, z_lo: Fraction, z_hi: Fraction) -> tuple[Fraction, Fraction]:
    """
    Compute certified enclosure [w_lo, w_hi] for w(x,y) = | hp(z) * dz | on z in [z_lo, z_hi].
    """
    # hp(z) = log2((1-z)/z), decreasing in z
    arg_lo = (1 - z_hi) / z_hi
    arg_hi = (1 - z_lo) / z_lo
    hp_lo, _ = certified_log2(arg_lo)
    _, hp_hi = certified_log2(arg_hi)
    
    # dz = (x-z)*(y-z) / ((1-x-y+2z) + t*(x+y-2z))
    # Numerator is positive and decreasing in z on z < min(x,y)
    num_lo = (x - z_hi) * (y - z_hi)
    num_hi = (x - z_lo) * (y - z_lo)
    
    # Denominator D(z) is linear with negative slope (2 - 2t < 0), decreasing in z
    def D(z: Fraction) -> Fraction:
        return (1 - x - y + 2 * z) + t * (x + y - 2 * z)
        
    den_lo = D(z_hi)
    den_hi = D(z_lo)
    if den_lo <= 0 or den_hi <= 0:
        raise ValueError("Denominator of dz is not strictly positive on the bracket")
        
    dz_lo = num_lo / den_hi
    dz_hi = num_hi / den_lo
    
    # Box product of hp and dz
    prods = [hp_lo * dz_lo, hp_lo * dz_hi, hp_hi * dz_lo, hp_hi * dz_hi]
    p_lo = min(prods)
    p_hi = max(prods)
    
    # Absolute value
    if p_lo >= 0:
        w_lo = p_lo
        w_hi = p_hi
    elif p_hi <= 0:
        w_lo = -p_hi
        w_hi = -p_lo
    else:
        w_lo = Fraction(0, 1)
        w_hi = max(-p_lo, p_hi)
        
    return w_lo, w_hi

def run_selftests():
    # 1. Product anchor
    n = 3
    t = Fraction(3, 2)
    u = {}
    for A in range(1 << n):
        pop = bin(A).count("1")
        u[A] = (Fraction(3, 10) ** pop) * (Fraction(7, 10) ** (3 - pop))
    
    # Build coupling
    pi = {}
    Z = Fraction(0, 1)
    for A in range(1 << n):
        for B in range(1 << n):
            pop = bin(A & B).count("1")
            val = u[A] * u[B] * (t ** pop)
            pi[(A, B)] = val
            Z += val
            
    # Check tables
    for i in range(1, n):
        m_i = (1 << (i - 1)) - 1
        c_i = 1 << (i - 1)
        for a in range(1 << (i - 1)):
            for b in range(1 << (i - 1)):
                p00 = sum(pi[(A, B)] for A in range(1 << n) for B in range(1 << n)
                          if (A & m_i) == a and (B & m_i) == b and (A & c_i) == 0 and (B & c_i) == 0) / Z
                p01 = sum(pi[(A, B)] for A in range(1 << n) for B in range(1 << n)
                          if (A & m_i) == a and (B & m_i) == b and (A & c_i) == 0 and (B & c_i) != 0) / Z
                p10 = sum(pi[(A, B)] for A in range(1 << n) for B in range(1 << n)
                          if (A & m_i) == a and (B & m_i) == b and (A & c_i) != 0 and (B & c_i) == 0) / Z
                p11 = sum(pi[(A, B)] for A in range(1 << n) for B in range(1 << n)
                          if (A & m_i) == a and (B & m_i) == b and (A & c_i) != 0 and (B & c_i) != 0) / Z
                if p00 > 0 and p01 > 0 and p10 > 0 and p11 > 0:
                    OR = (p00 * p11) / (p01 * p10)
                    assert OR == t, f"Product anchor failed: OR = {OR} != t = {t}"
                    d_lo, d_hi = certified_log2(OR / t)
                    assert d_lo <= 0 <= d_hi, "Product anchor dev does not contain 0"
    print("SELFTEST 1 (Product anchor): PASSED")

    # 2. log2 anchor
    l2_lo, l2_hi = certified_log2(Fraction(2, 1))
    assert l2_lo <= 1 <= l2_hi and (l2_hi - l2_lo) <= Fraction(1, 1 << 90), f"log2(2) anchor failed: [{l2_lo}, {l2_hi}]"
    
    l8_lo, l8_hi = certified_log2(Fraction(8, 1))
    assert l8_lo <= 3 <= l8_hi and (l8_hi - l8_lo) <= Fraction(1, 1 << 90), f"log2(8) anchor failed: [{l8_lo}, {l8_hi}]"
    
    lh_lo, lh_hi = certified_log2(Fraction(1, 2))
    assert lh_lo <= -1 <= lh_hi and (lh_hi - lh_lo) <= Fraction(1, 1 << 90), f"log2(1/2) anchor failed: [{lh_lo}, {lh_hi}]"
    print("SELFTEST 2 (log2 anchor): PASSED")

    # 3. Plackett anchor
    z_lo, z_hi = plackett_bracket_and_root(Fraction(1, 2), Fraction(1, 2), Fraction(9, 1))
    assert z_lo <= Fraction(3, 8) <= z_hi and (z_hi - z_lo) <= Fraction(1, 1 << 90), "Plackett anchor root 3/8 failed"
    print("SELFTEST 3 (Plackett anchor): PASSED")

def process_witness(w: dict) -> dict:
    name = w["name"]
    kind = w["kind"]
    n = w["n"]
    t = Fraction(w["t"][0], w["t"][1])
    
    u_raw = w["u"]
    u = {}
    for k_str, val in u_raw.items():
        A = int(k_str, 2)
        u[A] = Fraction(val[0], val[1])
        
    for A in range(1 << n):
        if A not in u:
            u[A] = Fraction(0, 1)

    # Compute unnormalized coupling
    pi = {}
    Z = Fraction(0, 1)
    for A in range(1 << n):
        uA = u[A]
        if uA == 0:
            continue
        for B in range(1 << n):
            uB = u[B]
            if uB == 0:
                continue
            pop = bin(A & B).count("1")
            val = uA * uB * (t ** pop)
            pi[(A, B)] = val
            Z += val
            
    assert Z > 0, "Z must be positive"

    # Marginals
    marginals = []
    for j in range(n):
        bit_j = 1 << j
        m_j = sum(pi[(A, B)] for (A, B) in pi if (A & bit_j) != 0) / Z
        marginals.append(m_j)
        
    max_marginal = max(marginals) if marginals else Fraction(0, 1)
    in_regime = all(m < Fraction(38271, 100000) for m in marginals)

    # Census tables
    agg_lo = Fraction(0, 1)
    agg_hi = Fraction(0, 1)
    num_lo = Fraction(0, 1)
    num_hi = Fraction(0, 1)
    
    n_tables = 0
    mixed_zero_tables = 0
    all_rplus = True

    for i in range(1, n):
        m_i = (1 << (i - 1)) - 1
        c_i = 1 << (i - 1)
        for a in range(1 << (i - 1)):
            for b in range(1 << (i - 1)):
                p00 = Fraction(0, 1)
                p01 = Fraction(0, 1)
                p10 = Fraction(0, 1)
                p11 = Fraction(0, 1)
                
                for (A, B), val in pi.items():
                    if (A & m_i) == a and (B & m_i) == b:
                        A_ci = (A & c_i) != 0
                        B_ci = (B & c_i) != 0
                        if not A_ci and not B_ci:
                            p00 += val
                        elif not A_ci and B_ci:
                            p01 += val
                        elif A_ci and not B_ci:
                            p10 += val
                        else:
                            p11 += val
                            
                p00 /= Z
                p01 /= Z
                p10 /= Z
                p11 /= Z
                
                mass = p00 + p01 + p10 + p11
                if mass == 0:
                    continue
                    
                is_all_pos = (p00 > 0 and p01 > 0 and p10 > 0 and p11 > 0)
                is_full_zero_row_col = (
                    (p00 == 0 and p01 == 0) or
                    (p10 == 0 and p11 == 0) or
                    (p00 == 0 and p10 == 0) or
                    (p01 == 0 and p11 == 0)
                )
                
                if not (is_all_pos or is_full_zero_row_col):
                    mixed_zero_tables += 1

                if is_all_pos:
                    n_tables += 1
                    x = (p00 + p01) / mass
                    y = (p00 + p10) / mass
                    OR = (p00 * p11) / (p01 * p10)
                    
                    dev_lo, dev_hi = certified_log2(OR / t)
                    
                    # Aggregate sum
                    agg_lo += mass * dev_lo
                    agg_hi += mass * dev_hi
                    
                    # Weight & root enclosure
                    z_lo, z_hi = plackett_bracket_and_root(x, y, t)
                    if z_hi >= Fraction(1, 2):
                        all_rplus = False
                        
                    w_lo, w_hi = compute_weight_enclosure(x, y, t, z_lo, z_hi)
                    
                    # Interval product mass * w(x,y) * dev
                    prods = [
                        w_lo * dev_lo,
                        w_lo * dev_hi,
                        w_hi * dev_lo,
                        w_hi * dev_hi
                    ]
                    term_lo = mass * min(prods)
                    term_hi = mass * max(prods)
                    
                    num_lo += term_lo
                    num_hi += term_hi

    out = {
        "name": name,
        "num_lo": str(num_lo) if kind != "aggregate" else None,
        "num_hi": str(num_hi) if kind != "aggregate" else None,
        "agg_lo": str(agg_lo),
        "agg_hi": str(agg_hi),
        "in_regime": in_regime,
        "max_marginal": float(max_marginal),
        "all_rplus": all_rplus if kind == "mm_abs_rplus" else None,
        "n_tables": n_tables,
        "mixed_zero_tables": mixed_zero_tables
    }
    return out

def main():
    run_selftests()
    
    if len(sys.argv) < 2:
        return
        
    input_path = sys.argv[1]
    with open(input_path, "r", encoding="utf-8") as f:
        data = json.load(f)
        
    for w in data.get("witnesses", []):
        res = process_witness(w)
        print(json.dumps(res))

if __name__ == "__main__":
    main()
