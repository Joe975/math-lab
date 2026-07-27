/* Fast kernel for Zaremba's conjecture (harness, tier 0).
 *
 * Canonical continued fractions via the Euclidean algorithm; the last
 * partial quotient is automatically >= 2 for reduced m/n with 0 < m < n.
 * Contract: problems/zaremba/PROBLEM.md. Reference semantics:
 * harness/zaremba/zaremba.py (this kernel must agree with it exactly).
 *
 * Build:  cc -O2 -o zscan zscan.c
 *
 * Modes:
 *   zscan check <A> <lo> <hi> [witness_file]
 *     For each n in [lo,hi], smallest m with gcd(m,n)=1 and K(m/n) <= A.
 *     Prints one summary line per million and any n with NO witness
 *     (a Zaremba counterexample for this A) as "FAIL n".
 *     With witness_file, writes "n m" per line for independent checking.
 *
 *   zscan exact <lo> <hi> [table_file]
 *     Exact z(n) for each n in [lo,hi] by exhaustive scan over m.
 *     Prints a histogram of z values and every n with z(n) >= 4.
 *     With table_file, writes "n z" per line.
 *
 * Everything is unsigned 64-bit integer arithmetic; no floating point.
 */

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef uint64_t u64;

/* K(m/n) if gcd(m,n)=1 and all quotients <= cap, else 0.
 * Early-aborts as soon as a quotient reaches `cap`+1 or the pair is seen
 * to be non-coprime, so cap doubles as a pruning bound. */
static inline unsigned kmax_capped(u64 n, u64 m, unsigned cap) {
    u64 p = n, q = m;
    unsigned worst = 0;
    while (q) {
        u64 a = p / q;
        u64 r = p % q;
        if (a > cap) return 0;
        if (a > worst) worst = (unsigned)a;
        p = q;
        q = r;
    }
    return p == 1 ? worst : 0;
}

static int mode_check(unsigned A, u64 lo, u64 hi, const char *wpath) {
    FILE *wf = NULL;
    if (wpath) {
        wf = fopen(wpath, "w");
        if (!wf) { perror(wpath); return 2; }
    }
    u64 fails = 0, slice_tries = 0, slice_n = 0;
    u64 worst_tries = 0, worst_tries_n = 0;
    for (u64 n = lo; n <= hi; n++) {
        u64 m, tries = 0, found = 0;
        /* Any m/n with all quotients <= A is at least the infimum of the
         * bounded-quotient set, [0; A,1,A,1,...]; for A = 5 that is
         * (3*sqrt(5)-5)/10 = 0.1708203932... Starting below it only burns
         * divisions on m that cannot qualify, so start just under the
         * infimum (17082/100000 < it); the first qualifying m is unchanged.
         * For A != 5 fall back to the coarse bound n/(A+1). */
        u64 start = A == 5 ? (n / 100000) * 17082 + (n % 100000) * 17082 / 100000
                           : n / (A + 1);
        for (m = start + 1; m < n; m++) {
            tries++;
            if (kmax_capped(n, m, A)) { found = 1; break; }
        }
        if (n == 1) { found = 1; m = 0; tries = 0; } /* z(1)=1 by convention */
        if (!found) {
            printf("FAIL %" PRIu64 "\n", n);
            fflush(stdout);
            fails++;
            continue;
        }
        if (wf) fprintf(wf, "%" PRIu64 " %" PRIu64 "\n", n, m);
        slice_tries += tries;
        slice_n++;
        if (tries > worst_tries) { worst_tries = tries; worst_tries_n = n; }
        if (n % 1000000 == 0) {
            printf("slice_to %" PRIu64 "  mean_tries %.3f  "
                   "worst_tries %" PRIu64 " at n=%" PRIu64 "\n",
                   n, (double)slice_tries / (double)slice_n,
                   worst_tries, worst_tries_n);
            fflush(stdout);
            slice_tries = slice_n = 0;
            worst_tries = worst_tries_n = 0;
        }
    }
    if (wf) fclose(wf);
    printf("done [%" PRIu64 ",%" PRIu64 "] A=%u failures=%" PRIu64 "\n",
           lo, hi, A, fails);
    return fails ? 1 : 0;
}

static int mode_exact(u64 lo, u64 hi, const char *tpath) {
    FILE *tf = NULL;
    if (tpath) {
        tf = fopen(tpath, "w");
        if (!tf) { perror(tpath); return 2; }
    }
    u64 hist[65];
    memset(hist, 0, sizeof hist);
    for (u64 n = lo; n <= hi; n++) {
        unsigned best;
        if (n < 2) {
            best = 1;
        } else {
            best = 64; /* scan below always improves on this for n >= 2 */
            for (u64 m = 1; m < n; m++) {
                /* quotients >= best can't improve, so prune at best-1 */
                unsigned k = kmax_capped(n, m, best - 1);
                if (k && k < best) {
                    best = k;
                    if (best == 2) break; /* canonical minimum for n >= 2 */
                }
            }
            if (best == 64) {
                /* unreachable: m=1 gives K=n; means n exceeded the cap */
                fprintf(stderr, "z(%" PRIu64 ") > 63; widen hist\n", n);
                return 2;
            }
        }
        hist[best]++;
        if (best >= 4)
            printf("z %" PRIu64 " = %u\n", n, best);
        if (tf) fprintf(tf, "%" PRIu64 " %u\n", n, best);
    }
    if (tf) fclose(tf);
    printf("histogram [%" PRIu64 ",%" PRIu64 "]:\n", lo, hi);
    for (int v = 1; v < 65; v++)
        if (hist[v])
            printf("  z=%d: %" PRIu64 "\n", v, hist[v]);
    return 0;
}

int main(int argc, char **argv) {
    if (argc >= 4 && strcmp(argv[1], "check") == 0) {
        unsigned A = (unsigned)strtoul(argv[2], NULL, 10);
        u64 lo = strtoull(argv[3], NULL, 10);
        u64 hi = strtoull(argv[4], NULL, 10);
        return mode_check(A, lo, hi, argc > 5 ? argv[5] : NULL);
    }
    if (argc >= 4 && strcmp(argv[1], "exact") == 0) {
        u64 lo = strtoull(argv[2], NULL, 10);
        u64 hi = strtoull(argv[3], NULL, 10);
        return mode_exact(lo, hi, argc > 4 ? argv[4] : NULL);
    }
    fprintf(stderr,
            "usage: zscan check <A> <lo> <hi> [witness_file]\n"
            "       zscan exact <lo> <hi> [table_file]\n");
    return 2;
}
