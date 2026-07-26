/* es_fcount.c -- exact Erdos-Straus representation counts f(p), parallel.
 *
 * f(p) = #{(x,y,z) : x <= y <= z, 4/p = 1/x + 1/y + 1/z}, p prime.
 *
 * Counting identity used, cross-checked against an independent pure-Python
 * counter:
 *   x in (p/4, 3p/4]; A = 4x - p, B = p x; solutions of 1/y + 1/z = A/B
 *   <-> divisors d | B^2, d <= B, dmin <= d, A | B + d, where
 *   dmin = max(0, 2x(2x-p)) enforces y >= x and d <= B enforces y <= z.
 *   Since p is prime and p > x, divisors of B^2 = p^2 x^2 with d <= B are
 *   d0 (d0 | x^2) and d0*p (d0 | x^2, d0 <= x).
 *
 * Differences from 001's helper: OpenMP parallel over primes, prime-range
 * chunking with a residue-class filter, and a prime-list mode for sampling.
 *
 * Usage:
 *   es_fcount pmin pmax mod res    f(p) for primes p in [pmin, pmax),
 *                                  p % mod == res  (mod 0 = all primes)
 *   es_fcount -l FILE pmax         f(p) for primes listed in FILE (one per
 *                                  line, each < pmax; primality not checked)
 * Output: lines "p f(p)", ascending in p (list mode: input order).
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

static int *spf;                 /* smallest prime factor, up to xmax */
static uint64_t **divs;          /* divisors of x^2, per x in [xlo, xmax] */
static int *ndivs;
static int xlo_g, xmax_g;

static void build_tables(int xlo, int xmax) {
    xlo_g = xlo; xmax_g = xmax;
    spf = calloc((size_t)xmax + 1, sizeof(int));
    if (!spf) { fprintf(stderr, "oom spf\n"); exit(1); }
    for (int i = 2; i <= xmax; i++)
        if (!spf[i]) for (long j = i; j <= xmax; j += i)
            if (!spf[j]) spf[j] = (int)i;
    divs  = malloc(((size_t)xmax + 1) * sizeof(uint64_t *));
    ndivs = calloc((size_t)xmax + 1, sizeof(int));
    if (!divs || !ndivs) { fprintf(stderr, "oom divs\n"); exit(1); }
    for (int x = xlo; x <= xmax; x++) {
        int pr[20], ex[20], np = 0, t = x;
        while (t > 1) {
            int q = spf[t], e = 0;
            while (t % q == 0) { t /= q; e++; }
            pr[np] = q; ex[np] = 2 * e; np++;
        }
        int cnt = 1;
        for (int i = 0; i < np; i++) cnt *= ex[i] + 1;
        uint64_t *D = malloc((size_t)cnt * sizeof(uint64_t));
        if (!D) { fprintf(stderr, "oom D\n"); exit(1); }
        D[0] = 1; int m = 1;
        for (int i = 0; i < np; i++) {
            int m0 = m; uint64_t pe = 1;
            for (int e = 1; e <= ex[i]; e++) {
                pe *= (uint64_t)pr[i];
                for (int j = 0; j < m0; j++) D[m++] = D[j] * pe;
            }
        }
        divs[x] = D; ndivs[x] = cnt;
    }
}

static long fcount(long p) {
    long count = 0;
    int xlow = (int)(p / 4 + 1), xhi = (int)((3 * p) / 4);
    for (int x = xlow; x <= xhi; x++) {
        uint64_t A = (uint64_t)(4L * x - p);
        uint64_t B = (uint64_t)p * (uint64_t)x;
        long t2 = 2L * x - p;
        uint64_t dmin = t2 > 0 ? (uint64_t)(2L * x * t2) : 0;
        uint64_t *D = divs[x];
        int nd = ndivs[x];
        for (int i = 0; i < nd; i++) {
            uint64_t d0 = D[i];
            /* d = d0: d0 <= x^2 < px = B always (x < p) */
            if (d0 >= dmin && (B + d0) % A == 0) count++;
            if (d0 <= (uint64_t)x) {          /* d = d0*p <= B */
                uint64_t d = d0 * (uint64_t)p;
                if (d >= dmin && (B + d) % A == 0) count++;
            }
        }
    }
    return count;
}

int main(int argc, char **argv) {
    long *plist = NULL; long np = 0;
    long pmin = 0, pmax = 0;
    if (argc >= 3 && !strcmp(argv[1], "-l")) {
        FILE *f = fopen(argv[2], "r");
        if (!f) { fprintf(stderr, "cannot open %s\n", argv[2]); return 1; }
        long cap = 1024; plist = malloc(cap * sizeof(long));
        long v;
        while (fscanf(f, "%ld", &v) == 1) {
            if (np == cap) { cap *= 2; plist = realloc(plist, cap * sizeof(long)); }
            plist[np++] = v;
            if (v > pmax) pmax = v;
            if (pmin == 0 || v < pmin) pmin = v;
        }
        fclose(f);
        pmax += 1;
    } else if (argc >= 5) {
        pmin = atol(argv[1]); pmax = atol(argv[2]);
        long mod = atol(argv[3]), res = atol(argv[4]);
        if (pmin < 2) pmin = 2;
        /* sieve primes in [pmin, pmax) */
        char *comp = calloc((size_t)pmax, 1);
        for (long i = 2; i * i < pmax; i++)
            if (!comp[i]) for (long j = i * i; j < pmax; j += i) comp[j] = 1;
        long cap = 1024; plist = malloc(cap * sizeof(long));
        for (long p = pmin; p < pmax; p++) {
            if (comp[p]) continue;
            if (mod > 0 && p % mod != res) continue;
            if (np == cap) { cap *= 2; plist = realloc(plist, cap * sizeof(long)); }
            plist[np++] = p;
        }
        free(comp);
    } else {
        fprintf(stderr, "usage: %s pmin pmax mod res | %s -l FILE\n",
                argv[0], argv[0]);
        return 1;
    }
    if (np == 0) return 0;
    int xlo = (int)(pmin / 4 + 1), xmax = (int)((3 * pmax) / 4 + 1);
    build_tables(xlo, xmax);
    long *out = malloc(np * sizeof(long));
    #pragma omp parallel for schedule(dynamic, 1)
    for (long i = 0; i < np; i++)
        out[i] = fcount(plist[i]);
    for (long i = 0; i < np; i++)
        printf("%ld %ld\n", plist[i], out[i]);
    return 0;
}
