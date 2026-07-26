/* graceful_core.c — exact graceful-labeling counter for trees (n <= 31).
 *
 * Protocol: reads trees from stdin, one per line:
 *     n p1 p2 ... p_{n-1}
 * where vertices are 0..n-1 in a connected order (each vertex i >= 1 has
 * parent p_i < i).  Writes one number per line to stdout: the number of
 * graceful labelings of the tree UP TO COMPLEMENTATION (f -> (n-1) - f).
 *
 * A graceful labeling is an injection f: V -> {0,...,n-1} with edge labels
 * |f(u)-f(v)| = {1,...,n-1}.  The complement of a graceful labeling is
 * graceful and always distinct from it, so the raw count is even; we count
 * raw/2 by restricting the root label:
 *     half = (#labelings with f(v0) <  ceil(m/2))
 *          + (#labelings with f(v0) == m/2, m even) / 2
 * with m = n-1.  (Complementation maps root label a -> m-a, pairing the
 * a < m/2 classes with a > m/2 and pairing labelings within a == m/2.)
 *
 * Compile: gcc -O3 -march=native -o graceful_core graceful_core.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, m;
static int par[32];
static int lab[32];
static unsigned long long cnt;

static void rec(int i, unsigned usedL, unsigned usedD)
{
    if (i == n) { cnt++; return; }
    const int pl = lab[par[i]];
    /* labels below parent's label: diff = pl - a */
    for (int a = 0; a < pl; a++) {
        if ((usedL >> a) & 1u) continue;
        const int d = pl - a;
        if ((usedD >> d) & 1u) continue;
        lab[i] = a;
        rec(i + 1, usedL | (1u << a), usedD | (1u << d));
    }
    /* labels above parent's label: diff = a - pl */
    for (int a = pl + 1; a <= m; a++) {
        if ((usedL >> a) & 1u) continue;
        const int d = a - pl;
        if ((usedD >> d) & 1u) continue;
        lab[i] = a;
        rec(i + 1, usedL | (1u << a), usedD | (1u << d));
    }
}

int main(void)
{
    char line[512];
    while (fgets(line, sizeof line, stdin)) {
        char *tok = strtok(line, " \t\n");
        if (!tok) continue;
        n = atoi(tok);
        if (n < 1 || n > 31) { fprintf(stderr, "bad n=%d\n", n); return 1; }
        for (int i = 1; i < n; i++) {
            tok = strtok(NULL, " \t\n");
            if (!tok) { fprintf(stderr, "short line\n"); return 1; }
            par[i] = atoi(tok);
            if (par[i] < 0 || par[i] >= i) { fprintf(stderr, "bad parent\n"); return 1; }
        }
        m = n - 1;
        if (n == 1) { printf("1\n"); continue; }

        unsigned long long half;
        /* root labels strictly below m/2 */
        cnt = 0;
        for (int a = 0; 2 * a < m; a++) {
            lab[0] = a;
            rec(1, 1u << a, 0u);
        }
        half = cnt;
        /* middle label when m even */
        if (m % 2 == 0) {
            cnt = 0;
            lab[0] = m / 2;
            rec(1, 1u << (m / 2), 0u);
            if (cnt % 2 != 0) { fprintf(stderr, "parity violation\n"); return 1; }
            half += cnt / 2;
        }
        printf("%llu\n", half);
    }
    return 0;
}
