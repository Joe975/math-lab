/* graceful_core.c — exact graceful-labeling counter for trees (n <= 31).
 *
 * Protocol: reads trees from stdin, one per line:
 *     n p1 p2 ... p_{n-1}
 * where vertices are 0..n-1 in a BFS order (each vertex i >= 1 has parent
 * p_i < i, and the children of any vertex are contiguous in the order).
 * Writes one number per line to stdout: the number of graceful labelings of
 * the tree UP TO COMPLEMENTATION (f -> (n-1) - f).
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
 * Search: vertices are labeled in the given BFS order, so each new vertex's
 * edge difference to its parent is determined immediately.  Pruning:
 *   - candidate labels iterated via bit tricks over the unused-label mask;
 *   - realizability of the LARGEST unused difference D: some future edge
 *     must be able to realize D.  Every future edge joins an unplaced vertex
 *     to either an "open" placed vertex (a placed vertex that still has
 *     unplaced children) or another unplaced vertex, so D is realizable only
 *     if (two unused labels D apart) or (an open label x with x+D or x-D
 *     unused) exists — three mask operations.
 *
 * Compile: gcc -O3 -march=native -o graceful_core graceful_core.c
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int n, m;
static int par[32];
static int lab[32];
static int haschild[32];
static int islast[32];      /* i is the last child of par[i] */
static unsigned LABMASK;    /* bits 0..m   */
static unsigned DIFFMASK;   /* bits 1..m   */
static unsigned long long cnt;

static void rec(int i, unsigned usedL, unsigned usedD, unsigned open)
{
    if (i == n) { cnt++; return; }

    const unsigned U = ~usedL & LABMASK;

    /* prune: largest unused difference must still be realizable */
    {
        const unsigned ud = ~usedD & DIFFMASK;   /* nonzero: i <= n-1 */
        const int D = 31 - __builtin_clz(ud);
        if (!((U & (U >> D)) | (open & (U >> D)) | (open & (U << D))))
            return;
    }

    const int pl = lab[par[i]];
    const unsigned childbit = haschild[i];       /* 0 or 1 */
    unsigned newopen_base = open;
    if (islast[i]) newopen_base &= ~(1u << pl);  /* parent closes */

    for (unsigned rest = U; rest; rest &= rest - 1) {
        const int a = __builtin_ctz(rest);
        const int d = a > pl ? a - pl : pl - a;
        if ((usedD >> d) & 1u) continue;
        lab[i] = a;
        rec(i + 1, usedL | (1u << a), usedD | (1u << d),
            newopen_base | (childbit << a));
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
        LABMASK = (m + 1 == 32) ? 0xffffffffu : ((1u << (m + 1)) - 1u);
        DIFFMASK = LABMASK & ~1u;

        memset(haschild, 0, sizeof haschild);
        for (int i = 1; i < n; i++) haschild[par[i]] = 1;
        int lastchild[32];
        memset(lastchild, -1, sizeof lastchild);
        for (int i = 1; i < n; i++) lastchild[par[i]] = i;
        for (int i = 1; i < n; i++) islast[i] = (lastchild[par[i]] == i);

        unsigned long long half;
        /* root labels strictly below m/2 */
        cnt = 0;
        for (int a = 0; 2 * a < m; a++) {
            lab[0] = a;
            rec(1, 1u << a, 0u, 1u << a);   /* root is open: n >= 2 */
        }
        half = cnt;
        /* middle label when m even */
        if (m % 2 == 0) {
            cnt = 0;
            lab[0] = m / 2;
            rec(1, 1u << (m / 2), 0u, 1u << (m / 2));
            if (cnt % 2 != 0) { fprintf(stderr, "parity violation\n"); return 1; }
            half += cnt / 2;
        }
        printf("%llu\n", half);
    }
    return 0;
}
