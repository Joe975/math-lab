/* cycle_filter.c - fast fixed-length-cycle presence filter for cubic graphs.
 *
 * Reads graph6 lines (n <= 62) on stdin.  For each graph computes:
 *   girth (exact), and presence of cycles of length exactly 4, 8, 16.
 * Prints notable graphs and a JSON summary line at the end.
 *
 * Output lines:
 *   SURV8 <g6> girth=<g> has16=<0|1>   -- C8-free graph (survivor)
 *   P8ONLY <g6> girth=<g>              -- pow2 intersection exactly {8}
 *                                         (no C4, has C8, no C16)
 *   CAND <g6>                          -- avoids C4, C8 and C16 entirely
 *                                         (counterexample candidate!)
 *   SUMMARY {...json...}               -- totals at EOF
 *
 * Flags:
 *   -v : per-graph verbose line "V <g6> girth=<g> c4=<b> c8=<b> c16=<b>"
 *        (for validation against an independent implementation)
 *
 * Algorithm: existence of a cycle of length exactly k is tested by DFS over
 * simple paths.  A cycle is rooted at its minimum vertex, so the DFS from
 * root r only visits vertices > r; pruning uses BFS distances back to the
 * root (computed in the subgraph induced by {v : v >= r}): a partial path at
 * vertex v with (k - depth) edges left to place cannot close unless
 * dist(v, root) <= k - depth.  Early exit on first cycle found.
 * Girth is computed as the smallest k in 3..n with a length-k cycle.
 *
 * Build: gcc -O2 -o cycle_filter cycle_filter.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define MAXN 62
#define MAXDEG 8

static int N;
static uint64_t adjmask[MAXN];
static int nbr[MAXN][MAXDEG], deg[MAXN];
static int distr[MAXN];          /* BFS distance from current root */

/* parse graph6 into globals; return 0 on failure */
static int parse_graph6(const char *s)
{
    int i, j, k, len = (int)strlen(s);
    if (len == 0) return 0;
    N = s[0] - 63;
    if (N < 1 || N > MAXN) return 0;
    memset(adjmask, 0, sizeof adjmask);
    memset(deg, 0, sizeof deg);
    k = 0;                       /* bit index into the edge bitstream */
    {
        int need = N * (N - 1) / 2;
        int have = (len - 1) * 6;
        if (have < need) return 0;
        for (j = 1; j < N; j++) {
            for (i = 0; i < j; i++) {
                int byte = s[1 + k / 6] - 63;
                if (byte < 0 || byte > 63) return 0;
                if ((byte >> (5 - k % 6)) & 1) {
                    adjmask[i] |= 1ULL << j;
                    adjmask[j] |= 1ULL << i;
                    if (deg[i] >= MAXDEG || deg[j] >= MAXDEG) return 0;
                    nbr[i][deg[i]++] = j;
                    nbr[j][deg[j]++] = i;
                }
                k++;
            }
        }
    }
    return 1;
}

/* DFS for a cycle of length exactly k through `root`, all vertices > root.
 * depth = number of edges from root to u along the current path. */
static int dfs_k(int u, uint64_t visited, int depth, int k, int root)
{
    int i;
    if (depth == k - 1)
        return (int)((adjmask[u] >> root) & 1ULL);
    for (i = 0; i < deg[u]; i++) {
        int v = nbr[u][i];
        if (v <= root) continue;             /* root handled at depth k-1 */
        if ((visited >> v) & 1ULL) continue;
        if (distr[v] > k - depth - 1) continue;   /* cannot get back */
        if (dfs_k(v, visited | (1ULL << v), depth + 1, k, root))
            return 1;
    }
    return 0;
}

/* BFS distances from root within the subgraph induced by {v >= root} */
static void bfs_from(int root)
{
    int queue[MAXN], qh = 0, qt = 0, i;
    for (i = 0; i < N; i++) distr[i] = 127;
    distr[root] = 0;
    queue[qt++] = root;
    while (qh < qt) {
        int u = queue[qh++];
        for (i = 0; i < deg[u]; i++) {
            int v = nbr[u][i];
            if (v >= root && distr[v] == 127) {
                distr[v] = distr[u] + 1;
                queue[qt++] = v;
            }
        }
    }
}

/* does the graph contain a cycle of length exactly k? */
static int has_cycle_len(int k)
{
    int root;
    if (k < 3 || k > N) return 0;
    for (root = 0; root < N; root++) {
        bfs_from(root);
        if (dfs_k(root, 1ULL << root, 0, k, root))
            return 1;
    }
    return 0;
}

static int girth(void)
{
    int k;
    for (k = 3; k <= N; k++)
        if (has_cycle_len(k))
            return k;
    return 0;                    /* forest (cannot happen for cubic) */
}

int main(int argc, char **argv)
{
    char line[256];
    int verbose = 0;
    long total = 0;
    long girth_hist[MAXN + 1];
    long combo[2][2][2];         /* [has4][has8][has16] */
    int i, j, k;

    for (i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-v")) verbose = 1;

    memset(girth_hist, 0, sizeof girth_hist);
    memset(combo, 0, sizeof combo);

    while (fgets(line, sizeof line, stdin)) {
        int g, c4, c8, c16;
        char *nl = strchr(line, '\n');
        if (nl) *nl = 0;
        if (!line[0]) continue;
        if (!parse_graph6(line)) {
            fprintf(stderr, "bad graph6: %s\n", line);
            return 1;
        }
        total++;
        g = girth();
        girth_hist[g]++;
        c4 = (g == 4) ? 1 : (g > 4 ? 0 : has_cycle_len(4));
        c8 = (g == 8) ? 1 : (g > 8 ? 0 : has_cycle_len(8));
        c16 = (g == 16) ? 1 : (g > 16 ? 0 : has_cycle_len(16));
        combo[c4][c8][c16]++;
        if (verbose)
            printf("V %s girth=%d c4=%d c8=%d c16=%d\n", line, g, c4, c8, c16);
        if (!c8)
            printf("SURV8 %s girth=%d c4=%d has16=%d\n", line, g, c4, c16);
        if (!c4 && c8 && !c16)
            printf("P8ONLY %s girth=%d\n", line, g);
        if (!c4 && !c8 && !c16)
            printf("CAND %s\n", line);
    }

    printf("SUMMARY {\"total\": %ld, \"girth_hist\": {", total);
    j = 0;
    for (k = 3; k <= MAXN; k++)
        if (girth_hist[k])
            printf("%s\"%d\": %ld", j++ ? ", " : "", k, girth_hist[k]);
    printf("}, \"combo_c4_c8_c16\": {");
    j = 0;
    for (i = 0; i < 2; i++)
        for (k = 0; k < 2; k++)
            for (int m = 0; m < 2; m++)
                if (combo[i][k][m])
                    printf("%s\"%d%d%d\": %ld", j++ ? ", " : "",
                           i, k, m, combo[i][k][m]);
    printf("}}\n");
    fflush(stdout);
    return 0;
}
