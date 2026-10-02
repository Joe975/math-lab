/* Multi-start local minimisation of the MUB defect
 *
 *   L(B_1..B_k) = sum_{a<b} sum_{i,j} (|<a_i|b_j>|^2 - 1/d)^2
 *
 * over k orthonormal bases of C^d, with B_1 fixed to the standard basis.
 * L = 0 iff the bases are mutually unbiased.
 *
 * Method: Riemannian-style projected gradient on each unitary U_b
 * (columns = basis vectors), retraction by modified Gram-Schmidt,
 * backtracking step control. Deterministic given the seed.
 *
 *   cc -O2 -o mub_search mub_search.c -lm
 *   ./mub_search d k starts seed iters [best.json]
 *
 * stdout: one line per start "start L", then a summary line "# best L".
 */
#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXD 8
#define MAXK 9
typedef double complex cx;

static int D, K;
static cx U[MAXK][MAXD][MAXD];   /* U[b][row][col]; U[0] = identity */
static cx G[MAXK][MAXD][MAXD];   /* Euclidean gradient wrt conj(U) */

static uint64_t rng_s;
static double urand(void) {      /* splitmix64 */
    uint64_t z = (rng_s += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return (z >> 11) * (1.0 / 9007199254740992.0);
}
static double gauss(void) {
    double u = urand(), v = urand();
    if (u < 1e-300) u = 1e-300;
    return sqrt(-2 * log(u)) * cos(2 * M_PI * v);
}

static void gram_schmidt(cx M[MAXD][MAXD]) {
    for (int c = 0; c < D; c++) {
        for (int p = 0; p < c; p++) {
            cx ip = 0;
            for (int r = 0; r < D; r++) ip += conj(M[r][p]) * M[r][c];
            for (int r = 0; r < D; r++) M[r][c] -= ip * M[r][p];
        }
        double n = 0;
        for (int r = 0; r < D; r++) n += creal(M[r][c] * conj(M[r][c]));
        n = sqrt(n);
        for (int r = 0; r < D; r++) M[r][c] /= n;
    }
}

/* loss, and if grad != 0 the Euclidean gradient into G */
static double loss(cx V[MAXK][MAXD][MAXD], int grad) {
    double L = 0, inv = 1.0 / D;
    if (grad) memset(G, 0, sizeof G);
    for (int a = 0; a < K; a++)
        for (int b = a + 1; b < K; b++) {
            cx R[MAXD][MAXD];
            for (int i = 0; i < D; i++)
                for (int j = 0; j < D; j++) {
                    cx g = 0;
                    for (int r = 0; r < D; r++) g += conj(V[a][r][i]) * V[b][r][j];
                    double e = creal(g * conj(g)) - inv;
                    L += e * e;
                    R[i][j] = 2 * e * g;   /* dL/d conj(g) */
                }
            if (!grad) continue;
            /* g_ij = sum_r conj(Va_ri) Vb_rj
             * dL/dconj(Vb_rj) = sum_i Va_ri R_ij ; dL/dconj(Va_ri) = sum_j Vb_rj conj(R_ij) */
            for (int r = 0; r < D; r++)
                for (int x = 0; x < D; x++) {
                    cx sb = 0, sa = 0;
                    for (int y = 0; y < D; y++) {
                        sb += V[a][r][y] * R[y][x];
                        sa += V[b][r][y] * conj(R[x][y]);
                    }
                    G[b][r][x] += sb;
                    G[a][r][x] += sa;
                }
        }
    return L;
}

static double optimise(int iters, double *gnorm_out) {
    static cx T[MAXK][MAXD][MAXD];
    double eta = 0.1, L = loss(U, 1), gn = 0;
    for (int it = 0; it < iters; it++) {
        /* tangent projection: P = G - U herm(U^H G) */
        static cx P[MAXK][MAXD][MAXD];
        gn = 0;
        for (int b = 1; b < K; b++) {
            cx A[MAXD][MAXD], H[MAXD][MAXD];
            for (int i = 0; i < D; i++)
                for (int j = 0; j < D; j++) {
                    cx s = 0;
                    for (int r = 0; r < D; r++) s += conj(U[b][r][i]) * G[b][r][j];
                    A[i][j] = s;
                }
            for (int i = 0; i < D; i++)
                for (int j = 0; j < D; j++) H[i][j] = 0.5 * (A[i][j] + conj(A[j][i]));
            for (int r = 0; r < D; r++)
                for (int j = 0; j < D; j++) {
                    cx s = 0;
                    for (int i = 0; i < D; i++) s += U[b][r][i] * H[i][j];
                    P[b][r][j] = G[b][r][j] - s;
                    gn += creal(P[b][r][j] * conj(P[b][r][j]));
                }
        }
        gn = sqrt(gn);
        if (gn < 1e-14 || L < 1e-28) break;
        for (;;) {
            memcpy(T, U, sizeof U);
            for (int b = 1; b < K; b++) {
                for (int r = 0; r < D; r++)
                    for (int j = 0; j < D; j++) T[b][r][j] -= eta * P[b][r][j];
                gram_schmidt(T[b]);
            }
            double Lt = loss(T, 0);
            if (Lt < L) { memcpy(U, T, sizeof U); L = Lt; eta *= 1.3; break; }
            eta *= 0.5;
            if (eta < 1e-18) { *gnorm_out = gn; return L; }
        }
        L = loss(U, 1);
    }
    *gnorm_out = gn;
    return L;
}

static void dump(const char *path, double L) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); exit(1); }
    fprintf(f, "{\"d\": %d, \"k\": %d, \"L\": %.17g,\n \"bases\": [\n", D, K, L);
    for (int b = 0; b < K; b++) {
        fprintf(f, "  [");
        for (int c = 0; c < D; c++) {
            fprintf(f, "%s[", c ? ", " : "");
            for (int r = 0; r < D; r++)
                fprintf(f, "%s[%.17g, %.17g]", r ? ", " : "", creal(U[b][r][c]), cimag(U[b][r][c]));
            fprintf(f, "]");
        }
        fprintf(f, "]%s\n", b + 1 < K ? "," : "");
    }
    fprintf(f, " ],\n \"convention\": \"bases[b][c] is basis vector c of basis b, entries [re, im]\"}\n");
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc < 6) {
        fprintf(stderr, "usage: %s d k starts seed iters [best.json]\n", argv[0]);
        return 2;
    }
    D = atoi(argv[1]); K = atoi(argv[2]);
    int starts = atoi(argv[3]), iters = atoi(argv[5]);
    rng_s = strtoull(argv[4], 0, 10);
    if (D < 2 || D > MAXD || K < 2 || K > MAXK) { fprintf(stderr, "bad d/k\n"); return 2; }
    static cx best[MAXK][MAXD][MAXD];
    double bestL = INFINITY;
    for (int s = 0; s < starts; s++) {
        memset(U, 0, sizeof U);
        for (int i = 0; i < D; i++) U[0][i][i] = 1;
        for (int b = 1; b < K; b++) {
            for (int r = 0; r < D; r++)
                for (int c = 0; c < D; c++) U[b][r][c] = gauss() + I * gauss();
            gram_schmidt(U[b]);
        }
        double gn, L = optimise(iters, &gn);
        printf("%d %.17g %.3g\n", s, L, gn);
        if (L < bestL) { bestL = L; memcpy(best, U, sizeof U); }
    }
    memcpy(U, best, sizeof U);
    printf("# d=%d k=%d starts=%d best %.17g\n", D, K, starts, bestL);
    if (argc > 6) dump(argv[6], bestL);
    return 0;
}
