/* census.c -- census of centrally symmetric polytopes with vertices in
 * {0,+-1}^4, enumerated as subsets of the 40 antipodal pairs of nonzero
 * points, up to the hyperoctahedral group B4 (signed coordinate
 * permutations, order 384).
 *
 * A subset of pairs is stored as a 40-bit mask over a FIXED pair ordering:
 * iterate p = (c0,c1,c2,c3), each ci in {-1,0,1} ascending, c0 outermost;
 * skip 0; keep p iff its first nonzero coordinate is +1.  That yields 40
 * representatives, indexed 0..39 in order of discovery.  The canonical
 * form of a mask is the minimum, as a uint64, over the 384 group images.
 *
 * Modes:
 *   census points                   print the 40 pair representatives
 *   census seed OUT                 canonical level-1 masks (the 4 orbits)
 *   census expand IN OUT            canonical level-(k+1) masks from level-k
 *   census canon                    stdin hex masks -> canonical hex masks
 *   census eval IN OUT START CNT    evaluate masks[START..START+CNT) of IN
 *
 * eval output, one line per mask:
 *   <hex> k=<pairs> st=<ok|degen|improper|punt> nv=- nf=- vol=- pvol=- P=-
 *
 *   degen    rank of the point set < 4 (not a body): no volume product
 *   improper some point of the mask is not a vertex of the hull; the hull
 *            equals that of a smaller mask, which the census reaches at its
 *            own level, so the mask is skipped rather than double-counted
 *   punt     an int64 overflow guard tripped; must be recomputed exactly
 *            (fractions) by the Python side.  Never silently approximated.
 *
 * Arithmetic contract: every combinatorial decision (facet identification,
 * sidedness, vertex tests, incidences) is EXACT int64 integer arithmetic
 * with overflow guards that downgrade the body to "punt".  Only the final
 * volume accumulation is double (error ~1e-13 relative, used purely as a
 * filter; every near-minimizer is re-verified exactly in Q by exact_check.py).
 *
 * Volume: divergence theorem recursing on dimension, as in the tier-0
 * harness: vol_d = (1/d) sum_F b_F * vol_{d-1}(proj_j F) / |a_{F,j}| with
 * a integer-primitive outward normals and j a coordinate with a_j != 0.
 * Polar volume needs no second vertex enumeration: the vertices of K* are
 * a_F/b_F over facets F of K, and the facet of K* supported by <v,y>=1
 * (v a vertex of K) has vertex set {a_F/b_F : F contains v}; those rational
 * points are scaled by L = lcm(b_F) to integers, guarded.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>

typedef int64_t i64;
typedef uint64_t u64;

/* ------------------------------------------------------------------ */
/* pair representatives and the B4 action                              */
/* ------------------------------------------------------------------ */

static int pairpt[40][4];      /* the 40 representatives */
static int code_of[81];        /* base-3 code -> pair index, -1 if none */
static int NPAIR = 0;

static int enc(const int *p) {
    return (p[0]+1) + 3*(p[1]+1) + 9*(p[2]+1) + 27*(p[3]+1);
}

static void build_pairs(void) {
    memset(code_of, -1, sizeof code_of);
    for (int c0 = -1; c0 <= 1; c0++)
    for (int c1 = -1; c1 <= 1; c1++)
    for (int c2 = -1; c2 <= 1; c2++)
    for (int c3 = -1; c3 <= 1; c3++) {
        int p[4] = {c0, c1, c2, c3};
        int fz = -1;
        for (int i = 0; i < 4; i++) if (p[i]) { fz = i; break; }
        if (fz < 0 || p[fz] != 1) continue;
        memcpy(pairpt[NPAIR], p, sizeof p);
        code_of[enc(p)] = NPAIR;
        NPAIR++;
    }
    if (NPAIR != 40) { fprintf(stderr, "pair count %d != 40\n", NPAIR); exit(2); }
}

static int pair_index(const int *p) {   /* normalize sign, look up */
    int q[4]; int fz = -1;
    for (int i = 0; i < 4; i++) if (p[i]) { fz = i; break; }
    if (fz < 0) return -1;
    for (int i = 0; i < 4; i++) q[i] = (p[fz] < 0) ? -p[i] : p[i];
    return code_of[enc(q)];
}

static int NG = 0;
static unsigned char ptab[384][40];    /* group element -> perm of pairs */

static void build_group(void) {
    int perm[24][4], np = 0;
    for (int a = 0; a < 4; a++) for (int b = 0; b < 4; b++)
    for (int c = 0; c < 4; c++) for (int d = 0; d < 4; d++) {
        if (a==b||a==c||a==d||b==c||b==d||c==d) continue;
        perm[np][0]=a; perm[np][1]=b; perm[np][2]=c; perm[np][3]=d; np++;
    }
    for (int pi = 0; pi < np; pi++)
    for (int s = 0; s < 16; s++) {
        for (int i = 0; i < 40; i++) {
            int q[4];
            for (int c = 0; c < 4; c++) {
                int v = pairpt[i][perm[pi][c]];
                q[c] = (s >> c & 1) ? -v : v;
            }
            int j = pair_index(q);
            if (j < 0) { fprintf(stderr, "group action broke\n"); exit(2); }
            ptab[NG][i] = (unsigned char)j;
        }
        NG++;
    }
    if (NG != 384) { fprintf(stderr, "group order %d != 384\n", NG); exit(2); }
}

static u64 canon(u64 m) {
    u64 best = m;
    int idx[12], k = 0;
    for (int i = 0; i < 40; i++) if (m >> i & 1) idx[k++] = i;
    for (int g = 0; g < NG; g++) {
        u64 t = 0;
        for (int j = 0; j < k; j++) t |= 1ULL << ptab[g][idx[j]];
        if (t < best) best = t;
    }
    return best;
}

/* ------------------------------------------------------------------ */
/* mask file IO                                                        */
/* ------------------------------------------------------------------ */

static u64 *read_masks(const char *path, long *n_out) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); exit(2); }
    long cap = 1 << 16, n = 0;
    u64 *v = malloc(cap * sizeof *v);
    char line[64];
    while (fgets(line, sizeof line, f)) {
        if (line[0] == '\n' || line[0] == '#') continue;
        if (n == cap) { cap *= 2; v = realloc(v, cap * sizeof *v); }
        v[n++] = strtoull(line, NULL, 16);
    }
    fclose(f);
    *n_out = n;
    return v;
}

static int cmp_u64(const void *a, const void *b) {
    u64 x = *(const u64 *)a, y = *(const u64 *)b;
    return x < y ? -1 : x > y;
}

/* ------------------------------------------------------------------ */
/* exact integer helpers                                               */
/* ------------------------------------------------------------------ */

static i64 gcd64(i64 a, i64 b) {
    if (a < 0) a = -a; if (b < 0) b = -b;
    while (b) { i64 t = a % b; a = b; b = t; }
    return a;
}

/* rank of an n x 4 integer matrix, fraction-free elimination.
 * Entries stay tiny here (inputs are 0,+-1), no guard needed. */
static int rank4(i64 m[][4], int n) {
    int r = 0;
    for (int c = 0; c < 4 && r < n; c++) {
        int p = -1;
        for (int i = r; i < n; i++) if (m[i][c]) { p = i; break; }
        if (p < 0) continue;
        for (int j = 0; j < 4; j++) { i64 t = m[r][j]; m[r][j] = m[p][j]; m[p][j] = t; }
        for (int i = r + 1; i < n; i++) {
            if (!m[i][c]) continue;
            i64 g = gcd64(m[i][c], m[r][c]);
            i64 fr = m[i][c] / g, fp = m[r][c] / g;
            for (int j = 0; j < 4; j++) m[i][j] = m[i][j] * fp - m[r][j] * fr;
        }
        r++;
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* generic brute-force hull volume, exact combinatorics, in 1..3 dims  */
/* ------------------------------------------------------------------ */

#define MAXP 512            /* max points passed to vol3/vol2 */
static int punt_flag;       /* set when an overflow guard trips */

#define GUARD ((i64)3000000000000000000LL)
static inline i64 gmul(i64 a, i64 b) {
    if (a == 0 || b == 0) return 0;
    i64 aa = a < 0 ? -a : a, bb = b < 0 ? -b : b;
    if (aa > GUARD / bb) { punt_flag = 1; return 0; }
    return a * b;
}

static double vol1(i64 p[][1], int n) {
    i64 lo = p[0][0], hi = p[0][0];
    for (int i = 1; i < n; i++) {
        if (p[i][0] < lo) lo = p[i][0];
        if (p[i][0] > hi) hi = p[i][0];
    }
    return (double)(hi - lo);
}

static double vol2(i64 p[][2], int n) {
    /* edges: pairs whose line has all points weakly on one side */
    i64 fa[128][2], fb[128]; int nf = 0;
    double acc = 0;
    for (int i = 0; i < n && !punt_flag; i++)
    for (int j = i + 1; j < n; j++) {
        i64 dx = p[j][0] - p[i][0], dy = p[j][1] - p[i][1];
        i64 a0 = dy, a1 = -dx;           /* normal to the edge */
        if (!a0 && !a1) continue;
        i64 b = gmul(a0, p[i][0]) + gmul(a1, p[i][1]);
        int pos = 0, neg = 0;
        for (int t = 0; t < n; t++) {
            i64 s = gmul(a0, p[t][0]) + gmul(a1, p[t][1]) - b;
            if (s > 0) pos = 1; else if (s < 0) neg = 1;
            if (pos && neg) break;
        }
        if (pos && neg) continue;
        if (pos) { a0 = -a0; a1 = -a1; b = -b; }
        i64 g = gcd64(a0, a1);
        a0 /= g; a1 /= g; b /= g;
        int dup = 0;
        for (int t = 0; t < nf; t++)
            if (fa[t][0]==a0 && fa[t][1]==a1 && fb[t]==b) { dup = 1; break; }
        if (dup) continue;
        if (nf == 128) { punt_flag = 1; break; }
        fa[nf][0]=a0; fa[nf][1]=a1; fb[nf]=b; nf++;
        int jc = a0 ? 0 : 1;
        i64 q[MAXP][1]; int nq = 0;
        for (int t = 0; t < n; t++)
            if (gmul(a0,p[t][0]) + gmul(a1,p[t][1]) == b)
                q[nq++][0] = p[t][1 - jc];
        i64 aj = jc ? a1 : a0; if (aj < 0) aj = -aj;
        acc += (double)b * vol1(q, nq) / (double)aj;
    }
    return acc / 2.0;
}

static double vol3(i64 p[][3], int n) {
    static i64 fa[2048][3], fb[2048]; int nf = 0;
    double acc = 0;
    if (n > MAXP) { punt_flag = 1; return 0; }
    for (int i = 0; i < n && !punt_flag; i++)
    for (int j = i + 1; j < n; j++)
    for (int k = j + 1; k < n; k++) {
        i64 d1[3], d2[3], a[3];
        for (int c = 0; c < 3; c++) { d1[c] = p[j][c]-p[i][c]; d2[c] = p[k][c]-p[i][c]; }
        a[0] = gmul(d1[1],d2[2]) - gmul(d1[2],d2[1]);
        a[1] = gmul(d1[2],d2[0]) - gmul(d1[0],d2[2]);
        a[2] = gmul(d1[0],d2[1]) - gmul(d1[1],d2[0]);
        if (!a[0] && !a[1] && !a[2]) continue;
        i64 b = gmul(a[0],p[i][0]) + gmul(a[1],p[i][1]) + gmul(a[2],p[i][2]);
        int pos = 0, neg = 0;
        for (int t = 0; t < n; t++) {
            i64 s = gmul(a[0],p[t][0]) + gmul(a[1],p[t][1]) + gmul(a[2],p[t][2]) - b;
            if (s > 0) pos = 1; else if (s < 0) neg = 1;
            if (pos && neg) break;
        }
        if (punt_flag) break;
        if (pos && neg) continue;
        if (pos) { a[0]=-a[0]; a[1]=-a[1]; a[2]=-a[2]; b=-b; }
        i64 g = gcd64(gcd64(a[0], a[1]), a[2]);
        a[0]/=g; a[1]/=g; a[2]/=g; b/=g;
        int dup = 0;
        for (int t = 0; t < nf; t++)
            if (fa[t][0]==a[0] && fa[t][1]==a[1] && fa[t][2]==a[2] && fb[t]==b) { dup = 1; break; }
        if (dup) continue;
        if (nf == 2048) { punt_flag = 1; break; }
        fa[nf][0]=a[0]; fa[nf][1]=a[1]; fa[nf][2]=a[2]; fb[nf]=b; nf++;
        int jc = a[0] ? 0 : (a[1] ? 1 : 2);
        i64 q[MAXP][2]; int nq = 0;
        for (int t = 0; t < n; t++)
            if (gmul(a[0],p[t][0]) + gmul(a[1],p[t][1]) + gmul(a[2],p[t][2]) == b) {
                int w = 0;
                for (int c = 0; c < 3; c++) if (c != jc) q[nq][w++] = p[t][c];
                nq++;
            }
        i64 aj = a[jc]; if (aj < 0) aj = -aj;
        acc += (double)b * vol2(q, nq) / (double)aj;
    }
    return acc / 3.0;
}

/* ------------------------------------------------------------------ */
/* the 4D body: facets, vertex test, volume, polar volume              */
/* ------------------------------------------------------------------ */

#define MAXV 40             /* 2 * max pairs evaluated */
#define MAXF 4096

static int npt;                     /* number of points (2k) */
static i64 pt[MAXV][4];             /* the points, entries 0,+-1 */
static int nfac;
static i64 fac_a[MAXF][4], fac_b[MAXF];
static u64 fac_on[MAXF];            /* incidence bitmask over points */

static int find_facets(void) {      /* returns 0 on internal error */
    nfac = 0;
    for (int i0 = 0; i0 < npt; i0++)
    for (int i1 = i0+1; i1 < npt; i1++)
    for (int i2 = i1+1; i2 < npt; i2++)
    for (int i3 = i2+1; i3 < npt; i3++) {
        i64 d[3][4], a[4];
        for (int c = 0; c < 4; c++) {
            d[0][c] = pt[i1][c] - pt[i0][c];
            d[1][c] = pt[i2][c] - pt[i0][c];
            d[2][c] = pt[i3][c] - pt[i0][c];
        }
        /* generalized cross product: a[l] = (-1)^l det(minor l) */
        for (int l = 0; l < 4; l++) {
            int cols[3], w = 0;
            for (int c = 0; c < 4; c++) if (c != l) cols[w++] = c;
            i64 det =
              d[0][cols[0]] * (d[1][cols[1]]*d[2][cols[2]] - d[1][cols[2]]*d[2][cols[1]])
            - d[0][cols[1]] * (d[1][cols[0]]*d[2][cols[2]] - d[1][cols[2]]*d[2][cols[0]])
            + d[0][cols[2]] * (d[1][cols[0]]*d[2][cols[1]] - d[1][cols[1]]*d[2][cols[0]]);
            a[l] = (l & 1) ? -det : det;
        }
        if (!a[0] && !a[1] && !a[2] && !a[3]) continue;
        i64 b = 0;
        for (int c = 0; c < 4; c++) b += a[c] * pt[i0][c];
        int pos = 0, neg = 0;
        for (int t = 0; t < npt; t++) {
            i64 s = -b;
            for (int c = 0; c < 4; c++) s += a[c] * pt[t][c];
            if (s > 0) pos = 1; else if (s < 0) neg = 1;
            if (pos && neg) break;
        }
        if (pos && neg) continue;
        if (pos) { for (int c = 0; c < 4; c++) a[c] = -a[c]; b = -b; }
        i64 g = gcd64(gcd64(a[0],a[1]), gcd64(a[2],a[3]));
        for (int c = 0; c < 4; c++) a[c] /= g;
        b /= g;
        int dup = 0;
        for (int t = 0; t < nfac; t++)
            if (fac_b[t]==b && fac_a[t][0]==a[0] && fac_a[t][1]==a[1]
                && fac_a[t][2]==a[2] && fac_a[t][3]==a[3]) { dup = 1; break; }
        if (dup) continue;
        if (b <= 0) return 0;       /* impossible for a symmetric body */
        if (nfac == MAXF) return 0;
        u64 on = 0;
        for (int t = 0; t < npt; t++) {
            i64 sv = -b;
            for (int c = 0; c < 4; c++) sv += a[c] * pt[t][c];
            if (sv == 0) on |= 1ULL << t;
        }
        for (int c = 0; c < 4; c++) fac_a[nfac][c] = a[c];
        fac_b[nfac] = b;
        fac_on[nfac] = on;
        nfac++;
    }
    return nfac >= 5;               /* a 4-polytope has >= 5 facets */
}

static int all_points_are_vertices(void) {
    for (int t = 0; t < npt; t++) {
        i64 m[MAXF][4]; int nm = 0;
        for (int f = 0; f < nfac; f++)
            if (fac_on[f] >> t & 1) {
                for (int c = 0; c < 4; c++) m[nm][c] = fac_a[f][c];
                nm++;
            }
        if (nm < 4 || rank4(m, nm) < 4) return 0;
    }
    return 1;
}

static double volume4(void) {
    double acc = 0;
    for (int f = 0; f < nfac; f++) {
        int jc = 0;
        while (!fac_a[f][jc]) jc++;
        i64 q[MAXV][3]; int nq = 0;
        for (int t = 0; t < npt; t++)
            if (fac_on[f] >> t & 1) {
                int w = 0;
                for (int c = 0; c < 4; c++) if (c != jc) q[nq][w++] = pt[t][c];
                nq++;
            }
        i64 aj = fac_a[f][jc]; if (aj < 0) aj = -aj;
        acc += (double)fac_b[f] * vol3(q, nq) / (double)aj;
    }
    return acc / 4.0;
}

/* vol(K*) = (1/4) sum over vertices v of vol3(proj_j {a_F/b_F : F>=v}) / |v_j|,
 * v_j = +-1 here.  Rational points scaled to integers by L = lcm(b_F). */
static double polar_volume4(void) {
    double acc = 0;
    for (int t = 0; t < npt && !punt_flag; t++) {
        int jc = 0;
        while (!pt[t][jc]) jc++;
        i64 L = 1;
        for (int f = 0; f < nfac; f++)
            if (fac_on[f] >> t & 1) {
                i64 b = fac_b[f], g = gcd64(L, b);
                L = gmul(L / g, b);
                if (punt_flag || L > 1000000) { punt_flag = 1; return 0; }
            }
        i64 q[MAXP][3]; int nq = 0;
        for (int f = 0; f < nfac; f++)
            if (fac_on[f] >> t & 1) {
                if (nq == MAXP) { punt_flag = 1; return 0; }
                i64 s = L / fac_b[f];
                int w = 0;
                for (int c = 0; c < 4; c++)
                    if (c != jc) {
                        q[nq][w] = gmul(fac_a[f][c], s);
                        if (q[nq][w] > 500000 || q[nq][w] < -500000) { punt_flag = 1; return 0; }
                        w++;
                    }
                nq++;
            }
        double v3 = vol3(q, nq);
        if (punt_flag) return 0;
        acc += v3 / ((double)L * (double)L * (double)L);
    }
    return acc / 4.0;
}

/* ------------------------------------------------------------------ */
/* modes                                                               */
/* ------------------------------------------------------------------ */

static void mode_points(void) {
    for (int i = 0; i < 40; i++)
        printf("%2d  %2d %2d %2d %2d\n", i,
               pairpt[i][0], pairpt[i][1], pairpt[i][2], pairpt[i][3]);
}

static void mode_seed(const char *out) {
    u64 v[40]; long n = 0;
    for (int i = 0; i < 40; i++) v[n++] = canon(1ULL << i);
    qsort(v, n, sizeof *v, cmp_u64);
    long m = 0;
    FILE *f = fopen(out, "w");
    for (long i = 0; i < n; i++)
        if (i == 0 || v[i] != v[i-1]) { fprintf(f, "%" PRIx64 "\n", v[i]); m++; }
    fclose(f);
    fprintf(stderr, "seed: %ld orbit(s)\n", m);
}

static void mode_expand(const char *in, const char *out) {
    long n; u64 *src = read_masks(in, &n);
    if (!n) { fprintf(stderr, "empty input\n"); exit(2); }
    int k = __builtin_popcountll(src[0]);
    long cap = n * (40 - k), m = 0;
    u64 *dst = malloc(cap * sizeof *dst);
    if (!dst) { perror("malloc"); exit(2); }
    for (long i = 0; i < n; i++)
        for (int j = 0; j < 40; j++)
            if (!(src[i] >> j & 1)) dst[m++] = canon(src[i] | 1ULL << j);
    qsort(dst, m, sizeof *dst, cmp_u64);
    FILE *f = fopen(out, "w");
    long uniq = 0;
    for (long i = 0; i < m; i++)
        if (i == 0 || dst[i] != dst[i-1]) { fprintf(f, "%" PRIx64 "\n", dst[i]); uniq++; }
    fclose(f);
    fprintf(stderr, "expand k=%d -> k=%d: %ld reps -> %ld orbit reps\n",
            k, k + 1, n, uniq);
    free(src); free(dst);
}

static void mode_canon(void) {
    char line[64];
    while (fgets(line, sizeof line, stdin)) {
        if (line[0] == '\n' || line[0] == '#') continue;
        printf("%" PRIx64 "\n", canon(strtoull(line, NULL, 16)));
    }
}

static void mode_eval(const char *in, const char *out, long start, long cnt) {
    long n; u64 *v = read_masks(in, &n);
    if (start < 0 || start > n) { fprintf(stderr, "bad start\n"); exit(2); }
    long end = start + cnt; if (end > n) end = n;
    FILE *f = fopen(out, "w");
    if (!f) { perror(out); exit(2); }
    for (long i = start; i < end; i++) {
        u64 m = v[i];
        int k = __builtin_popcountll(m);
        npt = 0;
        i64 rk[12][4]; int nr = 0;
        for (int j = 0; j < 40; j++)
            if (m >> j & 1) {
                for (int c = 0; c < 4; c++) {
                    pt[npt][c]   =  pairpt[j][c];
                    pt[npt+1][c] = -pairpt[j][c];
                    rk[nr][c]    =  pairpt[j][c];
                }
                npt += 2; nr++;
            }
        const char *st; int nv = 0, nf = 0;
        double vol = 0, pvol = 0, P = 0;
        punt_flag = 0;
        if (rank4(rk, nr) < 4) st = "degen";
        else if (!find_facets()) st = "punt";      /* structural surprise */
        else if (!all_points_are_vertices()) { st = "improper"; nf = nfac; }
        else {
            nv = npt; nf = nfac;
            vol = volume4();
            if (!punt_flag) pvol = polar_volume4();
            if (punt_flag) st = "punt";
            else if (vol <= 0 || pvol <= 0) st = "punt";
            else { st = "ok"; P = vol * pvol; }
        }
        fprintf(f, "%" PRIx64 " k=%d st=%s nv=%d nf=%d vol=%.17g pvol=%.17g P=%.17g\n",
                m, k, st, nv, nf, vol, pvol, P);
    }
    fclose(f);
    fprintf(stderr, "eval: wrote %ld..%ld of %ld from %s\n", start, end, n, in);
    free(v);
}

int main(int argc, char **argv) {
    build_pairs();
    build_group();
    if (argc >= 2 && !strcmp(argv[1], "points")) { mode_points(); return 0; }
    if (argc >= 3 && !strcmp(argv[1], "seed"))   { mode_seed(argv[2]); return 0; }
    if (argc >= 4 && !strcmp(argv[1], "expand")) { mode_expand(argv[2], argv[3]); return 0; }
    if (argc >= 2 && !strcmp(argv[1], "canon"))  { mode_canon(); return 0; }
    if (argc >= 6 && !strcmp(argv[1], "eval"))
        { mode_eval(argv[2], argv[3], atol(argv[4]), atol(argv[5])); return 0; }
    fprintf(stderr, "usage: census points|seed OUT|expand IN OUT|canon|eval IN OUT START CNT\n");
    return 2;
}
