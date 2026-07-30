/* SKEPTIC re-enumeration of the canonical translation-word universe.
 *
 * Independent re-implementation for the adversarial review: DFS over all
 * bounce words up to MAXLEN, tracking the linear part of the composed
 * unfolding isometry as integer triples (p, q, r) meaning p*alpha + q*beta
 * + r*pi.  A word is emitted when:
 *   (W1) no two cyclically adjacent letters equal,
 *   (W2) even length,
 *   (W3) linear part is a rotation with p = q = 0, r even,
 *   (W4) it is the lexicographic minimum over all rotations of itself and
 *        of its reversal,
 *   (W5) it is not u^k (k >= 2) with u itself satisfying (W1)-(W3).
 *
 * Derivation used for (W3), done by hand from the reflection composition:
 * reflecting a direction d in a mirror of direction t gives 2t - d; the
 * composed isometry alternates reflection/rotation, and its angle after
 * each step is 2t - angle.  Initial side directions: side 2 = AB has
 * direction 0; side 1 = CA has direction alpha (leaves A at angle alpha);
 * side 0 = BC has direction pi - beta.  These match the geometry of the
 * normalized triangle, and the whole scheme is cross-checked against the
 * tier-0 harness by skeptic_words_xcheck.py at short lengths.
 *
 * Usage: skeptic_enum MAXLEN OUTDIR   (writes words_LNN.txt per length)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXLEN 32

static int maxlen;
static long long counts[MAXLEN + 1];
static FILE *outf[MAXLEN + 1];
static char outdir[512];

typedef struct { long long p, q, r; } tri;

static int word[MAXLEN];

/* linear part of the composed isometry for word w of length n:
 * returns kind (0 = rotation, 1 = reflection) and angle in *ang. */
static int lin_part(const int *w, int n, tri *ang)
{
    tri dirs[3] = { {0, -1, 1}, {1, 0, 0}, {0, 0, 0} };
    tri a = {0, 0, 0};
    int kind = 0;
    for (int i = 0; i < n; i++) {
        tri t = dirs[w[i]];
        a.p = 2 * t.p - a.p; a.q = 2 * t.q - a.q; a.r = 2 * t.r - a.r;
        kind ^= 1;
        for (int s = 0; s < 3; s++) {
            dirs[s].p = 2 * t.p - dirs[s].p;
            dirs[s].q = 2 * t.q - dirs[s].q;
            dirs[s].r = 2 * t.r - dirs[s].r;
        }
    }
    *ang = a;
    return kind;
}

static int is_translation(const int *w, int n)
{
    tri a;
    if (n % 2) return 0;                    /* reflection, never */
    if (lin_part(w, n, &a)) return 0;
    return a.p == 0 && a.q == 0 && (a.r % 2) == 0;
}

/* lexicographic compare of rotation r of base vs current best word */
static int canonical(const int *w, int n)
{
    int rev[MAXLEN];
    for (int i = 0; i < n; i++) rev[i] = w[n - 1 - i];
    const int *bases[2] = { w, rev };
    for (int b = 0; b < 2; b++) {
        for (int r = 0; r < n; r++) {
            if (b == 0 && r == 0) continue;
            /* compare bases[b] rotated by r against w */
            for (int i = 0; i < n; i++) {
                int c = bases[b][(i + r) % n];
                if (c < w[i]) return 0;     /* something smaller exists */
                if (c > w[i]) break;
            }
        }
    }
    return 1;
}

static int power_of_translation(const int *w, int n)
{
    for (int p = 2; p < n; p++) {
        if (n % p) continue;
        int rep = 1;
        for (int i = p; i < n && rep; i++)
            if (w[i] != w[i % p]) rep = 0;
        if (!rep) continue;
        if (p % 2) continue;                /* odd u never a translation */
        if (w[0] == w[p - 1]) continue;     /* u not cyclically valid */
        if (is_translation(w, p)) return 1;
    }
    return 0;
}

static void emit(int n)
{
    if (!canonical(word, n)) return;
    if (power_of_translation(word, n)) return;
    counts[n]++;
    if (!outf[n]) {
        char path[600];
        snprintf(path, sizeof path, "%s/words_L%02d.txt", outdir, n);
        outf[n] = fopen(path, "w");
        if (!outf[n]) { perror(path); exit(1); }
    }
    for (int i = 0; i < n; i++) fputc('0' + word[i], outf[n]);
    fputc('\n', outf[n]);
}

static void dfs(int depth, tri ang, tri d0, tri d1, tri d2)
{
    if (depth >= 2 && depth % 2 == 0 && word[0] != word[depth - 1]
        && ang.p == 0 && ang.q == 0 && (ang.r % 2) == 0)
        emit(depth);
    if (depth == maxlen) return;
    int prev = depth ? word[depth - 1] : -1;
    tri dirs[3] = { d0, d1, d2 };
    for (int s = 0; s < 3; s++) {
        if (s == prev) continue;
        tri t = dirs[s];
        tri na = { 2 * t.p - ang.p, 2 * t.q - ang.q, 2 * t.r - ang.r };
        tri nd[3];
        for (int k = 0; k < 3; k++) {
            nd[k].p = 2 * t.p - dirs[k].p;
            nd[k].q = 2 * t.q - dirs[k].q;
            nd[k].r = 2 * t.r - dirs[k].r;
        }
        word[depth] = s;
        dfs(depth + 1, na, nd[0], nd[1], nd[2]);
    }
}

int main(int argc, char **argv)
{
    if (argc != 3) { fprintf(stderr, "usage: %s MAXLEN OUTDIR\n", argv[0]); return 2; }
    maxlen = atoi(argv[1]);
    if (maxlen < 2 || maxlen > MAXLEN) { fprintf(stderr, "bad MAXLEN\n"); return 2; }
    snprintf(outdir, sizeof outdir, "%s", argv[2]);
    tri ang = {0, 0, 0};
    tri d0 = {0, -1, 1}, d1 = {1, 0, 0}, d2 = {0, 0, 0};
    dfs(0, ang, d0, d1, d2);
    long long total = 0;
    for (int n = 2; n <= maxlen; n++) {
        if (counts[n]) printf("%d: %lld\n", n, counts[n]);
        total += counts[n];
        if (outf[n]) fclose(outf[n]);
    }
    printf("total: %lld\n", total);
    return 0;
}
