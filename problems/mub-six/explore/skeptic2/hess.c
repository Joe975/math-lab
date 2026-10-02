/* Independent Hessian analysis of 4 bases in C^6, basis 0 fixed.
 * L = sum_{a<b} sum_{i,j} (|<a_i|b_j>|^2 - 1/6)^2
 * Chart: M_k -> M_k C(A_k), C = Cayley (I-A/2)^{-1}(I+A/2), A anti-Hermitian,
 *   A = sum_m x_m T_m with Frobenius-orthonormal T_m:
 *   i E_jj (6), (E_kl - E_lk)/sqrt2 (15), i(E_kl + E_lk)/sqrt2 (15).
 * Hessian at x=0: exact via hyperdual numbers (forward-mode 2nd order AD),
 *   using C(A) = I + A + A^2/2 + O(A^3) (exact to 2nd order at A=0).
 * Cross-check: central second differences of L (true Cayley, long double), Richardson.
 * Polishing: Newton steps (pseudo-inverse dropping |lambda|<thr), retraction by true Cayley.
 * All arithmetic long double (64-bit mantissa).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
typedef long double R;
#define D 6
#define K 4
#define NP 36
#define N (3*NP)

typedef struct { R re, im; } C;
static C M[K][D][D]; /* M[k][row][col]; column = basis vector */

/* ---------- hyperdual ---------- */
typedef struct { R f, a, b, ab; } H;
typedef struct { H re, im; } HC;
static inline H hadd(H x, H y){ H r={x.f+y.f,x.a+y.a,x.b+y.b,x.ab+y.ab}; return r; }
static inline H hsub(H x, H y){ H r={x.f-y.f,x.a-y.a,x.b-y.b,x.ab-y.ab}; return r; }
static inline H hmul(H x, H y){ H r={x.f*y.f, x.f*y.a+x.a*y.f, x.f*y.b+x.b*y.f,
  x.f*y.ab+x.a*y.b+x.b*y.a+x.ab*y.f}; return r; }
static inline H hsc(H x, R s){ H r={x.f*s,x.a*s,x.b*s,x.ab*s}; return r; }
static inline HC hcadd(HC x, HC y){ HC r={hadd(x.re,y.re),hadd(x.im,y.im)}; return r; }
static inline HC hcmul(HC x, HC y){ HC r={hsub(hmul(x.re,y.re),hmul(x.im,y.im)),
  hadd(hmul(x.re,y.im),hmul(x.im,y.re))}; return r; }
static inline HC hcconjmul(HC x, HC y){ /* conj(x)*y */ HC r={hadd(hmul(x.re,y.re),hmul(x.im,y.im)),
  hsub(hmul(x.re,y.im),hmul(x.im,y.re))}; return r; }

/* generator m (0..35) -> anti-Hermitian matrix (complex, real) */
static void gen(int m, C T[D][D]){
  memset(T,0,sizeof(C)*D*D);
  R s=1.0L/sqrtl(2.0L);
  if(m<6){ T[m][m].im=1; return; }
  m-=6; int idx=0;
  for(int k=0;k<D;k++)for(int l=k+1;l<D;l++){
    if(idx==m){ T[k][l].re=s; T[l][k].re=-s; return; }
    if(idx+15==m){ T[k][l].im=s; T[l][k].im=s; return; }
    idx++;
  }
}
static C Tg[NP][D][D];

/* hyperdual L with x_i seeded in a, x_j seeded in b (i or j may be -1) */
static H L_hd(int i, int j){
  static HC U[K][D][D];
  for(int r=0;r<D;r++)for(int c=0;c<D;c++){ HC z; memset(&z,0,sizeof z);
    z.re.f=M[0][r][c].re; z.im.f=M[0][r][c].im; U[0][r][c]=z; }
  for(int k=1;k<K;k++){
    HC A[D][D]; memset(A,0,sizeof A);
    for(int p=0;p<2;p++){ int q = p? j : i; if(q<0) continue;
      if(q/NP != k-1) continue; int m=q%NP;
      for(int r=0;r<D;r++)for(int c=0;c<D;c++){
        if(p==0){ A[r][c].re.a+=Tg[m][r][c].re; A[r][c].im.a+=Tg[m][r][c].im; }
        else    { A[r][c].re.b+=Tg[m][r][c].re; A[r][c].im.b+=Tg[m][r][c].im; } } }
    /* G = I + A + A^2/2 */
    HC G[D][D];
    for(int r=0;r<D;r++)for(int c=0;c<D;c++){
      HC s; memset(&s,0,sizeof s);
      for(int t=0;t<D;t++) s=hcadd(s,hcmul(A[r][t],A[t][c]));
      s.re=hsc(s.re,0.5L); s.im=hsc(s.im,0.5L);
      s=hcadd(s,A[r][c]); if(r==c) s.re.f+=1; G[r][c]=s; }
    for(int r=0;r<D;r++)for(int c=0;c<D;c++){
      HC s; memset(&s,0,sizeof s);
      for(int t=0;t<D;t++){ HC m0; memset(&m0,0,sizeof m0); m0.re.f=M[k][r][t].re; m0.im.f=M[k][r][t].im;
        s=hcadd(s,hcmul(m0,G[t][c])); }
      U[k][r][c]=s; }
  }
  H L; memset(&L,0,sizeof L);
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)
    for(int p=0;p<D;p++)for(int q=0;q<D;q++){
      HC s; memset(&s,0,sizeof s);
      for(int r=0;r<D;r++) s=hcadd(s,hcconjmul(U[a][r][p],U[b][r][q]));
      H w=hadd(hmul(s.re,s.re),hmul(s.im,s.im)); w.f-=1.0L/D;
      L=hadd(L,hmul(w,w)); }
  return L;
}

/* ---------- plain long double ---------- */
static inline C cmul(C x, C y){ C r={x.re*y.re-x.im*y.im, x.re*y.im+x.im*y.re}; return r; }
static inline C cdiv(C x, C y){ R d=y.re*y.re+y.im*y.im; C r={(x.re*y.re+x.im*y.im)/d,(x.im*y.re-x.re*y.im)/d}; return r; }
/* Cayley: Q = (I - A/2)^{-1} (I + A/2) */
static void cayley(C A[D][D], C Q[D][D]){
  C P[D][D], Bm[D][D];
  for(int r=0;r<D;r++)for(int c=0;c<D;c++){
    P[r][c].re=(r==c)-A[r][c].re/2; P[r][c].im=-A[r][c].im/2;
    Bm[r][c].re=(r==c)+A[r][c].re/2; Bm[r][c].im=A[r][c].im/2; }
  for(int col=0;col<D;col++){ /* partial pivot */
    int pv=col; R best=-1;
    for(int r=col;r<D;r++){ R v=P[r][col].re*P[r][col].re+P[r][col].im*P[r][col].im; if(v>best){best=v;pv=r;} }
    if(pv!=col) for(int c=0;c<D;c++){ C t=P[col][c];P[col][c]=P[pv][c];P[pv][c]=t; t=Bm[col][c];Bm[col][c]=Bm[pv][c];Bm[pv][c]=t; }
    for(int r=0;r<D;r++) if(r!=col){
      C f=cdiv(P[r][col],P[col][col]);
      for(int c=0;c<D;c++){ C t=cmul(f,P[col][c]); P[r][c].re-=t.re; P[r][c].im-=t.im;
        t=cmul(f,Bm[col][c]); Bm[r][c].re-=t.re; Bm[r][c].im-=t.im; } } }
  for(int r=0;r<D;r++)for(int c=0;c<D;c++) Q[r][c]=cdiv(Bm[r][c],P[r][r]);
}
static void apply(const R *x, C out[K][D][D]){
  memcpy(out[0],M[0],sizeof(C)*D*D);
  for(int k=1;k<K;k++){
    C A[D][D]; memset(A,0,sizeof A);
    for(int m=0;m<NP;m++){ R v=x[(k-1)*NP+m]; if(v==0) continue;
      for(int r=0;r<D;r++)for(int c=0;c<D;c++){ A[r][c].re+=v*Tg[m][r][c].re; A[r][c].im+=v*Tg[m][r][c].im; } }
    C Q[D][D]; cayley(A,Q);
    for(int r=0;r<D;r++)for(int c=0;c<D;c++){ C s={0,0};
      for(int t=0;t<D;t++){ C u=cmul(M[k][r][t],Q[t][c]); s.re+=u.re; s.im+=u.im; } out[k][r][c]=s; } }
}
static R Lval(C U[K][D][D]){
  R L=0;
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int p=0;p<D;p++)for(int q=0;q<D;q++){
    R sr=0,si=0; for(int r=0;r<D;r++){ C x=U[a][r][p],y=U[b][r][q]; sr+=x.re*y.re+x.im*y.im; si+=x.re*y.im-x.im*y.re; }
    R w=sr*sr+si*si-1.0L/D; L+=w*w; }
  return L;
}
static R Lx(const R *x){ static C U[K][D][D]; apply(x,U); return Lval(U); }

/* Jacobi eigen for symmetric n x n */
static void jacobi(int n, R *A, R *ev, R *V){
  for(int i=0;i<n;i++)for(int j=0;j<n;j++) V[i*n+j]=(i==j);
  for(int sweep=0;sweep<100;sweep++){
    R off=0; for(int i=0;i<n;i++)for(int j=i+1;j<n;j++) off+=A[i*n+j]*A[i*n+j];
    if(off<1e-60L) break;
    for(int p=0;p<n;p++)for(int q=p+1;q<n;q++){
      R apq=A[p*n+q]; if(fabsl(apq)<1e-300L) continue;
      R th=(A[q*n+q]-A[p*n+p])/(2*apq);
      R t=(th>=0?1:-1)/(fabsl(th)+sqrtl(th*th+1)); R c=1/sqrtl(t*t+1), s=t*c;
      for(int k=0;k<n;k++){ R akp=A[k*n+p],akq=A[k*n+q]; A[k*n+p]=c*akp-s*akq; A[k*n+q]=s*akp+c*akq; }
      for(int k=0;k<n;k++){ R apk=A[p*n+k],aqk=A[q*n+k]; A[p*n+k]=c*apk-s*aqk; A[q*n+k]=s*apk+c*aqk; }
      for(int k=0;k<n;k++){ R vkp=V[k*n+p],vkq=V[k*n+q]; V[k*n+p]=c*vkp-s*vkq; V[k*n+q]=s*vkp+c*vkq; } } }
  for(int i=0;i<n;i++) ev[i]=A[i*n+i];
  /* sort ascending with vectors (columns) */
  for(int i=0;i<n;i++){ int b=i; for(int j=i+1;j<n;j++) if(ev[j]<ev[b]) b=j;
    if(b!=i){ R t=ev[i];ev[i]=ev[b];ev[b]=t; for(int k=0;k<n;k++){ t=V[k*n+i];V[k*n+i]=V[k*n+b];V[k*n+b]=t; } } }
}

static void hess_hd(R *Hm, R *g, R *L0){
  for(int i=0;i<N;i++)for(int j=i;j<N;j++){
    H r=L_hd(i,j); Hm[i*N+j]=Hm[j*N+i]=r.ab; if(i==j){ g[i]=r.a; *L0=r.f; } }
}

static void reorth(void){ /* modified Gram-Schmidt on columns of moving bases */
  for(int k=1;k<K;k++) for(int c=0;c<D;c++){
    for(int c2=0;c2<c;c2++){ C s={0,0};
      for(int r=0;r<D;r++){ s.re+=M[k][r][c2].re*M[k][r][c].re+M[k][r][c2].im*M[k][r][c].im;
        s.im+=M[k][r][c2].re*M[k][r][c].im-M[k][r][c2].im*M[k][r][c].re; }
      for(int r=0;r<D;r++){ C u=cmul(s,M[k][r][c2]); M[k][r][c].re-=u.re; M[k][r][c].im-=u.im; } }
    R nn=0; for(int r=0;r<D;r++) nn+=M[k][r][c].re*M[k][r][c].re+M[k][r][c].im*M[k][r][c].im;
    nn=sqrtl(nn); for(int r=0;r<D;r++){ M[k][r][c].re/=nn; M[k][r][c].im/=nn; } }
}

static R Hm[N*N], Hc[N*N], V[N*N], ev[N], g[N];

int main(int argc, char **argv){
  if(argc<2){ fprintf(stderr,"usage: hess in.json [out.json]\n"); return 1; }
  for(int m=0;m<NP;m++) gen(m,Tg[m]);
  /* parse: all numbers in order; layout [basis][col][row][re,im] */
  FILE *f=fopen(argv[1],"r"); fseek(f,0,SEEK_END); long n=ftell(f); rewind(f);
  char *buf=malloc(n+1); fread(buf,1,n,f); buf[n]=0; fclose(f);
  char *p=buf; int cnt=0; R vals[K*D*D*2];
  while(*p && cnt<K*D*D*2){ if((*p>='0'&&*p<='9')||*p=='-'||*p=='.'){ char *e; vals[cnt++]=strtold(p,&e); p=e; } else p++; }
  if(cnt!=K*D*D*2){ fprintf(stderr,"parse got %d\n",cnt); return 1; }
  for(int k=0;k<K;k++)for(int c=0;c<D;c++)for(int r=0;r<D;r++){
    int o=((k*D+c)*D+r)*2; M[k][r][c].re=vals[o]; M[k][r][c].im=vals[o+1]; }
  reorth();
  { R z[N]={0}; printf("L(start, after reorth) = %.19Lg\n", Lx(z)); }

  /* Newton polishing */
  R L0=0;
  for(int it=0; it<12; it++){
    hess_hd(Hm,g,&L0);
    R gn=0; for(int i=0;i<N;i++) gn+=g[i]*g[i]; gn=sqrtl(gn);
    memcpy(Hc,Hm,sizeof Hm); jacobi(N,Hc,ev,V);
    int nz=0; for(int i=0;i<N;i++) if(fabsl(ev[i])<1e-8L) nz++;
    printf("it %2d  L=%.19Lg  |g|=%.3Le  #|ev|<1e-8: %d  minev(nonzero-ish)=%.4Le\n",it,L0,gn,nz,ev[nz]);
    if(gn<1e-17L) break;
    R x[N]={0};
    for(int e=0;e<N;e++){ if(fabsl(ev[e])<1e-8L) continue; R c=0; for(int i=0;i<N;i++) c+=V[i*N+e]*g[i];
      for(int i=0;i<N;i++) x[i]-=c/ev[e]*V[i*N+e]; }
    static C U[K][D][D]; apply(x,U); memcpy(M,U,sizeof U); reorth();
  }
  hess_hd(Hm,g,&L0);
  { R gn=0; for(int i=0;i<N;i++) gn+=g[i]*g[i]; printf("FINAL L=%.19Lg |g|=%.3Le\n",L0,sqrtl(gn)); }

  /* symmetry generators in tangent coordinates */
  static R S[23][N]; int ns=0;
  for(int k=1;k<K;k++) for(int j=0;j<D;j++){ memset(S[ns],0,sizeof S[ns]); S[ns][(k-1)*NP+j]=1; ns++; }
  for(int j=0;j<D-1;j++){ /* left phase i E_jj on all bases (basis 0 = I: left phase == its own column phase, which is gauge) */
    memset(S[ns],0,sizeof S[ns]);
    for(int k=1;k<K;k++){ /* A = M^dag (i E_jj) M  : A[r][c] = i conj(M[j][r]) M[j][c] */
      C A[D][D]; for(int r=0;r<D;r++)for(int c=0;c<D;c++){ C a={M[k][j][r].re,-M[k][j][r].im}; C b=cmul(a,M[k][j][c]); A[r][c].re=-b.im; A[r][c].im=b.re; }
      for(int m=0;m<NP;m++){ R s=0; for(int r=0;r<D;r++)for(int c=0;c<D;c++) s+=Tg[m][r][c].re*A[r][c].re+Tg[m][r][c].im*A[r][c].im; S[ns][(k-1)*NP+m]=s; } }
    ns++; }
  /* Basis 0 is the identity? check */
  { R dev=0; for(int r=0;r<D;r++)for(int c=0;c<D;c++) dev+=fabsl(M[0][r][c].re-(r==c))+fabsl(M[0][r][c].im); printf("basis0 deviation from I: %.3Le\n",dev); }
  /* |H s| / |s| for each symmetry vector; rank of span via Gram-Schmidt */
  R maxres=0; static R Q[23][N]; int rank=0;
  for(int s=0;s<ns;s++){ R nn=0,hn=0; for(int i=0;i<N;i++){ nn+=S[s][i]*S[s][i]; R h=0; for(int j=0;j<N;j++) h+=Hm[i*N+j]*S[s][j]; hn+=h*h; }
    if(sqrtl(hn/nn)>maxres) maxres=sqrtl(hn/nn);
    R v[N]; memcpy(v,S[s],sizeof v);
    for(int t=0;t<rank;t++){ R d=0; for(int i=0;i<N;i++) d+=Q[t][i]*v[i]; for(int i=0;i<N;i++) v[i]-=d*Q[t][i]; }
    R vn=0; for(int i=0;i<N;i++) vn+=v[i]*v[i]; vn=sqrtl(vn);
    if(vn>1e-10L*sqrtl(nn)){ for(int i=0;i<N;i++) Q[rank][i]=v[i]/vn; rank++; } }
  printf("symmetry generators: %d, rank of span: %d, max |H s|/|s| = %.3Le\n",ns,rank,maxres);

  memcpy(Hc,Hm,sizeof Hm); jacobi(N,Hc,ev,V);
  printf("HD-Hessian spectrum (Frobenius-orthonormal Lie coords, A = sum x_m T_m, ||A||_F=|x|):\n");
  for(int i=0;i<N;i++){ printf("%3d % .10Le\n",i,ev[i]); }
  /* overlap of the near-null eigenvectors with symmetry span */
  { int nz=0; for(int i=0;i<N;i++) if(fabsl(ev[i])<1e-8L) nz++;
    R worst=1; for(int e=0;e<nz;e++){ R pr=0; for(int t=0;t<rank;t++){ R d=0; for(int i=0;i<N;i++) d+=Q[t][i]*V[i*N+e]; pr+=d*d; } if(pr<worst) worst=pr; }
    printf("#near-zero (|ev|<1e-8) = %d ; min squared projection of a null eigvec onto symmetry span = %.15Lf\n",nz,worst);
    /* Hessian restricted to orthogonal complement of symmetry span */
    static R P[N*N], Hr[N*N], Vr[N*N], evr[N];
    for(int i=0;i<N;i++)for(int j=0;j<N;j++){ R s=(i==j); for(int t=0;t<rank;t++) s-=Q[t][i]*Q[t][j]; P[i*N+j]=s; }
    static R T1[N*N];
    for(int i=0;i<N;i++)for(int j=0;j<N;j++){ R s=0; for(int k=0;k<N;k++) s+=P[i*N+k]*Hm[k*N+j]; T1[i*N+j]=s; }
    for(int i=0;i<N;i++)for(int j=0;j<N;j++){ R s=0; for(int k=0;k<N;k++) s+=T1[i*N+k]*P[k*N+j]; Hr[i*N+j]=s; }
    jacobi(N,Hr,evr,Vr);
    printf("projected-Hessian spectrum (symmetry span projected out), lowest 30 and top 3:\n");
    for(int i=0;i<30;i++) printf("  %3d % .10Le\n",i,evr[i]);
    for(int i=N-3;i<N;i++) printf("  %3d % .10Le\n",i,evr[i]);
  }

  /* finite-difference cross-check: central 2nd differences of L in true Cayley chart, Richardson */
  { R hs[2]={2e-3L,1e-3L}; static R F[2][N*N];
    R z[N]={0}; R l0=Lx(z);
    for(int s=0;s<2;s++){ R h=hs[s];
      for(int i=0;i<N;i++)for(int j=i;j<N;j++){ R x[N]={0}; R v;
        if(i==j){ x[i]=h; R lp=Lx(x); x[i]=-h; R lm=Lx(x); v=(lp-2*l0+lm)/(h*h); }
        else { x[i]=h;x[j]=h; R a=Lx(x); x[j]=-h; R b=Lx(x); x[i]=-h; R d=Lx(x); x[j]=h; R c=Lx(x); v=(a-b-c+d)/(4*h*h); }
        F[s][i*N+j]=F[s][j*N+i]=v; } }
    R md=0, mh=0, mdr=0; for(int i=0;i<N*N;i++){ R r=(4*F[1][i]-F[0][i])/3; R d=fabsl(r-Hm[i]); if(d>md) md=d; if(fabsl(Hm[i])>mh) mh=fabsl(Hm[i]); R d0=fabsl(F[1][i]-Hm[i]); if(d0>mdr) mdr=d0; Hc[i]=r; }
    printf("FD cross-check: max|H_FD(h=1e-3)-H_HD| = %.3Le ; max|H_FD(Richardson)-H_HD| = %.3Le ; max|H_HD| = %.3Le\n",mdr,md,mh);
    static R evf[N]; jacobi(N,Hc,evf,V);
    R me=0; for(int i=0;i<N;i++){ R d=fabsl(evf[i]-ev[i]); if(d>me) me=d; }
    printf("FD-Richardson spectrum vs HD spectrum: max eigenvalue diff = %.3Le\n",me);
    printf("FD-Richardson spectrum: ev[22]=%.6Le ev[23]=%.6Le ev[107]=%.6Le\n",evf[22],evf[23],evf[107]);
  }

  /* overlaps */
  { static R o[K*(K-1)/2*D*D]; int no=0;
    for(int a=0;a<K;a++)for(int b=a+1;b<K;b++)for(int pp=0;pp<D;pp++)for(int q=0;q<D;q++){
      R sr=0,si=0; for(int r=0;r<D;r++){ C x=M[a][r][pp],y=M[b][r][q]; sr+=x.re*y.re+x.im*y.im; si+=x.re*y.im-x.im*y.re; }
      o[no++]=sr*sr+si*si; }
    for(int i=0;i<no;i++)for(int j=i+1;j<no;j++) if(o[j]<o[i]){ R t=o[i];o[i]=o[j];o[j]=t; }
    printf("distinct |<a_i|b_j>|^2 values (cluster tol 1e-13), all 216 overlaps:\n");
    int i=0; while(i<no){ int j=i; R s=0; while(j<no && o[j]-o[i]<1e-13L){ s+=o[j]; j++; }
      printf("  %.18Lf  x%d   (spread %.2Le)\n", s/(j-i), j-i, o[j-1]-o[i]); i=j; }
    printf("per-pair distinct values:\n");
    for(int a=0;a<K;a++)for(int b=a+1;b<K;b++){ R v[36]; int nv=0;
      for(int pp=0;pp<D;pp++)for(int q=0;q<D;q++){ R sr=0,si=0; for(int r=0;r<D;r++){ C x=M[a][r][pp],y=M[b][r][q]; sr+=x.re*y.re+x.im*y.im; si+=x.re*y.im-x.im*y.re; } v[nv++]=sr*sr+si*si; }
      for(int x=0;x<36;x++)for(int y=x+1;y<36;y++) if(v[y]<v[x]){ R t=v[x];v[x]=v[y];v[y]=t; }
      printf("  pair (%d,%d):",a,b); int x=0; while(x<36){ int y=x; while(y<36&&v[y]-v[x]<1e-13L) y++; printf(" %.15Lf x%d;",v[x],y-x); x=y; } printf("\n"); }
  }
  if(argc>2){ FILE *g2=fopen(argv[2],"w"); fprintf(g2,"[");
    for(int k=0;k<K;k++){ fprintf(g2,"%s[",k?",":""); for(int c=0;c<D;c++){ fprintf(g2,"%s[",c?",":"");
      for(int r=0;r<D;r++) fprintf(g2,"%s[%.21Lg,%.21Lg]",r?",":"",M[k][r][c].re,M[k][r][c].im); fprintf(g2,"]"); } fprintf(g2,"]"); }
    fprintf(g2,"]\n"); fclose(g2); }
  return 0;
}
