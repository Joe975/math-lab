/* Skeptic re-implementation: MUB loss minimisation.
 * Parametrisation: body-frame exponential map  U_b <- U_b * exp(A_b),
 * A_b anti-Hermitian (A = iH, H Hermitian). Analytic gradient in the Lie
 * algebra; L-BFGS on the stacked anti-Hermitian coordinates (identity
 * transport in body frame), Armijo backtracking. Basis 0 fixed = identity.
 * exp computed by scaling-and-squaring Taylor (order 20, ||A||<=1/4).
 * usage: mubexp d k nstarts seed maxit [json_out]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

typedef double complex cplx;
static int D, K;
#define MAXD 8
#define MAXK 8
typedef cplx Mat[MAXD][MAXD];

static unsigned long long rs;
static double urand(void){ rs ^= rs>>12; rs ^= rs<<25; rs ^= rs>>27;
  return ((rs*2685821657736338717ULL)>>11)*(1.0/9007199254740992.0); }
static double gauss(void){ double u=urand(),v=urand(); if(u<1e-300)u=1e-300;
  return sqrt(-2*log(u))*cos(2*M_PI*v); }

static void matmul(Mat C, Mat A, Mat B){ Mat T; for(int i=0;i<D;i++)for(int j=0;j<D;j++){cplx s=0;for(int l=0;l<D;l++)s+=A[i][l]*B[l][j];T[i][j]=s;} memcpy(C,T,sizeof(Mat)); }
/* C = A^dagger B */
static void matmulH(Mat C, Mat A, Mat B){ Mat T; for(int i=0;i<D;i++)for(int j=0;j<D;j++){cplx s=0;for(int l=0;l<D;l++)s+=conj(A[l][i])*B[l][j];T[i][j]=s;} memcpy(C,T,sizeof(Mat)); }

/* Gram-Schmidt on columns: used for Haar init and as roundoff cleanup */
static void gs(Mat U){ for(int j=0;j<D;j++){ for(int p=0;p<2;p++) for(int q=0;q<j;q++){cplx s=0;for(int i=0;i<D;i++)s+=conj(U[i][q])*U[i][j];for(int i=0;i<D;i++)U[i][j]-=s*U[i][q];}
  double n=0;for(int i=0;i<D;i++)n+=creal(U[i][j]*conj(U[i][j]));n=sqrt(n);for(int i=0;i<D;i++)U[i][j]/=n;} }
static void haar(Mat U){ for(int i=0;i<D;i++)for(int j=0;j<D;j++)U[i][j]=gauss()+I*gauss(); gs(U); }

static void expm(Mat E, Mat A){ double nrm=0; for(int i=0;i<D;i++)for(int j=0;j<D;j++)nrm+=creal(A[i][j]*conj(A[i][j])); nrm=sqrt(nrm);
  int s=0; while(nrm>0.25){nrm/=2;s++;} double sc=ldexp(1.0,-s);
  Mat B,T; for(int i=0;i<D;i++)for(int j=0;j<D;j++)B[i][j]=A[i][j]*sc;
  for(int i=0;i<D;i++)for(int j=0;j<D;j++)E[i][j]=(i==j);
  memcpy(T,E,sizeof(Mat));
  for(int n=1;n<=20;n++){ matmul(T,T,B); for(int i=0;i<D;i++)for(int j=0;j<D;j++){T[i][j]/=n;E[i][j]+=T[i][j];} }
  for(int r=0;r<s;r++) matmul(E,E,E); }

static Mat U[MAXK];
/* loss and gradient (anti-Hermitian G[b], b>=1) */
static double lossgrad(Mat *Us, Mat *G){
  double L=0, invd=1.0/D;
  if(G) for(int b=0;b<K;b++) memset(G[b],0,sizeof(Mat));
  for(int a=0;a<K;a++)for(int b=a+1;b<K;b++){
    Mat M,W,Z; matmulH(M,Us[a],Us[b]);
    for(int i=0;i<D;i++)for(int j=0;j<D;j++){ double r=creal(M[i][j]*conj(M[i][j]))-invd; L+=r*r; W[i][j]=4*r*M[i][j]; }
    if(!G) continue;
    matmulH(Z,M,W); /* for b: skew(M^dag W) */
    for(int i=0;i<D;i++)for(int j=0;j<D;j++) G[b][i][j]+=0.5*(Z[i][j]-conj(Z[j][i]));
    /* for a: -skew(W M^dag) */
    for(int i=0;i<D;i++)for(int j=0;j<D;j++){cplx s=0;for(int l=0;l<D;l++)s+=W[i][l]*conj(M[j][l]);Z[i][j]=s;}
    for(int i=0;i<D;i++)for(int j=0;j<D;j++) G[a][i][j]-=0.5*(Z[i][j]-conj(Z[j][i]));
  }
  if(G) memset(G[0],0,sizeof(Mat));
  return L;
}
static double ip(Mat *X, Mat *Y){ double s=0; for(int b=1;b<K;b++)for(int i=0;i<D;i++)for(int j=0;j<D;j++) s+=creal(conj(X[b][i][j])*Y[b][i][j]); return s; }

static void step(Mat *Out, Mat *In, Mat *P, double t){ /* Out_b = In_b exp(t P_b) */
  for(int b=0;b<K;b++){ if(b==0){memcpy(Out[0],In[0],sizeof(Mat));continue;} Mat A,E; for(int i=0;i<D;i++)for(int j=0;j<D;j++)A[i][j]=t*P[b][i][j]; expm(E,A); matmul(Out[b],In[b],E);} }

#define MH 12
static double lastgn;
static double optimize(int maxit, int *iters){
  static Mat G[MAXK],Gn[MAXK],P[MAXK],Ut[MAXK],S[MH][MAXK],Y[MH][MAXK];
  double rho[MH],al[MH]; int h=0,hs=0;
  double L=lossgrad(U,G); int it, stall=0;
  for(it=0;it<maxit;it++){
    double gn=sqrt(ip(G,G)); if(gn<1e-12 || L<1e-28) break;
    /* two-loop recursion; P = -H g */
    memcpy(P,G,sizeof(Mat)*K);
    for(int m=0;m<hs;m++){int q=(h-1-m+MH)%MH; al[q]=rho[q]*ip(S[q],P); for(int b=1;b<K;b++)for(int i=0;i<D;i++)for(int j=0;j<D;j++)P[b][i][j]-=al[q]*Y[q][b][i][j];}
    if(hs>0){int q=(h-1+MH)%MH; double gam=ip(S[q],Y[q])/ip(Y[q],Y[q]); for(int b=1;b<K;b++)for(int i=0;i<D;i++)for(int j=0;j<D;j++)P[b][i][j]*=gam;}
    for(int m=hs-1;m>=0;m--){int q=(h-1-m+MH)%MH; double be=rho[q]*ip(Y[q],P); for(int b=1;b<K;b++)for(int i=0;i<D;i++)for(int j=0;j<D;j++)P[b][i][j]+=(al[q]-be)*S[q][b][i][j];}
    for(int b=1;b<K;b++)for(int i=0;i<D;i++)for(int j=0;j<D;j++)P[b][i][j]=-P[b][i][j];
    double dd=ip(P,G);
    if(!(dd<0)){ for(int b=1;b<K;b++)for(int i=0;i<D;i++)for(int j=0;j<D;j++)P[b][i][j]=-G[b][i][j]; dd=-gn*gn; hs=0; }
    double t = (hs==0)? fmin(1.0,0.1/gn) : 1.0, Ln=0; int ok=0;
    for(int ls=0;ls<60;ls++){ step(Ut,U,P,t); Ln=lossgrad(Ut,NULL); if(Ln<=L+1e-4*t*dd){ok=1;break;} t*=0.5; }
    if(!ok){ if(hs==0) break; hs=0; continue; }
    lossgrad(Ut,Gn);
    /* s = t P, y = Gn - G (body frame, identity transport) */
    for(int b=0;b<K;b++)for(int i=0;i<D;i++)for(int j=0;j<D;j++){S[h][b][i][j]=t*P[b][i][j]; Y[h][b][i][j]=Gn[b][i][j]-G[b][i][j];}
    double sy=ip(S[h],Y[h]); if(sy>1e-300){ rho[h]=1.0/sy; h=(h+1)%MH; if(hs<MH)hs++; }
    if(L-Ln < 1e-15*L) stall++; else stall=0;
    memcpy(U,Ut,sizeof(Mat)*K); memcpy(G,Gn,sizeof(Mat)*K); L=Ln;
    if(stall>=30) { it++; break; }
    if(it%500==499){ for(int b=1;b<K;b++) gs(U[b]); L=lossgrad(U,G); }
  }
  *iters=it; lastgn=sqrt(ip(G,G)); return L;
}
static double orthdev(Mat X){ Mat M; matmulH(M,X,X); double m=0; for(int i=0;i<D;i++)for(int j=0;j<D;j++){double e=cabs(M[i][j]-(i==j)); if(e>m)m=e;} return m; }

static int checkgrad(void){ /* finite-difference check of the analytic gradient */
  static Mat G[MAXK],P[MAXK],Up[MAXK],Um[MAXK];
  for(int b=1;b<K;b++) haar(U[b]); for(int i=0;i<D;i++)for(int j=0;j<D;j++)U[0][i][j]=(i==j);
  lossgrad(U,G);
  for(int b=0;b<K;b++)for(int i=0;i<D;i++)for(int j=0;j<=i;j++){cplx z=gauss()+I*gauss(); if(i==j) z=I*cimag(z); P[b][i][j]=z; P[b][j][i]=-conj(z);}
  double h=1e-5; step(Up,U,P,h); step(Um,U,P,-h);
  double fd=(lossgrad(Up,NULL)-lossgrad(Um,NULL))/(2*h), an=ip(G,P);
  fprintf(stderr,"gradcheck: analytic %.12g  fd %.12g  rel %.2e\n",an,fd,fabs(an-fd)/fabs(an)); return 0; }

static int cmpd(const void*a,const void*b){double x=*(double*)a,y=*(double*)b;return (x>y)-(x<y);}
int main(int argc,char**argv){
  if(argc<6){fprintf(stderr,"usage: d k nstarts seed maxit [json]\n");return 1;}
  D=atoi(argv[1]);K=atoi(argv[2]);int N=atoi(argv[3]);rs=strtoull(argv[4],0,10)*0x9E3779B97F4A7C15ULL+1;int maxit=atoi(argv[5]);
  for(int w=0;w<10;w++)urand();
  checkgrad();
  double *fin=malloc(sizeof(double)*N),best=1e300; static Mat Best[MAXK]; long totit=0; int capped=0;
  for(int s=0;s<N;s++){
    for(int i=0;i<D;i++)for(int j=0;j<D;j++)U[0][i][j]=(i==j);
    for(int b=1;b<K;b++) haar(U[b]);
    int its; double L=optimize(maxit,&its); totit+=its; if(its>=maxit)capped++;
    for(int b=1;b<K;b++) gs(U[b]); L=lossgrad(U,NULL); /* final roundoff cleanup, recompute */
    fin[s]=L; printf("start %d L %.15e iters %d gradnorm %.2e\n",s,L,its,lastgn); fflush(stdout);
    if(L<best){best=L;memcpy(Best,U,sizeof(Mat)*K);}
  }
  qsort(fin,N,sizeof(double),cmpd);
  double C2=K*(K-1)/2.0;
  printf("SUMMARY d=%d k=%d starts=%d best=%.15e  ASD_best=%.15f\n",D,K,N,best,1-best/(C2*(D-1)));
  printf("quantiles min %.10e q10 %.10e q25 %.10e med %.10e q75 %.10e q90 %.10e max %.10e\n",
    fin[0],fin[N/10],fin[N/4],fin[N/2],fin[3*N/4],fin[9*N/10],fin[N-1]);
  int nb=0; for(int s=0;s<N;s++) if(fin[s]<=best+1e-8) nb++;
  printf("starts within 1e-8 of best: %d ; within 1e-6: ",nb); nb=0; for(int s=0;s<N;s++) if(fin[s]<=best+1e-6) nb++; printf("%d ; hit maxit: %d ; mean iters %.1f\n",nb,capped,(double)totit/N);
  double od=0; for(int b=0;b<K;b++){double e=orthdev(Best[b]); if(e>od)od=e;} printf("best config max |U^dag U - I| = %.3e\n",od);
  if(argc>6){ FILE*f=fopen(argv[6],"w"); fprintf(f,"[");
    for(int b=0;b<K;b++){fprintf(f,"%s[",b?",":""); for(int j=0;j<D;j++){fprintf(f,"%s[",j?",":""); for(int i=0;i<D;i++) fprintf(f,"%s[%.17g,%.17g]",i?",":"",creal(Best[b][i][j]),cimag(Best[b][i][j])); fprintf(f,"]");} fprintf(f,"]");}
    fprintf(f,"]\n"); fclose(f); }
  return 0; }
