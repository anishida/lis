#define _POSIX_C_SOURCE 200809L
#include "pivot_ldlt.h"
#include "factor_cache.h"
#include "factor_private.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
static unsigned pass,fail;static size_t pairs,swaps;static double maxerr;static long double maxback;
static void check(int ok,const char *s){if(ok)pass++;else{fail++;fprintf(stderr,"FAIL %s\n",s);}}
static unsigned seed=12583;static double rnd(void){seed=1664525u*seed+1013904223u;return ((double)(seed>>8)/16777216.-.5)*2;}
static int solve_check(const p56_context *c,const double *a,unsigned n){
 double x[20],rhs[20],y[20],scratch[20];long double norm=0,bnorm=0,xnorm=0,res=0;double err=0;
 for(unsigned i=0;i<n;i++)x[i]=rnd();
 for(unsigned i=0;i<n;i++){long double sum=0,row=0;for(unsigned j=0;j<n;j++){sum+=(long double)a[i*n+j]*x[j];row+=fabsl(a[i*n+j]);}rhs[i]=(double)sum;norm=fmaxl(norm,row);bnorm=fmaxl(bnorm,fabsl(rhs[i]));}
 if(p56_apply(c,rhs,y,scratch,n))return 0;
 for(unsigned i=0;i<n;i++){long double sum=0;for(unsigned j=0;j<n;j++)sum+=(long double)a[i*n+j]*y[j];res=fmaxl(res,fabsl(sum-rhs[i]));xnorm=fmaxl(xnorm,fabsl(y[i]));err=fmax(err,fabs(x[i]-y[i]));}
 long double back=res/(norm*xnorm+bnorm);maxback=fmaxl(maxback,back);maxerr=fmax(maxerr,err);return err<1e-11&&back<1e-13L;
}
static int reconstruction(const p56_context *c,const double *a,unsigned n){
 double l[400]={0},d[400]={0};const p55_factor *f=c->factor;
 for(unsigned i=0;i<n;i++){l[i*n+i]=1;d[i*n+i]=f->diag[i];}
 for(size_t j=0;j<f->nb;j++){unsigned k=f->blocks[j].k;d[k*n+k+1]=d[(k+1)*n+k]=f->blocks[j].original_b;}
 for(size_t q=0;q<f->np;q++){const panel *p=f->p+q;size_t at=0;for(unsigned j=0;j<p->width;j++){unsigned k=p->first+j;for(unsigned r=j+1;r<p->width;r++)l[(p->first+r)*n+k]=p->values[at++];for(unsigned r=0;r<p->nrows;r++){unsigned i=p->short_index?p->first+p->width+((uint16_t*)p->rows)[r]:((uint32_t*)p->rows)[r];l[i*n+k]=p->values[at++];}}}
 long double err=0,norm=0;
 for(unsigned i=0;i<n;i++)for(unsigned j=0;j<n;j++){long double v=0;for(unsigned k=0;k<n;k++)for(unsigned r=0;r<n;r++)v+=(long double)l[i*n+k]*d[k*n+r]*l[j*n+r];unsigned u=c->perm?c->perm[i]:i,w=c->perm?c->perm[j]:j;err=fmaxl(err,fabsl(v-a[u*n+w]));norm=fmaxl(norm,fabsl(a[u*n+w]));}
 return err<=1e-13L*norm;
}
static void run_random(void){
 for(unsigned t=0;t<300;t++){
  /* Two panels: leading scalar pivots and trailing permuted 2x2 blocks.
     A = [D0, D0 V^T; V D0, S+V D0 V^T], known nonsingular. */
  enum{N=10};double a[N*N]={0},v[8][2];double scale=t%3==0?1e-100:t%3==1?1:1e100;unsigned perm[8];
  for(unsigned i=0;i<8;i++)perm[i]=i;
  for(unsigned i=7;i>0;i--){unsigned j=(unsigned)(fabs(rnd())*(i+1));if(j>i)j=i;unsigned z=perm[i];perm[i]=perm[j];perm[j]=z;}
  a[0]=4*scale;a[N+1]=-3*scale;
  for(unsigned i=0;i<8;i++)for(unsigned j=0;j<2;j++){v[i][j]=rnd()*.2;a[(i+2)*N+j]=a[j*N+i+2]=v[i][j]*a[j*N+j];}
  for(unsigned i=0;i<8;i++)for(unsigned j=0;j<=i;j++){double z=0;for(unsigned k=0;k<2;k++)z+=v[i][k]*a[k*N+k]*v[j][k];a[(i+2)*N+j+2]=a[(j+2)*N+i+2]=z;}
  for(unsigned i=0;i<8;i+=2){unsigned r=perm[i]+2,s=perm[i+1]+2;double z=(1+.3*rnd())*scale;a[r*N+s]+=z;a[s*N+r]+=z;}
  uint32_t rows[]={2,3,4,5,6,7,8,9};p58_symbolic sym[]={{0,2,8,rows},{2,8,0,NULL}};p56_context *c=NULL;p58_stats st={0};
  int rc=p58_factorize(&c,N,a,2,sym,.1,N,&st);int ok=!rc&&st.scalar==2&&st.pairs==4&&solve_check(c,a,N)&&reconstruction(c,a,N);check(ok,"random factorization/solve/reconstruct");if(!rc){pairs+=st.pairs;swaps+=st.swaps;}
  if(t%10==0&&c){char path[]="/tmp/lis58-cache-XXXXXX";int fd=mkstemp(path);uint8_t key[32]={58};p56_context *loaded=NULL;if(fd>=0)close(fd);int success=fd>=0&&!p57_save(c,path,key)&&!p57_load(&loaded,path,key,100000,N)&&reconstruction(loaded,a,N)&&solve_check(loaded,a,N);check(success,"cache roundtrip dynamic pivot");if(fd>=0)unlink(path);p56_destroy(loaded);}
  p56_destroy(c);
 }
}
static void run_coupled(void){
 /* Nonzero block multipliers induce genuine Schur updates in later panels.
    Build A=L D L^T, then reverse the order inside each panel. */
 for(unsigned t=0;t<120;t++){
  enum{N=8};double l[N*N]={0},d[N*N]={0},a[N*N]={0};
  double scale=t%3==0?1e-100:t%3==1?1:1e100;
  for(unsigned i=0;i<N;i++){l[i*N+i]=1;for(unsigned j=0;j<i;j++)if(!(i==j+1&&j%2==0))l[i*N+j]=rnd()*.08;}
  for(unsigned k=0;k<N;k+=2){d[k*N+k+1]=d[(k+1)*N+k]=(1+.2*rnd())*scale;}
  unsigned perm[]={0,2,1,3,4,6,5,7};
  for(unsigned i=0;i<N;i++)for(unsigned j=0;j<=i;j++){long double z=0;for(unsigned k=0;k<N;k++)for(unsigned r=0;r<N;r++)z+=(long double)l[perm[i]*N+k]*d[k*N+r]*l[perm[j]*N+r];a[i*N+j]=a[j*N+i]=(double)z;}
  uint32_t rows[]={4,5,6,7};p58_symbolic sym[]={{0,4,4,rows},{4,4,0,NULL}};p56_context *c=NULL;p58_stats st={0};
  int rc=p58_factorize(&c,N,a,2,sym,.1,N,&st);check(!rc&&st.pairs==4&&st.swaps==2&&reconstruction(c,a,N)&&solve_check(c,a,N),"coupled block Schur updates");
  if(!rc){pairs+=st.pairs;swaps+=st.swaps;}p56_destroy(c);
 }
}
int main(void){
 run_random();run_coupled();p56_context *c=NULL;p58_stats st={0};p58_symbolic full={0,4,0,NULL};
 double a[]={0,0,2,0, 0,0,0,3, 2,0,0,0, 0,3,0,0};
 int rc=p58_factorize(&c,4,a,1,&full,.1,4,&st);check(!rc&&st.swaps==1&&st.pairs==2&&solve_check(c,a,4)&&reconstruction(c,a,4),"nonadjacent pair");
 p56_context *old=c;p58_stats before=st;double singular[16]={0};check(p58_factorize(&c,4,singular,1,&full,.1,4,&st)==P58_NO_PIVOT&&c==old&&!memcmp(&st,&before,sizeof st)&&solve_check(c,a,4),"singular rollback");
 p58_symbolic split[]={{0,2,0,NULL},{2,2,0,NULL}};check(p58_factorize(&c,4,a,2,split,.1,4,NULL)==P58_NO_PIVOT&&c==old,"no cross-panel pivot");
 double fill[]={4,1,0,0,1,4,0,0,0,0,4,0,0,0,0,4};p58_symbolic tiny[]={{0,1,0,NULL},{1,3,0,NULL}};check(p58_factorize(&c,4,fill,2,tiny,.1,4,NULL)==P58_PATTERN&&c==old,"missing symbolic entry");
 double sp[]={4,0,1,0,0,0,0,2,1,0,0,0,0,2,0,0};uint32_t rr[]={2};p58_symbolic inc[]={{0,1,1,rr},{1,3,0,NULL}};
 /* Swap 2 with 3 would change membership in prior panel. */
 check(p58_factorize(&c,4,sp,2,inc,.1,4,NULL)==P58_NO_PIVOT&&c==old,"incompatible earlier row membership");
 double sparse[]={4,0,0,0,0,4,0,0,0,0,0,2,0,0,2,0};check(!p58_factorize(&c,4,sparse,2,split,.1,4,&st)&&st.pairs==1&&reconstruction(c,sparse,4)&&solve_check(c,sparse,4),"valid sparse pattern");old=c;
 double asym[16];memcpy(asym,sparse,sizeof asym);asym[1]=1;check(p58_factorize(&c,4,asym,1,&full,.1,4,NULL)==P55_INVALID&&c==old,"asymmetry");asym[1]=0;asym[0]=NAN;check(p58_factorize(&c,4,asym,1,&full,.1,4,NULL)==P55_INVALID&&c==old,"nonfinite input");
 check(p58_factorize(&c,4,sparse,1,&full,.1,3,NULL)==P55_INVALID&&c==old,"workspace limit");check(p58_factorize(&c,4,sparse,1,&full,0,4,NULL)==P55_INVALID&&c==old,"invalid threshold");
 double unstable[]={0,1,0,0,1,0,0,1000,0,0,1,0,0,1000,0,1};check(p58_factorize(&c,4,unstable,2,split,.1,4,NULL)==P58_NO_PIVOT&&c==old,"unstable pair rejected");
 p56_destroy(c);printf("pivot: %u/%u PASS; generated pairs=%zu swaps=%zu; max solution error=%.3e; max backward error=%.3Le\n",pass,pass+fail,pairs,swaps,maxerr,maxback);return fail?1:0;
}
