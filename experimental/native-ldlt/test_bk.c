#define _POSIX_C_SOURCE 200809L
#include "bk_ldlt.h"
#include "factor_cache.h"
#include "factor_private.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
static unsigned pass,fail;static size_t direct,secondary,alternate,pairs,swaps,old_ok;static long double growth,mult;static double maxerr;static long double maxback;
static void check(int ok,const char *s){if(ok)pass++;else{fail++;fprintf(stderr,"FAIL %s\n",s);}}
static unsigned seed=12583;static double rnd(void){seed=1664525u*seed+1013904223u;return ((double)(seed>>8)/16777216.-.5)*2;}
static int solve_check(const p56_context *c,const double *a,unsigned n){
 double x[20],rhs[20],y[20],scratch[20];long double norm=0,bnorm=0,xnorm=0,res=0;double err=0;
 for(unsigned i=0;i<n;i++)x[i]=rnd();
 for(unsigned i=0;i<n;i++){long double sum=0,row=0;for(unsigned j=0;j<n;j++){sum+=(long double)a[i*n+j]*x[j];row+=fabsl(a[i*n+j]);}rhs[i]=(double)sum;norm=fmaxl(norm,row);bnorm=fmaxl(bnorm,fabsl(rhs[i]));}
 if(p56_apply(c,rhs,y,scratch,n))return 0;
 for(unsigned i=0;i<n;i++){long double sum=0;for(unsigned j=0;j<n;j++)sum+=(long double)a[i*n+j]*y[j];res=fmaxl(res,fabsl(sum-rhs[i]));xnorm=fmaxl(xnorm,fabsl(y[i]));err=fmax(err,fabs(x[i]-y[i]));}
 long double back=res/(norm*xnorm+bnorm);maxback=fmaxl(maxback,back);maxerr=fmax(maxerr,err);return back<1e-13L;
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
static void accumulate(p59_stats s){direct+=s.direct;secondary+=s.secondary;alternate+=s.alternate;pairs+=s.pairs;swaps+=s.swaps;growth=fmaxl(growth,s.growth);mult=fmaxl(mult,s.max_multiplier);}
static void random_suite(void){
 for(unsigned t=0;t<600;t++){
  unsigned n=3+t%14;double a[400]={0};double scale=t%3==0?1e-100:t%3==1?1:1e100;
  for(unsigned i=0;i<n;i++)for(unsigned j=0;j<=i;j++){double z=rnd()*scale;if(i==j&&t%4==0)z=0;a[i*n+j]=a[j*n+i]=z;}
  p58_symbolic sym={0,n,0,NULL};p56_context *c=NULL,*old=NULL;p59_stats st={0};
  int rc=p59_factorize(&c,n,a,1,&sym,n,0,&st);
  check(!rc&&st.direct+st.secondary+st.alternate+2*st.pairs==n&&reconstruction(c,a,n)&&solve_check(c,a,n),"random dense BK");if(!rc)accumulate(st);
  int oldrc=p58_factorize(&old,n,a,1,&sym,.1,n,NULL);if(!oldrc)old_ok++;
  if(t%20==0&&c){char path[]="/tmp/lis59-cache-XXXXXX";int fd=mkstemp(path);uint8_t key[32]={59};p56_context *loaded=NULL;if(fd>=0)close(fd);check(fd>=0&&!p57_save(c,path,key)&&!p57_load(&loaded,path,key,100000,n)&&reconstruction(loaded,a,n)&&solve_check(loaded,a,n),"cache roundtrip BK");if(fd>=0)unlink(path);p56_destroy(loaded);}
  p56_destroy(c);p56_destroy(old);
 }
}
static void branch(const double *a,unsigned n,int which){p58_symbolic sym={0,n,0,NULL};p56_context *c=NULL;p59_stats st={0};int rc=p59_factorize(&c,n,a,1,&sym,n,0,&st);int selected=which==0?st.direct>0:which==1?st.secondary>0:which==2?st.alternate>0:st.pairs>0;check(!rc&&selected&&reconstruction(c,a,n)&&solve_check(c,a,n),"forced BK branch");p56_destroy(c);}
int main(void){
 random_suite();
 double a[]={4,1,1,3},b[]={.2,1,0,1,0,4,0,4,1},c[]={0,1,1,4},d[]={0,0,2,0,0,0,0,3,2,0,0,0,0,3,0,0};
 branch(a,2,0);branch(b,3,1);branch(c,2,2);branch(d,4,3);
 p56_context *ctx=NULL;p58_symbolic full={0,2,0,NULL};p59_stats st={0};check(!p59_factorize(&ctx,2,a,1,&full,2,0,&st),"initial factor");p56_context *old=ctx;
 double g[]={1,1,1,-1};p59_stats before=st;check(p59_factorize(&ctx,2,g,1,&full,2,1.5,&st)==P59_GROWTH&&ctx==old&&st.growth==before.growth&&st.direct==before.direct&&solve_check(ctx,a,2),"growth limit rollback");
 check(!p59_factorize(&ctx,2,g,1,&full,2,2,&st)&&st.growth==2&&reconstruction(ctx,g,2),"growth exact boundary");old=ctx;
 double zero[4]={0};check(p59_factorize(&ctx,2,zero,1,&full,2,0,NULL)==P58_NO_PIVOT&&ctx==old,"singular rollback");
 uint32_t rows[]={1};p58_symbolic split[]={{0,1,1,rows},{1,1,0,NULL}};check(p59_factorize(&ctx,2,c,2,split,2,0,NULL)==P58_NO_PIVOT&&ctx==old,"restricted alternate scalar");
 double pair[]={0,1,1,0};check(p59_factorize(&ctx,2,pair,2,split,2,0,NULL)==P58_NO_PIVOT&&ctx==old,"restricted pair");
 check(!p59_factorize(&ctx,2,a,2,split,2,0,&st)&&reconstruction(ctx,a,2),"valid multi-panel");old=ctx;
 check(p59_factorize(&ctx,2,a,1,&full,1,0,NULL)==P55_INVALID&&ctx==old,"workspace limit");
 check(p59_factorize(&ctx,2,a,1,&full,2,.5,NULL)==P55_INVALID&&ctx==old,"invalid growth limit");
 double bad[]={1,2,1,3};check(p59_factorize(&ctx,2,bad,1,&full,2,0,NULL)==P55_INVALID&&ctx==old,"asymmetry");
 /* Strong alternative diagonal: all eligible stage58 pairs violate its multiplier limit. */
 double hard[]={0,1,0,1,1000,100,0,100,1};p58_symbolic three={0,3,0,NULL};
 check(p58_factorize(&ctx,3,hard,1,&three,.1,3,NULL)==P58_NO_PIVOT&&ctx==old,"stage58 rejects hard matrix");
 check(!p59_factorize(&ctx,3,hard,1,&three,3,0,&st)&&st.alternate>0&&reconstruction(ctx,hard,3)&&solve_check(ctx,hard,3),"BK rescues hard matrix");
 /* Previously eliminated column has distinct nonzero entries: swapping rows of L is essential. */
 double coupled[]={4,1,2,1,.25,1.5,2,1.5,5};uint32_t rr[]={1,2};p58_symbolic two[]={{0,1,2,rr},{1,2,0,NULL}};
 check(!p59_factorize(&ctx,3,coupled,2,two,3,0,&st)&&st.alternate>0&&st.swaps>0&&reconstruction(ctx,coupled,3)&&solve_check(ctx,coupled,3),"alternate swap updates prior L");
 old=ctx;uint32_t r1[]={1};p58_symbolic incompatible[]={{0,1,1,r1},{1,2,0,NULL}};
 double structure[]={4,1,0,1,.25,1,0,1,4};
 check(p59_factorize(&ctx,3,structure,2,incompatible,3,0,NULL)==P58_NO_PIVOT&&ctx==old,"incompatible prior panel membership");
 p56_destroy(ctx);
 printf("BK: %u/%u PASS; random cases=600, stage58 success=%zu/600; direct=%zu secondary=%zu alternate=%zu pairs=%zu swaps=%zu; max growth=%.6Le max multiplier=%.6Le; max solution error=%.3e max backward error=%.3Le\n",pass,pass+fail,old_ok,direct,secondary,alternate,pairs,swaps,growth,mult,maxerr,maxback);return fail?1:0;
}
