#define _POSIX_C_SOURCE 200809L
#include "delayed_ldlt.h"
#include "factor_cache.h"
#include "factor_private.h"
#include <math.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
static unsigned pass,fail;static size_t delays,pairs,swaps,old_ok,rescued,peak_front,peak_pending;static long double growth;static double maxerr;static long double maxback;
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
static p60_limits limits(unsigned n){return (p60_limits){n,n,n*n,n,n*n,0};}
static void accumulate(p60_stats s){delays+=s.delays;pairs+=s.pairs;swaps+=s.swaps;growth=fmaxl(growth,s.growth);if(s.peak_front>peak_front)peak_front=s.peak_front;if(s.peak_pending>peak_pending)peak_pending=s.peak_pending;}
static int valid(const p56_context *c,const double *a,unsigned n,p60_stats st){return st.direct+st.secondary+st.alternate+2*st.pairs==n&&reconstruction(c,a,n)&&solve_check(c,a,n);}
static void roundtrip(const p56_context *c,const double *a,unsigned n){char path[]="/tmp/lis60-cache-XXXXXX";int fd=mkstemp(path);uint8_t key[32]={60};p56_context *loaded=NULL;if(fd>=0)close(fd);check(fd>=0&&!p57_save(c,path,key)&&!p57_load(&loaded,path,key,100000,n)&&reconstruction(loaded,a,n)&&solve_check(loaded,a,n),"cache roundtrip delayed");if(fd>=0)unlink(path);p56_destroy(loaded);}
static void random_suite(void){
 for(unsigned t=0;t<400;t++){
  unsigned n=3+t%14;double a[400]={0};double scale=t%3==0?1e-100:t%3==1?1:1e100;
  for(unsigned i=0;i<n;i++)for(unsigned j=0;j<=i;j++){double z=rnd()*scale;if(i==j&&t%3==0)z=0;a[i*n+j]=a[j*n+i]=z;}
  uint32_t rows[20];p58_symbolic s[20];size_t np=0;for(unsigned i=0;i<n;i++)rows[i]=i;
  for(unsigned i=0;i<n;){unsigned w=1+(t+i)%3;if(w>n-i)w=n-i;s[np++]=(p58_symbolic){i,w,n-i-w,rows+i+w};i+=w;}
  p56_context *c=NULL,*old=NULL;p60_stats st={0};p60_limits lim=limits(n);int rc=p60_factorize(&c,n,a,np,s,&lim,&st);check(!rc&&valid(c,a,n,st)&&st.factor_slots<=st.baseline_slots+lim.max_extra_slots,"random split panel");
  int oldrc=p59_factorize(&old,n,a,np,s,n,0,NULL);if(!oldrc)old_ok++;if(!rc){accumulate(st);if(oldrc==P58_NO_PIVOT)rescued++;if(t%20==0)roundtrip(c,a,n);}p56_destroy(c);p56_destroy(old);
 }
 for(unsigned t=0;t<120;t++){
  unsigned n=4+2*(t%7),m=n/2;double a[400]={0};uint32_t rows[20];p58_symbolic s[20];double scale=t%3==0?1e-100:t%3==1?1:1e100;
  for(unsigned i=0;i<m;i++){unsigned j=m+(i+t)%m;double b=(1+.25*rnd())*scale;a[i*n+j]=a[j*n+i]=b;rows[i]=j;s[i]=(p58_symbolic){i,1,1,rows+i};}for(unsigned i=m;i<n;i++)s[i]=(p58_symbolic){i,1,0,NULL};
  p56_context *c=NULL;p60_stats st={0};p60_limits lim=limits(n);int rc=p60_factorize(&c,n,a,n,s,&lim,&st);check(!rc&&st.delays>0&&st.pairs==m&&valid(c,a,n,st),"sparse long-distance pairs");if(!rc){accumulate(st);if(t%12==0)roundtrip(c,a,n);}p56_destroy(c);
 }
}
int main(void){
 random_suite();p56_context *c=NULL;uint32_t rr[]={1};p58_symbolic s[]={{0,1,1,rr},{1,1,0,NULL}};p60_limits lim=limits(2);p60_stats st={0};
 double pair[]={0,2,2,0},alt[]={0,1,1,4};
 check(p59_factorize(&c,2,pair,2,s,2,0,NULL)==P58_NO_PIVOT&&c==NULL,"stage59 rejects split pair");
 check(!p60_factorize(&c,2,pair,2,s,&lim,&st)&&st.delays==1&&st.pairs==1&&st.panels==1&&valid(c,pair,2,st),"stage60 rescues split pair");
 check(!p60_factorize(&c,2,alt,2,s,&lim,&st)&&st.alternate==1&&st.delays==1&&valid(c,alt,2,st),"delayed alternate");p56_context *old=c;p60_stats prior=st;
 lim.max_delays=0;check(p60_factorize(&c,2,pair,2,s,&lim,&st)==P60_LIMIT&&c==old&&st.delays==prior.delays,"delay cap rollback");lim=limits(2);
 lim.max_pending=0;check(p60_factorize(&c,2,pair,2,s,&lim,NULL)==P60_LIMIT&&c==old,"pending cap rollback");lim=limits(2);
 lim.max_front=1;check(p60_factorize(&c,2,pair,2,s,&lim,NULL)==P60_LIMIT&&c==old,"front cap rollback");lim=limits(2);
 lim.max_n=1;check(p60_factorize(&c,2,pair,2,s,&lim,NULL)==P55_INVALID&&c==old,"workspace cap rollback");lim=limits(2);
 double singular[]={0,0,0,0};check(p60_factorize(&c,2,singular,2,s,&lim,NULL)==P58_NO_PIVOT&&c==old,"singular rollback");
 p58_symbolic missing[]={{0,1,0,NULL},{1,1,0,NULL}};check(p60_factorize(&c,2,pair,2,missing,&lim,NULL)==P58_PATTERN&&c==old,"input pattern mismatch");
 double bad[]={1,2,1,3};check(p60_factorize(&c,2,bad,2,s,&lim,NULL)==P55_INVALID&&c==old,"asymmetry");
 double g[]={1,1,1,-1};lim.max_growth=1.5;check(p60_factorize(&c,2,g,2,s,&lim,NULL)==P59_GROWTH&&c==old,"growth cap rollback");
 uint32_t r2[]={2};p58_symbolic sparse[]={{0,1,1,r2},{1,1,1,r2},{2,1,0,NULL}};double fill[]={.1,0,1,0,0,2,1,2,0};lim=limits(3);lim.max_extra_slots=0;
 check(p60_factorize(&c,3,fill,3,sparse,&lim,NULL)==P60_LIMIT&&c==old&&solve_check(c,alt,2),"extra slots cap rollback");
 lim.max_extra_slots=1;check(!p60_factorize(&c,3,fill,3,sparse,&lim,&st)&&st.baseline_slots==2&&st.factor_slots==3&&valid(c,fill,3,st),"extra slots exact boundary");if(c)roundtrip(c,fill,3);
 p56_destroy(c);
 printf("DELAY: %u/%u PASS; random split=400, sparse pairs=120; stage59 split success=%zu/400; rescued=%zu/400; delays=%zu pairs=%zu swaps=%zu peak_front=%zu peak_pending=%zu; max growth=%.6Le; max solution error=%.3e max backward error=%.3Le\n",pass,pass+fail,old_ok,rescued,delays,pairs,swaps,peak_front,peak_pending,growth,maxerr,maxback);return fail?1:0;
}
