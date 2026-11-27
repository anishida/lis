#include "supernodal_ldlt.h"
#include "factor_cache.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#x);exit(1);}}while(0)
static unsigned tests;static uint32_t seed=7361;static uint32_t rnd(void){seed^=seed<<13;seed^=seed>>17;seed^=seed<<5;return seed;}
static p62_limits sl={128u*1024*1024,100000000};
static p61_limits nl={256,128u*1024*1024,10000000,1000000,0};
static long double worst;static double worst_x;
static void verify_solve(p56_context *ctx,uint32_t n,const size_t *col,const uint32_t *row,const double *a){
 double *b=calloc(n,sizeof *b),*x=malloc((size_t)n*sizeof *x),*y=malloc((size_t)n*sizeof *y),*w=malloc((size_t)n*sizeof *w);long double *r=calloc(n,sizeof *r),*d=calloc(n,sizeof *d);CHECK(b&&x&&y&&w&&r&&d);
 for(uint32_t i=0;i<n;i++)x[i]=((int)(i%17)-8)/9.;
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++){uint32_t i=row[z];b[i]+=a[z]*x[j];if(i!=j)b[j]+=a[z]*x[i];}
 CHECK(!p56_apply(ctx,b,y,w,n));
 for(uint32_t i=0;i<n;i++){r[i]=-b[i];d[i]=fabsl(b[i]);double e=fabs(y[i]-x[i]);CHECK(e<1e-9);if(e>worst_x)worst_x=e;}
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++){uint32_t i=row[z];long double t=(long double)a[z]*y[j];r[i]+=t;d[i]+=fabsl(t);if(i!=j){t=(long double)a[z]*y[i];r[j]+=t;d[j]+=fabsl(t);}}
 for(uint32_t i=0;i<n;i++){long double e=d[i]?fabsl(r[i])/d[i]:0;CHECK(e<1e-12L);if(e>worst)worst=e;}
 free(b);free(x);free(y);free(w);free(r);free(d);
}
static p64_limits lim={128u*1024*1024,100000000,500000000,8,64,1000000};
static size_t fronts1,fronts8,delays1,delays8,pairs8;
static void random_case(unsigned t){
 uint32_t n=2+rnd()%62,ri[4096];size_t cp[65],z=0;double a[4096],scale=t%3==0?1e-100:t%3==1?1e100:1.;
 for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=(j%2?1.:-1.)*scale*(2*n+1);for(uint32_t i=j+1;i<n;i++)if(rnd()%100<t%100){ri[z]=i;a[z++]=scale*((int)(rnd()%15)-7)/8.;}}cp[n]=z;
 p56_context *ctx=NULL;p61_stats st={0};p62_symbolic *sym=NULL;CHECK(!p62_analyze(&sym,n,z,cp,ri,&sl));
 uint32_t w=2+t%15;size_t expect=0;for(size_t g=0;g<sym->nsuper;g++)expect+=1+(sym->super[g+1]-sym->super[g]-1)/w;
 CHECK(!p65_factorize(&ctx,n,z,cp,ri,a,&sl,&nl,w,&st));CHECK(st.fronts==expect&&st.scalars+2*st.pairs==n);verify_solve(ctx,n,cp,ri,a);fronts1+=n;fronts8+=st.fronts;p56_destroy(ctx);p62_free(sym);tests++;
 if(t<150){p64_stats os={0};ctx=NULL;CHECK(!p65_ordered_factorize(&ctx,n,z,cp,ri,a,&lim,&sl,&nl,(p64_method)(t%4),w,&os));verify_solve(ctx,n,cp,ri,a);p56_destroy(ctx);tests++;}
}
static void paired(unsigned t){uint32_t n=2*(1+rnd()%30),ri[4096];size_t cp[65],z=0;double a[4096],scale=t%3==0?1e-100:t%3==1?1e100:1.;
 /* Weak nonzero inter-pair couplings on complete structural graph, strong pairs. */
 for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=0;for(uint32_t i=j+1;i<n;i++){ri[z]=i;a[z++]=scale*(i==j+1&&j%2==0?4.:((int)(rnd()%9)-4)*.001);}}cp[n]=z;
 p56_context *ctx=NULL;p61_stats st={0};CHECK(!p65_factorize(&ctx,n,z,cp,ri,a,&sl,&nl,1,&st));verify_solve(ctx,n,cp,ri,a);delays1+=st.delays;
 CHECK(!p65_factorize(&ctx,n,z,cp,ri,a,&sl,&nl,8,&st));verify_solve(ctx,n,cp,ri,a);delays8+=st.delays;pairs8+=st.pairs;CHECK(st.pairs>0);p56_destroy(ctx);tests++;
}
static void specials(void){size_t cp[]={0,2,3};uint32_t ri[]={0,1,1};double a[]={0,2,0};p56_context *ctx=NULL;p61_stats st={0};p61_limits low=nl;low.max_delays=0;
 CHECK(!p65_factorize(&ctx,2,3,cp,ri,a,&sl,&low,2,&st)&&st.fronts==1&&st.pairs==1&&st.delays==0);verify_solve(ctx,2,cp,ri,a);tests++;
 p56_context *old=ctx;p61_stats before=st;CHECK(p65_factorize(&ctx,2,3,cp,ri,a,&sl,&low,1,&st)==P60_LIMIT&&ctx==old&&!memcmp(&st,&before,sizeof st));tests++;
 low=nl;low.max_front=1;CHECK(p65_factorize(&ctx,2,3,cp,ri,a,&sl,&low,2,&st)==P60_LIMIT&&ctx==old);tests++;
 CHECK(p65_factorize(&ctx,2,3,cp,ri,a,&sl,&nl,0,&st)==P55_INVALID&&ctx==old);tests++;
 for(size_t b=1;b<before.peak_work_bytes;b+=17){low=nl;low.max_work_bytes=b;CHECK(p65_factorize(&ctx,2,3,cp,ri,a,&sl,&low,2,&st)==P60_LIMIT&&ctx==old&&!memcmp(&st,&before,sizeof st));tests++;}
 unsigned char key[32]={0};p56_context *restored=NULL;CHECK(!p57_save(ctx,"stage65-cache.bin",key));CHECK(!p57_load(&restored,"stage65-cache.bin",key,100000,2));verify_solve(restored,2,cp,ri,a);remove("stage65-cache.bin");p56_destroy(restored);p56_destroy(ctx);tests++;
 /* Two children point to different columns of one supernode, including delayed child. */
 size_t c[]={0,2,4,6,7};uint32_t r[]={0,2,1,3,2,3,3};double v[]={0,2,4,1,0,1,5};ctx=NULL;
 CHECK(!p65_factorize(&ctx,4,7,c,r,v,&sl,&nl,8,&st));CHECK(st.fronts==3&&st.delays>0);verify_solve(ctx,4,c,r,v);p56_destroy(ctx);tests++;
}
static void bench(uint32_t n,uint32_t side,uint32_t block){size_t *cp=malloc(((size_t)n+1)*sizeof *cp);uint32_t *ri=malloc((size_t)n*(block?block:3)*sizeof *ri);double *a=malloc((size_t)n*(block?block:3)*sizeof *a);CHECK(cp&&ri&&a);size_t z=0;
 for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=block?2*block:5;if(block){uint32_t end=(j/block+1)*block;if(end>n)end=n;for(uint32_t i=j+1;i<end;i++){ri[z]=i;a[z++]=-.25;}}else{if(j+1<n&&(!side||j%side+1<side)){ri[z]=j+1;a[z++]=-1;}if(side&&j+side<n){ri[z]=j+side;a[z++]=-1;}}}cp[n]=z;
 for(uint32_t w=1;w<=16;w*=4){double best=1e9; p64_stats st={0};for(int rep=0;rep<3;rep++){p56_context *ctx=NULL;clock_t start=clock();CHECK(!p65_ordered_factorize(&ctx,n,z,cp,ri,a,&lim,&sl,&nl,side?P64_ND:P64_NATURAL,w,&st));double tm=(double)(clock()-start)/CLOCKS_PER_SEC;if(tm<best)best=tm;verify_solve(ctx,n,cp,ri,a);p56_destroy(ctx);}printf("bench n=%u side=%u block=%u width=%u fronts=%zu max_front=%zu slots=%zu workspace=%zu cpu_min3=%.6f s\n",n,side,block,w,st.numeric.fronts,st.numeric.peak_front,st.numeric.factor_slots,st.numeric.peak_work_bytes,best);tests++;}
 free(cp);free(ri);free(a);
}
int main(void){for(unsigned i=0;i<400;i++)random_case(i);for(unsigned i=0;i<200;i++)paired(i);specials();bench(1600,40,0);bench(8192,0,16);bench(152592,0,0);printf("supernodal %u/%u PASS max_x=%.3e max_componentwise_backward=%.3Le random_fronts_legacy=%zu grouped=%zu paired_delays_legacy=%zu grouped=%zu grouped_pairs=%zu\n",tests,tests,worst_x,worst,fronts1,fronts8,delays1,delays8,pairs8);return 0;}
