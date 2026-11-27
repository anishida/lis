#include "minimum_degree.h"
#include "factor_private.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#define NONE UINT32_MAX
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
static p63_limits ol={128u*1024*1024,100000000,500000000};
static void oracle(uint32_t n,const size_t *col,const uint32_t *row,const p63_order *o){
 unsigned char g[64][64]={{0}},live[64];CHECK(n<=64);memset(live,1,n);size_t slots=0,front=0,pairs=0;
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++)if(row[z]!=j)g[row[z]][j]=g[j][row[z]]=1;
 for(uint32_t k=0;k<n;k++){uint32_t best=NONE;size_t d=SIZE_MAX;
  for(uint32_t j=0;j<n;j++)if(live[j]){size_t t=0;for(uint32_t i=0;i<n;i++)if(live[i]&&g[i][j])t++;if(t<d){best=j;d=t;}}
  CHECK(o->perm[k]==best&&o->inv[best]==k);slots+=d;pairs+=d*(d?d-1:0)/2;if(d+1>front)front=d+1;live[best]=0;
  for(uint32_t i=0;i<n;i++)if(live[i]&&g[i][best])for(uint32_t j=i+1;j<n;j++)if(live[j]&&g[j][best])g[i][j]=g[j][i]=1;
 }
 CHECK(o->slots==slots&&o->max_front==front&&o->pair_checks==pairs);
}
static void random_case(unsigned t){
 uint32_t n=1+rnd()%64,row[64*64];size_t col[65],z=0;double a[64*64],scale=t%3==0?1e-100:t%3==1?1e100:1.;
 for(uint32_t j=0;j<n;j++){col[j]=z;row[z]=j;a[z++]=(j%2?1.:-1.)*scale*(2*n+1);for(uint32_t i=j+1;i<n;i++)if(rnd()%100<t%100){row[z]=i;a[z++]=scale*((int)(rnd()%15)-7)/8.;}}
 col[n]=z;p63_order *o=NULL;CHECK(!p63_analyze(&o,n,z,col,row,&ol));oracle(n,col,row,o);
 if(t<300){p56_context *ctx=NULL;p63_stats st={0};CHECK(!p63_factorize(&ctx,n,z,col,row,a,&ol,&sl,&nl,&st));CHECK(st.numeric.factor_slots<=o->slots&&st.numeric.peak_front<=o->max_front&&st.summary.perm==NULL&&st.summary.inv==NULL);verify_solve(ctx,n,col,row,a);p56_destroy(ctx);}
 p63_free(o);tests++;
}
static void paired(unsigned t){uint32_t n=2*(1+rnd()%30),row[128];size_t col[65],z=0;double a[128],scale=t%2?1e100:1e-100;
 for(uint32_t j=0;j<n;j++){col[j]=z;if(t%2){row[z]=j;a[z++]=0;}if(j<n/2){row[z]=j+n/2;a[z++]=scale*(2+rnd()%7);}}col[n]=z;
 p63_order *o=NULL;CHECK(!p63_analyze(&o,n,z,col,row,&ol));oracle(n,col,row,o);p63_free(o);p56_context *ctx=NULL;p63_stats st={0};CHECK(!p63_factorize(&ctx,n,z,col,row,a,&ol,&sl,&nl,&st));CHECK(st.numeric.pairs==n/2&&st.numeric.delays==n/2);verify_solve(ctx,n,col,row,a);p56_destroy(ctx);tests++;
}
static void bad(void){size_t col[]={0,3,5,6};uint32_t row[]={0,1,2,1,2,2};double a[]={5,1,1,5,1,5};p63_order *o=NULL;CHECK(!p63_analyze(&o,3,6,col,row,&ol));p63_order *old=o;p63_limits l=ol;
 #define BAD(expr,code) do{CHECK((expr)==(code));CHECK(o==old);tests++;}while(0)
 BAD(p63_analyze(&o,3,5,col,row,&ol),P55_INVALID);BAD(p63_analyze(&o,3,6,col,NULL,&ol),P55_INVALID);row[1]=0;BAD(p63_analyze(&o,3,6,col,row,&ol),P55_INVALID);row[1]=1;
 l.max_slots=2;BAD(p63_analyze(&o,3,6,col,row,&l),P60_LIMIT);l=ol;l.max_pair_checks=0;BAD(p63_analyze(&o,3,6,col,row,&l),P60_LIMIT);
 l=ol;size_t peak=o->peak_work_bytes;for(size_t i=0;i<peak;i+=7){l.max_work_bytes=i;BAD(p63_analyze(&o,3,6,col,row,&l),P60_LIMIT);}l.max_work_bytes=peak;CHECK(!p63_analyze(&o,3,6,col,row,&l));p63_free(o);tests++;
 p56_context *ctx=NULL;p63_stats st={0};CHECK(!p63_factorize(&ctx,3,6,col,row,a,&ol,&sl,&nl,&st));p56_context *before=ctx;p63_stats saved=st;l=ol;l.max_slots=0;CHECK(p63_factorize(&ctx,3,6,col,row,a,&l,&sl,&nl,&st)==P60_LIMIT);CHECK(ctx==before&&!memcmp(&st,&saved,sizeof st));tests++;
 p61_limits num=nl;num.max_front=1;CHECK(p63_factorize(&ctx,3,6,col,row,a,&ol,&sl,&num,&st)==P60_LIMIT);CHECK(ctx==before&&!memcmp(&st,&saved,sizeof st));verify_solve(ctx,3,col,row,a);p56_destroy(ctx);tests++;
 size_t empty[]={0,0,0};o=NULL;CHECK(!p63_analyze(&o,2,0,empty,NULL,&ol));CHECK(o->slots==0&&o->max_front==1&&o->perm[0]==0);p63_free(o);tests++;
}
static void large(uint32_t n,uint32_t side,int numeric){
 size_t *col=malloc(((size_t)n+1)*sizeof *col);uint32_t *row=malloc(3*(size_t)n*sizeof *row);double *a=malloc(3*(size_t)n*sizeof *a);CHECK(col&&row&&a);size_t z=0;
 for(uint32_t j=0;j<n;j++){col[j]=z;row[z]=j;a[z++]=5;if(j+1<n&&(!side||j%side+1<side)){row[z]=j+1;a[z++]=-1;}if(side&&j+side<n){row[z]=j+side;a[z++]=-1;}}col[n]=z;
 p62_symbolic *s=NULL;CHECK(!p62_analyze(&s,n,z,col,row,&sl));p63_order *o=NULL;clock_t start=clock();CHECK(!p63_analyze(&o,n,z,col,row,&ol));double t=(double)(clock()-start)/CLOCKS_PER_SEC;
 printf("order n=%u side=%u natural_slots=%zu md_slots=%zu natural_front=%zu md_front=%zu order_peak_bytes=%zu pair_checks=%zu cpu=%.6f s\n",n,side,s->slots,o->slots,s->max_front,o->max_front,o->peak_work_bytes,o->pair_checks,t);
 if(side){CHECK(o->slots<s->slots);}else CHECK(o->slots==s->slots&&o->max_front==2);
 if(numeric){p56_context *ctx=NULL;p61_stats st={0};start=clock();CHECK(!p62_factorize(&ctx,n,z,col,row,a,&sl,&nl,&st));double natural=(double)(clock()-start)/CLOCKS_PER_SEC;verify_solve(ctx,n,col,row,a);p56_destroy(ctx);ctx=NULL;p63_stats os={0};start=clock();CHECK(!p63_factorize(&ctx,n,z,col,row,a,&ol,&sl,&nl,&os));double md=(double)(clock()-start)/CLOCKS_PER_SEC;CHECK(os.numeric.factor_slots==o->slots&&os.numeric.peak_front==o->max_front);verify_solve(ctx,n,col,row,a);p56_destroy(ctx);printf("factor n=%u natural_cpu=%.6f s md_cpu=%.6f s natural_numeric_bytes=%zu md_numeric_bytes=%zu reorder_bytes=%zu\n",n,natural,md,st.peak_work_bytes,os.numeric.peak_work_bytes,os.permutation_bytes);}
 p62_free(s);p63_free(o);free(col);free(row);free(a);tests++;
}
int main(void){for(unsigned i=0;i<700;i++)random_case(i);for(unsigned i=0;i<160;i++)paired(i);bad();large(152592,0,1);large(1600,40,1);large(16384,128,0);printf("ordering %u/%u PASS max_x=%.3e max_componentwise_backward=%.3Le\n",tests,tests,worst_x,worst);return 0;}
