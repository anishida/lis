#include "symbolic_ldlt.h"
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
/* Independent graph oracle: explicitly form all fill cliques. */
static void oracle(uint32_t n,const size_t *col,const uint32_t *row,const p62_symbolic *s){
 unsigned char g[64][64]={{0}};CHECK(n<=64);size_t slots=0,roots=0,front=0,ns=1;CHECK(s->super[0]==0);
 for(uint32_t j=0;j<n;j++)for(size_t z=col[j];z<col[j+1];z++)g[row[z]][j]=g[j][row[z]]=1;
 for(uint32_t j=0;j<n;j++){
  uint32_t p=NONE,c=1;for(uint32_t i=j+1;i<n;i++)if(g[i][j]){if(p==NONE)p=i;c++;}
  CHECK(s->parent[j]==p&&s->count[j]==c);slots+=c-1;roots+=p==NONE;if(c>front)front=c;
  for(uint32_t i=j+1;i<n;i++)if(g[i][j])for(uint32_t k=i+1;k<n;k++)if(g[k][j])g[k][i]=g[i][k]=1;
 }
 for(uint32_t j=0;j+1<n;j++){
  int merge=g[j+1][j];for(uint32_t i=j+2;i<n;i++)if(g[i][j]!=g[i][j+1])merge=0;
  if(!merge){CHECK(ns<s->nsuper&&s->super[ns]==j+1);ns++;}
 }
 CHECK(s->nsuper==ns&&s->super[ns]==n&&s->slots==slots&&s->roots==roots&&s->max_front==front);
}
static void random_case(unsigned t){
 uint32_t n=1+rnd()%64,row[64*64];size_t col[65],z=0;double a[64*64],scale=t%3==0?1e-100:t%3==1?1e100:1.;
 for(uint32_t j=0;j<n;j++){col[j]=z;row[z]=j;a[z++]=(j%2?1.:-1.)*scale*(2*n+1);for(uint32_t i=j+1;i<n;i++)if(rnd()%100<t%100){row[z]=i;a[z++]=scale*((int)(rnd()%15)-7)/8.;}}
 col[n]=z;p62_symbolic *s=NULL;CHECK(!p62_analyze(&s,n,z,col,row,&sl));oracle(n,col,row,s);p62_free(s);
 if(t<240){p56_context *ctx=NULL;p61_stats st={0};CHECK(!p62_factorize(&ctx,n,z,col,row,a,&sl,&nl,&st));verify_solve(ctx,n,col,row,a);p56_destroy(ctx);CHECK(st.fronts==n&&st.scalars==n);}
 tests++;
}
static void paired_case(unsigned t){uint32_t n=2*(1+rnd()%30),row[128];size_t col[65],z=0;double a[128],scale=t%2?1e100:1e-100;
 for(uint32_t j=0;j<n;j++){col[j]=z;if(t%2){row[z]=j;a[z++]=0;}if(j<n/2){row[z]=j+n/2;a[z++]=scale*(2+rnd()%7);}}col[n]=z;
 p62_symbolic *s=NULL;CHECK(!p62_analyze(&s,n,z,col,row,&sl));oracle(n,col,row,s);CHECK(s->roots==n/2&&s->slots==n/2);p62_free(s);
 p56_context *ctx=NULL;p61_stats st={0};CHECK(!p62_factorize(&ctx,n,z,col,row,a,&sl,&nl,&st));CHECK(st.pairs==n/2&&st.delays==n/2);verify_solve(ctx,n,col,row,a);p56_destroy(ctx);tests++;
}
static void bad(void){size_t col[]={0,2,3};uint32_t row[]={0,1,1};double a[]={0,2,0};p62_symbolic *s=NULL;CHECK(!p62_analyze(&s,2,3,col,row,&sl));p62_symbolic *old=s;p62_limits l=sl;
 #define BAD(expr,code) do{CHECK((expr)==(code));CHECK(s==old);tests++;}while(0)
 BAD(p62_analyze(&s,2,2,col,row,&sl),P55_INVALID);
 BAD(p62_analyze(&s,2,3,col,NULL,&sl),P55_INVALID);
 row[1]=0;BAD(p62_analyze(&s,2,3,col,row,&sl),P55_INVALID);row[1]=1;
 row[2]=0;BAD(p62_analyze(&s,2,3,col,row,&sl),P55_INVALID);row[2]=1;
 row[1]=2;BAD(p62_analyze(&s,2,3,col,row,&sl),P55_INVALID);row[1]=1;
 col[1]=4;BAD(p62_analyze(&s,2,3,col,row,&sl),P55_INVALID);col[1]=2;
 l.max_slots=0;BAD(p62_analyze(&s,2,3,col,row,&l),P60_LIMIT);
 l=sl;for(size_t i=0;i<old->peak_work_bytes;i+=7){l.max_work_bytes=i;BAD(p62_analyze(&s,2,3,col,row,&l),P60_LIMIT);}
 l.max_work_bytes=old->peak_work_bytes;CHECK(!p62_analyze(&s,2,3,col,row,&l));tests++;
 p62_free(s);s=NULL;size_t empty[]={0,0,0};CHECK(!p62_analyze(&s,2,0,empty,NULL,&sl));CHECK(s->roots==2&&s->slots==0&&s->nsuper==2);p62_free(s);tests++;
 p56_context *ctx=NULL;p61_stats st={0};CHECK(!p62_factorize(&ctx,2,3,col,row,a,&sl,&nl,&st));p56_context *before=ctx; p61_stats bs=st;
 l=sl;l.max_slots=0;CHECK(p62_factorize(&ctx,2,3,col,row,a,&l,&nl,&st)==P60_LIMIT);CHECK(ctx==before&&!memcmp(&st,&bs,sizeof st));tests++;
 p61_limits lim=nl;lim.max_front=1;CHECK(p62_factorize(&ctx,2,3,col,row,a,&sl,&lim,&st)==P60_LIMIT);CHECK(ctx==before&&!memcmp(&st,&bs,sizeof st));tests++;
 verify_solve(ctx,2,col,row,a);p56_destroy(ctx);tests++;
}
static void large(uint32_t n,uint32_t side,int numeric){
 size_t *col=malloc(((size_t)n+1)*sizeof *col);uint32_t *row=malloc(3*(size_t)n*sizeof *row);double *a=malloc(3*(size_t)n*sizeof *a);CHECK(col&&row&&a);size_t z=0;
 for(uint32_t j=0;j<n;j++){col[j]=z;row[z]=j;a[z++]=5;if(j+1<n&&(!side||j%side+1<side)){row[z]=j+1;a[z++]=-1;}if(side&&j+side<n){row[z]=j+side;a[z++]=-1;}}col[n]=z;
 clock_t start=clock();p62_symbolic *s=NULL;CHECK(!p62_analyze(&s,n,z,col,row,&sl));double t=(double)(clock()-start)/CLOCKS_PER_SEC;
 CHECK(s->roots==1&&s->max_front==(side?side+1:2));if(!side)CHECK(s->slots==n-1&&s->nsuper==n-1);
 printf("symbolic n=%u grid_side=%u nnz_lower=%zu slots=%zu supernodes=%zu max_front=%zu peak_bytes=%zu cpu=%.6f s\n",n,side,z,s->slots,s->nsuper,s->max_front,s->peak_work_bytes,t);
 if(numeric){p56_context *ctx=NULL;p61_stats st={0};start=clock();CHECK(!p62_factorize(&ctx,n,z,col,row,a,&sl,&nl,&st));double ft=(double)(clock()-start)/CLOCKS_PER_SEC;CHECK(st.peak_front==s->max_front&&st.factor_slots==s->slots);verify_solve(ctx,n,col,row,a);p56_destroy(ctx);printf("auto_factor n=%u cpu_symbolic_numeric=%.6f s numeric_peak_bytes=%zu\n",n,ft,st.peak_work_bytes);}
 p62_free(s);free(col);free(row);free(a);tests++;
}
int main(void){for(unsigned i=0;i<700;i++)random_case(i);for(unsigned i=0;i<160;i++)paired_case(i);bad();large(152592,0,1);large(1600,40,1);large(16384,128,0);printf("symbolic %u/%u PASS max_x=%.3e max_componentwise_backward=%.3Le\n",tests,tests,worst_x,worst);return 0;}
