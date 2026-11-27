#include "nested_dissection.h"
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
static void oracle(uint32_t n,const size_t *cp,const uint32_t *ri,const p63_order *o,const p64_estimate *e){
 unsigned char g[64][64]={{0}},seen[64]={0};size_t slots=0,front=0;long double work=0;CHECK(n<=64);
 for(uint32_t j=0;j<n;j++){CHECK(o->perm[j]<n&&!seen[o->perm[j]]&&o->inv[o->perm[j]]==j);seen[o->perm[j]]=1;for(size_t z=cp[j];z<cp[j+1];z++)if(ri[z]!=j)g[o->inv[j]][o->inv[ri[z]]]=g[o->inv[ri[z]]][o->inv[j]]=1;}
 for(uint32_t j=0;j<n;j++){size_t d=1;for(uint32_t i=j+1;i<n;i++)if(g[i][j]){d++;for(uint32_t k=i+1;k<n;k++)if(g[k][j])g[i][k]=g[k][i]=1;}slots+=d-1;if(d>front)front=d;work+=(long double)d*d;}
 CHECK(slots==e->slots&&front==e->front&&work==e->work&&o->slots==slots&&o->max_front==front);
}
static void random_case(unsigned t){
 uint32_t n=1+rnd()%64,ri[4096];size_t cp[65],z=0;double a[4096],scale=t%3==0?1e-100:t%3==1?1e100:1.;
 for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=(j%2?1.:-1.)*scale*(2*n+1);for(uint32_t i=j+1;i<n;i++)if(rnd()%100<t%100){ri[z]=i;a[z++]=scale*((int)(rnd()%15)-7)/8.;}}cp[n]=z;
 for(int m=P64_ND;m<=P64_AUTO;m++){p63_order *o=NULL;p64_stats st={0};CHECK(!p64_order(n,z,cp,ri,&lim,&sl,m,&o,&st));oracle(n,cp,ri,o,&st.candidates[st.selected]);for(int i=0;i<3;i++)if(!st.candidates[i].status)CHECK(st.candidates[st.selected].work<=st.candidates[i].work);p63_free(o);if(t<200){p56_context *ctx=NULL;CHECK(!p64_factorize(&ctx,n,z,cp,ri,a,&lim,&sl,&nl,m,&st));verify_solve(ctx,n,cp,ri,a);p56_destroy(ctx);}tests++;}
}
static void paired(unsigned t){uint32_t n=2*(1+rnd()%30),ri[128];size_t cp[65],z=0;double a[128],scale=t%2?1e100:1e-100;
 for(uint32_t j=0;j<n;j++){cp[j]=z;if(t%2){ri[z]=j;a[z++]=0;}if(j<n/2){ri[z]=j+n/2;a[z++]=scale*(2+rnd()%7);}}cp[n]=z;
 p56_context *ctx=NULL;p64_stats st={0};CHECK(!p64_factorize(&ctx,n,z,cp,ri,a,&lim,&sl,&nl,P64_ND,&st));CHECK(st.numeric.pairs==n/2&&st.numeric.delays==n/2);verify_solve(ctx,n,cp,ri,a);p56_destroy(ctx);tests++;
}
static void limits(void){size_t cp[]={0,3,5,6};uint32_t ri[]={0,1,2,1,2,2};double a[]={5,1,1,5,1,5};p63_order *o=NULL;p64_stats st={0};CHECK(!p64_order(3,6,cp,ri,&lim,&sl,P64_ND,&o,&st));p63_order *old=o;p64_stats saved=st;p64_limits l=lim;
 l.max_front=2;CHECK(p64_order(3,6,cp,ri,&l,&sl,P64_AUTO,&o,&st)==P60_LIMIT&&o==old&&!memcmp(&saved,&st,sizeof st));tests++;
 l=lim;l.max_work_bytes=0;CHECK(p64_order(3,6,cp,ri,&l,&sl,P64_ND,&o,&st)==P60_LIMIT&&o==old);tests++;
 l=lim;l.max_pair_checks=0;CHECK(!p64_order(3,6,cp,ri,&l,&sl,P64_AUTO,&o,&st));CHECK(st.candidates[P64_MD].status==P60_LIMIT&&st.selected==P64_NATURAL);tests++;
 CHECK(p64_order(3,5,cp,ri,&lim,&sl,P64_ND,&o,&st)==P55_INVALID);tests++;p63_free(o);
 p56_context *ctx=NULL;CHECK(!p64_factorize(&ctx,3,6,cp,ri,a,&lim,&sl,&nl,P64_AUTO,&st));p56_context *before=ctx;saved=st;p61_limits low=nl;low.max_front=1;CHECK(p64_factorize(&ctx,3,6,cp,ri,a,&lim,&sl,&low,P64_AUTO,&st)==P60_LIMIT&&ctx==before&&!memcmp(&saved,&st,sizeof st));verify_solve(ctx,3,cp,ri,a);p56_destroy(ctx);tests++;
 size_t empty[]={0,0,0};o=NULL;CHECK(!p64_order(2,0,empty,NULL,&lim,&sl,P64_ND,&o,&st));CHECK(o->slots==0&&o->max_front==1);p63_free(o);tests++;
}
static void large(uint32_t n,uint32_t side,int numeric){
 size_t *cp=malloc(((size_t)n+1)*sizeof *cp);uint32_t *ri=malloc(3*(size_t)n*sizeof *ri);double *a=malloc(3*(size_t)n*sizeof *a);CHECK(cp&&ri&&a);size_t z=0;
 for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=5;if(j+1<n&&(!side||j%side+1<side)){ri[z]=j+1;a[z++]=-1;}if(side&&j+side<n){ri[z]=j+side;a[z++]=-1;}}cp[n]=z;
 for(int m=0;m<=3;m++){p63_order *o=NULL;p64_stats st={0};clock_t start=clock();CHECK(!p64_order(n,z,cp,ri,&lim,&sl,m,&o,&st));double t=(double)(clock()-start)/CLOCKS_PER_SEC;p64_estimate e=st.candidates[st.selected];printf("order n=%u side=%u method=%d selected=%d slots=%zu front=%zu work=%.0Lf phase_bytes=%zu cpu=%.6f s\n",n,side,m,st.selected,e.slots,e.front,e.work,e.work_bytes,t);if(!side&&m==P64_AUTO)CHECK(st.selected==P64_NATURAL&&st.candidates[1].status==-1&&st.candidates[2].status==-1);p63_free(o);
 if(numeric&&(side||m==P64_AUTO)){p56_context *ctx=NULL;start=clock();CHECK(!p64_factorize(&ctx,n,z,cp,ri,a,&lim,&sl,&nl,m,&st));double f=(double)(clock()-start)/CLOCKS_PER_SEC;verify_solve(ctx,n,cp,ri,a);printf("factor n=%u method=%d cpu=%.6f s numeric_bytes=%zu front=%zu slots=%zu\n",n,m,f,st.numeric.peak_work_bytes,st.numeric.peak_front,st.numeric.factor_slots);p56_destroy(ctx);}tests++;}
 free(cp);free(ri);free(a);
}
int main(void){for(unsigned i=0;i<400;i++)random_case(i);for(unsigned i=0;i<120;i++)paired(i);limits();large(1600,40,1);large(16384,128,0);large(152592,0,1);printf("nested-dissection %u/%u PASS max_x=%.3e max_componentwise_backward=%.3Le\n",tests,tests,worst_x,worst);return 0;}
