#include "parallel_ldlt.h"
#include "schur_simd.h"
#include <omp.h>
#include <fenv.h>
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
static size_t paired_count,delay_count,parallel_runs;
static void identical(p56_context *x,p56_context *y,uint32_t n){
 double *b=malloc((size_t)n*sizeof *b),*a=malloc((size_t)n*sizeof *a),*c=malloc((size_t)n*sizeof *c),*w=malloc((size_t)n*sizeof *w);CHECK(b&&a&&c&&w);
 for(uint32_t j=0;j<n;j++)b[j]=((int)(j%19)-9)/7.;
 CHECK(!p56_apply(x,b,a,w,n)&&!p56_apply(y,b,c,w,n));CHECK(!memcmp(a,c,(size_t)n*sizeof *a));free(b);free(a);free(c);free(w);
}
static void random_case(unsigned t){
 uint32_t n=12+rnd()%49,ri[4096];size_t cp[65],z=0;double a[4096],scale=t%3==0?1e-100:t%3==1?1e100:1.;
 /* Sparse connected branches or disconnected components. Strict diagonal dominance. */
 for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=(j%2?1.:-1.)*scale*(2*n+1);for(uint32_t i=j+1;i<n;i++)if((t%2==0?i%4==j%4:i==n-1||i%4==j%4)&&rnd()%100<25+t%60){ri[z]=i;a[z++]=scale*((int)(rnd()%15)-7)/8.;}}cp[n]=z;
 p56_context *serial=NULL,*par=NULL;p61_stats ref={0};p67_stats st={0};uint32_t w=1+t%12;unsigned threads=2+t%3;
 CHECK(!p65_factorize(&serial,n,z,cp,ri,a,&sl,&nl,w,&ref));CHECK(!p67_factorize(&par,n,z,cp,ri,a,&sl,&nl,w,threads,&st));
 CHECK(st.numeric.pairs==ref.pairs&&st.numeric.scalars==ref.scalars&&st.numeric.factor_slots==ref.factor_slots&&st.numeric.delays==ref.delays&&st.numeric.fronts==ref.fronts);identical(serial,par,n);verify_solve(par,n,cp,ri,a);parallel_runs+=st.threads_used>1;tests++;
 if(t<100){p64_stats os;CHECK(!p65_ordered_factorize(&serial,n,z,cp,ri,a,&lim,&sl,&nl,(p64_method)(t%4),w,&os));CHECK(!p67_ordered_factorize(&par,n,z,cp,ri,a,&lim,&sl,&nl,(p64_method)(t%4),w,threads,&os,&st));identical(serial,par,n);verify_solve(par,n,cp,ri,a);tests++;}
 p56_destroy(serial);p56_destroy(par);
}
static void paired(unsigned t){
 uint32_t n=2*(4+rnd()%20),ri[4096];size_t cp[65],z=0;double a[4096],scale=t%3==0?1e-100:t%3==1?1e100:1.;uint32_t half=n/2;
 /* Match across distant columns, multiple leaves, delay into separate parents. */
 for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=0;if(j<half){ri[z]=j+half;a[z++]=scale*(1+j%5);}else if(t%2&&j+1<n){ri[z]=j+1;a[z++]=scale*.01;}}cp[n]=z;
 p56_context *serial=NULL,*par=NULL;p61_stats ref={0};p67_stats st={0};
 CHECK(!p65_factorize(&serial,n,z,cp,ri,a,&sl,&nl,1,&ref));CHECK(!p67_factorize(&par,n,z,cp,ri,a,&sl,&nl,1,4,&st));
 CHECK(st.numeric.pairs==ref.pairs&&st.numeric.delays==ref.delays&&st.numeric.pairs==half);identical(serial,par,n);verify_solve(par,n,cp,ri,a);paired_count+=st.numeric.pairs;delay_count+=st.numeric.delays;parallel_runs+=st.threads_used>1;tests++;p56_destroy(serial);p56_destroy(par);
}
static void special(void){
 size_t cp[]={0,2,4,5,6};uint32_t ri[]={0,2,1,3,2,3};double a[]={0,2,0,3,0,0};p56_context *p=NULL,*q=NULL;p67_stats st={0};p61_stats ref;
 CHECK(!p65_factorize(&q,4,6,cp,ri,a,&sl,&nl,1,&ref));CHECK(!p67_factorize(&p,4,6,cp,ri,a,&sl,&nl,1,4,&st));CHECK(st.threads_used>1&&st.numeric.delays==2&&st.numeric.pairs==2);identical(q,p,4);tests++;
 p56_context *old=p;p67_stats before=st;p61_limits low=nl;low.max_delays=1;CHECK(p67_factorize(&p,4,6,cp,ri,a,&sl,&low,1,4,&st)==P60_LIMIT&&p==old&&!memcmp(&st,&before,sizeof st));tests++;
 low=nl;low.max_factor_slots=1;CHECK(p67_factorize(&p,4,6,cp,ri,a,&sl,&low,1,4,&st)==P60_LIMIT&&p==old);tests++;
 low=nl;low.max_front=1;CHECK(p67_factorize(&p,4,6,cp,ri,a,&sl,&low,1,4,&st)==P60_LIMIT&&p==old);tests++;
 for(size_t b=1;b<before.numeric.peak_work_bytes;b+=23){low=nl;low.max_work_bytes=b;int rc=p67_factorize(&p,4,6,cp,ri,a,&sl,&low,1,4,&st);if(rc){CHECK(rc==P60_LIMIT&&p==old&&!memcmp(&st,&before,sizeof st));}else{CHECK(st.numeric.peak_work_bytes<=b);old=p;before=st;}tests++;}
 CHECK(p67_factorize(&p,4,6,cp,ri,a,&sl,&nl,1,0,&st)==P55_INVALID&&p==old);tests++;
 a[1]=0;CHECK(p67_factorize(&p,4,6,cp,ri,a,&sl,&nl,1,4,&st)==P58_NO_PIVOT&&p==old);tests++;a[1]=2;
 for(int m=1;m<=3;m++)if(p66_available((p66_mode)m)){CHECK(p66_set_mode((p66_mode)m));CHECK(!p65_factorize(&q,4,6,cp,ri,a,&sl,&nl,1,&ref));CHECK(!p67_factorize(&p,4,6,cp,ri,a,&sl,&nl,1,4,&st));identical(q,p,4);tests++;}
 CHECK(p66_set_mode(P66_AUTO));
 /* Concurrent solves on immutable factor, individual scratch and outputs. */
 int ok=1;
 #pragma omp parallel for num_threads(4) reduction(&:ok)
 for(int i=0;i<64;i++){double b[4]={2,3,2,3},x[4],w[4];ok&=!p56_apply(p,b,x,w,4)&&x[0]==1&&x[1]==1&&x[2]==1&&x[3]==1;}
 CHECK(ok);tests++;p56_destroy(q);p56_destroy(p);
}
static void bench(uint32_t n,uint32_t side,uint32_t block,int repetitions){
 size_t *cp=malloc(((size_t)n+1)*sizeof *cp);uint32_t *ri=malloc((size_t)n*(block?block:3)*sizeof *ri);double *a=malloc((size_t)n*(block?block:3)*sizeof *a);CHECK(cp&&ri&&a);size_t z=0;
 for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=block?2*block:5;if(block){uint32_t end=(j/block+1)*block;if(end>n)end=n;for(uint32_t i=j+1;i<end;i++){ri[z]=i;a[z++]=-.25;}}else{if(j+1<n&&(!side||j%side+1<side)){ri[z]=j+1;a[z++]=-1;}if(side&&j+side<n){ri[z]=j+side;a[z++]=-1;}}}cp[n]=z;
 p56_context *ref=NULL;for(unsigned th=1;th<=4;th*=2){double times[5];p67_stats st={0};p64_stats os;for(int rep=0;rep<repetitions;rep++){p56_context *ctx=NULL;double start=omp_get_wtime();CHECK(!p67_ordered_factorize(&ctx,n,z,cp,ri,a,&lim,&sl,&nl,side?P64_ND:P64_NATURAL,16,th,&os,&st));times[rep]=omp_get_wtime()-start;verify_solve(ctx,n,cp,ri,a);if(!ref){ref=ctx;ctx=NULL;}else identical(ref,ctx,n);p56_destroy(ctx);}for(int i=0;i<repetitions;i++)for(int j=i+1;j<repetitions;j++)if(times[j]<times[i]){double x=times[i];times[i]=times[j];times[j]=x;}
 printf("bench n=%u side=%u block=%u requested=%u used=%u levels=%zu ready=%zu workspace=%zu wall_median=%.6f s\n",n,side,block,th,st.threads_used,st.levels,st.max_ready,st.numeric.peak_work_bytes,times[repetitions/2]);tests++;}
 p56_destroy(ref);free(cp);free(ri);free(a);
}
int main(int argc,char **argv){omp_set_dynamic(0);if(argc>1&&!strcmp(argv[1],"bench")){bench(16384,128,0,3);bench(8192,0,64,3);bench(152592,0,0,3);return 0;}for(unsigned i=0;i<240;i++)random_case(i);for(unsigned i=0;i<120;i++)paired(i);special();bench(4096,64,0,1);bench(4096,0,64,1);bench(152592,0,0,1);printf("parallel %u/%u PASS max_x=%.3e max_componentwise_backward=%.3Le parallel_random_runs=%zu pairs=%zu delays=%zu\n",tests,tests,worst_x,worst,parallel_runs,paired_count,delay_count);return 0;}
