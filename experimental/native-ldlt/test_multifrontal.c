#include "multifrontal_ldlt.h"
#include "factor_cache.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#define ROOT UINT32_MAX
#define REQUIRE(x) do{if(!(x)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#x);exit(1);}}while(0)
static unsigned tests=0;static size_t delays=0,pairs=0;static long double worst=0;static double worst_x=0;
static uint32_t rng=19;static uint32_t rnd(void){rng=rng*1664525u+1013904223u;return rng;}
static p61_limits limits={256,128u*1024*1024,10000000,1000000,0};
static void solve(p56_context *ctx,uint32_t n,const size_t *col,const uint32_t *row,const double *a){
 double *b=calloc(n,sizeof *b),*x=malloc(n*sizeof *x),*y=malloc(n*sizeof *y),*scratch=malloc(n*sizeof *scratch);long double *r=calloc(n,sizeof *r),*den=calloc(n,sizeof *den);REQUIRE(b&&x&&y&&scratch&&r&&den);
 for(uint32_t i=0;i<n;i++)x[i]=((int)(i%13)-6)/7.0;
 for(uint32_t j=0;j<n;j++)for(size_t k=col[j];k<col[j+1];k++){uint32_t i=row[k];b[i]+=a[k]*x[j];if(i!=j)b[j]+=a[k]*x[i];}
 REQUIRE(p56_apply(ctx,b,y,scratch,n)==0);
 for(uint32_t i=0;i<n;i++){r[i]=-b[i];den[i]=fabsl(b[i]);double d=fabs(y[i]-x[i]);if(d>worst_x)worst_x=d;REQUIRE(d<1e-9);}
 for(uint32_t j=0;j<n;j++)for(size_t k=col[j];k<col[j+1];k++){uint32_t i=row[k];long double t=(long double)a[k]*y[j];r[i]+=t;den[i]+=fabsl(t);if(i!=j){t=(long double)a[k]*y[i];r[j]+=t;den[j]+=fabsl(t);}}
 for(uint32_t i=0;i<n;i++){long double d=den[i]?fabsl(r[i])/den[i]:0;if(d>worst)worst=d;REQUIRE(d<1e-12);}
 free(b);free(x);free(y);free(scratch);free(r);free(den);
}
static p61_stats run(uint32_t n,const size_t *col,const uint32_t *row,const double *a,const uint32_t *p,const p61_limits *lim){p56_context *c=NULL;p61_stats s={0};int rc=p61_factorize(&c,n,col,row,a,p,lim,&s);if(rc)fprintf(stderr,"rc=%d n=%u test=%u\n",rc,n,tests);REQUIRE(rc==0);REQUIRE(s.fronts==n&&s.scalars+2*s.pairs==n);solve(c,n,col,row,a);p56_destroy(c);tests++;delays+=s.delays;pairs+=s.pairs;return s;}
static void tree_case(unsigned t){
 uint32_t n=2+rnd()%47,p[48],row[48*48];size_t col[49],z=0;double a[48*48];double scale=t%3==0?1e-100:t%3==1?1e100:1;
 for(uint32_t j=0;j<n-1;j++){p[j]=j+1+rnd()%(n-j-1);}p[n-1]=ROOT;
 for(uint32_t j=0;j<n;j++){col[j]=z;row[z]=j;a[z++]=(j%2?1:-1)*scale*(double)(2*n);for(uint32_t i=p[j];i!=ROOT;i=p[i])if((rnd()>>16)%3!=0){row[z]=i;a[z++]=scale*((int)((rnd()>>16)%9)-4)/8.;}}col[n]=z;
 run(n,col,row,a,p,&limits);
}
/* Paired zero diagonals, a chain of fronts, additional weak couplings
   between adjacent pairs; exact solution obtained independently by CSC matvec. */
static void paired(unsigned t){uint32_t n=2*(2+rnd()%20),p[44],row[44*4];size_t col[45],z=0;double a[44*4],scale=t%3==0?1e-100:t%3==1?1e100:1;
 for(uint32_t j=0;j<n;j++){p[j]=j+1<n?j+1:ROOT;col[j]=z;row[z]=j;a[z++]=0;if(j%2==0){row[z]=j+1;a[z++]=scale*(2+(rnd()>>16)%5);}if(j+2<n){row[z]=j+2;a[z++]=scale*0.125;} }col[n]=z;
 p61_stats s=run(n,col,row,a,p,&limits);REQUIRE(s.pairs>0&&s.delays>0);
}
static void extra(void){
 /* Two delayed children arrive at a common parent; a separate root is
    processed independently. Alternative scalar pivot at the parent. */
 size_t col[]={0,2,4,5,6};uint32_t row[]={0,2,1,2,2,3},p[]={2,2,ROOT,ROOT};double a[]={0.01,2,0.02,3,20,-5};
 p61_stats st=run(4,col,row,a,p,&limits);REQUIRE(st.delays==2&&st.peak_front==3);
 /* Grandchild contribution skips the intermediate original column. */
 size_t cc[]={0,2,3,4};uint32_t rr[]={0,2,1,2},pp[]={1,2,ROOT};double aa[]={2,1,4,-3};run(3,cc,rr,aa,pp,&limits);
 /* No explicit diagonals: a nonsingular off-diagonal block. */
 size_t c2[]={0,1,1};uint32_t r2[]={1},p2[]={1,ROOT};double a2[]={2};run(2,c2,r2,a2,p2,&limits);
}
static void errors(void){size_t col[]={0,2,3};uint32_t row[]={0,1,1},p[]={1,ROOT};double a[]={0,2,0};p56_context *c=NULL;p61_stats s={0};REQUIRE(!p61_factorize(&c,2,col,row,a,p,&limits,&s));p56_context *old=c;p61_stats before=s;
 p61_limits l=limits;l.max_front=1;REQUIRE(p61_factorize(&c,2,col,row,a,p,&l,&s)==P60_LIMIT);tests++;
 l=limits;l.max_work_bytes=8;REQUIRE(p61_factorize(&c,2,col,row,a,p,&l,&s)==P60_LIMIT);tests++;
 l=limits;l.max_delays=0;REQUIRE(p61_factorize(&c,2,col,row,a,p,&l,&s)==P60_LIMIT);tests++;
 l=limits;l.max_factor_slots=0;REQUIRE(p61_factorize(&c,2,col,row,a,p,&l,&s)==P60_LIMIT);tests++;
 p[0]=ROOT;REQUIRE(p61_factorize(&c,2,col,row,a,p,&limits,&s)==P55_INVALID);tests++;p[0]=1;
 row[1]=0;REQUIRE(p61_factorize(&c,2,col,row,a,p,&limits,&s)==P55_INVALID);tests++;row[1]=1;
 a[1]=NAN;REQUIRE(p61_factorize(&c,2,col,row,a,p,&limits,&s)==P55_INVALID);tests++;
 a[1]=0;REQUIRE(p61_factorize(&c,2,col,row,a,p,&limits,&s)==P58_NO_PIVOT);tests++;a[1]=2;
 REQUIRE(c==old&&!memcmp(&before,&s,sizeof s));solve(c,2,col,row,a);tests++;
 /* Exhaust workspace budgets below the measured successful peak, covering
    allocation failures during initialization, assembly and packing. */
 for(size_t cap=1;cap<before.peak_work_bytes;cap+=13){l=limits;l.max_work_bytes=cap;REQUIRE(p61_factorize(&c,2,col,row,a,p,&l,&s)==P60_LIMIT);REQUIRE(c==old&&!memcmp(&before,&s,sizeof s));tests++;}
 unsigned char key[32]={0};p56_context *restored=NULL;const char *file="stage61-test-cache.bin";REQUIRE(!p57_save(c,file,key));REQUIRE(!p57_load(&restored,file,key,100000,2));solve(restored,2,col,row,a);remove(file);p56_destroy(restored);p56_destroy(c);tests++;
}
static void large(uint32_t n){size_t *col=malloc(((size_t)n+1)*sizeof *col);uint32_t *row=malloc(2*(size_t)n*sizeof *row),*p=malloc(n*sizeof *p);double *a=malloc(2*(size_t)n*sizeof *a);REQUIRE(col&&row&&p&&a);size_t z=0;
 for(uint32_t j=0;j<n;j++){col[j]=z;row[z]=j;a[z++]=4;p[j]=j+1<n?j+1:ROOT;if(j+1<n){row[z]=j+1;a[z++]=-1;}}col[n]=z;
 clock_t t=clock();p61_stats s={0};p56_context *ctx=NULL;REQUIRE(!p61_factorize(&ctx,n,col,row,a,p,&limits,&s));double ft=(double)(clock()-t)/CLOCKS_PER_SEC;t=clock();solve(ctx,n,col,row,a);double vt=(double)(clock()-t)/CLOCKS_PER_SEC;p56_destroy(ctx);tests++;
 printf("large n=%u nnz_lower=%zu front=%zu slots=%zu workspace=%zu cpu_factor=%.6f s cpu_solve_and_verify=%.6f s\n",n,z,s.peak_front,s.factor_slots,s.peak_work_bytes,ft,vt);REQUIRE(s.peak_front==2&&s.factor_slots==n-1);free(col);free(row);free(p);free(a);
}
static void grid_case(uint32_t side){uint32_t n=side*side;size_t *col=malloc(((size_t)n+1)*sizeof *col);uint32_t *row=malloc(3*(size_t)n*sizeof *row),*p=malloc(n*sizeof *p);double *a=malloc(3*(size_t)n*sizeof *a);REQUIRE(col&&row&&p&&a);size_t z=0;
 for(uint32_t j=0;j<n;j++){col[j]=z;row[z]=j;a[z++]=5;p[j]=j+1<n?j+1:ROOT;if(j%side+1<side){row[z]=j+1;a[z++]=-1;}if(j+side<n){row[z]=j+side;a[z++]=-1;}}col[n]=z;
 p61_stats s=run(n,col,row,a,p,&limits);printf("grid %ux%u n=%u nnz_lower=%zu front=%zu slots=%zu workspace=%zu\n",side,side,n,z,s.peak_front,s.factor_slots,s.peak_work_bytes);REQUIRE(s.peak_front==side+1);free(col);free(row);free(p);free(a);
}
static void large_indefinite(uint32_t n){size_t *col=malloc(((size_t)n+1)*sizeof *col);uint32_t *row=malloc(2*(size_t)n*sizeof *row),*p=malloc(n*sizeof *p);double *a=malloc(2*(size_t)n*sizeof *a);REQUIRE(col&&row&&p&&a);size_t z=0;
 for(uint32_t j=0;j<n;j++){col[j]=z;row[z]=j;a[z++]=0;p[j]=j%2?ROOT:j+1;if(j%2==0){row[z]=j+1;a[z++]=2;}}col[n]=z;
 p61_stats s=run(n,col,row,a,p,&limits);printf("large_indefinite n=%u nnz_lower=%zu front=%zu pairs=%zu delays=%zu slots=%zu workspace=%zu\n",n,z,s.peak_front,s.pairs,s.delays,s.factor_slots,s.peak_work_bytes);REQUIRE(s.peak_front==2&&s.pairs==n/2&&s.delays==n/2);free(col);free(row);free(p);free(a);
}
int main(void){for(unsigned i=0;i<300;i++)tree_case(i);for(unsigned i=0;i<200;i++)paired(i);extra();errors();large(10000);large(152592);grid_case(32);large_indefinite(152592);printf("multifrontal %u/%u PASS, pairs=%zu delays=%zu max_x=%.3e max_componentwise_backward=%.3Le\n",tests,tests,pairs,delays,worst_x,worst);return 0;}
