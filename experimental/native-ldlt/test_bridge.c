#include "buckle_bridge.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#define CHECK(x) do { if(!(x)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#x);exit(1);} }while(0)
static unsigned tests;static double max_error;static long double max_res;static uint32_t rng=69121;
static uint32_t rnd(void){rng^=rng<<13;rng^=rng>>17;rng^=rng<<5;return rng;}
typedef struct {p69_csr m;size_t *p;uint32_t *j;double *v;} owned;
static owned make(uint32_t n,const double *a,p69_storage s,unsigned base){
 owned x={0};x.p=malloc(((size_t)n+1)*sizeof *x.p);x.j=malloc((size_t)n*n*sizeof *x.j);x.v=malloc((size_t)n*n*sizeof *x.v);CHECK(x.p&&x.j&&x.v);size_t z=0;
 for(uint32_t i=0;i<n;i++){x.p[i]=z+base;for(uint32_t j=0;j<n;j++)if(a[(size_t)i*n+j]!=0&&(s==P69_FULL||(s==P69_LOWER&&j<=i)||(s==P69_UPPER&&j>=i))){x.j[z]=j+base;x.v[z++]=a[(size_t)i*n+j];}}x.p[n]=z+base;x.m=(p69_csr){n,z,x.p,x.j,x.v,base,s};return x;
}
static void drop(owned x){free(x.p);free(x.j);free(x.v);}
static void verify(p56_context *f,uint32_t n,const double *k,const double *g,double shift){
 double *b=calloc(n,sizeof *b),*x=malloc(n*sizeof *x),*w=malloc(n*sizeof *w),*y=malloc(n*sizeof *y);CHECK(b&&x&&w&&y);
 for(uint32_t i=0;i<n;i++)x[i]=((int)(i%13)-6)/7.;
 for(uint32_t i=0;i<n;i++){long double t=0;for(uint32_t j=0;j<n;j++)t+=((long double)k[(size_t)i*n+j]-(long double)shift*g[(size_t)i*n+j])*x[j];b[i]=(double)t;}
 CHECK(!p56_apply(f,b,y,w,n));for(uint32_t i=0;i<n;i++){double e=fabs(x[i]-y[i]);CHECK(e<1e-11);if(e>max_error)max_error=e;long double r=-b[i],d=fabsl(b[i]);for(uint32_t j=0;j<n;j++){long double t=((long double)k[(size_t)i*n+j]-(long double)shift*g[(size_t)i*n+j])*y[j];r+=t;d+=fabsl(t);}long double e2=d?fabsl(r)/d:0;CHECK(e2<1e-12L);if(e2>max_res)max_res=e2;}
 free(b);free(x);free(w);free(y);
}
static void random_case(unsigned t){uint32_t n=2+2*(rnd()%15);double k[32*32]={0},g[32*32]={0};double scale=t%3==0?1e-100:t%3==1?1e100:1,shift=t%4==0?0:t%4==1?.5:t%4==2?-2:2;
 for(uint32_t i=0;i<n;i++){k[(size_t)i*n+i]=(i%2?1:-1)*scale*(8*n);g[(size_t)i*n+i]=scale;for(uint32_t j=i+1;j<n;j++){if(rnd()%3==0)k[(size_t)i*n+j]=k[(size_t)j*n+i]=scale*((int)(rnd()%9)-4)/4.;if(rnd()%4==0)g[(size_t)i*n+j]=g[(size_t)j*n+i]=scale*((int)(rnd()%9)-4)/4.;}}
 owned a=make(n,k,(p69_storage)(t%3),t%2),b=make(n,g,(p69_storage)((t/3)%3),(t/2)%2);p69_options o;p69_defaults(&o);o.method=(p64_method)(t%4);p56_context *f=NULL;p69_stats st;CHECK(!p69_factorize(&f,&a.m,&b.m,shift,&o,&st));CHECK(st.parallel.threads_used==1&&st.shifted_nnz<=a.m.nnz+b.m.nnz);verify(f,n,k,g,shift);p56_destroy(f);drop(a);drop(b);tests++;}
static void pairs(void){for(unsigned t=0;t<60;t++){uint32_t n=2+2*(t%15);double k[32*32]={0},g[32*32]={0};for(uint32_t i=0;i<n;i++){k[i*n+i]=g[i*n+i]=1;}for(uint32_t i=0;i<n/2;i++){uint32_t j=i+n/2;g[i*n+j]=g[j*n+i]=2+i%4;}owned a=make(n,k,(p69_storage)(t%3),t%2),b=make(n,g,(p69_storage)((t+1)%3),(t+1)%2);p69_options o;p69_defaults(&o);o.method=(p64_method)(t%4);p56_context *f=NULL;p69_stats st;CHECK(!p69_factorize(&f,&a.m,&b.m,1,&o,&st));CHECK(p56_block_count(f)==n/2);verify(f,n,k,g,1);p56_destroy(f);drop(a);drop(b);tests++;}}
static void invalid(void){double a[]={4,1,1,5},z[]={0,0,0,0};owned k=make(2,a,P69_FULL,0),g=make(2,a,P69_UPPER,1);p69_options o;p69_defaults(&o);p56_context *f=NULL;p69_stats st;CHECK(!p69_factorize(&f,&k.m,NULL,0,&o,&st));p56_context *old=f;size_t saved=st.shifted_nnz;tests++;
 #define BAD(expr) do{CHECK((expr)!=0&&f==old&&st.shifted_nnz==saved);verify(f,2,a,z,0);tests++;}while(0)
 BAD(p69_factorize(&f,&k.m,NULL,1,&o,&st));BAD(p69_factorize(&f,&k.m,&g.m,NAN,&o,&st));
 k.v[1]=2;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));k.v[1]=1;
 k.j[1]=0;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));k.j[1]=1;
 k.j[1]=2;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));k.j[1]=1;
 k.p[1]=5;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));k.p[1]=2;
 k.m.base=2;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));k.m.base=0;
 k.m.storage=P69_LOWER;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));k.m.storage=P69_FULL;
 k.v[0]=INFINITY;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));k.v[0]=4;
 g.m.n=1;BAD(p69_factorize(&f,&k.m,&g.m,1,&o,&st));g.m.n=2;
 BAD(p69_factorize(&f,&k.m,&g.m,DBL_MAX,&o,&st));BAD(p69_factorize(&f,&k.m,&g.m,1,&o,&st));
 o.max_conversion_bytes=1;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));o.max_conversion_bytes=256u*1024*1024;
 o.numeric.max_work_bytes=1;BAD(p69_factorize(&f,&k.m,NULL,0,&o,&st));
 #undef BAD
 p56_destroy(f);drop(k);drop(g);
}
static void large(void){uint32_t n=152592;size_t *p=malloc(((size_t)n+1)*sizeof *p),*q=malloc(((size_t)n+1)*sizeof *q);uint32_t *j=malloc((size_t)3*n*sizeof *j),*r=malloc((size_t)n*sizeof *r);double *v=malloc((size_t)3*n*sizeof *v),*s=malloc((size_t)n*sizeof *s),*x=malloc((size_t)n*sizeof *x),*b=malloc((size_t)n*sizeof *b),*y=malloc((size_t)n*sizeof *y),*w=malloc((size_t)n*sizeof *w);CHECK(p&&q&&j&&r&&v&&s&&x&&b&&y&&w);size_t z=0;for(uint32_t i=0;i<n;i++){p[i]=z;if(i){j[z]=i-1;v[z++]=-1;}j[z]=i;v[z++]=5;if(i+1<n){j[z]=i+1;v[z++]=-1;}q[i]=i+1;r[i]=i+1;s[i]=2;x[i]=((int)(i%17)-8)/9.;}p[n]=z;q[n]=(size_t)n+1;for(uint32_t i=0;i<n;i++)b[i]=4*x[i]-(i?x[i-1]:0)-(i+1<n?x[i+1]:0);p69_csr a={n,z,p,j,v,0,P69_FULL},g={n,n,q,r,s,1,P69_LOWER};p69_options o;p69_defaults(&o);p56_context *f=NULL;p69_stats st;double t=omp_get_wtime();CHECK(!p69_factorize(&f,&a,&g,.5,&o,&st));double dt=omp_get_wtime()-t;t=omp_get_wtime();CHECK(!p56_apply(f,b,y,w,n));double solve=omp_get_wtime()-t;for(uint32_t i=0;i<n;i++){double e=fabs(x[i]-y[i]);CHECK(e<1e-12);if(e>max_error)max_error=e;long double r0=4.L*y[i]-(i?y[i-1]:0)-(i+1<n?y[i+1]:0)-b[i],d=4.L*fabsl(y[i])+(i?fabsl(y[i-1]):0)+(i+1<n?fabsl(y[i+1]):0)+fabsl(b[i]);long double er=d?fabsl(r0)/d:0;CHECK(er<1e-12L);if(er>max_res)max_res=er;}
 printf("large n=%u input_k=%zu input_g=%zu shifted_nnz=%zu conversion_bytes=%zu numeric_bytes=%zu threads=%u factorize_wall=%.6f s apply_wall=%.6f s\n",n,z,g.nnz,st.shifted_nnz,st.conversion_bytes,st.parallel.numeric.peak_work_bytes,st.parallel.threads_used,dt,solve);p56_destroy(f);free(p);free(q);free(j);free(r);free(v);free(s);free(x);free(b);free(y);free(w);tests++;
}
int main(void){for(unsigned t=0;t<360;t++)random_case(t);pairs();invalid();large();printf("bridge %u/%u PASS max_x=%.3e max_componentwise_backward=%.3Le\n",tests,tests,max_error,max_res);return 0;}
