#include "schur_simd.h"
#include "supernodal_ldlt.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#include <time.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#x);exit(1);}}while(0)
static unsigned tests; static uint32_t rng=73164; static long double max_kernel,max_back;static double max_x;
static uint32_t rnd(void){rng^=rng<<13;rng^=rng>>17;rng^=rng<<5;return rng;}
static double v(void){return ((int)(rnd()%2001)-1000)/997.;}
static p62_limits sl={128u*1024*1024,100000000};
static p61_limits nl={512,128u*1024*1024,10000000,1000000,0};
static p64_limits ol={128u*1024*1024,100000000,500000000,8,64,1000000};
static void compare(double *a,size_t m,size_t k,size_t w,int expected){
 size_t bytes=m*m*sizeof(double);double *ref=malloc(bytes),*res=malloc(bytes);CHECK(ref&&res);memcpy(ref,a,bytes);CHECK(p66_set_mode(P66_REFERENCE));long double peak=0;int ok=p66_select_kernel()(ref,m,k,w,&peak);CHECK(ok==expected);
 for(int mode=P66_AVX2;mode<=P66_AVX512;mode++){if(!p66_available(mode))continue;CHECK(p66_set_mode(mode)&&p66_effective_mode()==(p66_mode)mode);memcpy(res,a,bytes);long double p=0;CHECK(p66_select_kernel()(res,m,k,w,&p)==ok);
  if(ok){for(size_t i=0;i<m;i++)for(size_t j=0;j<m;j++){if(i<k+w||j<k+w){CHECK(!memcmp(res+i*m+j,a+i*m+j,sizeof(double)));continue;}CHECK(res[i*m+j]==res[j*m+i]);long double den=fabsl(ref[i*m+j]);long double d=fabsl((long double)res[i*m+j]-ref[i*m+j]);long double e=den?d/den:d;CHECK(e<3e-15L);if(e>max_kernel)max_kernel=e;}CHECK(fabsl(p-peak)<=3e-15L*peak);}
  tests++;
 }free(ref);free(res);
}
static void kernels(void){for(unsigned t=0;t<420;t++){size_t m=2+rnd()%96,w=1+t%2;if(m<=w)m=w+1;size_t k=rnd()%(m-w);double scale=t%3==0?1e-100:t%3==1?1e100:1;double *a=calloc(m*m,sizeof *a);CHECK(a);for(size_t i=k+w;i<m;i++){for(size_t j=k+w;j<=i;j++)a[i*m+j]=a[j*m+i]=v()*scale;for(size_t l=0;l<w;l++){a[(k+l)*m+i]=v();a[i*m+k+l]=v()*scale;}}compare(a,m,k,w,1);free(a);}
 /* rank-2 intermediate exceeds double range, final finite; true overflow;
    exact cancellation, subnormal outputs; vector-length tails. */
 for(unsigned t=0;t<4;t++){size_t m=19,w=2;double a[361]={0};for(size_t i=w;i<m;i++){for(size_t j=w;j<=i;j++)a[i*m+j]=a[j*m+i]=t==2?1: t==3?DBL_TRUE_MIN:DBL_MAX; a[i]=t==2?1:t==3?0:-1;a[i*m]=t==2?1:t==3?0:DBL_MAX;a[m+i]=t==0?1:0;a[i*m+1]=t==0?DBL_MAX:0;}compare(a,m,0,w,t!=1);}
 CHECK(!p66_set_mode((p66_mode)99));tests++;CHECK(p66_set_mode(P66_AUTO));CHECK(fesetround(FE_DOWNWARD)==0&&p66_effective_mode()==P66_REFERENCE);CHECK(fesetround(FE_TONEAREST)==0);tests++;
}
static void solve(p56_context *ctx,uint32_t n,const size_t *cp,const uint32_t *ri,const double *a){double *x=malloc(n*sizeof *x),*b=calloc(n,sizeof *b),*y=malloc(n*sizeof *y),*w=malloc(n*sizeof *w);long double *r=calloc(n,sizeof *r),*d=calloc(n,sizeof *d);CHECK(x&&b&&y&&w&&r&&d);for(uint32_t i=0;i<n;i++)x[i]=((int)(i%19)-9)/11.;for(uint32_t j=0;j<n;j++)for(size_t z=cp[j];z<cp[j+1];z++){uint32_t i=ri[z];b[i]+=a[z]*x[j];if(i!=j)b[j]+=a[z]*x[i];}CHECK(!p56_apply(ctx,b,y,w,n));for(uint32_t i=0;i<n;i++){double e=fabs(x[i]-y[i]);CHECK(e<1e-9);if(e>max_x)max_x=e;r[i]=-b[i];d[i]=fabsl(b[i]);}for(uint32_t j=0;j<n;j++)for(size_t z=cp[j];z<cp[j+1];z++){uint32_t i=ri[z];long double q=(long double)a[z]*y[j];r[i]+=q;d[i]+=fabsl(q);if(i!=j){q=(long double)a[z]*y[i];r[j]+=q;d[j]+=fabsl(q);}}for(uint32_t i=0;i<n;i++){long double e=d[i]?fabsl(r[i])/d[i]:0;CHECK(e<1e-12L);if(e>max_back)max_back=e;}free(x);free(b);free(y);free(w);free(r);free(d);}
static void factors(void){for(unsigned t=0;t<160;t++){uint32_t n=2*(5+rnd()%27),ri[4096];size_t cp[65],z=0;double a[4096],scale=t%3==0?1e-100:t%3==1?1e100:1;int pairs=t%2;for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=pairs?0:scale*(j%2?1:-1)*2*n;for(uint32_t i=j+1;i<n;i++){if(!pairs&&rnd()%100>=t%100)continue;ri[z]=i;a[z++]=scale*(pairs?(i==j+1&&j%2==0?4:v()*.001):v());}}cp[n]=z;
 for(int mode=P66_REFERENCE;mode<=P66_AVX512;mode++){if(!p66_available(mode))continue;CHECK(p66_set_mode(mode));p56_context *ctx=NULL;p64_stats st={0};CHECK(!p65_ordered_factorize(&ctx,n,z,cp,ri,a,&ol,&sl,&nl,t%4,1+t%16,&st));if(pairs)CHECK(st.numeric.pairs>0);solve(ctx,n,cp,ri,a);p56_destroy(ctx);tests++;}
 }}
static void bench(uint32_t n,uint32_t side,uint32_t block){size_t *cp=malloc(((size_t)n+1)*sizeof *cp);uint32_t *ri=malloc((size_t)n*(block?block:3)*sizeof *ri);double *a=malloc((size_t)n*(block?block:3)*sizeof *a);CHECK(cp&&ri&&a);size_t z=0;for(uint32_t j=0;j<n;j++){cp[j]=z;ri[z]=j;a[z++]=block?2*block:5;if(block){uint32_t end=(j/block+1)*block;if(end>n)end=n;for(uint32_t i=j+1;i<end;i++){ri[z]=i;a[z++]=-.25;}}else{if(j+1<n&&(!side||j%side+1<side)){ri[z]=j+1;a[z++]=-1;}if(side&&j+side<n){ri[z]=j+side;a[z++]=-1;}}}cp[n]=z;
 for(int mode=P66_REFERENCE;mode<=P66_AVX512;mode++){if(!p66_available(mode))continue;CHECK(p66_set_mode(mode));double tm[3];p64_stats st={0};for(int rep=0;rep<3;rep++){p56_context *ctx=NULL;clock_t start=clock();CHECK(!p65_ordered_factorize(&ctx,n,z,cp,ri,a,&ol,&sl,&nl,side?P64_ND:P64_NATURAL,16,&st));tm[rep]=(double)(clock()-start)/CLOCKS_PER_SEC;solve(ctx,n,cp,ri,a);p56_destroy(ctx);}double med=tm[0]+tm[1]+tm[2]-fmin(tm[0],fmin(tm[1],tm[2]))-fmax(tm[0],fmax(tm[1],tm[2]));printf("bench n=%u side=%u block=%u mode=%s cpu_median3=%.6f s fronts=%zu max_front=%zu slots=%zu\n",n,side,block,p66_mode_name(mode),med,st.numeric.fronts,st.numeric.peak_front,st.numeric.factor_slots);tests++;}free(cp);free(ri);free(a);}
int main(int argc,char **argv){(void)argv;CHECK(p66_set_mode(P66_AUTO));printf("AUTO=%s AVX2=%d AVX512=%d\n",p66_mode_name(p66_effective_mode()),!!p66_available(P66_AVX2),!!p66_available(P66_AVX512));if(argc==1){kernels();factors();}else{bench(1600,40,0);bench(4096,64,0);bench(4096,0,64);bench(152592,0,0);}printf("SIMD %u/%u PASS max_kernel_rel=%.3Le max_x=%.3e max_componentwise_backward=%.3Le\n",tests,tests,max_kernel,max_x,max_back);return 0;}
