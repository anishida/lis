#include "schur_simd.h"
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"FAIL %d %s\n",__LINE__,#x);exit(1);}}while(0)
static unsigned tests;
static double v(size_t i,size_t j){return ((int)((i*19+j*13)%37)-18)*0.03125;}
static void run(size_t n,size_t w,unsigned threads,int kind){
 size_t k=3,m=k+w+n,bytes=m*m*sizeof(double);
 double *a=calloc(m*m,sizeof(double)),*b=malloc(bytes);CHECK(a&&b);
 for(size_t i=k+w;i<m;i++){
  for(size_t j=k+w;j<=i;j++)a[i*m+j]=v(i,j);
  for(size_t j=i+1;j<m;j++)a[i*m+j]=NAN;
  for(size_t l=0;l<w;l++){a[(k+l)*m+i]=v(i,l);a[i*m+k+l]=v(l,i);}
 }
 if(kind==1){a[k*m+k+w]=DBL_MAX;a[(k+w)*m+k]=DBL_MAX;} /* range failure */
 if(kind==2){a[k*m+k+w]=-0.0;a[(k+w)*m+k]=DBL_TRUE_MIN;}
 memcpy(b,a,bytes);CHECK(p66_set_threads(1));long double pa=1,pb=1;int ra=p66_select_kernel_lower()(a,m,k,w,&pa);
 CHECK(p66_set_threads(threads));int rb=p66_select_kernel_lower()(b,m,k,w,&pb);CHECK(ra==rb);
 if(ra){CHECK(pa==pb);CHECK(!memcmp(a,b,bytes));}
 else CHECK(kind==1);
 free(a);free(b);tests++;
}
int main(void){
 CHECK(!p66_set_threads(0)&&!p66_set_threads(65));if(!p66_available(P66_AVX2)){puts("OpenMP AVX2 unavailable: SKIP");return 77;}CHECK(p66_set_mode(P66_AVX2));
 for(size_t w=1;w<=2;w++)for(size_t n=381;n<=389;n++)for(unsigned t=2;t<=4;t+=2)for(int c=0;c<=2;c++)run(n,w,t,c);
 for(size_t w=1;w<=2;w++)for(unsigned t=2;t<=4;t+=2)for(int c=0;c<=2;c++)run(769,w,t,c);
 /* Nested calls must stay serial; same result and calling thread mode. */
 #pragma omp parallel num_threads(2)
 {
  #pragma omp single
  {CHECK(p66_set_mode(P66_AVX2));run(400,2,4,0);}
 }
 CHECK(p66_set_threads(1));printf("Schur OpenMP %u/%u PASS\n",tests,tests);return 0;
}
