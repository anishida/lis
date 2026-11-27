#include "schur_simd.h"
#include <math.h>
#include <float.h>
#include <fenv.h>
#if (defined(__x86_64__) || defined(__i386__)) && (defined(__GNUC__) || defined(__clang__)) && !defined(P66_DISABLE_SIMD)
#define P66_X86 1
#include <immintrin.h>
#else
#define P66_X86 0
#endif
static _Thread_local p66_mode requested=P66_AUTO;
int p66_available(p66_mode mode){
 if(mode==P66_AUTO||mode==P66_REFERENCE)return 1;
#if P66_X86
 if(mode==P66_AVX2)return __builtin_cpu_supports("avx2")&&__builtin_cpu_supports("fma");
 if(mode==P66_AVX512)return __builtin_cpu_supports("avx512f")&&__builtin_cpu_supports("fma");
#endif
 return 0;
}
int p66_set_mode(p66_mode mode){if(!p66_available(mode))return 0;requested=mode;return 1;}
p66_mode p66_effective_mode(void){
 if(fegetround()!=FE_TONEAREST)return P66_REFERENCE;
#if P66_X86
 if(_mm_getcsr()&((1u<<15)|(1u<<6)))return P66_REFERENCE;
#endif
 if(requested!=P66_AUTO)return requested;
 /* AVX2 is the conservative default: the measured triangular tails
    and gather costs favor it on the development CPU. */
 if(p66_available(P66_AVX2))return P66_AVX2;
 if(p66_available(P66_AVX512))return P66_AVX512;
 return P66_REFERENCE;
}
const char *p66_mode_name(p66_mode mode){switch(mode){case P66_AUTO:return "AUTO";case P66_REFERENCE:return "REFERENCE";case P66_AVX2:return "AVX2";case P66_AVX512:return "AVX512";default:return "INVALID";}}
static inline int precise(double *a,size_t m,size_t k,size_t w,size_t i,size_t j,double *out){
 long double z=fmal(-(long double)a[k*m+i],a[j*m+k],a[i*m+j]);
 if(w==2)z=fmal(-(long double)a[(k+1)*m+i],a[j*m+k+1],z);
 if(!isfinite(z)||fabsl(z)>DBL_MAX)return 0;
 *out=(double)z;return 1;
}
static inline void commit(double *a,size_t m,size_t i,size_t j,double z,long double *peak){a[i*m+j]=a[j*m+i]=z;if(fabsl(z)>*peak)*peak=fabsl(z);}
static int reference(double *a,size_t m,size_t k,size_t w,long double *peak){
 for(size_t i=k+w;i<m;i++)for(size_t j=k+w;j<=i;j++){double z;if(!precise(a,m,k,w,i,j,&z))return 0;commit(a,m,i,j,z,peak);}
 return 1;
}
#if P66_X86
/* Rank-2 cancellation may amplify intermediate rounding. Recompute sensitive
   lanes and range extremes using original long-double precision. */
static inline int finish(double *a,size_t m,size_t k,size_t w,size_t i,size_t j,double t,double z,long double *peak){
 double az=fabs(z),at=fabs(t);
 if(!isfinite(z)||az>=DBL_MAX||az<DBL_MIN||(w==2&&(!isfinite(t)||at>=DBL_MAX||at<DBL_MIN||az<at*0.125))){
  if(!precise(a,m,k,w,i,j,&z))return 0;
 }
 commit(a,m,i,j,z,peak);return 1;
}
__attribute__((target("avx2,fma"))) static int avx2(double *a,size_t m,size_t k,size_t w,long double *peak){
 __m256i stride=_mm256_set_epi64x((long long)(3*m),(long long)(2*m),(long long)m,0);
 for(size_t i=k+w;i<m;i++){
  size_t j=k+w;__m256d x=_mm256_set1_pd(a[k*m+i]),y=_mm256_set1_pd(w==2?a[(k+1)*m+i]:0);
  for(;j<=i&&i-j+1>=4;j+=4){double t[4],z[4];__m256d c=_mm256_loadu_pd(a+i*m+j),u=_mm256_i64gather_pd(a+j*m+k,stride,8);
   __m256d v=_mm256_fnmadd_pd(x,u,c),r=v;
   if(w==2)r=_mm256_fnmadd_pd(y,_mm256_i64gather_pd(a+j*m+k+1,stride,8),v);
   _mm256_storeu_pd(t,v);_mm256_storeu_pd(z,r);
   for(size_t l=0;l<4;l++)if(!finish(a,m,k,w,i,j+l,t[l],z[l],peak))return 0;
  }
  for(;j<=i;j++){double z;if(!precise(a,m,k,w,i,j,&z))return 0;commit(a,m,i,j,z,peak);}
 }
 return 1;
}
__attribute__((target("avx512f,fma"))) static int avx512(double *a,size_t m,size_t k,size_t w,long double *peak){
 __m512i stride=_mm512_set_epi64((long long)(7*m),(long long)(6*m),(long long)(5*m),(long long)(4*m),(long long)(3*m),(long long)(2*m),(long long)m,0);
 for(size_t i=k+w;i<m;i++){
  size_t j=k+w;__m512d x=_mm512_set1_pd(a[k*m+i]),y=_mm512_set1_pd(w==2?a[(k+1)*m+i]:0);
  for(;j<=i&&i-j+1>=8;j+=8){double t[8],z[8];__m512d c=_mm512_loadu_pd(a+i*m+j),u=_mm512_i64gather_pd(stride,a+j*m+k,8);
   __m512d v=_mm512_fnmadd_pd(x,u,c),r=v;
   if(w==2)r=_mm512_fnmadd_pd(y,_mm512_i64gather_pd(stride,a+j*m+k+1,8),v);
   _mm512_storeu_pd(t,v);_mm512_storeu_pd(z,r);
   for(size_t l=0;l<8;l++)if(!finish(a,m,k,w,i,j+l,t[l],z[l],peak))return 0;
  }
  for(;j<=i;j++){double z;if(!precise(a,m,k,w,i,j,&z))return 0;commit(a,m,i,j,z,peak);}
 }
 return 1;
}
#endif
p66_kernel p66_select_kernel(void){
#if P66_X86
 switch(p66_effective_mode()){case P66_AVX512:return avx512;case P66_AVX2:return avx2;default:break;}
#endif
 return reference;
}
