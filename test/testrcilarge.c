/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
/* Native LIS bridge exercise beyond the 4096-element RCI OpenMP threshold. */
#include "lis_rci_bridge.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#define N 8192
#define CHECK(c) do{if(!(c)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c);rc=1;goto cleanup;}}while(0)
static void product(const double*x,double*y){for(size_t i=0;i<N;i++)y[i]=(4.0+(i%17)/17.0)*x[i]-(i?x[i-1]:0)-(i+1<N?x[i+1]:0);}
/* Drive all three RCI interfaces through the same native LIS matrix.
 * Build the RHS and validate the result independently of LIS matvec. */
static int cross_test(void)
{
 LIS_MATRIX A=NULL;
 LIS_VECTOR u=NULL,v=NULL;
 lis_rci_context *c=NULL;
 double *mem=NULL;
 int rc=0;
 CHECK(LIS_RCI_CHECK_CLIENT_ABI()==LIS_RCI_DONE);
 mem=calloc(6*N,sizeof(double));CHECK(mem);
 double *b=mem,*x0=b+N,*exact=x0+N,*ax=exact+N,*buf=ax+N,*reference=buf+N;
 for(int method=1;method<=4;method++){
  double lower=method==LIS_RCI_CG?-1.0:-0.6;
  double upper=method==LIS_RCI_CG?-1.0:-1.3;
  CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));
  CHECK(!lis_matrix_set_size(A,0,N));
  for(LIS_INT i=0;i<N;i++){
   CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0+(i%17)/17.0,A));
   if(i)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,lower,A));
   if(i+1<N)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,upper,A));
  }
  CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
  CHECK(!lis_vector_duplicate(A,&u));CHECK(!lis_vector_duplicate(A,&v));
  lis_rci_options o={N,200,12,1e-10,0,(lis_rci_method)method};
  c=lis_rci_create(&o);CHECK(c);
  for(int rep=0;rep<2;rep++){
   for(size_t i=0;i<N;i++){
    exact[i]=sin(.07*i+rep)+.3*cos(.19*(i+3));
    x0[i]=rep?.125:0;
   }
   for(size_t i=0;i<N;i++)
    b[i]=(4.0+(i%17)/17.0)*exact[i]+(i?lower*exact[i-1]:0)+(i+1<N?upper*exact[i+1]:0);
   for(int mode=0;mode<3;mode++){
    lis_rci_request q;lis_rci_request_ex ex;lis_rci_request_full f;
    int s=mode==0?lis_rci_start(c,b,x0,&q):mode==1?lis_rci_start_ex(c,b,x0,&ex):lis_rci_start_full(c,b,x0,&f);
    size_t count=0,precons=0;
    while(s>0){
     CHECK(++count<100000&&lis_rci_get_pending(c,&f)==s);
     const lis_rci_request_ex *p=&f.request;
     CHECK(p->n==N&&p->in&&p->token);
     long double scalar=0;
     if(s==LIS_RCI_MATVEC){
      for(size_t i=0;i<N;i++)buf[i]=p->in[i];
      CHECK(!lis_vector_set_values2(LIS_INS_VALUE,0,N,buf,u));
      CHECK(!lis_matvec(A,u,v));CHECK(!lis_vector_get_values(v,0,N,p->out));
     }else if(s==LIS_RCI_PRECON){
      precons++;
      for(size_t i=0;i<N;i++)
       p->out[i]=p->in[i]/((4.0+(i%17)/17.0)*(method==LIS_RCI_FGMRES?1+.025*((i+precons)%7):1));
     }else if(s==LIS_RCI_DOT||s==LIS_RCI_NORM){
      CHECK(mode!=0);
      for(size_t i=0;i<N;i++)scalar+=(long double)p->in[i]*(p->in2?p->in2[i]:p->in[i]);
      if(s==LIS_RCI_NORM)scalar=sqrtl(scalar);
     }else{
      CHECK(mode==2&&(s==LIS_RCI_AXPBY||s==LIS_RCI_DIVIDE));
      for(size_t i=0;i<N;i++){
       if(s==LIS_RCI_AXPBY)p->out[i]=(double)(f.alpha*p->in[i]+(p->in2?f.beta*p->in2[i]:0));
       else p->out[i]=(double)((long double)p->in[i]/f.alpha);
      }
     }
     uint64_t token=p->token;
     if(mode==0)s=lis_rci_advance(c,0,&q);
     else if(mode==1)s=lis_rci_advance_ex(c,token,0,scalar,&ex);
     else s=lis_rci_advance_full(c,token,0,scalar,&f);
    }
    const double *x=lis_rci_solution(c);lis_rci_info info;lis_rci_stats stats;
    CHECK(!lis_rci_get_info(c,&info)&&!lis_rci_get_stats(c,&stats));
    long double rn=0,bn=0,en=0,xn=0,dn=0;
    for(size_t i=0;i<N;i++){
     ax[i]=(4.0+(i%17)/17.0)*x[i]+(i?lower*x[i-1]:0)+(i+1<N?upper*x[i+1]:0);
     long double r=(long double)b[i]-ax[i],e=(long double)x[i]-exact[i];
     rn+=r*r;bn+=(long double)b[i]*b[i];en+=e*e;xn+=(long double)exact[i]*exact[i];
     if(mode==0)reference[i]=x[i];
     else {long double d=(long double)x[i]-reference[i];dn+=d*d;}
    }
    double rr=(double)sqrtl(rn/bn),er=(double)sqrtl(en/xn),dr=(double)sqrtl(dn/xn);
    int ok=s==LIS_RCI_DONE&&info.status==s&&rr<=1.05e-10&&er<=1e-8&&dr<=1e-8;
    printf("RCI CROSS method=%d mode=%d rhs=%d iter=%zu residual=%.9e error=%.9e delta=%.9e %s\n",
      method,mode,rep,info.iterations,rr,er,dr,ok?"PASS":"FAIL");
    CHECK(ok&&stats.issued[0]==info.matvecs&&stats.issued[1]==info.precons);
    for(size_t i=0;i<6;i++)CHECK(stats.issued[i]==stats.completed[i]);
    CHECK(mode==2?stats.issued[4]>0:stats.issued[4]==0);
   }
  }
  lis_rci_destroy(c);c=NULL;
  lis_vector_destroy(v);v=NULL;lis_vector_destroy(u);u=NULL;
  lis_matrix_destroy(A);A=NULL;
 }
cleanup:
 lis_rci_destroy(c);if(v)lis_vector_destroy(v);if(u)lis_vector_destroy(u);
 if(A)lis_matrix_destroy(A);free(mem);return rc;
}

int main(int argc,char**argv){
 LIS_MATRIX A=NULL;lis_rci_bridge*p=NULL;int rc=0;double*mem=NULL;
 if(lis_initialize(&argc,&argv))return 1;
#ifdef _OPENMP
 int threads=0;
#pragma omp parallel reduction(+:threads)
 {threads++;}
 printf("OpenMP actual team=%d\n",threads);
 const char*expect=getenv("OMP_NUM_THREADS");if(expect)CHECK(threads==atoi(expect));
#endif
 CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,N));
 for(LIS_INT i=0;i<N;i++){CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0+(i%17)/17.0,A));if(i)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1,A));if(i+1<N)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1,A));}
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
 mem=calloc(4*N,sizeof(double));CHECK(mem!=NULL);double*x=mem,*b=x+N,*exact=b+N,*ax=exact+N;
 for(int m=1;m<=4;m++){
  lis_rci_options o={N,200,30,1e-10,0,(lis_rci_method)m};CHECK(!lis_rci_bridge_create(A,&o,1,&p));
  for(int rep=0;rep<2;rep++){
   for(size_t i=0;i<N;i++){exact[i]=sin(.07*i+rep)+.3*cos(.19*(i+3));x[i]=rep?.125:0;}
   product(exact,b);lis_rci_info info;int s=lis_rci_bridge_solve(p,b,x,&info);product(x,ax);
   long double r=0,bn=0,e=0,en=0;for(size_t i=0;i<N;i++){long double d=(long double)b[i]-ax[i],t=(long double)x[i]-exact[i];r+=d*d;bn+=(long double)b[i]*b[i];e+=t*t;en+=(long double)exact[i]*exact[i];}
   double rr=(double)sqrtl(r/bn),er=(double)sqrtl(e/en);int ok=s==LIS_RCI_DONE&&rr<=1.05e-10&&er<=1e-8&&isfinite(rr)&&isfinite(er);
   printf("LARGE method=%d rhs=%d iter=%zu residual=%.9e error=%.9e %s\n",m,rep,info.iterations,rr,er,ok?"PASS":"FAIL");CHECK(ok);
  }
  lis_rci_bridge_destroy(p);p=NULL;
 }
 CHECK(!cross_test());
cleanup:
 free(mem);lis_rci_bridge_destroy(p);if(A)lis_matrix_destroy(A);lis_finalize();return rc;
}
