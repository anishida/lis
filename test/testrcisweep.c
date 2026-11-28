/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include "lis_rci.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#ifdef _OPENMP
#include <omp.h>
#endif
static unsigned cases,failures;
#define CHECK(e) do{cases++;if(!(e)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#e);failures++;}}while(0)
static void apply(size_t n,const double*x,double*y,int nonsym,int mode){
#ifdef _OPENMP
#pragma omp parallel for if(n>=4096)
#endif
 for(size_t i=0;i<n;i++){
  double d=4.+(double)(i%17)/17.;
  if(mode==1)y[i]=x[i];else if(mode==2)y[i]=0;else if(mode==3)y[i]=-x[i];
  else y[i]=d*x[i]-(i?x[i-1]:0)-(i+1<n?(nonsym?0.35:1.)*x[i+1]:0);
 }}
static int run(lis_rci_context*c,const double*b,const double*x0,int nonsym,int mode,int prec,int fault){lis_rci_request q;int st=lis_rci_start(c,b,x0,&q);size_t steps=0;while(st>0&&steps++<200000){if(fault==1)return lis_rci_advance(c,-77,&q);if(fault==2)return lis_rci_advance(c,0,&q);if(fault==3){for(size_t i=0;i<q.n;i++)q.out[i]=INFINITY;return lis_rci_advance(c,0,&q);}if(st==LIS_RCI_MATVEC)apply(q.n,q.in,q.out,nonsym,mode);else if(st==LIS_RCI_PRECON){for(size_t i=0;i<q.n;i++)q.out[i]=prec?q.in[i]/((4.+(double)(i%17)/17.)*(prec==2?(1.+.05*((steps+i)%7)):1.)):q.in[i];}else return -99;st=lis_rci_advance(c,0,&q);}return st;}
static void solve_case(lis_rci_method method,size_t n,size_t restart,int prec,int mode){
 lis_rci_options o={n,500,restart,1e-10,1e-13,method};lis_rci_context*c=lis_rci_create(&o);CHECK(c!=NULL);if(!c)return;
 double *buf=calloc(4*n,sizeof(double));if(!buf)exit(2);double *exact=buf,*b=buf+n,*x0=buf+2*n,*ax=buf+3*n;int nonsym=method!=LIS_RCI_CG;
 for(int rep=0;rep<3;rep++){
 for(size_t i=0;i<n;i++){exact[i]=sin((i+1)*.37+rep)+.3*cos((i+3)*.19);x0[i]=rep==0?0:(rep==1?.125:exact[i]);}
 apply(n,exact,b,nonsym,mode);int st=run(c,b,x0,nonsym,mode,prec,0);lis_rci_info info;lis_rci_get_info(c,&info);
 if(st!=LIS_RCI_DONE)fprintf(stderr,"method=%d n=%zu restart=%zu prec=%d rep=%d status=%d iter=%zu residual=%.17g\n",method,n,restart,prec,rep,st,info.iterations,info.residual);
 CHECK(st==LIS_RCI_DONE);const double*x=lis_rci_solution(c);apply(n,x,ax,nonsym,mode);long double r=0,bn=0,err=0,xn=0;
 for(size_t i=0;i<n;i++){long double d=(long double)b[i]-ax[i],e=(long double)x[i]-exact[i];r+=d*d;bn+=(long double)b[i]*b[i];err+=e*e;xn+=(long double)exact[i]*exact[i];}
 CHECK(sqrtl(r)<=fmaxl(o.atol,o.rtol*sqrtl(bn)));CHECK(sqrtl(err)<=1e-8L*fmaxl(1,sqrtl(xn)));CHECK(info.matvecs>=1);CHECK(info.iterations<=o.maxiter);CHECK(info.status==st);if(rep==2)CHECK(info.iterations==0);
 }
 free(buf);lis_rci_destroy(c);
}
static void edge(lis_rci_method m){lis_rci_options o={4,1,4,1e-13,1e-15,m};lis_rci_context*c=lis_rci_create(&o);double b[]={1,2,3,4},z[]={0,0,0,0},bad[]={NAN,0,0,0};lis_rci_request q;lis_rci_info info;CHECK(c!=NULL);if(!c)return;
 CHECK(lis_rci_advance(c,0,&q)==LIS_RCI_STATE);CHECK(run(c,b,NULL,0,1,0,0)==LIS_RCI_DONE);CHECK(run(c,z,NULL,0,1,0,0)==LIS_RCI_DONE);lis_rci_get_info(c,&info);CHECK(info.iterations==0);
 CHECK(run(c,b,NULL,0,2,0,0)==LIS_RCI_BREAKDOWN);CHECK(run(c,b,NULL,0,0,0,0)==LIS_RCI_MAXITER);CHECK(run(c,b,NULL,0,1,0,1)==LIS_RCI_CALLBACK);lis_rci_get_info(c,&info);CHECK(info.callback_code==-77);
 CHECK(run(c,b,NULL,0,1,0,2)==LIS_RCI_NONFINITE);CHECK(run(c,b,NULL,0,1,0,3)==LIS_RCI_NONFINITE);CHECK(lis_rci_start(c,bad,NULL,&q)==LIS_RCI_NONFINITE);
 CHECK(lis_rci_start(c,b,NULL,&q)==LIS_RCI_MATVEC);lis_rci_request saved=q;CHECK(lis_rci_start(c,b,NULL,&q)==LIS_RCI_STATE);CHECK(saved.in==q.in&&saved.out==q.out);CHECK(lis_rci_advance(c,9,&q)==LIS_RCI_CALLBACK);CHECK(run(c,b,NULL,0,1,0,0)==LIS_RCI_DONE);
 if(m==LIS_RCI_CG)CHECK(run(c,b,NULL,0,3,0,0)==LIS_RCI_BREAKDOWN);
 lis_rci_destroy(c);o.n=SIZE_MAX;CHECK(lis_rci_create(&o)==NULL);o.n=4;o.rtol=NAN;CHECK(lis_rci_create(&o)==NULL);o.rtol=1e-8;o.maxiter=0;CHECK(lis_rci_create(&o)==NULL);
}
int main(void){
#ifdef _OPENMP
 printf("OpenMP threads max=%d\n",omp_get_max_threads());
#else
 puts("Serial build");
#endif
 for(int m=1;m<=4;m++){for(size_t n=1;n<=127;n=2*n+1){for(int p=0;p<=(m==4?2:1);p++)solve_case((lis_rci_method)m,n,7,p,0);}solve_case((lis_rci_method)m,31,1,1,0);solve_case((lis_rci_method)m,17,16,0,1);solve_case((lis_rci_method)m,4096,12,1,0);if(m==4)solve_case((lis_rci_method)m,4096,12,2,0);edge((lis_rci_method)m);printf("method=%d cumulative_checks=%u failures=%u\n",m,cases,failures);}
 printf("RESULT checks=%u failures=%u\n",cases,failures);return failures?1:0;}
