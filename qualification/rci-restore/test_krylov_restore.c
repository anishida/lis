/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#include "lis_rci_krylov_restore.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#ifdef LIS_RCI_KRYLOV_NATIVE_LIBRARY
#define REFERENCE lis_rci_krylov_restore_reference64
#else
#define REFERENCE lis_krylov_symmetric_ex64
#endif
#define N 24
#define NEV 3
#define C(e) do{if(!(e)){fprintf(stderr,"KRYLOV FAIL %d: %s\n",__LINE__,#e);return 1;}}while(0)
typedef struct {uint64_t hash,count;int mode,psd;double sigma;uint64_t fail;} trace;
static void mix(trace *t,const void *p,size_t n)
{const unsigned char *a=p;for(size_t i=0;i<n;i++){t->hash^=a[i];t->hash*=UINT64_C(1099511628211);}}
static int apply(trace *t,uint32_t op,int64_t n,const double *x,double *y)
{
 t->count++;mix(t,&op,sizeof op);mix(t,x,(size_t)n*sizeof(double));
 if(t->count==t->fail)return 71;
 for(int64_t i=0;i<n;i++){
  double B=t->psd&&i>=n-4?0:1+.02*i;
  double A=(1+.11*i)*(1+.02*i);
  y[i]=op==LIS_RCI_OP_A?A*x[i]:op==LIS_RCI_OP_B?B*x[i]:x[i]/(A-t->sigma*B);
 }
 mix(t,y,(size_t)n*sizeof(double));return 0;
}
static int a(void *u,int64_t n,const double *x,double *y){return apply(u,LIS_RCI_OP_A,n,x,y);}
static int b(void *u,int64_t n,const double *x,double *y){return apply(u,LIS_RCI_OP_B,n,x,y);}
static int solve(void *u,int64_t n,const double *x,double *y){return apply(u,LIS_RCI_OP_SHIFTED_SOLVE,n,x,y);}
static int compare(int mode,int psd,int locking,int adaptive,const char *which,int reuse_flag,uint64_t fail)
{
 trace t1={.hash=UINT64_C(1469598103934665603),.mode=mode,.psd=psd,.sigma=.7,.fail=fail},t2=t1,t3=t1;
 lis_krylov_problem64 p={N,mode,t1.sigma,&t1,a,b,solve},p2=p,p3=p;p2.ctx=&t2;p3.ctx=&t3;
 lis_krylov_options o;lis_krylov_options_default(&o);
 o.nev=NEV;o.ncv=9;o.maxiter=120;o.tol=1e-9;o.which=which;o.hard_locking=locking;
 o.allow_semidefinite_metric=psd;o.inner_tol=mode?1e-12:0;
 lis_krylov_ext_options ext;lis_krylov_ext_options_default(&ext);
 ext.adaptive_ncv=adaptive;ext.ncv_max=adaptive?16:0;ext.restart_directions=adaptive?2:1;
 lis_krylov_workspace *w1=NULL,*w2=NULL,*w3=NULL;
 if(reuse_flag){C(!lis_krylov_workspace_create64(&p,&o,&ext,&w1));C(!lis_krylov_workspace_create64(&p2,&o,&ext,&w2));C(!lis_krylov_workspace_create64(&p3,&o,&ext,&w3));}
 double v1[NEV]={0},v2[NEV]={0},x1[NEV*N]={0},x2[NEV*N]={0},r1[NEV]={0},r2[NEV]={0};
 lis_krylov_stats s1,s2;lis_krylov_ext_stats e1,e2;lis_krylov_ext_stats_default(&e1);lis_krylov_ext_stats_default(&e2);
 int rc1=REFERENCE(&p,&o,&ext,w1,v1,x1,r1,&s1,&e1);
 lis_rci_krylov_restore *K=lis_rci_krylov_restore_create();C(K);
 lis_rci_krylov_restore_request q;
 int st=lis_rci_krylov_restore_start(K,&p2,&o,&ext,w2,NULL,v2,x2,r2,&s2,&e2,&q);
 uint64_t last=0;
 while(st>0){
  C(q.token>last&&q.n==N&&q.operation==(uint32_t)st&&q.sigma==p.sigma);last=q.token;
  C(lis_rci_krylov_restore_advance(K,q.token-1,0,&q)==LIS_RCI_PROTOCOL&&q.token==last);
  int err=apply(&t2,q.operation,q.n,q.in,q.out);
  st=lis_rci_krylov_restore_advance(K,q.token,err,&q);
 }
 int rc2=-1;lis_krylov_budget_stats counts;
 C(lis_rci_krylov_restore_status(K,&rc2,&counts)==st);
 if(rc1!=rc2||t1.hash!=t2.hash||t1.count!=t2.count)
  fprintf(stderr,"DIFF mode=%d psd=%d locking=%d adaptive=%d which=%s rc=%d/%d count=%llu/%llu hash=%llx/%llx\n",mode,psd,locking,adaptive,which,rc1,rc2,(unsigned long long)t1.count,(unsigned long long)t2.count,(unsigned long long)t1.hash,(unsigned long long)t2.hash);
 C(rc1==rc2&&t1.hash==t2.hash&&t1.count==t2.count&&counts.total==t2.count);
 C(s1.nconv==s2.nconv&&s1.nreturned==s2.nreturned&&s1.iterations==s2.iterations&&s1.restarts==s2.restarts&&s1.op_count==s2.op_count&&s1.metric_count==s2.metric_count&&s1.solve_count==s2.solve_count&&s1.max_residual==s2.max_residual);
 C(e1.ncv_initial==e2.ncv_initial&&e1.ncv_final==e2.ncv_final&&e1.ncv_growths==e2.ncv_growths&&e1.basis_vectors==e2.basis_vectors&&e1.projected_solves==e2.projected_solves&&e1.extra_seed_attempts==e2.extra_seed_attempts&&e1.extra_seed_accepted==e2.extra_seed_accepted);
 C(!memcmp(v1,v2,sizeof v1)&&!memcmp(x1,x2,sizeof x1)&&!memcmp(r1,r2,sizeof r1));
 double v3[NEV]={0},x3[NEV*N]={0},r3[NEV]={0};
 lis_krylov_stats s3;lis_krylov_ext_stats e3;lis_krylov_ext_stats_default(&e3);
 int rc3=lis_rci_krylov_restore_classic64(&p3,&o,&ext,w3,NULL,v3,x3,r3,&s3,&e3,&counts);
 C(rc1==rc3&&t1.hash==t3.hash&&t1.count==t3.count&&counts.total==t3.count);
 C(s1.nconv==s3.nconv&&s1.nreturned==s3.nreturned&&s1.iterations==s3.iterations&&s1.restarts==s3.restarts&&s1.op_count==s3.op_count&&s1.metric_count==s3.metric_count&&s1.solve_count==s3.solve_count&&s1.max_residual==s3.max_residual);
 C(e1.ncv_initial==e3.ncv_initial&&e1.ncv_final==e3.ncv_final&&e1.ncv_growths==e3.ncv_growths&&e1.basis_vectors==e3.basis_vectors&&e1.projected_solves==e3.projected_solves&&e1.extra_seed_attempts==e3.extra_seed_attempts&&e1.extra_seed_accepted==e3.extra_seed_accepted);
 C(!memcmp(v1,v3,sizeof v1)&&!memcmp(x1,x3,sizeof x1)&&!memcmp(r1,r3,sizeof r1));
 lis_rci_krylov_restore_destroy(K);lis_krylov_workspace_destroy(w1);lis_krylov_workspace_destroy(w2);lis_krylov_workspace_destroy(w3);return 0;
}
static int cancellation(void)
{
 trace t={.hash=1,.sigma=.7};lis_krylov_problem64 p={N,1,t.sigma,&t,a,b,solve};
 lis_krylov_options o;lis_krylov_options_default(&o);o.nev=NEV;o.ncv=9;o.inner_tol=1e-12;
 lis_krylov_workspace *w=NULL;C(!lis_krylov_workspace_create64(&p,&o,NULL,&w));
 lis_rci_krylov_restore *K=lis_rci_krylov_restore_create();C(K);
 double v[NEV],x[NEV*N],r[NEV];lis_krylov_stats s;
 lis_rci_krylov_restore_request q;uint64_t last=0;
 for(int stop=0;stop<20;stop++){
  int st=lis_rci_krylov_restore_start(K,&p,&o,NULL,w,NULL,v,x,r,&s,NULL,&q);
  C(st>0&&q.token>last);
  for(int i=0;i<stop;i++){C(st>0);int err=apply(&t,q.operation,q.n,q.in,q.out);st=lis_rci_krylov_restore_advance(K,q.token,err,&q);}
  C(st>0);last=q.token;
  C(lis_rci_krylov_restore_cancel(K,q.token-1)==LIS_RCI_PROTOCOL);
  C(lis_rci_krylov_restore_cancel(K,q.token)==LIS_RCI_CANCELLED&&s.nreturned==0);
 }
 lis_krylov_budget budget;lis_krylov_budget_default(&budget);budget.max_total=7;
 int st=lis_rci_krylov_restore_start(K,&p,&o,NULL,w,&budget,v,x,r,&s,NULL,&q);
 while(st>0){int err=apply(&t,q.operation,q.n,q.in,q.out);st=lis_rci_krylov_restore_advance(K,q.token,err,&q);}
 int native;lis_krylov_budget_stats counts;
 C(st==LIS_RCI_MAXITER&&lis_rci_krylov_restore_status(K,&native,&counts)==st);
 C(native==LIS_KRYLOV_BUDGET_EXHAUSTED&&counts.total==7&&counts.exhausted&&s.nreturned==0);
 st=lis_rci_krylov_restore_start(K,&p,&o,NULL,w,NULL,v,x,r,&s,NULL,&q);
 while(st>0){int err=apply(&t,q.operation,q.n,q.in,q.out);st=lis_rci_krylov_restore_advance(K,q.token,err,&q);}
 C(st==0||st==LIS_RCI_MAXITER);
 lis_rci_krylov_restore_destroy(K);lis_krylov_workspace_destroy(w);return 0;
}

/* Two live sessions interleave requests at different suspension depths.
 * Cancel one and ensure the other retains the original trace and results. */
static int interleaved(void)
{
 trace base={.hash=1,.sigma=.7},ta=base,tb=base;
 lis_krylov_problem64 pa={N,1,.7,&ta,a,b,solve},pb=pa,ref=pa;pb.ctx=&tb;ref.ctx=&base;
 lis_krylov_options o;lis_krylov_options_default(&o);o.nev=NEV;o.ncv=9;o.inner_tol=1e-12;o.maxiter=60;
 double va[NEV]={0},xa[NEV*N]={0},ra[NEV]={0},vb[NEV]={0},xb[NEV*N]={0},rb[NEV]={0},vr[NEV]={0},xr[NEV*N]={0},rr[NEV]={0};
 lis_krylov_stats sa,sb,sr;
 int r=REFERENCE(&ref,&o,NULL,NULL,vr,xr,rr,&sr,NULL);
 lis_rci_krylov_restore *ka=lis_rci_krylov_restore_create(),*kb=lis_rci_krylov_restore_create();C(ka&&kb);
 lis_rci_krylov_restore_request qa,qb;
 int aa=lis_rci_krylov_restore_start(ka,&pa,&o,NULL,NULL,NULL,va,xa,ra,&sa,NULL,&qa);
 int bb=lis_rci_krylov_restore_start(kb,&pb,&o,NULL,NULL,NULL,vb,xb,rb,&sb,NULL,&qb);
 for(int i=0;i<12;i++){
  C(aa>0&&bb>0);
  int e=apply(&tb,qb.operation,qb.n,qb.in,qb.out);bb=lis_rci_krylov_restore_advance(kb,qb.token,e,&qb);
  e=apply(&ta,qa.operation,qa.n,qa.in,qa.out);aa=lis_rci_krylov_restore_advance(ka,qa.token,e,&qa);
 }
 C(aa>0);C(lis_rci_krylov_restore_cancel(ka,qa.token)==LIS_RCI_CANCELLED);
 while(bb>0){int e=apply(&tb,qb.operation,qb.n,qb.in,qb.out);bb=lis_rci_krylov_restore_advance(kb,qb.token,e,&qb);}
 int native;lis_krylov_budget_stats counts;C(lis_rci_krylov_restore_status(kb,&native,&counts)==bb);
 C(native==r&&tb.hash==base.hash&&tb.count==base.count);
 C(!memcmp(vb,vr,sizeof vb)&&!memcmp(xb,xr,sizeof xb)&&!memcmp(rb,rr,sizeof rb));
 lis_rci_krylov_restore_destroy(ka);lis_rci_krylov_restore_destroy(kb);return 0;
}
/* n reaches the native parallel-kernel threshold. Compare in a fixed
 * thread configuration, not between different reduction partitions. */
static int large(void)
{
 const int64_t n=32768;const int nev=2;
 trace t1={.hash=1,.sigma=.7},t2=t1;
 lis_krylov_problem64 p={n,0,.7,&t1,a,b,solve},p2=p;p2.ctx=&t2;
 lis_krylov_options o;lis_krylov_options_default(&o);o.nev=nev;o.ncv=6;o.maxiter=5;
 double v1[2]={0},v2[2]={0},r1[2]={0},r2[2]={0};
 double *x1=calloc((size_t)n*nev,sizeof(double)),*x2=calloc((size_t)n*nev,sizeof(double));C(x1&&x2);
 lis_krylov_stats s1,s2;lis_krylov_budget_stats counts;
 int r=REFERENCE(&p,&o,NULL,NULL,v1,x1,r1,&s1,NULL);
 int r2code=lis_rci_krylov_restore_classic64(&p2,&o,NULL,NULL,NULL,v2,x2,r2,&s2,NULL,&counts);
 C(r==r2code&&t1.hash==t2.hash&&t1.count==t2.count&&counts.total==t1.count);
 C(!memcmp(v1,v2,sizeof v1)&&!memcmp(r1,r2,sizeof r1)&&!memcmp(x1,x2,(size_t)n*nev*sizeof(double)));
 C(s1.iterations==s2.iterations&&s1.op_count==s2.op_count&&s1.nconv==s2.nconv&&s1.nreturned==s2.nreturned);
 free(x1);free(x2);puts("KRYLOV LARGE PASS: n=32768 bitwise native differential");return 0;
}

int main(void)
{
 const char *which[]={"LM","SM","LA:ASC","SA:DESC","BE"};
 unsigned cases=0;
 for(int mode=0;mode<3;mode++)for(int lock=0;lock<2;lock++)for(int adapt=0;adapt<2;adapt++)for(int j=0;j<5;j++){
  C(!compare(mode,0,lock,adapt,which[j],(lock+adapt)%2,0));cases++;
 }
 for(int lock=0;lock<2;lock++)for(int adapt=0;adapt<2;adapt++){C(!compare(1,1,lock,adapt,"LM",1,0));cases++;}
 for(uint64_t f=1;f<=12;f++){C(!compare(1,0,0,1,"LM",1,f));cases++;}
 C(!cancellation());C(!interleaved());C(!large());
 printf("KRYLOV RESTORE PASS: %u three-way bitwise native differential cases, 20 cancellation boundaries, budget, independent contexts and workspace reuse\n",cases);
 return 0;
}
