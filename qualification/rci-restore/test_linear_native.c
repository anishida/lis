/* AI-assisted continuation reconstruction from the actual native LIS sources.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_rci_native_linear_restore.h"
#include "lislib.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define N 37
#define MAXH 201
#define C(e) do{if(!(e)){fprintf(stderr,"LINEAR DIFF FAIL %d: %s\n",__LINE__,#e);rc=1;goto cleanup;}}while(0)
typedef struct {uint64_t hash,calls,fail;} trace;
typedef struct {LIS_INT rc,iter; LIS_REAL residual,h[MAXH]; LIS_SCALAR x[N]; trace t;uint64_t a,p;} result;
static void mix(trace *t,const void *p,size_t n)
{const unsigned char *a=p;for(size_t i=0;i<n;i++){t->hash^=a[i];t->hash*=UINT64_C(1099511628211);}}
static LIS_INT matvec(void *u,const LIS_SCALAR *x,LIS_SCALAR *y)
{
 trace *t=u;t->calls++;mix(t,x,N*sizeof(*x));
 if(t->calls==t->fail)return 73;
 for(int i=0;i<N;i++)y[i]=(4+.01*i)*x[i]-(i?.7*x[i-1]:0)-(i+1<N?.7*x[i+1]:0);
 mix(t,y,N*sizeof(*y));return LIS_SUCCESS;
}
static uint64_t a_count,p_count;
static int force_cancel,force_incomplete;
static LIS_INT external(LIS_SOLVER solver,LIS_INT method)
{
 lis_rci_native_linear_restore *K=lis_rci_native_linear_restore_create();
 if(!K)return LIS_ERR_OUT_OF_MEMORY;
 lis_rci_native_linear_request q;
 int s=lis_rci_native_linear_restore_start_method(K,solver,method,1,&q);
 uint64_t last=0;int step=0;
 while(s>0){
  if(q.token<=last||lis_rci_native_linear_restore_advance(K,q.token-1,0,&q)!=LIS_RCI_PROTOCOL||q.token<=last){
   lis_rci_native_linear_restore_destroy(K);return 99;
  }
  last=q.token;++step;
  if(force_cancel&&step==2){
   if(lis_rci_native_linear_restore_cancel(K,q.token)!=LIS_RCI_CANCELLED){lis_rci_native_linear_restore_destroy(K);return 98;}
   s=LIS_RCI_CANCELLED;break;
  }
  LIS_INT e=0;
  if(force_incomplete&&step==2)q.out->value[0]=0;
  else e=q.operation==LIS_RCI_OP_A?lis_matvec(q.A,q.in,q.out):lis_psolve(q.solver,q.in,q.out);
  s=lis_rci_native_linear_restore_advance(K,q.token,e,&q);
 }
 LIS_INT native=0;
 int status=lis_rci_native_linear_restore_status(K,&native,&a_count,&p_count);
 if(status!=s||(force_cancel&&s!=LIS_RCI_CANCELLED)||(force_incomplete&&s!=LIS_RCI_NONFINITE))native=97;
 lis_rci_native_linear_restore_destroy(K);return native;
}
static int solve_case(const char *method,int mode,int zero,int maxiter,uint64_t fail,int user,result *r)
{
 LIS_MATRIX A=NULL;LIS_VECTOR b=NULL,x=NULL;LIS_SOLVER s=NULL;
 int rc=0;
 memset(r,0,sizeof(*r));r->t.hash=UINT64_C(1469598103934665603);r->t.fail=fail;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 if(user){C(!lis_matrix_set_user(A,&r->t,matvec,matvec));}
 else{
  for(LIS_INT i=0;i<N;i++){
   C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,4+.01*i,A));
   if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-.7,A));
   if(i+1<N)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-.7,A));
  }
  C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 }
 C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));C(!lis_solver_create(&s));
 double exact[N],rhs[N];for(int i=0;i<N;i++)exact[i]=sin(.17*i)+.3;
 trace init={0};C(!matvec(&init,exact,rhs));
 for(LIS_INT i=0;i<N;i++){C(!lis_vector_set_value(LIS_INS_VALUE,i,rhs[i],b));C(!lis_vector_set_value(LIS_INS_VALUE,i,.1,x));}
 char opts[200];snprintf(opts,sizeof(opts),"-i %s -p %s -tol 1e-11 -maxiter %d -restart 5 -initx_zeros %s -print mem",method,user?"none":"jacobi",maxiter,zero?"true":"false");
 C(!lis_solver_set_option(opts,s));
 a_count=p_count=0;lis_rci_native_linear_test_select(mode==0,mode==2?external:NULL);
 LIS_INT api=lis_solve(A,b,x,s);
 (void)api;
 C(!lis_solver_get_status(s,&r->rc));C(!lis_solver_get_iter(s,&r->iter));C(!lis_solver_get_residualnorm(s,&r->residual));
 C(!lis_vector_get_values(x,0,N,r->x));r->a=a_count;r->p=p_count;
 if(s->rhistory){
  LIS_INT count=r->iter;if(count>maxiter)count=maxiter;if(count<0)count=0;
  /* On callback failure the final history slot may never be written. */
  if(fail&&count>0)count--;
  if(count>0)memcpy(r->h,s->rhistory,(size_t)count*sizeof(LIS_REAL));
 }
 if(!fail&&!force_cancel&&!force_incomplete&&maxiter>2){
  C(r->rc==LIS_SUCCESS);
  for(int i=0;i<N;i++)C(fabs(r->x[i]-exact[i])<1e-8);
 }
cleanup:
 lis_rci_native_linear_test_select(0,NULL);
 if(s)lis_solver_destroy(s);if(x)lis_vector_destroy(x);if(b)lis_vector_destroy(b);if(A)lis_matrix_destroy(A);
 return rc;
}
static int compare(const char *method,int zero,int maxiter,uint64_t fail,int user)
{
 int rc=0;result a,b,c;
 C(!solve_case(method,0,zero,maxiter,fail,user,&a));
 C(!solve_case(method,1,zero,maxiter,fail,user,&b));
 C(!solve_case(method,2,zero,maxiter,fail,user,&c));
 if(a.rc!=b.rc||a.rc!=c.rc||a.iter!=b.iter||a.iter!=c.iter||memcmp(a.x,b.x,sizeof a.x)||memcmp(a.x,c.x,sizeof a.x))
  fprintf(stderr,"method=%s zero=%d maxiter=%d fail=%llu user=%d rc=%lld/%lld/%lld it=%lld/%lld/%lld\n",method,zero,maxiter,(unsigned long long)fail,user,(long long)a.rc,(long long)b.rc,(long long)c.rc,(long long)a.iter,(long long)b.iter,(long long)c.iter);
 C(a.rc==b.rc&&a.rc==c.rc&&a.iter==b.iter&&a.iter==c.iter);
 C((isnan(a.residual)&&isnan(b.residual)&&isnan(c.residual))||(a.residual==b.residual&&a.residual==c.residual));
 C(!memcmp(a.x,b.x,sizeof a.x)&&!memcmp(a.x,c.x,sizeof a.x));
 C(!memcmp(a.h,b.h,sizeof a.h)&&!memcmp(a.h,c.h,sizeof a.h));
 C(a.t.calls==b.t.calls&&a.t.calls==c.t.calls&&a.t.hash==b.t.hash&&a.t.hash==c.t.hash);
 C(c.a>0||(!fail&&a.rc==LIS_SUCCESS));C(c.p>0||fail==1);
cleanup:return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 const char *methods[]={"cg","bicgstab","gmres","fgmres"};
 int rc=0,tests=0;
 for(int m=0;m<4;m++)for(int z=0;z<2;z++){
  rc|=compare(methods[m],z,200,0,0);tests++;
  rc|=compare(methods[m],z,200,0,1);tests++;
  rc|=compare(methods[m],z,2,0,1);tests++;
  for(uint64_t f=1;f<=3;f++){rc|=compare(methods[m],z,200,f,1);tests++;}
 }
 for(int m=0;m<4;m++){
  result r;force_cancel=1;rc|=solve_case(methods[m],2,0,200,0,1,&r);force_cancel=0;
  if(r.rc!=LIS_BREAKDOWN)rc=1;
  force_incomplete=1;rc|=solve_case(methods[m],2,0,200,0,1,&r);force_incomplete=0;
  if(r.rc!=LIS_BREAKDOWN)rc=1;
 }
 lis_finalize();
 if(!rc)printf("NATIVE LINEAR PASS: %d three-way differential solves (CSR/Jacobi, USER/identity, nonzero/zero guess, maxiter, callback error), 8 cancel/incomplete cases\n",tests);
 return rc;
}
