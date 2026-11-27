/* AI-assisted continuation reconstruction from the actual native LIS sources.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_rci_native_gcro_restore.h"
#include "lislib.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#define N 48
#define K 3
#define C(e) do{if(!(e)){fprintf(stderr,"GCRO DIFF FAIL %d: %s\n",__LINE__,#e);rc=1;goto cleanup;}}while(0)
typedef struct {
 LIS_INT status,iter,dim,action,cross,builds,updates,rejects;
 LIS_REAL residual,quality,minquality;
 LIS_SCALAR x[N],basis[3*K*N];
 uint64_t mat,pre;
} answer;
static int equal_answer(const answer *a,const answer *b)
{
 if(a->status!=b->status||a->iter!=b->iter||a->dim!=b->dim||a->action!=b->action||a->cross!=b->cross||a->builds!=b->builds||a->updates!=b->updates||a->rejects!=b->rejects||a->residual!=b->residual||a->quality!=b->quality||a->minquality!=b->minquality||a->mat!=b->mat||a->pre!=b->pre)return 0;
 for(size_t i=0;i<N;i++)if(a->x[i]!=b->x[i])return 0;
 for(size_t i=0;i<3*K*N;i++)if(a->basis[i]!=b->basis[i])return 0;
 return 1;
}

static uint64_t mat,pre,stop;
static int fault_kind,forced_status;
static LIS_INT external(LIS_SOLVER s)
{
 lis_rci_native_gcro_restore *g=lis_rci_native_gcro_restore_create();
 if(!g)return LIS_ERR_OUT_OF_MEMORY;
 lis_rci_native_linear_request q;
 int st=lis_rci_native_gcro_restore_start(g,s,1,&q);
 uint64_t last=0,step=0;
 while(st>0){
  if(q.token<=last||lis_rci_native_gcro_restore_advance(g,q.token-1,0,&q)!=LIS_RCI_PROTOCOL){
   lis_rci_native_gcro_restore_destroy(g);return 99;
  }
  last=q.token;
  if(++step==stop){
   if(fault_kind==1)st=lis_rci_native_gcro_restore_cancel(g,q.token);
   else if(fault_kind==2){
    q.out->value[0]=0;
    st=lis_rci_native_gcro_restore_advance(g,q.token,0,&q);
   }else st=lis_rci_native_gcro_restore_advance(g,q.token,73,&q);
   forced_status=st;
   break;
  }
  LIS_INT e=q.operation==LIS_RCI_OP_A?lis_matvec(q.A,q.in,q.out):lis_psolve(q.solver,q.in,q.out);
  st=lis_rci_native_gcro_restore_advance(g,q.token,e,&q);
 }
 LIS_INT native=0;int status=lis_rci_native_gcro_restore_status(g,&native,&mat,&pre);
 if(status!=st)native=98;
 lis_rci_native_gcro_restore_destroy(g);return native;
}
static int sequence(int mode,answer out[3],int fault_rep,uint64_t fail_at,int kind)
{
 LIS_MATRIX A=NULL;LIS_VECTOR b=NULL,x=NULL;LIS_SOLVER s=NULL;
 int rc=0;memset(out,0,3*sizeof(*out));
 C(!lis_solver_create(&s));
 C(!lis_solver_set_option("-i gcrodr -p jacobi -restart 8 -gcrodr_keep 3 -maxiter 300 -tol 1e-11 -initx_zeros false",s));
 lis_rci_native_gcro_test_select(mode==0,mode==2?external:NULL);
 for(int rep=0;rep<3;rep++){
  double shift=.015*rep;
  C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
  for(LIS_INT i=0;i<N;i++){
   C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,2.6+.008*i+shift,A));
   if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-.7,A));
   if(i+1<N)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.1,A));
  }
  C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
  C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));
  LIS_SCALAR exact[N];for(int i=0;i<N;i++)exact[i]=sin(.19*i+.1*rep)+.4;
#ifdef _COMPLEX
  for(int j=0;j<N;j++)exact[j]+=I*(.15+.007*j);
#endif
  for(LIS_INT i=0;i<N;i++){
   LIS_SCALAR rhs=(2.6+.008*i+shift)*exact[i]-(i?.7*exact[i-1]:0)-(i+1<N?1.1*exact[i+1]:0);
   C(!lis_vector_set_value(LIS_INS_VALUE,i,rhs,b));C(!lis_vector_set_value(LIS_INS_VALUE,i,.1,x));
  }
  stop=rep==fault_rep?fail_at:0;fault_kind=kind;forced_status=0;mat=pre=0;
  (void)lis_solve(A,b,x,s);
  answer *a=out+rep;
  C(!lis_solver_get_status(s,&a->status));
  if(stop){
   C(forced_status==(kind==1?LIS_RCI_CANCELLED:kind==2?LIS_RCI_NONFINITE:LIS_RCI_CALLBACK));
   C(a->status==(kind==3?73:LIS_BREAKDOWN)&&!s->recycle_valid);
   goto cleanup;
  }
  C(a->status==LIS_SUCCESS);
  C(!lis_solver_get_iter(s,&a->iter));C(!lis_solver_get_residualnorm(s,&a->residual));
  C(!lis_solver_get_gcrodr_recycle(s,&a->dim,&a->action,&a->quality,&a->minquality));
  C(a->dim>0&&a->dim<=K);
  a->cross=s->recycle_cross_uses;a->builds=s->recycle_builds;a->updates=s->recycle_updates;a->rejects=s->recycle_rejects;
  C(!lis_vector_get_values(x,0,N,a->x));
  for(int i=0;i<N;i++)C(fabs(a->x[i]-exact[i])<1e-8);
  for(LIS_INT j=0;j<a->dim;j++){
   memcpy(a->basis+j*N,s->recycle_u[j]->value,N*sizeof(LIS_SCALAR));
   memcpy(a->basis+K*N+j*N,s->recycle_zu[j]->value,N*sizeof(LIS_SCALAR));
   memcpy(a->basis+2*K*N+j*N,s->recycle_c[j]->value,N*sizeof(LIS_SCALAR));
  }
  a->mat=mat;a->pre=pre;
  if(!fail_at)printf("GCRO mode=%d rhs=%d iter=%lld dim=%lld action=%lld cross=%lld A=%llu P=%llu residual=%.17g\n",mode,rep,(long long)a->iter,(long long)a->dim,(long long)a->action,(long long)a->cross,(unsigned long long)mat,(unsigned long long)pre,(double)a->residual);
  lis_vector_destroy(b);b=NULL;lis_vector_destroy(x);x=NULL;lis_matrix_destroy(A);A=NULL;
 }
cleanup:
 lis_rci_native_gcro_test_select(0,NULL);
 if(s)lis_solver_destroy(s);if(x)lis_vector_destroy(x);if(b)lis_vector_destroy(b);if(A)lis_matrix_destroy(A);
 return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 answer a[3],b[3],c[3],scratch[3];int rc=0;
 C(!sequence(0,a,-1,0,0));C(!sequence(1,b,-1,0,0));C(!sequence(2,c,-1,0,0));
 uint64_t requests[3];
 for(int i=0;i<3;i++){
  C(c[i].mat&&c[i].pre);requests[i]=c[i].mat+c[i].pre;c[i].mat=c[i].pre=0;
  C(equal_answer(a+i,b+i)&&equal_answer(a+i,c+i));
 }
 C(a[1].cross>a[0].cross&&a[2].cross>a[1].cross);
 /* All observed callbacks on both a cold and a warm solve, including
  * initial residual, Arnoldi, policy assessment, refresh and update.
  * Failures must use the original cleanup and invalidate partial recycle. */
 uint64_t faults=0;
 for(int rep=0;rep<2;rep++)for(uint64_t f=1;f<=requests[rep];f++)
  for(int kind=1;kind<=3;kind++){C(!sequence(2,scratch,rep,f,kind));faults++;}
 printf("NATIVE GCRO PASS: three-way exact-value sequence and recycle basis, %llu phase fault scenarios\n",(unsigned long long)faults);
cleanup:
 lis_finalize();return rc;
}
