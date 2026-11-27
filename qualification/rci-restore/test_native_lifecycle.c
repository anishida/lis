#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_rci_native_lifecycle.h"
#include <stdio.h>
#define C(x) do{if(!(x)){fprintf(stderr,"LIFECYCLE FAIL %d: %s\n",__LINE__,#x);return 1;}}while(0)
static int solve(const char *method)
{
 LIS_MATRIX A;LIS_VECTOR b,x,ref;LIS_SOLVER s,r;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,32));
 for(int i=0;i<32;i++){C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,3+.01*i,A));if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-.7,A));if(i<31)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-.7,A));}
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));C(!lis_vector_duplicate(A,&ref));C(!lis_vector_set_all(1,b));
 C(!lis_solver_create(&s));C(!lis_solver_create(&r));
 char opts[200];snprintf(opts,sizeof opts,"-i %s -p jacobi -tol 1e-11 -restart 8 -maxiter 300 -initx_zeros false",method);
 C(!lis_solver_set_option(opts,s));C(!lis_solver_set_option(opts,r));
 lis_rci_native_lifecycle *K=lis_rci_native_lifecycle_create();C(K);
 uint64_t last=0;
 for(int rep=0;rep<3;rep++){
  C(!lis_vector_set_all(.125,x));C(!lis_vector_set_all(.125,ref));
  C(!lis_solve(A,b,ref,r));C(r->retcode==LIS_SUCCESS);
  lis_rci_native_session_request q;
  int st=lis_rci_native_lifecycle_start(K,A,b,x,s,&q);C(st>0&&q.token>last);
  if(rep==0){
   last=q.token;C(lis_rci_native_lifecycle_cancel(K,last-1)==LIS_RCI_PROTOCOL);
   C(lis_rci_native_lifecycle_cancel(K,last)==LIS_RCI_CANCELLED);
   for(int i=0;i<32;i++)C(x->value[i]==.125);
   /* Keep classical warm-state history aligned with the cancellation. */
   lis_solver_destroy(r);C(!lis_solver_create(&r));C(!lis_solver_set_option(opts,r));continue;
  }
  int calls=0;
  while(st>0){
   C(q.token>last&&q.count==1&&q.representation==LIS_NR_VECTOR);last=q.token;
   C(lis_rci_native_lifecycle_advance(K,last-1,0,&q)==LIS_RCI_PROTOCOL);
   for(int i=0;i<32;i++)C(x->value[i]==.125);
   LIS_INT e=q.operation==LIS_RCI_OP_A?lis_matvec(q.A,q.in[0],q.out[0]):lis_psolve(q.solver,q.in[0],q.out[0]);calls++;
   st=lis_rci_native_lifecycle_advance(K,last,e,&q);
  }
  LIS_INT native;C(lis_rci_native_lifecycle_status(K,&native)==LIS_RCI_DONE&&st==LIS_RCI_DONE&&native==LIS_SUCCESS);
  C(s->iter==r->iter&&s->resid==r->resid);
  for(int i=0;i<32;i++)C(x->value[i]==ref->value[i]);
  printf("FULL LIFECYCLE PASS %s rep=%d requests=%d\n",method,rep,calls);
 }
 lis_rci_native_lifecycle_destroy(K);lis_solver_destroy(s);lis_solver_destroy(r);
 lis_vector_destroy(b);lis_vector_destroy(x);lis_vector_destroy(ref);lis_matrix_destroy(A);return 0;
}
int main(int argc,char **argv)
{
 C(!lis_initialize(&argc,&argv));
 const char *methods[]={"cg","bicgstab","gmres","fgmres","gcrodr","fgcrodr"};
 for(int i=0;i<6;i++)C(!solve(methods[i]));
 lis_finalize();return 0;
}
