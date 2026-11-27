/* Exercise the unified typed session through genuine native setup. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_esolver.h"
#include "lis_rci_native_session.h"
#include "lis_rci_native_linear_restore.h"
#include "lis_rci_native_gcro_restore.h"
#include "lis_rci_native_fgcro_restore.h"
#include "lis_rci_native_feast.h"
#include <stdio.h>
#define C(x) do{if(!(x)){fprintf(stderr,"UNIFIED FAIL %d: %s\n",__LINE__,#x);return 1;}}while(0)
static int seen,failed;
static LIS_INT run(lis_rci_native_session *s,int st,lis_rci_native_session_request q)
{
 uint64_t last=0,a,p;LIS_INT native=0;
 while(st>0){
  lis_rci_native_session_request pending;
  if(q.token<=last||lis_rci_native_session_pending(s,&pending)!=st||pending.token!=q.token||
     lis_rci_native_session_advance(s,q.token-1,0,&q)!=LIS_RCI_PROTOCOL||q.token!=pending.token){failed=1;break;}
  last=q.token;seen++;
  LIS_INT rc=0;
  if(q.representation==LIS_NR_SPLIT){
   if(q.operation!=LIS_RCI_OP_SHIFTED_SOLVE||q.shift_sign!=1){failed=1;break;}
   lis_vector_set_all(0,q.out_real);lis_vector_set_all(0,q.out_imag);
   rc=lis_efeast_shifted_solve_split(q.A,q.B,q.zr,q.zi,q.in_real,q.in_imag,q.out_real,q.out_imag,q.solver);
  }else if(q.representation==LIS_NR_VECTOR){
   if(q.count<=0||!q.in||!q.out){failed=1;break;}
   rc=q.operation==LIS_RCI_OP_PRECON?lis_psolve(q.solver,q.in[0],q.out[0]):
    q.block?lis_matvec_block(q.A,q.count,q.in,q.out):lis_matvec(q.A,q.in[0],q.out[0]);
  }else {failed=1;break;}
  st=lis_rci_native_session_advance(s,q.token,rc,&q);
 }
 lis_rci_native_session_status(s,&native,&a,&p,NULL);
 if(st!=LIS_RCI_DONE)failed=1;
 lis_rci_native_session_destroy(s);return failed?LIS_BREAKDOWN:native;
}
static LIS_INT linear(LIS_SOLVER solver,LIS_INT method)
{
 int kind=method==LIS_SOLVER_GCRODR?LIS_NR_GCRO:method==LIS_SOLVER_FGCRODR?LIS_NR_FGCRO:LIS_NR_LINEAR;
 lis_rci_native_session *s=lis_rci_native_session_create(kind);if(!s)return LIS_ERR_OUT_OF_MEMORY;
 lis_rci_native_session_request q;int st=lis_rci_native_session_start_linear(s,solver,method,&q);return run(s,st,q);
}
static LIS_INT gcro(LIS_SOLVER s){return linear(s,LIS_SOLVER_GCRODR);}
static LIS_INT fgcro(LIS_SOLVER s){return linear(s,LIS_SOLVER_FGCRODR);}
static LIS_INT feast(LIS_ESOLVER e)
{
 lis_rci_native_session *s=lis_rci_native_session_create(LIS_NR_FEAST);if(!s)return LIS_ERR_OUT_OF_MEMORY;
 lis_rci_native_session_request q;int st=lis_rci_native_session_start_feast(s,e,&q);return run(s,st,q);
}
static int diag(void *ctx,int64_t n,const double *x,double *y)
{(void)ctx;for(int64_t i=0;i<n;i++)y[i]=(i+1)*x[i];return 0;}
static int identity(void *ctx,int64_t n,const double *x,double *y)
{(void)ctx;for(int64_t i=0;i<n;i++)y[i]=x[i];return 0;}
static int inverse(void *ctx,int64_t n,const double *x,double *y)
{(void)ctx;for(int64_t i=0;i<n;i++)y[i]=x[i]/(i+1-.5);return 0;}
static int krylov(void)
{
 lis_krylov_problem64 p={12,LIS_KRYLOV_SHIFT_INVERT_B,.5,NULL,diag,identity,inverse};
 lis_krylov_options o;lis_krylov_options_default(&o);o.nev=2;o.ncv=6;o.inner_tol=1e-12;o.maxiter=100;
 double v[2],x[24],r[2];lis_krylov_stats stats;
 lis_rci_native_session *s=lis_rci_native_session_create(LIS_NR_KRYLOV);C(s);
 lis_rci_native_session_request q;uint64_t last=0;
 for(int rep=0;rep<2;rep++){
  int st=lis_rci_native_session_start_krylov(s,&p,&o,NULL,NULL,NULL,v,x,r,&stats,NULL,&q);
  if(st<=0){LIS_INT native;uint64_t a,p;lis_rci_native_session_status(s,&native,&a,&p,NULL);fprintf(stderr,"KRYLOV START rep=%d status=%d native=%lld\n",rep,st,(long long)native);}
  C(st>0&&q.token>last);
  if(!rep){last=q.token;C(lis_rci_native_session_cancel(s,last-1)==LIS_RCI_PROTOCOL);C(lis_rci_native_session_cancel(s,last)==LIS_RCI_CANCELLED);C(stats.nreturned==0);continue;}
  int shifted=0;
  while(st>0){
   C(q.representation==LIS_NR_RAW64&&q.n64==12&&q.token>last);
   last=q.token;C(lis_rci_native_session_advance(s,last-1,0,&q)==LIS_RCI_PROTOCOL);
   int rc;
   if(q.operation==LIS_RCI_OP_SHIFTED_SOLVE){C(q.shift_sign==-1&&q.sigma==.5);shifted++;rc=inverse(NULL,q.n64,q.raw_in,q.raw_out);}
   else if(q.operation==LIS_RCI_OP_B)rc=identity(NULL,q.n64,q.raw_in,q.raw_out);
   else{C(q.operation==LIS_RCI_OP_A);rc=diag(NULL,q.n64,q.raw_in,q.raw_out);}
   st=lis_rci_native_session_advance(s,last,rc,&q);
  }
  C(st==LIS_RCI_DONE&&shifted>0&&stats.nreturned==2);
 }
 lis_rci_native_session_destroy(s);return 0;
}
int main(int argc,char **argv)
{
 C(!lis_initialize(&argc,&argv));
 LIS_MATRIX A;LIS_VECTOR b,x;C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,12));
 for(int i=0;i<12;i++)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,i+1,A));
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));C(!lis_vector_duplicate(A,&b));C(!lis_vector_duplicate(A,&x));
 C(!lis_vector_set_all(1,b));
 lis_rci_native_linear_test_select(0,linear);lis_rci_native_gcro_test_select(0,gcro);lis_rci_native_fgcro_test_select(0,fgcro);
 const char *methods[]={"cg","bicgstab","gmres","fgmres","gcrodr","fgcrodr"};
 for(int i=0;i<6;i++){
  LIS_SOLVER s;C(!lis_solver_create(&s));char opts[160];
  snprintf(opts,sizeof opts,"-i %s -p none -restart 8 -gcrodr_keep 2 -fgcrodr_k 2 -tol 1e-10 -maxiter 300",methods[i]);
  C(!lis_solver_set_option(opts,s));C(!lis_vector_set_all(0,x));seen=failed=0;
  (void)lis_solve(A,b,x,s);C(!failed&&seen>0&&s->retcode==LIS_SUCCESS);
  for(int j=0;j<12;j++)C(fabs(x->value[j]-1.0/(j+1))<1e-8);
  printf("UNIFIED PASS %s requests=%d\n",methods[i],seen);lis_solver_destroy(s);
 }
 lis_rci_native_linear_test_select(0,NULL);lis_rci_native_gcro_test_select(0,NULL);lis_rci_native_fgcro_test_select(0,NULL);
 LIS_ESOLVER e;LIS_SCALAR ev;C(!lis_esolver_create(&e));C(!lis_vector_set_all(1,x));
 C(!lis_esolver_set_option("-e feast -ss 6 -emaxiter 20 -etol 1e-7 -feast_emin 1.5 -feast_emax 3.5 -feast_points 16",e));
 lis_rci_native_feast_test_select(0,feast);seen=failed=0;
 C(!lis_esolve(A,x,&ev,e));C(!failed&&seen>0&&e->retcode==LIS_SUCCESS&&e->nev==2);
 printf("UNIFIED PASS FEAST requests=%d\n",seen);lis_rci_native_feast_test_select(0,NULL);lis_esolver_destroy(e);
 C(!krylov());puts("UNIFIED PASS KRYLOV");
 lis_vector_destroy(b);lis_vector_destroy(x);lis_matrix_destroy(A);lis_finalize();return 0;
}
