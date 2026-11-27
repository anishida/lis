/* AI-assisted reconstruction from historical design descriptions, not recovered RC3 code.
 * Proposed for inclusion under LIS terms in COPYING, subject to contributor confirmation. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis_rci_native_restore.h"
#include "lis_precon.h"
#include <stdio.h>
#include <math.h>
#define N 64
#define C(e) do{if(!(e)){fprintf(stderr,"NATIVE FAIL %d: %s\n",__LINE__,#e);rc=1;goto cleanup;}}while(0)
static int run(int method,const char *precon)
{
 LIS_MATRIX A=NULL,B=NULL;LIS_SOLVER solver=NULL;LIS_PRECON pc=NULL;
 lis_rci_context *c=NULL;lis_rci_native_restore *s=NULL;
 int rc=0;
 double x[N],b[N],exact[N],bx[N],out[N];char opts[96];
 const char *names[]={"","cg","bicgstab","gmres","fgmres"};
 double lower=-.5,upper=method==LIS_RCI_CG?-.5:-1.1;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 C(!lis_matrix_create(LIS_COMM_WORLD,&B));C(!lis_matrix_set_size(B,0,N));
 for(LIS_INT i=0;i<N;i++){
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,4+.01*i,A));
  if(i)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,lower,A));
  if(i+1<N)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,upper,A));
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,1+.01*i,B));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(B));
 C(!lis_solver_create(&solver));
 snprintf(opts,sizeof opts,"-i %s -p %s -ssor_w 1.0",names[method],precon);
 C(!lis_solver_set_option(opts,solver));
 C(!lis_solver_set_matrix(A,solver));
 {LIS_PRECON temp=NULL;C(!lis_precon_create(solver,&temp));pc=temp;solver->precon=pc;}
 s=lis_rci_native_restore_create(A,B,solver);C(s);
 lis_rci_options o={N,200,12,1e-11,0,(lis_rci_method)method};
 c=lis_rci_create(&o);C(c);
 for(int rep=0;rep<2;rep++){
  for(int i=0;i<N;i++){exact[i]=sin(.1*i+rep)+.2;x[i]=.1;}
  for(int i=0;i<N;i++)b[i]=(4+.01*i)*exact[i]+(i?lower*exact[i-1]:0)+(i+1<N?upper*exact[i+1]:0);
  lis_rci_info info;
  C(!lis_rci_native_restore_solve(s,c,b,x,&info));C(info.precons>0);
  for(int i=0;i<N;i++)C(fabs(x[i]-exact[i])<1e-9);
  C(!lis_rci_native_restore_apply(s,LIS_RCI_OP_B,N,x,bx));
  for(int i=0;i<N;i++)C(fabs(bx[i]-(1+.01*i)*x[i])<1e-12);
 }
 C(!lis_matrix_mark_values_changed(A));
 for(int i=0;i<N;i++)out[i]=-3;
 C(lis_rci_native_restore_apply(s,LIS_RCI_OP_A,N,x,out)==LIS_RCI_PROTOCOL);
 lis_rci_info info;C(lis_rci_native_restore_solve(s,c,b,out,&info)==LIS_RCI_PROTOCOL);
 for(int i=0;i<N;i++)C(out[i]==-3);
cleanup:
 lis_rci_native_restore_destroy(s);lis_rci_destroy(c);
 if(solver)solver->precon=NULL;
 if(pc)lis_precon_destroy(pc);
 if(solver)lis_solver_destroy(solver);
 if(B)lis_matrix_destroy(B);
 if(A)lis_matrix_destroy(A);
 return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 int rc=0;for(int m=1;m<=4;m++)for(int p=0;p<2;p++)rc|=run(m,p?"ssor":"jacobi");
 lis_finalize();if(!rc)puts("NATIVE RESTORE PASS: 16 real LIS solves, native Jacobi/SSOR psolve, B matvec, generation invalidation");
 return rc;
}
