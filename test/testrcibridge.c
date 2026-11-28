/* AI-assisted contribution. Proposed for inclusion under the LIS terms in COPYING;
 * licensing is subject to contributor confirmation before merge. */
#include "lis_rci_bridge.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define CHECK(c) do {if(!(c)){fprintf(stderr,"FAIL line %d: %s\n",__LINE__,#c);rc=1;goto cleanup;}}while(0)
/* Native adaptation of the archived mock bridge failure checks.
 * Use distinct assembled matrices rather than invalid structural mutations. */
static int bridge_faults(void)
{
 LIS_MATRIX A=NULL;lis_rci_bridge *p=NULL;int rc=0;
 double b[]={2,3,4},bad[]={NAN,3,4},x[]={.5,.5,.5};
 lis_rci_options o={3,30,3,1e-10,0,LIS_RCI_CG};lis_rci_info info;
 for(int variant=0;variant<3;variant++){
  CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));
  CHECK(!lis_matrix_set_size(A,0,3));
  for(LIS_INT i=0;i<3;i++)
   CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i,i?2+i:variant==0?0:variant==1?-2:2,A));
  CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
  o.n=4;CHECK(lis_rci_bridge_create(A,&o,1,&p)==LIS_RCI_BADARG&&!p);
  o.n=3;CHECK(!lis_rci_bridge_create(A,&o,1,&p));
  if(variant<2){
   CHECK(lis_rci_bridge_solve(p,b,x,&info)==LIS_RCI_BADARG);
   CHECK(info.status==LIS_RCI_BADARG&&x[0]==.5&&x[1]==.5&&x[2]==.5);
  }else{
   CHECK(lis_rci_bridge_solve(p,bad,x,&info)==LIS_RCI_NONFINITE);
   CHECK(info.status==LIS_RCI_NONFINITE&&x[0]==.5&&x[1]==.5&&x[2]==.5);
   CHECK(lis_rci_bridge_solve(p,b,x,&info)==LIS_RCI_DONE);
   for(size_t i=0;i<3;i++)CHECK(fabs(x[i]-1)<1e-10);
  }
  lis_rci_bridge_destroy(p);p=NULL;lis_matrix_destroy(A);A=NULL;
 }
 puts("RCI NATIVE BRIDGE dimension/zero/negative diagonal/nonfinite/recovery PASS");
cleanup:
 lis_rci_bridge_destroy(p);if(A)lis_matrix_destroy(A);return rc;
}

int main(int argc,char **argv){
 LIS_MATRIX A=NULL; lis_rci_bridge *p=NULL; int rc=0; const LIS_INT n=64; double b[64],x[64];
 if(lis_initialize(&argc,&argv)!=0)return 1;
 CHECK(!bridge_faults());
 CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,n));
 for(LIS_INT i=0;i<n;i++){CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i,4.0,A));if(i)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1.0,A));if(i+1<n)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1.0,A));}
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
 for(int method=LIS_RCI_CG;method<=LIS_RCI_FGMRES;method++)for(int jac=0;jac<=1;jac++){
  lis_rci_options o={(size_t)n,200,12,1e-10,0,(lis_rci_method)method};CHECK(!lis_rci_bridge_create(A,&o,jac,&p));
  for(int rep=0;rep<2;rep++){
   for(LIS_INT i=0;i<n;i++){b[i]=(4-(i>0)-(i+1<n))*(rep+1.0);x[i]=0;}
   lis_rci_info info;CHECK(lis_rci_bridge_solve(p,b,x,&info)==LIS_RCI_DONE);
   for(LIS_INT i=0;i<n;i++)CHECK(fabs(x[i]-(rep+1.0))<1e-8);
   printf("method=%d jacobi=%d reuse=%d iterations=%zu PASS\n",method,jac,rep,info.iterations);
  }
  lis_rci_bridge_destroy(p);p=NULL;
 }
 cleanup:lis_rci_bridge_destroy(p);if(A)lis_matrix_destroy(A);lis_finalize();return rc;
}
