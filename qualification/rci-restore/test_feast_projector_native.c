/* Native projector differential qualification; real double non-MPI. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_rci_native_feast_projector.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
extern LIS_INT lis_efeast_projector_native_reference(LIS_MATRIX,LIS_MATRIX,LIS_REAL,LIS_REAL,LIS_REAL,LIS_INT,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER);
extern LIS_INT lis_efeast_shifted_solve_split(LIS_MATRIX,LIS_MATRIX,LIS_REAL,LIS_REAL,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_VECTOR,LIS_SOLVER);
#define N 4
#define C(x) do{if(!(x)){fprintf(stderr,"FEAST PROJECTOR FAIL %d: %s\n",__LINE__,#x);rc=1;goto cleanup;}}while(0)
static int probe(int generalized)
{
 LIS_MATRIX A=NULL,B=NULL;LIS_SOLVER s=NULL;LIS_VECTOR v[6]={0};
 lis_rci_native_feast_projector *f=NULL;int rc=0;
 C(!lis_matrix_create(LIS_COMM_WORLD,&A));C(!lis_matrix_set_size(A,0,N));
 if(generalized){C(!lis_matrix_create(LIS_COMM_WORLD,&B));C(!lis_matrix_set_size(B,0,N));}
 const double lambda[N]={1,2,4,7};
 for(LIS_INT i=0;i<N;i++){
  double b=generalized?1+.2*i:1;
  C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b*lambda[i],A));
  if(B)C(!lis_matrix_set_value(LIS_INS_VALUE,i,i,b,B));
 }
 C(!lis_matrix_set_type(A,LIS_MATRIX_CSR));C(!lis_matrix_assemble(A));
 if(B){C(!lis_matrix_set_type(B,LIS_MATRIX_CSR));C(!lis_matrix_assemble(B));}
 for(int i=0;i<6;i++)C(!lis_vector_duplicate(A,&v[i]));
 for(LIS_INT i=0;i<N;i++){
  C(!lis_vector_set_value(LIS_INS_VALUE,i,i+1,v[0]));
  C(!lis_vector_set_value(LIS_INS_VALUE,i,.2*(i+1),v[1]));
 }
 C(!lis_solver_create(&s));
 C(!lis_solver_set_option("-i gmres -p none -restart 16 -maxiter 200 -tol 1e-13 -initx_zeros true",s));
 C(!lis_efeast_projector_native_reference(A,B,1.5,4.5,.3,16,v[0],v[1],v[2],v[3],s));
 f=lis_rci_native_feast_projector_create();C(f);
 uint64_t last=0;
 for(int fail=0;fail<=32+2*generalized;fail++){
  lis_rci_native_feast_projector_request q;
  int st=lis_rci_native_feast_projector_start(f,A,B,1.5,4.5,.3,16,v[0],v[1],v[4],v[5],s,1,&q);
  int step=0;
  while(st>0){
   C(q.token>last);last=q.token;
   C(lis_rci_native_feast_projector_advance(f,q.token-1,0,&q)==FP_PROTOCOL);
   if(++step==fail){st=lis_rci_native_feast_projector_cancel(f,q.token);break;}
   LIS_INT e;
   if(q.operation==FP_B)e=lis_matvec(q.A,q.in_real,q.out_real);
   else{
    C(q.operation==FP_SHIFT);
    /* Write-only callback result; native solve is called with its usual zero guess. */
    C(!lis_vector_set_all(0,q.out_real));C(!lis_vector_set_all(0,q.out_imag));
    e=lis_efeast_shifted_solve_split(q.A,q.B,q.zr,q.zi,q.in_real,q.in_imag,q.out_real,q.out_imag,q.solver);
   }
   C(!e);
   st=lis_rci_native_feast_projector_advance(f,q.token,e,&q);
  }
  if(fail){C(st==FP_CANCELLED);continue;}
  C(st==FP_DONE&&step==32+2*generalized);
  C(!memcmp(v[2]->value,v[4]->value,N*sizeof(LIS_SCALAR)));
  C(!memcmp(v[3]->value,v[5]->value,N*sizeof(LIS_SCALAR)));
 }
 printf("NATIVE FEAST PROJECTOR PASS generalized=%d cancellations=%d\n",generalized,32+2*generalized);
cleanup:
 if(f)lis_rci_native_feast_projector_destroy(f);
 if(s)lis_solver_destroy(s);
 for(int i=0;i<6;i++)if(v[i])lis_vector_destroy(v[i]);
 if(B)lis_matrix_destroy(B);
 if(A)lis_matrix_destroy(A);
 return rc;
}
int main(int argc,char **argv){if(lis_initialize(&argc,&argv))return 1;int rc=probe(0)||probe(1);lis_finalize();return rc;}
