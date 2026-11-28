/* End-to-end classic dispatch benchmark; identical client for both libraries. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stddef.h>
#define CHECK(x) do{if(!(x)){fprintf(stderr,"BENCH FAIL %d %s\n",__LINE__,#x);rc=1;goto cleanup;}}while(0)
static int bench(const char *method,LIS_INT n,int repeat)
{
 LIS_MATRIX A=NULL;LIS_VECTOR b=NULL,x=NULL;LIS_SOLVER s=NULL;
 double *times=calloc((size_t)repeat,sizeof(*times));int rc=0;
 CHECK(times);
 CHECK(!lis_matrix_create(LIS_COMM_WORLD,&A));CHECK(!lis_matrix_set_size(A,0,n));
 for(LIS_INT i=0;i<n;i++){
  CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i,2.1+.01*(i%17),A));
  if(i)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i-1,-1,A));
  if(i+1<n)CHECK(!lis_matrix_set_value(LIS_INS_VALUE,i,i+1,-1,A));
 }
 CHECK(!lis_matrix_set_type(A,LIS_MATRIX_CSR));CHECK(!lis_matrix_assemble(A));
 CHECK(!lis_vector_duplicate(A,&b));CHECK(!lis_vector_duplicate(A,&x));
 CHECK(!lis_vector_set_all(1,x));CHECK(!lis_matvec(A,x,b));
 LIS_INT it=0;LIS_REAL res=0;
 for(int j=-1;j<repeat;j++){
  CHECK(!lis_solver_create(&s));
  char opts[256];snprintf(opts,sizeof opts,"-i %s -p jacobi -restart 24 -maxiter 500 -tol 1e-10 -initx_zeros true",method);
  CHECK(!lis_solver_set_option(opts,s));CHECK(!lis_vector_set_all(0,x));
  double start=lis_wtime();
  (void)lis_solve(A,b,x,s);
  double elapsed=lis_wtime()-start;
  LIS_INT status;CHECK(!lis_solver_get_status(s,&status));CHECK(status==LIS_SUCCESS);
  CHECK(!lis_solver_get_iter(s,&it));CHECK(!lis_solver_get_residualnorm(s,&res));
  for(LIS_INT i=0;i<n;i++)CHECK(fabs(x->value[i]-1)<1e-7);
  if(j>=0){times[j]=elapsed;printf("RESULT %s %lld %d %lld %a\n",method,(long long)n,j,(long long)it,(double)res);}
  lis_solver_destroy(s);s=NULL;
 }
 for(int i=0;i<repeat;i++)for(int j=i+1;j<repeat;j++)if(times[j]<times[i]){double t=times[i];times[i]=times[j];times[j]=t;}
 printf("TIME %s %lld %.9g\n",method,(long long)n,times[repeat/2]);
cleanup:
 if(s)lis_solver_destroy(s);if(x)lis_vector_destroy(x);if(b)lis_vector_destroy(b);if(A)lis_matrix_destroy(A);free(times);return rc;
}
int main(int argc,char **argv)
{
 if(lis_initialize(&argc,&argv))return 1;
 const char *methods[]={"cg","bicgstab","gmres","fgmres","gcrodr","fgcrodr"};
 int rc=0;for(int i=0;i<6;i++)for(int j=0;j<2;j++)rc|=bench(methods[i],j?32768:256,7);
 lis_finalize();return rc;
}
