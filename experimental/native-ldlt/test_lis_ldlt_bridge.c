#include "lis_ldlt_user_bridge.h"
#include <math.h>
#include <stdio.h>

#define OK(x) do{rc=(x);if(rc){fprintf(stderr,"line %d rc=%d\n",__LINE__,(int)rc);goto done;}}while(0)
/* Two independent indefinite 2x2 blocks with zero diagonal: require
 * 2x2 pivots, rather than a positive-definite surrogate. */
static int run(int mode,int method)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL,r=NULL;
    LIS_SOLVER s=NULL;
    LIS_INT rc=0,status,i;
    LIS_REAL residual,bnorm;
    int failed=1;
    OK(lis_matrix_create(LIS_COMM_WORLD,&A));
    OK(lis_matrix_set_size(A,0,4));
    for(i=0;i<4;i++){
        double v=(mode==2)?0.0:2.0+i/2;
        if(mode==1&&i==1)v=3.0;
        OK(lis_matrix_set_value(LIS_INS_VALUE,i,i,0.0,A));
        OK(lis_matrix_set_value(LIS_INS_VALUE,i,i^1,v,A));
    }
    OK(lis_matrix_set_type(A,LIS_MATRIX_CSR));
    OK(lis_matrix_assemble(A));
    OK(lis_vector_duplicate(A,&b));OK(lis_vector_duplicate(A,&x));OK(lis_vector_duplicate(A,&r));
    OK(lis_vector_set_all(1.0,b));OK(lis_vector_set_all(0.0,x));
    OK(lis_solver_create(&s));
    OK(lis_solver_set_option(method==0?"-i fgmres -p ldltexp -scale none -tol 1e-12 -maxiter 30":
       method==1?"-i gmres -p ldltexp -scale none -tol 1e-12 -maxiter 30":
       "-i cg -p ldltexp -scale none -tol 1e-12 -maxiter 30",s));
    rc=lis_solve(A,b,x,s);
    if(mode||method==2){
        if(!rc){OK(lis_solver_get_status(s,&status));if(status==LIS_SUCCESS)goto done;}
        failed=0;goto done;
    }
    if(rc)goto done;
    OK(lis_solver_get_status(s,&status));if(status!=LIS_SUCCESS)goto done;
    OK(lis_matvec(A,x,r));OK(lis_vector_axpy(-1.0,b,r));
    OK(lis_vector_nrm2(r,&residual));OK(lis_vector_nrm2(b,&bnorm));
    if(!isfinite(residual)||residual/bnorm>1e-12)goto done;
    printf("physical_relative_residual=%.17g\n",(double)(residual/bnorm));
    failed=0;
done:
    if(s)lis_solver_destroy(s);
    if(r)lis_vector_destroy(r);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    printf("LDLT_LIS mode=%d method=%d %s\n",mode,method,failed?"FAIL":"PASS");
    return failed;
}
int main(int argc,char **argv)
{
    int fail;
    if(lis_initialize(&argc,&argv))return 1;
    if(lis_ldlt_experimental_register()){lis_finalize();return 1;}
    fail=run(0,0);fail|=run(0,1);fail|=run(1,0);fail|=run(2,0);fail|=run(0,2);
    lis_finalize();
    if(!fail)puts("LDLT_LIS_BRIDGE=PASS");
    return fail;
}
