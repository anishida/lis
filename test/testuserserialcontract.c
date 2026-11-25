#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <stdio.h>
#include <math.h>
#if !defined(USE_MPI) && !defined(_COMPLEX)
typedef struct { LIS_INT n,calls; } CTX;
static LIS_INT ax(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)raw;
    LIS_INT i;
    c->calls++;
    for(i=0;i<c->n;i++)y[i]=2.0*x[i];
    return LIS_SUCCESS;
}
static int rejected(LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,CTX *ctx,
                    LIS_INT method,LIS_INT precision,LIS_INT use_at)
{
    LIS_SOLVER solver=NULL;
    LIS_INT err,before=ctx->calls;
    if(lis_solver_create(&solver))return 1;
    err=lis_solver_set_option("-i cg -p none -scale none -print none",solver);
    if(!err)
    {
        solver->options[LIS_OPTIONS_SOLVER]=method;
        solver->options[LIS_OPTIONS_PRECISION]=precision;
        solver->options[LIS_OPTIONS_USE_AT]=use_at;
        err=lis_solve(A,b,x,solver);
    }
    lis_solver_destroy(solver);
    if(err!=LIS_ERR_NOT_IMPLEMENTED || ctx->calls!=before)
    {
        fprintf(stderr,"unsupported USER configuration method=%d precision=%d use_at=%d err=%d\n",
                (int)method,(int)precision,(int)use_at,(int)err);
        return 1;
    }
    return 0;
}
static int exercise(void)
{
    CTX c={4,0};
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_INT status,i;
    double residual2=0.0;
    LIS_SOLVER solver=NULL;
    int failed=1;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,c.n) ||
       lis_matrix_set_user(A,&c,ax,NULL) || lis_vector_duplicate(A,&b) ||
       lis_vector_duplicate(A,&x) || lis_vector_set_all(1.0,b) ||
       lis_vector_set_all(0.0,x))goto done;
    if(rejected(A,b,x,&c,LIS_SOLVER_CG,LIS_PRECISION_QUAD,LIS_FALSE) ||
       rejected(A,b,x,&c,LIS_SOLVER_CG,LIS_PRECISION_SWITCH,LIS_FALSE) ||
       rejected(A,b,x,&c,LIS_SOLVER_JACOBI,LIS_PRECISION_DEFAULT,LIS_FALSE) ||
       rejected(A,b,x,&c,LIS_SOLVER_GS,LIS_PRECISION_DEFAULT,LIS_FALSE) ||
       rejected(A,b,x,&c,LIS_SOLVER_SOR,LIS_PRECISION_DEFAULT,LIS_FALSE) ||
       rejected(A,b,x,&c,LIS_SOLVER_BICG,LIS_PRECISION_DEFAULT,LIS_FALSE) ||
       rejected(A,b,x,&c,LIS_SOLVER_CG,LIS_PRECISION_DEFAULT,LIS_TRUE))goto done;
    if(lis_solver_create(&solver) ||
       lis_solver_set_option("-i cg -p none -tol 1e-10 -print none",solver) ||
       lis_solve(A,b,x,solver) || lis_solver_get_status(solver,&status) ||
       status!=LIS_SUCCESS || !c.calls)goto done;
    for(i=0;i<c.n;i++){
        double r=1.0-2.0*(double)x->value[i];
        if(!isfinite(r) || fabs((double)x->value[i]-0.5)>1e-10)goto done;
        residual2+=r*r;
    }
    if(!(sqrt(residual2/c.n)<=1e-10))goto done;
    failed=0;
done:
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}
#endif
int main(int argc,char **argv)
{
#if !defined(USE_MPI) && !defined(_COMPLEX)
    int failed;
    if(lis_initialize(&argc,&argv))return 1;
    failed=exercise();
    lis_finalize();
    printf("USER_SERIAL_CONTRACT %s\n",failed?"FAIL":"PASS");
    return failed?1:0;
#else
    (void)argc;(void)argv;
    puts("USER_SERIAL_CONTRACT SKIP: real serial test");
    return 77;
#endif
}
