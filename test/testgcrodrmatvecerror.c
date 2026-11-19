/* GCRO-DR must stop when a USER matvec callback fails in Arnoldi. */
#include "lis.h"
#include <math.h>
#include <stdio.h>

typedef struct {
    LIS_INT n,calls,fail_on;
} ERROR_OPERATOR;

static LIS_INT
matvec_fail_on_call(void *vctx,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    ERROR_OPERATOR *ctx=(ERROR_OPERATOR *)vctx;
    LIS_INT i;
    if(!ctx || !x || !y) return LIS_ERR_ILL_ARG;
    ctx->calls++;
    if(ctx->fail_on>0 && ctx->calls>=ctx->fail_on)
        return LIS_ERR_ILL_ARG;
    for(i=0;i<ctx->n;i++)
    {
        y[i]=4.0*x[i];
        if(i>0) y[i]-=x[i-1];
        if(i+1<ctx->n) y[i]-=x[i+1];
    }
    return LIS_SUCCESS;
}

int main(int argc,char **argv)
{
    enum { N=6 };
    ERROR_OPERATOR ctx={N,0,0};
    LIS_MATRIX A=NULL;
    LIS_SOLVER solver=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_INT i,err,status;
    int failed=1;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS) return 2;
    err=lis_matrix_create(LIS_COMM_WORLD,&A); if(err)goto done;
    err=lis_matrix_set_size(A,0,N); if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,matvec_fail_on_call,
                              matvec_fail_on_call);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<N;i++)
    {
        LIS_SCALAR exact=(LIS_SCALAR)(1+i);
        b->value[i]=4.0*exact;
        if(i>0)b->value[i]-=(LIS_SCALAR)i;
        if(i+1<N)b->value[i]-=(LIS_SCALAR)(i+2);
    }
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(
        "-i gcrodr -p none -scale none -restart 6 -gcrodr_k 1 "
        "-maxiter 50 -tol 1e-11 -print none",solver);
    if(err)goto done;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);if(err)goto done;
    err=lis_solve(A,b,x,solver);
    if(err || ctx.calls<1)goto done;
    for(i=0;i<N;i++)
    {
        if(fabs((double)x->value[i]-(double)(i+1))>1e-7)goto done;
    }
    /* Fail only AFTER at least one successful matrix-vector product;
     * previous code silently continued with stale Arnoldi vectors. */
    ctx.calls=0;
    ctx.fail_on=2;
    err=lis_vector_set_all((LIS_SCALAR)0.0,x);if(err)goto done;
    err=lis_solve(A,b,x,solver);
    if(err)goto done;
    err=lis_solver_get_status(solver,&status);
    if(err || status!=LIS_ERR_ILL_ARG || ctx.calls<2)goto done;
    failed=0;
done:
    if(failed)
        fprintf(stderr,"GCRODR_USER_MATVEC_ERROR FAIL err=%d calls=%d\n",
                (int)err,(int)ctx.calls);
    else
        puts("GCRODR_USER_MATVEC_ERROR PASS");
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
