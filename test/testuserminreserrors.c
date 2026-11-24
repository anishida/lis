/*
 * MINRES must honor matrix-free USER matvec and user preconditioner errors
 * on both the initial residual and subsequent Lanczos iterations.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

typedef struct {
    LIS_INT n, calls, fail_on;
} USER_CONTEXT;

static LIS_INT callback(void *data,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    USER_CONTEXT *ctx=(USER_CONTEXT *)data;
    LIS_INT i;
    ++ctx->calls;
    if(ctx->fail_on && ctx->calls==ctx->fail_on)
    {
        /* Simulate a callback writing only part of its result before
         * failing: consuming it would cause a spurious convergence status. */
        y[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<ctx->n;i++)
        y[i]=(LIS_SCALAR)(3.0+0.25*(double)i)*x[i];
    return LIS_SUCCESS;
}

static LIS_INT precon_calls=0,precon_error_at=0;

static LIS_INT error_precon_create(LIS_SOLVER solver,LIS_PRECON precon)
{
    (void)solver; (void)precon;
    return LIS_SUCCESS;
}

static LIS_INT error_precon_apply(LIS_SOLVER solver,LIS_VECTOR in,LIS_VECTOR out)
{
    (void)solver;
    ++precon_calls;
    if(precon_error_at && precon_calls==precon_error_at)
    {
        out->value[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    return lis_vector_copy(in,out);
}

static int run(LIS_INT fail_on,LIS_INT precon_fail_on)
{
    USER_CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT i,err=LIS_SUCCESS,status=LIS_SUCCESS;
    int failed=1;
    char opts[]="-i minres -p minerrpc -scale none -maxiter 100 -tol 1e-10 -print none";
    ctx.n=8;ctx.calls=0;ctx.fail_on=fail_on;
    precon_calls=0;precon_error_at=precon_fail_on;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,callback,callback);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<ctx.n;i++)
    {
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
        x->value[i]=(LIS_SCALAR)0.1;
    }
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(opts,solver);if(err)goto done;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_FALSE;

    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status))goto done;
    if(fail_on || precon_fail_on)
    {
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG ||
           (fail_on && ctx.calls<fail_on) ||
           (precon_fail_on && precon_calls<precon_fail_on))
        {
            fprintf(stderr,
                "MINRES USER error not propagated: matvec_at=%lld psolve_at=%lld ret=%lld status=%lld calls=%lld pcalls=%lld iter=%lld\n",
                (long long)fail_on,(long long)precon_fail_on,
                (long long)err,(long long)status,
                (long long)ctx.calls,(long long)precon_calls,
                (long long)solver->iter);
            goto done;
        }
    }
    else
    {
        if(err!=LIS_SUCCESS || status!=LIS_SUCCESS)goto done;
        for(i=0;i<ctx.n;i++)
        {
            double exact=(1.0+0.125*(double)i)/(3.0+0.25*(double)i);
            if(!isfinite((double)x->value[i]) ||
               fabs((double)x->value[i]-exact)>1e-7)goto done;
        }
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,
            "USER_MINRES_ERRORS FAIL matvec_at=%lld psolve_at=%lld status=%lld api=%lld\n",
            (long long)fail_on,(long long)precon_fail_on,
            (long long)status,(long long)err);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}

int main(int argc,char **argv)
{
#ifdef USE_MPI
    (void)argc;(void)argv;
    puts("USER_MINRES_ERRORS SKIP: serial callback test");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    if(lis_precon_register("minerrpc",error_precon_create,
                          error_precon_apply,error_precon_apply))
    {
        lis_finalize();
        return 1;
    }

    /* Initial and in-loop matrix-vector callback failures. */
    failed|=run(1,0);
    failed|=run(2,0);
    failed|=run(3,0);
    /* Initial and in-loop preconditioner failures. */
    failed|=run(0,1);
    failed|=run(0,2);
    /* No fault: verify the SPD diagonal system is actually solved. */
    failed|=run(0,0);
    lis_precon_register_free();
    lis_finalize();
    printf("USER_MINRES_ERRORS %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
