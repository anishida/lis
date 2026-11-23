#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

typedef struct { LIS_INT n,calls,fail_on; } CONTEXT;
static LIS_INT apply(void *data,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CONTEXT *ctx=(CONTEXT *)data;
    LIS_INT i;
    ctx->calls++;
    if(ctx->fail_on && ctx->calls==ctx->fail_on) return LIS_ERR_ILL_ARG;
    for(i=0;i<ctx->n;i++) y[i]=(LIS_SCALAR)(4.0+0.2*(double)i)*x[i];
    return LIS_SUCCESS;
}

/* A user-registered preconditioner probes checked psolve failures. */
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
    if(precon_error_at>0 && precon_calls==precon_error_at)
        return LIS_ERR_ILL_ARG;
    return lis_vector_copy(in,out);
}

static int run(const char *method,LIS_INT fail_on,LIS_INT precon_fail_on)
{
    CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=LIS_SUCCESS,status=LIS_SUCCESS,i;
    int failed=1;
    char options[160];
    ctx.n=8; ctx.calls=0; ctx.fail_on=fail_on;
    precon_calls=0; precon_error_at=precon_fail_on;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,apply,apply);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<ctx.n;i++)
    {
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
        x->value[i]=(LIS_SCALAR)0.125;
    }
    err=lis_solver_create(&solver);if(err)goto done;
    snprintf(options,sizeof(options),
        "-i %s -p %s -scale none -maxiter 200 -tol 1e-9 -print none",
        method,precon_fail_on ? "errpc" : "none");
    err=lis_solver_set_option(options,solver);if(err)goto done;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_FALSE;
    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status))goto done;
    if(fail_on || precon_fail_on)
    {
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG ||
           (fail_on && ctx.calls<fail_on) ||
           (precon_fail_on && precon_calls<precon_fail_on))
        {
            fprintf(stderr,"%s fail_on=%d err=%d status=%d calls=%d\n",
                    method,(int)fail_on,(int)err,(int)status,(int)ctx.calls);
            goto done;
        }
    }
    else
    {
        if(err!=LIS_SUCCESS || status!=LIS_SUCCESS)goto done;
        for(i=0;i<ctx.n;i++)
        {
            double exact=(1.0+0.125*(double)i)/(4.0+0.2*(double)i);
            if(fabs((double)x->value[i]-exact)>1e-7)goto done;
        }
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"USER_BICGSAFE_BICRSAFE FAIL method=%s fail_on=%d rc=%d\n",
                method,(int)fail_on,(int)err);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}
int main(int argc,char **argv)
{
#ifdef USE_MPI
    (void)argc; (void)argv;
    puts("testuserbicgsafe SKIP (serial-only USER callback regression)");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    if(lis_precon_register("errpc",error_precon_create,
                          error_precon_apply,error_precon_apply))
    {
        lis_finalize();
        return 1;
    }
    failed|=run("bicgsafe",2,0);
    failed|=run("bicgsafe",3,0);
    failed|=run("bicrsafe",2,0);
    failed|=run("bicrsafe",3,0);
    failed|=run("bicrsafe",4,0);
    failed|=run("bicgsafe",0,0);
    failed|=run("bicrsafe",0,0);
    failed|=run("bicgsafe",0,1);
    failed|=run("bicgsafe",0,2);
    failed|=run("bicrsafe",0,1);
    failed|=run("bicrsafe",0,2);
    lis_precon_register_free();
    lis_finalize();
    printf("USER_BICGSAFE_BICRSAFE %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
