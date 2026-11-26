#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

typedef struct { LIS_INT n,calls,fail_on; } USER_FAULT;
static LIS_INT diagonal(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    USER_FAULT *ctx=(USER_FAULT *)ptr;
    LIS_INT i;
    ++ctx->calls;
    if(ctx->fail_on>0 && ctx->calls==ctx->fail_on)
    {
        y[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<ctx->n;i++)
        y[i]=(LIS_SCALAR)(4.0+0.25*(double)i)*x[i];
    return LIS_SUCCESS;
}

static int exercise(const char *method,LIS_INT fail_on)
{
    USER_FAULT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT i,err=LIS_SUCCESS,status=LIS_SUCCESS;
    int failed=1;
    char opts[180];
    ctx.n=8;ctx.calls=0;ctx.fail_on=fail_on;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,diagonal,diagonal);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    for(i=0;i<ctx.n;i++)
    {
        x->value[i]=(LIS_SCALAR)0.1;
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
    }
    err=lis_solver_create(&solver);if(err)goto done;
    snprintf(opts,sizeof(opts),
        "-i %s -p none -scale none -maxiter 400 -tol 1e-8 -print none",
        method);
    err=lis_solver_set_option(opts,solver);if(err)goto done;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_FALSE;
    solver->options[LIS_OPTIONS_IDRS_RESTART]=2;
    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status))goto done;
    if(fail_on)
    {
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG ||
           ctx.calls!=fail_on)
        {
            fprintf(stderr,"%s injected=%d api=%d status=%d calls=%d\n",
                    method,(int)fail_on,(int)err,(int)status,(int)ctx.calls);
            goto done;
        }
    }
    else
    {
        if(err!=LIS_SUCCESS || status!=LIS_SUCCESS)goto done;
        for(i=0;i<ctx.n;i++)
        {
            double exact=(1.0+0.125*(double)i)/(4.0+0.25*(double)i);
            if(!(fabs((double)x->value[i]-exact)<=1e-6))goto done;
        }
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"USER_IDR_FAIL method=%s injected=%d api=%d\n",
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
    (void)argc;(void)argv;
    puts("SKIP: serial USER callback test");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=exercise("idr1",2);
    failed|=exercise("idr1",3);
    failed|=exercise("idrs",2);
    failed|=exercise("idrs",3);
    failed|=exercise("idrs",4);
    failed|=exercise("idr1",0);
    failed|=exercise("idrs",0);
    lis_finalize();
    printf("USER_IDR_CALLBACK %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
