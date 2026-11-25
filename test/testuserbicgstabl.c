#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

typedef struct { LIS_INT n,calls,fail_on; } USER_FAULT;
static LIS_INT diagonal(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    USER_FAULT *ctx=(USER_FAULT *)raw;
    LIS_INT i;
    ctx->calls++;
    if(ctx->fail_on>0 && ctx->calls==ctx->fail_on)
    {
        y[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<ctx->n;i++)
        y[i]=(LIS_SCALAR)(4.0+0.125*(double)i)*x[i];
    return LIS_SUCCESS;
}
static int exercise(LIS_INT fail_on,LIS_INT ell)
{
    USER_FAULT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT i,err=LIS_SUCCESS,status=LIS_SUCCESS;
    LIS_REAL resid=NAN;
    int failed=1;
    ctx.n=8;ctx.calls=0;ctx.fail_on=fail_on;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,diagonal,diagonal);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<ctx.n;i++)
    {
        /* Failure must restore x0; the separate numerical fix tests
         * healthy nonzero initial guesses. */
        x->value[i]=(LIS_SCALAR)(fail_on ? 0.15 : 0.0);
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
    }
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(
        "-i bicgstabl -p none -scale none -maxiter 120 -tol 1e-9 -print none",
        solver);
    if(err)goto done;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_FALSE;
    solver->options[LIS_OPTIONS_ELL]=ell;
    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status) || lis_solver_get_residualnorm(solver,&resid))goto done;
    if(fail_on)
    {
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG ||
           ctx.calls!=fail_on || !isnan((double)resid))
        {
            fprintf(stderr,"BiCGSTAB(l=%d) failed callback %d: api=%d status=%d calls=%d\n",
                    (int)ell,(int)fail_on,(int)err,(int)status,(int)ctx.calls);
            goto done;
        }
        for(i=0;i<ctx.n;i++)
            if(!(fabs((double)x->value[i]-0.15)<=1e-14))goto done;
    }
    else
    {
        /* The zero-initial-guess control must converge and return
         * an accurate iterate. Nonzero x0 is tested separately. */
        if(err!=LIS_SUCCESS || status!=LIS_SUCCESS)
        {
            fprintf(stderr,"BiCGSTAB(l=%d) control: api=%d status=%d calls=%d\n",
                    (int)ell,(int)err,(int)status,(int)ctx.calls);
            goto done;
        }
        for(i=0;i<ctx.n;i++)
        {
            const double exact=(1.0+0.125*(double)i)/(4.0+0.125*(double)i);
            if(!(fabs((double)x->value[i]-exact)<=1e-7))
            {
                fprintf(stderr,"BiCGSTAB(l=%d) output mismatch row=%d actual=%.17g exact=%.17g status=%d\n",
                        (int)ell,(int)i,(double)x->value[i],exact,(int)status);
                goto done;
            }
        }
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"USER_BICGSTABL FAIL ell=%d inject=%d api=%d\n",
                (int)ell,(int)fail_on,(int)err);
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
    failed|=exercise(2,1);
    failed|=exercise(3,1);
    failed|=exercise(2,2);
    failed|=exercise(3,2);
    failed|=exercise(4,2);
    failed|=exercise(0,1);
    failed|=exercise(0,2);
    lis_finalize();
    printf("USER_BICGSTABL %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
