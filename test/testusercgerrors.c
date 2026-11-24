/*
 * CG must propagate failed matrix-free USER matvec/preconditioner calls
 * without consuming their partially written outputs. The QUAD variant
 * is tested against a native CSR matrix: USER matrices are intentionally
 * restricted to default precision by LIS solver policy.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

typedef struct {
    LIS_INT n,calls,fail_on;
} CG_USER;

static LIS_INT cg_user_matvec(void *data,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CG_USER *ctx=(CG_USER *)data;
    LIS_INT i;
    ++ctx->calls;
    if(ctx->fail_on>0 && ctx->calls==ctx->fail_on)
    {
        y[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<ctx->n;i++)
        y[i]=(LIS_SCALAR)(3.0+0.25*(double)i)*x[i];
    return LIS_SUCCESS;
}

static LIS_INT precon_calls=0,precon_fail_on=0;
static LIS_INT cg_precon_create(LIS_SOLVER solver,LIS_PRECON precon)
{
    (void)solver;(void)precon;
    return LIS_SUCCESS;
}
static LIS_INT cg_precon_apply(LIS_SOLVER solver,LIS_VECTOR from,LIS_VECTOR to)
{
    (void)solver;
    ++precon_calls;
    if(precon_fail_on>0 && precon_calls==precon_fail_on)
    {
        to->value[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    return lis_vector_copy(from,to);
}

static int run(LIS_INT use_quad,LIS_INT matvec_fail_on,LIS_INT pc_fail_on)
{
    CG_USER ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT i,err=LIS_SUCCESS,status=0,completed=-1;
    LIS_REAL residual=NAN;
    int failed=1;
    char opt[]="-i cg -p cgerrpc -initx_zeros false -scale none -maxiter 100 -tol 1e-10 -print none";

    ctx.n=8;ctx.calls=0;ctx.fail_on=matvec_fail_on;
    precon_calls=0;precon_fail_on=pc_fail_on;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    if(use_quad)
    {
        for(i=0;i<ctx.n;i++)
        {
            err=lis_matrix_set_value(LIS_INS_VALUE,i,i,
                                     (LIS_SCALAR)(3.0+0.25*(double)i),A);
            if(err)goto done;
        }
        err=lis_matrix_set_type(A,LIS_MATRIX_CSR);if(err)goto done;
        err=lis_matrix_assemble(A);if(err)goto done;
    }
    else
    {
        err=lis_matrix_set_user(A,&ctx,cg_user_matvec,cg_user_matvec);
        if(err)goto done;
    }
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<ctx.n;i++)
    {
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
        x->value[i]=(LIS_SCALAR)0.1;
    }
    err=lis_solver_create(&solver);if(err)goto done;
    err=lis_solver_set_option(opt,solver);if(err)goto done;
#ifdef USE_QUAD_PRECISION
    if(use_quad)
        solver->options[LIS_OPTIONS_PRECISION]=LIS_PRECISION_QUAD;
#else
    if(use_quad)goto done;
#endif
    /* As with other LIS callback regressions, lis_solve() uses the
     * iterative retcode in solver status and returns API-level success. */
    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status) ||
       lis_solver_get_iter(solver,&completed) ||
       lis_solver_get_residualnorm(solver,&residual))goto done;
    if(matvec_fail_on || pc_fail_on)
    {
        LIS_INT expect_completed;
        if(use_quad)
            expect_completed=pc_fail_on-1;
        else if(matvec_fail_on)
            expect_completed=matvec_fail_on==3?1:0;
        else
            expect_completed=pc_fail_on-1;
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG ||
           completed!=expect_completed ||
           (matvec_fail_on && ctx.calls!=matvec_fail_on) ||
           (!use_quad && matvec_fail_on && precon_calls!=matvec_fail_on-1) ||
           (pc_fail_on && precon_calls!=pc_fail_on) ||
           (!use_quad && pc_fail_on && ctx.calls!=pc_fail_on))
        {
            fprintf(stderr,
                "CG %s error propagation: matvec_at=%lld pc_at=%lld api=%lld status=%lld iter=%lld calls=%lld pcalls=%lld\n",
                use_quad?"QUAD":"USER",(long long)matvec_fail_on,
                (long long)pc_fail_on,(long long)err,(long long)status,
                (long long)completed,(long long)ctx.calls,
                (long long)precon_calls);
            goto done;
        }
        for(i=0;i<ctx.n;i++)
        {
            double val=(double)x->value[i];
            if(!isfinite(val) ||
               (expect_completed==0 && fabs(val-0.1)>1e-15))
                goto done;
        }
    }
    else
    {
        if(err!=LIS_SUCCESS || status!=LIS_SUCCESS ||
           completed<1 || !isfinite((double)residual) ||
           residual>1e-10)goto done;
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
            "USER_CG_ERRORS FAIL quad=%lld matvec_at=%lld pc_at=%lld err=%lld status=%lld iter=%lld\n",
            (long long)use_quad,(long long)matvec_fail_on,
            (long long)pc_fail_on,(long long)err,(long long)status,
            (long long)completed);
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
    puts("USER_CG_ERRORS SKIP: serial callback regression");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    if(lis_precon_register("cgerrpc",cg_precon_create,
                          cg_precon_apply,cg_precon_apply))
    {
        lis_finalize();
        return 1;
    }
    failed|=run(0,1,0); /* initial residual helper (existing behavior) */
    failed|=run(0,2,0); /* in-loop USER matvec */
    failed|=run(0,3,0); /* later USER matvec */
    failed|=run(0,0,1); /* first preconditioner application */
    failed|=run(0,0,2); /* later preconditioner application */
    failed|=run(0,0,0); /* genuine convergence */
#ifdef USE_QUAD_PRECISION
    failed|=run(1,0,1); /* QUAD registered-preconditioner failure */
    failed|=run(1,0,2); /* later QUAD failure */
    failed|=run(1,0,0); /* QUAD convergence on native CSR */
#endif
    lis_precon_register_free();
    lis_finalize();
    printf("USER_CG_ERRORS %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
