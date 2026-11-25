#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <math.h>
#include <stdio.h>

typedef struct { LIS_INT n,calls,fail_on; } USER_FAULT;
static LIS_INT apply(void *raw,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    USER_FAULT *ctx=(USER_FAULT *)raw;
    LIS_INT i;
    ++ctx->calls;
    if(ctx->fail_on>0 && ctx->calls==ctx->fail_on)
    {
        y[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    for(i=0;i<ctx->n;i++)
        y[i]=(LIS_SCALAR)(4.0+0.2*(double)i)*x[i];
    return LIS_SUCCESS;
}

static LIS_INT pc_calls, pc_fail_on;
static LIS_INT pc_create(LIS_SOLVER solver,LIS_PRECON precon)
{
    (void)solver; (void)precon;
    return LIS_SUCCESS;
}
static LIS_INT pc_apply(LIS_SOLVER solver,LIS_VECTOR in,LIS_VECTOR out)
{
    (void)solver;
    ++pc_calls;
    if(pc_fail_on && pc_calls==pc_fail_on)
    {
        out->value[0]=(LIS_SCALAR)NAN;
        return LIS_ERR_ILL_ARG;
    }
    return lis_vector_copy(in,out);
}

static int run(const char *name,LIS_INT fail_on,LIS_INT precon_fail)
{
    USER_FAULT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,b=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT i,err=LIS_SUCCESS,status=LIS_SUCCESS;
    LIS_INT iterations=-1;
    LIS_REAL resid=NAN;
    int failed=1;
    char opts[200];
    pc_calls=0;pc_fail_on=precon_fail;
    ctx.n=8;ctx.calls=0;ctx.fail_on=fail_on;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,apply,apply);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    for(i=0;i<ctx.n;i++)
    {
        b->value[i]=(LIS_SCALAR)(1.0+0.125*(double)i);
        x->value[i]=(LIS_SCALAR)0.1;
    }
    err=lis_solver_create(&solver);if(err)goto done;
    snprintf(opts,sizeof(opts),"-i %s -p gmreserrpc -scale none -maxiter 200 -tol 1e-9 -print none -restart 8",name);
    if(fail_on==4 || precon_fail==3)
        snprintf(opts,sizeof(opts),"-i %s -p gmreserrpc -scale none -maxiter 200 -tol 1e-9 -print none -restart 2",name);
    err=lis_solver_set_option(opts,solver);
    if(err)goto done;
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS]=LIS_FALSE;
    err=lis_solve(A,b,x,solver);
    if(lis_solver_get_status(solver,&status) || lis_solver_get_iter(solver,&iterations) ||
       lis_solver_get_residualnorm(solver,&resid))goto done;
    if(fail_on || precon_fail)
    {
        if(err!=LIS_SUCCESS || status!=LIS_ERR_ILL_ARG ||
           (fail_on && ctx.calls!=fail_on) ||
           (precon_fail && pc_calls!=precon_fail) ||
           (fail_on==1 && iterations!=0) ||
           ((fail_on==1 || fail_on==4) ? !isnan((double)resid) : !isfinite((double)resid)))
        {
            fprintf(stderr,"GMRES_FGMRES fault %d api=%d status=%d calls=%d\n",
                   (int)fail_on,(int)err,(int)status,(int)ctx.calls);
            goto done;
        }
        for(i=0;i<ctx.n;i++)
            if(!isfinite((double)x->value[i]))goto done;

        /* Reuse the same solver after a callback failure with an exact
         * zero solution. A stale retcode must not turn it into failure. */
        ctx.fail_on=0; pc_fail_on=0;
        if(lis_vector_set_all(0.0,b) || lis_vector_set_all(0.0,x))goto done;
        err=lis_solve(A,b,x,solver);
        if(err || lis_solver_get_status(solver,&status) || status!=LIS_SUCCESS ||
           lis_solver_get_residualnorm(solver,&resid) || resid!=0.0)goto done;
        for(i=0;i<ctx.n;i++) if(x->value[i]!=0.0)goto done;
    }
    else
    {
        if(err!=LIS_SUCCESS || status!=LIS_SUCCESS || !isfinite((double)resid))
        {
            fprintf(stderr,"GMRES_FGMRES healthy api=%d status=%d calls=%d\n",
                    (int)err,(int)status,(int)ctx.calls);
            goto done;
        }
        for(i=0;i<ctx.n;i++)
        {
            double exact=(1.0+0.125*(double)i)/(4.0+0.2*(double)i);
            if(!(fabs((double)x->value[i]-exact)<=1e-7))goto done;
        }
    }
    failed=0;
done:
    if(failed)fprintf(stderr,"USER_GMRES_FGMRES FAIL inject=%d api=%d\n",
                      (int)fail_on,(int)err);
    if(solver)lis_solver_destroy(solver);
    if(b)lis_vector_destroy(b);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}
int main(int argc,char **argv)
{
#ifdef USE_MPI
    (void)argc;(void)argv;
    puts("USER_GMRES_FGMRES SKIP: serial callback test");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    if(lis_precon_register("gmreserrpc",pc_create,pc_apply,pc_apply))
    {
        lis_finalize();
        return 1;
    }
    failed|=run("gmres",1,0);
    failed|=run("gmres",2,0);
    failed|=run("gmres",3,0);
    failed|=run("gmres",4,0);
    failed|=run("gmres",0,3);
    failed|=run("gmres",0,1);
    failed|=run("gmres",0,2);
    failed|=run("gmres",0,0);
    failed|=run("fgmres",1,0);
    failed|=run("fgmres",2,0);
    failed|=run("fgmres",3,0);
    failed|=run("fgmres",4,0);
    failed|=run("fgmres",0,1);
    failed|=run("fgmres",0,2);
    failed|=run("fgmres",0,0);
    lis_precon_register_free();
    lis_finalize();
    printf("USER_GMRES_FGMRES %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
