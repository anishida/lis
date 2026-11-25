#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <math.h>
#include <stdio.h>

#if !defined(USE_MPI) && !defined(_COMPLEX)
typedef struct { LIS_INT n,calls,fail; LIS_SCALAR diag; LIS_INT code; } CTX;
static LIS_INT apply(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)ptr;
    LIS_INT i;
    ++c->calls;
    if(c->fail && c->calls==c->fail)
    {
        y[0]=(LIS_SCALAR)NAN;
        return c->code;
    }
    for(i=0;i<c->n;i++) y[i]=c->diag*x[i];
    return LIS_SUCCESS;
}
static LIS_INT create(CTX *c,LIS_MATRIX *A)
{
    LIS_INT err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err)return err;
    err=lis_matrix_set_size(*A,0,c->n);
    if(err)return err;
    return lis_matrix_set_user(*A,c,apply,apply);
}
static int run(int buckling,LIS_INT afail,LIS_INT bfail,LIS_INT code)
{
    CTX ac={4,0,afail,4.0,code},bc={4,0,bfail,2.0,code};
    LIS_MATRIX A=NULL,B=NULL;
    LIS_VECTOR x=NULL;
    LIS_ESOLVER e=NULL;
    LIS_SCALAR value=0;
    LIS_REAL residual=0;
    LIS_INT err,status,i,iter;
    int failed=1;
    if(create(&ac,&A))goto done;
    if(create(&bc,&B))goto done;
    if(lis_vector_duplicate(A,&x) || lis_esolver_create(&e))goto done;
    if(lis_esolver_set_option(buckling ?
        "-e gii -etransform buckling -shift 1 -etol 1e-10 -emaxiter 50 -eprint none" :
        "-e gii -etransform shift-invert -etol 1e-10 -emaxiter 50 -eprint none",e))goto done;
    if(lis_vector_set_all(1.0,x))goto done;
    err=lis_gesolve(A,B,x,&value,e);
    if(lis_esolver_get_status(e,&status) || lis_esolver_get_iter(e,&iter) ||
       lis_esolver_get_residualnorm(e,&residual))goto done;
    if(afail || bfail)
    {
        if(err!=LIS_SUCCESS || status!=code || iter!=0 ||
           !isnan((double)residual) ||
           (afail && ac.calls!=afail) || (bfail && bc.calls<bfail))
        {
            fprintf(stderr,"fault gen=%d A=%d B=%d api=%d status=%d iter=%d calls=%d/%d resid=%g\n",
                buckling,(int)afail,(int)bfail,(int)err,(int)status,(int)iter,
                (int)ac.calls,(int)bc.calls,(double)residual);
            goto done;
        }
    }
    else if(err!=LIS_SUCCESS || status!=LIS_SUCCESS ||
            !(fabs((double)value-2.0)<=1e-9) ||
            !(residual<=1e-10))goto done;
    if((code==LIS_MAXITER || code==LIS_BREAKDOWN) &&
       (!isnan((double)value) || (afail && ac.calls!=afail) ||
        (bfail && bc.calls!=bfail)))goto done;
    for(i=0;i<x->n;i++)if(!isfinite((double)x->value[i]))goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,"USER eigen inverse failure gen=%d A=%d B=%d\n",
                      buckling,(int)afail,(int)bfail);
    if(e)lis_esolver_destroy(e);
    if(x)lis_vector_destroy(x);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    return failed;
}

static LIS_INT varied_apply(void *ptr,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    CTX *c=(CTX *)ptr;
    LIS_INT i;
    ++c->calls;
    if(c->fail && c->calls==c->fail)return c->code;
    for(i=0;i<c->n;i++)y[i]=(i+1.0)*x[i];
    return LIS_SUCCESS;
}
static int origin(void)
{
    CTX c={4,0,0,0,LIS_MAXITER};
    LIS_CALLBACK_SCOPE scope;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_PRECON pc=NULL;
    LIS_INT status;
    int failed=1,active=0;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,c.n) ||
       lis_matrix_set_user(A,&c,varied_apply,varied_apply) ||
       lis_vector_duplicate(A,&b) || lis_vector_duplicate(A,&x) ||
       lis_vector_set_all(1.0,b) || lis_solver_create(&solver) ||
       lis_solver_set_option("-i cg -p none -maxiter 1 -tol 1e-14",solver))goto done;
    solver->A=A;
    if(lis_precon_create(solver,&pc)) { pc=NULL; goto done; }
    lis_callback_scope_begin(&scope); active=1;
    if(!lis_callback_scope_active()) { failed=0; goto done; }
    status=lis_esolver_inner_solve_attempt(A,b,x,solver,pc);
    if(status!=LIS_MAXITER || lis_callback_scope_status() ||
       !lis_esolver_inner_solve_can_retry(A,NULL,solver,status))goto done;
    lis_callback_scope_end(&scope); active=0;
    c.calls=0; c.fail=1;
    lis_vector_set_all(0.0,x);
    lis_callback_scope_begin(&scope); active=1;
    status=lis_esolver_inner_solve_attempt(A,b,x,solver,pc);
    if(status!=LIS_MAXITER || lis_callback_scope_status()!=status ||
       lis_esolver_inner_solve_can_retry(A,NULL,solver,status) || c.calls!=1)goto done;
    lis_callback_scope_end(&scope); active=0;
    if(lis_callback_scope_active() || lis_callback_scope_status())goto done;
    failed=0;
done:
    if(active)lis_callback_scope_end(&scope);
    if(pc)lis_precon_destroy(pc);
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    if(failed)fprintf(stderr,"USER eigen callback origin regression failed\n");
    return failed;
}
#endif
int main(int argc,char **argv)
{
#if defined(USE_MPI) || defined(_COMPLEX)
    (void)argc;(void)argv;
    puts("USER eigen inverse SKIP: real serial regression");
    return 77;
#else
    int failed=0;
    if(lis_initialize(&argc,&argv))return 1;
    failed|=run(0,0,1,LIS_ERR_ILL_ARG);
    failed|=run(1,1,0,LIS_ERR_ILL_ARG);
    failed|=run(0,0,0,LIS_ERR_ILL_ARG);
    failed|=run(1,0,0,LIS_ERR_ILL_ARG);
    /* A is first applied by the inner shift-invert solve. */
    failed|=run(0,1,0,LIS_MAXITER);
    failed|=run(0,1,0,LIS_BREAKDOWN);
    failed|=origin();
    lis_finalize();
    printf("USER_EIGEN_INVERSE %s\n",failed?"FAILED":"PASSED");
    return failed?1:0;
#endif
}
