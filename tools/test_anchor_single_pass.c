/* Native, serial real-double prototype regression. */
#include "lislib.h"
#include "lis_solver_options.h"
#include "lis_anchor_pilot.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct { LIS_INT stop_at, last, calls; } history;
static LIS_INT stop(LIS_INT iter,LIS_REAL res,void *data)
{
    history *h=data;
    if(iter!=h->last+1 || !isfinite((double)res))return LIS_TRUE;
    h->last=iter;h->calls++;
    return iter==h->stop_at;
}
static LIS_INT matrix(int kind,LIS_INT storage,LIS_MATRIX *A)
{
    LIS_INT i,rc;
    if((rc=lis_matrix_create(LIS_COMM_WORLD,A)) || (rc=lis_matrix_set_size(*A,0,80)))return rc;
    for(i=0;i<80;i++){
        double d=kind==0?1:kind==1?(double)(i+1):pow(10.0,-8.0+8.0*i/79.0);
        if((rc=lis_matrix_set_value(LIS_INS_VALUE,i,i,d,*A)))return rc;
    }
    if((rc=lis_matrix_set_type(*A,storage)))return rc;
    return lis_matrix_assemble(*A);
}
static int kernel_config(LIS_MATRIX P,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_INT maxiter,history *h,LIS_INT *status,LIS_INT *iter,LIS_REAL *res,int zero,LIS_REAL tol,int plain)
{
    LIS_SOLVER s=NULL;LIS_PRECON pc=NULL;LIS_INT rc;
    char options[]="-i fgmres -p ilut -scale none -print none -restart 60 -maxiter 100 -tol 1e-8 -ilutp 1 -ilutp_full 1 -ilutp_full_limit 60000000 -initx_zeros true";
    if((rc=lis_solver_create(&s)) || (rc=lis_solver_set_option(options,s)) ||
       (plain && (rc=lis_solver_set_option("-p none",s))) ||
       (rc=lis_solver_set_maxiter(s,maxiter)) || (rc=lis_solver_set_num_threads(s,1)) ||
       (rc=lis_solver_set_tolerance(s,tol)) ||
       (zero && (rc=lis_vector_set_all(0,x))) || (rc=lis_solver_sequence_prepare_precon(P,s,&pc)))goto done;
    s->options[LIS_OPTIONS_INITGUESS_ZEROS]=zero?LIS_TRUE:LIS_FALSE;
    rc=h?lis_solve_kernel_fgmres_stop(A,b,x,s,pc,stop,h):lis_solve_kernel(A,b,x,s,pc);
    if(!rc)rc=lis_solver_get_status(s,status);
    if(!rc)rc=lis_solver_get_iter(s,iter);
    if(!rc)rc=lis_solver_get_residualnorm(s,res);
done:
    if(s)lis_solver_destroy(s);
    return rc;
}
static int kernel(LIS_MATRIX P,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_INT maxiter,history *h,LIS_INT *status,LIS_INT *iter,LIS_REAL *res,int plain)
{
    return kernel_config(P,A,b,x,maxiter,h,status,iter,res,1,1e-8,plain);
}
static int compare(LIS_MATRIX P,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,LIS_VECTOR y,LIS_INT limit,LIS_INT at,int plain)
{
    LIS_INT a,bstatus,ia,ib;LIS_REAL ra,rb;
    history h={at,0,0};
    if(kernel(P,A,b,x,limit,NULL,&a,&ia,&ra,plain)||kernel(P,A,b,y,100,&h,&bstatus,&ib,&rb,plain))return 1;
    printf("KERNEL at=%lld ordinary_status=%lld stop_status=%lld ordinary_iter=%lld stop_iter=%lld hook_calls=%lld residual=%.15e / %.15e\n",
        (long long)at,(long long)a,(long long)bstatus,(long long)ia,(long long)ib,(long long)h.calls,(double)ra,(double)rb);
    if(a!=bstatus || !isfinite((double)ra) || ra!=rb || memcmp(x->value,y->value,(size_t)x->n*sizeof(LIS_SCALAR)))return 1;
    if(at==20 && (a!=LIS_MAXITER || ib!=20 || h.calls!=20))return 1;
    if(at==0 && (a!=LIS_SUCCESS || ia!=ib || h.calls!=ib))return 1;
    if(at==1 && (a!=LIS_SUCCESS || ia!=ib || h.calls!=1))return 1;
    return 0;
}
int main(int argc,char **argv)
{
    LIS_MATRIX P=NULL,A=NULL,B=NULL,C=NULL;LIS_VECTOR b=NULL,x=NULL,y=NULL;
    lis_anchor_pilot *ctx=NULL;lis_anchor_pilot_report r2;
    int fail=1;
    if(lis_initialize(&argc,&argv))return 1;
    if(matrix(0,LIS_MATRIX_CSR,&P)||matrix(1,LIS_MATRIX_CSR,&A)||
       matrix(2,LIS_MATRIX_CSR,&B)||matrix(2,LIS_MATRIX_CSC,&C)||lis_vector_duplicate(P,&b)||
       lis_vector_duplicate(P,&x)||lis_vector_duplicate(P,&y)||lis_vector_set_all(1,b))goto done;
    /* Stop at 20 must reconstruct exactly the ordinary capped iterate. */
    if(compare(P,B,b,x,y,20,20,0))goto done;
    /* A hook that never stops must be bitwise equivalent to unhooked FGMRES. */
    if(compare(P,A,b,x,y,100,0,0))goto done;
    /* Convergence at the same checkpoint must still be physically certified. */
    if(compare(P,P,b,x,y,100,1,0))goto done;
    /* The internal hook is not tied to the driver's CSR/ILUT policy.
     * A CSC solve with no preconditioner must match an ordinary capped
     * FGMRES solve, while the application driver still rejects CSC. */
    if(C->matrix_type!=LIS_MATRIX_CSC || compare(C,C,b,x,y,20,20,1))goto done;
    if(lis_vector_set_all(7,y) ||
       lis_anchor_fresh_solve(C,b,y,&r2)!=LIS_ERR_ILL_ARG || r2.accepted)goto done;
    for(LIS_INT i=0;i<y->n;i++)if(y->value[i]!=7)goto done;
    if(lis_anchor_pilot_create(C,&ctx)!=LIS_ERR_ILL_ARG || ctx)goto done;
    puts("GENERIC_FGMRES_STOP_CSC_NONE=PASS");
    {
        LIS_INT status,it;LIS_REAL r20,reported,nr,nb;
        LIS_VECTOR residual=NULL;
        if(kernel(P,A,b,x,20,NULL,&status,&it,&r20,0)||status!=LIS_MAXITER||!(r20>1e-8))goto done;
        if(kernel_config(P,A,b,x,80,NULL,&status,&it,&reported,0,1e-8/r20,0)||status!=LIS_SUCCESS)goto done;
        if(lis_vector_duplicate(A,&residual))goto done;
        if(lis_matvec(A,x,residual)||lis_vector_axpy(-1,b,residual)||lis_vector_nrm2(residual,&nr)||lis_vector_nrm2(b,&nb)){
            lis_vector_destroy(residual);goto done;
        }
        lis_vector_destroy(residual);
        printf("RESTART_MATCHED_PHYSICAL 20+%lld residual=%.15e\n",(long long)it,(double)(nr/nb));
        if(!(nr/nb<=1e-8))goto done;
    }
    if(lis_anchor_pilot_create(P,&ctx)||lis_anchor_pilot_solve(ctx,A,b,y,&r2))goto done;
    if(!r2.accepted||!r2.continued||r2.fallback||r2.factorizations||
       r2.pilot_iterations!=20||!(r2.relative_residual<=1e-8))goto done;
    if(lis_anchor_pilot_solve(ctx,B,b,x,&r2)||!r2.accepted||!r2.fallback||r2.continued||r2.pilot_iterations!=20||r2.factorizations!=1)goto done;
    if(lis_anchor_pilot_solve(ctx,P,b,x,&r2)||!r2.accepted||r2.fallback||r2.factorizations||r2.pilot_iterations!=1)goto done;
    fail=0;
done:
    lis_anchor_pilot_destroy(ctx);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(C)lis_matrix_destroy(C);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    if(P)lis_matrix_destroy(P);
    lis_finalize();
    printf("SINGLE_PASS_REGRESSION=%s\n",fail?"FAIL":"PASS");
    return fail;
}
