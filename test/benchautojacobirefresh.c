/* Optional timing probe: no speedup assertion, single process only.
 * Compare SAFE rebuild vs guarded JACOBI_NUMERIC on the same USER matrix.
 * Reproduce: autoreconf -fiv && mkdir build && cd build &&
 *            ../configure && make -j2 &&
 *            make -C test benchautojacobirefresh &&
 *            ./test/benchautojacobirefresh
 */
#include "lis.h"
#include "lis_solver.h"
#include "lis_matrix_materialize.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#ifdef _COMPLEX
#include <complex.h>
#endif

typedef struct {
    LIS_INT n;
    LIS_REAL scale;
    unsigned long long row_callbacks;
} BENCH_CTX;

static LIS_INT row_callback(void *ctx_,LIS_INT row,LIS_INT *nnz,
                            LIS_INT *cols,LIS_SCALAR *vals)
{
    BENCH_CTX *ctx=(BENCH_CTX *)ctx_;
    LIS_INT count=1,cap,k=0;
    if(ctx==NULL || nnz==NULL || row<0 || row>=ctx->n)
        return LIS_ERR_ILL_ARG;
    ctx->row_callbacks++;
    if(row>0) count++;
    if(row+1<ctx->n) count++;
    cap=*nnz;
    if(cols==NULL && vals==NULL) { *nnz=count; return LIS_SUCCESS; }
    if(cols==NULL || vals==NULL || cap<count) return LIS_ERR_ILL_ARG;
    if(row>0) {cols[k]=row-1; vals[k]=(LIS_SCALAR)(-ctx->scale); k++;}
    cols[k]=row; vals[k]=(LIS_SCALAR)(4*ctx->scale); k++;
    if(row+1<ctx->n) {cols[k]=row+1; vals[k]=(LIS_SCALAR)(-ctx->scale); k++;}
    *nnz=k;
    return LIS_SUCCESS;
}

static LIS_INT matvec_callback(void *ctx_,const LIS_SCALAR *x,LIS_SCALAR *y)
{
    BENCH_CTX *ctx=(BENCH_CTX *)ctx_;
    LIS_INT i,n;
    if(ctx==NULL || x==NULL || y==NULL) return LIS_ERR_ILL_ARG;
    n=ctx->n;
    for(i=0;i<n;i++)
    {
        y[i]=(LIS_SCALAR)(4*ctx->scale)*x[i];
        if(i>0) y[i]-=(LIS_SCALAR)ctx->scale*x[i-1];
        if(i+1<n) y[i]-=(LIS_SCALAR)ctx->scale*x[i+1];
    }
    return LIS_SUCCESS;
}

static int run(BENCH_CTX *ctx,LIS_MATRIX A,LIS_VECTOR b,LIS_VECTOR x,
               LIS_INT refresh_policy,int steps,double *seconds,
               unsigned long long *callback_count)
{
    LIS_MATRIX_AUTO_SESSION session=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err,policy,action,i,t;
    clock_t begin,end;
    int failed=1;
    err=lis_matrix_auto_session_create(&session);
    if(err)goto done;
    err=lis_matrix_auto_session_get_solver(session,&solver);
    if(err)goto done;
    {
        char opt[]="-i gmres -p jacobi -restart 10 -tol 1e-10 -maxiter 80 -scale none";
        err=lis_solver_set_option(opt,solver);
        if(err)goto done;
    }
    err=lis_solver_set_sequence_precon_refresh_policy(solver,refresh_policy);
    if(err)goto done;
    ctx->scale=(LIS_REAL)1;
    for(i=0;i<ctx->n;i++)
        b->value[i]=(LIS_SCALAR)((i==0 || i+1==ctx->n ? 3 : 2)*ctx->scale);
    err=lis_matrix_mark_values_changed(A); if(err)goto done;
    err=lis_vector_set_all((LIS_SCALAR)0,x); if(err)goto done;
    err=lis_matrix_auto_session_solve(session,A,b,x,&policy,&action);
    if(err || policy!=LIS_MATRIX_POLICY_CACHED_CSR)goto done;
    ctx->row_callbacks=0;
    begin=clock();
    for(t=0;t<steps;t++)
    {
        ctx->scale=(LIS_REAL)(1.0+0.02*(t+1));
        err=lis_matrix_mark_values_changed(A); if(err)goto done;
        for(i=0;i<ctx->n;i++)
            b->value[i]=(LIS_SCALAR)((i==0 || i+1==ctx->n ? 3 : 2)*ctx->scale);
        err=lis_vector_set_all((LIS_SCALAR)0,x); if(err)goto done;
        err=lis_matrix_auto_session_solve(session,A,b,x,&policy,&action);
        if(err || policy!=LIS_MATRIX_POLICY_CACHED_CSR ||
           action!=LIS_MATRIX_MATERIALIZE_VALUES)goto done;
        for(i=0;i<ctx->n;i++)
        {
#ifdef _COMPLEX
            if(cabs(x->value[i]-(LIS_SCALAR)1)>1.e-6)goto done;
#else
            if(fabs((double)(x->value[i]-(LIS_SCALAR)1))>1.e-6)goto done;
#endif
        }
    }
    end=clock();
    *seconds=(double)(end-begin)/CLOCKS_PER_SEC;
    *callback_count=ctx->row_callbacks;
    failed=0;
done:
    if(session)lis_matrix_auto_session_destroy(session);
    return failed;
}

int main(int argc,char **argv)
{
    const LIS_INT n=3000;
    const int steps=20;
    BENCH_CTX ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    double baseline=0,optimized=0;
    unsigned long long base_calls=0,opt_calls=0;
    LIS_INT err;
    int failed=1;
    err=lis_initialize(&argc,&argv);if(err)return 1;
    ctx.n=n;ctx.scale=(LIS_REAL)1;ctx.row_callbacks=0;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,n);if(err || A->n!=n)goto done;
    err=lis_matrix_set_user(A,&ctx,matvec_callback,matvec_callback);if(err)goto done;
    err=lis_matrix_set_user_row(A,row_callback);if(err)goto done;
    err=lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);if(err)goto done;
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    if(run(&ctx,A,b,x,LIS_SEQUENCE_PRECON_REFRESH_SAFE,
           steps,&baseline,&base_calls))goto done;
    if(run(&ctx,A,b,x,LIS_SEQUENCE_PRECON_REFRESH_JACOBI_NUMERIC,
           steps,&optimized,&opt_calls))goto done;
    printf("AUTO_JACOBI_BENCH n=%lld steps=%d safe_cpu_s=%.6f refresh_cpu_s=%.6f safe_row_callbacks=%llu refresh_row_callbacks=%llu\n",
           (long long)n,steps,baseline,optimized,base_calls,opt_calls);
    /* The verified fast path must never cause an extra callback scan.
     * SAFE performs the normal pattern verification and value refresh.
     */
    if(opt_calls!=base_calls)
    {
        fprintf(stderr,"Jacobi performed an extra CSR row scan (%llu vs %llu)\n",
                opt_calls,base_calls);
        goto done;
    }
    /* Timing varies across machines; do not impose a success threshold. */
    failed=0;
done:
    printf("AUTO_JACOBI_BENCH=%s\n",failed?"FAIL":"PASS");
    if(b)lis_vector_destroy(b);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    return failed;
}
