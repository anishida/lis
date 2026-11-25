/* USER-row callback microbenchmark: light/heavy callbacks, local-row order.
 * No speedup assertion: print timing and verify equal numerical checksums.
 * Intended as a reproducible complement to real FEA setup+solve traces.
 */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>

typedef struct {
    LIS_INT n;
    unsigned long scalar_calls,batch_calls;
    int heavy;
    volatile double sink;
} BENCH_CONTEXT;

static void simulate_dispatch_cost(BENCH_CONTEXT *c)
{
    int j;
    volatile double t=0.0;
    if(c->heavy)
    {
        /* Model expensive C/Fortran/FFI invocation setup, not row work.
         * Keep the loop observable to avoid compiler dead-code removal. */
        for(j=0;j<700;j++) t+=(double)(j+1)*0.00000001;
        c->sink+=t;
    }
}

static LIS_INT row_size(LIS_INT n,LIS_INT row)
{
    return (row>0)+1+(row+1<n);
}

static LIS_INT scalar_cb(void *opaque,LIS_INT row,LIS_INT *nnz,
                         LIS_INT *cols,LIS_SCALAR *values)
{
    BENCH_CONTEXT *c=(BENCH_CONTEXT *)opaque;
    LIS_INT k=0,needed;
    if(!c || row<0 || row>=c->n || nnz==NULL)return LIS_ERR_ILL_ARG;
    c->scalar_calls++;
    simulate_dispatch_cost(c);
    needed=row_size(c->n,row);
    if(!cols && !values){*nnz=needed;return LIS_SUCCESS;}
    if(!cols || !values || *nnz<needed)return LIS_ERR_ILL_ARG;
    if(row>0){cols[k]=row-1;values[k++]=(LIS_SCALAR)-1.0;}
    cols[k]=row;values[k++]=(LIS_SCALAR)4.0;
    if(row+1<c->n){cols[k]=row+1;values[k++]=(LIS_SCALAR)-1.0;}
    *nnz=needed;
    return LIS_SUCCESS;
}

static LIS_INT batch_cb(void *opaque,LIS_INT first,LIS_INT last,
                        LIS_INT *offsets,LIS_INT capacity,
                        LIS_INT *cols,LIS_SCALAR *values)
{
    BENCH_CONTEXT *c=(BENCH_CONTEXT *)opaque;
    LIS_INT row,k=0;
    if(!c || !offsets || first<0 || first>=last || last>c->n)
        return LIS_ERR_ILL_ARG;
    c->batch_calls++;
    simulate_dispatch_cost(c);
    offsets[0]=0;
    for(row=first;row<last;row++)
        offsets[row-first+1]=offsets[row-first]+row_size(c->n,row);
    if(!cols && !values)return LIS_SUCCESS;
    if(!cols || !values || capacity<offsets[last-first])
        return LIS_ERR_ILL_ARG;
    for(row=first;row<last;row++)
    {
        if(row>0){cols[k]=row-1;values[k++]=(LIS_SCALAR)-1.0;}
        cols[k]=row;values[k++]=(LIS_SCALAR)4.0;
        if(row+1<c->n){cols[k]=row+1;values[k++]=(LIS_SCALAR)-1.0;}
    }
    return LIS_SUCCESS;
}

static int run_one(int use_batch,int heavy,int random_access,
                   double *time_sec,double *sum_out,
                   unsigned long *dispatches)
{
    enum { N=2048, REPEATS=8 };
    BENCH_CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_PRECON_ROW_SOURCE src;
    LIS_INT err=0,i,repeat,k;
    clock_t start,stop;
    double checksum=0.0;
    int bad=1,initialized=0;
    memset(&ctx,0,sizeof(ctx));
    memset(&src,0,sizeof(src));
    ctx.n=N;ctx.heavy=heavy;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,N);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,scalar_cb);if(err)goto done;
    if(use_batch)
    {
        err=lis_matrix_set_user_rows_batch(A,batch_cb);
        if(err)goto done;
    }
    err=lis_precon_row_source_init(A,&src);if(err)goto done;
    initialized=1;
    start=clock();
    for(repeat=0;repeat<REPEATS;repeat++)
    {
        for(i=0;i<N;i++)
        {
            /* 1031 is coprime with 2048: full-cycle strided access. */
            LIS_INT row=random_access ? (i*1031)%N : i;
            err=lis_precon_row_source_get(&src,row);
            if(err)goto done;
            for(k=0;k<src.nnz;k++)
                checksum+=(double)src.value[k]*(double)(src.index[k]+1);
        }
    }
    stop=clock();
    *time_sec=(double)(stop-start)/CLOCKS_PER_SEC;
    *sum_out=checksum;
    *dispatches=ctx.scalar_calls+ctx.batch_calls;
    bad=0;
done:
    if(bad)fprintf(stderr,
        "USER_BATCH_BENCH FAIL batch=%d heavy=%d random=%d err=%d\n",
        use_batch,heavy,random_access,(int)err);
    if(initialized)lis_precon_row_source_destroy(&src);
    if(A)lis_matrix_destroy(A);
    return bad;
}

/* Benchmark the actual generic USER SpMV path, rather than only direct
 * row-source iteration. Repeats are fixed; timing is informational.
 */
static int run_spmv(int use_batch,int heavy,double *elapsed,
                    double *checksum,unsigned long *dispatches)
{
    enum { N=2048, REPEATS=16 };
    BENCH_CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_INT i,r,err=0;
    clock_t start,stop;
    double sum=0.0;
    int failed=1;
    memset(&ctx,0,sizeof(ctx));
    ctx.n=N;ctx.heavy=heavy;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,N);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,scalar_cb);if(err)goto done;
    if(use_batch)
    {
        err=lis_matrix_set_user_rows_batch(A,batch_cb);
        if(err)goto done;
    }
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&y);if(err)goto done;
    for(i=0;i<N;i++) x->value[i]=(LIS_SCALAR)(1.0+0.0001*i);
    start=clock();
    for(r=0;r<REPEATS;r++)
    {
        err=lis_matvec(A,x,y);
        if(err)goto done;
        for(i=0;i<N;i++)sum+=(double)y->value[i];
    }
    stop=clock();
    *elapsed=(double)(stop-start)/CLOCKS_PER_SEC;
    *checksum=sum;
    *dispatches=ctx.scalar_calls+ctx.batch_calls;
    failed=0;
done:
    if(failed)fprintf(stderr,
        "USER_BATCH_SPMV FAIL batch=%d heavy=%d err=%d\n",
        use_batch,heavy,(int)err);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* End-to-end solver microbenchmark. Unlike row-source measurements,
 * this exercises lis_solve's real USER matvec/CG path. No preconditioner
 * is used; true residuals are computed independently of the callback.
 */
static int run_cg_solve(int use_batch,int heavy,double *setup_sec,
                        double *solve_sec,LIS_INT *iterations,
                        double *true_relative,unsigned long *dispatches)
{
    enum { N=512, REPEATS=4 };
    BENCH_CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_SOLVER solver=NULL;
    LIS_INT err=0,status=0,iter=0,i,r;
    clock_t begin,end;
    double configure=0.0,solving=0.0,rsum=0.0,bsum=0.0;
    int failed=1;
    memset(&ctx,0,sizeof(ctx));
    ctx.n=N;ctx.heavy=heavy;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,N);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,scalar_cb);if(err)goto done;
    if(use_batch)
    {
        err=lis_matrix_set_user_rows_batch(A,batch_cb);
        if(err)goto done;
    }
    err=lis_vector_duplicate(A,&b);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    for(i=0;i<N;i++)
    {
        double exact_i=1.0+0.001*i;
        double rhs=4.0*exact_i;
        if(i>0) rhs-=(1.0+0.001*(i-1));
        if(i+1<N) rhs-=(1.0+0.001*(i+1));
        b->value[i]=(LIS_SCALAR)rhs;
        bsum+=rhs*rhs;
    }
    for(r=0;r<REPEATS;r++)
    {
        err=lis_vector_set_all((LIS_SCALAR)0.0,x);
        if(err)goto done;
        begin=clock();
        err=lis_solver_create(&solver);
        if(err)goto done;
        err=lis_solver_set_option(
            "-i cg -p none -scale none -tol 1e-11 "
            "-maxiter 150 -print none",solver);
        if(err)goto done;
        end=clock();
        configure+=(double)(end-begin)/CLOCKS_PER_SEC;
        begin=clock();
        err=lis_solve(A,b,x,solver);
        end=clock();
        solving+=(double)(end-begin)/CLOCKS_PER_SEC;
        if(err)goto done;
        err=lis_solver_get_status(solver,&status);
        if(err || status!=LIS_SUCCESS)goto done;
        err=lis_solver_get_iter(solver,&iter);
        if(err)goto done;
        lis_solver_destroy(solver);
        solver=NULL;
    }
    for(i=0;i<N;i++)
    {
        double ax=4.0*(double)x->value[i];
        if(i>0) ax-=(double)x->value[i-1];
        if(i+1<N) ax-=(double)x->value[i+1];
        ax-=(double)b->value[i];
        rsum+=ax*ax;
    }
    *true_relative=sqrt(rsum/bsum);
    if(!isfinite(*true_relative) || *true_relative>1e-9)
        goto done;
    *setup_sec=configure/REPEATS;
    *solve_sec=solving/REPEATS;
    *iterations=iter;
    *dispatches=ctx.scalar_calls+ctx.batch_calls;
    failed=0;
done:
    if(failed)
        fprintf(stderr,"USER_BATCH_SOLVE FAIL batch=%d heavy=%d "
                "err=%d status=%d true_rel=%.6e\n",
                use_batch,heavy,(int)err,(int)status,
                sqrt(bsum>0?rsum/bsum:0.0));
    if(solver)lis_solver_destroy(solver);
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    return failed;
}

int main(int argc,char **argv)
{
    int heavy,random_access,mode;
    int failed=0;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS)return 2;
    for(heavy=0;heavy<=1;heavy++)
    for(random_access=0;random_access<=1;random_access++)
    {
        double reference=0.0;
        for(mode=0;mode<=1;mode++)
        {
            double seconds=0.0,checksum=0.0;
            unsigned long calls=0;
            if(run_one(mode,heavy,random_access,&seconds,&checksum,&calls))
            {
                failed=1;
                goto done;
            }
            if(mode==0)reference=checksum;
            if(fabs(checksum-reference)>1e-10*(1.0+fabs(reference)))
            {
                failed=1;
                fprintf(stderr,"USER_BATCH_BENCH checksum mismatch\n");
                goto done;
            }
            printf("USER_BATCH_BENCH callback=%s rows=%s mode=%s "
                   "seconds=%.6f dispatches=%lu checksum=%.6f\n",
                   heavy?"heavy":"light",
                   random_access?"strided":"sequential",
                   mode?"batch":"scalar",seconds,calls,checksum);
        }
    }
    /* SpMV is the user-facing operation that benefits from this callback.
     * Compare both callback variants on the identical matrix and vectors. */
    for(heavy=0;heavy<=1;heavy++)
    {
        double reference=0.0;
        for(mode=0;mode<=1;mode++)
        {
            double seconds=0.0,checksum=0.0;
            unsigned long calls=0;
            if(run_spmv(mode,heavy,&seconds,&checksum,&calls))
            {
                failed=1;
                goto done;
            }
            if(mode==0)reference=checksum;
            if(fabs(checksum-reference)>1e-10*(1.0+fabs(reference)))
            {
                failed=1;
                fprintf(stderr,"USER_BATCH_SPMV checksum mismatch\n");
                goto done;
            }
            printf("USER_BATCH_SPMV callback=%s mode=%s "
                   "seconds=%.6f dispatches=%lu checksum=%.6f\n",
                   heavy?"heavy":"light",mode?"batch":"scalar",
                   seconds,calls,checksum);
        }
    }
    /* Real solver API: configuration and four full CG solves, with
     * independent analytical true-residual validation. */
    for(heavy=0;heavy<=1;heavy++)
    for(mode=0;mode<=1;mode++)
    {
        double setup=0.0,solve=0.0,true_rel=0.0;
        unsigned long calls=0;
        LIS_INT iter=0;
        if(run_cg_solve(mode,heavy,&setup,&solve,
                        &iter,&true_rel,&calls))
        {
            failed=1;
            goto done;
        }
        printf("USER_BATCH_SOLVE callback=%s mode=%s "
               "configure_s=%.6f solve_s=%.6f iterations=%d "
               "true_relative=%.3e dispatches=%lu\n",
               heavy?"heavy":"light",mode?"batch":"scalar",
               setup,solve,(int)iter,true_rel,calls);
    }
done:
    lis_finalize();
    return failed;
}
