/* Extended Krylov API contract and allocation-limit regression.
 * Build under a 32-bit target to test pre-multiplication bounds directly.
 * No large buffer is allocated by this test.
 */
#include "lis_krylov.h"
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

static int apply_diagonal(void *ctx, int n, const double *x, double *y)
{
    int i;
    (void)ctx;
    for (i = 0; i < n; ++i) y[i] = (double)(i + 1) * x[i];
    return 0;
}

static int apply_diagonal64(void *ctx,int64_t n,const double *x,double *y)
{
    int64_t i;
    (void)ctx;
    for(i=0;i<n;++i)y[i]=(double)(i+1)*x[i];
    return 0;
}

/* Synthetic inertia metadata exercises >32-bit counts without allocation. */
static int inertia64(void *ctx,int64_t n,double shift,int64_t *negative,int64_t *zero)
{
    (void)ctx;
    *negative=shift>0?n:0;
    *zero=0;
    return 0;
}

static int identity64(void *ctx,int64_t n,const double *x,double *y)
{
    (void)ctx;
    memcpy(y,x,(size_t)n*sizeof(double));
    return 0;
}
static int diagonal_inverse64(void *ctx,int64_t n,const double *x,double *y)
{
    int64_t i;
    const double sigma=*(const double*)ctx;
    for(i=0;i<n;++i)y[i]=x[i]/((double)(i+1)-sigma);
    return 0;
}
static int diagonal_inertia64(void *ctx,int64_t n,double shift,
                             int64_t *negative,int64_t *zero)
{
    int64_t i;
    (void)ctx;
    *negative=*zero=0;
    for(i=0;i<n;++i){
        double d=(double)(i+1)-shift;
        if(d<0)++*negative;
        else if(d==0)++*zero;
    }
    return 0;
}

static int fail_callback(void *ctx,int n,const double *x,double *y)
{
    (void)ctx;(void)n;(void)x;(void)y;
    return 17;
}

int main(void)
{
    lis_krylov_problem p = {0};
    lis_krylov_options options;
    lis_krylov_stats stats;
    lis_krylov_ext_options ext;
    lis_krylov_ext_stats extra;
    lis_krylov_workspace *w = NULL;
    double values[2], vectors[16], residuals[2];
    int rc;

    lis_krylov_workspace_destroy(NULL);
    lis_krylov_options_default(&options);
    lis_krylov_ext_options_default(&ext);
    lis_krylov_ext_stats_default(&extra);
    p.n = 8;
    p.apply_a = apply_diagonal;
    options.nev = 2;
    options.ncv = 5;
    options.tol = 1e-8;
    options.maxiter = 100;

    /* The output stats structure has an independent version/size
     * check. Invalid memory layouts must be rejected before use. */
    extra.version = LIS_KRYLOV_EXT_VERSION + 1u;
    rc = lis_krylov_symmetric_ex(&p, &options, &ext, NULL,
                                 values, vectors, residuals, &stats, &extra);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) {
        fprintf(stderr, "FAIL invalid statistics version: %d\n", rc);
        return 1;
    }
    lis_krylov_ext_stats_default(&extra);
    extra.struct_size--;
    rc = lis_krylov_symmetric_ex(&p, &options, &ext, NULL,
                                 values, vectors, residuals, &stats, &extra);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) {
        fprintf(stderr, "FAIL invalid statistics size: %d\n", rc);
        return 1;
    }
    lis_krylov_ext_stats_default(&extra);
    ext.version++;
    w = (lis_krylov_workspace *)(uintptr_t)1;
    rc = lis_krylov_workspace_create(&p, &options, &ext, &w);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT || w != NULL) {
        fprintf(stderr, "FAIL invalid workspace version: %d\n", rc);
        return 1;
    }
    lis_krylov_ext_options_default(&ext);
    ext.struct_size++;
    rc = lis_krylov_symmetric_ex(&p, &options, &ext, NULL,
                                 values, vectors, residuals, &stats, &extra);
    if (rc != LIS_KRYLOV_BAD_ARGUMENT) {
        fprintf(stderr, "FAIL invalid options size: %d\n", rc);
        return 1;
    }

    {
        lis_krylov_budget budget;
        lis_krylov_budget_stats used;
        int idx[5],count=-1;
        double v[5]={2,3,NAN,4,5},r[5]={0,1e-3,0,-1,1e-9};
        lis_krylov_ext_options_default(&ext);
        lis_krylov_budget_default(&budget);
        budget.max_apply_a=1;
        rc=lis_krylov_symmetric_budgeted(&p,&options,&ext,NULL,&budget,
                    values,vectors,residuals,&stats,&extra,&used);
        if(rc!=LIS_KRYLOV_BUDGET_EXHAUSTED || !used.exhausted ||
           used.apply_a!=1 || used.total!=1 || stats.nreturned || stats.nconv)
            return 1;
        p.apply_a=fail_callback;
        rc=lis_krylov_symmetric_budgeted(&p,&options,&ext,NULL,&budget,
                    values,vectors,residuals,&stats,&extra,&used);
        if(rc!=LIS_KRYLOV_CALLBACK_ERROR || used.exhausted || used.apply_a!=1)
            return 1;
        p.apply_a=apply_diagonal;
        lis_krylov_budget_default(&budget);
        /* A larger v1 prefix must be accepted without reading its tail. */
        budget.struct_size+=8;
        rc=lis_krylov_symmetric_budgeted(&p,&options,&ext,NULL,&budget,
                    values,vectors,residuals,&stats,&extra,&used);
        if(rc!=LIS_KRYLOV_SUCCESS || stats.nconv!=2 ||
           used.exhausted || used.total!=used.apply_a || !used.apply_a)
            return 1;
        rc=lis_krylov_select_converged(5,v,r,1,6,1e-8,idx,5,&count);
        if(rc || count!=2 || idx[0]!=0 || idx[1]!=4) return 1;
        count=-1;
        rc=lis_krylov_select_converged(5,v,r,1,6,1e-8,idx,1,&count);
        if(rc!=LIS_KRYLOV_BAD_ARGUMENT || count!=-1) return 1;
        puts("PASS callback budget and converged pair selection");
    }

    {
        lis_krylov_problem64 p64={0};
        double reference[2];
        lis_krylov_ext_options_default(&ext);
        lis_krylov_ext_stats_default(&extra);
        p64.n=p.n;p64.apply_a=apply_diagonal64;
        rc=lis_krylov_symmetric(&p,&options,reference,vectors,residuals,&stats);
        if(rc!=LIS_KRYLOV_SUCCESS)return 1;
        rc=lis_krylov_workspace_create64(&p64,&options,&ext,&w);
        if(rc || !w)return 1;
        rc=lis_krylov_symmetric_ex64(&p64,&options,&ext,w,values,vectors,
                                     residuals,&stats,&extra);
        if(rc || stats.nconv!=2 ||
           fabs(values[0]-reference[0])>1e-10 ||
           fabs(values[1]-reference[1])>1e-10)return 1;
        /* Same allocation works with the original callback ABI. */
        rc=lis_krylov_symmetric_ex(&p,&options,&ext,w,values,vectors,
                                   residuals,&stats,&extra);
        lis_krylov_workspace_destroy(w);w=NULL;
        if(rc || stats.nconv!=2)return 1;
        {
            lis_krylov_budget budget;
            lis_krylov_budget_stats used;
            int64_t count=0;
            lis_krylov_budget_default(&budget);
            budget.max_apply_a=1;
            rc=lis_krylov_symmetric_budgeted64(&p64,&options,&ext,NULL,&budget,
                        values,vectors,residuals,&stats,&extra,&used);
            if(rc!=LIS_KRYLOV_BUDGET_EXHAUSTED || used.apply_a!=1 ||
               stats.nreturned)return 1;
            rc=lis_krylov_count_interval64(INT64_C(5000000000),NULL,inertia64,
                                           -1,1,&count);
            if(rc || count!=INT64_C(5000000000))return 1;
        }
        p64.n=INT64_MAX;
        rc=lis_krylov_workspace_create64(&p64,&options,&ext,&w);
        if(rc!=LIS_KRYLOV_BAD_ARGUMENT || w)return 1;
        rc=lis_krylov_symmetric_ex64(&p64,&options,&ext,NULL,values,vectors,
                                     residuals,&stats,&extra);
        if(rc!=LIS_KRYLOV_BAD_ARGUMENT || stats.nreturned)return 1;
        {
            lis_krylov_options bad=options;
            const int64_t dims[]={INT64_MIN,-1,8,8};
            const int requested[]={1,1,INT_MIN,0};
            int i;
            bad.ncv=0;
            for(i=0;i<4;i++){
                p64.n=dims[i];bad.nev=requested[i];
                rc=lis_krylov_workspace_create64(&p64,&bad,&ext,&w);
                if(rc!=LIS_KRYLOV_BAD_ARGUMENT || w)return 1;
                rc=lis_krylov_symmetric_ex64(&p64,&bad,&ext,NULL,values,
                                             vectors,residuals,&stats,&extra);
                if(rc!=LIS_KRYLOV_BAD_ARGUMENT || stats.nreturned)return 1;
            }
        }
        {
            lis_krylov_options shift_options=options;
            lis_krylov_budget budget;
            lis_krylov_budget_stats used;
            int64_t count=0;
            double sigma=2.5;
            p64.n=8;p64.transform=LIS_KRYLOV_SHIFT_INVERT_B;
            p64.sigma=sigma;p64.ctx=&sigma;p64.apply_b=identity64;
            p64.solve=diagonal_inverse64;
            shift_options.ncv=8;shift_options.inner_tol=1e-10;
            rc=lis_krylov_solve_interval64(&p64,&shift_options,
                    diagonal_inertia64,1.5,3.5,values,vectors,residuals,
                    &stats,&count);
            if(rc || count!=2 || stats.nconv!=2 ||
               fabs(values[0]+values[1]-5)>1e-7 ||
               fabs(fabs(values[0]-values[1])-1)>1e-7)return 1;
            lis_krylov_budget_default(&budget);
            budget.max_apply_a=1;
            rc=lis_krylov_symmetric_budgeted64(&p64,&shift_options,&ext,NULL,
                    &budget,values,vectors,residuals,&stats,&extra,&used);
            if(rc!=LIS_KRYLOV_BUDGET_EXHAUSTED || !used.exhausted ||
               used.apply_a!=1 || !used.solve || stats.nreturned)return 1;
        }
        puts("PASS 64-bit vector API parity, shift-invert and overflow rejection");
    }

    /* On ILP32, n * ncv = 2^34 and n * (ncv/2) = 2^33.
     * Both products wrap to zero in size_t. Previously the solver
     * computed these products BEFORE its overflow check. The corrected
     * allocator must reject this shape without trying to allocate. */
    if (sizeof(size_t) == 4) {
        p.n = 131072;
        options.nev = 1;
        options.ncv = 131072;
        lis_krylov_ext_options_default(&ext);
        w = (lis_krylov_workspace *)(uintptr_t)1;
        rc = lis_krylov_workspace_create(&p, &options, &ext, &w);
        if (rc != LIS_KRYLOV_BAD_ARGUMENT || w != NULL) {
            fprintf(stderr, "FAIL ILP32 overflow detection: %d\n", rc);
            return 1;
        }
        puts("PASS ILP32 workspace multiplication overflow rejected");
    } else {
        puts("SKIP ILP32 overflow case (64-bit target)");
    }
    puts("PASS versioned API validation and NULL workspace cleanup");
    return 0;
}
