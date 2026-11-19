/* Optional batched USER row integration and cache lifecycle regression. */
#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include "lis_matrix_materialize.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

typedef struct {
    LIS_INT n,scalar_queries,batch_queries,fail_mismatch;
    LIS_INT fail_fill_once,negative_col_once,reverse_order;
    LIS_INT callback_oom_once,callback_repeated_error;
} CONTEXT;

static LIS_INT count_row(CONTEXT *ctx,LIS_INT row)
{
    return (row>0)+(LIS_INT)1+(row+1<ctx->n);
}

static LIS_INT get_row(void *opaque,LIS_INT row,LIS_INT *nnz,
                       LIS_INT *cols,LIS_SCALAR *vals)
{
    CONTEXT *ctx=(CONTEXT *)opaque;
    LIS_INT k=0,need;
    if(row<0 || row>=ctx->n || nnz==NULL) return LIS_ERR_ILL_ARG;
    need=count_row(ctx,row);
#ifdef _OPENMP
#pragma omp atomic update
#endif
    ctx->scalar_queries++;
    if(cols==NULL && vals==NULL){*nnz=need;return LIS_SUCCESS;}
    if(cols==NULL || vals==NULL || *nnz<need)return LIS_ERR_ILL_ARG;
    if(row>0){cols[k]=row-1;vals[k++]=-1.0;}
    cols[k]=row;vals[k++]=4.0;
    if(row+1<ctx->n){cols[k]=row+1;vals[k++]=-1.0;}
    *nnz=need;
    return LIS_SUCCESS;
}

static LIS_INT get_batch(void *opaque,LIS_INT first,LIS_INT last,
                         LIS_INT *offsets,LIS_INT capacity,
                         LIS_INT *cols,LIS_SCALAR *vals)
{
    CONTEXT *ctx=(CONTEXT *)opaque;
    LIS_INT row,k=0,i;
    if(!offsets || first<0 || last>ctx->n || first>=last)
        return LIS_ERR_ILL_ARG;
#ifdef _OPENMP
#pragma omp atomic update
#endif
    ctx->batch_queries++;
    offsets[0]=0;
    for(row=first;row<last;row++)
        offsets[row-first+1]=offsets[row-first]+count_row(ctx,row);
    if(!cols && !vals)return LIS_SUCCESS;
    if(!cols || !vals || capacity<offsets[last-first])
        return LIS_ERR_ILL_ARG;
    for(row=first;row<last;row++)
    {
        if(row>0){cols[k]=row-1;vals[k++]=-1.0;}
        cols[k]=row;vals[k++]=4.0;
        if(row+1<ctx->n){cols[k]=row+1;vals[k++]=-1.0;}
    }
    if(ctx->fail_fill_once)
    {
        ctx->fail_fill_once=0;
        return LIS_ERR_ILL_ARG;
    }
    if(ctx->callback_oom_once)
    {
        ctx->callback_oom_once=0;
        return LIS_OUT_OF_MEMORY;
    }
    if(ctx->callback_repeated_error>0)
    {
        ctx->callback_repeated_error--;
        return LIS_ERR_ILL_ARG;
    }
    if(ctx->negative_col_once)
    {
        ctx->negative_col_once=0;
        cols[0]=-1;
    }
    if(ctx->reverse_order)
    {
        for(row=first;row<last;row++)
        {
            LIS_INT lo=offsets[row-first],hi=offsets[row-first+1]-1;
            while(lo<hi)
            {
                LIS_INT tcol=cols[lo];
                LIS_SCALAR tval=vals[lo];
                cols[lo]=cols[hi]; vals[lo]=vals[hi];
                cols[hi]=tcol; vals[hi]=tval;
                lo++; hi--;
            }
        }
    }
    if(ctx->fail_mismatch)
    {
        i=last-first;
        offsets[i]++; /* Deliberate contract violation on fill. */
    }
    return LIS_SUCCESS;
}

static LIS_INT reset_batch_source(LIS_MATRIX A,LIS_PRECON_ROW_SOURCE *src)
{
    lis_precon_row_source_destroy(src);
    return lis_precon_row_source_init(A,src);
}

static int run_mode(int use_batch)
{
    enum { N=256 };
    CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_PRECON_ROW_SOURCE src;
    LIS_INT err=0,i,j,k,nnz;
    int failed=1,initialized=0;
    memset(&ctx,0,sizeof(ctx));
    memset(&src,0,sizeof(src));
    ctx.n=N;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,N);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,get_row);if(err)goto done;
    if(use_batch)
    {
        err=lis_matrix_set_user_rows_batch(A,get_batch);if(err)goto done;
    }
    err=lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);if(err)goto done;
    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&y);if(err)goto done;
    for(i=0;i<N;i++) x->value[i]=1.0+0.01*i;

    err=lis_matvec(A,x,y);if(err)goto done;
    for(i=0;i<N;i++)
    {
        double expected=4.0*x->value[i];
        if(i>0)expected-=x->value[i-1];
        if(i+1<N)expected-=x->value[i+1];
        if(fabs((double)y->value[i]-expected)>1e-11)goto done;
    }
    err=lis_matvech(A,x,y);if(err)goto done;
    for(i=0;i<N;i++)
    {
        double expected=4.0*x->value[i];
        if(i>0)expected-=x->value[i-1];
        if(i+1<N)expected-=x->value[i+1];
        if(fabs((double)y->value[i]-expected)>1e-11)goto done;
    }
    if(use_batch && ctx.scalar_queries!=0)
    {
        fprintf(stderr,"USER matvec unexpectedly used the scalar callback\n");
        goto done;
    }
    err=lis_precon_row_source_init(A,&src);if(err)goto done;
    initialized=1;
    /* Non-monotone order must trigger safe cache refresh. */
    for(k=0;k<3;k++)
    {
        LIS_INT row=k==0?97:(k==1?3:99);
        err=lis_precon_row_source_get(&src,row);if(err)goto done;
        nnz=count_row(&ctx,row);
        if(src.nnz!=nnz)goto done;
        for(j=0;j<src.nnz;j++)
        {
            if(src.index[j]<0 || src.index[j]>=N)goto done;
        }
    }


    if(use_batch && ctx.batch_queries>=2*N)
    {
        fprintf(stderr,"batch callback count too high %d\n",
                (int)ctx.batch_queries);
        goto done;
    }
    if(!use_batch && ctx.scalar_queries<2*N)goto done;

    if(use_batch)
    {
        /* A second distant miss favors one-row callbacks. Eight adjacent
         * requests restore batching after the random-access phase. */
        if(!src.batch_prefer_scalar)goto done;
        for(k=100;k<=107;k++)
        {
            err=lis_precon_row_source_get(&src,k);
            if(err || src.nnz!=3)goto done;
        }
        if(src.batch_prefer_scalar)goto done;
        err=reset_batch_source(A,&src);if(err)goto done;
        /* Explicit mismatch must be rejected, not produce stale row data. */
        ctx.fail_mismatch=1;
        err=lis_precon_row_source_get(&src,210);
        if(err!=LIS_ERR_ILL_ARG)goto done;
        ctx.fail_mismatch=0;
        err=lis_precon_row_source_get(&src,210);
        if(err || src.nnz!=3)goto done;
        err=reset_batch_source(A,&src);if(err)goto done;
        ctx.fail_fill_once=1;
        err=lis_precon_row_source_get(&src,230);
        if(err!=LIS_ERR_ILL_ARG)goto done;
        err=lis_precon_row_source_get(&src,230);
        if(err || src.nnz!=3)goto done;
        err=reset_batch_source(A,&src);if(err)goto done;
        ctx.negative_col_once=1;
        err=lis_precon_row_source_get(&src,247);
        if(err!=LIS_ERR_ILL_ARG)goto done;
        err=lis_precon_row_source_get(&src,247);
        if(err || src.nnz!=3)goto done;
        /* Column ordering has never been required for USER rows. */
        err=reset_batch_source(A,&src);if(err)goto done;
        ctx.reverse_order=1;
        err=lis_precon_row_source_get(&src,1);
        if(err || src.nnz!=3 || src.index[0]!=2)goto done;
        /* Application callback OOM is an error, NOT an internal
         * cache-allocation OOM: never silently downgrade it. */
        err=reset_batch_source(A,&src);if(err)goto done;
        ctx.scalar_queries=0; /* inspect only callback-error segment */
        ctx.callback_oom_once=1;
        err=lis_precon_row_source_get(&src,120);
        if(err!=LIS_OUT_OF_MEMORY || src.batch_disabled)goto done;
        err=lis_precon_row_source_get(&src,120);
        if(err || src.nnz!=3 || ctx.scalar_queries!=0)goto done;

        /* Four failures of the same callback/row must all propagate.
         * Locality adaptation is not allowed to hide a persistent error
         * by switching to the scalar callback on a later retry. */
        err=reset_batch_source(A,&src);if(err)goto done;
        ctx.callback_repeated_error=4;
        for(k=0;k<4;k++)
        {
            err=lis_precon_row_source_get(&src,160);
            if(err!=LIS_ERR_ILL_ARG || src.batch_prefer_scalar ||
               ctx.scalar_queries!=0)goto done;
        }
        err=lis_precon_row_source_get(&src,160);
        if(err || src.nnz!=3 || ctx.scalar_queries!=0)goto done;
    }
    failed=0;
done:
    if(failed)
        fprintf(stderr,"USER_BATCH FAIL mode=%d err=%d scalar=%d batch=%d\n",
                use_batch,(int)err,(int)ctx.scalar_queries,(int)ctx.batch_queries);
    else
        printf("USER_BATCH PASS mode=%d scalar_calls=%d batch_calls=%d\n",
               use_batch,(int)ctx.scalar_queries,(int)ctx.batch_queries);
    if(initialized)lis_precon_row_source_destroy(&src);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    return failed;
}

/* Verify registration reset and destruction without expanding LIS_MATRIX. */
static int check_registration_lifecycle(void)
{
    CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_INT err=0;
    int failed=1;
    memset(&ctx,0,sizeof(ctx));
    ctx.n=32;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,get_row);if(err)goto done;
    if(lis_matrix_user_batch_lookup(A)!=NULL)goto done;
    err=lis_matrix_set_user_rows_batch(A,get_batch);if(err)goto done;
    if(lis_matrix_user_batch_lookup(A)!=get_batch)goto done;
    err=lis_matrix_set_user_rows_batch(A,NULL);if(err)goto done;
    if(lis_matrix_user_batch_lookup(A)!=NULL)goto done;
    err=lis_matrix_set_user_rows_batch(A,get_batch);if(err)goto done;
    err=lis_matrix_set_user_row(A,get_row);if(err)goto done;
    if(lis_matrix_user_batch_lookup(A)!=NULL)goto done;
    /* Reinitializing an already configured USER descriptor is not supported
     * by LIS_MATRIX_CHECK_SET. Destroy must still forget its callback. */
    err=lis_matrix_set_user_rows_batch(A,get_batch);if(err)goto done;
    failed=0;
done:
    if(A)
    {
        LIS_MATRIX old=A;
        lis_matrix_destroy(A);
        if(lis_matrix_user_batch_lookup(old)!=NULL) failed=1;
    }
    if(failed)fprintf(stderr,"USER_BATCH lifecycle FAIL err=%d\n",(int)err);
    return failed;
}

/* Exercise the private binding registry concurrently on distinct matrices.
 * Setters still must not race with a solve for the SAME matrix. */
#ifdef _OPENMP
static int check_parallel_registrations(void)
{
    enum { MATRICES=8 };
    LIS_MATRIX matrices[MATRICES];
    CONTEXT contexts[MATRICES];
    LIS_INT i,err=0;
    int faults=0;
    memset(matrices,0,sizeof(matrices));
    memset(contexts,0,sizeof(contexts));
    for(i=0;i<MATRICES;i++)
    {
        contexts[i].n=32;
        err=lis_matrix_create(LIS_COMM_WORLD,&matrices[i]);
        if(err)goto done;
        err=lis_matrix_set_size(matrices[i],0,32);
        if(err)goto done;
        err=lis_matrix_set_user(matrices[i],&contexts[i],NULL,NULL);
        if(err)goto done;
        err=lis_matrix_set_user_row(matrices[i],get_row);
        if(err)goto done;
    }
#pragma omp parallel for num_threads(8) reduction(+:faults)
    for(i=0;i<MATRICES;i++)
    {
        if(lis_matrix_set_user_rows_batch(matrices[i],get_batch)
           !=LIS_SUCCESS) faults++;
        if(lis_matrix_user_batch_lookup(matrices[i])!=get_batch)
            faults++;
        if(lis_matrix_set_user_rows_batch(matrices[i],NULL)
           !=LIS_SUCCESS) faults++;
        if(lis_matrix_user_batch_lookup(matrices[i])!=NULL)
            faults++;
        if(lis_matrix_set_user_rows_batch(matrices[i],get_batch)
           !=LIS_SUCCESS) faults++;
    }
done:
    for(i=0;i<MATRICES;i++)
        if(matrices[i])lis_matrix_destroy(matrices[i]);
    if(err || faults)
    {
        fprintf(stderr,"USER_BATCH parallel registry FAIL err=%d faults=%d\n",
                (int)err,faults);
        return 1;
    }
    return 0;
}
#endif

#ifdef LIS_TEST_OOM
static int check_batch_cache_oom(void)
{
    CONTEXT ctx;
    LIS_MATRIX A=NULL;
    LIS_PRECON_ROW_SOURCE src;
    LIS_INT err=0;
    int bad=1,initialized=0;
    memset(&ctx,0,sizeof(ctx));memset(&src,0,sizeof(src));
    ctx.n=64;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,ctx.n);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,get_row);if(err)goto done;
    lis_test_oom_begin(1);
    err=lis_matrix_set_user_rows_batch(A,get_batch);
    lis_test_oom_end();
    if(err!=LIS_OUT_OF_MEMORY || lis_matrix_user_batch_lookup(A)!=NULL)
        goto done;
    err=lis_matrix_set_user_rows_batch(A,get_batch);if(err)goto done;
    /* Registration is transactional: a failed replacement keeps old data. */
    lis_test_oom_begin(1);
    err=lis_matrix_set_user_rows_batch(A,get_batch);
    lis_test_oom_end();
    if(err!=LIS_OUT_OF_MEMORY ||
       lis_matrix_user_batch_lookup(A)!=get_batch)goto done;
    err=lis_precon_row_source_init(A,&src);if(err)goto done;
    initialized=1;
    lis_test_oom_begin(1);
    err=lis_precon_row_source_get(&src,3);
    lis_test_oom_end();
    if(err || !src.batch_disabled || src.nnz!=3 ||
       src.index[1]!=3 || ctx.scalar_queries==0)
        goto done;
    err=lis_precon_row_source_get(&src,4);
    if(err || src.nnz!=3)goto done;
    bad=0;
done:
    lis_test_oom_end();
    if(initialized)lis_precon_row_source_destroy(&src);
    if(A)lis_matrix_destroy(A);
    if(bad)fprintf(stderr,"USER_BATCH OOM fallback FAIL err=%d\n",(int)err);
    return bad;
}
#endif

/* Materialization is a core consumer of USER row sources. Verify its
 * entire CSR representation, not just matvec output, in both modes. */
static int check_materialization(int use_batch)
{
    enum { N=256 };
    CONTEXT ctx;
    LIS_MATRIX A=NULL,B=NULL;
    LIS_INT err=0,i,k,expected,cols[3];
    LIS_SCALAR vals[3];
    int failed=1;
    memset(&ctx,0,sizeof(ctx));
    ctx.n=N;
    err=lis_matrix_create(LIS_COMM_WORLD,&A);if(err)goto done;
    err=lis_matrix_set_size(A,0,N);if(err)goto done;
    err=lis_matrix_set_user(A,&ctx,NULL,NULL);if(err)goto done;
    err=lis_matrix_set_user_row(A,get_row);if(err)goto done;
    if(use_batch)
    {
        err=lis_matrix_set_user_rows_batch(A,get_batch);if(err)goto done;
    }
    err=lis_matrix_set_user_row_thread_safe(A,LIS_TRUE);if(err)goto done;
    err=lis_matrix_materialize_csr_owned(A,&B);if(err)goto done;
    if(B==NULL || B->matrix_type!=LIS_MATRIX_CSR ||
       B->ptr==NULL || B->index==NULL || B->value==NULL ||
       B->ptr[N]!=3*N-2)goto done;
    for(i=0;i<N;i++)
    {
        expected=count_row(&ctx,i);
        if(B->ptr[i+1]-B->ptr[i]!=expected)goto done;
        k=0;
        if(i>0){cols[k]=i-1;vals[k++]=-1.0;}
        cols[k]=i; vals[k++]=4.0;
        if(i+1<N){cols[k]=i+1;vals[k++]=-1.0;}
        for(k=0;k<expected;k++)
            if(B->index[B->ptr[i]+k]!=cols[k] ||
               fabs((double)B->value[B->ptr[i]+k]-
                    (double)vals[k])>1e-12)goto done;
    }
    if(use_batch && (ctx.scalar_queries!=0 || ctx.batch_queries==0))
        goto done;
    if(!use_batch && ctx.scalar_queries==0)goto done;
    failed=0;
done:
    if(failed)fprintf(stderr,
        "USER_BATCH materialization FAIL mode=%d err=%d\\n",
        use_batch,(int)err);
    if(B)lis_matrix_destroy(B);
    if(A)lis_matrix_destroy(A);
    return failed;
}

int main(int argc,char **argv)
{
    int failed;
    if(lis_initialize(&argc,&argv)!=LIS_SUCCESS)return 2;
    failed=run_mode(0)|run_mode(1);
    failed|=check_materialization(0);
    failed|=check_materialization(1);
    failed|=check_registration_lifecycle();
#ifdef _OPENMP
    failed|=check_parallel_registrations();
#endif
#ifdef LIS_TEST_OOM
    failed|=check_batch_cache_oom();
#endif
    lis_finalize();
    return failed;
}
