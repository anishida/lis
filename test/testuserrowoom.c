#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lislib.h"
#include <stdio.h>
#include <math.h>

#if defined(LIS_TEST_OOM) && !defined(USE_MPI) && !defined(_COMPLEX)
static LIS_INT row(void *ctx,LIS_INT i,LIS_INT *nnz,LIS_INT *cols,LIS_SCALAR *vals)
{
    LIS_INT k,need=i==0?1:3;
    (void)ctx;
    if(!cols && !vals){*nnz=need;return LIS_SUCCESS;}
    if(!cols || !vals || *nnz<need)return LIS_ERR_ILL_ARG;
    for(k=0;k<need;k++){cols[k]=k;vals[k]=1.0+i+k;}
    *nnz=need;
    return LIS_SUCCESS;
}
static LIS_INT batch(void *ctx,LIS_INT begin,LIS_INT end,LIS_INT *offsets,
                    LIS_INT capacity,LIS_INT *cols,LIS_SCALAR *vals)
{
    LIS_INT i,k=0,n;
    offsets[0]=0;
    for(i=begin;i<end;i++)
    {
        n=i==0?1:3;
        offsets[i-begin+1]=offsets[i-begin]+n;
    }
    if(!cols && !vals)return LIS_SUCCESS;
    if(!cols || !vals || capacity<offsets[end-begin])return LIS_ERR_ILL_ARG;
    for(i=begin;i<end;i++){
        n=capacity-k;
        if(row(ctx,i,&n,cols+k,vals+k))return LIS_ERR_ILL_ARG;
        k+=n;
    }
    return LIS_SUCCESS;
}
static int valid(LIS_PRECON_ROW_SOURCE *s,LIS_INT i)
{
    LIS_INT k,need=i==0?1:3;
    if(s->nnz!=need)return 0;
    for(k=0;k<need;k++)
        if(s->index[k]!=k || s->value[k]!=(LIS_SCALAR)(1.0+i+k))return 0;
    return 1;
}
static int exercise(LIS_MATRIX A,int batched,size_t fail_at)
{
    LIS_PRECON_ROW_SOURCE s;
    size_t live=lis_test_oom_live(),before,calls;
    LIS_INT *old_index=NULL;
    LIS_SCALAR *old_value=NULL;
    LIS_INT err;
    int failed=1,initialized=0;
    if(lis_precon_row_source_init(A,&s))goto done;
    initialized=1;
    if(!batched && (lis_precon_row_source_get(&s,0) || !valid(&s,0)))goto done;
    old_index=s.index_buffer;
    old_value=s.value_buffer;
    before=lis_test_oom_live();
    lis_test_oom_begin(fail_at);
    err=lis_precon_row_source_get(&s,1);
    calls=lis_test_oom_calls();
    lis_test_oom_end();
    if(batched)
    {
        if(err || !s.batch_disabled || calls<fail_at || !valid(&s,1))goto done;
    }
    else
    {
        if(err!=LIS_OUT_OF_MEMORY || calls!=fail_at ||
           lis_test_oom_live()!=before || s.capacity!=1 ||
           s.index_buffer!=old_index || s.value_buffer!=old_value ||
           !old_index || !old_value || old_index[0]!=0 || old_value[0]!=1.0 ||
           s.index!=NULL || s.value!=NULL || s.nnz!=0)goto done;
        /* A failed get invalidates the published view, not owned scratch. */
    }
    /* A fresh request after failure must still return the correct row. */
    if(lis_precon_row_source_get(&s,1) || !valid(&s,1) ||
       lis_precon_row_source_get(&s,0) || !valid(&s,0))goto done;
    failed=0;
done:
    lis_test_oom_end();
    if(initialized)lis_precon_row_source_destroy(&s);
    if(lis_test_oom_live()!=live)failed=1;
    if(failed)fprintf(stderr,"USER_ROW_OOM failed batch=%d allocation=%zu\n",batched,fail_at);
    return failed;
}

/* Seed a short batch at the end of the matrix, then request a larger
 * batch from row zero. Failed growth must retain the old owned buffers
 * and fall back to scalar access without exposing a stale batch view. */
static int exercise_batch_growth(LIS_MATRIX A,size_t fail_at)
{
    LIS_PRECON_ROW_SOURCE s;
    size_t live=lis_test_oom_live(),before,calls;
    LIS_INT *old_index=NULL;
    LIS_SCALAR *old_value=NULL;
    LIS_INT old_capacity=0,k,err;
    int failed=1,initialized=0;
    if(lis_precon_row_source_init(A,&s))goto done;
    initialized=1;
    if(lis_precon_row_source_get(&s,3) || !valid(&s,3) ||
       s.batch_disabled || s.batch_capacity!=3 ||
       s.batch_first!=3 || s.batch_last!=4)goto done;
    old_index=s.batch_index; old_value=s.batch_value;
    old_capacity=s.batch_capacity;
    if(!old_index || !old_value)goto done;
    before=lis_test_oom_live();
    lis_test_oom_begin(fail_at);
    err=lis_precon_row_source_get(&s,0);
    calls=lis_test_oom_calls();
    lis_test_oom_end();
    if(err || calls<fail_at || !s.batch_disabled ||
       s.batch_first!=0 || s.batch_last!=0 ||
       s.batch_index!=old_index || s.batch_value!=old_value ||
       s.batch_capacity!=old_capacity ||
       !valid(&s,0) || s.index!=s.index_buffer || s.value!=s.value_buffer ||
       lis_test_oom_live()!=before+2)goto done;
    for(k=0;k<old_capacity;k++)
        if(old_index[k]!=k || old_value[k]!=(LIS_SCALAR)(4.0+k))goto done;
    if(lis_precon_row_source_get(&s,1) || !valid(&s,1) ||
       lis_precon_row_source_get(&s,3) || !valid(&s,3) ||
       !s.batch_disabled || s.batch_index!=old_index ||
       s.batch_value!=old_value || lis_test_oom_live()!=before+2)goto done;
    failed=0;
done:
    lis_test_oom_end();
    if(initialized)lis_precon_row_source_destroy(&s);
    if(lis_test_oom_live()!=live)failed=1;
    if(failed)fprintf(stderr,"USER_ROW_OOM batch growth failed allocation=%zu\n",fail_at);
    return failed;
}
#endif

int main(int argc,char **argv)
{
#if defined(LIS_TEST_OOM) && !defined(USE_MPI) && !defined(_COMPLEX)
    LIS_MATRIX A=NULL;
    int failed=1;
    if(lis_initialize(&argc,&argv))return 1;
    if(lis_matrix_create(LIS_COMM_WORLD,&A) || lis_matrix_set_size(A,0,4) ||
       lis_matrix_set_user(A,NULL,NULL,NULL) || lis_matrix_set_user_row(A,row))goto done;
    failed=exercise(A,0,1)|exercise(A,0,2);
    if(lis_matrix_set_user_rows_batch(A,batch)){failed=1;goto done;}
    failed|=exercise(A,1,1)|exercise(A,1,2);
    failed|=exercise_batch_growth(A,1)|exercise_batch_growth(A,2);
done:
    if(A)lis_matrix_destroy(A);
    lis_finalize();
    printf("USER_ROW_OOM %s\n",failed?"FAIL":"PASS");
    return failed?1:0;
#else
    (void)argc;(void)argv;
    puts("USER_ROW_OOM SKIP: real serial OOM build required");
    return 77;
#endif
}
