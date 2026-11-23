#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <limits.h>
#include <stdio.h>

/* 0: valid, 1: signed NNZ overflow, 2: fill/query mismatch,
 * 3: negative NNZ returned by the query callback. */
typedef struct { int mode; } USER_TO_CSR_CTX;

static LIS_INT user_get_row(void *data,LIS_INT row,LIS_INT *nnz,
                             LIS_INT *cols,LIS_SCALAR *vals)
{
    USER_TO_CSR_CTX *ctx=(USER_TO_CSR_CTX *)data;
    LIS_INT max_nnz=(sizeof(LIS_INT)==sizeof(int))
                   ? (LIS_INT)INT_MAX : (LIS_INT)LLONG_MAX;
    if(row<0 || row>=2 || nnz==NULL) return LIS_ERR_ILL_ARG;
    if(cols==NULL && vals==NULL)
    {
        *nnz=(ctx->mode==1 && row==0) ? max_nnz :
             (ctx->mode==3 && row==0) ? -1 : 1;
        return LIS_SUCCESS;
    }
    if(cols==NULL || vals==NULL || *nnz<1) return LIS_ERR_ILL_ARG;
    cols[0]=row;
    vals[0]=(LIS_SCALAR)(2+row);
    /* The callback writes no more than the supplied capacity but
     * reports a changed row size, which conversion must reject. */
    *nnz=(ctx->mode==2 && row==0) ? 2 : 1;
    return LIS_SUCCESS;
}

static int run_case(int mode)
{
    USER_TO_CSR_CTX ctx;
    LIS_MATRIX user=NULL,csr=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_INT err=LIS_SUCCESS;
    int failed=1;
    ctx.mode=mode;
    err=lis_matrix_create(LIS_COMM_WORLD,&user);
    if(err) goto done;
    err=lis_matrix_set_size(user,0,2);
    if(err) goto done;
    err=lis_matrix_set_user(user,&ctx,NULL,NULL);
    if(err) goto done;
    err=lis_matrix_set_user_row(user,user_get_row);
    if(err) goto done;
    err=lis_matrix_duplicate(user,&csr);
    if(err) goto done;
    err=lis_matrix_set_type(csr,LIS_MATRIX_CSR);
    if(err) goto done;
    err=lis_matrix_convert(user,csr);
    if(mode!=0)
    {
        if(err==LIS_ERR_ILL_ARG) failed=0;
        goto done;
    }
    if(err) goto done;
    err=lis_vector_duplicate(csr,&x);
    if(err) goto done;
    err=lis_vector_duplicate(csr,&y);
    if(err) goto done;
    lis_vector_set_all((LIS_SCALAR)1,x);
    err=lis_matvec(csr,x,y);
    if(err) goto done;
    if(y->value[0]!=(LIS_SCALAR)2 ||
       y->value[1]!=(LIS_SCALAR)3) goto done;
    failed=0;
done:
    if(failed) fprintf(stderr,"USER_TO_CSR FAIL mode=%d rc=%d\n",
                       mode,(int)err);
    if(y) lis_vector_destroy(y);
    if(x) lis_vector_destroy(x);
    if(csr) lis_matrix_destroy(csr);
    if(user) lis_matrix_destroy(user);
    return failed;
}

int main(int argc,char **argv)
{
#ifdef USE_MPI
    /* This regression exercises the serial USER-to-CSR path.  MPI USER
     * duplication without communication metadata is covered separately
     * by testusermpiduplicate (#247).  The pre-existing MPI conversion
     * path would otherwise segfault before any bound check is tested. */
    (void)argc;
    (void)argv;
    puts("USER_TO_CSR_BOUNDS SKIP (serial-only conversion regression)");
    return 77; /* Automake's standard SKIP result. */
#else
    int failed=0;
    if(lis_initialize(&argc,&argv)) return 1;
    failed |= run_case(1);
    failed |= run_case(2);
    failed |= run_case(3);
    failed |= run_case(0);
    lis_finalize();
    printf("USER_TO_CSR_BOUNDS %s\n",failed ? "FAILED" : "PASSED");
    return failed ? 1 : 0;
#endif
}
