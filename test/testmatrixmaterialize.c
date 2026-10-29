#include "lis.h"
#include "lis_matrix_materialize.h"
#include <stdio.h>

typedef struct
{
    LIS_INT variant;
    LIS_REAL scale;
} TEST_ROW_CTX;

static void
scalar_set_real(LIS_SCALAR *value,LIS_REAL x)
{
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
    (*value)[0]=x;
    (*value)[1]=(LIS_REAL)0.0;
#else
    *value=(LIS_SCALAR)x;
#endif
}

static int
scalar_equal_real(const LIS_SCALAR *value,LIS_REAL expected)
{
    LIS_SCALAR target;
    scalar_set_real(&target,expected);
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
    return (*value)[0]==target[0] && (*value)[1]==target[1];
#else
    return *value==target;
#endif
}

static LIS_INT
test_get_row(
    void *user_data,
    LIS_INT local_row,
    LIS_INT *nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    TEST_ROW_CTX *ctx=(TEST_ROW_CTX *)user_data;
    LIS_INT required;
    LIS_REAL s;

    if(ctx==NULL || nnz==NULL || local_row<0 || local_row>=3)
        return LIS_ERR_ILL_ARG;

    required=(local_row==1) ? 3 : 2;
    if(ctx->variant!=0 && local_row!=1)
        required=3;

    if(cols==NULL && vals==NULL)
    {
        *nnz=required;
        return LIS_SUCCESS;
    }
    if(cols==NULL || vals==NULL || *nnz<required)
    {
        *nnz=required;
        return LIS_ERR_ILL_ARG;
    }

    s=ctx->scale;
    if(local_row==0)
    {
        cols[0]=0; scalar_set_real(&vals[0],(LIS_REAL)2.0*s);
        cols[1]=1; scalar_set_real(&vals[1],(LIS_REAL)-1.0*s);
        if(ctx->variant!=0)
        {
            cols[2]=2; scalar_set_real(&vals[2],(LIS_REAL)0.5*s);
        }
    }
    else if(local_row==1)
    {
        cols[0]=0; scalar_set_real(&vals[0],(LIS_REAL)-1.0*s);
        cols[1]=1; scalar_set_real(&vals[1],(LIS_REAL)2.0*s);
        cols[2]=2; scalar_set_real(&vals[2],(LIS_REAL)-1.0*s);
    }
    else
    {
        if(ctx->variant!=0)
        {
            cols[0]=0; scalar_set_real(&vals[0],(LIS_REAL)0.5*s);
            cols[1]=1; scalar_set_real(&vals[1],(LIS_REAL)-1.0*s);
            cols[2]=2; scalar_set_real(&vals[2],(LIS_REAL)2.0*s);
        }
        else
        {
            cols[0]=1; scalar_set_real(&vals[0],(LIS_REAL)-1.0*s);
            cols[1]=2; scalar_set_real(&vals[1],(LIS_REAL)2.0*s);
        }
    }

    *nnz=required;
    return LIS_SUCCESS;
}

static LIS_INT
create_user_matrix(TEST_ROW_CTX *ctx,LIS_MATRIX *A)
{
    LIS_INT err;

    *A=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,A);
    if(err) return err;
    err=lis_matrix_set_size(*A,0,3);
    if(err) goto fail;
    err=lis_matrix_set_user(*A,ctx,NULL,NULL);
    if(err) goto fail;
    err=lis_matrix_set_user_row(*A,test_get_row);
    if(err) goto fail;
    err=lis_matrix_set_user_row_thread_safe(*A,LIS_TRUE);
    if(err) goto fail;
    return LIS_SUCCESS;

fail:
    lis_matrix_destroy(*A);
    *A=NULL;
    return err;
}

static int
check_csr_matches_context(LIS_MATRIX A,const TEST_ROW_CTX *ctx)
{
    LIS_INT cols[3];
    LIS_SCALAR vals[3];
    LIS_INT i,k,nnz,p0,p1,err;
    TEST_ROW_CTX copy=*ctx;

    if(A==NULL || A->matrix_type!=LIS_MATRIX_CSR)
        return 1;

    for(i=0;i<3;i++)
    {
        nnz=3;
        err=test_get_row(&copy,i,&nnz,cols,vals);
        if(err) return 1;
        p0=A->ptr[i];
        p1=A->ptr[i+1];
        if(p1-p0!=nnz) return 1;
        for(k=0;k<nnz;k++)
        {
            if(A->index[p0+k]!=cols[k]) return 1;
#if defined(_COMPLEX) && !defined(HAVE_COMPLEX_H)
            if(A->value[p0+k][0]!=vals[k][0] ||
               A->value[p0+k][1]!=vals[k][1]) return 1;
#else
            if(A->value[p0+k]!=vals[k]) return 1;
#endif
        }
    }
    return 0;
}

static int
check_operator_values(LIS_MATRIX A,LIS_REAL scale)
{
    /* Both operands use variant 1 here, so their sum is the same pattern. */
    TEST_ROW_CTX expected;
    expected.variant=1;
    expected.scale=scale;
    return check_csr_matches_context(A,&expected);
}

int main(int argc,char **argv)
{
    TEST_ROW_CTX actx={0,(LIS_REAL)1.0};
    TEST_ROW_CTX bctx={1,(LIS_REAL)1.0};
    LIS_MATRIX A=NULL,B=NULL,op=NULL,csr=NULL,previous=NULL;
    LIS_MATRIX_CSR_CACHE acache=NULL,opcache=NULL;
    LIS_INT action=-1,err;
    int failed=0;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=create_user_matrix(&actx,&A);
    if(err) { failed=1; goto cleanup; }

    err=lis_matrix_materialize_csr_cached(A,&acache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD ||
       check_csr_matches_context(csr,&actx))
    {
        fprintf(stderr,"initial cached CSR build failed\n");
        failed=1; goto cleanup;
    }
    previous=csr;

    err=lis_matrix_materialize_csr_cached(A,&acache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REUSE || csr!=previous)
    {
        fprintf(stderr,"unchanged cached CSR was not reused\n");
        failed=1; goto cleanup;
    }

    actx.scale=(LIS_REAL)2.0;
    err=lis_matrix_mark_values_changed(A);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_materialize_csr_cached(A,&acache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES || csr!=previous ||
       check_csr_matches_context(csr,&actx))
    {
        fprintf(stderr,"value-only cached CSR refresh failed\n");
        failed=1; goto cleanup;
    }

    actx.variant=1;
    err=lis_matrix_mark_structure_changed(A);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_materialize_csr_cached(A,&acache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD || csr==previous ||
       check_csr_matches_context(csr,&actx))
    {
        fprintf(stderr,"structural cached CSR rebuild failed\n");
        failed=1; goto cleanup;
    }

    err=create_user_matrix(&bctx,&B);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_create_operator(
        (LIS_SCALAR)1.0,A,(LIS_SCALAR)1.0,B,&op);
    if(err) { failed=1; goto cleanup; }

    err=lis_matrix_materialize_csr_cached(op,&opcache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD ||
       check_operator_values(csr,(LIS_REAL)3.0))
    {
        fprintf(stderr,"operator cached CSR build failed\n");
        failed=1; goto cleanup;
    }
    previous=csr;

    /* The operator object is unchanged; only one nested operand is updated. */
    bctx.scale=(LIS_REAL)4.0;
    err=lis_matrix_mark_values_changed(B);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_materialize_csr_cached(op,&opcache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES || csr!=previous ||
       check_operator_values(csr,(LIS_REAL)6.0))
    {
        fprintf(stderr,"nested operator value refresh failed\n");
        failed=1; goto cleanup;
    }

    /* A nested structural generation change conservatively rebuilds. */
    bctx.variant=0;
    err=lis_matrix_mark_structure_changed(B);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_materialize_csr_cached(op,&opcache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD || csr==previous)
    {
        fprintf(stderr,"nested operator structural rebuild failed\n");
        failed=1; goto cleanup;
    }

    printf("LIS_MATRIX_MATERIALIZE_CACHE PASSED\n");

cleanup:
    lis_matrix_csr_cache_destroy(opcache);
    lis_matrix_csr_cache_destroy(acache);
    if(op) lis_matrix_destroy(op);
    if(B) lis_matrix_destroy(B);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed ? 1 : 0;
}
