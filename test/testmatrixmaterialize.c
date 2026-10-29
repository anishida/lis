#include "lis.h"
#include "lis_matrix.h"
#include "lis_matrix_materialize.h"
#include "lis_matrix_operator.h"
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
materialize_test_matvec(
    void *user_data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    LIS_INT i,k,nnz,err,cols[3];
    LIS_SCALAR vals[3];
    TEST_ROW_CTX *ctx=(TEST_ROW_CTX *)user_data;

    if(ctx==NULL || x==NULL || y==NULL)
        return LIS_ERR_ILL_ARG;

    for(i=0;i<3;i++)
    {
        nnz=3;
        err=test_get_row(ctx,i,&nnz,cols,vals);
        if(err) return err;
        y[i]=(LIS_SCALAR)0.0;
        for(k=0;k<nnz;k++)
            y[i]+=vals[k]*x[cols[k]];
    }
    return LIS_SUCCESS;
}

static LIS_INT
materialize_test_matvech(
    void *user_data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    /* The test matrix has real-valued entries in real and complex builds. */
    return materialize_test_matvec(user_data,x,y);
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
    err=lis_matrix_set_user(*A,ctx,materialize_test_matvec,materialize_test_matvech);
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

static int
test_context_entry(
    const TEST_ROW_CTX *ctx,
    LIS_INT row,
    LIS_INT col,
    LIS_REAL *value)
{
    LIS_REAL s;

    if(ctx==NULL || value==NULL || row<0 || row>=3 || col<0 || col>=3)
        return 0;

    s=ctx->scale;
    if(row==0)
    {
        if(col==0) { *value=(LIS_REAL)2.0*s; return 1; }
        if(col==1) { *value=(LIS_REAL)-1.0*s; return 1; }
        if(col==2 && ctx->variant!=0)
        {
            *value=(LIS_REAL)0.5*s;
            return 1;
        }
        return 0;
    }

    if(row==1)
    {
        if(col==0) { *value=(LIS_REAL)-1.0*s; return 1; }
        if(col==1) { *value=(LIS_REAL)2.0*s; return 1; }
        if(col==2) { *value=(LIS_REAL)-1.0*s; return 1; }
        return 0;
    }

    if(col==0 && ctx->variant!=0)
    {
        *value=(LIS_REAL)0.5*s;
        return 1;
    }
    if(col==1) { *value=(LIS_REAL)-1.0*s; return 1; }
    if(col==2) { *value=(LIS_REAL)2.0*s; return 1; }
    return 0;
}

static int
check_scaled_nested_csr(
    LIS_MATRIX csr,
    const TEST_ROW_CTX *actx,
    const TEST_ROW_CTX *bctx,
    const LIS_REAL d[3],
    LIS_REAL alpha,
    LIS_REAL beta)
{
    LIS_INT row,col,p,end;

    if(csr==NULL || csr->matrix_type!=LIS_MATRIX_CSR)
        return 1;

    for(row=0;row<3;row++)
    {
        p=csr->ptr[row];
        end=csr->ptr[row+1];

        for(col=0;col<3;col++)
        {
            LIS_REAL av=(LIS_REAL)0.0;
            LIS_REAL bv=(LIS_REAL)0.0;
            LIS_REAL expected=(LIS_REAL)0.0;
            int have_a=test_context_entry(actx,row,col,&av);
            int have_b=test_context_entry(bctx,row,col,&bv);

            if(!have_a && !have_b)
                continue;

            if(p>=end || csr->index[p]!=col)
                return 1;

            if(have_a)
                expected+=alpha*d[row]*av*d[col];
            if(have_b)
                expected+=beta*bv;

            if(!scalar_equal_real(&csr->value[p],expected))
                return 1;
            p++;
        }

        if(p!=end)
            return 1;
    }

    return 0;
}


static int
test_auto_matrix_policy(LIS_MATRIX A)
{
    LIS_SOLVER direct_solver=NULL,cached_solver=NULL;
    LIS_MATRIX_CSR_CACHE cache=NULL;
    LIS_MATRIX selected=NULL,first=NULL;
    LIS_INT policy=-1,action=-99,err;
    int failed=0;

    err=lis_solver_create(&direct_solver);
    if(err) { failed=1; goto done; }
    err=lis_solver_set_option("-i gmres -p ilut",direct_solver);
    if(err) { failed=1; goto done; }

    err=lis_matrix_prepare_solver_matrix(
        A,direct_solver,&cache,&selected,&policy,&action);
    if(err || selected!=A || policy!=LIS_MATRIX_POLICY_DIRECT ||
       action!=LIS_MATRIX_POLICY_ACTION_DIRECT)
    {
        fprintf(stderr,"AUTO policy did not retain a row-backed USER matrix for ILUT\n");
        failed=1; goto done;
    }

    err=lis_solver_create(&cached_solver);
    if(err) { failed=1; goto done; }
    err=lis_solver_set_option("-i gmres -p ssor",cached_solver);
    if(err) { failed=1; goto done; }

    err=lis_matrix_prepare_solver_matrix(
        A,cached_solver,&cache,&selected,&policy,&action);
    if(err || selected==A || selected==NULL ||
       selected->matrix_type!=LIS_MATRIX_CSR ||
       policy!=LIS_MATRIX_POLICY_CACHED_CSR ||
       action!=LIS_MATRIX_MATERIALIZE_REBUILD)
    {
        fprintf(stderr,"AUTO policy did not materialize CSR for SSOR\n");
        failed=1; goto done;
    }
    first=selected;

    err=lis_matrix_prepare_solver_matrix(
        A,cached_solver,&cache,&selected,&policy,&action);
    if(err || selected!=first || action!=LIS_MATRIX_MATERIALIZE_REUSE)
    {
        fprintf(stderr,"AUTO policy did not reuse unchanged CSR\n");
        failed=1; goto done;
    }

    err=lis_matrix_mark_values_changed(A);
    if(err) { failed=1; goto done; }
    err=lis_matrix_prepare_solver_matrix(
        A,cached_solver,&cache,&selected,&policy,&action);
    if(err || selected!=first || action!=LIS_MATRIX_MATERIALIZE_VALUES)
    {
        fprintf(stderr,"AUTO policy did not refresh CSR values in place\n");
        failed=1; goto done;
    }

    err=lis_matrix_mark_structure_changed(A);
    if(err) { failed=1; goto done; }
    err=lis_matrix_prepare_solver_matrix(
        A,cached_solver,&cache,&selected,&policy,&action);
    if(err || selected==first || action!=LIS_MATRIX_MATERIALIZE_REBUILD)
    {
        fprintf(stderr,"AUTO policy did not rebuild CSR after structure generation change\n");
        failed=1; goto done;
    }

done:
    lis_matrix_csr_cache_destroy(cache);
    if(cached_solver) lis_solver_destroy(cached_solver);
    if(direct_solver) lis_solver_destroy(direct_solver);
    return failed;
}

int main(int argc,char **argv)
{
    TEST_ROW_CTX actx={0,(LIS_REAL)1.0};
    TEST_ROW_CTX bctx={1,(LIS_REAL)1.0};
    LIS_REAL d[3]={(LIS_REAL)2.0,(LIS_REAL)3.0,(LIS_REAL)4.0};
    LIS_MATRIX A=NULL,B=NULL,op=NULL,scaled=NULL,nested=NULL;
    LIS_MATRIX csr=NULL,previous=NULL;
    LIS_MATRIX_CSR_CACHE acache=NULL,opcache=NULL,nestedcache=NULL;
    LIS_VECTOR scalevec=NULL;
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

    /* Build D*A*D inside another linear combination. */
    err=lis_vector_duplicate(A,&scalevec); if(err) { failed=1; goto cleanup; }
    scalar_set_real(&scalevec->value[0],d[0]);
    scalar_set_real(&scalevec->value[1],d[1]);
    scalar_set_real(&scalevec->value[2],d[2]);
    err=lis_matrix_create_scaled_operator(A,scalevec,&scaled);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_create_operator(
        (LIS_SCALAR)0.5,scaled,(LIS_SCALAR)2.0,B,&nested);
    if(err) { failed=1; goto cleanup; }

    err=lis_matrix_materialize_csr_cached(
        nested,&nestedcache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD ||
       check_scaled_nested_csr(
           csr,&actx,&bctx,d,(LIS_REAL)0.5,(LIS_REAL)2.0))
    {
        fprintf(stderr,"scaled nested operator materialization failed\n");
        failed=1; goto cleanup;
    }
    previous=csr;

    /* A grandchild value generation must refresh values without reallocating. */
    actx.scale=(LIS_REAL)3.0;
    err=lis_matrix_mark_values_changed(A);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_materialize_csr_cached(
        nested,&nestedcache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES || csr!=previous ||
       check_scaled_nested_csr(
           csr,&actx,&bctx,d,(LIS_REAL)0.5,(LIS_REAL)2.0))
    {
        fprintf(stderr,"nested grandchild value refresh failed\n");
        failed=1; goto cleanup;
    }

    err=lis_matrix_operator_set_coefficients(
        nested,(LIS_SCALAR)1.0,(LIS_SCALAR)-1.0);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_materialize_csr_cached(
        nested,&nestedcache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES || csr!=previous ||
       check_scaled_nested_csr(
           csr,&actx,&bctx,d,(LIS_REAL)1.0,(LIS_REAL)-1.0))
    {
        fprintf(stderr,"operator coefficient cache refresh failed\n");
        failed=1; goto cleanup;
    }

    d[1]=(LIS_REAL)5.0;
    scalar_set_real(&scalevec->value[1],d[1]);
    err=lis_matrix_mark_values_changed(scaled);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_materialize_csr_cached(
        nested,&nestedcache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES || csr!=previous ||
       check_scaled_nested_csr(
           csr,&actx,&bctx,d,(LIS_REAL)1.0,(LIS_REAL)-1.0))
    {
        fprintf(stderr,"scaled operator cache refresh failed\n");
        failed=1; goto cleanup;
    }

    /* Structural change below the scaled child must force a rebuild. */
    actx.variant=0;
    err=lis_matrix_mark_structure_changed(A);
    if(err) { failed=1; goto cleanup; }
    err=lis_matrix_materialize_csr_cached(
        nested,&nestedcache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD || csr==previous ||
       check_scaled_nested_csr(
           csr,&actx,&bctx,d,(LIS_REAL)1.0,(LIS_REAL)-1.0))
    {
        fprintf(stderr,"scaled nested structural rebuild failed\n");
        failed=1; goto cleanup;
    }

    if(test_auto_matrix_policy(A))
    {
        failed=1;
        goto cleanup;
    }

    printf("LIS_MATRIX_MATERIALIZE_CACHE PASSED\n");

cleanup:
    if(nested) lis_matrix_destroy(nested);
    if(scaled) lis_matrix_destroy(scaled);
    if(scalevec) lis_vector_destroy(scalevec);
    lis_matrix_csr_cache_destroy(nestedcache);
    lis_matrix_csr_cache_destroy(opcache);
    lis_matrix_csr_cache_destroy(acache);
    if(op) lis_matrix_destroy(op);
    if(B) lis_matrix_destroy(B);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed ? 1 : 0;
}
