#include "lis.h"
#include "lis_matrix_split_view.h"
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#ifdef _COMPLEX
#include <complex.h>
#endif

static int near_value(LIS_SCALAR value,double expected)
{
#ifdef _COMPLEX
    return cabs(value-(LIS_SCALAR)expected)<1.0e-12;
#else
    return fabs((double)value-expected)<1.0e-12;
#endif
}

int main(int argc,char **argv)
{
    const int64_t ptr[4]={1,2,3,3};
    const int64_t index[2]={2,3};
    const LIS_SCALAR diag[3]={
        (LIS_SCALAR)4.0,(LIS_SCALAR)5.0,(LIS_SCALAR)6.0};
    const LIS_SCALAR off[2]={
        (LIS_SCALAR)-1.0,(LIS_SCALAR)-2.0};
    const int64_t ptr2[4]={0,0,1,1};
    const int64_t index2[1]={2};
    const LIS_SCALAR diag2[3]={
        (LIS_SCALAR)2.0,(LIS_SCALAR)3.0,(LIS_SCALAR)4.0};
    const LIS_SCALAR lower2[1]={(LIS_SCALAR)5.0};
    const LIS_SCALAR upper2[1]={(LIS_SCALAR)-2.0};
    const int64_t bad_ptr[4]={0,1,1,1};
    const int64_t bad_index[1]={INT64_MAX};
    const int64_t nonmonotone_ptr[4]={0,1,0,1};
    const int64_t wrong_end_ptr[4]={0,0,1,2};
    const int64_t wrong_index[1]={0};
    LIS_MATRIX A=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_UNSIGNED_INT s0=0,v0=0,s1=0,v1=0,s2=0,v2=0;
    LIS_INT err;
    int failed=0;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=lis_matrix_create_split_triangular_view_i64(
        LIS_COMM_WORLD,3,2,ptr,index,diag,off,NULL,1,&A,&view);
    if(err)
    {
        fprintf(stderr,"64-bit split-view import failed\n");
        failed=1;
        goto cleanup;
    }

    err=lis_vector_duplicate(A,&x); if(err) { failed=1; goto cleanup; }
    err=lis_vector_duplicate(A,&y); if(err) { failed=1; goto cleanup; }

    x->value[0]=(LIS_SCALAR)1.0;
    x->value[1]=(LIS_SCALAR)2.0;
    x->value[2]=(LIS_SCALAR)3.0;

    err=lis_matvec(A,x,y);
    if(err || !near_value(y->value[0],2.0) ||
       !near_value(y->value[1],3.0) ||
       !near_value(y->value[2],14.0))
    {
        fprintf(stderr,"64-bit split-view matvec mismatch\n");
        failed=1;
        goto cleanup;
    }

    err=lis_matrix_get_generation(A,&s0,&v0);
    if(err) { failed=1; goto cleanup; }

    /* Change nnz, origin and symmetry mode through the fixed-width API. */
    err=lis_matrix_split_triangular_view_rebind_i64(
        view,1,ptr2,index2,diag2,lower2,upper2,0);
    if(err)
    {
        fprintf(stderr,"64-bit split-view rebind failed\n");
        failed=1;
        goto cleanup;
    }

    err=lis_matrix_get_generation(A,&s1,&v1);
    if(err || s1!=s0+1 || v1!=v0+1)
    {
        fprintf(stderr,"64-bit split-view rebind generation mismatch\n");
        failed=1;
        goto cleanup;
    }

    err=lis_matvec(A,x,y);
    if(err || !near_value(y->value[0],2.0) ||
       !near_value(y->value[1],0.0) ||
       !near_value(y->value[2],22.0))
    {
        fprintf(stderr,"64-bit split-view rebind matvec mismatch\n");
        failed=1;
        goto cleanup;
    }

    /* Invalid fixed-width structure must roll back without a generation bump. */
    err=lis_matrix_split_triangular_view_rebind_i64(
        view,1,bad_ptr,bad_index,diag2,lower2,upper2,0);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"64-bit split-view invalid rebind was accepted\n");
        failed=1;
        goto cleanup;
    }

    if(lis_matrix_split_triangular_view_rebind_i64(
           view,1,nonmonotone_ptr,index2,diag2,lower2,upper2,0)==LIS_SUCCESS ||
       lis_matrix_split_triangular_view_rebind_i64(
           view,1,wrong_end_ptr,index2,diag2,lower2,upper2,0)==LIS_SUCCESS ||
       lis_matrix_split_triangular_view_rebind_i64(
           view,1,ptr2,wrong_index,diag2,lower2,upper2,0)==LIS_SUCCESS)
    {
        fprintf(stderr,"invalid 64-bit column structure was accepted\n");
        failed=1; goto cleanup;
    }

    err=lis_matrix_get_generation(A,&s2,&v2);
    if(err || s2!=s1 || v2!=v1)
    {
        fprintf(stderr,"64-bit split-view failed rebind changed generation\n");
        failed=1;
        goto cleanup;
    }

    err=lis_matvec(A,x,y);
    if(err || !near_value(y->value[0],2.0) ||
       !near_value(y->value[1],0.0) ||
       !near_value(y->value[2],22.0))
    {
        fprintf(stderr,"64-bit split-view rollback mismatch\n");
        failed=1;
        goto cleanup;
    }

    /* Null upper values select the documented symmetric representation. */
    err=lis_matrix_split_triangular_view_rebind_i64(
        view,1,ptr2,index2,diag2,lower2,NULL,0);
    if(err || lis_matvec(A,x,y) ||
       !near_value(y->value[0],2.0) ||
       !near_value(y->value[1],21.0) ||
       !near_value(y->value[2],22.0))
    {
        fprintf(stderr,"symmetric 64-bit rebind with NULL upper values failed\n");
        failed=1; goto cleanup;
    }

    /* Reject invalid API inputs without changing the active factor view. */
    if(lis_matrix_split_triangular_view_rebind_i64(
           view,-1,ptr2,index2,diag2,lower2,NULL,0)==LIS_SUCCESS ||
       lis_matrix_split_triangular_view_rebind_i64(
           view,1,NULL,index2,diag2,lower2,NULL,0)==LIS_SUCCESS ||
       lis_matrix_split_triangular_view_rebind_i64(
           view,1,ptr2,NULL,diag2,lower2,NULL,0)==LIS_SUCCESS ||
       lis_matrix_split_triangular_view_rebind_i64(
           view,1,ptr2,index2,NULL,lower2,NULL,0)==LIS_SUCCESS ||
       lis_matrix_split_triangular_view_rebind_i64(
           view,1,ptr2,index2,diag2,lower2,NULL,2)==LIS_SUCCESS)
    {
        fprintf(stderr,"invalid 64-bit API arguments were accepted\n");
        failed=1; goto cleanup;
    }
    if(lis_matvec(A,x,y) ||
       !near_value(y->value[1],21.0) ||
       !near_value(y->value[2],22.0))
    {
        fprintf(stderr,"invalid 64-bit rebind damaged active view\n");
        failed=1; goto cleanup;
    }

    if(sizeof(LIS_INT)<sizeof(int64_t))
    {
        const int64_t dummy_ptr[1]={0};
        const LIS_SCALAR dummy_diag[1]={(LIS_SCALAR)1.0};
        LIS_MATRIX too_large_A=NULL;
        LIS_MATRIX_SPLIT_TRIANGULAR_VIEW too_large_view=NULL;
        int64_t too_large_n=(int64_t)INT_MAX+INT64_C(1);

        err=lis_matrix_create_split_triangular_view_i64(
            LIS_COMM_WORLD,too_large_n,0,dummy_ptr,NULL,
            dummy_diag,NULL,NULL,0,&too_large_A,&too_large_view);
        if(err==LIS_SUCCESS || too_large_A!=NULL || too_large_view!=NULL)
        {
            fprintf(stderr,"64-bit split-view range rejection failed\n");
            if(too_large_A) lis_matrix_destroy(too_large_A);
            lis_matrix_split_triangular_view_destroy(too_large_view);
            failed=1;
            goto cleanup;
        }
    }

    printf("LIS_SPLIT_TRIANGULAR_VIEW_I64 PASSED\n");

cleanup:
    if(y) lis_vector_destroy(y);
    if(x) lis_vector_destroy(x);
    if(A) lis_matrix_destroy(A);
    lis_matrix_split_triangular_view_destroy(view);
    lis_finalize();
    return failed ? 1 : 0;
}
