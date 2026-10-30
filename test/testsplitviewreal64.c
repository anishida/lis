#include "lis.h"
#include "lis_matrix.h"
#include "lis_matrix_split_view.h"
#include "lis_matrix_materialize.h"
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#ifdef _COMPLEX
#include <complex.h>
#endif

static int near(LIS_SCALAR x,double y)
{
#ifdef _COMPLEX
    return cabs(x-(LIS_SCALAR)y)<1.0e-11;
#else
    return fabs((double)x-y)<1.0e-11;
#endif
}

static int check_vector(LIS_VECTOR x,const double want[3])
{
    return near(x->value[0],want[0]) &&
           near(x->value[1],want[1]) &&
           near(x->value[2],want[2]);
}

int main(int argc,char **argv)
{
    const int64_t ptr[4]={0,2,3,3};
    const int64_t idx[3]={1,2,2};
    double diag[3]={4,5,6};
    double lower[3]={1,2,3};
    double upper[3]={7,8,9};
    const LIS_SCALAR typed_diag[3]={
        (LIS_SCALAR)4,(LIS_SCALAR)5,(LIS_SCALAR)6};
    const LIS_SCALAR typed_lower[3]={
        (LIS_SCALAR)1,(LIS_SCALAR)2,(LIS_SCALAR)3};
    const LIS_SCALAR typed_upper[3]={
        (LIS_SCALAR)7,(LIS_SCALAR)8,(LIS_SCALAR)9};
    const double y_normal[3]={42,38,26};
    const double y_transpose[3]={12,26,44};
    const double y_symmetric[3]={12,20,26};
    const double y_changed[3]={48,38,26};
    const double y_restored[3]={43,38,26};
    LIS_MATRIX_SPLIT_TRIANGULAR_TOPOLOGY topology=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_SHARED_VIEW view=NULL,native_view=NULL,sym_view=NULL;
    LIS_MATRIX A=NULL,native=NULL,sym=NULL,csr=NULL,old_csr=NULL;
    LIS_MATRIX_CSR_CACHE cache=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_UNSIGNED_INT sg0=0,vg0=0,sg1=0,vg1=0;
    LIS_SCALAR diag_read[3],vals[3];
    LIS_INT cols[3],count,action,err;
    int failed=1;

    err=lis_initialize(&argc,&argv);
    if(err)return 1;

    err=lis_matrix_create_split_triangular_topology_i64(
        3,3,ptr,idx,0,&topology);
    if(err)goto done;

    err=lis_matrix_create_split_triangular_shared_view_real64(
        LIS_COMM_WORLD,topology,diag,lower,upper,&A,&view);
    if(err)goto done;
    err=lis_matrix_create_split_triangular_shared_view(
        LIS_COMM_WORLD,topology,typed_diag,typed_lower,typed_upper,
        &native,&native_view);
    if(err)goto done;
    err=lis_matrix_create_split_triangular_shared_view_real64(
        LIS_COMM_WORLD,topology,diag,lower,NULL,&sym,&sym_view);
    if(err)goto done;
    lis_matrix_split_triangular_topology_destroy(topology);
    topology=NULL;

    err=lis_vector_duplicate(A,&x);if(err)goto done;
    err=lis_vector_duplicate(A,&y);if(err)goto done;
    x->value[0]=(LIS_SCALAR)1;
    x->value[1]=(LIS_SCALAR)2;
    x->value[2]=(LIS_SCALAR)3;

    if(lis_matvec(A,x,y) || !check_vector(y,y_normal))
        goto done;
    if(lis_matvec(native,x,y) || !check_vector(y,y_normal))
        goto done;
    if(lis_matvech(A,x,y) || !check_vector(y,y_transpose))
        goto done;
    if(lis_matvech(native,x,y) || !check_vector(y,y_transpose))
        goto done;
    if(lis_matvec(sym,x,y) || !check_vector(y,y_symmetric))
        goto done;
    if(lis_matvech(sym,x,y) || !check_vector(y,y_symmetric))
        goto done;

    count=3;
    err=lis_matrix_get_user_row(A,0,&count,cols,vals);
    if(err || count!=3 ||
       cols[0]!=0 || cols[1]!=1 || cols[2]!=2 ||
       !near(vals[0],4) || !near(vals[1],7) || !near(vals[2],8))
        goto done;
    err=lis_matrix_get_user_row(A,1,&count,NULL,NULL);
    if(err || count!=3)goto done;
    err=A->user_get_diagonal(A->user_data,diag_read);
    if(err || !near(diag_read[0],4) || !near(diag_read[1],5) ||
       !near(diag_read[2],6))
        goto done;

    err=lis_matrix_materialize_csr_cached(A,&cache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REBUILD)goto done;
    old_csr=csr;
    err=lis_matrix_materialize_csr_cached(A,&cache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_REUSE || old_csr!=csr)
        goto done;

    err=lis_matrix_get_generation(A,&sg0,&vg0);
    if(err)goto done;
    diag[0]=10;
    if(lis_matrix_split_triangular_shared_view_set_values_real64(
           view,diag,lower,upper))goto done;
    err=lis_matrix_get_generation(A,&sg1,&vg1);
    if(err || sg1!=sg0 || vg1<=vg0)goto done;
    if(lis_matvec(A,x,y) || !check_vector(y,y_changed))goto done;
    err=lis_matrix_materialize_csr_cached(A,&cache,&csr,&action);
    if(err || action!=LIS_MATRIX_MATERIALIZE_VALUES || csr!=old_csr)
        goto done;
    if(lis_matvec(csr,x,y) || !check_vector(y,y_changed))goto done;

    /* Rebind to a different external array without allocating scalar copies. */
    {
        const double d2[3]={5,5,6};
        if(lis_matrix_split_triangular_shared_view_set_values_real64(
               view,d2,lower,upper))goto done;
        if(lis_matvec(A,x,y) || !check_vector(y,y_restored))goto done;
        if(lis_matrix_split_triangular_shared_view_set_values_real64(
               view,d2,lower,NULL)==LIS_SUCCESS)goto done;
        if(lis_matrix_split_triangular_shared_view_set_values(
               view,typed_diag,typed_lower,typed_upper)==LIS_SUCCESS)
            goto done;
        if(lis_matvec(A,x,y) || !check_vector(y,y_restored))goto done;
    }

    /* The native setter must also reject mismatched scalar storage. */
    if(lis_matrix_split_triangular_shared_view_set_values_real64(
           native_view,diag,lower,upper)==LIS_SUCCESS)goto done;

    failed=0;
done:
    if(failed)fprintf(stderr,"real64 shared view regression failed (err=%d)\n",(int)err);
    if(cache)lis_matrix_csr_cache_destroy(cache);
    if(y)lis_vector_destroy(y);
    if(x)lis_vector_destroy(x);
    if(A)lis_matrix_destroy(A);
    if(native)lis_matrix_destroy(native);
    if(sym)lis_matrix_destroy(sym);
    if(view)lis_matrix_split_triangular_shared_view_destroy(view);
    if(native_view)lis_matrix_split_triangular_shared_view_destroy(native_view);
    if(sym_view)lis_matrix_split_triangular_shared_view_destroy(sym_view);
    if(topology)lis_matrix_split_triangular_topology_destroy(topology);
    lis_finalize();
    printf("LIS_SPLIT_VIEW_REAL64 %s\n",failed?"FAILED":"PASSED");
    return failed;
}
