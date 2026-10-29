#include "lis.h"
#include "lis_matrix_split_view.h"
#include <math.h>
#include <stdio.h>
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

static int near_scalar(LIS_SCALAR value,LIS_SCALAR expected)
{
#ifdef _COMPLEX
    return cabs(value-expected)<1.0e-12;
#else
    return fabs((double)(value-expected))<1.0e-12;
#endif
}

static int set_x123(LIS_VECTOR x)
{
    if(x==NULL || x->n!=3) return 1;
    x->value[0]=(LIS_SCALAR)1.0;
    x->value[1]=(LIS_SCALAR)2.0;
    x->value[2]=(LIS_SCALAR)3.0;
    return 0;
}

static int check_y(LIS_VECTOR y,double y0,double y1,double y2)
{
    return !near_value(y->value[0],y0) ||
           !near_value(y->value[1],y1) ||
           !near_value(y->value[2],y2);
}

static int check_middle_row_sorted(LIS_MATRIX A)
{
    LIS_INT nnz=3;
    LIS_INT cols[3]={-1,-1,-1};
    LIS_SCALAR vals[3];
    LIS_INT err;

    err=lis_matrix_get_user_row(A,1,&nnz,cols,vals);
    if(err || nnz!=3) return 1;
    if(cols[0]!=0 || cols[1]!=1 || cols[2]!=2) return 1;
    return 0;
}

static int run_symmetric_origin1(void)
{
    const LIS_INT ptr[4]={1,2,3,3};
    const LIS_INT index[2]={2,3};
    const LIS_SCALAR diag0[3]={(LIS_SCALAR)4.0,(LIS_SCALAR)5.0,(LIS_SCALAR)6.0};
    const LIS_SCALAR off0[2]={(LIS_SCALAR)-1.0,(LIS_SCALAR)-2.0};
    const LIS_SCALAR diag1[3]={(LIS_SCALAR)5.0,(LIS_SCALAR)6.0,(LIS_SCALAR)7.0};
    const LIS_SCALAR off1[2]={(LIS_SCALAR)-0.5,(LIS_SCALAR)-1.0};
    LIS_MATRIX A=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_INT err;
    int rc=1;

    err=lis_matrix_create_split_triangular_view(
        LIS_COMM_WORLD,3,2,ptr,index,diag0,off0,NULL,1,&A,&view);
    if(err) goto done;
    if(check_middle_row_sorted(A)) goto done;

    err=lis_vector_duplicate(A,&x); if(err) goto done;
    err=lis_vector_duplicate(A,&y); if(err) goto done;
    if(set_x123(x)) goto done;

    err=lis_matvec(A,x,y); if(err) goto done;
    if(check_y(y,2.0,3.0,14.0)) goto done;
    err=lis_matvech(A,x,y); if(err) goto done;
    if(check_y(y,2.0,3.0,14.0)) goto done;

    err=lis_matrix_split_triangular_view_set_values(view,diag1,off1,NULL);
    if(err) goto done;
    err=lis_matvec(A,x,y); if(err) goto done;
    if(check_y(y,4.0,8.5,19.0)) goto done;
    err=lis_matvech(A,x,y); if(err) goto done;
    if(check_y(y,4.0,8.5,19.0)) goto done;

    rc=0;
done:
    if(y) lis_vector_destroy(y);
    if(x) lis_vector_destroy(x);
    if(A) lis_matrix_destroy(A);
    lis_matrix_split_triangular_view_destroy(view);
    return rc;
}

static int run_nonsymmetric_origin0(void)
{
    const LIS_INT ptr[4]={0,1,2,2};
    const LIS_INT index[2]={1,2};
    const LIS_SCALAR diag0[3]={(LIS_SCALAR)4.0,(LIS_SCALAR)5.0,(LIS_SCALAR)6.0};
    const LIS_SCALAR lower0[2]={(LIS_SCALAR)-1.0,(LIS_SCALAR)-2.0};
    const LIS_SCALAR upper0[2]={(LIS_SCALAR)-3.0,(LIS_SCALAR)-4.0};
    const LIS_SCALAR diag1[3]={(LIS_SCALAR)5.0,(LIS_SCALAR)6.0,(LIS_SCALAR)7.0};
    const LIS_SCALAR lower1[2]={(LIS_SCALAR)-0.5,(LIS_SCALAR)-1.0};
    const LIS_SCALAR upper1[2]={(LIS_SCALAR)-1.5,(LIS_SCALAR)-2.0};
    LIS_MATRIX A=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_INT err;
    int rc=1;

    err=lis_matrix_create_split_triangular_view(
        LIS_COMM_WORLD,3,2,ptr,index,diag0,lower0,upper0,0,&A,&view);
    if(err) goto done;
    if(check_middle_row_sorted(A)) goto done;

    err=lis_vector_duplicate(A,&x); if(err) goto done;
    err=lis_vector_duplicate(A,&y); if(err) goto done;
    if(set_x123(x)) goto done;

    err=lis_matvec(A,x,y); if(err) goto done;
    if(check_y(y,-2.0,-3.0,14.0)) goto done;
    err=lis_matvech(A,x,y); if(err) goto done;
    if(check_y(y,2.0,1.0,10.0)) goto done;

    err=lis_matrix_split_triangular_view_set_values(
        view,diag1,lower1,upper1);
    if(err) goto done;
    err=lis_matvec(A,x,y); if(err) goto done;
    if(check_y(y,2.0,5.5,19.0)) goto done;
    err=lis_matvech(A,x,y); if(err) goto done;
    if(check_y(y,4.0,7.5,17.0)) goto done;

    /* Symmetry mode is structural and cannot change during a value update. */
    err=lis_matrix_split_triangular_view_set_values(view,diag1,lower1,NULL);
    if(err==LIS_SUCCESS) goto done;

    rc=0;
done:
    if(y) lis_vector_destroy(y);
    if(x) lis_vector_destroy(x);
    if(A) lis_matrix_destroy(A);
    lis_matrix_split_triangular_view_destroy(view);
    return rc;
}

#ifdef _COMPLEX
static int run_complex_adjoint(void)
{
    const LIS_INT ptr[3]={0,1,1};
    const LIS_INT index[1]={1};
    const LIS_SCALAR diag[2]={1.0+2.0*I,2.0-1.0*I};
    const LIS_SCALAR lower[1]={3.0+2.0*I};
    const LIS_SCALAR upper[1]={-1.0+4.0*I};
    LIS_MATRIX A=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_SCALAR e0,e1;
    LIS_INT err;
    int rc=1;

    err=lis_matrix_create_split_triangular_view(
        LIS_COMM_WORLD,2,1,ptr,index,diag,lower,upper,0,&A,&view);
    if(err) goto done;
    err=lis_vector_duplicate(A,&x); if(err) goto done;
    err=lis_vector_duplicate(A,&y); if(err) goto done;

    x->value[0]=1.0+2.0*I;
    x->value[1]=-2.0+1.0*I;

    e0=conj(diag[0])*x->value[0]+conj(lower[0])*x->value[1];
    e1=conj(upper[0])*x->value[0]+conj(diag[1])*x->value[1];

    err=lis_matvech(A,x,y); if(err) goto done;
    if(!near_scalar(y->value[0],e0) || !near_scalar(y->value[1],e1))
        goto done;

    rc=0;
done:
    if(y) lis_vector_destroy(y);
    if(x) lis_vector_destroy(x);
    if(A) lis_matrix_destroy(A);
    lis_matrix_split_triangular_view_destroy(view);
    return rc;
}
#endif

static int reject_invalid_storage(void)
{
    const LIS_SCALAR diag[3]={(LIS_SCALAR)1.0,(LIS_SCALAR)1.0,(LIS_SCALAR)1.0};
    const LIS_SCALAR lower[2]={(LIS_SCALAR)1.0,(LIS_SCALAR)1.0};
    const LIS_INT ptr_diag[4]={0,1,1,1};
    const LIS_INT index_diag[1]={0};
    const LIS_INT ptr_unsorted[4]={0,2,2,2};
    const LIS_INT index_unsorted[2]={2,1};
    LIS_MATRIX A=NULL;
    LIS_MATRIX_SPLIT_TRIANGULAR_VIEW view=NULL;
    LIS_INT err;

    err=lis_matrix_create_split_triangular_view(
        LIS_COMM_WORLD,3,1,ptr_diag,index_diag,diag,lower,NULL,0,&A,&view);
    if(err==LIS_SUCCESS || A!=NULL || view!=NULL) return 1;

    err=lis_matrix_create_split_triangular_view(
        LIS_COMM_WORLD,3,2,ptr_unsorted,index_unsorted,diag,lower,NULL,0,&A,&view);
    if(err==LIS_SUCCESS || A!=NULL || view!=NULL) return 1;

    return 0;
}

int main(int argc,char **argv)
{
    LIS_INT err;
    int rc=1;

    err=lis_initialize(&argc,&argv); if(err) return 1;

    if(run_symmetric_origin1()) goto done;
    if(run_nonsymmetric_origin0()) goto done;
#ifdef _COMPLEX
    if(run_complex_adjoint()) goto done;
#endif
    if(reject_invalid_storage()) goto done;

    printf("LIS_SPLIT_TRIANGULAR_VIEW PASSED\n");
    rc=0;

done:
    lis_finalize();
    return rc;
}
