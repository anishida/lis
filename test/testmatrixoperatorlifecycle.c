#include "lis.h"
#include "lis_matrix_operator.h"
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

static LIS_INT make_matrix(
    LIS_INT n,
    LIS_SCALAR a00,LIS_SCALAR a01,
    LIS_SCALAR a10,LIS_SCALAR a11,
    LIS_MATRIX *A)
{
    LIS_INT err;

    *A=NULL;
    err=lis_matrix_create(LIS_COMM_WORLD,A); if(err) return err;
    err=lis_matrix_set_size(*A,0,n); if(err) goto fail;

    if(n==2)
    {
        err=lis_matrix_set_value(LIS_INS_VALUE,0,0,a00,*A); if(err) goto fail;
        err=lis_matrix_set_value(LIS_INS_VALUE,0,1,a01,*A); if(err) goto fail;
        err=lis_matrix_set_value(LIS_INS_VALUE,1,0,a10,*A); if(err) goto fail;
        err=lis_matrix_set_value(LIS_INS_VALUE,1,1,a11,*A); if(err) goto fail;
    }
    else
    {
        LIS_INT i;
        for(i=0;i<n;i++)
        {
            err=lis_matrix_set_value(
                LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,*A);
            if(err) goto fail;
        }
    }

    err=lis_matrix_assemble(*A); if(err) goto fail;
    return LIS_SUCCESS;

fail:
    if(*A) lis_matrix_destroy(*A);
    *A=NULL;
    return err;
}

static int check_y(LIS_VECTOR y,double e0,double e1)
{
    return !near_value(y->value[0],e0) || !near_value(y->value[1],e1);
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL,B=NULL,C=NULL,D=NULL,E=NULL,op=NULL,nested=NULL;
    LIS_VECTOR x=NULL,y=NULL;
    LIS_UNSIGNED_INT s0=0,v0=0,s1=0,v1=0,s2=0,v2=0,s3=0,v3=0;
    LIS_INT err;
    int failed=0;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=make_matrix(2,2,0,0,3,&A); if(err) { failed=1; goto done; }
    err=make_matrix(2,1,1,0,2,&B); if(err) { failed=1; goto done; }
    err=make_matrix(2,4,0,1,1,&C); if(err) { failed=1; goto done; }
    err=make_matrix(2,1,0,0,2,&D); if(err) { failed=1; goto done; }
    err=make_matrix(3,0,0,0,0,&E); if(err) { failed=1; goto done; }

    err=lis_matrix_create_operator(
        (LIS_SCALAR)1.0,A,(LIS_SCALAR)1.0,B,&op);
    if(err) { failed=1; goto done; }
    err=lis_matrix_create_operator(
        (LIS_SCALAR)0.5,op,(LIS_SCALAR)1.0,B,&nested);
    if(err) { failed=1; goto done; }

    err=lis_vector_duplicate(op,&x); if(err) { failed=1; goto done; }
    err=lis_vector_duplicate(op,&y); if(err) { failed=1; goto done; }
    x->value[0]=(LIS_SCALAR)1.0;
    x->value[1]=(LIS_SCALAR)2.0;

    err=lis_matvec(op,x,y);
    if(err || check_y(y,5.0,10.0))
    {
        fprintf(stderr,"initial operator matvec mismatch\n");
        failed=1; goto done;
    }

    err=lis_matrix_get_generation(op,&s0,&v0);
    if(err) { failed=1; goto done; }
    err=lis_matrix_operator_set_coefficients(
        op,(LIS_SCALAR)2.0,(LIS_SCALAR)-1.0);
    if(err) { failed=1; goto done; }
    err=lis_matrix_get_generation(op,&s1,&v1);
    if(err || s1!=s0 || v1!=v0+1)
    {
        fprintf(stderr,"operator coefficient generation mismatch\n");
        failed=1; goto done;
    }

    err=lis_matvec(op,x,y);
    if(err || check_y(y,1.0,8.0))
    {
        fprintf(stderr,"updated operator matvec mismatch\n");
        failed=1; goto done;
    }

    /* A parent operator must observe the mutable child without reconstruction. */
    err=lis_matvec(nested,x,y);
    if(err || check_y(y,3.5,8.0))
    {
        fprintf(stderr,"nested coefficient update mismatch\n");
        failed=1; goto done;
    }

    err=lis_matrix_operator_rebind(op,C,D);
    if(err) { failed=1; goto done; }
    err=lis_matrix_get_generation(op,&s2,&v2);
    if(err || s2!=s1+1 || v2!=v1+1)
    {
        fprintf(stderr,"operator rebind generation mismatch\n");
        failed=1; goto done;
    }

    err=lis_matvec(op,x,y);
    if(err || check_y(y,7.0,2.0))
    {
        fprintf(stderr,"rebound operator matvec mismatch\n");
        failed=1; goto done;
    }

    err=lis_matvec(nested,x,y);
    if(err || check_y(y,6.5,5.0))
    {
        fprintf(stderr,"nested rebind update mismatch\n");
        failed=1; goto done;
    }

    /* Self-reference and indirect cycles must be rejected transactionally. */
    err=lis_matrix_operator_rebind(op,op,D);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"self-referential rebind was accepted\n");
        failed=1; goto done;
    }
    err=lis_matrix_operator_rebind(op,nested,D);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"indirect cyclic rebind was accepted\n");
        failed=1; goto done;
    }
    err=lis_matrix_get_generation(op,&s3,&v3);
    if(err || s3!=s2 || v3!=v2)
    {
        fprintf(stderr,"cyclic rebind changed matrix generation\n");
        failed=1; goto done;
    }

    /* Incompatible dimensions must leave the active operator untouched. */
    err=lis_matrix_operator_rebind(op,E,D);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"invalid operator rebind was accepted\n");
        failed=1; goto done;
    }
    err=lis_matrix_get_generation(op,&s3,&v3);
    if(err || s3!=s2 || v3!=v2)
    {
        fprintf(stderr,"failed operator rebind changed generation\n");
        failed=1; goto done;
    }

    err=lis_matvec(op,x,y);
    if(err || check_y(y,7.0,2.0))
    {
        fprintf(stderr,"operator rollback mismatch\n");
        failed=1; goto done;
    }
    err=lis_matvec(nested,x,y);
    if(err || check_y(y,6.5,5.0))
    {
        fprintf(stderr,"nested rollback mismatch\n");
        failed=1; goto done;
    }

    printf("LIS_MATRIX_OPERATOR_LIFECYCLE PASSED\n");

done:
    if(y) lis_vector_destroy(y);
    if(x) lis_vector_destroy(x);
    if(nested) lis_matrix_destroy(nested);
    if(op) lis_matrix_destroy(op);
    if(E) lis_matrix_destroy(E);
    if(D) lis_matrix_destroy(D);
    if(C) lis_matrix_destroy(C);
    if(B) lis_matrix_destroy(B);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return failed ? 1 : 0;
}
