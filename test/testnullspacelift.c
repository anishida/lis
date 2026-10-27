#include "lis.h"
#include "lis_nullspace_lift.h"
#include <math.h>
#include <stdio.h>
#ifdef _COMPLEX
#include <complex.h>
#endif

static int set_global_pair(LIS_VECTOR v,LIS_SCALAR value)
{
    LIS_INT is,ie,i,err;

    err=lis_vector_get_range(v,&is,&ie);
    if(err) return 1;

    for(i=is;i<ie;i++)
    {
        if(i==0 || i==1)
            v->value[i-is]=value;
    }
    return 0;
}

static int subtract_global_pair(LIS_VECTOR v,LIS_SCALAR value)
{
    LIS_INT is,ie,i,err;

    err=lis_vector_get_range(v,&is,&ie);
    if(err) return 1;

    for(i=is;i<ie;i++)
    {
        if(i==0 || i==1)
            v->value[i-is]-=value;
    }
    return 0;
}

static int build_basis(LIS_MATRIX A,LIS_VECTOR *q0,LIS_VECTOR *q1)
{
    LIS_INT is,ie,i,err;

    err=lis_vector_duplicate(A,q0); if(err) return 1;
    err=lis_vector_duplicate(A,q1); if(err) return 1;
    err=lis_vector_set_all((LIS_SCALAR)0.0,*q0); if(err) return 1;
    err=lis_vector_set_all((LIS_SCALAR)0.0,*q1); if(err) return 1;

    err=lis_vector_get_range(*q0,&is,&ie); if(err) return 1;
    for(i=is;i<ie;i++)
    {
        if(i==0 || i==1)
            (*q0)->value[i-is]=(LIS_SCALAR)1.0;
        if(i==2 || i==3)
            (*q1)->value[i-is]=(LIS_SCALAR)1.0;
    }
    return 0;
}

static int check_lift(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR lift,
    LIS_SCALAR rhs_value,
    LIS_SCALAR expected_value)
{
    LIS_NULLSPACE_LIFT_DIAGNOSTICS diag;
    LIS_REAL nrm;
    LIS_INT err;

    err=lis_vector_set_all((LIS_SCALAR)0.0,b); if(err) return 1;
    if(set_global_pair(b,rhs_value)) return 1;

    err=lis_solver_build_nullspace_lift(
        solver,b,(LIS_REAL)2.0,(LIS_REAL)1.0e-12,lift,&diag);
    if(err) return 1;

    if(diag.nullspace_dim!=2 || diag.active_count!=1)
        return 1;

    if(subtract_global_pair(lift,expected_value)) return 1;
    err=lis_vector_nrm2(lift,&nrm);
    if(err || nrm>(LIS_REAL)1.0e-11)
        return 1;

    return 0;
}

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR q0=NULL,q1=NULL,b=NULL,lift=NULL;
    LIS_VECTOR basis[2];
    LIS_SOLVER solver=NULL;
    LIS_INT err;
    int rc=1;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,4);
    if(err) goto done;

    if(build_basis(A,&q0,&q1)) goto done;
    err=lis_vector_duplicate(A,&b); if(err) goto done;
    err=lis_vector_duplicate(A,&lift); if(err) goto done;

    basis[0]=q0; basis[1]=q1;
    err=lis_solver_create(&solver); if(err) goto done;
    err=lis_solver_set_nullspace(solver,2,basis); if(err) goto done;

    if(check_lift(
        solver,b,lift,
        (LIS_SCALAR)2.0,
        (LIS_SCALAR)1.0))
        goto done;

#ifdef _COMPLEX
    /* Exercise the Hermitian q^H b path with a purely imaginary component. */
    if(check_lift(
        solver,b,lift,
        (LIS_SCALAR)(2.0*I),
        (LIS_SCALAR)(1.0*I)))
        goto done;
#endif

    printf("LIS_NULLSPACE_LIFT PASSED\n");
    rc=0;

done:
    if(solver) lis_solver_destroy(solver);
    if(lift) lis_vector_destroy(lift);
    if(b) lis_vector_destroy(b);
    if(q1) lis_vector_destroy(q1);
    if(q0) lis_vector_destroy(q0);
    if(A) lis_matrix_destroy(A);
    lis_finalize();
    return rc;
}
