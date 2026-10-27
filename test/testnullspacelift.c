#include "lis.h"
#include "lis_nullspace_lift.h"
#include <math.h>
#include <stdio.h>

int main(int argc,char **argv)
{
    LIS_MATRIX A=NULL;
    LIS_VECTOR q0=NULL,q1=NULL,b=NULL,lift=NULL;
    LIS_VECTOR basis[2];
    LIS_SOLVER solver=NULL;
    LIS_NULLSPACE_LIFT_DIAGNOSTICS diag;
    LIS_REAL nrm;
    LIS_INT err;
    int rc=1;

    err=lis_initialize(&argc,&argv);
    if(err) return 1;

    err=lis_matrix_create(LIS_COMM_WORLD,&A);
    if(err) goto done;
    err=lis_matrix_set_size(A,0,4);
    if(err) goto done;

    err=lis_vector_duplicate(A,&q0); if(err) goto done;
    err=lis_vector_duplicate(A,&q1); if(err) goto done;
    err=lis_vector_duplicate(A,&b); if(err) goto done;
    err=lis_vector_duplicate(A,&lift); if(err) goto done;

    lis_vector_set_all((LIS_SCALAR)0.0,q0);
    lis_vector_set_all((LIS_SCALAR)0.0,q1);
    lis_vector_set_all((LIS_SCALAR)0.0,b);

    q0->value[0]=(LIS_SCALAR)1.0;
    q0->value[1]=(LIS_SCALAR)1.0;
    q1->value[2]=(LIS_SCALAR)1.0;
    q1->value[3]=(LIS_SCALAR)1.0;
    b->value[0]=(LIS_SCALAR)2.0;
    b->value[1]=(LIS_SCALAR)2.0;

    basis[0]=q0; basis[1]=q1;
    err=lis_solver_create(&solver); if(err) goto done;
    err=lis_solver_set_nullspace(solver,2,basis); if(err) goto done;

    err=lis_solver_build_nullspace_lift(
        solver,b,(LIS_REAL)2.0,(LIS_REAL)1.0e-12,lift,&diag);
    if(err) goto done;

    if(diag.nullspace_dim!=2 || diag.active_count!=1)
        goto done;

    lift->value[0]-=(LIS_SCALAR)1.0;
    lift->value[1]-=(LIS_SCALAR)1.0;
    err=lis_vector_nrm2(lift,&nrm);
    if(err || nrm>(LIS_REAL)1.0e-11)
        goto done;

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
