#include <math.h>
#include <stdio.h>
#ifdef _COMPLEX
#include <complex.h>
#endif
#include "lis.h"
#include "lis_nullspace_lift.h"

#define Q0_INDEX0 0
#define Q0_INDEX1 2
#define Q1_INDEX0 1
#define Q1_INDEX1 3

static int real_close(LIS_REAL actual,LIS_REAL expected,LIS_REAL tolerance)
{
    return fabs(actual-expected) <= tolerance*(1.0+fabs(expected));
}

static int set_global_entry(LIS_VECTOR v,LIS_INT index,LIS_SCALAR value)
{
    LIS_INT is,ie,err;

    err=lis_vector_get_range(v,&is,&ie);
    if(err) return 1;

    if(index>=is && index<ie)
        v->value[index-is]=value;
    return 0;
}

static int set_global_pair(
    LIS_VECTOR v,LIS_INT index0,LIS_INT index1,LIS_SCALAR value)
{
    if(set_global_entry(v,index0,value)) return 1;
    if(set_global_entry(v,index1,value)) return 1;
    return 0;
}

static int subtract_global_pair(
    LIS_VECTOR v,LIS_INT index0,LIS_INT index1,LIS_SCALAR value)
{
    LIS_INT is,ie,i,err;

    err=lis_vector_get_range(v,&is,&ie);
    if(err) return 1;

    for(i=is;i<ie;i++)
    {
        if(i==index0 || i==index1)
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
        if(i==Q0_INDEX0 || i==Q0_INDEX1)
        {
#ifdef _COMPLEX
            (*q0)->value[i-is]=(LIS_SCALAR)(1.0+1.0*I);
#else
            (*q0)->value[i-is]=(LIS_SCALAR)1.0;
#endif
        }
        if(i==Q1_INDEX0 || i==Q1_INDEX1)
            (*q1)->value[i-is]=(LIS_SCALAR)1.0;
    }
    return 0;
}

#ifdef _COMPLEX
static int check_complex_dot(LIS_VECTOR q0,LIS_VECTOR b)
{
    LIS_SCALAR coefficient;
    LIS_INT err;

    err=lis_vector_set_all((LIS_SCALAR)0.0,b); if(err) return 1;
    if(set_global_pair(
        b,Q0_INDEX0,Q0_INDEX1,(LIS_SCALAR)2.0)) return 1;

    err=lis_vector_dot(q0,b,&coefficient);
    if(err) return 1;

    if(fabs(creal(coefficient)-4.0)>(LIS_REAL)1.0e-11)
        return 1;
    if(fabs(cimag(coefficient)+4.0)>(LIS_REAL)1.0e-11)
        return 1;

    return 0;
}
#endif

static int check_lift(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR lift,
    LIS_SCALAR rhs_value,
    LIS_SCALAR expected_value)
{
    const LIS_REAL sqrt_two=(LIS_REAL)1.4142135623730950488;
    LIS_NULLSPACE_LIFT_DIAGNOSTICS diag;
    LIS_REAL nrm;
    LIS_REAL expected_rhs_norm;
    LIS_REAL expected_lift_norm;
    LIS_INT err;

    err=lis_vector_set_all((LIS_SCALAR)0.0,b); if(err) return 1;
    if(set_global_pair(b,Q0_INDEX0,Q0_INDEX1,rhs_value)) return 1;

    err=lis_solver_build_nullspace_lift(
        solver,b,(LIS_REAL)2.0,(LIS_REAL)1.0e-12,lift,&diag);
    if(err) return 1;

    expected_rhs_norm=sqrt_two*(LIS_REAL)fabs(rhs_value);
    expected_lift_norm=sqrt_two*(LIS_REAL)fabs(expected_value);

    if(diag.nullspace_dim!=2 || diag.active_count!=1)
        return 1;
    if(!real_close(diag.rhs_norm,expected_rhs_norm,(LIS_REAL)1.0e-11))
        return 1;
    if(!real_close(diag.max_relative_component,(LIS_REAL)1.0,
                   (LIS_REAL)1.0e-11))
        return 1;
    if(!real_close(diag.lift_norm,expected_lift_norm,(LIS_REAL)1.0e-11))
        return 1;

    if(subtract_global_pair(
        lift,Q0_INDEX0,Q0_INDEX1,expected_value)) return 1;
    err=lis_vector_nrm2(lift,&nrm);
    if(err || nrm>(LIS_REAL)1.0e-11)
        return 1;

    return 0;
}

static int check_zero_rhs(
    LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR lift)
{
    LIS_NULLSPACE_LIFT_DIAGNOSTICS diag;
    LIS_REAL nrm;
    LIS_INT err;

    err=lis_vector_set_all((LIS_SCALAR)0.0,b); if(err) return 1;
    err=lis_vector_set_all((LIS_SCALAR)3.0,lift); if(err) return 1;

    err=lis_solver_build_nullspace_lift(
        solver,b,(LIS_REAL)2.0,(LIS_REAL)0.5,lift,&diag);
    if(err) return 1;

    if(diag.nullspace_dim!=2 || diag.active_count!=0)
        return 1;
    if(diag.rhs_norm!=(LIS_REAL)0.0 ||
       diag.max_relative_component!=(LIS_REAL)0.0 ||
       diag.lift_norm!=(LIS_REAL)0.0)
        return 1;

    err=lis_vector_nrm2(lift,&nrm);
    if(err || nrm!=(LIS_REAL)0.0)
        return 1;

    return 0;
}

static int check_threshold(
    LIS_SOLVER solver,LIS_VECTOR b,LIS_VECTOR lift)
{
    const LIS_REAL one_over_sqrt_two=
        (LIS_REAL)0.7071067811865475244;
    LIS_NULLSPACE_LIFT_DIAGNOSTICS diag;
    LIS_REAL boundary,nrm;
    LIS_INT err;

    err=lis_vector_set_all((LIS_SCALAR)0.0,b); if(err) return 1;
    if(set_global_entry(b,Q0_INDEX0,(LIS_SCALAR)2.0)) return 1;

    err=lis_solver_build_nullspace_lift(
        solver,b,(LIS_REAL)2.0,(LIS_REAL)2.0,lift,&diag);
    if(err) return 1;

    if(diag.nullspace_dim!=2 || diag.active_count!=0)
        return 1;
    if(!real_close(diag.rhs_norm,(LIS_REAL)2.0,(LIS_REAL)1.0e-11))
        return 1;
    if(!real_close(diag.max_relative_component,one_over_sqrt_two,
                   (LIS_REAL)1.0e-11))
        return 1;
    if(diag.lift_norm!=(LIS_REAL)0.0)
        return 1;

    boundary=diag.max_relative_component;

    err=lis_solver_build_nullspace_lift(
        solver,b,(LIS_REAL)2.0,boundary,lift,&diag);
    if(err) return 1;

    if(diag.active_count!=0 || diag.lift_norm!=(LIS_REAL)0.0)
        return 1;
    if(!real_close(diag.max_relative_component,boundary,
                   (LIS_REAL)1.0e-14))
        return 1;

    err=lis_solver_build_nullspace_lift(
        solver,b,(LIS_REAL)2.0,
        boundary*(LIS_REAL)(1.0-1.0e-8),lift,&diag);
    if(err) return 1;

    if(diag.active_count!=1)
        return 1;
    if(!real_close(diag.max_relative_component,boundary,
                   (LIS_REAL)1.0e-14))
        return 1;
    if(!real_close(diag.lift_norm,one_over_sqrt_two,
                   (LIS_REAL)1.0e-11))
        return 1;

    if(subtract_global_pair(
        lift,Q0_INDEX0,Q0_INDEX1,(LIS_SCALAR)0.5)) return 1;
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

#ifdef _COMPLEX
    if(check_complex_dot(q0,b)) goto done;
#endif

    basis[0]=q0; basis[1]=q1;
    err=lis_solver_create(&solver); if(err) goto done;
    err=lis_solver_set_nullspace(solver,2,basis); if(err) goto done;

    if(check_lift(
        solver,b,lift,
        (LIS_SCALAR)2.0,
        (LIS_SCALAR)1.0))
        goto done;

    if(check_zero_rhs(solver,b,lift))
        goto done;

    if(check_threshold(solver,b,lift))
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
