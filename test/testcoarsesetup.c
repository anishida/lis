#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include <stdio.h>

#include "lis.h"
#include "lis_solver.h"

#define TEST_N 8
#define EPS_DIAG 1.0e-14


static LIS_REAL
test_scalar_real(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return creal(value);
#else
    return value;
#endif
}


static LIS_INT
create_diagonal_matrix(
    LIS_SCALAR first,
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;

    LIS_INT is,ie,i,err;


    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if(err)
        return err;


    err = lis_matrix_set_size(
        A,
        0,
        TEST_N);

    if(err)
        goto fail;


    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);

    if(err)
        goto fail;


    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if(err)
        goto fail;


    for(i=is;i<ie;i++)
    {
        LIS_SCALAR value;

        if(i==0)
            value = first;
        else if(i==1)
            value = (LIS_SCALAR)4.0;
        else
            value = (LIS_SCALAR)(i+3);


        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            value,
            A);

        if(err)
            goto fail;
    }


    err = lis_matrix_assemble(A);

    if(err)
        goto fail;


    *Aout = A;

    return LIS_SUCCESS;


fail:

    if(A)
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
create_basis(
    LIS_MATRIX A,
    LIS_VECTOR *z0,
    LIS_VECTOR *z1)
{
    LIS_INT err;


    *z0 = NULL;
    *z1 = NULL;


    err = lis_vector_duplicate(
        A,
        z0);

    if(err)
        return err;


    err = lis_vector_duplicate(
        A,
        z1);

    if(err)
    {
        lis_vector_destroy(*z0);
        *z0 = NULL;

        return err;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        *z0);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        *z1);


    err = lis_vector_set_value(
        LIS_INS_VALUE,
        0,
        (LIS_SCALAR)1.0,
        *z0);

    if(err)
        return err;


    err = lis_vector_set_value(
        LIS_INS_VALUE,
        1,
        (LIS_SCALAR)1.0,
        *z1);

    return err;
}


static int
near_scalar(
    LIS_SCALAR actual,
    LIS_SCALAR expected,
    LIS_REAL rel)
{
    LIS_REAL scale;

    scale = fabs(expected);

    /*
     * Preserve relative sensitivity for small but nonzero
     * coarse eigenvalues. Only exact zero uses an absolute
     * reference scale.
     */
    if(scale==0.0)
        scale = 1.0;


    return
        fabs(actual-expected)
        <=
        rel*scale;
}


int
main(
    int argc,
    char **argv)
{
    LIS_MATRIX A = NULL;
    LIS_MATRIX Asing = NULL;

    LIS_VECTOR z0 = NULL;
    LIS_VECTOR z1 = NULL;

    LIS_VECTOR basis[2];

    LIS_SOLVER solver = NULL;

    LIS_SCALAR rhs[2];
    LIS_SCALAR coeff[2];

    LIS_NEAR_NULLSPACE_COARSE coarse;

    LIS_INT err;

    int failed = 0;


    err = lis_initialize(
        &argc,
        &argv);

    if(err)
        return 1;


    err = create_diagonal_matrix(
        (LIS_SCALAR)EPS_DIAG,
        &A);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = create_diagonal_matrix(
        (LIS_SCALAR)0.0,
        &Asing);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = create_basis(
        A,
        &z0,
        &z1);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_solver_create(
        &solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    basis[0] = z0;
    basis[1] = z1;


    err = lis_solver_set_near_nullspace(
        solver,
        2,
        basis);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    /*
     * Caller vectors are no longer required after set().
     */
    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z0);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z1);


    err = lis_solver_near_nullspace_coarse_setup(
        solver,
        A,
        LIS_SCALE_NONE);


    printf(
        "SETUP rc=%d\n",
        (int)err);


    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    coarse =
        solver->near_nullspace_coarse;


    if(
        coarse==NULL
        ||
        !coarse->ready
        ||
        coarse->dim!=2)
    {
        fprintf(
            stderr,
            "FAIL coarse state not ready\n");

        failed = 1;
        goto cleanup;
    }


    printf(
        "COARSE dim=%d ready=%d\n",
        (int)coarse->dim,
        (int)coarse->ready);


    printf(
        "E00=%.16e E01=%.16e E10=%.16e E11=%.16e\n",
        (double)test_scalar_real(coarse->E[0]),
        (double)test_scalar_real(coarse->E[1]),
        (double)test_scalar_real(coarse->E[2]),
        (double)test_scalar_real(coarse->E[3]));


    if(
        !near_scalar(
            coarse->E[0],
            (LIS_SCALAR)EPS_DIAG,
            1.0e-12)
        ||
        !near_scalar(
            coarse->E[1],
            (LIS_SCALAR)0.0,
            1.0e-12)
        ||
        !near_scalar(
            coarse->E[2],
            (LIS_SCALAR)0.0,
            1.0e-12)
        ||
        !near_scalar(
            coarse->E[3],
            (LIS_SCALAR)4.0,
            1.0e-12))
    {
        fprintf(
            stderr,
            "FAIL unexpected coarse matrix\n");

        failed = 1;
    }


    rhs[0] =
        (LIS_SCALAR)(2.0*EPS_DIAG);

    rhs[1] =
        (LIS_SCALAR)12.0;


    err = lis_solver_near_nullspace_coarse_solve(
        solver,
        rhs,
        coeff);


    printf(
        "COARSE_SOLVE rc=%d c0=%.16e c1=%.16e\n",
        (int)err,
        (double)test_scalar_real(coeff[0]),
        (double)test_scalar_real(coeff[1]));


    if(
        err
        ||
        !near_scalar(
            coeff[0],
            (LIS_SCALAR)2.0,
            1.0e-11)
        ||
        !near_scalar(
            coeff[1],
            (LIS_SCALAR)3.0,
            1.0e-11))
    {
        fprintf(
            stderr,
            "FAIL reusable coarse solve\n");

        failed = 1;
    }


    /*
     * Rebuild from a singular effective operator.
     * This must fail cleanly and invalidate old state.
     */
    err = lis_solver_near_nullspace_coarse_setup(
        solver,
        Asing,
        LIS_SCALE_NONE);


    printf(
        "SINGULAR rc=%d expected=%d coarse=%p\n",
        (int)err,
        (int)LIS_BREAKDOWN,
        (void *)solver->near_nullspace_coarse);


    if(
        err!=LIS_BREAKDOWN
        ||
        solver->near_nullspace_coarse!=NULL)
    {
        fprintf(
            stderr,
            "FAIL singular coarse operator handling\n");

        failed = 1;
    }


    /*
     * Recovery after failed setup.
     */
    err = lis_solver_near_nullspace_coarse_setup(
        solver,
        A,
        LIS_SCALE_NONE);


    printf(
        "RECOVERY rc=%d ready=%d\n",
        (int)err,
        solver->near_nullspace_coarse
        ?
        (int)solver->near_nullspace_coarse->ready
        :
        -1);


    if(
        err
        ||
        solver->near_nullspace_coarse==NULL
        ||
        !solver->near_nullspace_coarse->ready)
    {
        fprintf(
            stderr,
            "FAIL recovery after coarse setup breakdown\n");

        failed = 1;
    }


    /*
     * Clearing the user basis must also clear
     * the derived coarse state.
     */
    err = lis_solver_clear_near_nullspace(
        solver);


    printf(
        "CLEAR rc=%d dim=%d coarse=%p\n",
        (int)err,
        (int)solver->near_nullspace_dim,
        (void *)solver->near_nullspace_coarse);


    if(
        err
        ||
        solver->near_nullspace_dim!=0
        ||
        solver->near_nullspace_coarse!=NULL)
    {
        fprintf(
            stderr,
            "FAIL clear did not invalidate coarse state\n");

        failed = 1;
    }


cleanup:

    if(solver)
        lis_solver_destroy(solver);

    if(z1)
        lis_vector_destroy(z1);

    if(z0)
        lis_vector_destroy(z0);

    if(Asing)
        lis_matrix_destroy(Asing);

    if(A)
        lis_matrix_destroy(A);


    printf(
        "LIS_COARSE_SETUP_STAGE7B %s\n",
        failed ? "FAILED" : "PASSED");


    lis_finalize();


    return failed ? 1 : 0;
}
