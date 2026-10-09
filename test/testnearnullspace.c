#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lis.h"

#include <math.h>
#include <stdio.h>

#define TEST_N 8


static int
check_owned_copy(
    LIS_SOLVER solver,
    LIS_INT expected_dim,
    LIS_REAL expected_norm)
{
    LIS_REAL nrm = 0.0;

    if( solver->near_nullspace_dim!=expected_dim )
    {
        fprintf(
            stderr,
            "FAIL dimension expected=%d actual=%d\n",
            (int)expected_dim,
            (int)solver->near_nullspace_dim);

        return 1;
    }


    if( expected_dim==0 )
    {
        if( solver->near_nullspace!=NULL )
        {
            fprintf(
                stderr,
                "FAIL cleared basis pointer is not NULL\n");

            return 1;
        }

        return 0;
    }


    if(
        solver->near_nullspace==NULL
        ||
        solver->near_nullspace[0]==NULL)
    {
        fprintf(
            stderr,
            "FAIL owned basis is NULL\n");

        return 1;
    }


    if(
        lis_vector_nrm2(
            solver->near_nullspace[0],
            &nrm))
    {
        fprintf(
            stderr,
            "FAIL cannot evaluate owned basis norm\n");

        return 1;
    }


    if(
        fabs(nrm-expected_norm)
        >
        1.0e-12*(1.0+expected_norm))
    {
        fprintf(
            stderr,
            "FAIL norm expected=%.16e actual=%.16e\n",
            (double)expected_norm,
            (double)nrm);

        return 1;
    }


    return 0;
}


int
main(
    int argc,
    char **argv)
{
    LIS_SOLVER solver = NULL;

    LIS_VECTOR v0 = NULL;
    LIS_VECTOR v1 = NULL;
    LIS_VECTOR vbad = NULL;

    LIS_VECTOR basis2[2];
    LIS_VECTOR basis1[1];
    LIS_VECTOR null_basis[1];

    LIS_INT err;

    LIS_REAL expected;

    int failed = 0;


    err = lis_initialize(
        &argc,
        &argv);

    if(err)
        return 1;


    err = lis_solver_create(
        &solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_create(
        LIS_COMM_WORLD,
        &v0);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_set_size(
        v0,
        0,
        TEST_N);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(
        v0,
        &v1);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_create(
        LIS_COMM_WORLD,
        &vbad);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_set_size(
        vbad,
        0,
        TEST_N+1);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        v0);

    lis_vector_set_all(
        (LIS_SCALAR)3.0,
        v1);


    basis2[0] = v0;
    basis2[1] = v1;


    err = lis_solver_set_near_nullspace(
        solver,
        2,
        basis2);


    printf(
        "SET2 rc=%d dim=%d\n",
        (int)err,
        (int)solver->near_nullspace_dim);


    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    if(
        solver->near_nullspace[0]==v0
        ||
        solver->near_nullspace[1]==v1)
    {
        fprintf(
            stderr,
            "FAIL solver did not make owned copies\n");

        failed = 1;
    }


    /*
     * Mutating caller-owned vectors must not modify
     * the solver-owned copies.
     */
    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        v0);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        v1);


    expected =
        sqrt((LIS_REAL)TEST_N);


    failed |=
        check_owned_copy(
            solver,
            2,
            expected);


    /*
     * Replacement must be transactional and release
     * the previously owned basis only after the new
     * basis has been copied successfully.
     */
    lis_vector_set_all(
        (LIS_SCALAR)2.0,
        v0);

    basis1[0] = v0;


    err = lis_solver_set_near_nullspace(
        solver,
        1,
        basis1);


    printf(
        "REPLACE1 rc=%d dim=%d\n",
        (int)err,
        (int)solver->near_nullspace_dim);


    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        v0);


    expected =
        2.0*sqrt((LIS_REAL)TEST_N);


    failed |=
        check_owned_copy(
            solver,
            1,
            expected);


    /*
     * Invalid replacement must leave the existing
     * valid basis unchanged.
     */
    null_basis[0] = NULL;

    err = lis_solver_set_near_nullspace(
        solver,
        1,
        null_basis);


    printf(
        "NULL_VECTOR rc=%d expected=%d dim=%d\n",
        (int)err,
        (int)LIS_ERR_ILL_ARG,
        (int)solver->near_nullspace_dim);


    if(err!=LIS_ERR_ILL_ARG)
        failed = 1;


    failed |=
        check_owned_copy(
            solver,
            1,
            expected);


    /*
     * Different layouts are rejected.
     */
    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        v0);

    basis2[0] = v0;
    basis2[1] = vbad;


    err = lis_solver_set_near_nullspace(
        solver,
        2,
        basis2);


    printf(
        "BAD_LAYOUT rc=%d expected=%d dim=%d\n",
        (int)err,
        (int)LIS_ERR_ILL_ARG,
        (int)solver->near_nullspace_dim);


    if(err!=LIS_ERR_ILL_ARG)
        failed = 1;


    failed |=
        check_owned_copy(
            solver,
            1,
            expected);


    err = lis_solver_clear_near_nullspace(
        solver);


    printf(
        "CLEAR rc=%d dim=%d ptr=%p\n",
        (int)err,
        (int)solver->near_nullspace_dim,
        (void *)solver->near_nullspace);


    if(err)
        failed = 1;


    failed |=
        check_owned_copy(
            solver,
            0,
            0.0);


    /*
     * nvec == 0 is also defined as clear.
     */
    err = lis_solver_set_near_nullspace(
        solver,
        0,
        NULL);


    printf(
        "SET0 rc=%d dim=%d\n",
        (int)err,
        (int)solver->near_nullspace_dim);


    if(err)
        failed = 1;


    err = lis_solver_set_near_nullspace(
        solver,
        -1,
        NULL);


    printf(
        "NEGATIVE_DIM rc=%d expected=%d\n",
        (int)err,
        (int)LIS_ERR_ILL_ARG);


    if(err!=LIS_ERR_ILL_ARG)
        failed = 1;


    /*
     * Leave one owned basis installed so solver_destroy()
     * exercises the ownership cleanup path.
     */
    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        v1);

    basis1[0] = v1;


    err = lis_solver_set_near_nullspace(
        solver,
        1,
        basis1);


    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    printf(
        "DESTROY_WITH_BASIS dim=%d\n",
        (int)solver->near_nullspace_dim);


cleanup:

    if(solver)
    {
        lis_solver_destroy(
            solver);

        solver = NULL;
    }


    if(vbad)
        lis_vector_destroy(vbad);

    if(v1)
        lis_vector_destroy(v1);

    if(v0)
        lis_vector_destroy(v0);


    printf(
        "LIS_NEAR_NULLSPACE_STAGE7A %s\n",
        failed ? "FAILED" : "PASSED");


    lis_finalize();


    return failed ? 1 : 0;
}
