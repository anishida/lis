#include "lis.h"
#include <stdio.h>

static LIS_INT
create_matrix(LIS_MATRIX *A)
{
    LIS_INT err;

    err = lis_matrix_create(LIS_COMM_WORLD, A);
    if (err) return err;

    err = lis_matrix_set_size(*A, 0, 2);
    if (err) return err;

    err = lis_matrix_set_value(LIS_INS_VALUE, 0, 0, 2.0, *A);
    if (err) return err;

    err = lis_matrix_set_value(LIS_INS_VALUE, 1, 1, 3.0, *A);
    if (err) return err;

    return lis_matrix_assemble(*A);
}

int
main(int argc, char **argv)
{
    LIS_INT err, status;
    LIS_MATRIX A = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_SOLVER solver = NULL;
    int failed = 0;

    err = lis_initialize(&argc, &argv);
    if (err) return 1;

    err = create_matrix(&A);
    if (err) goto cleanup;

    err = lis_vector_duplicate(A, &b);
    if (err) goto cleanup;

    err = lis_vector_duplicate(A, &x);
    if (err) goto cleanup;

    lis_vector_set_all(1.0, b);
    lis_vector_set_all(0.0, x);

    err = lis_solver_create(&solver);
    if (err) goto cleanup;

    err = lis_solver_set_option("-i cg -p none -maxiter 20", solver);
    if (err) goto cleanup;

    /*
     * First solve creates a valid solver-owned residual history.
     */
    err = lis_solve(A, b, x, solver);
    if (err)
    {
        printf("first lis_solve returned %d\n", (int)err);
        failed = 1;
        goto cleanup;
    }

    if (solver->rhistory == NULL)
    {
        printf("first solve did not create rhistory\n");
        failed = 1;
        goto cleanup;
    }


    /*
     * Force matrix conversion setup to fail during a repeated solve.
     * The failed solve must not leave stale residual-history state.
     */
    solver->options[LIS_OPTIONS_STORAGE] = 999;

    err = lis_solve(A, b, x, solver);


    if (err == LIS_SUCCESS)
    {
        printf("expected second solve to fail\n");
        failed = 1;
    }

    if (solver->rhistory != NULL)
    {
        printf("solver->rhistory is non-NULL after failed setup\n");

        /*
         * Keep cleanup safe if this regression ever reappears.
         */
        solver->rhistory = NULL;
        failed = 1;
    }

    err = lis_solver_get_status(solver, &status);
    if (err || status != LIS_ERR_ILL_ARG)
    {
        printf("status after failed solve = %d, expected %d\n",
               (int)status, (int)LIS_ERR_ILL_ARG);
        failed = 1;
    }

    /*
     * The same solver must remain reusable after the setup error.
     */
    solver->options[LIS_OPTIONS_STORAGE] = 0;
    lis_vector_set_all(0.0, x);

    err = lis_solve(A, b, x, solver);
    if (err)
    {
        printf("third lis_solve returned %d\n", (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_get_status(solver, &status);
    if (err || status != LIS_SUCCESS)
    {
        printf("status after recovery solve = %d, expected LIS_SUCCESS\n",
               (int)status);
        failed = 1;
    }

    if (solver->rhistory == NULL)
    {
        printf("recovery solve did not create rhistory\n");
        failed = 1;
    }

cleanup:
    if (solver) lis_solver_destroy(solver);
    if (x) lis_vector_destroy(x);
    if (b) lis_vector_destroy(b);
    if (A) lis_matrix_destroy(A);

    lis_finalize();

    if (failed)
    {
        printf("LIS_RHISTORY_ERROR FAILED\n");
        return 1;
    }

    printf("LIS_RHISTORY_ERROR PASSED\n");
    return 0;
}
