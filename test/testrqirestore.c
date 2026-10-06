#include "lis.h"
#include <stdio.h>

static LIS_INT
create_diagonal_matrix(LIS_MATRIX *A, LIS_SCALAR a0, LIS_SCALAR a1)
{
    LIS_INT err;

    err = lis_matrix_create(LIS_COMM_WORLD, A);
    if (err) return err;

    err = lis_matrix_set_size(*A, 0, 2);
    if (err) return err;

    err = lis_matrix_set_value(LIS_INS_VALUE, 0, 0, a0, *A);
    if (err) return err;

    err = lis_matrix_set_value(LIS_INS_VALUE, 1, 1, a1, *A);
    if (err) return err;

    return lis_matrix_assemble(*A);
}

static LIS_INT
check_diagonal_unchanged(LIS_MATRIX A, LIS_VECTOR before)
{
    LIS_INT err;
    LIS_VECTOR after;
    LIS_REAL nrm2;

    err = lis_vector_duplicate(A, &after);
    if (err) return err;

    err = lis_matrix_get_diagonal(A, after);
    if (err)
    {
        lis_vector_destroy(after);
        return err;
    }

    lis_vector_axpy(-1.0, before, after);
    lis_vector_nrm2(after, &nrm2);
    lis_vector_destroy(after);

    if (nrm2 != 0.0)
    {
        printf("matrix diagonal was not restored: ||delta||_2 = %e\n",
               (double)nrm2);
        return LIS_FAILS;
    }

    return LIS_SUCCESS;
}

static LIS_INT
test_rqi_restore(void)
{
    LIS_INT err, check, status;
    LIS_MATRIX A;
    LIS_VECTOR x, diag;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    A = NULL;
    x = NULL;
    diag = NULL;
    esolver = NULL;

    err = create_diagonal_matrix(&A, 2.0, 3.0);
    if (err) goto cleanup;

    err = lis_vector_duplicate(A, &x);
    if (err) goto cleanup;

    err = lis_vector_duplicate(A, &diag);
    if (err) goto cleanup;

    lis_vector_set_all(1.0, x);

    err = lis_matrix_get_diagonal(A, diag);
    if (err) goto cleanup;

    err = lis_esolver_create(&esolver);
    if (err) goto cleanup;

    /*
     * -maxiter_noimp -1 is supplied on the command line.
     * The inner RQI linear solver accepts the option via
     * lis_solver_set_optionC() and rejects it in lis_solve_kernel(),
     * after the temporary matrix shift has already been applied.
     */
    lis_esolver_set_option("-e rqi -emaxiter 2", esolver);

    err = lis_esolve(A, x, &evalue, esolver);

    if (err != LIS_SUCCESS)
    {
        printf("lis_esolve returned API error %d\n", (int)err);
        err = LIS_FAILS;
        goto cleanup;
    }

    lis_esolver_get_status(esolver, &status);

    if (status != LIS_ERR_ILL_ARG)
    {
        printf("RQI status %d, expected LIS_ERR_ILL_ARG (%d)\n",
               (int)status, (int)LIS_ERR_ILL_ARG);
        err = LIS_FAILS;
        goto cleanup;
    }

    check = check_diagonal_unchanged(A, diag);
    if (check)
    {
        printf("RQI did not restore A after solve-kernel error\n");
        err = LIS_FAILS;
        goto cleanup;
    }

    err = LIS_SUCCESS;

cleanup:
    if (esolver) lis_esolver_destroy(esolver);
    if (diag) lis_vector_destroy(diag);
    if (x) lis_vector_destroy(x);
    if (A) lis_matrix_destroy(A);

    return err;
}

static LIS_INT
test_grqi_restore(void)
{
    LIS_INT err, check, status;
    LIS_MATRIX A, B;
    LIS_VECTOR x, diag;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    A = NULL;
    B = NULL;
    x = NULL;
    diag = NULL;
    esolver = NULL;

    err = create_diagonal_matrix(&A, 2.0, 3.0);
    if (err) goto cleanup;

    err = create_diagonal_matrix(&B, 1.0, 1.0);
    if (err) goto cleanup;

    err = lis_vector_duplicate(A, &x);
    if (err) goto cleanup;

    err = lis_vector_duplicate(A, &diag);
    if (err) goto cleanup;

    lis_vector_set_all(1.0, x);

    err = lis_matrix_get_diagonal(A, diag);
    if (err) goto cleanup;

    err = lis_esolver_create(&esolver);
    if (err) goto cleanup;

    lis_esolver_set_option("-e grqi -emaxiter 2", esolver);

    err = lis_gesolve(A, B, x, &evalue, esolver);

    if (err != LIS_SUCCESS)
    {
        printf("lis_gesolve returned API error %d\n", (int)err);
        err = LIS_FAILS;
        goto cleanup;
    }

    lis_esolver_get_status(esolver, &status);

    if (status != LIS_ERR_ILL_ARG)
    {
        printf("GRQI status %d, expected LIS_ERR_ILL_ARG (%d)\n",
               (int)status, (int)LIS_ERR_ILL_ARG);
        err = LIS_FAILS;
        goto cleanup;
    }

    check = check_diagonal_unchanged(A, diag);
    if (check)
    {
        printf("GRQI did not restore A after solve-kernel error\n");
        err = LIS_FAILS;
        goto cleanup;
    }

    err = LIS_SUCCESS;

cleanup:
    if (esolver) lis_esolver_destroy(esolver);
    if (diag) lis_vector_destroy(diag);
    if (x) lis_vector_destroy(x);
    if (B) lis_matrix_destroy(B);
    if (A) lis_matrix_destroy(A);

    return err;
}

int
main(int argc, char **argv)
{
    LIS_INT err;
    int test_argc;
    char arg0[] = "testrqirestore";
    char arg1[] = "-maxiter_noimp";
    char arg2[] = "-1";
    char *test_argv[] = {arg0, arg1, arg2, NULL};
    char **test_argvp = test_argv;

    (void)argc;
    (void)argv;

    test_argc = 3;
    err = lis_initialize(&test_argc, &test_argvp);
    if (err) return 1;

    err = test_rqi_restore();
    if (err)
    {
        printf("LIS_RQI_RESTORE FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_grqi_restore();
    if (err)
    {
        printf("LIS_GRQI_RESTORE FAILED\n");
        lis_finalize();
        return 1;
    }

    printf("LIS_RQI_RESTORE PASSED\n");

    lis_finalize();
    return 0;
}
