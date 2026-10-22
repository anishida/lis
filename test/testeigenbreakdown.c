#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include "lis.h"
#include <stdio.h>


static LIS_INT
create_diagonal_matrix(
    LIS_MATRIX *A,
    LIS_SCALAR a0,
    LIS_SCALAR a1)
{
    LIS_INT err;

    *A = NULL;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        A);
    if( err )
        return err;

    err = lis_matrix_set_size(
        *A,
        0,
        2);
    if( err )
        return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,
        0,
        0,
        a0,
        *A);
    if( err )
        return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,
        1,
        1,
        a1,
        *A);
    if( err )
        return err;

    return lis_matrix_assemble(*A);
}


static LIS_INT
check_standard_breakdown(
    const char *name,
    const char *options,
    LIS_MATRIX A,
    LIS_SCALAR initial_value)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    x = NULL;
    esolver = NULL;
    evalue = 0.0;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err )
        goto cleanup;

    err = lis_vector_set_all(
        initial_value,
        x);
    if( err )
        goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    err = lis_esolver_set_option(
        (char *)options,
        esolver);
    if( err )
        goto cleanup;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS &&
        api!=LIS_BREAKDOWN )
    {
        printf(
            "%s: unexpected API return %d\n",
            name,
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,
        &status);
    if( err )
        goto cleanup;

    if( status!=LIS_BREAKDOWN )
    {
        printf(
            "%s: status %d, expected LIS_BREAKDOWN (%d)\n",
            name,
            (int)status,
            (int)LIS_BREAKDOWN);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: %s -> LIS_BREAKDOWN\n",
        name);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
check_generalized_breakdown(
    const char *name,
    const char *options,
    LIS_MATRIX A,
    LIS_MATRIX B)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    x = NULL;
    esolver = NULL;
    evalue = 0.0;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err )
        goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x);
    if( err )
        goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    err = lis_esolver_set_option(
        (char *)options,
        esolver);
    if( err )
        goto cleanup;

    api = lis_gesolve(
        A,
        B,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS &&
        api!=LIS_BREAKDOWN )
    {
        printf(
            "%s: unexpected API return %d\n",
            name,
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,
        &status);
    if( err )
        goto cleanup;

    if( status!=LIS_BREAKDOWN )
    {
        printf(
            "%s: status %d, expected LIS_BREAKDOWN (%d)\n",
            name,
            (int)status,
            (int)LIS_BREAKDOWN);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: %s -> LIS_BREAKDOWN\n",
        name);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
check_arnoldi_happy_breakdown(
    LIS_MATRIX A)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    x = NULL;
    esolver = NULL;
    evalue = 0.0;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err )
        goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x);
    if( err )
        goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    /*
     * For A = I, the first Arnoldi vector is already an
     * invariant subspace.  Orthogonalization therefore
     * produces an exactly zero next Hessenberg norm.
     *
     * This is a finite happy breakdown and must not be
     * reported as LIS_BREAKDOWN.
     */
    err = lis_esolver_set_option(
        "-e ai -ss 2 -emaxiter 4 -etol 1.0e-12",
        esolver);
    if( err )
        goto cleanup;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Arnoldi happy breakdown: API return %d\n",
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,
        &status);
    if( err )
        goto cleanup;

    if( status==LIS_BREAKDOWN )
    {
        printf(
            "Arnoldi happy breakdown was reported as LIS_BREAKDOWN\n");
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Arnoldi finite happy breakdown -> status %d\n",
        (int)status);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
test_power_zero_operator(void)
{
    LIS_INT err;
    LIS_MATRIX A;

    A = NULL;

    err = create_diagonal_matrix(
        &A,
        0.0,
        0.0);
    if( err )
        goto cleanup;

    err = check_standard_breakdown(
        "Power zero operator",
        "-e pi -emaxiter 4 -etol 1.0e-12",
        A,
        1.0);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_power_degenerate_metric(void)
{
    LIS_INT err;
    LIS_MATRIX A;
    LIS_MATRIX B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,
        1.0,
        2.0);
    if( err )
        goto cleanup;

    err = create_diagonal_matrix(
        &B,
        0.0,
        0.0);
    if( err )
        goto cleanup;

    err = check_generalized_breakdown(
        "Generalized Power degenerate metric",
        "-e gpi -emaxiter 4 -etol 1.0e-12",
        A,
        B);

cleanup:
    if( B )
        lis_matrix_destroy(B);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_inverse_degenerate_metric(void)
{
    LIS_INT err;
    LIS_MATRIX A;
    LIS_MATRIX B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,
        2.0,
        3.0);
    if( err )
        goto cleanup;

    err = create_diagonal_matrix(
        &B,
        0.0,
        0.0);
    if( err )
        goto cleanup;

    err = check_generalized_breakdown(
        "Generalized Inverse degenerate metric",
        "-e gii -emaxiter 4 -etol 1.0e-12",
        A,
        B);

cleanup:
    if( B )
        lis_matrix_destroy(B);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_grqi_degenerate_metric(void)
{
    LIS_INT err;
    LIS_MATRIX A;
    LIS_MATRIX B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,
        2.0,
        3.0);
    if( err )
        goto cleanup;

    err = create_diagonal_matrix(
        &B,
        0.0,
        0.0);
    if( err )
        goto cleanup;

    err = check_generalized_breakdown(
        "Generalized RQI degenerate metric",
        "-e grqi -emaxiter 4 -etol 1.0e-12",
        A,
        B);

cleanup:
    if( B )
        lis_matrix_destroy(B);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_cg_zero_initial_vector(void)
{
    LIS_INT err;
    LIS_MATRIX A;

    A = NULL;

    err = create_diagonal_matrix(
        &A,
        1.0,
        2.0);
    if( err )
        goto cleanup;

    err = check_standard_breakdown(
        "Conjugate Gradient zero initial vector",
        "-e cg -initx_ones false -emaxiter 4 -etol 1.0e-12",
        A,
        0.0);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_cg_degenerate_metric(void)
{
    LIS_INT err;
    LIS_MATRIX A;
    LIS_MATRIX B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,
        1.0,
        2.0);
    if( err )
        goto cleanup;

    err = create_diagonal_matrix(
        &B,
        0.0,
        0.0);
    if( err )
        goto cleanup;

    err = check_generalized_breakdown(
        "Generalized Conjugate Gradient degenerate metric",
        "-e gcg -emaxiter 4 -etol 1.0e-12",
        A,
        B);

cleanup:
    if( B )
        lis_matrix_destroy(B);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_cr_zero_rayleigh(void)
{
    LIS_INT err;
    LIS_MATRIX A;

    A = NULL;

    /*
     * With the default all-ones initial vector and
     * A = diag(1,-1), the initial Rayleigh quotient is zero.
     */
    err = create_diagonal_matrix(
        &A,
        1.0,
        -1.0);
    if( err )
        goto cleanup;

    err = check_standard_breakdown(
        "Conjugate Residual zero Rayleigh denominator",
        "-e cr -emaxiter 4 -etol 1.0e-12",
        A,
        1.0);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_cr_degenerate_metric(void)
{
    LIS_INT err;
    LIS_MATRIX A;
    LIS_MATRIX B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,
        1.0,
        2.0);
    if( err )
        goto cleanup;

    err = create_diagonal_matrix(
        &B,
        0.0,
        0.0);
    if( err )
        goto cleanup;

    err = check_generalized_breakdown(
        "Generalized Conjugate Residual degenerate metric",
        "-e gcr -emaxiter 4 -etol 1.0e-12",
        A,
        B);

cleanup:
    if( B )
        lis_matrix_destroy(B);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_arnoldi_happy_breakdown(void)
{
    LIS_INT err;
    LIS_MATRIX A;

    A = NULL;

    err = create_diagonal_matrix(
        &A,
        1.0,
        1.0);
    if( err )
        goto cleanup;

    err = check_arnoldi_happy_breakdown(
        A);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


int
main(
    int argc,
    char **argv)
{
    LIS_INT err;

    err = lis_initialize(
        &argc,
        &argv);
    if( err )
        return 1;

    err = test_power_zero_operator();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN POWER FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_power_degenerate_metric();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN GPOWER FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_inverse_degenerate_metric();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN GINVERSE FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_grqi_degenerate_metric();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN GRQI FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_cg_zero_initial_vector();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN CG FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_cg_degenerate_metric();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN GCG FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_cr_zero_rayleigh();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN CR FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_cr_degenerate_metric();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN GCR FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_arnoldi_happy_breakdown();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN ARNOLDI FAILED\n");
        lis_finalize();
        return 1;
    }

    printf(
        "LIS_EIGEN_BREAKDOWN PASSED\n");

    lis_finalize();

    return 0;
}
