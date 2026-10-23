#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <stdio.h>
#include <math.h>

#include "lis.h"
#include "lis_esolver.h"


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
create_diagonal_matrix3(
    LIS_MATRIX *A,
    LIS_SCALAR a0,
    LIS_SCALAR a1,
    LIS_SCALAR a2)
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
        3);
    if( err )
        return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,0,0,a0,*A);
    if( err )
        return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,1,1,a1,*A);
    if( err )
        return err;

    err = lis_matrix_set_value(
        LIS_INS_VALUE,2,2,a2,*A);
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
check_standard_certification_rejects(
    const char *name,
    const char *options,
    LIS_MATRIX A,
    LIS_REAL tol)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_REAL residual;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    residual = 0.0;

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

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    /*
     * The public solve call reports completed API execution with
     * LIS_SUCCESS.  The eigensolver result is read through retcode.
     */
    if( api!=LIS_SUCCESS )
    {
        printf(
            "%s: API return %d, expected LIS_SUCCESS (%d)\n",
            name,
            (int)api,
            (int)LIS_SUCCESS);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,
        &status);
    if( err )
        goto cleanup;

    if( status!=LIS_MAXITER )
    {
        printf(
            "%s: status %d, expected LIS_MAXITER (%d)\n",
            name,
            (int)status,
            (int)LIS_MAXITER);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_residualnorm(
        esolver,
        &residual);
    if( err )
        goto cleanup;

    if( !(residual>tol) )
    {
        printf(
            "%s: physical residual %.17e, expected > %.17e\n",
            name,
            (double)residual,
            (double)tol);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: %s -> LIS_MAXITER, physical residual %.6e\n",
        name,
        (double)residual);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
check_generalized_certification_rejects(
    const char *name,
    const char *options,
    LIS_MATRIX A,
    LIS_MATRIX B,
    LIS_REAL tol)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_REAL residual;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    residual = 0.0;

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

    /*
     * The public solve call reports completed API execution with
     * LIS_SUCCESS.  The eigensolver result is read through retcode.
     */
    if( api!=LIS_SUCCESS )
    {
        printf(
            "%s: API return %d, expected LIS_SUCCESS (%d)\n",
            name,
            (int)api,
            (int)LIS_SUCCESS);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,
        &status);
    if( err )
        goto cleanup;

    if( status!=LIS_MAXITER )
    {
        printf(
            "%s: status %d, expected LIS_MAXITER (%d)\n",
            name,
            (int)status,
            (int)LIS_MAXITER);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_residualnorm(
        esolver,
        &residual);
    if( err )
        goto cleanup;

    if( !(residual>tol) )
    {
        printf(
            "%s: physical residual %.17e, expected > %.17e\n",
            name,
            (double)residual,
            (double)tol);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: %s -> LIS_MAXITER, physical residual %.6e\n",
        name,
        (double)residual);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
check_standard_recovery(
    const char *name,
    const char *options,
    LIS_MATRIX A,
    LIS_SCALAR expected,
    LIS_REAL tol)
{
    LIS_INT err,api,status;
    LIS_REAL residual,error;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    residual = 0.0;

    err = lis_vector_duplicate(A,&x);
    if( err ) goto cleanup;

    err = lis_vector_set_all(1.0,x);
    if( err ) goto cleanup;

    err = lis_esolver_create(&esolver);
    if( err ) goto cleanup;

    err = lis_esolver_set_option(
        (char *)options,esolver);
    if( err ) goto cleanup;

    api = lis_esolve(
        A,x,&evalue,esolver);
    if( api!=LIS_SUCCESS )
    {
        printf(
            "%s: API return %d, expected LIS_SUCCESS\n",
            name,(int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,&status);
    if( err ) goto cleanup;

    if( status!=LIS_SUCCESS )
    {
        printf(
            "%s: status %d, expected LIS_SUCCESS\n",
            name,(int)status);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_residualnorm(
        esolver,&residual);
    if( err ) goto cleanup;

#ifdef _COMPLEX
    error = cabs(evalue-expected);
#else
    error = fabs(evalue-expected);
#endif

    if( !isfinite((double)residual) ||
        residual>tol ||
        !isfinite((double)error) ||
        error>1.0e-5 )
    {
        printf(
            "%s: error %.17e, residual %.17e\n",
            name,
            (double)error,
            (double)residual);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: %s -> residual %.6e\n",
        name,
        (double)residual);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);
    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
check_generalized_recovery(
    const char *name,
    const char *options,
    LIS_MATRIX A,
    LIS_MATRIX B,
    LIS_SCALAR expected,
    LIS_REAL tol)
{
    LIS_INT err,api,status;
    LIS_REAL residual,error;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    residual = 0.0;

    err = lis_vector_duplicate(A,&x);
    if( err ) goto cleanup;

    err = lis_vector_set_all(1.0,x);
    if( err ) goto cleanup;

    err = lis_esolver_create(&esolver);
    if( err ) goto cleanup;

    err = lis_esolver_set_option(
        (char *)options,esolver);
    if( err ) goto cleanup;

    api = lis_gesolve(
        A,B,x,&evalue,esolver);
    if( api!=LIS_SUCCESS )
    {
        printf(
            "%s: API return %d, expected LIS_SUCCESS\n",
            name,(int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,&status);
    if( err ) goto cleanup;

    if( status!=LIS_SUCCESS )
    {
        printf(
            "%s: status %d, expected LIS_SUCCESS\n",
            name,(int)status);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_residualnorm(
        esolver,&residual);
    if( err ) goto cleanup;

#ifdef _COMPLEX
    error = cabs(evalue-expected);
#else
    error = fabs(evalue-expected);
#endif

    if( !isfinite((double)residual) ||
        residual>tol ||
        !isfinite((double)error) ||
        error>1.0e-5 )
    {
        printf(
            "%s: error %.17e, residual %.17e\n",
            name,
            (double)error,
            (double)residual);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: %s -> residual %.6e\n",
        name,
        (double)residual);

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
test_power_physical_certification(void)
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

    /*
     * A large spectral shift makes the shifted residual estimator
     * artificially small after one iteration, although the candidate
     * is still far from an eigenpair of the physical problem.
     */
    err = check_standard_certification_rejects(
        "Power shifted false convergence",
        "-e pi -shift 1.0e6 -emaxiter 1 -etol 1.0e-6",
        A,
        1.0e-6);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_power_physical_certification(void)
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
        1.0,
        1.0);
    if( err )
        goto cleanup;

    err = check_generalized_certification_rejects(
        "Generalized Power shifted false convergence",
        "-e gpi -shift -1.0e6 -emaxiter 1 -etol 1.0e-6",
        A,
        B,
        1.0e-6);

cleanup:
    if( B )
        lis_matrix_destroy(B);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_subspace_physical_certification(void)
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

    err = check_standard_certification_rejects(
        "Subspace shifted false convergence",
        "-e si -ie pi -ss 1 -shift 1.0e6 "
        "-emaxiter 1 -etol 1.0e-6",
        A,
        1.0e-6);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_subspace_physical_certification(void)
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
        1.0,
        1.0);
    if( err )
        goto cleanup;

    err = check_generalized_certification_rejects(
        "Generalized Subspace shifted false convergence",
        "-e gsi -ige gpi -ss 1 -shift -1.0e6 "
        "-emaxiter 1 -etol 1.0e-6",
        A,
        B,
        1.0e-6);

cleanup:
    if( B )
        lis_matrix_destroy(B);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_inverse_shift_recovery(void)
{
    LIS_INT err;
    LIS_MATRIX A;

    A = NULL;

    err = create_diagonal_matrix(
        &A,1.0,2.0);
    if( err ) goto cleanup;

    err = check_standard_recovery(
        "Inverse singular-shift recovery",
        "-e ii -shift 1.0 -emaxiter 8 -etol 1.0e-8",
        A,1.0,1.0e-8);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_inverse_shift_recovery(void)
{
    LIS_INT err;
    LIS_MATRIX A,B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,1.0,2.0);
    if( err ) goto cleanup;

    err = create_diagonal_matrix(
        &B,1.0,1.0);
    if( err ) goto cleanup;

    err = check_generalized_recovery(
        "Generalized Inverse singular-shift recovery",
        "-e gii -shift 1.0 -emaxiter 8 -etol 1.0e-8",
        A,B,1.0,1.0e-8);

cleanup:
    if( B )
        lis_matrix_destroy(B);
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_rqi_shift_recovery(void)
{
    LIS_INT err;
    LIS_MATRIX A;

    A = NULL;

    err = create_diagonal_matrix3(
        &A,1.0,2.0,3.0);
    if( err ) goto cleanup;

    err = check_standard_recovery(
        "RQI singular-shift recovery",
        "-e rqi -emaxiter 8 -etol 1.0e-8",
        A,2.0,1.0e-8);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_grqi_shift_recovery(void)
{
    LIS_INT err;
    LIS_MATRIX A,B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix3(
        &A,1.0,2.0,3.0);
    if( err ) goto cleanup;

    err = create_diagonal_matrix3(
        &B,1.0,1.0,1.0);
    if( err ) goto cleanup;

    err = check_generalized_recovery(
        "Generalized RQI singular-shift recovery",
        "-e grqi -emaxiter 8 -etol 1.0e-8",
        A,B,2.0,1.0e-8);

cleanup:
    if( B )
        lis_matrix_destroy(B);
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_subspace_inverse_shift_recovery(void)
{
    LIS_INT err;
    LIS_MATRIX A;

    A = NULL;

    err = create_diagonal_matrix(
        &A,1.0,2.0);
    if( err ) goto cleanup;

    err = check_standard_recovery(
        "Subspace/II singular-shift recovery",
        "-e si -ie ii -ss 1 -shift 1.0 -emaxiter 8 -etol 1.0e-8",
        A,1.0,1.0e-8);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_subspace_inverse_shift_recovery(void)
{
    LIS_INT err;
    LIS_MATRIX A,B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,1.0,2.0);
    if( err ) goto cleanup;

    err = create_diagonal_matrix(
        &B,1.0,1.0);
    if( err ) goto cleanup;

    err = check_generalized_recovery(
        "Generalized Subspace/GII singular-shift recovery",
        "-e gsi -ige gii -ss 1 -shift 1.0 -emaxiter 8 -etol 1.0e-8",
        A,B,1.0,1.0e-8);

cleanup:
    if( B )
        lis_matrix_destroy(B);
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
check_transformed_gli_recovery(
    const char *name,
    const char *options,
    LIS_MATRIX A,
    LIS_MATRIX B,
    LIS_SCALAR expected,
    LIS_REAL tol)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_REAL residual;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    residual = 0.0;

    err = lis_vector_duplicate(A,&x);
    if( err ) goto cleanup;

    err = lis_vector_set_all(1.0,x);
    if( err ) goto cleanup;

    err = lis_esolver_create(&esolver);
    if( err ) goto cleanup;

    err = lis_esolver_set_option(
        (char *)options,
        esolver);
    if( err ) goto cleanup;

    api = lis_gesolve(
        A,B,x,&evalue,esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "%s: API return %d, expected LIS_SUCCESS (%d)\n",
            name,
            (int)api,
            (int)LIS_SUCCESS);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,
        &status);
    if( err ) goto cleanup;

    if( status!=LIS_SUCCESS )
    {
        printf(
            "%s: status %d, expected LIS_SUCCESS (%d)\n",
            name,
            (int)status,
            (int)LIS_SUCCESS);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_residualnorm(
        esolver,
        &residual);
    if( err ) goto cleanup;

    if( !isfinite((double)residual) || residual>tol )
    {
        printf(
            "%s: physical residual %.17e, expected <= %.17e\n",
            name,
            (double)residual,
            (double)tol);
        err = LIS_FAILS;
        goto cleanup;
    }

#ifdef _COMPLEX
    if( cabs(evalue-expected)>1.0e-6 )
#else
    if( fabs(evalue-expected)>1.0e-6 )
#endif
    {
        printf(
            "%s: eigenvalue mismatch\n",
            name);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: %s -> residual %.6e\n",
        name,
        (double)residual);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
test_gli_shift_invert_singular_shift_recovery(void)
{
    LIS_INT err;
    LIS_MATRIX A,B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,1.0,2.0);
    if( err ) goto cleanup;

    err = create_diagonal_matrix(
        &B,1.0,1.0);
    if( err ) goto cleanup;

    err = check_transformed_gli_recovery(
        "GLI shift-invert singular-shift recovery",
        "-e gli -etransform shift-invert -shift 1.0 "
        "-ss 2 -nev 1 -emaxiter 8 -etol 1.0e-8",
        A,B,1.0,1.0e-8);

cleanup:
    if( B )
        lis_matrix_destroy(B);
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_gli_buckling_singular_shift_recovery(void)
{
    LIS_INT err;
    LIS_MATRIX A,B;

    A = NULL;
    B = NULL;

    err = create_diagonal_matrix(
        &A,1.0,2.0);
    if( err ) goto cleanup;

    err = create_diagonal_matrix(
        &B,1.0,1.0);
    if( err ) goto cleanup;

    err = check_transformed_gli_recovery(
        "GLI buckling singular-shift recovery",
        "-e gli -etransform buckling -shift 1.0 "
        "-ss 2 -nev 1 -emaxiter 8 -etol 1.0e-8",
        A,B,1.0,1.0e-8);

cleanup:
    if( B )
        lis_matrix_destroy(B);
    if( A )
        lis_matrix_destroy(A);

    return err;
}




static LIS_INT
test_inverse_shift_recovery_exhaustion(void)
{
    LIS_INT err;
    LIS_MATRIX A;
    LIS_SCALAR shift;
    LIS_SCALAR delta;

    A = NULL;
    shift = 1.0;

    /*
     * Construct one exact singular mode for each bounded recovery
     * shift: sigma, sigma+delta and sigma-delta.  With an all-ones
     * right-hand side every trial is inconsistent in at least one
     * component, so the finite recovery sequence must exhaust and
     * propagate a breakdown rather than looping or reporting success.
     */
    delta = lis_esolver_shift_perturbation(shift);

    if( delta==0.0 )
    {
        printf(
            "Inverse recovery exhaustion: zero perturbation\n");
        return LIS_FAILS;
    }

    err = create_diagonal_matrix3(
        &A,
        shift,
        shift+delta,
        shift-delta);
    if( err )
        goto cleanup;

    err = check_standard_breakdown(
        "Inverse bounded-shift recovery exhaustion",
        "-e ii -shift 1.0 -emaxiter 2 -etol 1.0e-12",
        A,
        1.0);

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_subspace_gpi_retry_exhaustion(void)
{
    LIS_INT err;
    LIS_MATRIX A;
    LIS_MATRIX B;

    A = NULL;
    B = NULL;

    /*
     * The GPI normalization uses A and is well-defined for A=I.
     * The subsequent B^-1*w solve is deliberately inconsistent:
     * B=diag(1,0) while w has a nonzero second component.
     *
     * This reaches the GPI same-operator retry/rebuild path without
     * introducing a spectral-shift perturbation.
     */
    err = create_diagonal_matrix(
        &A,
        1.0,
        1.0);
    if( err )
        goto cleanup;

    err = create_diagonal_matrix(
        &B,
        1.0,
        0.0);
    if( err )
        goto cleanup;

    err = check_generalized_breakdown(
        "Generalized Subspace/GPI retry exhaustion",
        "-e gsi -ige gpi -ss 1 -emaxiter 2 -etol 1.0e-12",
        A,
        B);

cleanup:
    if( B )
        lis_matrix_destroy(B);

    if( A )
        lis_matrix_destroy(A);

    return err;
}




#define LIS_EIGEN_ORTHO_N 30
#define LIS_EIGEN_ORTHO_SS 20
#define LIS_EIGEN_ORTHO_TOL 1.0e-8


static LIS_REAL
eigen_scalar_magnitude(
    LIS_SCALAR value)
{
#ifdef _COMPLEX
    return cabs(value);
#else
    return fabs(value);
#endif
}


static LIS_INT
create_eigen_orthogonality_matrix(
    LIS_MATRIX *A)
{
    LIS_INT err;
    LIS_INT i;
    LIS_SCALAR value;
    const LIS_SCALAR ratio =
        (LIS_SCALAR)1.6102620275609394;

    *A = NULL;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        A);
    if( err )
        return err;

    err = lis_matrix_set_size(
        *A,
        0,
        LIS_EIGEN_ORTHO_N);
    if( err )
        return err;

    value = 1.0;

    for(i=0;i<LIS_EIGEN_ORTHO_N;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            value,
            *A);
        if( err )
            return err;

        value *= ratio;
    }

    return lis_matrix_assemble(*A);
}


static LIS_INT
create_eigen_identity_matrix(
    LIS_MATRIX *B)
{
    LIS_INT err;
    LIS_INT i;

    *B = NULL;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        B);
    if( err )
        return err;

    err = lis_matrix_set_size(
        *B,
        0,
        LIS_EIGEN_ORTHO_N);
    if( err )
        return err;

    for(i=0;i<LIS_EIGEN_ORTHO_N;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            1.0,
            *B);
        if( err )
            return err;
    }

    return lis_matrix_assemble(*B);
}


static LIS_INT
standard_lanczos_basis_defect(
    LIS_ESOLVER esolver,
    LIS_REAL *defect)
{
    LIS_INT err;
    LIS_INT i,j;
    LIS_SCALAR dot;
    LIS_REAL current;
    LIS_VECTOR *v;

    if( esolver==NULL ||
        esolver->work==NULL ||
        esolver->worklen<
            LIS_EIGEN_ORTHO_SS+2 )
        return LIS_FAILS;

    v = &esolver->work[1];
    *defect = 0.0;

    for(i=1;i<LIS_EIGEN_ORTHO_SS;i++)
    {
        err = lis_vector_dot(
            v[i],
            v[i],
            &dot);
        if( err )
            return err;

        current =
            eigen_scalar_magnitude(
                dot-1.0);

        if( current>*defect )
            *defect = current;

        for(j=1;j<i;j++)
        {
            err = lis_vector_dot(
                v[j],
                v[i],
                &dot);
            if( err )
                return err;

            current =
                eigen_scalar_magnitude(dot);

            if( current>*defect )
                *defect = current;
        }
    }

    return LIS_SUCCESS;
}


static LIS_INT
generalized_lanczos_basis_defect(
    LIS_ESOLVER esolver,
    LIS_REAL *defect)
{
    LIS_INT err;
    LIS_INT i,j;
    LIS_SCALAR dot;
    LIS_SCALAR target;
    LIS_REAL current;
    LIS_VECTOR *w;
    LIS_VECTOR *v;

    if( esolver==NULL ||
        esolver->work==NULL ||
        esolver->worklen<
            4+2*LIS_EIGEN_ORTHO_SS )
        return LIS_FAILS;

    w = &esolver->work[2];
    v = &esolver->work[
        2+LIS_EIGEN_ORTHO_SS];

    *defect = 0.0;

    for(i=1;i<LIS_EIGEN_ORTHO_SS;i++)
    {
        for(j=1;j<LIS_EIGEN_ORTHO_SS;j++)
        {
            err = lis_vector_dot(
                v[i],
                w[j],
                &dot);
            if( err )
                return err;

            target =
                i==j ? 1.0 : 0.0;

            current =
                eigen_scalar_magnitude(
                    dot-target);

            if( current>*defect )
                *defect = current;
        }
    }

    return LIS_SUCCESS;
}




#define LIS_EIGEN_ARNOLDI_ORTHO_TOL 1.0e-8


static LIS_INT
create_arnoldi_orthogonality_matrix(
    LIS_MATRIX *A)
{
    LIS_INT err;
    LIS_INT i;
    LIS_SCALAR value;
    const LIS_SCALAR ratio =
        (LIS_SCALAR)1.05;

    *A = NULL;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        A);
    if( err )
        return err;

    err = lis_matrix_set_size(
        *A,
        0,
        LIS_EIGEN_ORTHO_N);
    if( err )
        return err;

    value = 1.0;

    for(i=0;i<LIS_EIGEN_ORTHO_N;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            value,
            *A);
        if( err )
            return err;

        value *= ratio;
    }

    return lis_matrix_assemble(*A);
}


static LIS_INT
arnoldi_basis_defect(
    LIS_VECTOR *v,
    LIS_INT count,
    LIS_REAL *defect)
{
    LIS_INT err;
    LIS_INT i,j;
    LIS_SCALAR dot;
    LIS_REAL current;

    *defect = 0.0;

    for(i=0;i<count;i++)
    {
        err = lis_vector_dot(
            v[i],
            v[i],
            &dot);
        if( err )
            return err;

        current =
            eigen_scalar_magnitude(
                dot-1.0);

        if( current>*defect )
            *defect = current;

        for(j=0;j<i;j++)
        {
            err = lis_vector_dot(
                v[j],
                v[i],
                &dot);
            if( err )
                return err;

            current =
                eigen_scalar_magnitude(dot);

            if( current>*defect )
                *defect = current;
        }
    }

    return LIS_SUCCESS;
}


static LIS_INT
test_standard_arnoldi_orthogonality(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_REAL defect;
    LIS_VECTOR *v;

    A = NULL;
    x = NULL;
    esolver = NULL;
    defect = 0.0;

    err = create_arnoldi_orthogonality_matrix(
        &A);
    if( err )
        goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err )
        goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    err = lis_esolver_set_option(
        "-e ai -ss 20 -emaxiter 4 "
        "-etol 1.0e-12",
        esolver);
    if( err )
        goto cleanup;

    esolver->options[
        LIS_EOPTIONS_RVAL] = 1;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Standard Arnoldi orthogonality: "
            "API return %d\n",
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    if( esolver->work==NULL ||
        esolver->worklen<
            LIS_EIGEN_ORTHO_SS+2 )
    {
        err = LIS_FAILS;
        goto cleanup;
    }

    v = &esolver->work[1];

    err = arnoldi_basis_defect(
        v,
        LIS_EIGEN_ORTHO_SS,
        &defect);
    if( err )
        goto cleanup;

    if( !isfinite((double)defect) ||
        defect>LIS_EIGEN_ARNOLDI_ORTHO_TOL )
    {
        printf(
            "Standard Arnoldi orthogonality "
            "defect %.17e, expected <= %.17e\n",
            (double)defect,
            (double)LIS_EIGEN_ARNOLDI_ORTHO_TOL);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Standard Arnoldi "
        "orthogonality defect %.6e\n",
        (double)defect);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);
    if( x )
        lis_vector_destroy(x);
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_arnoldi_orthogonality(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_MATRIX A;
    LIS_MATRIX B;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_REAL defect;
    LIS_VECTOR *v;

    A = NULL;
    B = NULL;
    x = NULL;
    esolver = NULL;
    defect = 0.0;

    err = create_arnoldi_orthogonality_matrix(
        &A);
    if( err )
        goto cleanup;

    err = create_eigen_identity_matrix(
        &B);
    if( err )
        goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err )
        goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    err = lis_esolver_set_option(
        "-e gai -ss 20 -emaxiter 4 "
        "-etol 1.0e-12",
        esolver);
    if( err )
        goto cleanup;

    esolver->options[
        LIS_EOPTIONS_RVAL] = 1;

    api = lis_gesolve(
        A,
        B,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Generalized Arnoldi orthogonality: "
            "API return %d\n",
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    if( esolver->work==NULL ||
        esolver->worklen<
            LIS_EIGEN_ORTHO_SS+3 )
    {
        err = LIS_FAILS;
        goto cleanup;
    }

    v = &esolver->work[2];

    err = arnoldi_basis_defect(
        v,
        LIS_EIGEN_ORTHO_SS,
        &defect);
    if( err )
        goto cleanup;

    if( !isfinite((double)defect) ||
        defect>LIS_EIGEN_ARNOLDI_ORTHO_TOL )
    {
        printf(
            "Generalized Arnoldi orthogonality "
            "defect %.17e, expected <= %.17e\n",
            (double)defect,
            (double)LIS_EIGEN_ARNOLDI_ORTHO_TOL);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Generalized Arnoldi "
        "orthogonality defect %.6e\n",
        (double)defect);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);
    if( x )
        lis_vector_destroy(x);
    if( B )
        lis_matrix_destroy(B);
    if( A )
        lis_matrix_destroy(A);

    return err;
}




#define LIS_EIGEN_LOCK_N 8
#define LIS_EIGEN_LOCK_SS 3
#define LIS_EIGEN_LOCK_TOL 1.0e-6
#define LIS_EIGEN_LOCK_OVERLAP_TOL 1.0e-10
#define LIS_EIGEN_LOCK_EVALUE_TOL 1.0e-5


static const LIS_REAL lis_eigen_lock_power_expected[
    LIS_EIGEN_LOCK_SS] = {
        2.00,
        1.80,
        1.60
    };

static const LIS_REAL lis_eigen_lock_inverse_expected[
    LIS_EIGEN_LOCK_SS] = {
        0.50,
        0.60,
        0.70
    };


static LIS_INT
locking_evalue_sequence(
    const char *name,
    LIS_ESOLVER esolver,
    const LIS_REAL *expected)
{
    LIS_INT i;
    LIS_REAL error;

    for(i=0;i<LIS_EIGEN_LOCK_SS;i++)
    {
        error = eigen_scalar_magnitude(
            esolver->evalue[i] -
            (LIS_SCALAR)expected[i]);

        if( !isfinite((double)error) ||
            error>LIS_EIGEN_LOCK_EVALUE_TOL )
        {
            printf(
                "%s: eigenvalue[%d] error %.17e, expected %.17e\n",
                name,
                (int)i,
                (double)error,
                (double)expected[i]);
            return LIS_FAILS;
        }
    }

    printf(
        "PASS: %s eigenvalue sequence\n",
        name);

    return LIS_SUCCESS;
}


static LIS_INT
create_locking_matrix(
    LIS_MATRIX *A)
{
    static const LIS_REAL diag[LIS_EIGEN_LOCK_N] = {
        0.50,
        0.60,
        0.70,
        1.00,
        1.30,
        1.60,
        1.80,
        2.00
    };

    LIS_INT err;
    LIS_INT i;

    *A = NULL;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        A);
    if( err )
        return err;

    err = lis_matrix_set_size(
        *A,
        0,
        LIS_EIGEN_LOCK_N);
    if( err )
        return err;

    for(i=0;i<LIS_EIGEN_LOCK_N;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            diag[i],
            *A);
        if( err )
            return err;
    }

    return lis_matrix_assemble(*A);
}


static LIS_INT
create_locking_identity(
    LIS_MATRIX *B)
{
    LIS_INT err;
    LIS_INT i;

    *B = NULL;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        B);
    if( err )
        return err;

    err = lis_matrix_set_size(
        *B,
        0,
        LIS_EIGEN_LOCK_N);
    if( err )
        return err;

    for(i=0;i<LIS_EIGEN_LOCK_N;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            1.0,
            *B);
        if( err )
            return err;
    }

    return lis_matrix_assemble(*B);
}


static LIS_INT
locking_euclidean_overlap(
    LIS_ESOLVER esolver,
    LIS_REAL *max_overlap)
{
    LIS_INT err;
    LIS_INT i,j;
    LIS_REAL ni,nj;
    LIS_REAL value;
    LIS_SCALAR dot;

    *max_overlap = 0.0;

    for(i=0;i<LIS_EIGEN_LOCK_SS;i++)
    {
        err = lis_vector_nrm2(
            esolver->evector[i],
            &ni);
        if( err )
            return err;

        if( !(ni>0.0) )
            return LIS_BREAKDOWN;

        for(j=0;j<i;j++)
        {
            err = lis_vector_nrm2(
                esolver->evector[j],
                &nj);
            if( err )
                return err;

            if( !(nj>0.0) )
                return LIS_BREAKDOWN;

            err = lis_vector_dot(
                esolver->evector[j],
                esolver->evector[i],
                &dot);
            if( err )
                return err;

            value =
                eigen_scalar_magnitude(dot) /
                (ni*nj);

            if( value>*max_overlap )
                *max_overlap = value;
        }
    }

    return LIS_SUCCESS;
}


static LIS_INT
locking_b_overlap(
    LIS_MATRIX B,
    LIS_ESOLVER esolver,
    LIS_REAL *max_overlap)
{
    LIS_INT err;
    LIS_INT i,j;
    LIS_REAL value;
    LIS_SCALAR numerator;
    LIS_SCALAR di,dj;
    LIS_VECTOR Bi;
    LIS_VECTOR Bj;

    Bi = NULL;
    Bj = NULL;
    *max_overlap = 0.0;

    err = lis_vector_duplicate(
        B,
        &Bi);
    if( err )
        goto cleanup;

    err = lis_vector_duplicate(
        B,
        &Bj);
    if( err )
        goto cleanup;

    for(i=0;i<LIS_EIGEN_LOCK_SS;i++)
    {
        err = lis_matvec(
            B,
            esolver->evector[i],
            Bi);
        if( err )
            goto cleanup;

        err = lis_vector_dot(
            esolver->evector[i],
            Bi,
            &di);
        if( err )
            goto cleanup;

        if( !lis_esolver_is_valid_denominator(di) )
        {
            err = LIS_BREAKDOWN;
            goto cleanup;
        }

        for(j=0;j<i;j++)
        {
            err = lis_matvec(
                B,
                esolver->evector[j],
                Bj);
            if( err )
                goto cleanup;

            err = lis_vector_dot(
                esolver->evector[j],
                Bj,
                &dj);
            if( err )
                goto cleanup;

            if( !lis_esolver_is_valid_denominator(dj) )
            {
                err = LIS_BREAKDOWN;
                goto cleanup;
            }

            err = lis_vector_dot(
                esolver->evector[j],
                Bi,
                &numerator);
            if( err )
                goto cleanup;

            value =
                eigen_scalar_magnitude(numerator) /
                sqrt(
                    eigen_scalar_magnitude(di) *
                    eigen_scalar_magnitude(dj));

            if( value>*max_overlap )
                *max_overlap = value;
        }
    }

    err = LIS_SUCCESS;

cleanup:
    if( Bj )
        lis_vector_destroy(Bj);
    if( Bi )
        lis_vector_destroy(Bi);

    return err;
}


static LIS_INT
check_locked_standard_subspace(
    const char *name,
    const char *inner,
    const LIS_REAL *expected,
    LIS_MATRIX A)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_REAL overlap;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    char options[256];

    x = NULL;
    esolver = NULL;
    overlap = 0.0;

    err = lis_vector_duplicate(A,&x);
    if( err )
        goto cleanup;

    err = lis_vector_set_all(1.0,x);
    if( err )
        goto cleanup;

    err = lis_esolver_create(&esolver);
    if( err )
        goto cleanup;

    snprintf(
        options,
        sizeof(options),
        "-e si -ie %s -ss %d "
        "-emaxiter 600 -etol %.17e",
        inner,
        LIS_EIGEN_LOCK_SS,
        (double)LIS_EIGEN_LOCK_TOL);

    err = lis_esolver_set_option(
        options,
        esolver);
    if( err )
        goto cleanup;

    api = lis_esolve(
        A,x,&evalue,esolver);
    if( api!=LIS_SUCCESS )
    {
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,&status);
    if( err )
        goto cleanup;

    if( status!=LIS_SUCCESS )
    {
        printf(
            "%s: status %d, expected LIS_SUCCESS\n",
            name,(int)status);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = locking_euclidean_overlap(
        esolver,&overlap);
    if( err )
        goto cleanup;

    if( !isfinite((double)overlap) ||
        overlap>LIS_EIGEN_LOCK_OVERLAP_TOL )
    {
        printf(
            "%s: overlap %.17e, expected <= %.17e\n",
            name,
            (double)overlap,
            (double)LIS_EIGEN_LOCK_OVERLAP_TOL);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = locking_evalue_sequence(
        name,
        esolver,
        expected);
    if( err )
        goto cleanup;

    printf(
        "PASS: %s locking overlap %.6e\n",
        name,
        (double)overlap);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);
    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
check_locked_generalized_subspace(
    const char *name,
    const char *inner,
    const LIS_REAL *expected,
    LIS_MATRIX A,
    LIS_MATRIX B)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_REAL overlap;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    char options[256];

    x = NULL;
    esolver = NULL;
    overlap = 0.0;

    err = lis_vector_duplicate(A,&x);
    if( err )
        goto cleanup;

    err = lis_vector_set_all(1.0,x);
    if( err )
        goto cleanup;

    err = lis_esolver_create(&esolver);
    if( err )
        goto cleanup;

    snprintf(
        options,
        sizeof(options),
        "-e gsi -ige %s -ss %d "
        "-emaxiter 600 -etol %.17e",
        inner,
        LIS_EIGEN_LOCK_SS,
        (double)LIS_EIGEN_LOCK_TOL);

    err = lis_esolver_set_option(
        options,
        esolver);
    if( err )
        goto cleanup;

    api = lis_gesolve(
        A,B,x,&evalue,esolver);
    if( api!=LIS_SUCCESS )
    {
        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,&status);
    if( err )
        goto cleanup;

    if( status!=LIS_SUCCESS )
    {
        printf(
            "%s: status %d, expected LIS_SUCCESS\n",
            name,(int)status);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = locking_b_overlap(
        B,esolver,&overlap);
    if( err )
        goto cleanup;

    if( !isfinite((double)overlap) ||
        overlap>LIS_EIGEN_LOCK_OVERLAP_TOL )
    {
        printf(
            "%s: B-overlap %.17e, expected <= %.17e\n",
            name,
            (double)overlap,
            (double)LIS_EIGEN_LOCK_OVERLAP_TOL);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = locking_evalue_sequence(
        name,
        esolver,
        expected);
    if( err )
        goto cleanup;

    printf(
        "PASS: %s B-locking overlap %.6e\n",
        name,
        (double)overlap);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);
    if( x )
        lis_vector_destroy(x);

    return err;
}


static LIS_INT
test_subspace_locked_modes(void)
{
    LIS_INT err;
    LIS_MATRIX A;
    LIS_MATRIX B;

    A = NULL;
    B = NULL;

    err = create_locking_matrix(&A);
    if( err )
        goto cleanup;

    err = create_locking_identity(&B);
    if( err )
        goto cleanup;

    err = check_locked_standard_subspace(
        "Subspace/PI",
        "pi",
        lis_eigen_lock_power_expected,
        A);
    if( err )
        goto cleanup;

    err = check_locked_standard_subspace(
        "Subspace/II",
        "ii",
        lis_eigen_lock_inverse_expected,
        A);
    if( err )
        goto cleanup;

    err = check_locked_generalized_subspace(
        "Generalized Subspace/GPI",
        "gpi",
        lis_eigen_lock_power_expected,
        A,
        B);
    if( err )
        goto cleanup;

    err = check_locked_generalized_subspace(
        "Generalized Subspace/GII",
        "gii",
        lis_eigen_lock_inverse_expected,
        A,
        B);
    if( err )
        goto cleanup;

    printf(
        "PASS: Subspace locked-mode regressions\n");

    err = LIS_SUCCESS;

cleanup:
    if( B )
        lis_matrix_destroy(B);
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_standard_lanczos_orthogonality(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_REAL defect;

    A = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    defect = 0.0;

    err = create_eigen_orthogonality_matrix(
        &A);
    if( err )
        goto cleanup;

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
        "-e li -ss 20 -emaxiter 4 "
        "-etol 1.0e-12",
        esolver);
    if( err )
        goto cleanup;

    esolver->options[
        LIS_EOPTIONS_RVAL] = 1;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Standard Lanczos orthogonality: "
            "API return %d\n",
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = standard_lanczos_basis_defect(
        esolver,
        &defect);
    if( err )
        goto cleanup;

    if( !isfinite((double)defect) ||
        defect>LIS_EIGEN_ORTHO_TOL )
    {
        printf(
            "Standard Lanczos orthogonality "
            "defect %.17e, expected <= %.17e\n",
            (double)defect,
            (double)LIS_EIGEN_ORTHO_TOL);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Standard Lanczos "
        "orthogonality defect %.6e\n",
        (double)defect);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);
    if( x )
        lis_vector_destroy(x);
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_generalized_lanczos_orthogonality(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_MATRIX A;
    LIS_MATRIX B;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_REAL defect;

    A = NULL;
    B = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    defect = 0.0;

    err = create_eigen_orthogonality_matrix(
        &A);
    if( err )
        goto cleanup;

    err = create_eigen_identity_matrix(
        &B);
    if( err )
        goto cleanup;

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
        "-e gli -ss 20 -emaxiter 4 "
        "-etol 1.0e-12",
        esolver);
    if( err )
        goto cleanup;

    esolver->options[
        LIS_EOPTIONS_RVAL] = 1;

    api = lis_gesolve(
        A,
        B,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Generalized Lanczos orthogonality: "
            "API return %d\n",
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = generalized_lanczos_basis_defect(
        esolver,
        &defect);
    if( err )
        goto cleanup;

    if( !isfinite((double)defect) ||
        defect>LIS_EIGEN_ORTHO_TOL )
    {
        printf(
            "Generalized Lanczos B-orthogonality "
            "defect %.17e, expected <= %.17e\n",
            (double)defect,
            (double)LIS_EIGEN_ORTHO_TOL);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Generalized Lanczos "
        "B-orthogonality defect %.6e\n",
        (double)defect);

    err = LIS_SUCCESS;

cleanup:
    if( esolver )
        lis_esolver_destroy(esolver);
    if( x )
        lis_vector_destroy(x);
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


    err = test_power_physical_certification();
    if( err )
    {
        printf(
            "LIS_EIGEN_CERTIFICATION POWER FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_power_physical_certification();
    if( err )
    {
        printf(
            "LIS_EIGEN_CERTIFICATION GPOWER FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_subspace_physical_certification();
    if( err )
    {
        printf(
            "LIS_EIGEN_CERTIFICATION SUBSPACE FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_subspace_physical_certification();
    if( err )
    {
        printf(
            "LIS_EIGEN_CERTIFICATION GSUBSPACE FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_inverse_shift_recovery();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY INVERSE FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_inverse_shift_recovery();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY GINVERSE FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_rqi_shift_recovery();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY RQI FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_grqi_shift_recovery();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY GRQI FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_subspace_inverse_shift_recovery();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY SUBSPACE_II FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_subspace_inverse_shift_recovery();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY GSUBSPACE_GII FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_gli_shift_invert_singular_shift_recovery();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY GLI_SHIFT_INVERT FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_gli_buckling_singular_shift_recovery();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY GLI_BUCKLING FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_inverse_shift_recovery_exhaustion();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY EXHAUSTION FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_subspace_gpi_retry_exhaustion();
    if( err )
    {
        printf(
            "LIS_EIGEN_RECOVERY GSUBSPACE_GPI_RETRY FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_subspace_locked_modes();
    if( err )
    {
        printf(
            "LIS_EIGEN_LOCKING SUBSPACE FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_standard_lanczos_orthogonality();
    if( err )
    {
        printf(
            "LIS_EIGEN_ORTHOGONALITY LANCZOS FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_lanczos_orthogonality();
    if( err )
    {
        printf(
            "LIS_EIGEN_ORTHOGONALITY "
            "GENERALIZED LANCZOS FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_standard_arnoldi_orthogonality();
    if( err )
    {
        printf(
            "LIS_EIGEN_ORTHOGONALITY ARNOLDI FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_generalized_arnoldi_orthogonality();
    if( err )
    {
        printf(
            "LIS_EIGEN_ORTHOGONALITY "
            "GENERALIZED ARNOLDI FAILED\n");
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
