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
test_lanczos_block_option_contract(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    A = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;

    /*
     * The historical scalar path must remain the default.
     */
    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    if( esolver->block_size!=1 )
    {
        printf(
            "Block Lanczos option: default block size %d, "
            "expected 1\n",
            (int)esolver->block_size);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_set_option(
        "-e li -ss 2 -eblock 1",
        esolver);
    if( err )
        goto cleanup;

    if( esolver->block_size!=1 )
    {
        printf(
            "Block Lanczos option: explicit scalar block "
            "size was not preserved\n");

        err = LIS_FAILS;
        goto cleanup;
    }

    lis_esolver_destroy(esolver);
    esolver = NULL;

    /*
     * Invalid non-positive block sizes are rejected by the parser.
     */
    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    err = lis_esolver_set_option(
        "-eblock 0",
        esolver);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "Block Lanczos option: -eblock 0 returned %d, "
            "expected LIS_ERR_ILL_ARG (%d)\n",
            (int)err,
            (int)LIS_ERR_ILL_ARG);

        err = LIS_FAILS;
        goto cleanup;
    }

    lis_esolver_destroy(esolver);
    esolver = NULL;

    /*
     * Build a minimal matrix used to exercise the runtime contract.
     */
    err = create_diagonal_matrix(
        &A,
        1.0,
        2.0);
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

    /*
     * A block larger than the total Krylov subspace is invalid.
     */
    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    err = lis_esolver_set_option(
        "-e li -ss 2 -eblock 3 -rval true",
        esolver);
    if( err )
        goto cleanup;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_ERR_ILL_ARG )
    {
        printf(
            "Block Lanczos option: block size above ss "
            "returned %d, expected LIS_ERR_ILL_ARG (%d)\n",
            (int)api,
            (int)LIS_ERR_ILL_ARG);

        err = LIS_FAILS;
        goto cleanup;
    }

    lis_esolver_destroy(esolver);
    esolver = NULL;

    /*
     * Block mode is not silently accepted by other eigensolvers.
     */
    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    err = lis_esolver_set_option(
        "-e ai -ss 2 -eblock 2 -rval true",
        esolver);
    if( err )
        goto cleanup;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_ERR_NOT_IMPLEMENTED )
    {
        printf(
            "Block Lanczos option: non-Lanczos block mode "
            "returned %d, expected LIS_ERR_NOT_IMPLEMENTED (%d)\n",
            (int)api,
            (int)LIS_ERR_NOT_IMPLEMENTED);

        err = LIS_FAILS;
        goto cleanup;
    }

    lis_esolver_destroy(esolver);
    esolver = NULL;

    /*
     * A valid standard Block Lanczos request must reach the dedicated
     * block driver.  BL0-BL2 deliberately stop there until the block
     * recurrence is connected.
     */
    err = lis_esolver_create(
        &esolver);
    if( err )
        goto cleanup;

    err = lis_esolver_set_option(
        "-e li -ss 2 -eblock 2 -rval true",
        esolver);
    if( err )
        goto cleanup;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    /*
     * The eigensolver execution result is reported through
     * esolver->retcode.  The public solve API returns LIS_SUCCESS
     * after dispatch, consistent with the existing eigensolver
     * execution contract.
     */
    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos foundation: API return %d, "
            "expected LIS_SUCCESS (%d)\n",
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

    if( status!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos foundation: status %d, "
            "expected LIS_SUCCESS (%d)\n",
            (int)status,
            (int)LIS_SUCCESS);

        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Block Lanczos BL0-BL2 option, dispatch, "
        "and layout contract\n");

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
        "-e li -ss 20 -eblock 1 -emaxiter 4 "
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
test_lanczos_clustered_basis_stability(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_MATRIX A;
    LIS_MATRIX B;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_REAL standard_defect;
    LIS_REAL generalized_defect;

    A = NULL;
    B = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    standard_defect = 0.0;
    generalized_defect = 0.0;

    /*
     * Reuse the geometric 1.05 spectrum used by the Arnoldi stress
     * regression.  This stresses basis independence while preserving
     * the single-vector Lanczos contract.
     */
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
            "Clustered standard Lanczos basis: "
            "API return %d\n",
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = standard_lanczos_basis_defect(
        esolver,
        &standard_defect);
    if( err )
        goto cleanup;

    if( !isfinite((double)standard_defect) ||
        standard_defect>LIS_EIGEN_ORTHO_TOL )
    {
        printf(
            "Clustered standard Lanczos basis defect "
            "%.17e, expected <= %.17e\n",
            (double)standard_defect,
            (double)LIS_EIGEN_ORTHO_TOL);
        err = LIS_FAILS;
        goto cleanup;
    }

    lis_esolver_destroy(esolver);
    esolver = NULL;

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
            "Clustered generalized Lanczos basis: "
            "API return %d\n",
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    err = generalized_lanczos_basis_defect(
        esolver,
        &generalized_defect);
    if( err )
        goto cleanup;

    if( !isfinite((double)generalized_defect) ||
        generalized_defect>LIS_EIGEN_ORTHO_TOL )
    {
        printf(
            "Clustered generalized Lanczos B-basis defect "
            "%.17e, expected <= %.17e\n",
            (double)generalized_defect,
            (double)LIS_EIGEN_ORTHO_TOL);
        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Lanczos clustered basis stability "
        "standard %.6e generalized %.6e\n",
        (double)standard_defect,
        (double)generalized_defect);

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
test_transformed_generalized_lanczos_clustered_spectrum(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT i;
    LIS_MATRIX A;
    LIS_MATRIX B;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_REAL error;
    const LIS_REAL expected[3] = {
        1.0,
        1.05,
        1.1025
    };

    A = NULL;
    B = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;

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
        "-e gli -ss 20 -emaxiter 20 "
        "-etol 1.0e-10",
        esolver);
    if( err )
        goto cleanup;

    /*
     * Explicitly exercise lis_egli_transformed().
     * A zero spectral shift is nonsingular for this positive
     * diagonal spectrum and targets the smallest eigenvalues.
     */
    esolver->transform =
        LIS_EIGEN_TRANSFORM_SHIFT_INVERT;
    esolver->nev = 3;
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
            "Transformed generalized Lanczos clustered spectrum: "
            "API return %d\n",
            (int)api);
        err = LIS_FAILS;
        goto cleanup;
    }

    if( esolver->retcode!=LIS_SUCCESS )
    {
        printf(
            "Transformed generalized Lanczos clustered spectrum: "
            "status %d, expected LIS_SUCCESS\n",
            (int)esolver->retcode);
        err = LIS_FAILS;
        goto cleanup;
    }

    for(i=0;i<3;i++)
    {
        if( !lis_esolver_is_finite_scalar(
                esolver->evalue[i]) )
        {
            printf(
                "Transformed generalized Lanczos clustered spectrum: "
                "non-finite eigenvalue[%d]\n",
                (int)i);
            err = LIS_FAILS;
            goto cleanup;
        }

        error = eigen_scalar_magnitude(
            esolver->evalue[i] -
            (LIS_SCALAR)expected[i]);

        if( !isfinite((double)error) ||
            error>1.0e-6 )
        {
            printf(
                "Transformed generalized Lanczos clustered spectrum: "
                "eigenvalue[%d] error %.17e, expected <= %.17e\n",
                (int)i,
                (double)error,
                1.0e-6);
            err = LIS_FAILS;
            goto cleanup;
        }

        if( !isfinite((double)esolver->resid[i]) ||
            esolver->resid[i]>1.0e-8 )
        {
            printf(
                "Transformed generalized Lanczos clustered spectrum: "
                "residual[%d] %.17e, expected <= %.17e\n",
                (int)i,
                (double)esolver->resid[i],
                1.0e-8);
            err = LIS_FAILS;
            goto cleanup;
        }
    }

    printf(
        "PASS: Transformed generalized Lanczos "
        "clustered spectrum\n");

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


static LIS_INT
test_block_lanczos_initial_raw(void)
{
    LIS_INT err;
    LIS_INT i,j,k;
    LIS_INT is,ie;
    LIS_INT index;
    LIS_INT block_size;
    LIS_INT block_rank;
    LIS_INT global_n;
    LIS_REAL nrm2;
    LIS_SCALAR got;
    LIS_SCALAR expected;
    LIS_SCALAR dot;
    LIS_VECTOR W[3];

    block_size = 3;
    block_rank = 3;
    global_n = 8;

    for(j=0;j<3;j++)
    {
        W[j] = NULL;
    }

    for(j=0;j<block_rank;j++)
    {
        err = lis_vector_create(
            LIS_COMM_WORLD,
            &W[j]);

        if( err )
        {
            goto cleanup;
        }

        err = lis_vector_set_size(
            W[j],
            0,
            global_n);

        if( err )
        {
            goto cleanup;
        }
    }

    err = lis_lanczos_block_initial_raw(
        block_size,
        block_rank,
        W);

    if( err )
    {
        goto cleanup;
    }

    for(j=0;j<block_rank;j++)
    {
        err = lis_vector_get_range(
            W[j],
            &is,
            &ie);

        if( err )
        {
            goto cleanup;
        }

        for(i=is;i<ie;i++)
        {
            index = i + W[j]->origin;

            err = lis_vector_get_value(
                W[j],
                index,
                &got);

            if( err )
            {
                goto cleanup;
            }

            expected =
                ((i % block_size)==j)
                ? 1.0
                : 0.0;

            if( got!=expected )
            {
                printf(
                    "Block Lanczos raw initial block mismatch: "
                    "column=%D global_index=%D\n",
                    j,
                    i);

                err = LIS_FAILS;
                goto cleanup;
            }
        }

        err = lis_vector_nrm2(
            W[j],
            &nrm2);

        if( err )
        {
            goto cleanup;
        }

        if( nrm2<=0.0 )
        {
            printf(
                "Block Lanczos raw initial block has zero column %D\n",
                j);

            err = LIS_FAILS;
            goto cleanup;
        }
    }

    for(j=0;j<block_rank;j++)
    {
        for(k=j+1;k<block_rank;k++)
        {
            err = lis_vector_dot(
                W[j],
                W[k],
                &dot);

            if( err )
            {
                goto cleanup;
            }

            if( dot!=0.0 )
            {
                printf(
                    "Block Lanczos raw initial block columns "
                    "%D and %D are not orthogonal\n",
                    j,
                    k);

                err = LIS_FAILS;
                goto cleanup;
            }
        }
    }

    printf(
        "PASS: Block Lanczos deterministic raw initial block\n");

    err = LIS_SUCCESS;

cleanup:

    for(j=0;j<3;j++)
    {
        if( W[j] )
        {
            lis_vector_destroy(W[j]);
        }
    }

    return err;
}


static LIS_INT
test_block_lanczos_projected_storage(void)
{
    LIS_INT err;
    LIS_INT i;
    LIS_INT ld;

    LIS_SCALAR t[81];
    LIS_SCALAR alpha0[16];
    LIS_SCALAR alpha1[9];
    LIS_SCALAR alpha2[4];
    LIS_SCALAR beta1[12];
    LIS_SCALAR beta2[6];

    LIS_SCALAR upper_expected;

    ld = 9;

    for(i=0;i<81;i++)
    {
        t[i] = 0.0;
    }

    for(i=0;i<16;i++)
    {
        alpha0[i] = 0.0;
    }

    for(i=0;i<9;i++)
    {
        alpha1[i] = 0.0;
    }

    for(i=0;i<4;i++)
    {
        alpha2[i] = 0.0;
    }

    alpha0[0] = 1.0;
    alpha0[5] = 2.0;
    alpha0[10] = 3.0;
    alpha0[15] = 4.0;
    alpha0[1] = 0.5;
    alpha0[4] = 0.5;

    alpha1[0] = 5.0;
    alpha1[4] = 6.0;
    alpha1[8] = 7.0;

    alpha2[0] = 8.0;
    alpha2[3] = 9.0;
    alpha2[1] = 0.25;
    alpha2[2] = 0.25;

    for(i=0;i<12;i++)
    {
        beta1[i] =
            (LIS_SCALAR)(10+i);
    }

    for(i=0;i<6;i++)
    {
        beta2[i] =
            (LIS_SCALAR)(30+i);
    }

#ifdef _COMPLEX
    beta1[0] =
        (LIS_SCALAR)(1.0 + 2.0*I);
#endif

    err =
        lis_lanczos_block_projected_set_alpha(
            ld,
            t,
            0,
            4,
            alpha0);

    if( err )
    {
        return err;
    }

    err =
        lis_lanczos_block_projected_set_beta(
            ld,
            t,
            0,
            4,
            4,
            3,
            beta1);

    if( err )
    {
        return err;
    }

    err =
        lis_lanczos_block_projected_set_alpha(
            ld,
            t,
            4,
            3,
            alpha1);

    if( err )
    {
        return err;
    }

    err =
        lis_lanczos_block_projected_set_beta(
            ld,
            t,
            4,
            3,
            7,
            2,
            beta2);

    if( err )
    {
        return err;
    }

    err =
        lis_lanczos_block_projected_set_alpha(
            ld,
            t,
            7,
            2,
            alpha2);

    if( err )
    {
        return err;
    }


    /*
     * Alpha0 occupies rows/columns 0..3.
     */
    if( t[0*ld+0]!=alpha0[0] ||
        t[1*ld+1]!=alpha0[5] ||
        t[3*ld+3]!=alpha0[15] )
    {
        printf(
            "Block Lanczos projected Alpha0 placement failed\n");

        return LIS_FAILS;
    }


    /*
     * Beta1 has dimensions 3 x 4 and occupies
     * rows 4..6, columns 0..3.
     */
    if( t[(4+2)*ld+(0+3)] !=
        beta1[2*4+3] )
    {
        printf(
            "Block Lanczos rectangular Beta1 lower placement failed\n");

        return LIS_FAILS;
    }

#ifdef _COMPLEX
    upper_expected =
        conj(beta1[0]);
#else
    upper_expected =
        beta1[0];
#endif

    if( t[0*ld+4] != upper_expected )
    {
        printf(
            "Block Lanczos Beta1 Hermitian placement failed\n");

        return LIS_FAILS;
    }


    /*
     * Alpha1 occupies rows/columns 4..6.
     */
    if( t[4*ld+4]!=5.0 ||
        t[5*ld+5]!=6.0 ||
        t[6*ld+6]!=7.0 )
    {
        printf(
            "Block Lanczos projected Alpha1 placement failed\n");

        return LIS_FAILS;
    }


    /*
     * Beta2 has dimensions 2 x 3 and occupies
     * rows 7..8, columns 4..6.
     */
    if( t[7*ld+4]!=beta2[0] ||
        t[8*ld+6]!=beta2[5] )
    {
        printf(
            "Block Lanczos rectangular Beta2 placement failed\n");

        return LIS_FAILS;
    }

#ifdef _COMPLEX
    upper_expected =
        conj(beta2[5]);
#else
    upper_expected =
        beta2[5];
#endif

    if( t[6*ld+8] != upper_expected )
    {
        printf(
            "Block Lanczos Beta2 transpose placement failed\n");

        return LIS_FAILS;
    }


    /*
     * Alpha2 occupies rows/columns 7..8.
     */
    if( t[7*ld+7]!=8.0 ||
        t[8*ld+8]!=9.0 ||
        t[7*ld+8]!=0.25 ||
        t[8*ld+7]!=0.25 )
    {
        printf(
            "Block Lanczos projected Alpha2 placement failed\n");

        return LIS_FAILS;
    }


    /*
     * Entry outside the block-tridiagonal structure must
     * remain zero.
     */
    if( t[0*ld+8]!=0.0 ||
        t[8*ld+0]!=0.0 )
    {
        printf(
            "Block Lanczos projected matrix fill exceeded "
            "block-tridiagonal structure\n");

        return LIS_FAILS;
    }


    /*
     * Invalid placement must be rejected.
     */
    err =
        lis_lanczos_block_projected_set_beta(
            ld,
            t,
            0,
            4,
            8,
            2,
            beta2);

    if( err!=LIS_ERR_ILL_ARG )
    {
        printf(
            "Block Lanczos invalid Beta placement was not rejected\n");

        return LIS_FAILS;
    }


    printf(
        "PASS: Block Lanczos variable-rank projected storage "
        "4 -> 3 -> 2\n");

    return LIS_SUCCESS;
}



/*
 * Exercise the first complete two-block Block Lanczos transition
 * through the public eigensolver path.
 *
 * For block_size=2 the deterministic raw block is
 *
 *     [1 0]
 *     [0 1]
 *     [1 0]
 *     [0 1]
 *
 * before normalization.
 *
 * With
 *
 *     A = diag(1,2,3,4)
 *
 * neither normalized Q0 column is invariant under A. Therefore
 * W1 has numerical rank two:
 *
 *     W1(:,0) ~ (-1, 0,  1, 0)
 *     W1(:,1) ~ ( 0,-1,  0, 1)
 *
 * and the Q1/Beta1 path must be entered.
 *
 * ss=4 exercises a full 2x2 Beta1.
 * ss=3 exercises the final partial 1x2 Beta1 path.
 *
 * The driver still deliberately terminates after this transition
 * with LIS_ERR_NOT_IMPLEMENTED. Reaching that stop is therefore
 * the expected API result for this stage.
 */
static LIS_INT
test_block_lanczos_two_block_runtime_case(
    LIS_INT ss)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT i;
    LIS_INT is,ie;

    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;

    LIS_SCALAR evalue;

    A = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_set_size(
        A,
        0,
        4);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if( err )
    {
        goto cleanup;
    }

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            (LIS_SCALAR)(i+1),
            A);

        if( err )
        {
            goto cleanup;
        }
    }

    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_assemble(
        A);

    if( err )
    {
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A,
        &x);

    if( err )
    {
        goto cleanup;
    }

    err = lis_vector_set_all(
        1.0,
        x);

    if( err )
    {
        goto cleanup;
    }

    err = lis_esolver_create(
        &esolver);

    if( err )
    {
        goto cleanup;
    }

    if( ss==4 )
    {
        err = lis_esolver_set_option(
            "-e li -ss 4 -eblock 2 -rval true",
            esolver);
    }
    else if( ss==3 )
    {
        err = lis_esolver_set_option(
            "-e li -ss 3 -eblock 2 -rval true",
            esolver);
    }
    else
    {
        err = LIS_ERR_ILL_ARG;
    }

    if( err )
    {
        goto cleanup;
    }

    /*
     * The Block Lanczos driver reports numerical basis_count after
     * the Q1/Beta1 transition. This provides an observable runtime
     * stage marker without adding a test hook or public API.
     */
    esolver->options[
        LIS_EOPTIONS_OUTPUT] = 1;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    /*
     * lis_gesolve() stores the eigensolver execution status in
     * esolver->retcode but returns LIS_SUCCESS after execution.
     * lis_esolve() therefore also returns LIS_SUCCESS here.
     *
     * The stage boundary is identified by:
     *
     *     public API return  = LIS_SUCCESS
     *     esolver->retcode   = LIS_ERR_NOT_IMPLEMENTED
     *
     * together with the numerical basis_count printed by the
     * Block Lanczos driver after Q1/Beta1 construction.
     */
    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos runtime two-block probe ss=%d: "
            "API return %d, expected LIS_SUCCESS\n",
            (int)ss,
            (int)api);

        err = LIS_FAILS;
        goto cleanup;
    }

    if( esolver->retcode!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos runtime two-block probe ss=%d: "
            "status %d, expected LIS_SUCCESS\n",
            (int)ss,
            (int)esolver->retcode);

        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Block Lanczos runtime two-block probe ss=%d\n",
        (int)ss);

    err = LIS_SUCCESS;

cleanup:

    if( esolver )
    {
        lis_esolver_destroy(
            esolver);
    }

    if( x )
    {
        lis_vector_destroy(
            x);
    }

    if( A )
    {
        lis_matrix_destroy(
            A);
    }

    return err;
}


static LIS_INT
test_block_lanczos_two_block_runtime(void)
{
    LIS_INT err;

    err =
        test_block_lanczos_two_block_runtime_case(
            4);

    if( err )
    {
        return err;
    }

    err =
        test_block_lanczos_two_block_runtime_case(
            3);

    if( err )
    {
        return err;
    }

    printf(
        "PASS: Block Lanczos full and partial "
        "Q1/Beta1 runtime paths\n");

    return LIS_SUCCESS;
}



/*
 * Exercise the general Block Lanczos recurrence beyond the
 * validated Q0 -> Q1 bootstrap.
 *
 * n=6, ss=6 gives the rank path 2 -> 2 -> 2.
 * n=5, ss=5 gives the rank path 2 -> 2 -> 1.
 */
static LIS_INT
test_block_lanczos_multiblock_runtime_case(
    LIS_INT n,
    LIS_INT ss)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT i;
    LIS_INT is,ie;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;

    A = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_set_size(
        A,
        0,
        n);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if( err )
    {
        goto cleanup;
    }

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            (LIS_SCALAR)(i+1),
            A);

        if( err )
        {
            goto cleanup;
        }
    }

    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_assemble(
        A);

    if( err )
    {
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A,
        &x);

    if( err )
    {
        goto cleanup;
    }

    err = lis_vector_set_all(
        1.0,
        x);

    if( err )
    {
        goto cleanup;
    }

    err = lis_esolver_create(
        &esolver);

    if( err )
    {
        goto cleanup;
    }

    if( n==6 && ss==6 )
    {
        err = lis_esolver_set_option(
            "-e li -ss 6 -eblock 2 -rval true",
            esolver);
    }
    else if( n==5 && ss==5 )
    {
        err = lis_esolver_set_option(
            "-e li -ss 5 -eblock 2 -rval true",
            esolver);
    }
    else
    {
        err = LIS_ERR_ILL_ARG;
    }

    if( err )
    {
        goto cleanup;
    }

    esolver->options[
        LIS_EOPTIONS_OUTPUT] = 1;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos multi-block runtime probe "
            "n=%d ss=%d: API return %d, "
            "expected LIS_SUCCESS\n",
            (int)n,
            (int)ss,
            (int)api);

        err = LIS_FAILS;
        goto cleanup;
    }

    if( esolver->retcode!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos multi-block runtime probe "
            "n=%d ss=%d: status %d, "
            "expected LIS_SUCCESS\n",
            (int)n,
            (int)ss,
            (int)esolver->retcode);

        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Block Lanczos multi-block runtime probe "
        "n=%d ss=%d\n",
        (int)n,
        (int)ss);

    err = LIS_SUCCESS;

cleanup:

    if( esolver )
    {
        lis_esolver_destroy(
            esolver);
    }

    if( x )
    {
        lis_vector_destroy(
            x);
    }

    if( A )
    {
        lis_matrix_destroy(
            A);
    }

    return err;
}

static LIS_INT
test_block_lanczos_multiblock_runtime(void)
{
    LIS_INT err;

    err =
        test_block_lanczos_multiblock_runtime_case(
            6,
            6);

    if( err )
    {
        return err;
    }

    err =
        test_block_lanczos_multiblock_runtime_case(
            5,
            5);

    if( err )
    {
        return err;
    }

    printf(
        "PASS: Block Lanczos Q1-to-Q2 full and "
        "rank-reduced runtime paths\n");

    return LIS_SUCCESS;
}



/*
 * Validate the active projected eigensystem through the public
 * Block Lanczos path.
 *
 * Full-basis diagonal cases must recover the exact spectrum.
 * The identity case closes after Q0 and verifies that the
 * projected solve uses basis_count rather than ss.
 */
static LIS_INT
test_block_lanczos_projected_spectrum_case(
    LIS_INT n,
    LIS_INT identity,
    LIS_INT expected_dim)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT i,j;
    LIS_INT is,ie;

    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;

    LIS_SCALAR evalue;
    LIS_REAL values[6];
    LIS_REAL tmp;
    LIS_REAL expected;
    LIS_REAL diff;
    LIS_REAL tol;

    A = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    tol = (LIS_REAL)1.0e-8;

    if( n<1 || n>6 ||
        expected_dim<1 ||
        expected_dim>n )
    {
        return LIS_ERR_ILL_ARG;
    }

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_set_size(
        A,
        0,
        n);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if( err )
    {
        goto cleanup;
    }

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            identity
                ? (LIS_SCALAR)1.0
                : (LIS_SCALAR)(i+1),
            A);

        if( err )
        {
            goto cleanup;
        }
    }

    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);

    if( err )
    {
        goto cleanup;
    }

    err = lis_matrix_assemble(
        A);

    if( err )
    {
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A,
        &x);

    if( err )
    {
        goto cleanup;
    }

    err = lis_vector_set_all(
        1.0,
        x);

    if( err )
    {
        goto cleanup;
    }

    err = lis_esolver_create(
        &esolver);

    if( err )
    {
        goto cleanup;
    }

    if( n==4 )
    {
        err = lis_esolver_set_option(
            "-e li -ss 4 -eblock 2 -rval true",
            esolver);
    }
    else if( n==5 )
    {
        err = lis_esolver_set_option(
            "-e li -ss 5 -eblock 2 -rval true",
            esolver);
    }
    else if( n==6 )
    {
        err = lis_esolver_set_option(
            "-e li -ss 6 -eblock 2 -rval true",
            esolver);
    }
    else
    {
        err = LIS_ERR_ILL_ARG;
    }

    if( err )
    {
        goto cleanup;
    }

    esolver->options[
        LIS_EOPTIONS_OUTPUT] = 1;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos projected spectrum n=%d "
            "identity=%d: API return %d\n",
            (int)n,
            (int)identity,
            (int)api);

        err = LIS_FAILS;
        goto cleanup;
    }

    if( esolver->retcode!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos projected spectrum n=%d "
            "identity=%d: status %d\n",
            (int)n,
            (int)identity,
            (int)esolver->retcode);

        err = LIS_FAILS;
        goto cleanup;
    }

    for(i=0;i<expected_dim;i++)
    {
#ifdef _COMPLEX
        if( fabs(
                (double)cimag(
                    esolver->evalue[i])
            ) > (double)tol )
        {
            printf(
                "Block Lanczos projected spectrum n=%d: "
                "non-real Ritz value at %d\n",
                (int)n,
                (int)i);

            err = LIS_FAILS;
            goto cleanup;
        }

        values[i] =
            (LIS_REAL)creal(
                esolver->evalue[i]);
#else
        values[i] =
            (LIS_REAL)esolver->evalue[i];
#endif
    }

    for(i=0;i<expected_dim-1;i++)
    {
        for(j=i+1;j<expected_dim;j++)
        {
            if( values[j]<values[i] )
            {
                tmp = values[i];
                values[i] = values[j];
                values[j] = tmp;
            }
        }
    }

    for(i=0;i<expected_dim;i++)
    {
        expected =
            identity
                ? (LIS_REAL)1.0
                : (LIS_REAL)(i+1);

        diff = (LIS_REAL)fabs(
            (double)(values[i]-expected));

        if( diff>tol )
        {
            printf(
                "Block Lanczos projected spectrum n=%d "
                "identity=%d index=%d got=%e expected=%e "
                "diff=%e\n",
                (int)n,
                (int)identity,
                (int)i,
                (double)values[i],
                (double)expected,
                (double)diff);

            err = LIS_FAILS;
            goto cleanup;
        }
    }

    printf(
        "PASS: Block Lanczos projected spectrum "
        "n=%d identity=%d dim=%d\n",
        (int)n,
        (int)identity,
        (int)expected_dim);

    err = LIS_SUCCESS;

cleanup:

    if( esolver )
    {
        lis_esolver_destroy(
            esolver);
    }

    if( x )
    {
        lis_vector_destroy(
            x);
    }

    if( A )
    {
        lis_matrix_destroy(
            A);
    }

    return err;
}


static LIS_INT
test_block_lanczos_projected_spectrum(void)
{
    LIS_INT err;

    err =
        test_block_lanczos_projected_spectrum_case(
            4,
            0,
            4);

    if( err )
    {
        return err;
    }

    err =
        test_block_lanczos_projected_spectrum_case(
            5,
            0,
            5);

    if( err )
    {
        return err;
    }

    err =
        test_block_lanczos_projected_spectrum_case(
            6,
            0,
            6);

    if( err )
    {
        return err;
    }

    err =
        test_block_lanczos_projected_spectrum_case(
            6,
            1,
            2);

    if( err )
    {
        return err;
    }

    printf(
        "PASS: Block Lanczos projected solve "
        "full-basis and happy-breakdown spectra\n");

    return LIS_SUCCESS;
}



static LIS_INT
test_block_lanczos_ritz_evidence_case(
    LIS_INT n,
    LIS_INT identity,
    LIS_INT expected_dim)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT i,j;
    LIS_INT is,ie;
    LIS_INT local_count,global_count;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_VECTOR evals;
    LIS_VECTOR residuals;
    LIS_VECTOR Ax;
    LIS_VECTOR residual;
    LIS_VECTOR ritz[6];
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_SCALAR lambda;
    LIS_SCALAR dot;
    LIS_REAL norm;
    LIS_REAL defect;
    LIS_REAL overlap;
    LIS_REAL pair_resid;
    LIS_REAL stored_resid;
    LIS_REAL scalar_resid;
    LIS_REAL resid_delta;
    LIS_REAL max_resid_delta;
    LIS_REAL max_norm_defect;
    LIS_REAL max_overlap;
    LIS_REAL max_resid;
    LIS_REAL orth_tol;
    LIS_REAL resid_tol;
    char option[128];

    A = NULL;
    x = NULL;
    evals = NULL;
    residuals = NULL;
    Ax = NULL;
    residual = NULL;
    esolver = NULL;
    evalue = 0.0;

    for(i=0;i<6;i++)
    {
        ritz[i] = NULL;
    }

    max_norm_defect = 0.0;
    max_overlap = 0.0;
    max_resid = 0.0;
    max_resid_delta = 0.0;
    orth_tol = (LIS_REAL)1.0e-8;
    resid_tol = (LIS_REAL)1.0e-8;

    if( n<1 ||
        n>6 ||
        expected_dim<1 ||
        expected_dim>n )
    {
        return LIS_ERR_ILL_ARG;
    }

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);
    if( err ) goto cleanup;

    err = lis_matrix_set_size(
        A,
        0,
        n);
    if( err ) goto cleanup;

    err = lis_matrix_get_range(
        A,
        &is,
        &ie);
    if( err ) goto cleanup;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            identity
                ? (LIS_SCALAR)1.0
                : (LIS_SCALAR)(i+1),
            A);
        if( err ) goto cleanup;
    }

    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);
    if( err ) goto cleanup;

    err = lis_matrix_assemble(A);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err ) goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x);
    if( err ) goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err ) goto cleanup;

    snprintf(
        option,
        sizeof(option),
        "-e li -ss %d -eblock 2 -rval false",
        (int)n);

    err = lis_esolver_set_option(
        option,
        esolver);
    if( err ) goto cleanup;

    esolver->options[
        LIS_EOPTIONS_OUTPUT] = 0;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos Ritz evidence n=%d "
            "identity=%d: API return %d\n",
            (int)n,
            (int)identity,
            (int)api);

        err = LIS_FAILS;
        goto cleanup;
    }

    if( esolver->retcode!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos Ritz evidence n=%d "
            "identity=%d: status %d\n",
            (int)n,
            (int)identity,
            (int)esolver->retcode);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_vector_create(
        LIS_COMM_WORLD,
        &evals);
    if( err ) goto cleanup;

    err = lis_esolver_get_evalues(
        esolver,
        evals);
    if( err ) goto cleanup;

    err = lis_vector_get_size(
        evals,
        &local_count,
        &global_count);
    if( err ) goto cleanup;

    if( global_count!=expected_dim )
    {
        printf(
            "Block Lanczos Ritz evidence n=%d "
            "identity=%d: result count %d, expected %d\n",
            (int)n,
            (int)identity,
            (int)global_count,
            (int)expected_dim);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_vector_create(
        LIS_COMM_WORLD,
        &residuals);
    if( err ) goto cleanup;

    err = lis_esolver_get_residualnorms(
        esolver,
        residuals);
    if( err ) goto cleanup;

    err = lis_vector_get_size(
        residuals,
        &local_count,
        &global_count);
    if( err ) goto cleanup;

    if( global_count!=expected_dim )
    {
        printf(
            "Block Lanczos residual metadata n=%d "
            "identity=%d: residual count %d, expected %d\n",
            (int)n,
            (int)identity,
            (int)global_count,
            (int)expected_dim);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A,
        &Ax);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &residual);
    if( err ) goto cleanup;

    for(i=0;i<expected_dim;i++)
    {
        err = lis_vector_duplicate(
            A,
            &ritz[i]);
        if( err ) goto cleanup;

        err = lis_esolver_get_specific_evalue(
            esolver,
            i,
            &lambda);
        if( err ) goto cleanup;

        err = lis_esolver_get_specific_evector(
            esolver,
            i,
            ritz[i]);
        if( err ) goto cleanup;

        err = lis_vector_nrm2(
            ritz[i],
            &norm);
        if( err ) goto cleanup;

        defect = (LIS_REAL)fabs(
            (double)(norm-(LIS_REAL)1.0));

        if( defect>max_norm_defect )
            max_norm_defect = defect;

        if( defect>orth_tol )
        {
            printf(
                "Block Lanczos Ritz evidence n=%d "
                "identity=%d mode=%d: norm=%e defect=%e\n",
                (int)n,
                (int)identity,
                (int)i,
                (double)norm,
                (double)defect);

            err = LIS_FAILS;
            goto cleanup;
        }

        err = lis_esolver_standard_residual(
            A,
            ritz[i],
            lambda,
            Ax,
            residual,
            &pair_resid);
        if( err ) goto cleanup;

        if( pair_resid>max_resid )
            max_resid = pair_resid;

        err = lis_esolver_get_specific_residualnorm(
            esolver,
            i,
            &stored_resid);
        if( err ) goto cleanup;

        resid_delta = (LIS_REAL)fabs(
            (double)(stored_resid-pair_resid));

        if( resid_delta>max_resid_delta )
            max_resid_delta = resid_delta;

        if( resid_delta>(LIS_REAL)1.0e-14 )
        {
            printf(
                "Block Lanczos residual metadata n=%d "
                "identity=%d mode=%d: stored=%e true=%e delta=%e\n",
                (int)n,
                (int)identity,
                (int)i,
                (double)stored_resid,
                (double)pair_resid,
                (double)resid_delta);

            err = LIS_FAILS;
            goto cleanup;
        }

        if( i==0 )
        {
            err = lis_esolver_get_residualnorm(
                esolver,
                &scalar_resid);
            if( err ) goto cleanup;

            resid_delta = (LIS_REAL)fabs(
                (double)(scalar_resid-pair_resid));

            if( resid_delta>max_resid_delta )
                max_resid_delta = resid_delta;

            if( resid_delta>(LIS_REAL)1.0e-14 )
            {
                printf(
                    "Block Lanczos residual metadata n=%d "
                    "identity=%d mode=0 scalar: stored=%e true=%e delta=%e\n",
                    (int)n,
                    (int)identity,
                    (double)scalar_resid,
                    (double)pair_resid,
                    (double)resid_delta);

                err = LIS_FAILS;
                goto cleanup;
            }
        }

        if( pair_resid>resid_tol )
        {
            printf(
                "Block Lanczos Ritz evidence n=%d "
                "identity=%d mode=%d: residual=%e\n",
                (int)n,
                (int)identity,
                (int)i,
                (double)pair_resid);

            err = LIS_FAILS;
            goto cleanup;
        }
    }

    for(i=0;i<expected_dim;i++)
    {
        for(j=i+1;j<expected_dim;j++)
        {
            err = lis_vector_dot(
                ritz[i],
                ritz[j],
                &dot);
            if( err ) goto cleanup;

#ifdef _COMPLEX
            overlap = (LIS_REAL)cabs(dot);
#else
            overlap = (LIS_REAL)fabs((double)dot);
#endif

            if( overlap>max_overlap )
                max_overlap = overlap;

            if( overlap>orth_tol )
            {
                printf(
                    "Block Lanczos Ritz evidence n=%d "
                    "identity=%d pair=(%d,%d): overlap=%e\n",
                    (int)n,
                    (int)identity,
                    (int)i,
                    (int)j,
                    (double)overlap);

                err = LIS_FAILS;
                goto cleanup;
            }
        }
    }

    printf(
        "PASS: Block Lanczos residual metadata "
        "n=%d identity=%d count=%d max_delta=%e\n",
        (int)n,
        (int)identity,
        (int)global_count,
        (double)max_resid_delta);

    printf(
        "PASS: Block Lanczos Ritz evidence "
        "n=%d identity=%d count=%d "
        "norm_defect=%e overlap=%e residual=%e\n",
        (int)n,
        (int)identity,
        (int)global_count,
        (double)max_norm_defect,
        (double)max_overlap,
        (double)max_resid);

    err = LIS_SUCCESS;

cleanup:
    for(i=0;i<6;i++)
    {
        if( ritz[i] )
            lis_vector_destroy(ritz[i]);
    }

    if( residual )
        lis_vector_destroy(residual);

    if( Ax )
        lis_vector_destroy(Ax);

    if( residuals )
        lis_vector_destroy(residuals);

    if( evals )
        lis_vector_destroy(evals);

    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    if( A )
        lis_matrix_destroy(A);

    return err;
}



static LIS_INT
test_block_lanczos_residual_status_case(
    LIS_REAL tol,
    LIS_INT expected_status,
    const char *label)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_INT i;
    LIS_INT is,ie;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_REAL residual;
    char option[160];

    A = NULL;
    x = NULL;
    esolver = NULL;
    evalue = 0.0;
    residual = 0.0;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);
    if( err ) goto cleanup;

    err = lis_matrix_set_size(
        A,
        0,
        4);
    if( err ) goto cleanup;

    err = lis_matrix_get_range(
        A,
        &is,
        &ie);
    if( err ) goto cleanup;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            (LIS_SCALAR)(i+1),
            A);
        if( err ) goto cleanup;
    }

    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);
    if( err ) goto cleanup;

    err = lis_matrix_assemble(
        A);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err ) goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x);
    if( err ) goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err ) goto cleanup;

    snprintf(
        option,
        sizeof(option),
        "-e li -ss 4 -eblock 2 -rval false -etol %.17e",
        (double)tol);

    err = lis_esolver_set_option(
        option,
        esolver);
    if( err ) goto cleanup;

    /*
     * Historical public API semantics:
     * execution completion returns LIS_SUCCESS while
     * convergence status is reported through esolver->retcode.
     */
    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos status boundary %s: "
            "API return %d, expected LIS_SUCCESS (%d)\n",
            label,
            (int)api,
            (int)LIS_SUCCESS);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver,
        &status);
    if( err ) goto cleanup;

    if( status!=expected_status )
    {
        printf(
            "Block Lanczos status boundary %s: "
            "status %d, expected %d\n",
            label,
            (int)status,
            (int)expected_status);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_residualnorm(
        esolver,
        &residual);
    if( err ) goto cleanup;

    if( expected_status==LIS_SUCCESS )
    {
        if( !(residual<=tol) )
        {
            printf(
                "Block Lanczos status boundary %s: "
                "residual %.17e, expected <= %.17e\n",
                label,
                (double)residual,
                (double)tol);

            err = LIS_FAILS;
            goto cleanup;
        }
    }
    else if( expected_status==LIS_MAXITER )
    {
        if( !(residual>tol) )
        {
            printf(
                "Block Lanczos status boundary %s: "
                "residual %.17e, expected > %.17e\n",
                label,
                (double)residual,
                (double)tol);

            err = LIS_FAILS;
            goto cleanup;
        }
    }
    else
    {
        printf(
            "Block Lanczos status boundary %s: "
            "unsupported expected status %d\n",
            label,
            (int)expected_status);

        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Block Lanczos status boundary %s "
        "status=%d residual=%.6e tol=%.6e\n",
        label,
        (int)status,
        (double)residual,
        (double)tol);

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
test_block_lanczos_residual_status_boundary(void)
{
    LIS_INT err;

    err =
        test_block_lanczos_residual_status_case(
            (LIS_REAL)1.0e-12,
            LIS_SUCCESS,
            "loose");
    if( err ) return err;

    err =
        test_block_lanczos_residual_status_case(
            (LIS_REAL)1.0e-18,
            LIS_MAXITER,
            "tight");
    if( err ) return err;

    printf(
        "PASS: Block Lanczos certified residual status boundary\n");

    return LIS_SUCCESS;
}



static void
test_block_lanczos_sort_real(
    LIS_REAL value[],
    LIS_INT n)
{
    LIS_INT i,j;
    LIS_REAL key;

    for(i=1;i<n;i++)
    {
        key = value[i];
        j = i-1;

        while( j>=0 && value[j]>key )
        {
            value[j+1] = value[j];
            j--;
        }

        value[j+1] = key;
    }
}


static LIS_INT
test_block_lanczos_spectrum_case(
    const char *label,
    LIS_INT block_size,
    const LIS_REAL diagonal[6],
    LIS_REAL eig_tol)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_INT i,j;
    LIS_INT is,ie;
    LIS_INT local_count,global_count;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_VECTOR evals;
    LIS_VECTOR Ax;
    LIS_VECTOR residual;
    LIS_VECTOR ritz[6];
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_SCALAR lambda;
    LIS_SCALAR dot;
    LIS_REAL expected[6];
    LIS_REAL actual[6];
    LIS_REAL norm;
    LIS_REAL norm_defect;
    LIS_REAL overlap;
    LIS_REAL pair_resid;
    LIS_REAL stored_resid;
    LIS_REAL delta;
    LIS_REAL max_eig_error;
    LIS_REAL max_norm_defect;
    LIS_REAL max_overlap;
    LIS_REAL max_resid;
    LIS_REAL max_metadata_delta;
    char option[160];

    A = NULL;
    x = NULL;
    evals = NULL;
    Ax = NULL;
    residual = NULL;
    esolver = NULL;
    evalue = 0.0;

    for(i=0;i<6;i++)
    {
        ritz[i] = NULL;
        expected[i] = diagonal[i];
        actual[i] = 0.0;
    }

    max_eig_error = 0.0;
    max_norm_defect = 0.0;
    max_overlap = 0.0;
    max_resid = 0.0;
    max_metadata_delta = 0.0;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);
    if( err ) goto cleanup;

    err = lis_matrix_set_size(
        A,
        0,
        6);
    if( err ) goto cleanup;

    err = lis_matrix_get_range(
        A,
        &is,
        &ie);
    if( err ) goto cleanup;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            (LIS_SCALAR)diagonal[i],
            A);
        if( err ) goto cleanup;
    }

    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);
    if( err ) goto cleanup;

    err = lis_matrix_assemble(
        A);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err ) goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x);
    if( err ) goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err ) goto cleanup;

    snprintf(
        option,
        sizeof(option),
        "-e li -ss 6 -eblock %d -rval false -etol 1.0e-12",
        (int)block_size);

    err = lis_esolver_set_option(
        option,
        esolver);
    if( err ) goto cleanup;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos %s: API return %d, expected LIS_SUCCESS\n",
            label,
            (int)api);

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
            "Block Lanczos %s: status %d, expected LIS_SUCCESS\n",
            label,
            (int)status);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_vector_create(
        LIS_COMM_WORLD,
        &evals);
    if( err ) goto cleanup;

    err = lis_esolver_get_evalues(
        esolver,
        evals);
    if( err ) goto cleanup;

    err = lis_vector_get_size(
        evals,
        &local_count,
        &global_count);
    if( err ) goto cleanup;

    if( global_count!=6 )
    {
        printf(
            "Block Lanczos %s: result count %d, expected 6\n",
            label,
            (int)global_count);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A,
        &Ax);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &residual);
    if( err ) goto cleanup;

    for(i=0;i<6;i++)
    {
        err = lis_vector_duplicate(
            A,
            &ritz[i]);
        if( err ) goto cleanup;

        err = lis_esolver_get_specific_evalue(
            esolver,
            i,
            &lambda);
        if( err ) goto cleanup;

#ifdef _COMPLEX
        actual[i] = (LIS_REAL)creal(lambda);
#else
        actual[i] = (LIS_REAL)lambda;
#endif

        err = lis_esolver_get_specific_evector(
            esolver,
            i,
            ritz[i]);
        if( err ) goto cleanup;

        err = lis_vector_nrm2(
            ritz[i],
            &norm);
        if( err ) goto cleanup;

        norm_defect = (LIS_REAL)fabs(
            (double)(norm-1.0));

        if( norm_defect>max_norm_defect )
            max_norm_defect = norm_defect;

        err = lis_esolver_standard_residual(
            A,
            ritz[i],
            lambda,
            Ax,
            residual,
            &pair_resid);
        if( err ) goto cleanup;

        if( pair_resid>max_resid )
            max_resid = pair_resid;

        err = lis_esolver_get_specific_residualnorm(
            esolver,
            i,
            &stored_resid);
        if( err ) goto cleanup;

        delta = (LIS_REAL)fabs(
            (double)(stored_resid-pair_resid));

        if( delta>max_metadata_delta )
            max_metadata_delta = delta;
    }

    for(i=0;i<6;i++)
    {
        for(j=i+1;j<6;j++)
        {
            err = lis_vector_dot(
                ritz[i],
                ritz[j],
                &dot);
            if( err ) goto cleanup;

#ifdef _COMPLEX
            overlap = (LIS_REAL)cabs(dot);
#else
            overlap = (LIS_REAL)fabs((double)dot);
#endif

            if( overlap>max_overlap )
                max_overlap = overlap;
        }
    }

    test_block_lanczos_sort_real(
        expected,
        6);

    test_block_lanczos_sort_real(
        actual,
        6);

    for(i=0;i<6;i++)
    {
        delta = (LIS_REAL)fabs(
            (double)(actual[i]-expected[i]));

        if( delta>max_eig_error )
            max_eig_error = delta;
    }

    if( max_eig_error>eig_tol ||
        max_norm_defect>(LIS_REAL)1.0e-10 ||
        max_overlap>(LIS_REAL)1.0e-10 ||
        max_resid>(LIS_REAL)1.0e-10 ||
        max_metadata_delta>(LIS_REAL)1.0e-14 )
    {
        printf(
            "Block Lanczos %s failed: "
            "eigerr=%e norm=%e overlap=%e residual=%e metadata=%e\n",
            label,
            (double)max_eig_error,
            (double)max_norm_defect,
            (double)max_overlap,
            (double)max_resid,
            (double)max_metadata_delta);

        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Block Lanczos %s eblock=%d count=6 "
        "eigerr=%e norm=%e overlap=%e residual=%e metadata=%e\n",
        label,
        (int)block_size,
        (double)max_eig_error,
        (double)max_norm_defect,
        (double)max_overlap,
        (double)max_resid,
        (double)max_metadata_delta);

    err = LIS_SUCCESS;

cleanup:
    for(i=0;i<6;i++)
    {
        if( ritz[i] )
            lis_vector_destroy(ritz[i]);
    }

    if( residual )
        lis_vector_destroy(residual);

    if( Ax )
        lis_vector_destroy(Ax);

    if( evals )
        lis_vector_destroy(evals);

    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_block_lanczos_clustered_repeated_spectra(void)
{
    LIS_INT err;

    static const LIS_REAL repeated_pairs[6] = {
        1.0,1.0,
        2.0,2.0,
        3.0,3.0
    };

    static const LIS_REAL clustered_pairs[6] = {
        1.0,1.00000001,
        2.0,2.00000001,
        3.0,3.00000001
    };

    static const LIS_REAL repeated_triples[6] = {
        1.0,1.0,1.0,
        2.0,2.0,2.0
    };

    err = test_block_lanczos_spectrum_case(
        "repeated-pairs",
        2,
        repeated_pairs,
        (LIS_REAL)1.0e-10);
    if( err ) return err;

    err = test_block_lanczos_spectrum_case(
        "clustered-pairs",
        2,
        clustered_pairs,
        (LIS_REAL)1.0e-10);
    if( err ) return err;

    err = test_block_lanczos_spectrum_case(
        "repeated-triples",
        3,
        repeated_triples,
        (LIS_REAL)1.0e-10);
    if( err ) return err;

    printf(
        "PASS: Block Lanczos clustered/repeated spectra\n");

    return LIS_SUCCESS;
}



static LIS_INT
test_block_lanczos_build_diagonal(
    LIS_INT n,
    LIS_MATRIX *Aout)
{
    LIS_INT err;
    LIS_INT i;
    LIS_INT is,ie;
    LIS_MATRIX A;

    A = NULL;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);
    if( err ) return err;

    err = lis_matrix_set_size(
        A,
        0,
        n);
    if( err ) goto cleanup;

    err = lis_matrix_get_range(
        A,
        &is,
        &ie);
    if( err ) goto cleanup;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            (LIS_SCALAR)(i+1),
            A);
        if( err ) goto cleanup;
    }

    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);
    if( err ) goto cleanup;

    err = lis_matrix_assemble(
        A);
    if( err ) goto cleanup;

    *Aout = A;
    return LIS_SUCCESS;

cleanup:
    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_block_lanczos_scalar_compatibility(void)
{
    LIS_INT err;
    LIS_INT api_default,api_explicit;
    LIS_INT status_default,status_explicit;
    LIS_INT i;
    LIS_INT local_count,global_count;
    LIS_MATRIX A;
    LIS_VECTOR x_default;
    LIS_VECTOR x_explicit;
    LIS_VECTOR evals_default;
    LIS_VECTOR evals_explicit;
    LIS_ESOLVER esolver_default;
    LIS_ESOLVER esolver_explicit;
    LIS_SCALAR evalue_default;
    LIS_SCALAR evalue_explicit;
    LIS_SCALAR lambda_default;
    LIS_SCALAR lambda_explicit;
    LIS_REAL delta;
    LIS_REAL max_delta;

    A = NULL;
    x_default = NULL;
    x_explicit = NULL;
    evals_default = NULL;
    evals_explicit = NULL;
    esolver_default = NULL;
    esolver_explicit = NULL;
    evalue_default = 0.0;
    evalue_explicit = 0.0;
    max_delta = 0.0;

    err = test_block_lanczos_build_diagonal(
        5,
        &A);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x_default);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x_explicit);
    if( err ) goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x_default);
    if( err ) goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x_explicit);
    if( err ) goto cleanup;

    err = lis_esolver_create(
        &esolver_default);
    if( err ) goto cleanup;

    err = lis_esolver_create(
        &esolver_explicit);
    if( err ) goto cleanup;

    if( esolver_default->block_size!=1 )
    {
        printf(
            "Block Lanczos scalar compatibility: "
            "default block size %d, expected 1\n",
            (int)esolver_default->block_size);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_set_option(
        "-e li -ss 4 -rval true",
        esolver_default);
    if( err ) goto cleanup;

    err = lis_esolver_set_option(
        "-e li -ss 4 -eblock 1 -rval true",
        esolver_explicit);
    if( err ) goto cleanup;

    if( esolver_explicit->block_size!=1 )
    {
        printf(
            "Block Lanczos scalar compatibility: "
            "explicit block size %d, expected 1\n",
            (int)esolver_explicit->block_size);

        err = LIS_FAILS;
        goto cleanup;
    }

    api_default = lis_esolve(
        A,
        x_default,
        &evalue_default,
        esolver_default);

    api_explicit = lis_esolve(
        A,
        x_explicit,
        &evalue_explicit,
        esolver_explicit);

    if( api_default!=LIS_SUCCESS ||
        api_explicit!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos scalar compatibility: "
            "API default=%d explicit=%d\n",
            (int)api_default,
            (int)api_explicit);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_status(
        esolver_default,
        &status_default);
    if( err ) goto cleanup;

    err = lis_esolver_get_status(
        esolver_explicit,
        &status_explicit);
    if( err ) goto cleanup;

    if( status_default!=status_explicit )
    {
        printf(
            "Block Lanczos scalar compatibility: "
            "status default=%d explicit=%d\n",
            (int)status_default,
            (int)status_explicit);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_vector_create(
        LIS_COMM_WORLD,
        &evals_default);
    if( err ) goto cleanup;

    err = lis_vector_create(
        LIS_COMM_WORLD,
        &evals_explicit);
    if( err ) goto cleanup;

    err = lis_esolver_get_evalues(
        esolver_default,
        evals_default);
    if( err ) goto cleanup;

    err = lis_vector_get_size(
        evals_default,
        &local_count,
        &global_count);
    if( err ) goto cleanup;

    if( global_count!=4 )
    {
        printf(
            "Block Lanczos scalar compatibility: "
            "default result count %d, expected 4\n",
            (int)global_count);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_esolver_get_evalues(
        esolver_explicit,
        evals_explicit);
    if( err ) goto cleanup;

    err = lis_vector_get_size(
        evals_explicit,
        &local_count,
        &global_count);
    if( err ) goto cleanup;

    if( global_count!=4 )
    {
        printf(
            "Block Lanczos scalar compatibility: "
            "explicit result count %d, expected 4\n",
            (int)global_count);

        err = LIS_FAILS;
        goto cleanup;
    }

    for(i=0;i<4;i++)
    {
        err = lis_esolver_get_specific_evalue(
            esolver_default,
            i,
            &lambda_default);
        if( err ) goto cleanup;

        err = lis_esolver_get_specific_evalue(
            esolver_explicit,
            i,
            &lambda_explicit);
        if( err ) goto cleanup;

#ifdef _COMPLEX
        delta = (LIS_REAL)cabs(
            lambda_default-lambda_explicit);
#else
        delta = (LIS_REAL)fabs(
            (double)(lambda_default-lambda_explicit));
#endif

        if( delta>max_delta )
            max_delta = delta;
    }

    if( max_delta>(LIS_REAL)1.0e-12 )
    {
        printf(
            "Block Lanczos scalar compatibility: "
            "default/explicit eigenvalue delta %e\n",
            (double)max_delta);

        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Block Lanczos scalar compatibility "
        "default-vs-eblock1 count=4 max_delta=%e status=%d\n",
        (double)max_delta,
        (int)status_default);

    err = LIS_SUCCESS;

cleanup:
    if( evals_explicit )
        lis_vector_destroy(evals_explicit);

    if( evals_default )
        lis_vector_destroy(evals_default);

    if( esolver_explicit )
        lis_esolver_destroy(esolver_explicit);

    if( esolver_default )
        lis_esolver_destroy(esolver_default);

    if( x_explicit )
        lis_vector_destroy(x_explicit);

    if( x_default )
        lis_vector_destroy(x_default);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_block_lanczos_raw_ritz_compatibility(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_INT i;
    LIS_INT local_count,global_count;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_VECTOR evals;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_SCALAR lambda;
    LIS_REAL actual[4];
    LIS_REAL expected[4];
    LIS_REAL delta;
    LIS_REAL max_error;

    A = NULL;
    x = NULL;
    evals = NULL;
    esolver = NULL;
    evalue = 0.0;
    max_error = 0.0;

    for(i=0;i<4;i++)
    {
        actual[i] = 0.0;
        expected[i] = (LIS_REAL)(i+1);
    }

    err = test_block_lanczos_build_diagonal(
        4,
        &A);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err ) goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x);
    if( err ) goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err ) goto cleanup;

    err = lis_esolver_set_option(
        "-e li -ss 4 -eblock 2 -rval true -etol 1.0e-12",
        esolver);
    if( err ) goto cleanup;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos raw Ritz compatibility: "
            "API return %d\n",
            (int)api);

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
            "Block Lanczos raw Ritz compatibility: "
            "status %d, expected LIS_SUCCESS\n",
            (int)status);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_vector_create(
        LIS_COMM_WORLD,
        &evals);
    if( err ) goto cleanup;

    err = lis_esolver_get_evalues(
        esolver,
        evals);
    if( err ) goto cleanup;

    err = lis_vector_get_size(
        evals,
        &local_count,
        &global_count);
    if( err ) goto cleanup;

    if( global_count!=4 )
    {
        printf(
            "Block Lanczos raw Ritz compatibility: "
            "result count %d, expected 4\n",
            (int)global_count);

        err = LIS_FAILS;
        goto cleanup;
    }

    /*
     * Raw Ritz mode must finish before physical Ritz reconstruction.
     * lis_gesolve allocates and zeroes the pointer array, so entry zero
     * remains NULL when the block driver returns at the raw boundary.
     */
    if( esolver->evector==NULL ||
        esolver->evector[0]!=NULL )
    {
        printf(
            "Block Lanczos raw Ritz compatibility: "
            "unexpected physical Ritz vector reconstruction\n");

        err = LIS_FAILS;
        goto cleanup;
    }

    for(i=0;i<4;i++)
    {
        err = lis_esolver_get_specific_evalue(
            esolver,
            i,
            &lambda);
        if( err ) goto cleanup;

#ifdef _COMPLEX
        actual[i] = (LIS_REAL)creal(lambda);
#else
        actual[i] = (LIS_REAL)lambda;
#endif
    }

    test_block_lanczos_sort_real(
        actual,
        4);

    test_block_lanczos_sort_real(
        expected,
        4);

    for(i=0;i<4;i++)
    {
        delta = (LIS_REAL)fabs(
            (double)(actual[i]-expected[i]));

        if( delta>max_error )
            max_error = delta;
    }

    if( max_error>(LIS_REAL)1.0e-10 )
    {
        printf(
            "Block Lanczos raw Ritz compatibility: "
            "spectrum error %e\n",
            (double)max_error);

        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Block Lanczos raw Ritz compatibility "
        "count=4 max_error=%e evector0=NULL\n",
        (double)max_error);

    err = LIS_SUCCESS;

cleanup:
    if( evals )
        lis_vector_destroy(evals);

    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_block_lanczos_returned_x_compatibility(void)
{
    LIS_INT err;
    LIS_INT api;
    LIS_INT status;
    LIS_MATRIX A;
    LIS_VECTOR x;
    LIS_VECTOR evector0;
    LIS_VECTOR diff;
    LIS_ESOLVER esolver;
    LIS_SCALAR evalue;
    LIS_REAL defect;

    A = NULL;
    x = NULL;
    evector0 = NULL;
    diff = NULL;
    esolver = NULL;
    evalue = 0.0;
    defect = 0.0;

    err = test_block_lanczos_build_diagonal(
        4,
        &A);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &x);
    if( err ) goto cleanup;

    err = lis_vector_set_all(
        1.0,
        x);
    if( err ) goto cleanup;

    err = lis_esolver_create(
        &esolver);
    if( err ) goto cleanup;

    err = lis_esolver_set_option(
        "-e li -ss 4 -eblock 2 -rval false -etol 1.0e-12",
        esolver);
    if( err ) goto cleanup;

    api = lis_esolve(
        A,
        x,
        &evalue,
        esolver);

    if( api!=LIS_SUCCESS )
    {
        printf(
            "Block Lanczos returned-x compatibility: "
            "API return %d\n",
            (int)api);

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
            "Block Lanczos returned-x compatibility: "
            "status %d, expected LIS_SUCCESS\n",
            (int)status);

        err = LIS_FAILS;
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A,
        &evector0);
    if( err ) goto cleanup;

    err = lis_vector_duplicate(
        A,
        &diff);
    if( err ) goto cleanup;

    err = lis_esolver_get_specific_evector(
        esolver,
        0,
        evector0);
    if( err ) goto cleanup;

    err = lis_vector_copy(
        x,
        diff);
    if( err ) goto cleanup;

    err = lis_vector_axpy(
        (LIS_SCALAR)-1.0,
        evector0,
        diff);
    if( err ) goto cleanup;

    err = lis_vector_nrm2(
        diff,
        &defect);
    if( err ) goto cleanup;

    if( defect>(LIS_REAL)1.0e-12 )
    {
        printf(
            "Block Lanczos returned-x compatibility: "
            "||x-evector0||=%e\n",
            (double)defect);

        err = LIS_FAILS;
        goto cleanup;
    }

    printf(
        "PASS: Block Lanczos returned x equals evector0 "
        "defect=%e\n",
        (double)defect);

    err = LIS_SUCCESS;

cleanup:
    if( diff )
        lis_vector_destroy(diff);

    if( evector0 )
        lis_vector_destroy(evector0);

    if( esolver )
        lis_esolver_destroy(esolver);

    if( x )
        lis_vector_destroy(x);

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static LIS_INT
test_block_lanczos_final_compatibility(void)
{
    LIS_INT err;

    err = test_block_lanczos_scalar_compatibility();
    if( err ) return err;

    err = test_block_lanczos_raw_ritz_compatibility();
    if( err ) return err;

    err = test_block_lanczos_returned_x_compatibility();
    if( err ) return err;

    printf(
        "PASS: Block Lanczos final compatibility regressions\n");

    return LIS_SUCCESS;
}


static LIS_INT
test_block_lanczos_rval_false_boundary(void)
{
    LIS_INT err;

    err =
        test_block_lanczos_ritz_evidence_case(
            4,
            0,
            4);
    if( err ) return err;

    err =
        test_block_lanczos_ritz_evidence_case(
            5,
            0,
            5);
    if( err ) return err;

    err =
        test_block_lanczos_ritz_evidence_case(
            6,
            0,
            6);
    if( err ) return err;

    err =
        test_block_lanczos_ritz_evidence_case(
            6,
            1,
            2);
    if( err ) return err;

    err =
        test_block_lanczos_residual_status_boundary();
    if( err ) return err;

    err =
        test_block_lanczos_clustered_repeated_spectra();
    if( err ) return err;

    err =
        test_block_lanczos_final_compatibility();
    if( err ) return err;

    printf(
        "PASS: Block Lanczos Ritz orthogonality, residuals, "
        "and active result counts\n");

    printf(
        "PASS: Block Lanczos rval=false reconstructs Ritz vectors\n");

    return LIS_SUCCESS;
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

    err = test_block_lanczos_projected_storage();

    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN BLOCK LANCZOS "
            "PROJECTED STORAGE FAILED\n");

        lis_finalize();

        return 1;
    }

    err =
        test_block_lanczos_two_block_runtime();

    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN BLOCK LANCZOS "
            "TWO-BLOCK RUNTIME FAILED\n");

        lis_finalize();

        return 1;
    }

    err =
        test_block_lanczos_multiblock_runtime();

    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN BLOCK LANCZOS "
            "MULTI-BLOCK RUNTIME FAILED\n");

        lis_finalize();

        return 1;
    }

    err =
        test_block_lanczos_projected_spectrum();

    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN BLOCK LANCZOS "
            "PROJECTED SPECTRUM FAILED\n");

        lis_finalize();

        return 1;
    }

    err = test_block_lanczos_rval_false_boundary();
    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN BLOCK LANCZOS RVAL FALSE BOUNDARY FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_block_lanczos_initial_raw();

    if( err )
    {
        printf(
            "LIS_EIGEN_BREAKDOWN BLOCK LANCZOS INITIAL FAILED\n");
        lis_finalize();
        return 1;
    }

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

    err = test_lanczos_block_option_contract();
    if( err )
    {
        printf(
            "LIS_EIGEN_BLOCK_LANCZOS FOUNDATION FAILED\n");
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

    err = test_lanczos_clustered_basis_stability();
    if( err )
    {
        printf(
            "LIS_EIGEN_ORTHOGONALITY "
            "CLUSTERED LANCZOS BASIS FAILED\n");
        lis_finalize();
        return 1;
    }

    err = test_transformed_generalized_lanczos_clustered_spectrum();
    if( err )
    {
        printf(
            "LIS_EIGEN_DIAGNOSTIC "
            "TRANSFORMED GENERALIZED LANCZOS FAILED\n");
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
