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
