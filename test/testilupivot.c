#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lislib.h"

#include <math.h>
#include <stdio.h>

#define TEST_N       8
#define TEST_BLOCKS  4

static const LIS_REAL block_scale[TEST_BLOCKS] = {
        (LIS_REAL)1.0e-12,
        (LIS_REAL)1.0,
        (LIS_REAL)1.0e12,
        (LIS_REAL)1.0e6
};

static const LIS_REAL relative_delta =
        (LIS_REAL)1.0e-12;

static const LIS_REAL pivot_tol =
        (LIS_REAL)1.0e-6;


static LIS_REAL
real_abs(LIS_REAL x)
{
        return x < (LIS_REAL)0.0 ? -x : x;
}


static LIS_REAL
scalar_abs(LIS_SCALAR x)
{
        return (LIS_REAL)fabs(x);
}


/*
 * Build four independent 2x2 blocks
 *
 *              [ 1       1       ]
 *      A_b = s [                   ]
 *              [ 1   1 + delta    ]
 *
 * The second elimination pivot is approximately s*delta,
 * while the original row infinity norm is approximately s.
 *
 * Therefore a relative pivot tolerance tau must regularize
 * the second pivot to approximately
 *
 *      tau * ||row||_inf
 *
 * independently of the absolute value of s.
 */
static LIS_INT
create_test_matrix(LIS_MATRIX *A)
{
        LIS_INT err;
        LIS_INT is,ie,i,b;
        LIS_SCALAR s;

        *A = NULL;

        err = lis_matrix_create(LIS_COMM_WORLD,A);
        if( err ) return err;

        err = lis_matrix_set_size(*A,0,TEST_N);
        if( err ) return err;

        err = lis_matrix_get_range(*A,&is,&ie);
        if( err ) return err;

        for(i=is;i<ie;i++)
        {
                b = i/2;
                s = (LIS_SCALAR)block_scale[b];

                if( (i%2)==0 )
                {
                        err = lis_matrix_set_value(
                                LIS_INS_VALUE,
                                i,i,
                                s,
                                *A);
                        if( err ) return err;

                        err = lis_matrix_set_value(
                                LIS_INS_VALUE,
                                i,i+1,
                                s,
                                *A);
                        if( err ) return err;
                }
                else
                {
                        err = lis_matrix_set_value(
                                LIS_INS_VALUE,
                                i,i-1,
                                s,
                                *A);
                        if( err ) return err;

                        err = lis_matrix_set_value(
                                LIS_INS_VALUE,
                                i,i,
                                s*(LIS_SCALAR)
                                  ((LIS_REAL)1.0 + relative_delta),
                                *A);
                        if( err ) return err;
                }
        }

        err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
        if( err ) return err;

        return lis_matrix_assemble(*A);
}


static int
check_regularized_diagonal(const char *name,
                           LIS_PRECON precon)
{
        LIS_INT b,i;
        LIS_SCALAR inverse_pivot;
        LIS_SCALAR pivot;
        LIS_REAL actual;
        LIS_REAL row_scale;
        LIS_REAL expected;
        LIS_REAL relerr;

        if( precon==NULL || precon->D==NULL )
        {
                fprintf(stderr,
                        "%s: missing inverse pivot vector\n",
                        name);
                return 1;
        }

        if( precon->D->n!=TEST_N )
        {
                fprintf(stderr,
                        "%s: test requires one process "
                        "(local n=%d, expected %d)\n",
                        name,
                        (int)precon->D->n,
                        TEST_N);
                return 1;
        }

        for(b=0;b<TEST_BLOCKS;b++)
        {
                i = 2*b + 1;

                inverse_pivot = precon->D->value[i];

                if( !(scalar_abs(inverse_pivot) >
                      (LIS_REAL)0.0) )
                {
                        fprintf(stderr,
                                "%s block %d: invalid inverse pivot\n",
                                name,(int)b);
                        return 1;
                }

                pivot = (LIS_SCALAR)1.0 / inverse_pivot;
                actual = scalar_abs(pivot);

                row_scale =
                        block_scale[b] *
                        ((LIS_REAL)1.0 + relative_delta);

                expected = pivot_tol * row_scale;

                if( !(actual > (LIS_REAL)0.0) )
                {
                        fprintf(stderr,
                                "%s block %d: "
                                "non-finite or zero pivot\n",
                                name,(int)b);
                        return 1;
                }

                relerr =
                        real_abs(actual-expected) /
                        expected;

                if( !(relerr <= (LIS_REAL)5.0e-6) )
                {
                        fprintf(stderr,
                                "%s block %d: "
                                "pivot=% .15e expected=% .15e "
                                "relerr=% .15e\n",
                                name,
                                (int)b,
                                (double)actual,
                                (double)expected,
                                (double)relerr);
                        return 1;
                }
        }

        return 0;
}


static int
run_preconditioner(LIS_MATRIX A,
                   const char *name)
{
        LIS_SOLVER solver = NULL;
        LIS_PRECON precon = NULL;
        LIS_INT err;
        int failed = 0;
        char options[256];

        err = lis_solver_create(&solver);
        if( err )
        {
                fprintf(stderr,
                        "%s: lis_solver_create failed: %d\n",
                        name,(int)err);
                return 1;
        }

        solver->A = A;

        snprintf(
                options,sizeof(options),
                "-p %s "
                "-ilu_fill 0 "
                "-iluc_drop 0 "
                "-iluc_rate 5 "
                "-scale none "
                "-ilu_pivot_tol %.17g",
                name,
                (double)pivot_tol);

        err = lis_solver_set_option(options,solver);

        if( err )
        {
                fprintf(stderr,
                        "%s: -ilu_pivot_tol rejected: %d\n",
                        name,(int)err);
                failed = 1;
                goto cleanup;
        }

        precon = NULL;
        err = lis_precon_create(solver,&precon);

        if( err )
        {
                fprintf(stderr,
                        "%s: preconditioner creation failed: %d\n",
                        name,(int)err);
                failed = 1;
                goto cleanup;
        }

        failed |=
                check_regularized_diagonal(name,precon);

cleanup:

        if( precon )
                lis_precon_destroy(precon);

        if( solver )
                lis_solver_destroy(solver);

        return failed;
}


static LIS_INT
set_test_options(LIS_SOLVER solver,
                 const char *name,
                 const char *pivot_value)
{
        char options[256];

        if( pivot_value )
        {
                snprintf(
                        options,sizeof(options),
                        "-p %s "
                        "-ilu_fill 0 "
                        "-iluc_drop 0 "
                        "-iluc_rate 5 "
                        "-scale none "
                        "-ilu_pivot_tol %s",
                        name,
                        pivot_value);
        }
        else
        {
                snprintf(
                        options,sizeof(options),
                        "-p %s "
                        "-ilu_fill 0 "
                        "-iluc_drop 0 "
                        "-iluc_rate 5 "
                        "-scale none",
                        name);
        }

        return lis_solver_set_option(options,solver);
}


static int
check_default_parameter(void)
{
        LIS_SOLVER solver = NULL;
        LIS_INT err;
        LIS_REAL value;

        err = lis_solver_create(&solver);

        if( err )
        {
                fprintf(stderr,
                        "default parameter: "
                        "lis_solver_create failed: %d\n",
                        (int)err);
                return 1;
        }

        value = (LIS_REAL)solver->params[
                LIS_PARAMS_ILU_PIVOT_TOL-LIS_OPTIONS_LEN];

        if( value!=(LIS_REAL)0.0 )
        {
                fprintf(stderr,
                        "default pivot tolerance is %.17g, "
                        "expected 0\n",
                        (double)value);

                lis_solver_destroy(solver);
                return 1;
        }

        lis_solver_destroy(solver);

        return 0;
}


static int
capture_inverse_diagonal(const char *name,
                         const char *pivot_value,
                         LIS_SCALAR values[TEST_N])
{
        LIS_MATRIX A = NULL;
        LIS_SOLVER solver = NULL;
        LIS_PRECON precon = NULL;
        LIS_INT err;
        LIS_INT i;
        int failed = 0;

        err = create_test_matrix(&A);

        if( err )
        {
                fprintf(stderr,
                        "%s compatibility: "
                        "matrix creation failed: %d\n",
                        name,(int)err);
                return 1;
        }

        if( A->n!=TEST_N )
        {
                fprintf(stderr,
                        "%s compatibility: "
                        "test requires one MPI process\n",
                        name);
                failed = 1;
                goto cleanup;
        }

        err = lis_solver_create(&solver);

        if( err )
        {
                fprintf(stderr,
                        "%s compatibility: "
                        "solver creation failed: %d\n",
                        name,(int)err);
                failed = 1;
                goto cleanup;
        }

        solver->A = A;

        err = set_test_options(
                solver,
                name,
                pivot_value);

        if( err )
        {
                fprintf(stderr,
                        "%s compatibility: "
                        "option setup failed: %d\n",
                        name,(int)err);
                failed = 1;
                goto cleanup;
        }

        err = lis_precon_create(
                solver,
                &precon);

        if( err )
        {
                fprintf(stderr,
                        "%s compatibility: "
                        "preconditioner creation failed: %d\n",
                        name,(int)err);

                /*
                 * lis_precon_create destroys the partially created
                 * preconditioner on failure.
                 */
                precon = NULL;

                failed = 1;
                goto cleanup;
        }

        if( precon->D==NULL ||
            precon->D->n!=TEST_N )
        {
                fprintf(stderr,
                        "%s compatibility: "
                        "invalid inverse diagonal\n",
                        name);
                failed = 1;
                goto cleanup;
        }

        for(i=0;i<TEST_N;i++)
        {
                values[i] = precon->D->value[i];
        }

cleanup:

        if( precon )
                lis_precon_destroy(precon);

        if( solver )
                lis_solver_destroy(solver);

        if( A )
                lis_matrix_destroy(A);

        return failed;
}


static int
check_default_zero_compatibility(const char *name)
{
        LIS_SCALAR default_diag[TEST_N];
        LIS_SCALAR zero_diag[TEST_N];
        LIS_REAL a,b,diff,scale;
        LIS_INT i;
        int failed;

        failed = capture_inverse_diagonal(
                name,
                NULL,
                default_diag);

        if( failed )
                return 1;

        failed = capture_inverse_diagonal(
                name,
                "0",
                zero_diag);

        if( failed )
                return 1;

        for(i=0;i<TEST_N;i++)
        {
                a = scalar_abs(default_diag[i]);
                b = scalar_abs(zero_diag[i]);

                scale = a>b ? a : b;

                if( scale<(LIS_REAL)1.0 )
                        scale = (LIS_REAL)1.0;

                diff = scalar_abs(
                        default_diag[i] -
                        zero_diag[i]);

                if( diff >
                    (LIS_REAL)1.0e-12*scale )
                {
                        fprintf(stderr,
                                "%s compatibility row %d: "
                                "default and explicit zero differ "
                                "(rel=% .15e)\n",
                                name,
                                (int)i,
                                (double)(diff/scale));

                        return 1;
                }
        }

        return 0;
}


static int
check_negative_tolerance(const char *name)
{
        LIS_MATRIX A = NULL;
        LIS_SOLVER solver = NULL;
        LIS_PRECON precon = NULL;
        LIS_INT err;
        int failed = 0;

        err = create_test_matrix(&A);

        if( err )
        {
                fprintf(stderr,
                        "%s negative tolerance: "
                        "matrix creation failed: %d\n",
                        name,(int)err);
                return 1;
        }

        err = lis_solver_create(&solver);

        if( err )
        {
                fprintf(stderr,
                        "%s negative tolerance: "
                        "solver creation failed: %d\n",
                        name,(int)err);
                failed = 1;
                goto cleanup;
        }

        solver->A = A;

        err = set_test_options(
                solver,
                name,
                "-1e-6");

        /*
         * It is valid for option validation to reject the value
         * immediately.  Otherwise preconditioner creation must
         * reject it with LIS_ERR_ILL_ARG.
         */
        if( err==LIS_ERR_ILL_ARG )
        {
                goto cleanup;
        }

        if( err )
        {
                fprintf(stderr,
                        "%s negative tolerance: "
                        "unexpected option error: %d\n",
                        name,(int)err);
                failed = 1;
                goto cleanup;
        }

        err = lis_precon_create(
                solver,
                &precon);

        if( err==LIS_ERR_ILL_ARG )
        {
                precon = NULL;
                goto cleanup;
        }

        if( err )
        {
                fprintf(stderr,
                        "%s negative tolerance: "
                        "expected LIS_ERR_ILL_ARG, got %d\n",
                        name,(int)err);

                precon = NULL;
                failed = 1;
                goto cleanup;
        }

        fprintf(stderr,
                "%s negative tolerance: "
                "invalid value was accepted\n",
                name);

        failed = 1;

cleanup:

        if( precon )
                lis_precon_destroy(precon);

        if( solver )
                lis_solver_destroy(solver);

        if( A )
                lis_matrix_destroy(A);

        return failed;
}


static LIS_INT
create_zero_row_matrix(LIS_MATRIX *A)
{
        LIS_INT err;
        LIS_INT is,ie,i;
        LIS_SCALAR value;

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

        err = lis_matrix_get_range(
                *A,
                &is,
                &ie);

        if( err )
                return err;

        /*
         * Row 0 has scale 1.
         * Row 1 is structurally present but has zero scale.
         */
        for(i=is;i<ie;i++)
        {
                value =
                        i==0
                        ? (LIS_SCALAR)1.0
                        : (LIS_SCALAR)0.0;

                err = lis_matrix_set_value(
                        LIS_INS_VALUE,
                        i,
                        i,
                        value,
                        *A);

                if( err )
                        return err;
        }

        err = lis_matrix_set_type(
                *A,
                LIS_MATRIX_CSR);

        if( err )
                return err;

        return lis_matrix_assemble(*A);
}


static int
check_zero_row_breakdown(const char *name)
{
        LIS_MATRIX A = NULL;
        LIS_SOLVER solver = NULL;
        LIS_PRECON precon = NULL;
        LIS_INT err;
        int failed = 0;

        err = create_zero_row_matrix(&A);

        if( err )
        {
                fprintf(stderr,
                        "%s zero row: "
                        "matrix creation failed: %d\n",
                        name,(int)err);
                return 1;
        }

        if( A->n!=2 )
        {
                fprintf(stderr,
                        "%s zero row: "
                        "test requires one MPI process\n",
                        name);
                failed = 1;
                goto cleanup;
        }

        err = lis_solver_create(&solver);

        if( err )
        {
                fprintf(stderr,
                        "%s zero row: "
                        "solver creation failed: %d\n",
                        name,(int)err);
                failed = 1;
                goto cleanup;
        }

        solver->A = A;

        err = set_test_options(
                solver,
                name,
                "1e-6");

        if( err )
        {
                fprintf(stderr,
                        "%s zero row: "
                        "option setup failed: %d\n",
                        name,(int)err);
                failed = 1;
                goto cleanup;
        }

        err = lis_precon_create(
                solver,
                &precon);

        if( err==LIS_BREAKDOWN )
        {
                precon = NULL;
                goto cleanup;
        }

        if( err )
        {
                fprintf(stderr,
                        "%s zero row: "
                        "expected LIS_BREAKDOWN, got %d\n",
                        name,(int)err);

                precon = NULL;
                failed = 1;
                goto cleanup;
        }

        fprintf(stderr,
                "%s zero row: "
                "positive pivot tolerance accepted "
                "a zero-scale row\n",
                name);

        failed = 1;

cleanup:

        if( precon )
                lis_precon_destroy(precon);

        if( solver )
                lis_solver_destroy(solver);

        if( A )
                lis_matrix_destroy(A);

        return failed;
}


int
main(int argc, char *argv[])
{
        LIS_MATRIX A = NULL;
        LIS_INT err;
        int failed = 0;

        err = lis_initialize(&argc,&argv);

        if( err )
                return 1;

        /*
         * Public option contract:
         * omission means exactly zero.
         */
        failed |= check_default_parameter();

        /*
         * Positive relative pivot tolerance:
         * ILU(k), ILUT and ILUC must regularize relative
         * to the original row infinity norm.
         */
        err = create_test_matrix(&A);

        if( err )
        {
                fprintf(stderr,
                        "matrix creation failed: %d\n",
                        (int)err);

                lis_finalize();
                return 1;
        }

        if( A->n!=TEST_N )
        {
                fprintf(stderr,
                        "testilupivot currently requires "
                        "one MPI process\n");

                lis_matrix_destroy(A);
                lis_finalize();
                return 1;
        }

        failed |= run_preconditioner(A,"ilu");
        failed |= run_preconditioner(A,"ilut");
        failed |= run_preconditioner(A,"iluc");

        lis_matrix_destroy(A);
        A = NULL;

        /*
         * Backward compatibility:
         * omitted tolerance and explicit zero must produce
         * equivalent inverse pivots.
         */
        failed |=
                check_default_zero_compatibility("ilu");

        failed |=
                check_default_zero_compatibility("ilut");

        failed |=
                check_default_zero_compatibility("iluc");

        /*
         * Negative tolerance is invalid.
         */
        failed |=
                check_negative_tolerance("ilu");

        failed |=
                check_negative_tolerance("ilut");

        failed |=
                check_negative_tolerance("iluc");

        /*
         * A positive relative tolerance cannot be defined
         * for a row whose original infinity norm is zero.
         */
        failed |=
                check_zero_row_breakdown("ilu");

        failed |=
                check_zero_row_breakdown("ilut");

        failed |=
                check_zero_row_breakdown("iluc");

        if( failed )
        {
                fprintf(stderr,
                        "LIS_ILU_PIVOT_STAGE6B RED\n");

                lis_finalize();
                return 1;
        }

        printf(
                "LIS_ILU_PIVOT_STAGE6B PASSED\n");

        lis_finalize();

        return 0;
}
