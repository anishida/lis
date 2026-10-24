#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include <stdio.h>
#include <math.h>

#include "lis.h"
#include "lis_esolver.h"


#define TEST_N 3


static LIS_REAL scalar_real_part(
        LIS_SCALAR value)
{
#ifdef _COMPLEX
        return (LIS_REAL)creal(value);
#else
        return (LIS_REAL)value;
#endif
}


static LIS_REAL scalar_imag_part(
        LIS_SCALAR value)
{
#ifdef _COMPLEX
        return (LIS_REAL)cimag(value);
#else
        (void)value;
        return (LIS_REAL)0.0;
#endif
}


static int close_real(
        LIS_REAL a,
        LIS_REAL b,
        LIS_REAL tol)
{
        LIS_REAL scale;

        scale =
                (LIS_REAL)1.0 +
                fabs(b);

        return
                fabs(a-b) <=
                tol*scale;
}


typedef struct
{
        LIS_INT n;
        LIS_INT kind;
} FEAST_USER_CONTEXT;


#define FEAST_USER_A 1
#define FEAST_USER_B 2


static LIS_INT feast_user_matvec(
        void *user_data,
        const LIS_SCALAR *x,
        LIS_SCALAR *y)
{
        FEAST_USER_CONTEXT *ctx;

        ctx =
                (FEAST_USER_CONTEXT *)user_data;

        if( ctx==NULL ||
            x==NULL ||
            y==NULL ||
            ctx->n!=TEST_N )
        {
                return LIS_ERR_ILL_ARG;
        }

        if( ctx->kind==FEAST_USER_A )
        {
                y[0] =
                        (LIS_SCALAR)4.0*x[0] +
                        x[1];

                y[1] =
                        x[0] +
                        (LIS_SCALAR)3.0*x[1] +
                        (LIS_SCALAR)2.0*x[2];

                y[2] =
                        (LIS_SCALAR)2.0*x[1] +
                        (LIS_SCALAR)5.0*x[2];

                return LIS_SUCCESS;
        }

        if( ctx->kind==FEAST_USER_B )
        {
                y[0] =
                        (LIS_SCALAR)2.0*x[0];

                y[1] =
                        x[1];

                y[2] =
                        (LIS_SCALAR)3.0*x[2];

                return LIS_SUCCESS;
        }

        return LIS_ERR_ILL_ARG;
}


static LIS_INT create_user_matrix(
        LIS_INT kind,
        FEAST_USER_CONTEXT *ctx,
        LIS_MATRIX *matrix)
{
        LIS_INT err;

        if( ctx==NULL ||
            matrix==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        *matrix = NULL;

        err = lis_matrix_create(
                LIS_COMM_WORLD,
                matrix);

        if( err )
        {
                return err;
        }

        err = lis_matrix_set_size(
                *matrix,
                0,
                TEST_N);

        if( err )
        {
                return err;
        }

        ctx->n =
                (*matrix)->n;

        ctx->kind =
                kind;

        return lis_matrix_set_user(
                *matrix,
                ctx,
                feast_user_matvec,
                feast_user_matvec);
}


static int check_matrix_type(
        const char *name,
        LIS_MATRIX matrix,
        LIS_INT expected)
{
        if( matrix==NULL )
        {
                fprintf(
                        stderr,
                        "FAIL: %s is NULL\n",
                        name);

                return 1;
        }

        if( matrix->matrix_type!=expected )
        {
                fprintf(
                        stderr,
                        "FAIL: %s matrix type changed: got %d expected %d\n",
                        name,
                        (int)matrix->matrix_type,
                        (int)expected);

                return 1;
        }

        return 0;
}


static LIS_INT create_matrix_a(
        LIS_MATRIX *A)
{
        static const LIS_REAL values[TEST_N][TEST_N] = {
                {
                        (LIS_REAL)4.0,
                        (LIS_REAL)1.0,
                        (LIS_REAL)0.0
                },
                {
                        (LIS_REAL)1.0,
                        (LIS_REAL)3.0,
                        (LIS_REAL)2.0
                },
                {
                        (LIS_REAL)0.0,
                        (LIS_REAL)2.0,
                        (LIS_REAL)5.0
                }
        };

        LIS_INT i;
        LIS_INT j;
        LIS_INT err;

        *A = NULL;

        err = lis_matrix_create(
                LIS_COMM_WORLD,
                A);

        if( err )
        {
                return err;
        }

        err = lis_matrix_set_size(
                *A,
                0,
                TEST_N);

        if( err )
        {
                return err;
        }

        for(i=0;i<TEST_N;i++)
        {
                for(j=0;j<TEST_N;j++)
                {
                        if( values[i][j]==0.0 )
                        {
                                continue;
                        }

                        err = lis_matrix_set_value(
                                LIS_INS_VALUE,
                                i,
                                j,
                                (LIS_SCALAR)values[i][j],
                                *A);

                        if( err )
                        {
                                return err;
                        }
                }
        }

        return lis_matrix_assemble(
                *A);
}


static LIS_INT create_matrix_b(
        LIS_MATRIX *B)
{
        static const LIS_REAL diag[TEST_N] = {
                (LIS_REAL)2.0,
                (LIS_REAL)1.0,
                (LIS_REAL)3.0
        };

        LIS_INT i;
        LIS_INT err;

        *B = NULL;

        err = lis_matrix_create(
                LIS_COMM_WORLD,
                B);

        if( err )
        {
                return err;
        }

        err = lis_matrix_set_size(
                *B,
                0,
                TEST_N);

        if( err )
        {
                return err;
        }

        for(i=0;i<TEST_N;i++)
        {
                err = lis_matrix_set_value(
                        LIS_INS_VALUE,
                        i,
                        i,
                        (LIS_SCALAR)diag[i],
                        *B);

                if( err )
                {
                        return err;
                }
        }

        return lis_matrix_assemble(
                *B);
}


static LIS_INT set_input_vectors(
        LIS_VECTOR x_real,
        LIS_VECTOR x_imag)
{
        static const LIS_REAL xr[TEST_N] = {
                (LIS_REAL)1.0,
                (LIS_REAL)-2.0,
                (LIS_REAL)0.5
        };

        static const LIS_REAL xi[TEST_N] = {
                (LIS_REAL)0.25,
                (LIS_REAL)1.0,
                (LIS_REAL)-1.5
        };

        LIS_INT i;
        LIS_INT err;

        for(i=0;i<TEST_N;i++)
        {
                err = lis_vector_set_value(
                        LIS_INS_VALUE,
                        i,
                        (LIS_SCALAR)xr[i],
                        x_real);

                if( err )
                {
                        return err;
                }

                err = lis_vector_set_value(
                        LIS_INS_VALUE,
                        i,
                        (LIS_SCALAR)xi[i],
                        x_imag);

                if( err )
                {
                        return err;
                }
        }

        return LIS_SUCCESS;
}


static int check_vector(
        const char *name,
        LIS_VECTOR vector,
        const LIS_REAL expected[TEST_N])
{
        LIS_SCALAR value;

        LIS_REAL real_value;
        LIS_REAL imag_value;

        LIS_INT i;
        LIS_INT err;

        for(i=0;i<TEST_N;i++)
        {
                err = lis_vector_get_value(
                        vector,
                        i,
                        &value);

                if( err )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s get value %d returned %d\n",
                                name,
                                (int)i,
                                (int)err);

                        return 1;
                }

                real_value =
                        scalar_real_part(
                                value);

                imag_value =
                        scalar_imag_part(
                                value);

                if( !close_real(
                            real_value,
                            expected[i],
                            (LIS_REAL)5.0e-13) )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s[%d] = %.17g expected %.17g\n",
                                name,
                                (int)i,
                                (double)real_value,
                                (double)expected[i]);

                        return 1;
                }

                /*
                 * Split vectors themselves must remain purely real,
                 * including when LIS is compiled with complex scalars.
                 */
                if( !close_real(
                            imag_value,
                            (LIS_REAL)0.0,
                            (LIS_REAL)5.0e-13) )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s[%d] has scalar imaginary part %.17g\n",
                                name,
                                (int)i,
                                (double)imag_value);

                        return 1;
                }
        }

        return 0;
}


static int run_shift_case(
        const char *name,
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_VECTOR x_real,
        LIS_VECTOR x_imag,
        LIS_VECTOR y_real,
        LIS_VECTOR y_imag,
        LIS_VECTOR work_real,
        LIS_VECTOR work_imag,
        const LIS_REAL expected_real[TEST_N],
        const LIS_REAL expected_imag[TEST_N])
{
        LIS_INT err;

        err = lis_efeast_shifted_matvec_split(
                A,
                B,
                (LIS_REAL)2.0,
                (LIS_REAL)0.75,
                x_real,
                x_imag,
                y_real,
                y_imag,
                work_real,
                work_imag);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s shifted matvec returned %d\n",
                        name,
                        (int)err);

                return 1;
        }

        if( check_vector(
                    "y_real",
                    y_real,
                    expected_real) )
        {
                return 1;
        }

        if( check_vector(
                    "y_imag",
                    y_imag,
                    expected_imag) )
        {
                return 1;
        }

        return 0;
}


static int check_vector_match(
        const char *name,
        LIS_VECTOR actual,
        LIS_VECTOR expected,
        LIS_REAL tol)
{
        LIS_SCALAR actual_value;
        LIS_SCALAR expected_value;

        LIS_REAL actual_real;
        LIS_REAL actual_imag;
        LIS_REAL expected_real;
        LIS_REAL expected_imag;

        LIS_INT i;
        LIS_INT err;

        for(i=0;i<TEST_N;i++)
        {
                err = lis_vector_get_value(
                        actual,
                        i,
                        &actual_value);

                if( err )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s actual get value %d returned %d\n",
                                name,
                                (int)i,
                                (int)err);

                        return 1;
                }

                err = lis_vector_get_value(
                        expected,
                        i,
                        &expected_value);

                if( err )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s expected get value %d returned %d\n",
                                name,
                                (int)i,
                                (int)err);

                        return 1;
                }

                actual_real =
                        scalar_real_part(
                                actual_value);

                actual_imag =
                        scalar_imag_part(
                                actual_value);

                expected_real =
                        scalar_real_part(
                                expected_value);

                expected_imag =
                        scalar_imag_part(
                                expected_value);

                if( !close_real(
                            actual_real,
                            expected_real,
                            tol) ||
                    !close_real(
                            actual_imag,
                            expected_imag,
                            tol) )
                {
                        fprintf(
                                stderr,
                                "FAIL: %s[%d] = "
                                "(%.17g, %.17g), expected "
                                "(%.17g, %.17g)\n",
                                name,
                                (int)i,
                                (double)actual_real,
                                (double)actual_imag,
                                (double)expected_real,
                                (double)expected_imag);

                        return 1;
                }
        }

        return 0;
}


static int run_shift_solve_case(
        const char *name,
        LIS_MATRIX A,
        LIS_MATRIX B,
        LIS_VECTOR exact_real,
        LIS_VECTOR exact_imag,
        LIS_VECTOR rhs_real,
        LIS_VECTOR rhs_imag,
        LIS_VECTOR work_real,
        LIS_VECTOR work_imag)
{
        LIS_VECTOR solve_real;
        LIS_VECTOR solve_imag;

        LIS_SOLVER solver;

        LIS_INT iter;
        LIS_INT err;

        LIS_REAL residual;

        solve_real = NULL;
        solve_imag = NULL;

        solver = NULL;

        iter = 0;
        residual = 0.0;

        /*
         * Build an exact RHS from the known split-complex solution.
         */
        err = lis_efeast_shifted_matvec_split(
                A,
                B,
                (LIS_REAL)2.0,
                (LIS_REAL)0.75,
                exact_real,
                exact_imag,
                rhs_real,
                rhs_imag,
                work_real,
                work_imag);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s RHS construction returned %d\n",
                        name,
                        (int)err);

                return 1;
        }

        err = lis_vector_duplicate(
                A,
                &solve_real);

        if( err ) goto fail;

        err = lis_vector_duplicate(
                A,
                &solve_imag);

        if( err ) goto fail;

        err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                solve_real);

        if( err ) goto fail;

        err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                solve_imag);

        if( err ) goto fail;

        err = lis_solver_create(
                &solver);

        if( err ) goto fail;

        /*
         * The F3-2 bridge is deliberately matrix-free.
         * Use a nonsymmetric Krylov method that requires matvec only.
         */
        solver->options[
                LIS_OPTIONS_SOLVER] =
                LIS_SOLVER_GMRES;

        solver->options[
                LIS_OPTIONS_PRECON] =
                LIS_PRECON_TYPE_NONE;

        solver->options[
                LIS_OPTIONS_SCALE] =
                LIS_SCALE_NONE;

        solver->options[
                LIS_OPTIONS_STORAGE] =
                0;

        solver->options[
                LIS_OPTIONS_USE_AT] =
                LIS_FALSE;

        solver->options[
                LIS_OPTIONS_PRECISION] =
                LIS_PRECISION_DOUBLE;

        solver->options[
                LIS_OPTIONS_MAXITER] =
                200;

        solver->options[
                LIS_OPTIONS_RESTART] =
                12;

        solver->options[
                LIS_OPTIONS_INITGUESS_ZEROS] =
                LIS_TRUE;

        solver->params[
                LIS_PARAMS_RESID
                -
                LIS_OPTIONS_LEN] =
                (LIS_SCALAR)1.0e-12;

        err = lis_efeast_shifted_solve_split(
                A,
                B,
                (LIS_REAL)2.0,
                (LIS_REAL)0.75,
                rhs_real,
                rhs_imag,
                solve_real,
                solve_imag,
                solver);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: %s shifted solve returned %d\n",
                        name,
                        (int)err);

                goto fail;
        }

        err = lis_solver_get_iter(
                solver,
                &iter);

        if( err ) goto fail;

        err = lis_solver_get_residualnorm(
                solver,
                &residual);

        if( err ) goto fail;

        if( iter<=0 ||
            !isfinite((double)residual) ||
            residual>(LIS_REAL)1.0e-9 )
        {
                fprintf(
                        stderr,
                        "FAIL: %s solve diagnostics "
                        "iter=%d residual=%.17g\n",
                        name,
                        (int)iter,
                        (double)residual);

                goto fail;
        }

        if( check_vector_match(
                    "solve_real",
                    solve_real,
                    exact_real,
                    (LIS_REAL)2.0e-9) )
        {
                goto fail;
        }

        if( check_vector_match(
                    "solve_imag",
                    solve_imag,
                    exact_imag,
                    (LIS_REAL)2.0e-9) )
        {
                goto fail;
        }

        printf(
                "PASS: %s shifted solve "
                "iter=%d residual=%.6e\n",
                name,
                (int)iter,
                (double)residual);

        if( solver!=NULL )
        {
                lis_solver_destroy(
                        solver);
        }

        if( solve_imag!=NULL )
        {
                lis_vector_destroy(
                        solve_imag);
        }

        if( solve_real!=NULL )
        {
                lis_vector_destroy(
                        solve_real);
        }

        return 0;

fail:

        if( solver!=NULL )
        {
                lis_solver_destroy(
                        solver);
        }

        if( solve_imag!=NULL )
        {
                lis_vector_destroy(
                        solve_imag);
        }

        if( solve_real!=NULL )
        {
                lis_vector_destroy(
                        solve_real);
        }

        return 1;
}


int main(
        int argc,
        char **argv)
{
        static const LIS_REAL generalized_real[TEST_N] = {
                (LIS_REAL)1.625,
                (LIS_REAL)-0.75,
                (LIS_REAL)7.875
        };

        static const LIS_REAL generalized_imag[TEST_N] = {
                (LIS_REAL)0.5,
                (LIS_REAL)0.25,
                (LIS_REAL)-2.375
        };

        static const LIS_REAL standard_real[TEST_N] = {
                (LIS_REAL)-0.1875,
                (LIS_REAL)-0.75,
                (LIS_REAL)3.625
        };

        static const LIS_REAL standard_imag[TEST_N] = {
                (LIS_REAL)-0.75,
                (LIS_REAL)0.25,
                (LIS_REAL)2.875
        };

        LIS_MATRIX A;
        LIS_MATRIX B;

        LIS_MATRIX user_a;
        LIS_MATRIX user_b;

        LIS_MATRIX operator_a;
        LIS_MATRIX operator_b;

        FEAST_USER_CONTEXT user_a_ctx;
        FEAST_USER_CONTEXT user_b_ctx;

        LIS_VECTOR x_real;
        LIS_VECTOR x_imag;
        LIS_VECTOR y_real;
        LIS_VECTOR y_imag;
        LIS_VECTOR work_real;
        LIS_VECTOR work_imag;

        LIS_INT err;

        int failed;

        A = NULL;
        B = NULL;

        user_a = NULL;
        user_b = NULL;

        operator_a = NULL;
        operator_b = NULL;

        user_a_ctx.n = 0;
        user_a_ctx.kind = 0;

        user_b_ctx.n = 0;
        user_b_ctx.kind = 0;

        x_real = NULL;
        x_imag = NULL;
        y_real = NULL;
        y_imag = NULL;
        work_real = NULL;
        work_imag = NULL;

        failed = 0;

        err = lis_initialize(
                &argc,
                &argv);

        if( err )
        {
                return (int)err;
        }

        err = create_matrix_a(
                &A);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: create A returned %d\n",
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        err = create_matrix_b(
                &B);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: create B returned %d\n",
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        err = create_user_matrix(
                FEAST_USER_A,
                &user_a_ctx,
                &user_a);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: create USER A returned %d\n",
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        err = create_user_matrix(
                FEAST_USER_B,
                &user_b_ctx,
                &user_b);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: create USER B returned %d\n",
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        /*
         * Build exact operator aliases of A and B:
         *
         * operator_a = 1*A + 0*B
         * operator_b = 0*A + 1*B
         *
         * They exercise the LIS_MATRIX_OPERATOR matvec path without
         * changing the analytical shifted system.
         */
        err = lis_matrix_create_operator(
                (LIS_SCALAR)1.0,
                A,
                (LIS_SCALAR)0.0,
                B,
                &operator_a);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: create OPERATOR A returned %d\n",
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        err = lis_matrix_create_operator(
                (LIS_SCALAR)0.0,
                A,
                (LIS_SCALAR)1.0,
                B,
                &operator_b);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: create OPERATOR B returned %d\n",
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        if( check_matrix_type(
                    "USER A",
                    user_a,
                    LIS_MATRIX_USER) ||
            check_matrix_type(
                    "USER B",
                    user_b,
                    LIS_MATRIX_USER) ||
            check_matrix_type(
                    "OPERATOR A",
                    operator_a,
                    LIS_MATRIX_OPERATOR) ||
            check_matrix_type(
                    "OPERATOR B",
                    operator_b,
                    LIS_MATRIX_OPERATOR) )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &x_real);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &x_imag);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &y_real);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &y_imag);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &work_real);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &work_imag);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = set_input_vectors(
                x_real,
                x_imag);

        if( err )
        {
                fprintf(
                        stderr,
                        "FAIL: set input vectors returned %d\n",
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        /*
         * Generalized problem:
         *
         *     (zB-A)x
         */
        if( run_shift_case(
                    "generalized",
                    A,
                    B,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag,
                    generalized_real,
                    generalized_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * Standard problem:
         *
         *     (zI-A)x
         */
        if( run_shift_case(
                    "standard",
                    A,
                    NULL,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag,
                    standard_real,
                    standard_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * USER + USER
         */
        if( run_shift_case(
                    "USER+USER",
                    user_a,
                    user_b,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag,
                    generalized_real,
                    generalized_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * USER + CSR
         */
        if( run_shift_case(
                    "USER+CSR",
                    user_a,
                    B,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag,
                    generalized_real,
                    generalized_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * CSR + USER
         */
        if( run_shift_case(
                    "CSR+USER",
                    A,
                    user_b,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag,
                    generalized_real,
                    generalized_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * OPERATOR + OPERATOR
         */
        if( run_shift_case(
                    "OPERATOR+OPERATOR",
                    operator_a,
                    operator_b,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag,
                    generalized_real,
                    generalized_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * USER standard problem, B = I.
         */
        if( run_shift_case(
                    "USER standard",
                    user_a,
                    NULL,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag,
                    standard_real,
                    standard_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * OPERATOR standard problem, B = I.
         */
        if( run_shift_case(
                    "OPERATOR standard",
                    operator_a,
                    NULL,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag,
                    standard_real,
                    standard_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * F3-2 shifted linear solves.
         *
         * The exact solution is the existing x_real/x_imag pair.
         * RHS vectors are generated through the F2 shifted matvec,
         * then recovered through the matrix-free real 2N bridge.
         */

        if( run_shift_solve_case(
                    "CSR+CSR",
                    A,
                    B,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        if( run_shift_solve_case(
                    "USER+USER",
                    user_a,
                    user_b,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        if( run_shift_solve_case(
                    "OPERATOR+OPERATOR",
                    operator_a,
                    operator_b,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        if( run_shift_solve_case(
                    "CSR standard",
                    A,
                    NULL,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        if( run_shift_solve_case(
                    "USER standard",
                    user_a,
                    NULL,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        if( run_shift_solve_case(
                    "OPERATOR standard",
                    operator_a,
                    NULL,
                    x_real,
                    x_imag,
                    y_real,
                    y_imag,
                    work_real,
                    work_imag) )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * The FEAST shifted path must not materialize or convert
         * USER/OPERATOR matrices to CSR as a side effect.
         */
        if( check_matrix_type(
                    "USER A after shifted matvec",
                    user_a,
                    LIS_MATRIX_USER) ||
            check_matrix_type(
                    "USER B after shifted matvec",
                    user_b,
                    LIS_MATRIX_USER) ||
            check_matrix_type(
                    "OPERATOR A after shifted matvec",
                    operator_a,
                    LIS_MATRIX_OPERATOR) ||
            check_matrix_type(
                    "OPERATOR B after shifted matvec",
                    operator_b,
                    LIS_MATRIX_OPERATOR) )
        {
                failed = 1;
                goto cleanup;
        }

cleanup:

        if( work_imag!=NULL )
        {
                lis_vector_destroy(
                        work_imag);
        }

        if( work_real!=NULL )
        {
                lis_vector_destroy(
                        work_real);
        }

        if( y_imag!=NULL )
        {
                lis_vector_destroy(
                        y_imag);
        }

        if( y_real!=NULL )
        {
                lis_vector_destroy(
                        y_real);
        }

        if( x_imag!=NULL )
        {
                lis_vector_destroy(
                        x_imag);
        }

        if( x_real!=NULL )
        {
                lis_vector_destroy(
                        x_real);
        }

        if( operator_b!=NULL )
        {
                lis_matrix_destroy(
                        operator_b);
        }

        if( operator_a!=NULL )
        {
                lis_matrix_destroy(
                        operator_a);
        }

        if( user_b!=NULL )
        {
                lis_matrix_destroy(
                        user_b);
        }

        if( user_a!=NULL )
        {
                lis_matrix_destroy(
                        user_a);
        }

        if( B!=NULL )
        {
                lis_matrix_destroy(
                        B);
        }

        if( A!=NULL )
        {
                lis_matrix_destroy(
                        A);
        }

        lis_finalize();

        if( failed )
        {
                return 1;
        }

        printf(
                "PASS: native FEAST split shifted operator/solve "
                "CSR/USER/OPERATOR\n");

        return 0;
}
