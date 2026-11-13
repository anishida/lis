#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include "lis.h"
#include <stdio.h>

typedef struct
{
        LIS_INT n;
        LIS_INT gn;
        LIS_INT is;

        LIS_SCALAR diag;
        LIS_SCALAR offdiag;

        unsigned long matvec_calls;
        unsigned long diagonal_calls;
        unsigned long row_query_calls;
        unsigned long row_fill_calls;
} USER_PRECON_CTX;


static LIS_INT
user_matvec(
        void *user_data,
        const LIS_SCALAR *x,
        LIS_SCALAR *y)
{
        USER_PRECON_CTX *ctx;
        LIS_INT i;

        ctx = (USER_PRECON_CTX *)user_data;

        ctx->matvec_calls++;

        if( ctx->n==1 )
        {
                y[0] = ctx->diag*x[0];
                return LIS_SUCCESS;
        }

        y[0] =
                ctx->diag*x[0]
                +
                ctx->offdiag*x[1];

        for(i=1;i<ctx->n-1;i++)
        {
                y[i] =
                        ctx->offdiag*x[i-1]
                        +
                        ctx->diag*x[i]
                        +
                        ctx->offdiag*x[i+1];
        }

        y[ctx->n-1] =
                ctx->offdiag*x[ctx->n-2]
                +
                ctx->diag*x[ctx->n-1];

        return LIS_SUCCESS;
}


static LIS_INT
user_get_diagonal(
        void *user_data,
        LIS_SCALAR *d)
{
        USER_PRECON_CTX *ctx;
        LIS_INT i;

        ctx = (USER_PRECON_CTX *)user_data;

        ctx->diagonal_calls++;

        for(i=0;i<ctx->n;i++)
        {
                d[i] = ctx->diag;
        }

        return LIS_SUCCESS;
}


static LIS_INT
user_get_row(
        void *user_data,
        LIS_INT local_row,
        LIS_INT *nnz,
        LIS_INT *cols,
        LIS_SCALAR *vals)
{
        USER_PRECON_CTX *ctx;
        LIS_INT global_row;
        LIS_INT required;
        LIS_INT capacity;
        LIS_INT k;

        ctx = (USER_PRECON_CTX *)user_data;

        if( nnz==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        if( local_row<0 || local_row>=ctx->n )
        {
                return LIS_ERR_ILL_ARG;
        }

        global_row = ctx->is + local_row;

        required = 1;

        if( global_row>0 )
        {
                required++;
        }

        if( global_row+1<ctx->gn )
        {
                required++;
        }

        if( cols==NULL && vals==NULL )
        {
                ctx->row_query_calls++;

                *nnz = required;

                return LIS_SUCCESS;
        }

        if( cols==NULL || vals==NULL )
        {
                return LIS_ERR_ILL_ARG;
        }

        ctx->row_fill_calls++;

        capacity = *nnz;

        if( capacity<required )
        {
                *nnz = required;
                return LIS_ERR_ILL_ARG;
        }

        k = 0;

        if( global_row>0 )
        {
                cols[k] = global_row-1;
                vals[k] = ctx->offdiag;
                k++;
        }

        cols[k] = global_row;
        vals[k] = ctx->diag;
        k++;

        if( global_row+1<ctx->gn )
        {
                cols[k] = global_row+1;
                vals[k] = ctx->offdiag;
                k++;
        }

        *nnz = k;

        return LIS_SUCCESS;
}


static void
reset_counters(
        USER_PRECON_CTX *ctx)
{
        ctx->matvec_calls      = 0;
        ctx->diagonal_calls    = 0;
        ctx->row_query_calls   = 0;
        ctx->row_fill_calls    = 0;
}


static int
run_case(
        const char *name,
        const char *options,
        int use_diagonal,
        int use_row,
        int expect_success,
        int require_row_calls,
        int require_diagonal_calls)
{
        const LIS_INT requested_n = 64;

        USER_PRECON_CTX ctx;

        LIS_MATRIX A;
        LIS_VECTOR b;
        LIS_VECTOR x;
        LIS_VECTOR exact;
        LIS_SOLVER solver;

        LIS_INT err;
        LIS_INT status;
        LIS_INT local_n;
        LIS_INT global_n;
        LIS_INT is;
        LIS_INT ie;

        LIS_REAL error_norm;

        int failed;

        A      = NULL;
        b      = NULL;
        x      = NULL;
        exact  = NULL;
        solver = NULL;

        failed = 0;

        err = lis_matrix_create(
                LIS_COMM_WORLD,
                &A);

        if( err )
        {
                fprintf(
                        stderr,
                        "%s: lis_matrix_create failed: %d\n",
                        name,
                        (int)err);

                return 1;
        }

        err = lis_matrix_set_size(
                A,
                0,
                requested_n);

        if( err )
        {
                fprintf(
                        stderr,
                        "%s: lis_matrix_set_size failed: %d\n",
                        name,
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        err = lis_matrix_get_size(
                A,
                &local_n,
                &global_n);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_matrix_get_range(
                A,
                &is,
                &ie);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        (void)ie;

        ctx.n  = local_n;
        ctx.gn = global_n;
        ctx.is = is;

        ctx.diag    = (LIS_SCALAR)4.0;
        ctx.offdiag = (LIS_SCALAR)-1.0;

        reset_counters(&ctx);

        err = lis_matrix_set_user(
                A,
                &ctx,
                user_matvec,
                NULL);

        if( err )
        {
                fprintf(
                        stderr,
                        "%s: lis_matrix_set_user failed: %d\n",
                        name,
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        if( use_diagonal )
        {
                err = lis_matrix_set_user_diagonal(
                        A,
                        user_get_diagonal);

                if( err )
                {
                        fprintf(
                                stderr,
                                "%s: lis_matrix_set_user_diagonal failed: %d\n",
                                name,
                                (int)err);

                        failed = 1;
                        goto cleanup;
                }
        }

        if( use_row )
        {
                err = lis_matrix_set_user_row(
                        A,
                        user_get_row);

                if( err )
                {
                        fprintf(
                                stderr,
                                "%s: lis_matrix_set_user_row failed: %d\n",
                                name,
                                (int)err);

                        failed = 1;
                        goto cleanup;
                }
        }

        err = lis_vector_duplicate(
                A,
                &b);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &x);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_duplicate(
                A,
                &exact);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_set_all(
                (LIS_SCALAR)1.0,
                exact);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        /*
         * Build b with the same USER operator.
         * The exact solution is x = 1.
         */
        err = lis_matvec(
                A,
                exact,
                b);

        if( err )
        {
                fprintf(
                        stderr,
                        "%s: RHS matvec failed: %d\n",
                        name,
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        reset_counters(&ctx);

        err = lis_vector_set_all(
                (LIS_SCALAR)0.0,
                x);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_solver_create(
                &solver);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_solver_set_option(
                (char *)options,
                solver);

        if( err )
        {
                fprintf(
                        stderr,
                        "%s: lis_solver_set_option failed: %d\n",
                        name,
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        err = lis_solve(
                A,
                b,
                x,
                solver);

        if( !expect_success )
        {
                if( err!=LIS_ERR_NOT_IMPLEMENTED )
                {
                        fprintf(
                                stderr,
                                "%s: returned %d, expected LIS_ERR_NOT_IMPLEMENTED=%d\n",
                                name,
                                (int)err,
                                (int)LIS_ERR_NOT_IMPLEMENTED);

                        failed = 1;
                        goto cleanup;
                }

                if( ctx.row_query_calls!=0 ||
                    ctx.row_fill_calls!=0 )
                {
                        fprintf(
                                stderr,
                                "%s: row callback was reached despite missing ROW_ACCESS\n",
                                name);

                        failed = 1;
                        goto cleanup;
                }

                printf(
                        "CASE=%s PASS expected_err=%d "
                        "matvec=%lu diagonal=%lu "
                        "row_query=%lu row_fill=%lu\n",
                        name,
                        (int)err,
                        ctx.matvec_calls,
                        ctx.diagonal_calls,
                        ctx.row_query_calls,
                        ctx.row_fill_calls);

                goto cleanup;
        }

        if( err )
        {
                fprintf(
                        stderr,
                        "%s: lis_solve failed: %d\n",
                        name,
                        (int)err);

                failed = 1;
                goto cleanup;
        }

        err = lis_solver_get_status(
                solver,
                &status);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        if( status!=LIS_SUCCESS )
        {
                fprintf(
                        stderr,
                        "%s: solver status=%d\n",
                        name,
                        (int)status);

                failed = 1;
                goto cleanup;
        }

        if( ctx.matvec_calls==0 )
        {
                fprintf(
                        stderr,
                        "%s: USER matvec was not used\n",
                        name);

                failed = 1;
                goto cleanup;
        }

        if( require_row_calls )
        {
                if( ctx.row_query_calls==0 ||
                    ctx.row_fill_calls==0 )
                {
                        fprintf(
                                stderr,
                                "%s: USER row callback was not materialized "
                                "(query=%lu fill=%lu)\n",
                                name,
                                ctx.row_query_calls,
                                ctx.row_fill_calls);

                        failed = 1;
                        goto cleanup;
                }
        }
        else
        {
                if( ctx.row_query_calls!=0 ||
                    ctx.row_fill_calls!=0 )
                {
                        fprintf(
                                stderr,
                                "%s: unexpected USER row access "
                                "(query=%lu fill=%lu)\n",
                                name,
                                ctx.row_query_calls,
                                ctx.row_fill_calls);

                        failed = 1;
                        goto cleanup;
                }
        }

        if( require_diagonal_calls &&
            ctx.diagonal_calls==0 )
        {
                fprintf(
                        stderr,
                        "%s: USER diagonal callback was not used\n",
                        name);

                failed = 1;
                goto cleanup;
        }

        err = lis_vector_axpy(
                (LIS_SCALAR)-1.0,
                exact,
                x);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        err = lis_vector_nrm2(
                x,
                &error_norm);

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        if( error_norm>(LIS_REAL)1.0e-8 )
        {
                fprintf(
                        stderr,
                        "%s: solution error %.15e\n",
                        name,
                        (double)error_norm);

                failed = 1;
                goto cleanup;
        }

        printf(
                "CASE=%s PASS error=%.6e "
                "matvec=%lu diagonal=%lu "
                "row_query=%lu row_fill=%lu\n",
                name,
                (double)error_norm,
                ctx.matvec_calls,
                ctx.diagonal_calls,
                ctx.row_query_calls,
                ctx.row_fill_calls);


cleanup:

        if( solver )
        {
                lis_solver_destroy(solver);
        }

        if( exact )
        {
                lis_vector_destroy(exact);
        }

        if( x )
        {
                lis_vector_destroy(x);
        }

        if( b )
        {
                lis_vector_destroy(b);
        }

        if( A )
        {
                lis_matrix_destroy(A);
        }

        return failed;
}


int
main(
        int argc,
        char *argv[])
{
        int failed;

        failed = 0;

        if( lis_initialize(&argc,&argv) )
        {
                return 1;
        }

#ifdef USE_MPI
        {
                int nprocs;

                MPI_Comm_size(
                        LIS_COMM_WORLD,
                        &nprocs);

                if( nprocs!=1 )
                {
                        printf(
                                "LIS_USER_PRECON_GATE SKIPPED: "
                                "serial regression requires one MPI rank\n");

                        lis_finalize();

                        return 0;
                }
        }
#endif

        failed += run_case(
                "USER_ILUT",
                "-i gmres -p ilut "
                "-scale none "
                "-tol 1.0e-12 "
                "-maxiter 200 "
                "-maxiter_noimp 0 "
                "-print none",
                1,
                1,
                1,
                1,
                0);

        failed += run_case(
                "USER_ILU",
                "-i gmres -p ilu "
                "-fill 0 "
                "-scale none "
                "-tol 1.0e-12 "
                "-maxiter 200 "
                "-maxiter_noimp 0 "
                "-print none",
                0,
                1,
                1,
                1,
                0);

        failed += run_case(
                "USER_ILUC",
                "-i gmres -p iluc "
                "-iluc_drop 0 "
                "-iluc_rate 5 "
                "-scale none "
                "-tol 1.0e-12 "
                "-maxiter 200 "
                "-maxiter_noimp 0 "
                "-print none",
                0,
                1,
                1,
                1,
                0);

        failed += run_case(
                "USER_JACOBI",
                "-i gmres -p jacobi "
                "-scale none "
                "-tol 1.0e-12 "
                "-maxiter 200 "
                "-maxiter_noimp 0 "
                "-print none",
                1,
                0,
                1,
                0,
                1);

        /*
         * Hybrid is a composition container.  A USER matrix is valid when
         * both its inner Krylov solver and inner preconditioner can consume
         * the capabilities provided by the USER matrix.
         */
        failed += run_case(
                "USER_HYBRID_ILUT",
                "-i gmres -p hybrid "
                "-hybrid_i gmres "
                "-hybrid_p ilut "
                "-hybrid_maxiter 50 "
                "-hybrid_restart 20 "
                "-hybrid_tol 1.0e-10 "
                "-restart 20 "
                "-scale none "
                "-tol 1.0e-12 "
                "-maxiter 200 "
                "-maxiter_noimp 0 "
                "-print none",
                1,
                1,
                1,
                1,
                0);

        /*
         * The same composition must fail before preconditioner setup when
         * the USER matrix does not provide row access required by ILUT.
         */
        failed += run_case(
                "USER_HYBRID_ILUT_NO_ROW",
                "-i gmres -p hybrid "
                "-hybrid_i gmres "
                "-hybrid_p ilut "
                "-hybrid_maxiter 50 "
                "-hybrid_restart 20 "
                "-hybrid_tol 1.0e-10 "
                "-restart 20 "
                "-scale none "
                "-tol 1.0e-12 "
                "-maxiter 200 "
                "-maxiter_noimp 0 "
                "-print none",
                1,
                0,
                0,
                0,
                0);

        /*
         * Validate inner solver requirements as well as inner PC
         * requirements.  BiCG requires the USER adjoint matvec callback.
         */
        failed += run_case(
                "USER_HYBRID_BICG_NO_MATVECH",
                "-i gmres -p hybrid "
                "-hybrid_i bicg "
                "-hybrid_p none "
                "-hybrid_maxiter 20 "
                "-hybrid_tol 1.0e-10 "
                "-restart 20 "
                "-scale none "
                "-tol 1.0e-12 "
                "-maxiter 200 "
                "-maxiter_noimp 0 "
                "-print none",
                0,
                0,
                0,
                0,
                0);

        failed += run_case(
                "USER_ILUT_NO_ROW",
                "-i gmres -p ilut "
                "-scale none "
                "-tol 1.0e-12 "
                "-maxiter 200 "
                "-maxiter_noimp 0 "
                "-print none",
                1,
                0,
                0,
                0,
                0);

        if( failed )
        {
                fprintf(
                        stderr,
                        "LIS_USER_PRECON_GATE FAILED cases=%d\n",
                        failed);

                lis_finalize();

                return 1;
        }

        printf(
                "LIS_USER_PRECON_GATE PASSED\n");

        lis_finalize();

        return 0;
}
