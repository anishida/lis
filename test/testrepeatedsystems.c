#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"

#include <stdio.h>
#include <stdlib.h>

#define REPEAT_MAGIC 0x34454452

typedef struct
{
    int magic;
    int updates;
    LIS_MATRIX matrix;
} REPEAT_CTX;

static REPEAT_CTX *live_ctx = NULL;

static int psd_create_count = 0;
static int psd_update_count = 0;
static int psd_destroy_count = 0;
static int psolve_count = 0;
static int bad_matrix = 0;


/* ------------------------------------------------------------ */
/* Matrix helpers                                               */
/* ------------------------------------------------------------ */

static LIS_INT
create_diag_matrix(LIS_SCALAR base, LIS_MATRIX *A)
{
    LIS_INT err;
    LIS_INT i,is,ie;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,4);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            base + (LIS_SCALAR)i,
            *A);
        if( err ) return err;
    }

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if( err ) return err;

    return lis_matrix_assemble(*A);
}


static LIS_INT
check_solution(LIS_MATRIX A,
               LIS_VECTOR x,
               LIS_VECTOR exact,
               LIS_REAL *error)
{
    LIS_VECTOR diff = NULL;
    LIS_INT err;

    err = lis_vector_duplicate(A,&diff);
    if( err ) return err;

    err = lis_vector_copy(x,diff);
    if( err )
    {
        lis_vector_destroy(diff);
        return err;
    }

    err = lis_vector_axpy((LIS_SCALAR)-1.0,exact,diff);
    if( err )
    {
        lis_vector_destroy(diff);
        return err;
    }

    err = lis_vector_nrm2(diff,error);

    lis_vector_destroy(diff);

    return err;
}


/* ------------------------------------------------------------ */
/* USERDEF / PSD preconditioner                                 */
/* ------------------------------------------------------------ */

static LIS_INT
ordinary_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    (void)solver;
    (void)precon;

    return LIS_SUCCESS;
}


static LIS_INT
repeat_psd_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    REPEAT_CTX *ctx;
    LIS_INT err;

    ctx = (REPEAT_CTX *)malloc(sizeof(REPEAT_CTX));
    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic = REPEAT_MAGIC;
    ctx->updates = 0;
    ctx->matrix = solver->A;

    err = lis_precon_set_user_data(precon,ctx);
    if( err )
    {
        free(ctx);
        return err;
    }

    live_ctx = ctx;
    psd_create_count++;

    return LIS_SUCCESS;
}


static LIS_INT
repeat_psd_update(LIS_SOLVER solver, LIS_PRECON precon)
{
    REPEAT_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (REPEAT_CTX *)user_data;

    if( ctx==NULL ||
        ctx!=live_ctx ||
        ctx->magic!=REPEAT_MAGIC )
    {
        return LIS_FAILS;
    }

    /*
     * The update must see the matrix installed in the persistent
     * solver before the update call.
     */
    ctx->matrix = solver->A;
    ctx->updates++;

    psd_update_count++;

    return LIS_SUCCESS;
}


static LIS_INT
repeat_psolve(LIS_SOLVER solver,
              LIS_VECTOR b,
              LIS_VECTOR x)
{
    if( live_ctx==NULL ||
        live_ctx->magic!=REPEAT_MAGIC )
    {
        return LIS_FAILS;
    }

    /*
     * A persistent preconditioner must always correspond to the
     * matrix currently used by the solver.
     */
    if( solver->A!=live_ctx->matrix )
    {
        bad_matrix = 1;
        return LIS_FAILS;
    }

    psolve_count++;

    /*
     * Identity action is enough here.  This test checks lifecycle,
     * updates, matrix ownership and warm-start behaviour rather than
     * preconditioner quality.
     */
    return lis_vector_copy(b,x);
}


static LIS_INT
repeat_psolveh(LIS_SOLVER solver,
               LIS_VECTOR b,
               LIS_VECTOR x)
{
    return repeat_psolve(solver,b,x);
}


static LIS_INT
repeat_destroy(LIS_PRECON precon)
{
    REPEAT_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (REPEAT_CTX *)user_data;

    if( ctx==NULL )
    {
        return LIS_SUCCESS;
    }

    if( ctx!=live_ctx ||
        ctx->magic!=REPEAT_MAGIC )
    {
        return LIS_FAILS;
    }

    ctx->magic = 0;

    free(ctx);

    live_ctx = NULL;
    psd_destroy_count++;

    return lis_precon_set_user_data(precon,NULL);
}


/* ------------------------------------------------------------ */
/* Main regression                                              */
/* ------------------------------------------------------------ */

int
main(int argc, char **argv)
{
    LIS_MATRIX A0 = NULL;
    LIS_MATRIX A1 = NULL;
    LIS_MATRIX A2 = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR exact = NULL;

    LIS_SOLVER solver = NULL;
    LIS_PRECON precon = NULL;

    REPEAT_CTX *ctx0 = NULL;
    REPEAT_CTX *ctx = NULL;
    void *user_data = NULL;

    LIS_INT err;
    LIS_INT status;
    LIS_INT iter;

    LIS_REAL error;

    int psolve_before;
    int failed = 0;
    int precon_alive = 0;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;


    /* -------------------------------------------------------- */
    /* Register persistent USERDEF + PSD callbacks              */
    /* -------------------------------------------------------- */

    err = lis_precon_register_ex(
        "repeat4ed",
        ordinary_create,
        repeat_psolve,
        repeat_psolveh,
        repeat_destroy);

    if( err )
    {
        fprintf(stderr,
                "lis_precon_register_ex failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_precon_register_psd(
        "repeat4ed",
        repeat_psd_create,
        repeat_psd_update);

    if( err )
    {
        fprintf(stderr,
                "lis_precon_register_psd failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* A0, A1, A2: same dimensions, changing coefficients      */
    /* -------------------------------------------------------- */

    err = create_diag_matrix((LIS_SCALAR)2.0,&A0);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = create_diag_matrix((LIS_SCALAR)3.0,&A1);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = create_diag_matrix((LIS_SCALAR)4.0,&A2);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(A0,&b);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_duplicate(A0,&x);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_duplicate(A0,&exact);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    /*
     * Exact solution for all three systems.
     *
     * RHS changes together with the matrix:
     *
     *     b_k = A_k * exact
     *
     * This lets A1 verify warm-start deterministically:
     * solution from A0 is already the exact solution of A1.
     */
    err = lis_vector_set_all((LIS_SCALAR)1.0,exact);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_set_all((LIS_SCALAR)0.0,x);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Persistent solver                                       */
    /* -------------------------------------------------------- */

    err = lis_solver_create(&solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_option(
        "-i cg -p repeat4ed "
        "-scale none -print none "
        "-maxiter 50 -tol 1e-12",
        solver);

    if( err )
    {
        fprintf(stderr,
                "lis_solver_set_option failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    /*
     * Reuse the solution supplied in x instead of resetting every
     * repeated system to a zero initial guess.
     */
    solver->options[LIS_OPTIONS_INITGUESS_ZEROS] = LIS_FALSE;


    /* ======================================================== */
    /* SYSTEM A0                                               */
    /* ======================================================== */

    err = lis_matvec(A0,exact,b);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(A0,solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_create(solver,&precon);
    if( err )
    {
        fprintf(stderr,
                "initial PSD create failed: %d\n",
                (int)err);
        precon = NULL;
        failed = 1;
        goto cleanup;
    }

    precon_alive = 1;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    ctx0 = (REPEAT_CTX *)user_data;

    if( ctx0==NULL ||
        ctx0!=live_ctx ||
        ctx0->matrix!=A0 ||
        ctx0->updates!=0 )
    {
        fprintf(stderr,
                "invalid persistent context after create\n");
        failed = 1;
        goto cleanup;
    }

    psolve_before = psolve_count;

    err = lis_solve_kernel(A0,b,x,solver,precon);
    if( err )
    {
        fprintf(stderr,
                "A0 solve_kernel failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_get_status(solver,&status);
    if( err || status!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "A0 status=%d\n",
                (int)status);
        failed = 1;
        goto cleanup;
    }

    err = check_solution(A0,x,exact,&error);
    if( err || error>1.0e-10 )
    {
        fprintf(stderr,
                "A0 solution error=%.15e\n",
                (double)error);
        failed = 1;
        goto cleanup;
    }

    if( psolve_count<=psolve_before )
    {
        fprintf(stderr,
                "A0 did not use persistent preconditioner\n");
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* SYSTEM A1                                               */
    /* ======================================================== */

    /*
     * x still contains the exact A0 solution, which is also the
     * exact A1 solution.  With warm start enabled, A1 should
     * converge immediately from the supplied x.
     */
    err = lis_matvec(A1,exact,b);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(A1,solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_update(solver,precon);
    if( err )
    {
        fprintf(stderr,
                "A1 PSD update failed: %d\n",
                (int)err);

        /*
         * lis_precon_psd_update() owns destruction on update error.
         */
        precon_alive = 0;
        precon = NULL;

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_get_user_data(precon,&user_data);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    ctx = (REPEAT_CTX *)user_data;

    if( ctx!=ctx0 ||
        ctx->matrix!=A1 ||
        ctx->updates!=1 )
    {
        fprintf(stderr,
                "persistent context changed after A1 update\n");
        failed = 1;
        goto cleanup;
    }

    err = lis_solve_kernel(A1,b,x,solver,precon);
    if( err )
    {
        fprintf(stderr,
                "A1 solve_kernel failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_get_status(solver,&status);
    if( err || status!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "A1 status=%d\n",
                (int)status);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_get_iter(solver,&iter);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    /*
     * LIS reports immediate convergence of the initial residual
     * as iteration 1.  If the kernel silently reset x to zero,
     * this nonuniform diagonal system would require more work.
     */
    if( iter!=1 )
    {
        fprintf(stderr,
                "A1 warm start iter=%d, expected 1\n",
                (int)iter);
        failed = 1;
        goto cleanup;
    }

    err = check_solution(A1,x,exact,&error);
    if( err || error>1.0e-10 )
    {
        fprintf(stderr,
                "A1 solution error=%.15e\n",
                (double)error);
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* SYSTEM A2                                               */
    /* ======================================================== */

    err = lis_matvec(A2,exact,b);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    /*
     * Keep a nonzero warm start but perturb the previous solution.
     * A2 therefore performs real iterations and exercises the
     * preconditioner after its second update.
     */
    err = lis_vector_scale((LIS_SCALAR)0.9,x);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(A2,solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_update(solver,precon);
    if( err )
    {
        fprintf(stderr,
                "A2 PSD update failed: %d\n",
                (int)err);

        precon_alive = 0;
        precon = NULL;

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_get_user_data(precon,&user_data);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    ctx = (REPEAT_CTX *)user_data;

    if( ctx!=ctx0 ||
        ctx->matrix!=A2 ||
        ctx->updates!=2 )
    {
        fprintf(stderr,
                "persistent context changed after A2 update\n");
        failed = 1;
        goto cleanup;
    }

    psolve_before = psolve_count;

    err = lis_solve_kernel(A2,b,x,solver,precon);
    if( err )
    {
        fprintf(stderr,
                "A2 solve_kernel failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_get_status(solver,&status);
    if( err || status!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "A2 status=%d\n",
                (int)status);
        failed = 1;
        goto cleanup;
    }

    if( psolve_count<=psolve_before )
    {
        fprintf(stderr,
                "A2 did not use updated persistent preconditioner\n");
        failed = 1;
        goto cleanup;
    }

    err = check_solution(A2,x,exact,&error);
    if( err || error>1.0e-10 )
    {
        fprintf(stderr,
                "A2 solution error=%.15e\n",
                (double)error);
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Final lifecycle checks                                   */
    /* -------------------------------------------------------- */

    if( psd_create_count!=1 ||
        psd_update_count!=2 ||
        psd_destroy_count!=0 ||
        bad_matrix!=0 )
    {
        fprintf(stderr,
                "counts before destroy create=%d update=%d "
                "destroy=%d bad_matrix=%d\n",
                psd_create_count,
                psd_update_count,
                psd_destroy_count,
                bad_matrix);
        failed = 1;
        goto cleanup;
    }


cleanup:

    if( precon_alive && precon )
    {
        err = lis_precon_destroy(precon);

        if( err )
        {
            fprintf(stderr,
                    "lis_precon_destroy failed: %d\n",
                    (int)err);
            failed = 1;
        }

        precon_alive = 0;
        precon = NULL;
    }

    if( solver )
    {
        lis_solver_destroy(solver);
        solver = NULL;
    }

    if( exact )
    {
        lis_vector_destroy(exact);
        exact = NULL;
    }

    if( x )
    {
        lis_vector_destroy(x);
        x = NULL;
    }

    if( b )
    {
        lis_vector_destroy(b);
        b = NULL;
    }

    if( A2 )
    {
        lis_matrix_destroy(A2);
        A2 = NULL;
    }

    if( A1 )
    {
        lis_matrix_destroy(A1);
        A1 = NULL;
    }

    if( A0 )
    {
        lis_matrix_destroy(A0);
        A0 = NULL;
    }

    lis_precon_register_free();
    lis_finalize();

    if( !failed &&
        psd_destroy_count!=1 )
    {
        fprintf(stderr,
                "destroy count=%d, expected 1\n",
                psd_destroy_count);
        failed = 1;
    }

    if( failed )
    {
        fprintf(stderr,
                "LIS_REPEATED_SYSTEMS_STAGE4E FAILED\n");
        return 1;
    }

    printf(
        "LIS_REPEATED_SYSTEMS_STAGE4E PASSED "
        "create=%d update=%d destroy=%d psolve=%d\n",
        psd_create_count,
        psd_update_count,
        psd_destroy_count,
        psolve_count);

    return 0;
}
