#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PARAM_MAGIC 0x36435055
#define TEST_N 4

static const LIS_REAL tol0  = 1.0e-6;
static const LIS_REAL tol1  = 1.0e-10;
static const LIS_REAL tol2  = 1.0e-4;

static const LIS_REAL drop0 = 0.125;
static const LIS_REAL drop1 = 0.250;
static const LIS_REAL drop2 = 0.375;

typedef struct
{
    int magic;
    int updates;

    LIS_SOLVER solver;
    LIS_PRECON precon;

    LIS_REAL create_tol;
    LIS_REAL create_drop;

    LIS_REAL update_tol;
    LIS_REAL update_drop;
} PARAM_CTX;


static PARAM_CTX *live_ctx = NULL;

static int create_count  = 0;
static int update_count  = 0;
static int destroy_count = 0;
static int psolve_count  = 0;


static LIS_REAL
read_param(LIS_SOLVER solver, LIS_INT param)
{
    return fabs(
        solver->params[
            param-LIS_OPTIONS_LEN]);
}


static int
real_close(LIS_REAL a, LIS_REAL b)
{
    LIS_REAL diff;
    LIS_REAL scale;

    diff  = fabs(a-b);
    scale = fabs(a) + fabs(b) + 1.0;

    return diff <= 1.0e-12*scale;
}


static LIS_INT
create_test_matrix(LIS_MATRIX *A)
{
    LIS_INT err;
    LIS_INT i,is,ie;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,TEST_N);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            (LIS_SCALAR)(2.0+i),
            *A);

        if( err ) return err;
    }

    err = lis_matrix_set_type(
        *A,
        LIS_MATRIX_CSR);

    if( err ) return err;

    return lis_matrix_assemble(*A);
}


/* ------------------------------------------------------------ */
/* Registered USERDEF preconditioner                             */
/* ------------------------------------------------------------ */

static LIS_INT
ordinary_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    (void)solver;
    (void)precon;

    return LIS_SUCCESS;
}


static LIS_INT
param_psd_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    PARAM_CTX *ctx;
    LIS_INT err;

    ctx = (PARAM_CTX *)malloc(
        sizeof(PARAM_CTX));

    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic       = PARAM_MAGIC;
    ctx->updates     = 0;
    ctx->solver      = solver;
    ctx->precon      = precon;

    ctx->create_tol =
        read_param(
            solver,
            LIS_PARAMS_RESID);

    ctx->create_drop =
        read_param(
            solver,
            LIS_PARAMS_DROP);

    ctx->update_tol  = -1.0;
    ctx->update_drop = -1.0;

    err = lis_precon_set_user_data(
        precon,
        ctx);

    if( err )
    {
        free(ctx);
        return err;
    }

    live_ctx = ctx;
    create_count++;

    return LIS_SUCCESS;
}


static LIS_INT
param_psd_update(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    PARAM_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err ) return err;

    ctx = (PARAM_CTX *)user_data;

    if( ctx==NULL ||
        ctx!=live_ctx ||
        ctx->magic!=PARAM_MAGIC ||
        ctx->solver!=solver ||
        ctx->precon!=precon )
    {
        return LIS_FAILS;
    }

    ctx->update_tol =
        read_param(
            solver,
            LIS_PARAMS_RESID);

    ctx->update_drop =
        read_param(
            solver,
            LIS_PARAMS_DROP);

    ctx->updates++;
    update_count++;

    return LIS_SUCCESS;
}


static LIS_INT
param_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    if( live_ctx==NULL ||
        live_ctx->magic!=PARAM_MAGIC ||
        live_ctx->solver!=solver )
    {
        return LIS_FAILS;
    }

    psolve_count++;

    return lis_vector_copy(b,x);
}


static LIS_INT
param_psolveh(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    return param_psolve(
        solver,
        b,
        x);
}


static LIS_INT
param_destroy(
    LIS_PRECON precon)
{
    PARAM_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;
    LIS_INT result;

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err ) return err;

    ctx = (PARAM_CTX *)user_data;

    if( ctx==NULL )
    {
        return LIS_SUCCESS;
    }

    result = LIS_SUCCESS;

    if( ctx!=live_ctx ||
        ctx->magic!=PARAM_MAGIC ||
        ctx->precon!=precon )
    {
        result = LIS_FAILS;
    }

    ctx->magic = 0;
    free(ctx);

    live_ctx = NULL;
    destroy_count++;

    err = lis_precon_set_user_data(
        precon,
        NULL);

    if( err ) return err;

    return result;
}


/* ------------------------------------------------------------ */
/* Solve contract helper                                        */
/* ------------------------------------------------------------ */

static int
run_solve_and_check_tol(
    const char *label,
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_SOLVER solver,
    LIS_PRECON precon,
    LIS_REAL expected_tol)
{
    LIS_INT err;
    LIS_INT status;
    int psolve_before;

    err = lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    if( err )
    {
        fprintf(
            stderr,
            "%s: reset x failed: %d\n",
            label,
            (int)err);

        return 1;
    }

    psolve_before = psolve_count;

    err = lis_solve_kernel(
        A,
        b,
        x,
        solver,
        precon);

    if( err )
    {
        fprintf(
            stderr,
            "%s: solve_kernel returned %d\n",
            label,
            (int)err);

        return 1;
    }

    err = lis_solver_get_status(
        solver,
        &status);

    if( err )
    {
        fprintf(
            stderr,
            "%s: get_status failed: %d\n",
            label,
            (int)err);

        return 1;
    }

    if( status!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "%s: solver status=%d\n",
            label,
            (int)status);

        return 1;
    }

    /*
     * With nrm2_r convergence, lis_solver_get_initial_residual()
     * copies the current LIS_PARAMS_RESID value into solver->tol.
     * This therefore verifies that the tolerance is refreshed for
     * each solve rather than retained from an earlier solve.
     */
    if( !real_close(
            solver->tol,
            expected_tol) )
    {
        fprintf(
            stderr,
            "%s: solver tol=%.17g expected %.17g\n",
            label,
            (double)solver->tol,
            (double)expected_tol);

        return 1;
    }

    if( psolve_count<=psolve_before )
    {
        fprintf(
            stderr,
            "%s: persistent preconditioner was not used\n",
            label);

        return 1;
    }

    return 0;
}


/* ------------------------------------------------------------ */
/* Main regression                                              */
/* ------------------------------------------------------------ */

int
main(int argc, char **argv)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR exact = NULL;

    LIS_SOLVER solver = NULL;
    LIS_PRECON precon = NULL;

    PARAM_CTX *ctx0 = NULL;
    PARAM_CTX *ctx = NULL;
    void *user_data = NULL;

    LIS_INT err;

    int failed = 0;
    int precon_alive = 0;
    int initialized = 0;

    err = lis_initialize(
        &argc,
        &argv);

    if( err )
    {
        return 1;
    }

    initialized = 1;


    /* -------------------------------------------------------- */
    /* Register reusable USERDEF / PSD preconditioner           */
    /* -------------------------------------------------------- */

    err = lis_precon_register_ex(
        "param6c",
        ordinary_create,
        param_psolve,
        param_psolveh,
        param_destroy);

    if( err )
    {
        fprintf(
            stderr,
            "lis_precon_register_ex failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_register_psd(
        "param6c",
        param_psd_create,
        param_psd_update);

    if( err )
    {
        fprintf(
            stderr,
            "lis_precon_register_psd failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Matrix and vectors                                       */
    /* -------------------------------------------------------- */

    err = create_test_matrix(&A);
    if( err )
    {
        fprintf(
            stderr,
            "create_test_matrix failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_vector_duplicate(A,&b);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_duplicate(A,&x);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_duplicate(A,&exact);
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

    err = lis_matvec(
        A,
        exact,
        b);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Persistent solver                                        */
    /* -------------------------------------------------------- */

    err = lis_solver_create(&solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_option(
        "-i cg "
        "-p param6c "
        "-scale none "
        "-print none "
        "-conv_cond nrm2_r "
        "-maxiter 50 "
        "-tol 1.0e-6 "
        "-iluc_drop 0.125",
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "initial solver options failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(
        A,
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* CREATE CONTRACT                                          */
    /* ======================================================== */

    err = lis_precon_psd_create(
        solver,
        &precon);

    if( err )
    {
        fprintf(
            stderr,
            "PSD create failed: %d\n",
            (int)err);

        precon = NULL;
        failed = 1;
        goto cleanup;
    }

    precon_alive = 1;

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    ctx0 = (PARAM_CTX *)user_data;

    if( ctx0==NULL ||
        ctx0!=live_ctx ||
        create_count!=1 ||
        update_count!=0 ||
        destroy_count!=0 ||
        !real_close(ctx0->create_tol,tol0) ||
        !real_close(ctx0->create_drop,drop0) )
    {
        fprintf(
            stderr,
            "create contract mismatch "
            "create=%d update=%d destroy=%d "
            "tol=%.17g drop=%.17g\n",
            create_count,
            update_count,
            destroy_count,
            ctx0 ? (double)ctx0->create_tol : -1.0,
            ctx0 ? (double)ctx0->create_drop : -1.0);

        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* SOLVE #1: initial tolerance                              */
    /* ======================================================== */

    if( run_solve_and_check_tol(
            "solve0",
            A,
            b,
            x,
            solver,
            precon,
            tol0) )
    {
        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=0 ||
        destroy_count!=0 )
    {
        fprintf(
            stderr,
            "solve0 changed preconditioner lifecycle\n");

        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* UPDATE: change solver tolerance and preconditioner param  */
    /* ======================================================== */

    err = lis_solver_set_option(
        "-tol 1.0e-10 "
        "-iluc_drop 0.250",
        solver);

    if( err )
    {
        fprintf(
            stderr,
            "parameter update failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( !real_close(
            read_param(
                solver,
                LIS_PARAMS_RESID),
            tol1) ||
        !real_close(
            read_param(
                solver,
                LIS_PARAMS_DROP),
            drop1) )
    {
        fprintf(
            stderr,
            "updated solver parameters not stored\n");

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_update(
        solver,
        precon);

    if( err )
    {
        /*
         * Existing PSD contract: update failure destroys the
         * caller-owned preconditioner.
         */
        precon = NULL;
        precon_alive = 0;

        fprintf(
            stderr,
            "PSD update failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    ctx = (PARAM_CTX *)user_data;

    if( ctx!=ctx0 ||
        ctx!=live_ctx ||
        create_count!=1 ||
        update_count!=1 ||
        destroy_count!=0 ||
        ctx->updates!=1 ||
        !real_close(ctx->update_tol,tol1) ||
        !real_close(ctx->update_drop,drop1) )
    {
        fprintf(
            stderr,
            "PSD update did not observe current parameters "
            "create=%d update=%d destroy=%d "
            "tol=%.17g drop=%.17g\n",
            create_count,
            update_count,
            destroy_count,
            ctx ? (double)ctx->update_tol : -1.0,
            ctx ? (double)ctx->update_drop : -1.0);

        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* SOLVE #2: updated tolerance                              */
    /* ======================================================== */

    if( run_solve_and_check_tol(
            "solve1",
            A,
            b,
            x,
            solver,
            precon,
            tol1) )
    {
        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=1 ||
        destroy_count!=0 )
    {
        fprintf(
            stderr,
            "solve1 changed preconditioner lifecycle\n");

        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* SOLVE #3: tolerance-only change                          */
    /*                                                         */
    /* No PSD update is required when only a solve-time         */
    /* tolerance changes.                                      */
    /* ======================================================== */

    err = lis_solver_set_option(
        "-tol 1.0e-4",
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    if( run_solve_and_check_tol(
            "solve2",
            A,
            b,
            x,
            solver,
            precon,
            tol2) )
    {
        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=1 ||
        destroy_count!=0 ||
        ctx0->updates!=1 ||
        !real_close(ctx0->update_drop,drop1) )
    {
        fprintf(
            stderr,
            "tolerance-only solve unexpectedly updated "
            "the persistent preconditioner\n");

        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* SECOND PSD UPDATE: preconditioner parameter only          */
    /* ======================================================== */

    err = lis_solver_set_option(
        "-iluc_drop 0.375",
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    /*
     * Changing a preconditioner parameter does not implicitly
     * update a caller-owned persistent preconditioner.  Its
     * existing context must remain unchanged until the explicit
     * PSD update call below.
     */
    if( update_count!=1 ||
        !real_close(ctx0->update_drop,drop1) )
    {
        fprintf(
            stderr,
            "preconditioner changed before explicit PSD update\n");

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_update(
        solver,
        precon);

    if( err )
    {
        precon = NULL;
        precon_alive = 0;

        fprintf(
            stderr,
            "second PSD update failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=2 ||
        destroy_count!=0 ||
        ctx0->updates!=2 ||
        !real_close(ctx0->update_tol,tol2) ||
        !real_close(ctx0->update_drop,drop2) )
    {
        fprintf(
            stderr,
            "second PSD update contract mismatch "
            "create=%d update=%d destroy=%d "
            "tol=%.17g drop=%.17g\n",
            create_count,
            update_count,
            destroy_count,
            (double)ctx0->update_tol,
            (double)ctx0->update_drop);

        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* SOLVE #4: same object remains usable after second update  */
    /* ======================================================== */

    if( run_solve_and_check_tol(
            "solve3",
            A,
            b,
            x,
            solver,
            precon,
            tol2) )
    {
        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=2 ||
        destroy_count!=0 )
    {
        fprintf(
            stderr,
            "solve3 changed preconditioner lifecycle\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Explicit cleanup contract                                */
    /* -------------------------------------------------------- */

    err = lis_precon_destroy(precon);
    precon = NULL;
    precon_alive = 0;

    if( err )
    {
        fprintf(
            stderr,
            "preconditioner destroy failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=2 ||
        destroy_count!=1 )
    {
        fprintf(
            stderr,
            "final lifecycle counts "
            "create=%d update=%d destroy=%d\n",
            create_count,
            update_count,
            destroy_count);

        failed = 1;
        goto cleanup;
    }


cleanup:

    if( precon_alive && precon )
    {
        lis_precon_destroy(precon);
        precon = NULL;
        precon_alive = 0;
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

    if( A )
    {
        lis_matrix_destroy(A);
        A = NULL;
    }

    lis_precon_register_free();

    if( initialized )
    {
        lis_finalize();
    }

    if( failed )
    {
        fprintf(
            stderr,
            "LIS_PARAM_UPDATE_STAGE6C FAILED\n");

        return 1;
    }

    printf(
        "LIS_PARAM_UPDATE_STAGE6C PASSED "
        "create=%d update=%d destroy=%d psolve=%d\n",
        create_count,
        update_count,
        destroy_count,
        psolve_count);

    return 0;
}
