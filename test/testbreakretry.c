#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"

#include <stdio.h>
#include <stdlib.h>

#define BREAK_RETRY_MAGIC 0x36445242
#define TEST_N 4

typedef struct
{
    int magic;
    int updates;
    LIS_MATRIX matrix;
} BREAK_RETRY_CTX;

static BREAK_RETRY_CTX *active_ctx = NULL;

static int create_count = 0;
static int update_count = 0;
static int destroy_count = 0;
static int psolve_count = 0;
static int bad_matrix = 0;


/* ------------------------------------------------------------ */
/* Matrix helpers                                                */
/* ------------------------------------------------------------ */

static LIS_INT
create_zero_matrix(LIS_MATRIX *A)
{
    LIS_INT i,is,ie,err;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,TEST_N);
    if( err ) return err;

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    /*
     * Keep explicit zero diagonal entries so the matrix is a
     * normally assembled CSR operator while A*p remains zero.
     * CG must therefore encounter dot(p,A*p)==0.
     */
    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            (LIS_SCALAR)0.0,
            *A);

        if( err ) return err;
    }

    return lis_matrix_assemble(*A);
}


static LIS_INT
create_good_matrix(LIS_MATRIX *A)
{
    LIS_INT i,is,ie,err;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,TEST_N);
    if( err ) return err;

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
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

    return lis_matrix_assemble(*A);
}


/* ------------------------------------------------------------ */
/* USERDEF / PSD preconditioner                                  */
/* ------------------------------------------------------------ */

static LIS_INT
ordinary_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    (void)solver;
    (void)precon;

    return LIS_SUCCESS;
}


static LIS_INT
break_retry_psd_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    BREAK_RETRY_CTX *ctx;
    LIS_INT err;

    ctx = (BREAK_RETRY_CTX *)malloc(
        sizeof(BREAK_RETRY_CTX));

    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic = BREAK_RETRY_MAGIC;
    ctx->updates = 0;
    ctx->matrix = solver->A;

    err = lis_precon_set_user_data(
        precon,
        ctx);

    if( err )
    {
        free(ctx);
        return err;
    }

    active_ctx = ctx;
    create_count++;

    return LIS_SUCCESS;
}


static LIS_INT
break_retry_psd_update(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    BREAK_RETRY_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err ) return err;

    ctx = (BREAK_RETRY_CTX *)user_data;

    if( ctx==NULL ||
        ctx!=active_ctx ||
        ctx->magic!=BREAK_RETRY_MAGIC )
    {
        return LIS_FAILS;
    }

    /*
     * The explicit recovery update must observe the new system
     * installed in the solver after the breakdown.
     */
    ctx->matrix = solver->A;
    ctx->updates++;

    update_count++;

    return LIS_SUCCESS;
}


static LIS_INT
break_retry_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    BREAK_RETRY_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    if( solver->precon==NULL )
    {
        return LIS_FAILS;
    }

    err = lis_precon_get_user_data(
        solver->precon,
        &user_data);

    if( err ) return err;

    ctx = (BREAK_RETRY_CTX *)user_data;

    if( ctx==NULL ||
        ctx!=active_ctx ||
        ctx->magic!=BREAK_RETRY_MAGIC )
    {
        return LIS_FAILS;
    }

    if( solver->A!=ctx->matrix )
    {
        bad_matrix = 1;
        return LIS_FAILS;
    }

    psolve_count++;

    return lis_vector_copy(b,x);
}


static LIS_INT
break_retry_psolveh(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    return break_retry_psolve(
        solver,
        b,
        x);
}


static LIS_INT
break_retry_destroy(
    LIS_PRECON precon)
{
    BREAK_RETRY_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err ) return err;

    ctx = (BREAK_RETRY_CTX *)user_data;

    if( ctx==NULL )
    {
        return LIS_SUCCESS;
    }

    if( ctx!=active_ctx ||
        ctx->magic!=BREAK_RETRY_MAGIC )
    {
        return LIS_FAILS;
    }

    ctx->magic = 0;
    free(ctx);

    active_ctx = NULL;
    destroy_count++;

    return lis_precon_set_user_data(
        precon,
        NULL);
}


/* ------------------------------------------------------------ */
/* Vector / result helpers                                       */
/* ------------------------------------------------------------ */

static LIS_INT
prepare_vectors(
    LIS_MATRIX A,
    LIS_VECTOR *b,
    LIS_VECTOR *x,
    LIS_VECTOR *exact,
    LIS_VECTOR *work)
{
    LIS_INT err;

    err = lis_vector_duplicate(A,b);
    if( err ) return err;

    err = lis_vector_duplicate(A,x);
    if( err ) return err;

    err = lis_vector_duplicate(A,exact);
    if( err ) return err;

    err = lis_vector_duplicate(A,work);
    if( err ) return err;

    err = lis_vector_set_all(
        (LIS_SCALAR)1.0,
        *exact);

    if( err ) return err;

    err = lis_matvec(
        A,
        *exact,
        *b);

    if( err ) return err;

    return lis_vector_set_all(
        (LIS_SCALAR)0.0,
        *x);
}


static LIS_INT
solution_error(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_VECTOR work,
    LIS_REAL *nrm2)
{
    LIS_INT err;

    err = lis_matvec(
        A,
        x,
        work);

    if( err ) return err;

    err = lis_vector_axpy(
        (LIS_SCALAR)-1.0,
        b,
        work);

    if( err ) return err;

    return lis_vector_nrm2(
        work,
        nrm2);
}


/* ------------------------------------------------------------ */
/* Main contract test                                            */
/* ------------------------------------------------------------ */

int
main(int argc, char **argv)
{
    LIS_MATRIX A_bad = NULL;
    LIS_MATRIX A_good = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR exact = NULL;
    LIS_VECTOR work = NULL;

    LIS_SOLVER solver = NULL;

    LIS_PRECON precon = NULL;
    LIS_PRECON precon0 = NULL;

    BREAK_RETRY_CTX *ctx0 = NULL;
    void *user_data = NULL;

    LIS_INT err;
    LIS_INT kernel_bad;
    LIS_INT status_bad = -999;
    LIS_INT kernel_good;
    LIS_INT status_good = -999;

    LIS_REAL error = -1.0;

    int psolve_after_bad = 0;
    int precon_alive = 0;
    int failed = 0;


    err = lis_initialize(
        &argc,
        &argv);

    if( err )
    {
        return 1;
    }


    /* -------------------------------------------------------- */
    /* Registration                                             */
    /* -------------------------------------------------------- */

    err = lis_precon_register_ex(
        "retry6d",
        ordinary_create,
        break_retry_psolve,
        break_retry_psolveh,
        break_retry_destroy);

    if( err )
    {
        fprintf(stderr,
            "register_ex failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_register_psd(
        "retry6d",
        break_retry_psd_create,
        break_retry_psd_update);

    if( err )
    {
        fprintf(stderr,
            "register_psd failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Matrices / vectors                                       */
    /* -------------------------------------------------------- */

    err = create_zero_matrix(
        &A_bad);

    if( err )
    {
        fprintf(stderr,
            "create_zero_matrix failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = create_good_matrix(
        &A_good);

    if( err )
    {
        fprintf(stderr,
            "create_good_matrix failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = prepare_vectors(
        A_good,
        &b,
        &x,
        &exact,
        &work);

    if( err )
    {
        fprintf(stderr,
            "prepare_vectors failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Solver                                                   */
    /* -------------------------------------------------------- */

    err = lis_solver_create(
        &solver);

    if( err )
    {
        fprintf(stderr,
            "solver create failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_option(
        "-i cg "
        "-p retry6d "
        "-scale none "
        "-print none "
        "-tol 1.0e-12 "
        "-maxiter 20 "
        "-maxiter_noimp 0",
        solver);

    if( err )
    {
        fprintf(stderr,
            "solver options failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Persistent PSD create                                    */
    /* -------------------------------------------------------- */

    err = lis_solver_set_matrix(
        A_bad,
        solver);

    if( err )
    {
        fprintf(stderr,
            "set A_bad failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_create(
        solver,
        &precon);

    if( err )
    {
        fprintf(stderr,
            "PSD create failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    precon_alive = 1;
    precon0 = precon;

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    ctx0 = (BREAK_RETRY_CTX *)user_data;

    if( ctx0==NULL ||
        ctx0!=active_ctx ||
        ctx0->magic!=BREAK_RETRY_MAGIC ||
        ctx0->matrix!=A_bad ||
        ctx0->updates!=0 ||
        create_count!=1 ||
        update_count!=0 ||
        destroy_count!=0 )
    {
        fprintf(stderr,
            "initial PSD lifecycle mismatch\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* BREAKDOWN                                                */
    /* -------------------------------------------------------- */

    err = lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    kernel_bad = lis_solve_kernel(
        A_bad,
        b,
        x,
        solver,
        precon);

    err = lis_solver_get_status(
        solver,
        &status_bad);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    /*
     * Iterative breakdown is a solver result, not a kernel API
     * failure. The kernel call succeeds and the result is exposed
     * through solver status.
     */
    if( kernel_bad!=LIS_SUCCESS ||
        status_bad!=LIS_BREAKDOWN )
    {
        fprintf(stderr,
            "breakdown result kernel=%d status=%d\n",
            (int)kernel_bad,
            (int)status_bad);

        failed = 1;
        goto cleanup;
    }

    psolve_after_bad = psolve_count;

    if( psolve_after_bad<=0 )
    {
        fprintf(stderr,
            "breakdown solve did not use preconditioner\n");

        failed = 1;
        goto cleanup;
    }

    /*
     * A solver breakdown must not destroy or replace a
     * caller-owned persistent preconditioner.
     */
    if( precon!=precon0 ||
        active_ctx!=ctx0 ||
        ctx0->magic!=BREAK_RETRY_MAGIC ||
        ctx0->matrix!=A_bad ||
        ctx0->updates!=0 ||
        create_count!=1 ||
        update_count!=0 ||
        destroy_count!=0 ||
        bad_matrix!=0 )
    {
        fprintf(stderr,
            "persistent state changed after breakdown\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* RECOVERY UPDATE                                          */
    /* -------------------------------------------------------- */

    err = lis_solver_set_matrix(
        A_good,
        solver);

    if( err )
    {
        fprintf(stderr,
            "set A_good failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_update(
        solver,
        precon);

    if( err )
    {
        /*
         * lis_precon_psd_update() destroys the preconditioner
         * on callback failure. Never reuse that pointer.
         */
        precon_alive = 0;
        precon = NULL;

        fprintf(stderr,
            "recovery PSD update failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( precon!=precon0 ||
        active_ctx!=ctx0 ||
        ctx0->magic!=BREAK_RETRY_MAGIC ||
        ctx0->matrix!=A_good ||
        ctx0->updates!=1 ||
        create_count!=1 ||
        update_count!=1 ||
        destroy_count!=0 ||
        bad_matrix!=0 )
    {
        fprintf(stderr,
            "recovery update lifecycle mismatch\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* RETRY                                                    */
    /* -------------------------------------------------------- */

    err = lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    kernel_good = lis_solve_kernel(
        A_good,
        b,
        x,
        solver,
        precon);

    err = lis_solver_get_status(
        solver,
        &status_good);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = solution_error(
        A_good,
        b,
        x,
        work,
        &error);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    if( kernel_good!=LIS_SUCCESS ||
        status_good!=LIS_SUCCESS ||
        error>1.0e-10 )
    {
        fprintf(stderr,
            "retry result kernel=%d status=%d error=%.15e\n",
            (int)kernel_good,
            (int)status_good,
            (double)error);

        failed = 1;
        goto cleanup;
    }

    if( psolve_count<=psolve_after_bad )
    {
        fprintf(stderr,
            "retry did not use persistent preconditioner\n");

        failed = 1;
        goto cleanup;
    }

    if( precon!=precon0 ||
        active_ctx!=ctx0 ||
        ctx0->magic!=BREAK_RETRY_MAGIC ||
        ctx0->matrix!=A_good ||
        ctx0->updates!=1 ||
        create_count!=1 ||
        update_count!=1 ||
        destroy_count!=0 ||
        bad_matrix!=0 )
    {
        fprintf(stderr,
            "persistent state changed during retry\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* FINAL DESTROY                                            */
    /* -------------------------------------------------------- */

    err = lis_precon_destroy(
        precon);

    precon_alive = 0;
    precon = NULL;

    if( err )
    {
        fprintf(stderr,
            "final destroy failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=1 ||
        destroy_count!=1 ||
        active_ctx!=NULL ||
        bad_matrix!=0 )
    {
        fprintf(stderr,
            "final lifecycle create=%d update=%d destroy=%d "
            "psolve=%d bad_matrix=%d\n",
            create_count,
            update_count,
            destroy_count,
            psolve_count,
            bad_matrix);

        failed = 1;
        goto cleanup;
    }


cleanup:

    if( precon_alive && precon )
    {
        lis_precon_destroy(precon);
        precon = NULL;
    }

    if( solver )
        lis_solver_destroy(solver);

    if( work )
        lis_vector_destroy(work);

    if( exact )
        lis_vector_destroy(exact);

    if( x )
        lis_vector_destroy(x);

    if( b )
        lis_vector_destroy(b);

    if( A_good )
        lis_matrix_destroy(A_good);

    if( A_bad )
        lis_matrix_destroy(A_bad);

    lis_precon_register_free();
    lis_finalize();

    if( failed )
    {
        fprintf(stderr,
            "LIS_BREAK_RETRY_STAGE6D FAILED\n");

        return 1;
    }

    printf(
        "LIS_BREAK_RETRY_STAGE6D PASSED "
        "create=%d update=%d destroy=%d psolve=%d\n",
        create_count,
        update_count,
        destroy_count,
        psolve_count);

    return 0;
}
