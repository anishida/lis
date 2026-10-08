#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"

#include <stdio.h>
#include <stdlib.h>

#define SYSTEM6E_MAGIC 0x36455359
#define SYSTEM6E_N 24


typedef struct
{
    int magic;
    int updates;
    LIS_MATRIX matrix;
} SYSTEM6E_CTX;


static SYSTEM6E_CTX *active_ctx = NULL;

static int create_count = 0;
static int update_count = 0;
static int destroy_count = 0;
static int psolve_count = 0;
static int bad_matrix = 0;


/* ------------------------------------------------------------ */
/* Matrix helpers                                                */
/* ------------------------------------------------------------ */

static LIS_INT
create_tridiag(
    LIS_SCALAR diag,
    LIS_SCALAR offdiag,
    LIS_MATRIX *A)
{
    LIS_INT i,is,ie,err;

    err = lis_matrix_create(
        LIS_COMM_WORLD,
        A);

    if( err ) return err;

    err = lis_matrix_set_size(
        *A,
        0,
        SYSTEM6E_N);

    if( err ) return err;

    err = lis_matrix_set_type(
        *A,
        LIS_MATRIX_CSR);

    if( err ) return err;

    err = lis_matrix_get_range(
        *A,
        &is,
        &ie);

    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        if( i>0 )
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i-1,
                offdiag,
                *A);

            if( err ) return err;
        }

        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            diag,
            *A);

        if( err ) return err;

        if( i<SYSTEM6E_N-1 )
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i+1,
                offdiag,
                *A);

            if( err ) return err;
        }
    }

    return lis_matrix_assemble(*A);
}


/* ------------------------------------------------------------ */
/* Changing exact solutions                                      */
/* ------------------------------------------------------------ */

static LIS_SCALAR
exact_value(
    LIS_INT i,
    int mode)
{
    LIS_REAL x;

    if( mode==0 )
    {
        x = 0.75
          + 0.025*(LIS_REAL)(i+1)
          + 0.004*(LIS_REAL)((i*7)%11);
    }
    else if( mode==1 )
    {
        x = 0.35
          + 0.017*(LIS_REAL)((i+1)*(i+1))
          - 0.021*(LIS_REAL)((i*5)%9);
    }
    else
    {
        x = 1.10
          - 0.019*(LIS_REAL)(i+1)
          + 0.031*(LIS_REAL)((i*3)%7);
    }

    return (LIS_SCALAR)x;
}


static LIS_INT
set_exact(
    LIS_MATRIX A,
    LIS_VECTOR exact,
    int mode)
{
    LIS_INT i,is,ie,err;

    err = lis_vector_set_all(
        (LIS_SCALAR)0.0,
        exact);

    if( err ) return err;

    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        err = lis_vector_set_value(
            LIS_INS_VALUE,
            i,
            exact_value(i,mode),
            exact);

        if( err ) return err;
    }

    return LIS_SUCCESS;
}


/* ------------------------------------------------------------ */
/* USERDEF persistent PSD preconditioner                          */
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
system_psd_create(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    SYSTEM6E_CTX *ctx;
    LIS_INT err;

    ctx = (SYSTEM6E_CTX *)malloc(
        sizeof(SYSTEM6E_CTX));

    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic = SYSTEM6E_MAGIC;
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
system_psd_update(
    LIS_SOLVER solver,
    LIS_PRECON precon)
{
    SYSTEM6E_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err ) return err;

    ctx = (SYSTEM6E_CTX *)user_data;

    if( ctx==NULL ||
        ctx!=active_ctx ||
        ctx->magic!=SYSTEM6E_MAGIC )
    {
        return LIS_FAILS;
    }

    ctx->matrix = solver->A;
    ctx->updates++;

    update_count++;

    return LIS_SUCCESS;
}


static LIS_INT
system_psolve(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    SYSTEM6E_CTX *ctx;
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

    ctx = (SYSTEM6E_CTX *)user_data;

    if( ctx==NULL ||
        ctx!=active_ctx ||
        ctx->magic!=SYSTEM6E_MAGIC )
    {
        return LIS_FAILS;
    }

    /*
     * Every psolve in the changing-system sequence must observe
     * the matrix installed by the most recent PSD create/update.
     */
    if( solver->A!=ctx->matrix )
    {
        bad_matrix = 1;
        return LIS_FAILS;
    }

    psolve_count++;

    /*
     * Identity preconditioner.  Keeping it as a registered USERDEF
     * preconditioner lets this test exercise persistent PSD
     * ownership without making convergence depend on a particular
     * built-in factorization.
     */
    return lis_vector_copy(
        b,
        x);
}


static LIS_INT
system_psolveh(
    LIS_SOLVER solver,
    LIS_VECTOR b,
    LIS_VECTOR x)
{
    return system_psolve(
        solver,
        b,
        x);
}


static LIS_INT
system_destroy(
    LIS_PRECON precon)
{
    SYSTEM6E_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(
        precon,
        &user_data);

    if( err ) return err;

    ctx = (SYSTEM6E_CTX *)user_data;

    if( ctx==NULL )
    {
        return LIS_SUCCESS;
    }

    if( ctx!=active_ctx ||
        ctx->magic!=SYSTEM6E_MAGIC )
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
/* Result helpers                                                */
/* ------------------------------------------------------------ */

static LIS_INT
solution_error(
    LIS_VECTOR x,
    LIS_VECTOR exact,
    LIS_VECTOR work,
    LIS_REAL *error)
{
    LIS_INT err;

    err = lis_vector_copy(
        x,
        work);

    if( err ) return err;

    err = lis_vector_axpy(
        (LIS_SCALAR)-1.0,
        exact,
        work);

    if( err ) return err;

    return lis_vector_nrm2(
        work,
        error);
}


static LIS_INT
run_system(
    int system_id,
    int exact_mode,
    LIS_REAL requested_tol,
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_VECTOR exact,
    LIS_VECTOR work,
    LIS_SOLVER solver,
    LIS_PRECON precon,
    LIS_PRECON precon0,
    SYSTEM6E_CTX *ctx0,
    LIS_INT expected_updates,
    LIS_INT *iter_out,
    LIS_REAL *resid_out,
    LIS_REAL *error_out)
{
    LIS_INT err;
    LIS_INT kernel;
    LIS_INT status;
    LIS_INT iter;

    LIS_REAL resid;
    LIS_REAL error;

    int psolve_before;


    /*
     * Both the exact solution and RHS change from one system to
     * the next:
     *
     *     b_k = A_k * exact_k
     */
    err = set_exact(
        A,
        exact,
        exact_mode);

    if( err )
    {
        fprintf(stderr,
            "SYSTEM%d set exact failed: %d\n",
            system_id,
            (int)err);

        return err;
    }

    err = lis_matvec(
        A,
        exact,
        b);

    if( err )
    {
        fprintf(stderr,
            "SYSTEM%d RHS creation failed: %d\n",
            system_id,
            (int)err);

        return err;
    }

    /*
     * Use a deterministic zero initial guess for every system.
     * This keeps the tolerance effect independent of warm-start
     * behavior already covered by testrepeatedsystems.
     */
    err = lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    if( err ) return err;

    psolve_before = psolve_count;

    kernel = lis_solve_kernel(
        A,
        b,
        x,
        solver,
        precon);

    err = lis_solver_get_status(
        solver,
        &status);

    if( err ) return err;

    err = lis_solver_get_iter(
        solver,
        &iter);

    if( err ) return err;

    err = lis_solver_get_residualnorm(
        solver,
        &resid);

    if( err ) return err;

    err = solution_error(
        x,
        exact,
        work,
        &error);

    if( err ) return err;


    if( kernel!=LIS_SUCCESS ||
        status!=LIS_SUCCESS )
    {
        fprintf(stderr,
            "SYSTEM%d solve kernel=%d status=%d\n",
            system_id,
            (int)kernel,
            (int)status);

        return LIS_FAILS;
    }


    /*
     * The same caller-owned preconditioner and its user context
     * must survive the entire changing-system sequence.
     */
    if( precon!=precon0 ||
        active_ctx!=ctx0 ||
        active_ctx==NULL ||
        active_ctx->magic!=SYSTEM6E_MAGIC ||
        active_ctx->matrix!=A ||
        active_ctx->updates!=expected_updates ||
        create_count!=1 ||
        update_count!=expected_updates ||
        destroy_count!=0 ||
        bad_matrix!=0 )
    {
        fprintf(stderr,
            "SYSTEM%d persistent lifecycle mismatch\n",
            system_id);

        return LIS_FAILS;
    }


    if( psolve_count<=psolve_before )
    {
        fprintf(stderr,
            "SYSTEM%d did not use persistent preconditioner\n",
            system_id);

        return LIS_FAILS;
    }


    /*
     * Do not inspect solver->params here.  The integration
     * contract is observable behavior: the current solve must
     * satisfy the tolerance most recently supplied through the
     * public option API.
     */
    if( resid >
        (LIS_REAL)(1.05*requested_tol) )
    {
        fprintf(stderr,
            "SYSTEM%d residual %.15e exceeds tolerance %.15e\n",
            system_id,
            (double)resid,
            (double)requested_tol);

        return LIS_FAILS;
    }


    /*
     * These systems are well conditioned.  A loose condition-
     * number allowance verifies the actual solution without
     * coupling the test to exact iteration counts.
     */
    if( error >
        (LIS_REAL)(20.0*requested_tol) )
    {
        fprintf(stderr,
            "SYSTEM%d solution error %.15e exceeds bound %.15e\n",
            system_id,
            (double)error,
            (double)(20.0*requested_tol));

        return LIS_FAILS;
    }


    printf(
        "SYSTEM%d tol=%.3e iter=%d "
        "resid=%.15e error=%.15e "
        "update=%d psolve=%d\n",
        system_id,
        (double)requested_tol,
        (int)iter,
        (double)resid,
        (double)error,
        update_count,
        psolve_count);


    *iter_out = iter;
    *resid_out = resid;
    *error_out = error;

    return LIS_SUCCESS;
}


/* ------------------------------------------------------------ */
/* Main                                                          */
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
    LIS_VECTOR work = NULL;

    LIS_SOLVER solver = NULL;

    LIS_PRECON precon = NULL;
    LIS_PRECON precon0 = NULL;

    SYSTEM6E_CTX *ctx0 = NULL;
    void *user_data = NULL;

    LIS_INT iter0 = -1;
    LIS_INT iter1 = -1;
    LIS_INT iter2 = -1;

    LIS_REAL resid0 = -1.0;
    LIS_REAL resid1 = -1.0;
    LIS_REAL resid2 = -1.0;

    LIS_REAL error0 = -1.0;
    LIS_REAL error1 = -1.0;
    LIS_REAL error2 = -1.0;

    LIS_INT err;

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
    /* Register generic persistent PSD preconditioner            */
    /* -------------------------------------------------------- */

    err = lis_precon_register_ex(
        "system6e",
        ordinary_create,
        system_psolve,
        system_psolveh,
        system_destroy);

    if( err )
    {
        fprintf(stderr,
            "lis_precon_register_ex failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = lis_precon_register_psd(
        "system6e",
        system_psd_create,
        system_psd_update);

    if( err )
    {
        fprintf(stderr,
            "lis_precon_register_psd failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Same dimensions, changing matrix coefficients             */
    /* -------------------------------------------------------- */

    err = create_tridiag(
        (LIS_SCALAR)4.0,
        (LIS_SCALAR)-1.0,
        &A0);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = create_tridiag(
        (LIS_SCALAR)5.0,
        (LIS_SCALAR)-1.25,
        &A1);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = create_tridiag(
        (LIS_SCALAR)6.0,
        (LIS_SCALAR)-1.5,
        &A2);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Reusable vectors                                          */
    /* -------------------------------------------------------- */

    err = lis_vector_duplicate(
        A0,
        &b);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A0,
        &x);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A0,
        &exact);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_duplicate(
        A0,
        &work);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* One solver for all systems                                */
    /* -------------------------------------------------------- */

    err = lis_solver_create(
        &solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_option(
        "-i cg "
        "-p system6e "
        "-scale none "
        "-print none "
        "-maxiter 100 "
        "-maxiter_noimp 0",
        solver);

    if( err )
    {
        fprintf(stderr,
            "initial solver options failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* SYSTEM 0                                                  */
    /* A0 + b0 + exact0 + tol0                                  */
    /* -------------------------------------------------------- */

    err = lis_solver_set_option(
        "-tol 1.0e-2",
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(
        A0,
        solver);

    if( err )
    {
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

    ctx0 = (SYSTEM6E_CTX *)user_data;

    if( ctx0==NULL ||
        ctx0!=active_ctx ||
        ctx0->magic!=SYSTEM6E_MAGIC ||
        ctx0->matrix!=A0 ||
        ctx0->updates!=0 ||
        create_count!=1 ||
        update_count!=0 ||
        destroy_count!=0 )
    {
        fprintf(stderr,
            "initial persistent lifecycle mismatch\n");

        failed = 1;
        goto cleanup;
    }

    err = run_system(
        0,
        0,
        (LIS_REAL)1.0e-2,
        A0,
        b,
        x,
        exact,
        work,
        solver,
        precon,
        precon0,
        ctx0,
        0,
        &iter0,
        &resid0,
        &error0);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* SYSTEM 1                                                  */
    /* new A + new b + new exact + tighter solve tolerance       */
    /* -------------------------------------------------------- */

    err = lis_solver_set_option(
        "-tol 1.0e-6",
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(
        A1,
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_update(
        solver,
        precon);

    if( err )
    {
        /*
         * Failed PSD update destroys the preconditioner.
         * Never reuse the pointer after such an error.
         */
        precon_alive = 0;
        precon = NULL;

        fprintf(stderr,
            "SYSTEM1 PSD update failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = run_system(
        1,
        1,
        (LIS_REAL)1.0e-6,
        A1,
        b,
        x,
        exact,
        work,
        solver,
        precon,
        precon0,
        ctx0,
        1,
        &iter1,
        &resid1,
        &error1);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* SYSTEM 2                                                  */
    /* new A + new b + new exact + tighter solve tolerance       */
    /* -------------------------------------------------------- */

    err = lis_solver_set_option(
        "-tol 1.0e-10",
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(
        A2,
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_precon_psd_update(
        solver,
        precon);

    if( err )
    {
        precon_alive = 0;
        precon = NULL;

        fprintf(stderr,
            "SYSTEM2 PSD update failed: %d\n",
            (int)err);

        failed = 1;
        goto cleanup;
    }

    err = run_system(
        2,
        2,
        (LIS_REAL)1.0e-10,
        A2,
        b,
        x,
        exact,
        work,
        solver,
        precon,
        precon0,
        ctx0,
        2,
        &iter2,
        &resid2,
        &error2);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Integrated sequence checks                                */
    /* -------------------------------------------------------- */

    /*
     * Do not require exact iteration counts.  They are diagnostic
     * only and may vary slightly between supported platforms.
     */
    if( iter0<=0 ||
        iter1<=0 ||
        iter2<=0 )
    {
        fprintf(stderr,
            "invalid iteration counts: %d %d %d\n",
            (int)iter0,
            (int)iter1,
            (int)iter2);

        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=2 ||
        destroy_count!=0 ||
        active_ctx!=ctx0 ||
        ctx0->updates!=2 ||
        ctx0->matrix!=A2 ||
        bad_matrix!=0 )
    {
        fprintf(stderr,
            "pre-destroy lifecycle mismatch\n");

        failed = 1;
        goto cleanup;
    }


    /* -------------------------------------------------------- */
    /* Explicit final destroy                                    */
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
        update_count!=2 ||
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

    if( A2 )
        lis_matrix_destroy(A2);

    if( A1 )
        lis_matrix_destroy(A1);

    if( A0 )
        lis_matrix_destroy(A0);

    lis_precon_register_free();
    lis_finalize();


    if( failed )
    {
        fprintf(stderr,
            "LIS_CHANGING_SYSTEM_STAGE6E FAILED\n");

        return 1;
    }


    printf(
        "LIS_CHANGING_SYSTEM_STAGE6E PASSED "
        "iter=%d/%d/%d "
        "resid=%.3e/%.3e/%.3e "
        "create=%d update=%d destroy=%d psolve=%d\n",
        (int)iter0,
        (int)iter1,
        (int)iter2,
        (double)resid0,
        (double)resid1,
        (double)resid2,
        create_count,
        update_count,
        destroy_count,
        psolve_count);

    return 0;
}
