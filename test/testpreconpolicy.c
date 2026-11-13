#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"

#include <stdio.h>
#include <stdlib.h>

#define POLICY_MAGIC 0x34465043

typedef struct
{
    int magic;
    int generation;
    int updates;
    LIS_MATRIX matrix;
} POLICY_CTX;

static POLICY_CTX *active_ctx = NULL;

static int create_count = 0;
static int update_count = 0;
static int destroy_count = 0;
static int psolve_count = 0;
static int bad_matrix = 0;


/* ------------------------------------------------------------ */
/* Matrix helpers                                               */
/* ------------------------------------------------------------ */

static LIS_INT
create_scaled_identity(LIS_SCALAR scale, LIS_MATRIX *A)
{
    LIS_INT err;
    LIS_INT i,is,ie;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,2);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            scale,
            *A);
        if( err ) return err;
    }

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if( err ) return err;

    return lis_matrix_assemble(*A);
}


static LIS_INT
create_stagnation_matrix(LIS_MATRIX *A)
{
    LIS_INT err;
    LIS_INT i,is,ie;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,2);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        if( i==0 )
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,
                0,
                1,
                (LIS_SCALAR)1.0,
                *A);
            if( err ) return err;
        }
        else if( i==1 )
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,
                1,
                0,
                (LIS_SCALAR)-1.0,
                *A);
            if( err ) return err;
        }
    }

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if( err ) return err;

    return lis_matrix_assemble(*A);
}


static LIS_INT
set_stagnation_rhs(LIS_VECTOR b)
{
    LIS_INT err;
    LIS_INT is,ie;

    err = lis_vector_set_all((LIS_SCALAR)0.0,b);
    if( err ) return err;

    err = lis_vector_get_range(b,&is,&ie);
    if( err ) return err;

    if( is<=0 && 0<ie )
    {
        err = lis_vector_set_value(
            LIS_INS_VALUE,
            0,
            (LIS_SCALAR)1.0,
            b);
        if( err ) return err;
    }

    return LIS_SUCCESS;
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

    err = lis_vector_axpy(
        (LIS_SCALAR)-1.0,
        exact,
        diff);

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
policy_psd_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    POLICY_CTX *ctx;
    LIS_INT err;

    ctx = (POLICY_CTX *)malloc(sizeof(POLICY_CTX));
    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic = POLICY_MAGIC;
    ctx->generation = create_count + 1;
    ctx->updates = 0;
    ctx->matrix = solver->A;

    err = lis_precon_set_user_data(precon,ctx);
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
policy_psd_update(LIS_SOLVER solver, LIS_PRECON precon)
{
    POLICY_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (POLICY_CTX *)user_data;

    if( ctx==NULL ||
        ctx!=active_ctx ||
        ctx->magic!=POLICY_MAGIC )
    {
        return LIS_FAILS;
    }

    ctx->matrix = solver->A;
    ctx->updates++;

    update_count++;

    return LIS_SUCCESS;
}


static LIS_INT
policy_psolve(LIS_SOLVER solver,
              LIS_VECTOR b,
              LIS_VECTOR x)
{
    if( active_ctx==NULL ||
        active_ctx->magic!=POLICY_MAGIC )
    {
        return LIS_FAILS;
    }

    if( solver->A!=active_ctx->matrix )
    {
        bad_matrix = 1;
        return LIS_FAILS;
    }

    psolve_count++;

    return lis_vector_copy(b,x);
}


static LIS_INT
policy_psolveh(LIS_SOLVER solver,
               LIS_VECTOR b,
               LIS_VECTOR x)
{
    return policy_psolve(solver,b,x);
}


static LIS_INT
policy_destroy(LIS_PRECON precon)
{
    POLICY_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (POLICY_CTX *)user_data;

    if( ctx==NULL )
    {
        return LIS_SUCCESS;
    }

    if( ctx!=active_ctx ||
        ctx->magic!=POLICY_MAGIC )
    {
        return LIS_FAILS;
    }

    ctx->magic = 0;

    free(ctx);

    active_ctx = NULL;
    destroy_count++;

    return lis_precon_set_user_data(precon,NULL);
}


/* ------------------------------------------------------------ */
/* Checks                                                       */
/* ------------------------------------------------------------ */

static int
check_context(LIS_PRECON precon,
              LIS_MATRIX matrix,
              int generation,
              int updates,
              const char *label)
{
    POLICY_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err )
    {
        fprintf(stderr,
                "%s: get_user_data failed: %d\n",
                label,
                (int)err);
        return 1;
    }

    ctx = (POLICY_CTX *)user_data;

    if( ctx==NULL ||
        ctx!=active_ctx ||
        ctx->magic!=POLICY_MAGIC ||
        ctx->generation!=generation ||
        ctx->updates!=updates ||
        ctx->matrix!=matrix )
    {
        fprintf(stderr,
                "%s: invalid context "
                "gen=%d updates=%d\n",
                label,
                ctx ? ctx->generation : -1,
                ctx ? ctx->updates : -1);
        return 1;
    }

    return 0;
}


static int
run_good_system(const char *label,
                LIS_MATRIX A,
                LIS_VECTOR b,
                LIS_VECTOR x,
                LIS_VECTOR exact,
                LIS_SOLVER solver,
                LIS_PRECON precon)
{
    LIS_INT err;
    LIS_INT status;
    LIS_REAL error;
    int before;

    err = lis_matvec(A,exact,b);
    if( err )
    {
        fprintf(stderr,
                "%s: lis_matvec failed: %d\n",
                label,
                (int)err);
        return 1;
    }

    err = lis_vector_set_all((LIS_SCALAR)0.0,x);
    if( err ) return 1;

    before = psolve_count;

    err = lis_solve_kernel(
        A,
        b,
        x,
        solver,
        precon);

    if( err!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "%s: solve_kernel returned %d\n",
                label,
                (int)err);
        return 1;
    }

    err = lis_solver_get_status(solver,&status);
    if( err ) return 1;

    if( status!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "%s: solver status=%d\n",
                label,
                (int)status);
        return 1;
    }

    if( psolve_count<=before )
    {
        fprintf(stderr,
                "%s: preconditioner was not used\n",
                label);
        return 1;
    }

    err = check_solution(
        A,
        x,
        exact,
        &error);

    if( err || error>1.0e-10 )
    {
        fprintf(stderr,
                "%s: solution error=%.15e\n",
                label,
                (double)error);
        return 1;
    }

    return 0;
}


/* ------------------------------------------------------------ */
/* Main                                                         */
/* ------------------------------------------------------------ */

int
main(int argc, char **argv)
{
    LIS_MATRIX A0 = NULL;
    LIS_MATRIX A1 = NULL;
    LIS_MATRIX A2 = NULL;
    LIS_MATRIX A_retry = NULL;
    LIS_MATRIX A_stag = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR exact = NULL;

    LIS_SOLVER solver = NULL;
    LIS_PRECON precon = NULL;

    POLICY_CTX *generation1_ctx = NULL;
    POLICY_CTX *generation2_ctx = NULL;

    LIS_INT err;
    LIS_INT status;
    LIS_INT iter;

    int psolve_before;
    int failed = 0;
    int precon_alive = 0;


    err = lis_initialize(&argc,&argv);
    if( err ) return 1;


    err = lis_precon_register_ex(
        "policy4f",
        ordinary_create,
        policy_psolve,
        policy_psolveh,
        policy_destroy);

    if( err )
    {
        fprintf(stderr,
                "register_ex failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    err = lis_precon_register_psd(
        "policy4f",
        policy_psd_create,
        policy_psd_update);

    if( err )
    {
        fprintf(stderr,
                "register_psd failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    err = create_scaled_identity(
        (LIS_SCALAR)2.0,
        &A0);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = create_scaled_identity(
        (LIS_SCALAR)3.0,
        &A1);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = create_scaled_identity(
        (LIS_SCALAR)4.0,
        &A2);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = create_scaled_identity(
        (LIS_SCALAR)5.0,
        &A_retry);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = create_stagnation_matrix(&A_stag);
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

    err = lis_vector_set_all(
        (LIS_SCALAR)1.0,
        exact);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_solver_create(&solver);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_option(
        "-i gmres -p policy4f "
        "-scale none -print none "
        "-restart 1 -maxiter 40 "
        "-maxiter_noimp 0 "
        "-tol 1e-12",
        solver);

    if( err )
    {
        fprintf(stderr,
                "solver options failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* CREATE                                                   */
    /* ======================================================== */

    err = lis_solver_set_matrix(A0,solver);
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
                "CREATE failed: %d\n",
                (int)err);
        precon = NULL;
        failed = 1;
        goto cleanup;
    }

    precon_alive = 1;

    generation1_ctx = active_ctx;

    if( check_context(
            precon,
            A0,
            1,
            0,
            "CREATE") )
    {
        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=0 ||
        destroy_count!=0 )
    {
        fprintf(stderr,
                "CREATE counts=%d/%d/%d\n",
                create_count,
                update_count,
                destroy_count);
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* KEEP                                                     */
    /* ======================================================== */

    if( run_good_system(
            "KEEP-1",
            A0,
            b,
            x,
            exact,
            solver,
            precon) )
    {
        failed = 1;
        goto cleanup;
    }

    if( run_good_system(
            "KEEP-2",
            A0,
            b,
            x,
            exact,
            solver,
            precon) )
    {
        failed = 1;
        goto cleanup;
    }

    if( active_ctx!=generation1_ctx ||
        check_context(
            precon,
            A0,
            1,
            0,
            "KEEP") )
    {
        fprintf(stderr,
                "KEEP changed persistent context\n");
        failed = 1;
        goto cleanup;
    }

    if( create_count!=1 ||
        update_count!=0 ||
        destroy_count!=0 )
    {
        fprintf(stderr,
                "KEEP counts=%d/%d/%d\n",
                create_count,
                update_count,
                destroy_count);
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* UPDATE                                                   */
    /* ======================================================== */

    err = lis_solver_set_matrix(A1,solver);
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
         * update() destroys the preconditioner on failure.
         */
        precon = NULL;
        precon_alive = 0;

        fprintf(stderr,
                "UPDATE failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( active_ctx!=generation1_ctx ||
        check_context(
            precon,
            A1,
            1,
            1,
            "UPDATE") )
    {
        fprintf(stderr,
                "UPDATE replaced context\n");
        failed = 1;
        goto cleanup;
    }

    if( run_good_system(
            "UPDATE",
            A1,
            b,
            x,
            exact,
            solver,
            precon) )
    {
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* REBUILD                                                  */
    /* ======================================================== */

    err = lis_precon_destroy(precon);

    /*
     * lis_precon_destroy() frees the object even when a user
     * destroy callback reports an error, so never reuse pointer.
     */
    precon = NULL;
    precon_alive = 0;

    if( err )
    {
        fprintf(stderr,
                "REBUILD destroy failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    if( destroy_count!=1 ||
        active_ctx!=NULL )
    {
        fprintf(stderr,
                "REBUILD did not destroy generation 1\n");
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(A2,solver);
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
                "REBUILD create failed: %d\n",
                (int)err);
        precon = NULL;
        failed = 1;
        goto cleanup;
    }

    precon_alive = 1;
    generation2_ctx = active_ctx;

    if( check_context(
            precon,
            A2,
            2,
            0,
            "REBUILD") )
    {
        failed = 1;
        goto cleanup;
    }

    if( create_count!=2 ||
        update_count!=1 ||
        destroy_count!=1 )
    {
        fprintf(stderr,
                "REBUILD counts=%d/%d/%d\n",
                create_count,
                update_count,
                destroy_count);
        failed = 1;
        goto cleanup;
    }

    if( run_good_system(
            "REBUILD",
            A2,
            b,
            x,
            exact,
            solver,
            precon) )
    {
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* STAGNATION                                               */
    /* ======================================================== */

    err = lis_solver_set_matrix(
        A_stag,
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
        precon = NULL;
        precon_alive = 0;

        fprintf(stderr,
                "STAGNATION update failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( active_ctx!=generation2_ctx ||
        check_context(
            precon,
            A_stag,
            2,
            1,
            "STAGNATION") )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_option(
        "-maxiter_noimp 1",
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = set_stagnation_rhs(b);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    psolve_before = psolve_count;

    err = lis_solve_kernel(
        A_stag,
        b,
        x,
        solver,
        precon);

    /*
     * lis_solve_kernel() itself succeeds.  The iterative result
     * is reported through solver status.
     */
    if( err!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "STAGNATION solve_kernel returned %d\n",
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

    err = lis_solver_get_iter(
        solver,
        &iter);
    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    if( status!=LIS_MAXITER )
    {
        fprintf(stderr,
                "STAGNATION status=%d, expected %d\n",
                (int)status,
                (int)LIS_MAXITER);
        failed = 1;
        goto cleanup;
    }

    if( iter<=0 || iter>=40 )
    {
        fprintf(stderr,
                "STAGNATION iter=%d\n",
                (int)iter);
        failed = 1;
        goto cleanup;
    }

    if( psolve_count<=psolve_before )
    {
        fprintf(stderr,
                "STAGNATION did not use preconditioner\n");
        failed = 1;
        goto cleanup;
    }

    /*
     * Stagnation must not destroy caller-owned preconditioner.
     */
    if( active_ctx!=generation2_ctx ||
        destroy_count!=1 )
    {
        fprintf(stderr,
                "STAGNATION destroyed persistent preconditioner\n");
        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* RETRY via UPDATE                                         */
    /* ======================================================== */

    err = lis_solver_set_option(
        "-maxiter_noimp 0",
        solver);

    if( err )
    {
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_matrix(
        A_retry,
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
        precon = NULL;
        precon_alive = 0;

        fprintf(stderr,
                "RETRY update failed: %d\n",
                (int)err);

        failed = 1;
        goto cleanup;
    }

    if( active_ctx!=generation2_ctx ||
        check_context(
            precon,
            A_retry,
            2,
            2,
            "RETRY") )
    {
        failed = 1;
        goto cleanup;
    }

    if( run_good_system(
            "RETRY",
            A_retry,
            b,
            x,
            exact,
            solver,
            precon) )
    {
        failed = 1;
        goto cleanup;
    }


    if( create_count!=2 ||
        update_count!=3 ||
        destroy_count!=1 ||
        bad_matrix!=0 )
    {
        fprintf(stderr,
                "counts before final destroy "
                "create=%d update=%d destroy=%d bad_matrix=%d\n",
                create_count,
                update_count,
                destroy_count,
                bad_matrix);

        failed = 1;
        goto cleanup;
    }


    /* ======================================================== */
    /* FINAL DESTROY                                            */
    /* ======================================================== */

    err = lis_precon_destroy(precon);

    precon = NULL;
    precon_alive = 0;

    if( err )
    {
        fprintf(stderr,
                "final destroy failed: %d\n",
                (int)err);
        failed = 1;
        goto cleanup;
    }

    if( destroy_count!=2 ||
        active_ctx!=NULL )
    {
        fprintf(stderr,
                "final destroy count=%d\n",
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

    if( A_stag )
    {
        lis_matrix_destroy(A_stag);
        A_stag = NULL;
    }

    if( A_retry )
    {
        lis_matrix_destroy(A_retry);
        A_retry = NULL;
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


    if( failed )
    {
        fprintf(stderr,
                "LIS_PRECON_POLICY_STAGE4F FAILED\n");
        return 1;
    }

    printf(
        "LIS_PRECON_POLICY_STAGE4F PASSED "
        "create=%d update=%d destroy=%d psolve=%d\n",
        create_count,
        update_count,
        destroy_count,
        psolve_count);

    return 0;
}
