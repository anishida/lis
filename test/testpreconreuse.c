#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"
#include <stdio.h>
#include <stdlib.h>

#define REUSE_MAGIC 0x34445245

#define CTX_ORDINARY 1
#define CTX_PSD      2

typedef struct
{
    int magic;
    int kind;
    int updates;
} REUSE_CTX;

static int ordinary_create_count = 0;
static int ordinary_destroy_count = 0;

static int psd_create_count = 0;
static int psd_update_count = 0;
static int psd_destroy_count = 0;

static int bad_context = 0;

static LIS_MATRIX expected_update_matrix = NULL;

static LIS_INT create_context(LIS_PRECON precon, int kind)
{
    REUSE_CTX *ctx;
    LIS_INT err;

    ctx = (REUSE_CTX *)malloc(sizeof(REUSE_CTX));
    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic = REUSE_MAGIC;
    ctx->kind = kind;
    ctx->updates = 0;

    err = lis_precon_set_user_data(precon,ctx);
    if( err )
    {
        free(ctx);
        return err;
    }

    return LIS_SUCCESS;
}

static LIS_INT ordinary_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    (void)solver;

    ordinary_create_count++;

    return create_context(precon,CTX_ORDINARY);
}

static LIS_INT identity_psolve(LIS_SOLVER solver,
                               LIS_VECTOR b,
                               LIS_VECTOR x)
{
    (void)solver;

    return lis_vector_copy(b,x);
}

static LIS_INT identity_psolveh(LIS_SOLVER solver,
                                LIS_VECTOR b,
                                LIS_VECTOR x)
{
    (void)solver;

    return lis_vector_copy(b,x);
}

static LIS_INT psd_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    (void)solver;

    psd_create_count++;

    return create_context(precon,CTX_PSD);
}

static LIS_INT psd_update(LIS_SOLVER solver, LIS_PRECON precon)
{
    REUSE_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    psd_update_count++;

    if( expected_update_matrix!=NULL &&
        solver->A!=expected_update_matrix )
    {
        fprintf(stderr,
                "PSD update sees the wrong solver matrix\n");
        return LIS_FAILS;
    }

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (REUSE_CTX *)user_data;

    if( ctx==NULL ||
        ctx->magic!=REUSE_MAGIC ||
        ctx->kind!=CTX_PSD )
    {
        return LIS_FAILS;
    }

    ctx->updates++;

    return LIS_SUCCESS;
}

static LIS_INT user_destroy(LIS_PRECON precon)
{
    REUSE_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (REUSE_CTX *)user_data;

    if( ctx==NULL )
    {
        return LIS_SUCCESS;
    }

    if( ctx->magic!=REUSE_MAGIC )
    {
        bad_context = 1;
        return LIS_FAILS;
    }

    if( ctx->kind==CTX_ORDINARY )
    {
        ordinary_destroy_count++;
    }
    else if( ctx->kind==CTX_PSD )
    {
        psd_destroy_count++;
    }
    else
    {
        bad_context = 1;
        return LIS_FAILS;
    }

    ctx->magic = 0;
    free(ctx);

    return lis_precon_set_user_data(precon,NULL);
}

static LIS_INT create_identity_matrix(LIS_MATRIX *A)
{
    LIS_INT i,is,ie,err;
    const LIS_INT n = 2;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,n);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,i,i,(LIS_SCALAR)1.0,*A);
        if( err ) return err;
    }

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if( err ) return err;

    return lis_matrix_assemble(*A);
}

/*
 * Skew-symmetric 2x2 matrix:
 *
 *     [ 0  1 ]
 *     [-1  0 ]
 *
 * With b=e_0 and restarted GMRES(1), the one-dimensional
 * Krylov correction cannot reduce the residual.  This gives a
 * deterministic no-improvement path for -maxiter_noimp.
 */
static LIS_INT create_stagnation_matrix(LIS_MATRIX *A)
{
    LIS_INT i,is,ie,err;
    const LIS_INT n = 2;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if( err ) return err;

    err = lis_matrix_set_size(*A,0,n);
    if( err ) return err;

    err = lis_matrix_get_range(*A,&is,&ie);
    if( err ) return err;

    for(i=is;i<ie;i++)
    {
        if( i==0 )
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,0,1,(LIS_SCALAR)1.0,*A);
            if( err ) return err;
        }
        else if( i==1 )
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,1,0,(LIS_SCALAR)-1.0,*A);
            if( err ) return err;
        }
    }

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if( err ) return err;

    return lis_matrix_assemble(*A);
}

static LIS_INT create_vectors(LIS_MATRIX A,
                              LIS_VECTOR *b,
                              LIS_VECTOR *x)
{
    LIS_INT err,is,ie;

    err = lis_vector_duplicate(A,b);
    if( err ) return err;

    err = lis_vector_duplicate(A,x);
    if( err ) return err;

    err = lis_vector_set_all((LIS_SCALAR)0.0,*b);
    if( err ) return err;

    err = lis_vector_get_range(*b,&is,&ie);
    if( err ) return err;

    /*
     * Use b=e_0.  For the skew-symmetric stagnation matrix,
     *
     *     A*b = -e_1,
     *
     * so restarted GMRES(1) has an exactly orthogonal search
     * direction and cannot reduce the residual.
     *
     * The range check also keeps this valid for MPI.
     */
    if( is<=0 && 0<ie )
    {
        err = lis_vector_set_value(
            LIS_INS_VALUE,0,(LIS_SCALAR)1.0,*b);
        if( err ) return err;
    }

    return lis_vector_set_all((LIS_SCALAR)0.0,*x);
}

static int test_reusable_psd(LIS_MATRIX A_stag,
                             LIS_MATRIX A_good)
{
    LIS_SOLVER solver = NULL;
    LIS_PRECON precon = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    REUSE_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err,status,iter;

    err = create_vectors(A_stag,&b,&x);
    if( err ) return 1;

    err = lis_solver_create(&solver);
    if( err ) return 1;

    err = lis_solver_set_option(
        "-i gmres -p reuse4d "
        "-scale none -print none "
        "-restart 1 -maxiter 20 "
        "-maxiter_noimp 1 -tol 1e-14",
        solver);
    if( err ) return 1;

    err = lis_solver_set_matrix(A_stag,solver);
    if( err ) return 1;

    err = lis_precon_psd_create(solver,&precon);
    if( err )
    {
        fprintf(stderr,
                "PSD create failed: %d\n",(int)err);
        return 1;
    }

    if( psd_create_count!=1 || psd_destroy_count!=0 )
    {
        fprintf(stderr,
                "unexpected initial PSD counts %d/%d\n",
                psd_create_count,psd_destroy_count);
        return 1;
    }

    /*
     * First solve deliberately stagnates.
     *
     * lis_solve_kernel() itself returns LIS_SUCCESS; the
     * iterative-solver result is available through solver status.
     */
    err = lis_solve_kernel(A_stag,b,x,solver,precon);
    if( err!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "stagnating solve_kernel returned %d\n",(int)err);
        return 1;
    }

    err = lis_solver_get_status(solver,&status);
    if( err ) return 1;

    err = lis_solver_get_iter(solver,&iter);
    if( err ) return 1;

    if( status!=LIS_MAXITER )
    {
        fprintf(stderr,
                "stagnating solve status=%d, expected LIS_MAXITER=%d\n",
                (int)status,(int)LIS_MAXITER);
        return 1;
    }

    /*
     * With maxiter_noimp=1 this must stop well before maxiter=20.
     */
    if( iter>=20 )
    {
        fprintf(stderr,
                "maxiter_noimp did not stop the solve early; iter=%d\n",
                (int)iter);
        return 1;
    }

    /*
     * Early stopping must not destroy caller-owned PSD state.
     */
    if( psd_destroy_count!=0 )
    {
        fprintf(stderr,
                "preconditioner destroyed after early stop\n");
        return 1;
    }

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return 1;

    ctx = (REUSE_CTX *)user_data;
    if( ctx==NULL ||
        ctx->magic!=REUSE_MAGIC ||
        ctx->kind!=CTX_PSD ||
        ctx->updates!=0 )
    {
        fprintf(stderr,
                "PSD context invalid after early stop\n");
        return 1;
    }

    /*
     * A new tangent/system matrix is installed before PSD update.
     * The update callback must see A_good.
     */
    expected_update_matrix = A_good;

    err = lis_solver_set_matrix(A_good,solver);
    if( err ) return 1;

    err = lis_precon_psd_update(solver,precon);
    if( err )
    {
        fprintf(stderr,
                "PSD update after early stop failed: %d\n",(int)err);
        return 1;
    }

    expected_update_matrix = NULL;

    if( psd_update_count!=1 || psd_destroy_count!=0 )
    {
        fprintf(stderr,
                "unexpected PSD update counts update=%d destroy=%d\n",
                psd_update_count,psd_destroy_count);
        return 1;
    }

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return 1;

    ctx = (REUSE_CTX *)user_data;
    if( ctx==NULL ||
        ctx->magic!=REUSE_MAGIC ||
        ctx->kind!=CTX_PSD ||
        ctx->updates!=1 )
    {
        fprintf(stderr,
                "PSD context invalid after update\n");
        return 1;
    }

    /*
     * Retry with a well-conditioned matrix using the same
     * preconditioner object.
     */
    lis_vector_set_all((LIS_SCALAR)0.0,x);

    err = lis_solve_kernel(A_good,b,x,solver,precon);
    if( err!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "retry solve_kernel returned %d\n",(int)err);
        return 1;
    }

    err = lis_solver_get_status(solver,&status);
    if( err ) return 1;

    if( status!=LIS_SUCCESS )
    {
        fprintf(stderr,
                "retry status=%d, expected LIS_SUCCESS\n",
                (int)status);
        return 1;
    }

    if( psd_destroy_count!=0 )
    {
        fprintf(stderr,
                "preconditioner destroyed before explicit cleanup\n");
        return 1;
    }

    err = lis_precon_destroy(precon);
    if( err ) return 1;
    precon = NULL;

    if( psd_destroy_count!=1 )
    {
        fprintf(stderr,
                "PSD destroy count=%d, expected 1\n",
                psd_destroy_count);
        return 1;
    }

    lis_solver_destroy(solver);
    lis_vector_destroy(b);
    lis_vector_destroy(x);

    return 0;
}

static int test_lis_solve_error_cleanup(LIS_MATRIX A)
{
    LIS_SOLVER solver = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_INT err;

    err = create_vectors(A,&b,&x);
    if( err ) return 1;

    err = lis_solver_create(&solver);
    if( err ) return 1;

    /*
     * -maxiter_noimp=-1 is rejected inside lis_solve_kernel(),
     * after lis_solve() has already created the preconditioner.
     *
     * lis_solve() owns that preconditioner and must destroy it
     * before returning the kernel error.
     */
    err = lis_solver_set_option(
        "-i gmres -p reuse4d "
        "-scale none -print none "
        "-maxiter_noimp -1",
        solver);
    if( err ) return 1;

    err = lis_solve(A,b,x,solver);

    if( err!=LIS_ERR_ILL_ARG )
    {
        fprintf(stderr,
                "invalid maxiter_noimp returned %d, expected %d\n",
                (int)err,(int)LIS_ERR_ILL_ARG);
        return 1;
    }

    if( ordinary_create_count!=1 )
    {
        fprintf(stderr,
                "ordinary create count=%d, expected 1\n",
                ordinary_create_count);
        return 1;
    }

    /*
     * lis_solve() owns the preconditioner it creates and must
     * destroy it when lis_solve_kernel() returns an error.
     */
    if( ordinary_destroy_count!=1 )
    {
        fprintf(stderr,
                "ordinary destroy count=%d, expected 1\n",
                ordinary_destroy_count);
        return 1;
    }

    lis_solver_destroy(solver);
    lis_vector_destroy(b);
    lis_vector_destroy(x);

    return 0;
}

int main(int argc, char *argv[])
{
    LIS_MATRIX A_stag = NULL;
    LIS_MATRIX A_good = NULL;
    LIS_INT err;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;

    err = lis_precon_register_ex("reuse4d",
                                 ordinary_create,
                                 identity_psolve,
                                 identity_psolveh,
                                 user_destroy);
    if( err )
    {
        fprintf(stderr,
                "lis_precon_register_ex failed: %d\n",(int)err);
        return 1;
    }

    err = lis_precon_register_psd("reuse4d",
                                  psd_create,
                                  psd_update);
    if( err )
    {
        fprintf(stderr,
                "lis_precon_register_psd failed: %d\n",(int)err);
        return 1;
    }

    err = create_stagnation_matrix(&A_stag);
    if( err )
    {
        fprintf(stderr,
                "create_stagnation_matrix failed: %d\n",(int)err);
        return 1;
    }

    err = create_identity_matrix(&A_good);
    if( err )
    {
        fprintf(stderr,
                "create_identity_matrix failed: %d\n",(int)err);
        return 1;
    }

    if( test_reusable_psd(A_stag,A_good) )
    {
        return 1;
    }

    if( test_lis_solve_error_cleanup(A_good) )
    {
        return 1;
    }

    if( bad_context )
    {
        fprintf(stderr,"bad user context detected\n");
        return 1;
    }

    lis_matrix_destroy(A_stag);
    lis_matrix_destroy(A_good);

    lis_precon_register_free();
    lis_finalize();

    printf("LIS_PRECON_REUSE_STAGE4D PASSED\n");

    return 0;
}
