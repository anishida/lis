#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif
#include "lis.h"
#include <stdio.h>
#include <stdlib.h>

#define CTX_MAGIC 0x4c495350

typedef struct
{
    LIS_INT n;
} USER_IDENTITY;

typedef struct
{
    int magic;
} USER_PRECON_CTX;

static int create_count = 0;
static int psolve_count = 0;
static int psolveh_count = 0;
static int destroy_count = 0;
static int destroy_bad_context = 0;

static int fail_create_count = 0;
static int fail_destroy_count = 0;
static int fail_destroy_bad_context = 0;

static LIS_INT identity_matvec(void *user_data,
                               const LIS_SCALAR *x,
                               LIS_SCALAR *y)
{
    USER_IDENTITY *A = (USER_IDENTITY *)user_data;
    LIS_INT i;

    for(i=0;i<A->n;i++)
    {
        y[i] = x[i];
    }

    return LIS_SUCCESS;
}

static LIS_INT identity_matvech(void *user_data,
                                const LIS_SCALAR *x,
                                LIS_SCALAR *y)
{
    return identity_matvec(user_data,x,y);
}

static LIS_INT lifecycle_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    USER_PRECON_CTX *ctx;
    LIS_INT err;

    (void)solver;

    create_count++;

    ctx = (USER_PRECON_CTX *)malloc(sizeof(USER_PRECON_CTX));
    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic = CTX_MAGIC;

    err = lis_precon_set_user_data(precon,ctx);
    if( err )
    {
        free(ctx);
        return err;
    }

    return LIS_SUCCESS;
}

static LIS_INT lifecycle_psolve(LIS_SOLVER solver,
                                LIS_VECTOR b,
                                LIS_VECTOR x)
{
    USER_PRECON_CTX *ctx;
    void *user_data = NULL;
    LIS_INT i, err;

    psolve_count++;

    err = lis_precon_get_user_data(solver->precon,&user_data);
    if( err ) return err;

    ctx = (USER_PRECON_CTX *)user_data;
    if( ctx==NULL || ctx->magic!=CTX_MAGIC )
    {
        return LIS_FAILS;
    }

    for(i=0;i<b->n;i++)
    {
        x->value[i] = b->value[i];
    }

    return LIS_SUCCESS;
}

static LIS_INT lifecycle_psolveh(LIS_SOLVER solver,
                                 LIS_VECTOR b,
                                 LIS_VECTOR x)
{
    USER_PRECON_CTX *ctx;
    void *user_data = NULL;
    LIS_INT i, err;

    psolveh_count++;

    err = lis_precon_get_user_data(solver->precon,&user_data);
    if( err ) return err;

    ctx = (USER_PRECON_CTX *)user_data;
    if( ctx==NULL || ctx->magic!=CTX_MAGIC )
    {
        return LIS_FAILS;
    }

    for(i=0;i<b->n;i++)
    {
        x->value[i] = b->value[i];
    }

    return LIS_SUCCESS;
}

static LIS_INT lifecycle_destroy(LIS_PRECON precon)
{
    USER_PRECON_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    destroy_count++;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (USER_PRECON_CTX *)user_data;
    if( ctx==NULL || ctx->magic!=CTX_MAGIC )
    {
        destroy_bad_context = 1;
        return LIS_FAILS;
    }

    ctx->magic = 0;
    free(ctx);

    return lis_precon_set_user_data(precon,NULL);
}

static LIS_INT failing_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    USER_PRECON_CTX *ctx;
    LIS_INT err;

    (void)solver;

    fail_create_count++;

    ctx = (USER_PRECON_CTX *)malloc(sizeof(USER_PRECON_CTX));
    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic = CTX_MAGIC;

    err = lis_precon_set_user_data(precon,ctx);
    if( err )
    {
        free(ctx);
        return err;
    }

    return LIS_FAILS;
}

static LIS_INT failing_destroy(LIS_PRECON precon)
{
    USER_PRECON_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    fail_destroy_count++;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (USER_PRECON_CTX *)user_data;
    if( ctx==NULL || ctx->magic!=CTX_MAGIC )
    {
        fail_destroy_bad_context = 1;
        return LIS_FAILS;
    }

    ctx->magic = 0;
    free(ctx);

    return lis_precon_set_user_data(precon,NULL);
}

int main(int argc, char *argv[])
{
    LIS_MATRIX A;
    LIS_VECTOR b, x, exact;
    LIS_SOLVER solver;
    USER_IDENTITY shell;
    LIS_INT i, err;
    LIS_INT local_n, global_n;
    const LIS_INT n = 32;

    /*
     * Compile-time part of the RED test:
     * this typedef does not exist before Stage 4B2.
     */
    LIS_PRECON_DESTROY_XXX destroy_callback = lifecycle_destroy;

    (void)destroy_callback;

    lis_initialize(&argc,&argv);

    err = lis_precon_register_ex("lifeok",
                                 lifecycle_create,
                                 lifecycle_psolve,
                                 lifecycle_psolveh,
                                 lifecycle_destroy);
    if( err )
    {
        fprintf(stderr,"lis_precon_register_ex lifeok failed: %d\n",
                (int)err);
        lis_finalize();
        return 1;
    }

    err = lis_precon_register_ex("lifefail",
                                 failing_create,
                                 lifecycle_psolve,
                                 lifecycle_psolveh,
                                 failing_destroy);
    if( err )
    {
        fprintf(stderr,"lis_precon_register_ex lifefail failed: %d\n",
                (int)err);
        lis_precon_register_free();
        lis_finalize();
        return 1;
    }

    lis_matrix_create(LIS_COMM_WORLD,&A);
    lis_matrix_set_size(A,0,n);
    lis_matrix_get_size(A,&local_n,&global_n);
    (void)global_n;

    shell.n = local_n;

    err = lis_matrix_set_user(A,&shell,
                              identity_matvec,
                              identity_matvech);
    if( err )
    {
        fprintf(stderr,"lis_matrix_set_user failed: %d\n",(int)err);
        lis_matrix_destroy(A);
        lis_precon_register_free();
        lis_finalize();
        return 1;
    }

    lis_vector_duplicate(A,&b);
    lis_vector_duplicate(A,&x);
    lis_vector_duplicate(A,&exact);

    for(i=0;i<exact->n;i++)
    {
        exact->value[i] = 1.0;
    }

    err = lis_matvec(A,exact,b);
    if( err )
    {
        fprintf(stderr,"lis_matvec failed: %d\n",(int)err);
        return 1;
    }

    /*
     * Normal lifecycle.
     *
     * BiCG is intentional: it exercises both psolve and psolveh.
     */
    lis_vector_set_all(0.0,x);

    lis_solver_create(&solver);
    err = lis_solver_set_option(
        "-i bicg -p lifeok -tol 1.0e-12 -maxiter 100 -print none",
        solver);
    if( err )
    {
        fprintf(stderr,"lis_solver_set_option lifeok failed: %d\n",
                (int)err);
        return 1;
    }

    err = lis_solve(A,b,x,solver);
    if( err )
    {
        fprintf(stderr,"lis_solve lifeok failed: %d\n",(int)err);
        return 1;
    }

    lis_solver_destroy(solver);

    if( create_count!=1 )
    {
        fprintf(stderr,"normal create count = %d, expected 1\n",
                create_count);
        return 1;
    }

    if( psolve_count<1 )
    {
        fprintf(stderr,"normal psolve was not called\n");
        return 1;
    }

    if( psolveh_count<1 )
    {
        fprintf(stderr,"normal psolveh was not called\n");
        return 1;
    }

    if( destroy_count!=1 )
    {
        fprintf(stderr,"normal destroy count = %d, expected 1\n",
                destroy_count);
        return 1;
    }

    if( destroy_bad_context )
    {
        fprintf(stderr,"normal destroy received invalid user_data\n");
        return 1;
    }

    /*
     * Error lifecycle.
     *
     * failing_create stores user_data and then returns LIS_FAILS.
     * lis_precon_create() must route the partially-created object through
     * lis_precon_destroy(), which must invoke failing_destroy exactly once.
     */
    lis_vector_set_all(0.0,x);

    lis_solver_create(&solver);
    err = lis_solver_set_option(
        "-i bicg -p lifefail -tol 1.0e-12 -maxiter 100 -print none",
        solver);
    if( err )
    {
        fprintf(stderr,"lis_solver_set_option lifefail failed: %d\n",
                (int)err);
        return 1;
    }

    err = lis_solve(A,b,x,solver);

    if( err!=LIS_FAILS )
    {
        fprintf(stderr,
                "failing create returned solve status %d, expected %d\n",
                (int)err,(int)LIS_FAILS);
        return 1;
    }

    lis_solver_destroy(solver);

    if( fail_create_count!=1 )
    {
        fprintf(stderr,"error create count = %d, expected 1\n",
                fail_create_count);
        return 1;
    }

    if( fail_destroy_count!=1 )
    {
        fprintf(stderr,"error destroy count = %d, expected 1\n",
                fail_destroy_count);
        return 1;
    }

    if( fail_destroy_bad_context )
    {
        fprintf(stderr,"error destroy received invalid user_data\n");
        return 1;
    }

    lis_vector_destroy(exact);
    lis_vector_destroy(x);
    lis_vector_destroy(b);
    lis_matrix_destroy(A);

    lis_precon_register_free();
    lis_finalize();

    printf("PASSED\n");
    return 0;
}
