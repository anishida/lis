#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"
#include <stdio.h>
#include <stdlib.h>

#define PSD_MAGIC 0x50534443

#define CTX_NORMAL      1
#define CTX_CREATE_FAIL 2
#define CTX_UPDATE_FAIL 3

typedef struct
{
    int magic;
    int kind;
    int updates;
} USER_PSD_CTX;

static int normal_create_count = 0;
static int normal_update_count = 0;
static int normal_destroy_count = 0;

static int fail_create_count = 0;
static int fail_create_destroy_count = 0;

static int fail_update_create_count = 0;
static int fail_update_count = 0;
static int fail_update_destroy_count = 0;

static int destroy_bad_context = 0;

static LIS_INT ordinary_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    (void)solver;
    (void)precon;

    return LIS_SUCCESS;
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

static LIS_INT create_context(LIS_PRECON precon, int kind)
{
    USER_PSD_CTX *ctx;
    LIS_INT err;

    ctx = (USER_PSD_CTX *)malloc(sizeof(USER_PSD_CTX));
    if( ctx==NULL )
    {
        return LIS_OUT_OF_MEMORY;
    }

    ctx->magic = PSD_MAGIC;
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

static LIS_INT normal_psd_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    (void)solver;

    normal_create_count++;

    return create_context(precon,CTX_NORMAL);
}

static LIS_INT normal_psd_update(LIS_SOLVER solver, LIS_PRECON precon)
{
    USER_PSD_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    (void)solver;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (USER_PSD_CTX *)user_data;
    if( ctx==NULL ||
        ctx->magic!=PSD_MAGIC ||
        ctx->kind!=CTX_NORMAL )
    {
        return LIS_FAILS;
    }

    ctx->updates++;
    normal_update_count++;

    return LIS_SUCCESS;
}

static LIS_INT failing_psd_create(LIS_SOLVER solver, LIS_PRECON precon)
{
    LIS_INT err;

    (void)solver;

    fail_create_count++;

    err = create_context(precon,CTX_CREATE_FAIL);
    if( err ) return err;

    return LIS_FAILS;
}

static LIS_INT update_fail_psd_create(LIS_SOLVER solver,
                                      LIS_PRECON precon)
{
    (void)solver;

    fail_update_create_count++;

    return create_context(precon,CTX_UPDATE_FAIL);
}

static LIS_INT failing_psd_update(LIS_SOLVER solver, LIS_PRECON precon)
{
    USER_PSD_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    (void)solver;

    fail_update_count++;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (USER_PSD_CTX *)user_data;
    if( ctx==NULL ||
        ctx->magic!=PSD_MAGIC ||
        ctx->kind!=CTX_UPDATE_FAIL )
    {
        return LIS_FAILS;
    }

    return LIS_FAILS;
}

static LIS_INT user_destroy(LIS_PRECON precon)
{
    USER_PSD_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return err;

    ctx = (USER_PSD_CTX *)user_data;

    /*
     * A USERDEF object may have no PSD state at all.
     * In that case there is nothing for this callback to release.
     */
    if( ctx==NULL )
    {
        return LIS_SUCCESS;
    }

    if( ctx->magic!=PSD_MAGIC )
    {
        destroy_bad_context = 1;
        return LIS_FAILS;
    }

    switch( ctx->kind )
    {
    case CTX_NORMAL:
        normal_destroy_count++;
        break;
    case CTX_CREATE_FAIL:
        fail_create_destroy_count++;
        break;
    case CTX_UPDATE_FAIL:
        fail_update_destroy_count++;
        break;
    default:
        destroy_bad_context = 1;
        return LIS_FAILS;
    }

    ctx->magic = 0;
    free(ctx);

    return lis_precon_set_user_data(precon,NULL);
}

static LIS_INT create_test_matrix(LIS_MATRIX *A)
{
    LIS_INT i, is, ie, err;
    const LIS_INT n = 16;

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

static LIS_INT create_solver_for_precon(LIS_MATRIX A,
                                        const char *name,
                                        LIS_SOLVER *solver)
{
    char options[128];
    LIS_INT err;

    err = lis_solver_create(solver);
    if( err ) return err;

    snprintf(options,sizeof(options),
             "-i gmres -p %s -scale none -print none",
             name);

    err = lis_solver_set_option(options,*solver);
    if( err ) return err;

    return lis_solver_set_matrix(A,*solver);
}

int main(int argc, char *argv[])
{
    LIS_MATRIX A;
    LIS_SOLVER solver;
    LIS_PRECON precon;
    USER_PSD_CTX *ctx;
    void *user_data = NULL;
    LIS_INT err;

    err = lis_initialize(&argc,&argv);
    if( err ) return 1;

    /*
     * Normal PSD USERDEF.
     */
    err = lis_precon_register_ex("psdok",
                                 ordinary_create,
                                 identity_psolve,
                                 identity_psolveh,
                                 user_destroy);
    if( err )
    {
        fprintf(stderr,"register psdok failed: %d\n",(int)err);
        return 1;
    }

    err = lis_precon_register_psd("psdok",
                                  normal_psd_create,
                                  normal_psd_update);
    if( err )
    {
        fprintf(stderr,"register PSD psdok failed: %d\n",(int)err);
        return 1;
    }

    /*
     * USERDEF whose PSD create fails after allocating its context.
     */
    err = lis_precon_register_ex("psdcfail",
                                 ordinary_create,
                                 identity_psolve,
                                 identity_psolveh,
                                 user_destroy);
    if( err ) return 1;

    err = lis_precon_register_psd("psdcfail",
                                  failing_psd_create,
                                  normal_psd_update);
    if( err ) return 1;

    /*
     * USERDEF whose PSD update fails.
     */
    err = lis_precon_register_ex("psdufail",
                                 ordinary_create,
                                 identity_psolve,
                                 identity_psolveh,
                                 user_destroy);
    if( err ) return 1;

    err = lis_precon_register_psd("psdufail",
                                  update_fail_psd_create,
                                  failing_psd_update);
    if( err ) return 1;

    /*
     * Registered USERDEF without PSD callbacks.
     */
    err = lis_precon_register_ex("psdnone",
                                 ordinary_create,
                                 identity_psolve,
                                 identity_psolveh,
                                 user_destroy);
    if( err ) return 1;

    err = create_test_matrix(&A);
    if( err )
    {
        fprintf(stderr,"create_test_matrix failed: %d\n",(int)err);
        return 1;
    }

    /*
     * Normal create + repeated updates must retain the same context.
     */
    err = create_solver_for_precon(A,"psdok",&solver);
    if( err ) return 1;

    err = lis_precon_psd_create(solver,&precon);
    if( err )
    {
        fprintf(stderr,"normal psd_create failed: %d\n",(int)err);
        return 1;
    }

    err = lis_precon_psd_update(solver,precon);
    if( err )
    {
        fprintf(stderr,"normal first psd_update failed: %d\n",(int)err);
        return 1;
    }

    err = lis_precon_psd_update(solver,precon);
    if( err )
    {
        fprintf(stderr,"normal second psd_update failed: %d\n",(int)err);
        return 1;
    }

    err = lis_precon_get_user_data(precon,&user_data);
    if( err ) return 1;

    ctx = (USER_PSD_CTX *)user_data;
    if( ctx==NULL ||
        ctx->magic!=PSD_MAGIC ||
        ctx->kind!=CTX_NORMAL ||
        ctx->updates!=2 )
    {
        fprintf(stderr,"normal PSD context is invalid\n");
        return 1;
    }

    if( normal_create_count!=1 || normal_update_count!=2 )
    {
        fprintf(stderr,
                "normal counts create=%d update=%d, expected 1/2\n",
                normal_create_count,normal_update_count);
        return 1;
    }

    err = lis_precon_destroy(precon);
    if( err ) return 1;

    lis_solver_destroy(solver);

    if( normal_destroy_count!=1 )
    {
        fprintf(stderr,"normal destroy count=%d, expected 1\n",
                normal_destroy_count);
        return 1;
    }

    /*
     * A failing PSD create must propagate its error and destroy
     * partially-created USERDEF state exactly once.
     */
    err = create_solver_for_precon(A,"psdcfail",&solver);
    if( err ) return 1;

    err = lis_precon_psd_create(solver,&precon);

    if( err!=LIS_FAILS )
    {
        fprintf(stderr,
                "failing psd_create returned %d, expected %d\n",
                (int)err,(int)LIS_FAILS);
        return 1;
    }

    lis_solver_destroy(solver);

    if( fail_create_count!=1 ||
        fail_create_destroy_count!=1 )
    {
        fprintf(stderr,
                "create-fail counts create=%d destroy=%d, expected 1/1\n",
                fail_create_count,fail_create_destroy_count);
        return 1;
    }

    /*
     * A failing PSD update must propagate its error and destroy
     * the USERDEF preconditioner exactly once.
     */
    err = create_solver_for_precon(A,"psdufail",&solver);
    if( err ) return 1;

    err = lis_precon_psd_create(solver,&precon);
    if( err ) return 1;

    err = lis_precon_psd_update(solver,precon);

    if( err!=LIS_FAILS )
    {
        fprintf(stderr,
                "failing psd_update returned %d, expected %d\n",
                (int)err,(int)LIS_FAILS);
        return 1;
    }

    /*
     * precon was destroyed by lis_precon_psd_update() on error.
     * Do not destroy it a second time.
     */
    lis_solver_destroy(solver);

    if( fail_update_create_count!=1 ||
        fail_update_count!=1 ||
        fail_update_destroy_count!=1 )
    {
        fprintf(stderr,
                "update-fail counts create=%d update=%d destroy=%d, expected 1/1/1\n",
                fail_update_create_count,
                fail_update_count,
                fail_update_destroy_count);
        return 1;
    }

    /*
     * USERDEF without registered PSD callbacks remains unsupported.
     */
    err = create_solver_for_precon(A,"psdnone",&solver);
    if( err ) return 1;

    err = lis_precon_psd_create(solver,&precon);

    if( err!=LIS_ERR_NOT_IMPLEMENTED )
    {
        fprintf(stderr,
                "missing PSD callback returned %d, expected %d\n",
                (int)err,(int)LIS_ERR_NOT_IMPLEMENTED);
        return 1;
    }

    lis_solver_destroy(solver);

    if( destroy_bad_context )
    {
        fprintf(stderr,"destroy callback received invalid context\n");
        return 1;
    }

    lis_matrix_destroy(A);
    lis_precon_register_free();
    lis_finalize();

    printf("PASSED\n");
    return 0;
}
