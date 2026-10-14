#include <math.h>
#include <stdio.h>

#include <string.h>

#include "lislib.h"


extern LIS_INT
lis_precon_mfamg_get_smoother_omega(
    LIS_MATRIX A,
    LIS_REAL *omega);


LIS_INT
lis_precon_mfamg_get_hierarchy_info(
    LIS_PRECON precon,
    LIS_INT max_levels,
    LIS_INT *nlevels,
    LIS_INT *coarse_n,
    LIS_INT *level_n,
    LIS_INT *owns_A);

extern LIS_INT
lis_precon_mfamg_get_global_row(
    LIS_MATRIX A,
    LIS_INT local_row,
    LIS_INT *nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals);

extern LIS_INT
lis_precon_mfamg_get_distributed_aggregate_info(
    LIS_MATRIX A,
    LIS_INT map_capacity,
    LIS_INT *local_nc,
    LIS_INT *global_nc,
    LIS_INT *coarse_is,
    LIS_INT *coarse_ie,
    LIS_INT *map_n,
    LIS_INT *fine_to_coarse);


extern LIS_INT
lis_precon_mfamg_get_distributed_prolongation_info(
    LIS_MATRIX A,
    LIS_VECTOR z,
    LIS_INT fine_capacity,
    LIS_INT coarse_capacity,
    LIS_INT *map_n,
    LIS_INT *local_nc,
    LIS_INT *global_nc,
    LIS_SCALAR *fine_prolong,
    LIS_SCALAR *coarse_z);


extern LIS_INT
lis_precon_mfamg_get_distributed_hierarchy_info(
    LIS_MATRIX A,
    LIS_VECTOR z,
    LIS_INT level_capacity,
    LIS_INT terminal_row_capacity,
    LIS_INT terminal_nnz_capacity,
    LIS_INT *nlevels,
    LIS_INT *local_n,
    LIS_INT *global_n,
    LIS_INT *owns_A,
    LIS_INT *terminal_np,
    LIS_INT *terminal_nnz,
    LIS_INT *terminal_row,
    LIS_INT *terminal_global_cols,
    LIS_SCALAR *terminal_values,
    LIS_INT terminal_dense_capacity,
    LIS_INT terminal_vector_capacity,
    LIS_SCALAR *terminal_dense,
    LIS_SCALAR *terminal_solution);


extern LIS_INT
lis_precon_mfamg_get_distributed_coarse_info(
    LIS_MATRIX A,
    LIS_VECTOR z,
    LIS_INT row_capacity,
    LIS_INT nnz_capacity,
    LIS_INT *local_nc,
    LIS_INT *global_nc,
    LIS_INT *coarse_is,
    LIS_INT *coarse_ie,
    LIS_INT *coarse_np,
    LIS_INT *local_nnz,
    LIS_INT *row,
    LIS_INT *global_cols,
    LIS_SCALAR *values);



LIS_INT
lis_precon_mfamg_build_multivector_basis(
    LIS_INT n,
    LIS_INT aggregate_count,
    LIS_INT *aggregate,
    LIS_INT mode_count,
    LIS_VECTOR basis[],
    LIS_INT **aggregate_rank_out,
    LIS_INT **aggregate_offset_out,
    LIS_INT **row_ptr_out,
    LIS_SCALAR **prolong_out,
    LIS_SCALAR **coarse_basis_out,
    LIS_INT *coarse_n_out);


#define TEST_N 24
#define HARD_N 31
#define STRENGTH_N 8
#define SOLVE_TOL 1.0e-8
#define MATCH_TOL 1.0e-12

#define EXPECTED_MFAMG_LEVELS 3
#define EXPECTED_MFAMG_LEVEL0_N 24
#define EXPECTED_MFAMG_LEVEL1_N 12
#define EXPECTED_MFAMG_LEVEL2_N 6
#define EXPECTED_MFAMG_COARSE_N 6


typedef struct
{
    LIS_INT n;
    LIS_INT gn;
    LIS_INT is;
} USER_CTX;


/* ============================================================
 * USER matrix callbacks
 * ============================================================ */

static LIS_INT
user_matvec(
    void *user_data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    USER_CTX *ctx;
    LIS_INT i;


    ctx =
        (USER_CTX *)user_data;


    for(i=0;i<ctx->n;i++)
    {
        y[i] =
            (LIS_SCALAR)2.0
            *
            x[i];

        if(i>0)
        {
            y[i] -=
                x[i-1];
        }

        if(i+1<ctx->n)
        {
            y[i] -=
                x[i+1];
        }
    }


    return LIS_SUCCESS;
}


static LIS_INT
user_diagonal(
    void *user_data,
    LIS_SCALAR *d)
{
    USER_CTX *ctx;
    LIS_INT i;


    ctx =
        (USER_CTX *)user_data;


    for(i=0;i<ctx->n;i++)
    {
        d[i] =
            (LIS_SCALAR)2.0;
    }


    return LIS_SUCCESS;
}


static LIS_INT
user_row(
    void *user_data,
    LIS_INT local_row,
    LIS_INT *nnz,
    LIS_INT *cols,
    LIS_SCALAR *vals)
{
    USER_CTX *ctx;

    LIS_INT global_row;
    LIS_INT need;
    LIS_INT k;


    ctx =
        (USER_CTX *)user_data;


    if(
        ctx==NULL
        ||
        nnz==NULL
        ||
        local_row<0
        ||
        local_row>=ctx->n)
    {
        return LIS_ERR_ILL_ARG;
    }


    global_row =
        ctx->is
        +
        local_row;


    if(
        global_row<0
        ||
        global_row>=ctx->gn)
    {
        return LIS_ERR_ILL_ARG;
    }


    need = 1;


    if(global_row>0)
    {
        need++;
    }


    if(global_row+1<ctx->gn)
    {
        need++;
    }


    if(
        cols==NULL
        &&
        vals==NULL)
    {
        *nnz = need;

        return LIS_SUCCESS;
    }


    if(
        cols==NULL
        ||
        vals==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    if(*nnz<need)
    {
        *nnz = need;

        return LIS_ERR_ILL_ARG;
    }


    k = 0;


    if(global_row>0)
    {
        cols[k] =
            global_row-1;

        vals[k] =
            (LIS_SCALAR)-1.0;

        k++;
    }


    cols[k] =
        global_row;

    vals[k] =
        (LIS_SCALAR)2.0;

    k++;


    if(global_row+1<ctx->gn)
    {
        cols[k] =
            global_row+1;

        vals[k] =
            (LIS_SCALAR)-1.0;

        k++;
    }


    *nnz = k;


    return LIS_SUCCESS;
}

/* ============================================================
 * Matrix creation
 * ============================================================ */

static LIS_INT
create_csr(
    LIS_MATRIX *A)
{
    LIS_INT i;
    LIS_INT err;


    *A = NULL;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            A);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            *A,
            0,
            TEST_N);

    if(err)
        return err;


    for(i=0;i<TEST_N;i++)
    {
        if(i>0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-1.0,
                    *A);

            if(err)
                return err;
        }


        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)2.0,
                *A);

        if(err)
            return err;


        if(i+1<TEST_N)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+1,
                    (LIS_SCALAR)-1.0,
                    *A);

            if(err)
                return err;
        }
    }


    return
        lis_matrix_assemble(
            *A);
}


static LIS_INT
create_user(
    USER_CTX *ctx,
    LIS_MATRIX *A)
{
    LIS_INT err;


    if(
        ctx==NULL
        ||
        A==NULL)
    {
        return LIS_ERR_ILL_ARG;
    }


    *A = NULL;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            A);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            *A,
            0,
            TEST_N);

    if(err)
        return err;


    ctx->n =
        (*A)->n;

    ctx->gn =
        (*A)->gn;

    ctx->is =
        (*A)->is;


    err =
        lis_matrix_set_user(
            *A,
            ctx,
            user_matvec,
            user_matvec);

    if(err)
        return err;


    err =
        lis_matrix_set_user_diagonal(
            *A,
            user_diagonal);

    if(err)
        return err;


    return
        lis_matrix_set_user_row(
            *A,
            user_row);
}

/* ============================================================
 * Helpers
 * ============================================================ */

static void
build_manufactured_rhs(
    LIS_VECTOR b)
{
    LIS_SCALAR exact[TEST_N];
    LIS_INT i;


    for(i=0;i<TEST_N;i++)
    {
        exact[i] =
            (LIS_SCALAR)
            (
                1.0
                +
                0.1
                *
                (double)i
            );
    }


    for(i=0;i<TEST_N;i++)
    {
        b->value[i] =
            (LIS_SCALAR)2.0
            *
            exact[i];

        if(i>0)
        {
            b->value[i] -=
                exact[i-1];
        }

        if(i+1<TEST_N)
        {
            b->value[i] -=
                exact[i+1];
        }
    }
}


static LIS_INT
relative_residual(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_REAL *rel)
{
    LIS_VECTOR ax;
    LIS_VECTOR r;

    LIS_REAL bnrm;
    LIS_REAL rnrm;

    LIS_INT err;


    ax = NULL;
    r = NULL;


    err =
        lis_vector_duplicate(
            A,
            &ax);

    if(err)
        return err;


    err =
        lis_vector_duplicate(
            A,
            &r);

    if(err)
    {
        lis_vector_destroy(
            ax);

        return err;
    }


    err =
        lis_matvec(
            A,
            x,
            ax);

    if(err)
        goto cleanup;


    err =
        lis_vector_copy(
            b,
            r);

    if(err)
        goto cleanup;


    err =
        lis_vector_axpy(
            (LIS_SCALAR)-1.0,
            ax,
            r);

    if(err)
        goto cleanup;


    err =
        lis_vector_nrm2(
            b,
            &bnrm);

    if(err)
        goto cleanup;


    err =
        lis_vector_nrm2(
            r,
            &rnrm);

    if(err)
        goto cleanup;


    if(!(bnrm>0.0))
    {
        err =
            LIS_ERR_ILL_ARG;

        goto cleanup;
    }


    *rel =
        rnrm
        /
        bnrm;


cleanup:

    lis_vector_destroy(
        r);

    lis_vector_destroy(
        ax);


    return err;
}


static LIS_REAL
max_relative_difference(
    const LIS_SCALAR *a,
    const LIS_SCALAR *b)
{
    LIS_REAL diff;
    LIS_REAL scale;
    LIS_REAL max_diff;

    LIS_INT i;


    max_diff =
        0.0;


    for(i=0;i<TEST_N;i++)
    {
        diff =
            fabs(
                a[i]
                -
                b[i]);

        scale =
            1.0
            +
            fabs(
                a[i]);

        diff /=
            scale;


        if(diff>max_diff)
        {
            max_diff =
                diff;
        }
    }


    return max_diff;
}


/* ============================================================
 * Direct preconditioner application
 * ============================================================ */

static int
run_direct_case(
    const char *name,
    LIS_MATRIX A,
    LIS_SCALAR *result)
{
    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR z;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver;
    LIS_PRECON precon;

    LIS_REAL rel;
    LIS_REAL mag;

    LIS_INT i;
    LIS_INT err;
    LIS_INT destroy_err;

    LIS_INT hierarchy_nlevels;
    LIS_INT hierarchy_coarse_n;

    LIS_INT hierarchy_level_n[
        EXPECTED_MFAMG_LEVELS];

    LIS_INT hierarchy_owns_A[
        EXPECTED_MFAMG_LEVELS];

    int failed;


    b = NULL;
    x = NULL;
    z = NULL;

    solver = NULL;
    precon = NULL;

    hierarchy_nlevels = -1;
    hierarchy_coarse_n = -1;

    for(
        i=0;
        i<EXPECTED_MFAMG_LEVELS;
        i++)
    {
        hierarchy_level_n[i] = -1;
        hierarchy_owns_A[i] = -1;
    }

    failed = 0;


    err =
        lis_vector_duplicate(
            A,
            &b);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
        goto fail;


    build_manufactured_rhs(
        b);


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    err =
        lis_solver_create(
            &solver);

    if(err)
        goto fail;


    err =
        lis_solver_set_option(
            "-i gmres "
            "-p mfamg "
            "-scale none "
            "-print none",
            solver);

    if(err)
        goto fail;


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
        goto fail;


    err =
        lis_solver_set_matrix(
            A,
            solver);

    if(err)
        goto fail;


    err =
        lis_precon_create(
            solver,
            &precon);


    printf(
        "%s DIRECT_CREATE rc=%d\n",
        name,
        (int)err);


    if(err)
    {
        /*
         * Failed creation is cleaned up inside lis_precon_create().
         */
        precon = NULL;

        failed = 1;

        goto cleanup;
    }


    solver->precon =
        precon;


    err =
        lis_precon_mfamg_get_hierarchy_info(
            precon,
            EXPECTED_MFAMG_LEVELS,
            &hierarchy_nlevels,
            &hierarchy_coarse_n,
            hierarchy_level_n,
            hierarchy_owns_A);


    printf(
        "%s HIERARCHY rc=%d "
        "nlevels=%d coarse_n=%d "
        "levels=%d,%d,%d "
        "owns=%d,%d,%d\n",
        name,
        (int)err,
        (int)hierarchy_nlevels,
        (int)hierarchy_coarse_n,
        (int)hierarchy_level_n[0],
        (int)hierarchy_level_n[1],
        (int)hierarchy_level_n[2],
        (int)hierarchy_owns_A[0],
        (int)hierarchy_owns_A[1],
        (int)hierarchy_owns_A[2]);


    if(err)
    {
        fprintf(
            stderr,
            "%s MFAMG hierarchy query failed rc=%d\n",
            name,
            (int)err);

        failed = 1;

        goto cleanup;
    }


    if(
        hierarchy_nlevels
            !=
        EXPECTED_MFAMG_LEVELS
        ||
        hierarchy_coarse_n
            !=
        EXPECTED_MFAMG_COARSE_N
        ||
        hierarchy_level_n[0]
            !=
        EXPECTED_MFAMG_LEVEL0_N
        ||
        hierarchy_level_n[1]
            !=
        EXPECTED_MFAMG_LEVEL1_N
        ||
        hierarchy_level_n[2]
            !=
        EXPECTED_MFAMG_LEVEL2_N
        ||
        hierarchy_owns_A[0]
            !=
        LIS_FALSE
        ||
        hierarchy_owns_A[1]
            !=
        LIS_TRUE
        ||
        hierarchy_owns_A[2]
            !=
        LIS_TRUE)
    {
        fprintf(
            stderr,
            "%s unexpected MFAMG hierarchy\n",
            name);

        failed = 1;
    }


    if(
        hierarchy_nlevels<3
        ||
        hierarchy_coarse_n
            !=
        hierarchy_level_n[
            EXPECTED_MFAMG_LEVELS-1])
    {
        fprintf(
            stderr,
            "%s MFAMG hierarchy is not genuinely multilevel "
            "or terminal size is inconsistent\n",
            name);

        failed = 1;
    }


    if(
        solver->near_nullspace_coarse
        !=
        NULL)
    {
        fprintf(
            stderr,
            "%s unexpected legacy coarse state before psolve\n",
            name);

        failed = 1;
    }


    err =
        lis_psolve(
            solver,
            b,
            x);


    printf(
        "%s DIRECT_PSOLVE rc=%d\n",
        name,
        (int)err);


    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    if(
        solver->near_nullspace_coarse
        !=
        NULL)
    {
        fprintf(
            stderr,
            "%s unexpected legacy coarse state after psolve\n",
            name);

        failed = 1;
    }


    err =
        relative_residual(
            A,
            b,
            x,
            &rel);

    if(err)
        goto fail;


    printf(
        "%s DIRECT_RELRES=%.12e\n",
        name,
        (double)rel);


    if(!(rel<1.0))
    {
        fprintf(
            stderr,
            "%s direct MFAMG did not improve residual\n",
            name);

        failed = 1;
    }


    for(i=0;i<TEST_N;i++)
    {
        mag =
            fabs(
                x->value[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            fprintf(
                stderr,
                "%s direct result is non-finite\n",
                name);

            failed = 1;

            break;
        }


        result[i] =
            x->value[i];
    }


cleanup:

    if(precon)
    {
        if(solver)
        {
            solver->precon =
                NULL;
        }


        destroy_err =
            lis_precon_destroy(
                precon);


        printf(
            "%s DIRECT_DESTROY rc=%d\n",
            name,
            (int)destroy_err);


        if(destroy_err)
        {
            failed = 1;
        }


        precon = NULL;
    }


    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(z)
    {
        lis_vector_destroy(
            z);
    }

    if(x)
    {
        lis_vector_destroy(
            x);
    }

    if(b)
    {
        lis_vector_destroy(
            b);
    }


    return failed;


fail:

    fprintf(
        stderr,
        "%s direct setup failed rc=%d\n",
        name,
        (int)err);

    failed = 1;

    goto cleanup;
}


/* ============================================================
 * Full Krylov solve
 * ============================================================ */

static int
run_solve_case(
    const char *name,
    LIS_MATRIX A,
    LIS_SCALAR *result)
{
    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR z;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver;

    LIS_REAL rel;
    LIS_REAL mag;

    LIS_INT i;
    LIS_INT err;

    int failed;


    b = NULL;
    x = NULL;
    z = NULL;
    solver = NULL;

    failed = 0;


    err =
        lis_vector_duplicate(
            A,
            &b);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
        goto fail;


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    err =
        lis_solver_create(
            &solver);

    if(err)
        goto fail;


    err =
        lis_solver_set_option(
            "-i gmres "
            "-p mfamg "
            "-scale none "
            "-tol 1.0e-10 "
            "-maxiter 200 "
            "-print none",
            solver);

    if(err)
        goto fail;


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
        goto fail;


    err =
        lis_solve(
            A,
            b,
            x,
            solver);


    printf(
        "%s SOLVE rc=%d\n",
        name,
        (int)err);


    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    if(
        solver->near_nullspace_coarse
        !=
        NULL)
    {
        fprintf(
            stderr,
            "%s full solve used legacy coarse state\n",
            name);

        failed = 1;
    }


    err =
        relative_residual(
            A,
            b,
            x,
            &rel);

    if(err)
        goto fail;


    printf(
        "%s SOLVE_RELRES=%.12e\n",
        name,
        (double)rel);


    if(
        !(
            rel
            <=
            SOLVE_TOL
        ))
    {
        fprintf(
            stderr,
            "%s full solve residual is too large\n",
            name);

        failed = 1;
    }


    for(i=0;i<TEST_N;i++)
    {
        mag =
            fabs(
                x->value[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            fprintf(
                stderr,
                "%s solve result is non-finite\n",
                name);

            failed = 1;

            break;
        }


        result[i] =
            x->value[i];
    }


cleanup:

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(z)
    {
        lis_vector_destroy(
            z);
    }

    if(x)
    {
        lis_vector_destroy(
            x);
    }

    if(b)
    {
        lis_vector_destroy(
            b);
    }


    return failed;


fail:

    fprintf(
        stderr,
        "%s solve setup failed rc=%d\n",
        name,
        (int)err);

    failed = 1;

    goto cleanup;
}


/* ============================================================
 * Zero near-nullspace rejection
 * ============================================================ */

static int
test_zero_mode(
    LIS_MATRIX A)
{
    LIS_VECTOR z;
    LIS_VECTOR basis[1];

    LIS_SOLVER solver;
    LIS_PRECON precon;

    LIS_INT err;

    int failed;


    z = NULL;
    solver = NULL;
    precon = NULL;

    failed = 0;


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
        return 1;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z);


    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    err =
        lis_solver_set_option(
            "-i gmres "
            "-p mfamg "
            "-scale none "
            "-print none",
            solver);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    err =
        lis_solver_set_matrix(
            A,
            solver);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    err =
        lis_precon_create(
            solver,
            &precon);


    if(err)
    {
        /*
         * lis_precon_create() already released the failed object.
         */
        precon = NULL;
    }


    printf(
        "ZERO_MODE rc=%d expected=%d\n",
        (int)err,
        (int)LIS_ERR_ILL_ARG);


    if(
        err
        !=
        LIS_ERR_ILL_ARG)
    {
        failed = 1;
    }


cleanup:

    if(precon)
    {
        lis_precon_destroy(
            precon);
    }

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(z)
    {
        lis_vector_destroy(
            z);
    }


    return failed;
}



static int
run_multivector_direct_case(
    const char *name,
    LIS_MATRIX A,
    LIS_SCALAR *result)
{
    enum
    {
        MODE_COUNT = 2,
        MAX_LEVELS = 32
    };

    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR z0;
    LIS_VECTOR z1;

    LIS_VECTOR basis[MODE_COUNT];

    LIS_SOLVER solver;
    LIS_PRECON precon;

    LIS_REAL rel;
    LIS_REAL mag;

    LIS_INT hierarchy_nlevels;
    LIS_INT hierarchy_coarse_n;

    LIS_INT hierarchy_level_n[MAX_LEVELS];
    LIS_INT hierarchy_owns_A[MAX_LEVELS];

    LIS_INT i;
    LIS_INT level_id;
    LIS_INT err;
    LIS_INT destroy_err;

    int failed;


    b = NULL;
    x = NULL;
    z0 = NULL;
    z1 = NULL;

    solver = NULL;
    precon = NULL;

    hierarchy_nlevels = -1;
    hierarchy_coarse_n = -1;

    failed = 0;


    for(i=0;i<MAX_LEVELS;i++)
    {
        hierarchy_level_n[i] = -1;
        hierarchy_owns_A[i] = -1;
    }


    err =
        lis_vector_duplicate(
            A,
            &b);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &z0);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &z1);

    if(err)
        goto fail;


    build_manufactured_rhs(
        b);


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z0);


    for(i=0;i<TEST_N;i++)
    {
        if(i<10)
        {
            z1->value[i] =
                (LIS_SCALAR)(i%2);
        }
        else
        {
            z1->value[i] =
                (LIS_SCALAR)0.0;
        }
    }


    err =
        lis_solver_create(
            &solver);

    if(err)
        goto fail;


    err =
        lis_solver_set_option(
            "-i gmres "
            "-p mfamg "
            "-scale none "
            "-print none",
            solver);

    if(err)
        goto fail;


    basis[0] =
        z0;

    basis[1] =
        z1;


    err =
        lis_solver_set_near_nullspace(
            solver,
            MODE_COUNT,
            basis);

    if(err)
        goto fail;


    err =
        lis_solver_set_matrix(
            A,
            solver);

    if(err)
        goto fail;


    err =
        lis_precon_create(
            solver,
            &precon);


    printf(
        "%s MULTIVECTOR_CREATE rc=%d\n",
        name,
        (int)err);


    if(err)
    {
        precon = NULL;
        failed = 1;
        goto cleanup;
    }


    solver->precon =
        precon;


    err =
        lis_precon_mfamg_get_hierarchy_info(
            precon,
            MAX_LEVELS,
            &hierarchy_nlevels,
            &hierarchy_coarse_n,
            hierarchy_level_n,
            hierarchy_owns_A);


    printf(
        "%s MULTIVECTOR_HIERARCHY rc=%d "
        "nlevels=%d coarse_n=%d\n",
        name,
        (int)err,
        (int)hierarchy_nlevels,
        (int)hierarchy_coarse_n);


    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    if(
        hierarchy_nlevels!=3
        ||
        hierarchy_nlevels>MAX_LEVELS
        ||
        hierarchy_level_n[0]!=TEST_N
        ||
        hierarchy_level_n[1]!=17
        ||
        hierarchy_level_n[2]<=0
        ||
        hierarchy_level_n[2]>=17
        ||
        hierarchy_level_n[2]>16
        ||
        hierarchy_coarse_n
            !=
        hierarchy_level_n[
            hierarchy_nlevels-1]
        ||
        hierarchy_owns_A[0]!=LIS_FALSE)
    {
        failed = 1;
    }


    for(
        level_id=0;
        level_id<hierarchy_nlevels;
        level_id++)
    {
        printf(
            "%s MULTIVECTOR_LEVEL "
            "level=%d n=%d owns=%d\n",
            name,
            (int)level_id,
            (int)hierarchy_level_n[level_id],
            (int)hierarchy_owns_A[level_id]);


        if(
            hierarchy_level_n[level_id]<=0
            ||
            (
                level_id>0
                &&
                hierarchy_owns_A[level_id]
                    !=
                LIS_TRUE
            ))
        {
            failed = 1;
        }


        if(
            level_id+1<hierarchy_nlevels
            &&
            hierarchy_level_n[level_id+1]
                >=
            hierarchy_level_n[level_id])
        {
            failed = 1;
        }
    }


    if(failed)
        goto cleanup;


    printf(
        "%s MULTIVECTOR_MULTILEVEL=PASSED "
        "nlevels=%d levels=%d,%d,%d terminal=%d\n",
        name,
        (int)hierarchy_nlevels,
        (int)hierarchy_level_n[0],
        (int)hierarchy_level_n[1],
        (int)hierarchy_level_n[2],
        (int)hierarchy_coarse_n);


    err =
        lis_psolve(
            solver,
            b,
            x);


    printf(
        "%s MULTIVECTOR_PSOLVE rc=%d\n",
        name,
        (int)err);


    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        relative_residual(
            A,
            b,
            x,
            &rel);

    if(err)
        goto fail;


    printf(
        "%s MULTIVECTOR_RELRES=%.12e\n",
        name,
        (double)rel);


    if(
        rel!=rel
        ||
        !(rel<1.0))
    {
        failed = 1;
    }


    for(i=0;i<TEST_N;i++)
    {
        mag =
            fabs(
                x->value[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            failed = 1;
            break;
        }


        result[i] =
            x->value[i];
    }


cleanup:

    if(precon)
    {
        if(solver)
        {
            solver->precon =
                NULL;
        }


        destroy_err =
            lis_precon_destroy(
                precon);

        printf(
            "%s MULTIVECTOR_DESTROY rc=%d\n",
            name,
            (int)destroy_err);


        if(destroy_err)
            failed = 1;

        precon = NULL;
    }


    if(solver)
        lis_solver_destroy(solver);

    if(z1)
        lis_vector_destroy(z1);

    if(z0)
        lis_vector_destroy(z0);

    if(x)
        lis_vector_destroy(x);

    if(b)
        lis_vector_destroy(b);


    return failed;


fail:

    fprintf(
        stderr,
        "%s multivector setup failed rc=%d\n",
        name,
        (int)err);

    failed = 1;

    goto cleanup;
}


static int
check_serial_multivector(
    LIS_MATRIX A,
    LIS_MATRIX U)
{
    LIS_SCALAR csr_result[TEST_N];
    LIS_SCALAR user_result[TEST_N];

    LIS_REAL diff;

    int failed;


    failed = 0;


    failed |=
        run_multivector_direct_case(
            "CSR",
            A,
            csr_result);


    failed |=
        run_multivector_direct_case(
            "USER",
            U,
            user_result);


    if(failed)
    {
        fprintf(
            stderr,
            "MFAMG SERIAL MULTIVECTOR TEST FAILED\n");

        return 1;
    }


    diff =
        max_relative_difference(
            csr_result,
            user_result);


    printf(
        "MULTIVECTOR_CSR_USER_MAX_REL_DIFF=%.12e\n",
        (double)diff);


    if(diff>MATCH_TOL)
    {
        fprintf(
            stderr,
            "MFAMG serial multivector CSR/USER mismatch\n");

        return 1;
    }


    printf(
        "MFAMG SERIAL MULTIVECTOR TEST PASSED\n");


    return 0;
}


/* ============================================================
 * Serial contract
 * ============================================================ */




static int
check_mfamg_adaptive_smoother(void)
{
    LIS_MATRIX A;

    LIS_INT i;
    LIS_INT err;

    LIS_REAL omega;
    LIS_REAL expected;
    LIS_REAL diff;

    int failed;


    A = NULL;

    omega = 0.0;

    expected =
        (LIS_REAL)(20.0/33.0);

    failed = 0;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        lis_matrix_set_size(
            A,
            0,
            4);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    for(i=0;i<4;i++)
    {
        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)1.0,
                A);

        if(err)
        {
            failed = 1;
            goto cleanup;
        }


        if(i>0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-0.6,
                    A);

            if(err)
            {
                failed = 1;
                goto cleanup;
            }
        }


        if(i<3)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+1,
                    (LIS_SCALAR)-0.6,
                    A);

            if(err)
            {
                failed = 1;
                goto cleanup;
            }
        }
    }


    err =
        lis_matrix_assemble(
            A);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        lis_precon_mfamg_get_smoother_omega(
            A,
            &omega);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    diff =
        fabs(
            omega
            -
            expected);


    if(
        omega!=omega
        ||
        diff>1.0e-12)
    {
        failed = 1;
    }


cleanup:

    if(A)
    {
        lis_matrix_destroy(
            A);
    }


    printf(
        "SMOOTHER_OMEGA=%s "
        "omega=%.12e expected=%.12e\n",
        failed
        ?
        "FAILED"
        :
        "PASSED",
        (double)omega,
        (double)expected);


    return failed;
}


static int
check_mfamg_multivector_basis_foundation(void)
{
    enum
    {
        LOCAL_N = 6,
        LOCAL_AGGREGATES = 2,
        LOCAL_MODES = 2
    };

    LIS_VECTOR z0;
    LIS_VECTOR z1;

    LIS_VECTOR basis[LOCAL_MODES];

    LIS_INT aggregate[LOCAL_N] =
    {
        0, 0, 0,
        1, 1, 1
    };

    LIS_INT *aggregate_rank;
    LIS_INT *aggregate_offset;
    LIS_INT *row_ptr;

    LIS_SCALAR *prolong;
    LIS_SCALAR *coarse_basis;

    LIS_INT coarse_n;

    LIS_INT i;
    LIS_INT a;
    LIS_INT m;
    LIS_INT q;
    LIS_INT r;
    LIS_INT err;

    LIS_SCALAR recovered;
    LIS_SCALAR gram;
    LIS_SCALAR target;

    LIS_REAL reconstruction_error;
    LIS_REAL orthogonality_error;
    LIS_REAL value_error;

    int failed;


    z0 = NULL;
    z1 = NULL;

    aggregate_rank = NULL;
    aggregate_offset = NULL;
    row_ptr = NULL;

    prolong = NULL;
    coarse_basis = NULL;

    coarse_n = -1;

    reconstruction_error = 0.0;
    orthogonality_error = 0.0;

    failed = 0;


    err =
        lis_vector_create(
            LIS_COMM_WORLD,
            &z0);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_set_size(
            z0,
            0,
            LOCAL_N);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            z0,
            &z1);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    /*
     * Aggregate 0 has local rank 2.
     */
    z0->value[0] = (LIS_SCALAR)1.0;
    z0->value[1] = (LIS_SCALAR)1.0;
    z0->value[2] = (LIS_SCALAR)1.0;

    z1->value[0] = (LIS_SCALAR)0.0;
    z1->value[1] = (LIS_SCALAR)1.0;
    z1->value[2] = (LIS_SCALAR)2.0;


    /*
     * Aggregate 1 has local rank 1 because z1=2*z0.
     */
    z0->value[3] = (LIS_SCALAR)1.0;
    z0->value[4] = (LIS_SCALAR)1.0;
    z0->value[5] = (LIS_SCALAR)1.0;

    z1->value[3] = (LIS_SCALAR)2.0;
    z1->value[4] = (LIS_SCALAR)2.0;
    z1->value[5] = (LIS_SCALAR)2.0;


    basis[0] =
        z0;

    basis[1] =
        z1;


    err =
        lis_precon_mfamg_build_multivector_basis(
            LOCAL_N,
            LOCAL_AGGREGATES,
            aggregate,
            LOCAL_MODES,
            basis,
            &aggregate_rank,
            &aggregate_offset,
            &row_ptr,
            &prolong,
            &coarse_basis,
            &coarse_n);

    if(err)
    {
        failed = 1;
        goto report;
    }


    if(
        aggregate_rank==NULL
        ||
        aggregate_offset==NULL
        ||
        row_ptr==NULL
        ||
        prolong==NULL
        ||
        coarse_basis==NULL
        ||
        coarse_n!=3
        ||
        aggregate_rank[0]!=2
        ||
        aggregate_rank[1]!=1
        ||
        aggregate_offset[0]!=0
        ||
        aggregate_offset[1]!=2
        ||
        aggregate_offset[2]!=3
        ||
        row_ptr[0]!=0
        ||
        row_ptr[1]!=2
        ||
        row_ptr[2]!=4
        ||
        row_ptr[3]!=6
        ||
        row_ptr[4]!=7
        ||
        row_ptr[5]!=8
        ||
        row_ptr[6]!=9)
    {
        failed = 1;
        goto report;
    }


    /*
     * Verify B = P * Bc.
     */
    for(m=0;m<LOCAL_MODES;m++)
    {
        for(i=0;i<LOCAL_N;i++)
        {
            a =
                aggregate[i];

            recovered =
                (LIS_SCALAR)0.0;


            for(
                q=0;
                q<aggregate_rank[a];
                q++)
            {
                recovered +=
                    prolong[
                        row_ptr[i]
                        +
                        q]
                    *
                    coarse_basis[
                        (
                            aggregate_offset[a]
                            +
                            q
                        )
                        *
                        LOCAL_MODES
                        +
                        m];
            }


            value_error =
                fabs(
                    recovered
                    -
                    basis[m]->value[i]);


            if(
                value_error
                >
                reconstruction_error)
            {
                reconstruction_error =
                    value_error;
            }
        }
    }


    /*
     * Verify Q_a^H Q_a = I.
     */
    for(a=0;a<LOCAL_AGGREGATES;a++)
    {
        for(
            q=0;
            q<aggregate_rank[a];
            q++)
        {
            for(
                r=0;
                r<aggregate_rank[a];
                r++)
            {
                gram =
                    (LIS_SCALAR)0.0;


                for(i=0;i<LOCAL_N;i++)
                {
                    if(aggregate[i]!=a)
                        continue;


                    gram +=
                        conj(
                            prolong[
                                row_ptr[i]
                                +
                                q])
                        *
                        prolong[
                            row_ptr[i]
                            +
                            r];
                }


                target =
                    (q==r)
                    ?
                    (LIS_SCALAR)1.0
                    :
                    (LIS_SCALAR)0.0;


                value_error =
                    fabs(
                        gram
                        -
                        target);


                if(
                    value_error
                    >
                    orthogonality_error)
                {
                    orthogonality_error =
                        value_error;
                }
            }
        }
    }


    if(
        reconstruction_error>1.0e-11
        ||
        orthogonality_error>1.0e-11)
    {
        failed = 1;
    }


report:

    printf(
        "MULTIVECTOR_BASIS=%s "
        "coarse_n=%d "
        "ranks=%d,%d "
        "entries=%d "
        "reconstruction=% .12e "
        "orthogonality=% .12e\n",
        failed
        ?
        "FAILED"
        :
        "PASSED",
        (int)coarse_n,
        aggregate_rank
        ?
        (int)aggregate_rank[0]
        :
        -1,
        aggregate_rank
        ?
        (int)aggregate_rank[1]
        :
        -1,
        row_ptr
        ?
        (int)row_ptr[LOCAL_N]
        :
        -1,
        (double)reconstruction_error,
        (double)orthogonality_error);


cleanup:

    if(coarse_basis)
        lis_free(coarse_basis);

    if(prolong)
        lis_free(prolong);

    if(row_ptr)
        lis_free(row_ptr);

    if(aggregate_offset)
        lis_free(aggregate_offset);

    if(aggregate_rank)
        lis_free(aggregate_rank);

    if(z1)
        lis_vector_destroy(z1);

    if(z0)
        lis_vector_destroy(z0);


    return failed;
}


/* ============================================================
 * Persistent MFAMG regression helpers
 * ============================================================ */

static LIS_INT
create_scaled_csr(
    LIS_SCALAR scale,
    LIS_MATRIX *A)
{
    LIS_INT i;
    LIS_INT err;


    if(
        A==NULL
        ||
        !(fabs(scale)>0.0))
    {
        return LIS_ERR_ILL_ARG;
    }


    *A = NULL;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            A);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            *A,
            0,
            TEST_N);

    if(err)
        goto fail;


    for(i=0;i<TEST_N;i++)
    {
        if(i>0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-1.0
                    *
                    scale,
                    *A);

            if(err)
                goto fail;
        }


        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)2.0
                *
                scale,
                *A);

        if(err)
            goto fail;


        if(i+1<TEST_N)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+1,
                    (LIS_SCALAR)-1.0
                    *
                    scale,
                    *A);

            if(err)
                goto fail;
        }
    }


    err =
        lis_matrix_assemble(
            *A);

    if(err)
        goto fail;


    return LIS_SUCCESS;


fail:

    if(*A)
    {
        lis_matrix_destroy(
            *A);

        *A = NULL;
    }


    return err;
}


/* ============================================================
 * In-place preconditioner application
 * ============================================================ */

static int
check_mfamg_alias_safety(void)
{
    LIS_MATRIX A;

    LIS_VECTOR b;
    LIS_VECTOR out;
    LIS_VECTOR alias;
    LIS_VECTOR z;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver;
    LIS_PRECON precon;

    LIS_REAL diff;

    LIS_INT err;

    int failed;


    A = NULL;

    b = NULL;
    out = NULL;
    alias = NULL;
    z = NULL;

    solver = NULL;
    precon = NULL;

    failed = 0;


    err =
        create_csr(
            &A);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &out);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &alias);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
        goto fail;


    build_manufactured_rhs(
        b);


    err =
        lis_vector_copy(
            b,
            alias);

    if(err)
        goto fail;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        out);


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    err =
        lis_solver_create(
            &solver);

    if(err)
        goto fail;


    err =
        lis_solver_set_option(
            "-i gmres "
            "-p mfamg "
            "-scale none "
            "-print none",
            solver);

    if(err)
        goto fail;


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
        goto fail;


    err =
        lis_solver_set_matrix(
            A,
            solver);

    if(err)
        goto fail;


    err =
        lis_precon_create(
            solver,
            &precon);

    if(err)
        goto fail;


    solver->precon =
        precon;


    err =
        lis_psolve(
            solver,
            b,
            out);

    if(err)
        goto fail;


    err =
        lis_psolve(
            solver,
            alias,
            alias);

    if(err)
        goto fail;


    diff =
        max_relative_difference(
            out->value,
            alias->value);


    printf(
        "MFAMG_ALIAS_MAX_REL_DIFF=%.12e\n",
        (double)diff);


    if(
        diff!=diff
        ||
        diff>(LIS_REAL)1.0e-12)
    {
        fprintf(
            stderr,
            "MFAMG alias-safe psolve mismatch\n");

        failed = 1;
    }


cleanup:

    if(solver)
    {
        solver->precon =
            NULL;
    }


    if(precon)
    {
        lis_precon_destroy(
            precon);
    }


    if(solver)
    {
        lis_solver_destroy(
            solver);
    }


    if(z)
    {
        lis_vector_destroy(
            z);
    }


    if(alias)
    {
        lis_vector_destroy(
            alias);
    }


    if(out)
    {
        lis_vector_destroy(
            out);
    }


    if(b)
    {
        lis_vector_destroy(
            b);
    }


    if(A)
    {
        lis_matrix_destroy(
            A);
    }


    if(failed)
    {
        fprintf(
            stderr,
            "MFAMG ALIAS SAFETY TEST FAILED\n");

        return 1;
    }


    printf(
        "MFAMG ALIAS SAFETY TEST PASSED\n");


    return 0;


fail:

    fprintf(
        stderr,
        "MFAMG alias safety setup failed rc=%d\n",
        (int)err);

    failed = 1;

    goto cleanup;
}


/* ============================================================
 * Persistent MFAMG lifecycle
 * ============================================================ */

static int
check_mfamg_psd_lifecycle(void)
{
    LIS_MATRIX matrices[3];

    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR z;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver;

    LIS_PRECON precon;
    LIS_PRECON precon0;

    void *mfamg0;

    LIS_REAL scales[3];
    LIS_REAL direct_rel;
    LIS_REAL solve_rel;

    LIS_INT hierarchy_nlevels;
    LIS_INT hierarchy_coarse_n;

    LIS_INT hierarchy_level_n[
        EXPECTED_MFAMG_LEVELS];

    LIS_INT hierarchy_owns_A[
        EXPECTED_MFAMG_LEVELS];

    LIS_INT system_id;
    LIS_INT i;
    LIS_INT err;
    LIS_INT destroy_err;

    int precon_alive;
    int failed;


    matrices[0] = NULL;
    matrices[1] = NULL;
    matrices[2] = NULL;

    b = NULL;
    x = NULL;
    z = NULL;

    solver = NULL;

    precon = NULL;
    precon0 = NULL;

    mfamg0 = NULL;

    scales[0] =
        (LIS_REAL)1.0;

    scales[1] =
        (LIS_REAL)1.25;

    scales[2] =
        (LIS_REAL)0.75;

    precon_alive = 0;
    failed = 0;


    err =
        create_csr(
            &matrices[0]);

    if(err)
        goto fail;


    err =
        create_scaled_csr(
            (LIS_SCALAR)scales[1],
            &matrices[1]);

    if(err)
        goto fail;


    err =
        create_scaled_csr(
            (LIS_SCALAR)scales[2],
            &matrices[2]);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            matrices[0],
            &b);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            matrices[0],
            &x);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            matrices[0],
            &z);

    if(err)
        goto fail;


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    err =
        lis_solver_create(
            &solver);

    if(err)
        goto fail;


    err =
        lis_solver_set_option(
            "-i fgmres "
            "-p mfamg "
            "-scale none "
            "-tol 1.0e-10 "
            "-maxiter 200 "
            "-print none",
            solver);

    if(err)
        goto fail;


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
        goto fail;


    for(
        system_id=0;
        system_id<3;
        system_id++)
    {
        build_manufactured_rhs(
            b);


        if(scales[system_id]!=(LIS_REAL)1.0)
        {
            err =
                lis_vector_scale(
                    (LIS_SCALAR)
                    scales[system_id],
                    b);

            if(err)
                goto fail;
        }


        lis_vector_set_all(
            (LIS_SCALAR)0.0,
            x);


        err =
            lis_solver_set_matrix(
                matrices[system_id],
                solver);

        if(err)
            goto fail;


        if(system_id==0)
        {
            err =
                lis_precon_psd_create(
                    solver,
                    &precon);

            printf(
                "MFAMG_PSD_CREATE rc=%d\n",
                (int)err);


            if(err)
            {
                precon = NULL;
                failed = 1;

                goto cleanup;
            }


            precon_alive = 1;
            precon0 = precon;

            mfamg0 =
                (void *)precon->mfamg;

            if(mfamg0==NULL)
            {
                fprintf(
                    stderr,
                    "MFAMG PSD initial state is null\n");

                failed = 1;
                goto cleanup;
            }
        }
        else
        {
            err =
                lis_precon_psd_update(
                    solver,
                    precon);

            printf(
                "MFAMG_PSD_UPDATE%d rc=%d\n",
                (int)system_id,
                (int)err);

            if(
                !err
                &&
                (void *)precon->mfamg!=mfamg0)
            {
                fprintf(
                    stderr,
                    "MFAMG PSD same-pattern update replaced state system=%d\n",
                    (int)system_id);

                failed = 1;
                goto cleanup;
            }


            if(err)
            {
                /*
                 * The PSD dispatcher owns destruction on update
                 * failure.
                 */
                precon_alive = 0;
                precon = NULL;

                failed = 1;

                goto cleanup;
            }
        }


        if(
            precon==NULL
            ||
            precon!=precon0
            ||
            precon->A!=matrices[system_id]
            ||
            precon->precon_type
                !=
                LIS_PRECON_TYPE_MFAMG)
        {
            fprintf(
                stderr,
                "MFAMG PSD identity mismatch "
                "system=%d\n",
                (int)system_id);

            failed = 1;

            goto cleanup;
        }


        hierarchy_nlevels = -1;
        hierarchy_coarse_n = -1;

        for(
            i=0;
            i<EXPECTED_MFAMG_LEVELS;
            i++)
        {
            hierarchy_level_n[i] = -1;
            hierarchy_owns_A[i] = -1;
        }


        err =
            lis_precon_mfamg_get_hierarchy_info(
                precon,
                EXPECTED_MFAMG_LEVELS,
                &hierarchy_nlevels,
                &hierarchy_coarse_n,
                hierarchy_level_n,
                hierarchy_owns_A);

        if(err)
            goto fail;


        if(
            hierarchy_nlevels
                !=
                EXPECTED_MFAMG_LEVELS
            ||
            hierarchy_coarse_n
                !=
                EXPECTED_MFAMG_COARSE_N
            ||
            hierarchy_level_n[0]
                !=
                EXPECTED_MFAMG_LEVEL0_N
            ||
            hierarchy_level_n[1]
                !=
                EXPECTED_MFAMG_LEVEL1_N
            ||
            hierarchy_level_n[2]
                !=
                EXPECTED_MFAMG_LEVEL2_N
            ||
            hierarchy_owns_A[0]
                !=
                LIS_FALSE
            ||
            hierarchy_owns_A[1]
                !=
                LIS_TRUE
            ||
            hierarchy_owns_A[2]
                !=
                LIS_TRUE)
        {
            fprintf(
                stderr,
                "MFAMG PSD hierarchy mismatch "
                "system=%d levels=%d "
                "sizes=%d,%d,%d coarse=%d\n",
                (int)system_id,
                (int)hierarchy_nlevels,
                (int)hierarchy_level_n[0],
                (int)hierarchy_level_n[1],
                (int)hierarchy_level_n[2],
                (int)hierarchy_coarse_n);

            failed = 1;

            goto cleanup;
        }


        /*
         * Exercise the persistent hierarchy directly first.
         */
        solver->precon =
            precon;


        err =
            lis_psolve(
                solver,
                b,
                x);


        solver->precon =
            NULL;


        if(err)
        {
            fprintf(
                stderr,
                "MFAMG PSD direct apply failed "
                "system=%d rc=%d\n",
                (int)system_id,
                (int)err);

            failed = 1;

            goto cleanup;
        }


        err =
            relative_residual(
                matrices[system_id],
                b,
                x,
                &direct_rel);

        if(err)
            goto fail;


        if(
            direct_rel!=direct_rel
            ||
            !(direct_rel<1.0))
        {
            fprintf(
                stderr,
                "MFAMG PSD direct residual failed "
                "system=%d rel=%.12e\n",
                (int)system_id,
                (double)direct_rel);

            failed = 1;

            goto cleanup;
        }


        /*
         * The full repeated solve uses the same caller-owned
         * persistent preconditioner.
         */
        lis_vector_set_all(
            (LIS_SCALAR)0.0,
            x);


        err =
            lis_solve_kernel(
                matrices[system_id],
                b,
                x,
                solver,
                precon);

        if(err)
        {
            fprintf(
                stderr,
                "MFAMG PSD solve failed "
                "system=%d rc=%d\n",
                (int)system_id,
                (int)err);

            failed = 1;

            goto cleanup;
        }


        err =
            relative_residual(
                matrices[system_id],
                b,
                x,
                &solve_rel);

        if(err)
            goto fail;


        printf(
            "MFAMG_PSD_SYSTEM%d "
            "scale=%.6e "
            "state_reused=%d "
            "direct_relres=%.12e "
            "solve_relres=%.12e "
            "levels=%d,%d,%d\n",
            (int)system_id,
            (double)scales[system_id],
            (int)
            (
                (void *)precon->mfamg
                ==
                mfamg0
            ),
            (double)direct_rel,
            (double)solve_rel,
            (int)hierarchy_level_n[0],
            (int)hierarchy_level_n[1],
            (int)hierarchy_level_n[2]);


        if(
            solve_rel!=solve_rel
            ||
            !(
                solve_rel
                <=
                SOLVE_TOL
            ))
        {
            fprintf(
                stderr,
                "MFAMG PSD solve residual failed "
                "system=%d rel=%.12e\n",
                (int)system_id,
                (double)solve_rel);

            failed = 1;

            goto cleanup;
        }
    }


cleanup:

    if(precon_alive && precon)
    {
        destroy_err =
            lis_precon_destroy(
                precon);

        printf(
            "MFAMG_PSD_DESTROY rc=%d\n",
            (int)destroy_err);


        if(destroy_err)
        {
            failed = 1;
        }


        precon = NULL;
        precon_alive = 0;
    }


    if(solver)
    {
        solver->precon =
            NULL;

        lis_solver_destroy(
            solver);
    }


    if(z)
    {
        lis_vector_destroy(
            z);
    }


    if(x)
    {
        lis_vector_destroy(
            x);
    }


    if(b)
    {
        lis_vector_destroy(
            b);
    }


    for(
        i=2;
        i>=0;
        i--)
    {
        if(matrices[i])
        {
            lis_matrix_destroy(
                matrices[i]);
        }
    }


    if(failed)
    {
        fprintf(
            stderr,
            "MFAMG PSD LIFECYCLE TEST FAILED\n");

        return 1;
    }


    printf(
        "MFAMG PSD LIFECYCLE TEST PASSED\n");


    return 0;


fail:

    fprintf(
        stderr,
        "MFAMG PSD lifecycle setup failed rc=%d\n",
        (int)err);

    failed = 1;

    goto cleanup;
}


static LIS_INT
create_rewired_csr(
    LIS_MATRIX *A)
{
    LIS_INT i;
    LIS_INT err;


    if(A==NULL)
        return LIS_ERR_ILL_ARG;


    *A = NULL;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            A);

    if(err)
        return err;


    err =
        lis_matrix_set_size(
            *A,
            0,
            TEST_N);

    if(err)
        return err;


    /*
     * Keep the same number of entries per row as create_csr(),
     * but change the graph connectivity.
     *
     * Removed:
     *
     *     5 -- 6
     *    10 -- 11
     *
     * Added:
     *
     *     5 -- 10
     *     6 -- 11
     *
     * Therefore n, nnz, and CSR ptr[] remain unchanged while
     * CSR index[] differs.
     */
    for(i=0;i<TEST_N;i++)
    {
        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)4.0,
                *A);

        if(err)
            return err;
    }


    for(i=0;i<TEST_N-1;i++)
    {
        if(
            i==5
            ||
            i==10)
        {
            continue;
        }


        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i+1,
                (LIS_SCALAR)-1.0,
                *A);

        if(err)
            return err;


        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i+1,
                i,
                (LIS_SCALAR)-1.0,
                *A);

        if(err)
            return err;
    }


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            5,
            10,
            (LIS_SCALAR)-1.0,
            *A);

    if(err)
        return err;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            10,
            5,
            (LIS_SCALAR)-1.0,
            *A);

    if(err)
        return err;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            6,
            11,
            (LIS_SCALAR)-1.0,
            *A);

    if(err)
        return err;


    err =
        lis_matrix_set_value(
            LIS_INS_VALUE,
            11,
            6,
            (LIS_SCALAR)-1.0,
            *A);

    if(err)
        return err;


    return
        lis_matrix_assemble(
            *A);
}


static int
check_mfamg_psd_pattern_fallback(void)
{
    LIS_MATRIX A0;
    LIS_MATRIX A1;

    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR exact;
    LIS_VECTOR z;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver;

    LIS_PRECON precon;
    LIS_PRECON precon0;

    void *mfamg0;

    LIS_REAL relres;

    LIS_INT i;
    LIS_INT err;

    LIS_INT same_nnz;
    LIS_INT same_ptr;
    LIS_INT different_index;

    int precon_alive;
    int failed;


    A0 = NULL;
    A1 = NULL;

    b = NULL;
    x = NULL;
    exact = NULL;
    z = NULL;

    solver = NULL;

    precon = NULL;
    precon0 = NULL;

    mfamg0 = NULL;

    relres =
        (LIS_REAL)-1.0;

    precon_alive = 0;
    failed = 0;


    err =
        create_csr(
            &A0);

    if(err)
        goto fail;


    err =
        create_rewired_csr(
            &A1);

    if(err)
        goto fail;


    if(
        A0->matrix_type!=LIS_MATRIX_CSR
        ||
        A1->matrix_type!=LIS_MATRIX_CSR
        ||
        A0->n!=A1->n
        ||
        A0->n!=TEST_N
        ||
        A0->ptr==NULL
        ||
        A1->ptr==NULL
        ||
        A0->index==NULL
        ||
        A1->index==NULL)
    {
        fprintf(
            stderr,
            "MFAMG pattern fallback matrix setup is invalid\n");

        failed = 1;
        goto cleanup;
    }


    same_nnz =
        (
            A0->nnz
            ==
            A1->nnz
        );


    same_ptr =
        LIS_TRUE;

    for(
        i=0;
        i<=A0->n;
        i++)
    {
        if(A0->ptr[i]!=A1->ptr[i])
        {
            same_ptr =
                LIS_FALSE;

            break;
        }
    }


    different_index =
        LIS_FALSE;

    if(same_nnz)
    {
        for(
            i=0;
            i<A0->nnz;
            i++)
        {
            if(A0->index[i]!=A1->index[i])
            {
                different_index =
                    LIS_TRUE;

                break;
            }
        }
    }


    printf(
        "MFAMG_PATTERN_VARIANT "
        "same_nnz=%d "
        "same_ptr=%d "
        "different_index=%d\n",
        (int)same_nnz,
        (int)same_ptr,
        (int)different_index);


    if(
        !same_nnz
        ||
        !same_ptr
        ||
        !different_index)
    {
        fprintf(
            stderr,
            "MFAMG pattern fallback test did not isolate CSR index change\n");

        failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            A1,
            &b);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A1,
            &x);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A1,
            &exact);

    if(err)
        goto fail;


    err =
        lis_vector_duplicate(
            A1,
            &z);

    if(err)
        goto fail;


    for(
        i=0;
        i<TEST_N;
        i++)
    {
        exact->value[i] =
            (LIS_SCALAR)
            (
                1.0
                +
                0.05
                *
                (double)i
            );
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    err =
        lis_matvec(
            A1,
            exact,
            b);

    if(err)
        goto fail;


    err =
        lis_solver_create(
            &solver);

    if(err)
        goto fail;


    err =
        lis_solver_set_option(
            "-i fgmres "
            "-p mfamg "
            "-scale none "
            "-tol 1.0e-10 "
            "-maxiter 200 "
            "-print none",
            solver);

    if(err)
        goto fail;


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
        goto fail;


    err =
        lis_solver_set_matrix(
            A0,
            solver);

    if(err)
        goto fail;


    err =
        lis_precon_psd_create(
            solver,
            &precon);

    printf(
        "MFAMG_PATTERN_FALLBACK_CREATE rc=%d\n",
        (int)err);


    if(err)
    {
        precon = NULL;
        failed = 1;

        goto cleanup;
    }


    precon_alive = 1;

    precon0 =
        precon;

    mfamg0 =
        (void *)precon->mfamg;


    if(mfamg0==NULL)
    {
        fprintf(
            stderr,
            "MFAMG pattern fallback initial state is null\n");

        failed = 1;
        goto cleanup;
    }


    err =
        lis_solver_set_matrix(
            A1,
            solver);

    if(err)
        goto fail;


    err =
        lis_precon_psd_update(
            solver,
            precon);

    printf(
        "MFAMG_PATTERN_FALLBACK_UPDATE rc=%d\n",
        (int)err);


    if(err)
    {
        /*
         * PSD update dispatcher destroys the persistent
         * preconditioner when its update callback fails.
         */
        precon_alive = 0;
        precon = NULL;

        failed = 1;
        goto cleanup;
    }


    if(precon!=precon0)
    {
        fprintf(
            stderr,
            "MFAMG pattern fallback replaced persistent preconditioner object\n");

        failed = 1;
        goto cleanup;
    }


    if(
        precon->mfamg==NULL
        ||
        (void *)precon->mfamg
            ==
        mfamg0)
    {
        fprintf(
            stderr,
            "MFAMG changed-pattern update did not replace state\n");

        failed = 1;
        goto cleanup;
    }


    if(precon->A!=A1)
    {
        fprintf(
            stderr,
            "MFAMG changed-pattern update did not install new matrix\n");

        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    err =
        lis_solve_kernel(
            A1,
            b,
            x,
            solver,
            precon);

    if(err)
        goto fail;


    err =
        relative_residual(
            A1,
            b,
            x,
            &relres);

    if(err)
        goto fail;


    printf(
        "MFAMG_PATTERN_FALLBACK "
        "precon_preserved=%d "
        "state_replaced=%d "
        "matrix_updated=%d "
        "relres=%.12e\n",
        (int)
        (
            precon
            ==
            precon0
        ),
        (int)
        (
            (void *)precon->mfamg
            !=
            mfamg0
        ),
        (int)
        (
            precon->A
            ==
            A1
        ),
        (double)relres);


    if(
        relres!=relres
        ||
        relres>(LIS_REAL)1.0e-8)
    {
        fprintf(
            stderr,
            "MFAMG pattern fallback solve residual is too large\n");

        failed = 1;
        goto cleanup;
    }


cleanup:

    if(solver)
    {
        solver->precon =
            NULL;
    }


    if(
        precon_alive
        &&
        precon)
    {
        lis_precon_destroy(
            precon);

        precon = NULL;
    }


    if(solver)
    {
        lis_solver_destroy(
            solver);

        solver = NULL;
    }


    if(z)
    {
        lis_vector_destroy(
            z);

        z = NULL;
    }


    if(exact)
    {
        lis_vector_destroy(
            exact);

        exact = NULL;
    }


    if(x)
    {
        lis_vector_destroy(
            x);

        x = NULL;
    }


    if(b)
    {
        lis_vector_destroy(
            b);

        b = NULL;
    }


    if(A1)
    {
        lis_matrix_destroy(
            A1);

        A1 = NULL;
    }


    if(A0)
    {
        lis_matrix_destroy(
            A0);

        A0 = NULL;
    }


    if(failed)
    {
        fprintf(
            stderr,
            "MFAMG PATTERN FALLBACK TEST FAILED\n");

        return 1;
    }


    printf(
        "MFAMG PATTERN FALLBACK TEST PASSED\n");


    return 0;


fail:

    fprintf(
        stderr,
        "MFAMG pattern fallback operation failed rc=%d\n",
        (int)err);

    failed = 1;

    goto cleanup;
}


static int
run_serial(void)
{
    LIS_MATRIX A;
    LIS_MATRIX U;

    USER_CTX ctx;

    LIS_SCALAR csr_direct[TEST_N];
    LIS_SCALAR user_direct[TEST_N];

    LIS_SCALAR csr_solve[TEST_N];
    LIS_SCALAR user_solve[TEST_N];

    LIS_REAL diff;

    LIS_INT err;

    int failed;


    A = NULL;
    U = NULL;

    failed = 0;

    failed |=
        check_mfamg_multivector_basis_foundation();


    failed |=
        check_mfamg_adaptive_smoother();


    failed |=
        check_mfamg_alias_safety();


    failed |=
        check_mfamg_psd_lifecycle();


    failed |=
        check_mfamg_psd_pattern_fallback();


    err =
        create_csr(
            &A);

    if(err)
    {
        fprintf(
            stderr,
            "CSR matrix creation failed rc=%d\n",
            (int)err);

        failed = 1;

        goto cleanup;
    }


    err =
        create_user(
            &ctx,
            &U);

    if(err)
    {
        fprintf(
            stderr,
            "USER matrix creation failed rc=%d\n",
            (int)err);

        failed = 1;

        goto cleanup;
    }


    failed |=
        run_direct_case(
            "CSR",
            A,
            csr_direct);


    failed |=
        run_direct_case(
            "USER",
            U,
            user_direct);


    diff =
        max_relative_difference(
            csr_direct,
            user_direct);


    printf(
        "DIRECT_CSR_USER_MAX_REL_DIFF=%.12e\n",
        (double)diff);


    if(diff>MATCH_TOL)
    {
        failed = 1;
    }


    failed |=
        run_solve_case(
            "CSR",
            A,
            csr_solve);


    failed |=
        run_solve_case(
            "USER",
            U,
            user_solve);


    diff =
        max_relative_difference(
            csr_solve,
            user_solve);


    printf(
        "SOLVE_CSR_USER_MAX_REL_DIFF=%.12e\n",
        (double)diff);


    if(diff>MATCH_TOL)
    {
        failed = 1;
    }


    failed |=
        check_serial_multivector(
            A,
            U);


    failed |=
        test_zero_mode(
            A);


cleanup:

    if(U)
    {
        lis_matrix_destroy(
            U);
    }

    if(A)
    {
        lis_matrix_destroy(
            A);
    }


    if(failed)
    {
        fprintf(
            stderr,
            "MFAMG SERIAL TEST FAILED\n");

        return 1;
    }


    printf(
        "MFAMG SERIAL TEST PASSED\n");


    return 0;
}


/* ============================================================
 * MPI global-row regression and serial-only rejection
 * ============================================================ */

static int
check_mfamg_global_rows(
    const char *name,
    LIS_MATRIX A,
    int rank)
{
    LIS_INT local_row;
    LIS_INT global_row;
    LIS_INT expected_nnz;
    LIS_INT nnz;
    LIS_INT k;
    LIS_INT err;

    LIS_INT cols[3];
    LIS_SCALAR vals[3];

    int seen_lower;
    int seen_diag;
    int seen_upper;
    int failed;


    if(
        name==NULL
        ||
        A==NULL)
    {
        return 1;
    }


    failed = 0;


    for(
        local_row=0;
        local_row<A->n;
        local_row++)
    {
        global_row =
            A->is
            +
            local_row;


        expected_nnz = 1;

        if(global_row>0)
        {
            expected_nnz++;
        }

        if(global_row+1<A->gn)
        {
            expected_nnz++;
        }


        nnz = 0;


        err =
            lis_precon_mfamg_get_global_row(
                A,
                local_row,
                &nnz,
                NULL,
                NULL);


        if(
            err
            ||
            nnz!=expected_nnz)
        {
            fprintf(
                stderr,
                "%s MPI row query failed on rank %d row %d: rc=%d nnz=%d expected=%d\n",
                name,
                rank,
                (int)global_row,
                (int)err,
                (int)nnz,
                (int)expected_nnz);

            failed = 1;

            continue;
        }


        nnz = 3;


        err =
            lis_precon_mfamg_get_global_row(
                A,
                local_row,
                &nnz,
                cols,
                vals);


        if(
            err
            ||
            nnz!=expected_nnz)
        {
            fprintf(
                stderr,
                "%s MPI row fill failed on rank %d row %d: rc=%d nnz=%d expected=%d\n",
                name,
                rank,
                (int)global_row,
                (int)err,
                (int)nnz,
                (int)expected_nnz);

            failed = 1;

            continue;
        }


        seen_lower = 0;
        seen_diag = 0;
        seen_upper = 0;


        for(k=0;k<nnz;k++)
        {
            if(
                global_row>0
                &&
                cols[k]==global_row-1
                &&
                fabs(
                    vals[k]
                    -
                    (LIS_SCALAR)-1.0)
                    <=
                    MATCH_TOL
                &&
                !seen_lower)
            {
                seen_lower = 1;

                continue;
            }


            if(
                cols[k]==global_row
                &&
                fabs(
                    vals[k]
                    -
                    (LIS_SCALAR)2.0)
                    <=
                    MATCH_TOL
                &&
                !seen_diag)
            {
                seen_diag = 1;

                continue;
            }


            if(
                global_row+1<A->gn
                &&
                cols[k]==global_row+1
                &&
                fabs(
                    vals[k]
                    -
                    (LIS_SCALAR)-1.0)
                    <=
                    MATCH_TOL
                &&
                !seen_upper)
            {
                seen_upper = 1;

                continue;
            }


            fprintf(
                stderr,
                "%s MPI row entry failed on rank %d row %d entry %d col=%d\n",
                name,
                rank,
                (int)global_row,
                (int)k,
                (int)cols[k]);

            failed = 1;
        }


        if(
            (global_row>0 && !seen_lower)
            ||
            !seen_diag
            ||
            (global_row+1<A->gn && !seen_upper))
        {
            fprintf(
                stderr,
                "%s MPI row coverage failed on rank %d row %d\n",
                name,
                rank,
                (int)global_row);

            failed = 1;
        }
    }


    return failed;
}


static int
check_mfamg_distributed_aggregates(
    LIS_MATRIX A,
    int rank,
    int nprocs)
{
    LIS_INT local_nc;
    LIS_INT global_nc;
    LIS_INT coarse_is;
    LIS_INT coarse_ie;
    LIS_INT map_n;

    LIS_INT expected_local_nc;
    LIS_INT expected_global_nc;
    LIS_INT expected_coarse_is;
    LIS_INT expected_coarse_ie;

    LIS_INT fine_to_coarse[TEST_N];

    LIS_INT i;
    LIS_INT global_fine;
    LIS_INT expected_gid;
    LIS_INT err;

    int local_failed;
    int global_failed;


    local_nc = -1;
    global_nc = -1;
    coarse_is = -1;
    coarse_ie = -1;
    map_n = -1;

    local_failed = 0;
    global_failed = 0;


    if(
        A==NULL
        ||
        nprocs!=2
        ||
        TEST_N!=24
        ||
        A->n!=12
        ||
        A->is%2!=0)
    {
        local_failed = 1;
    }
    else
    {
        err =
            lis_precon_mfamg_get_distributed_aggregate_info(
                A,
                TEST_N,
                &local_nc,
                &global_nc,
                &coarse_is,
                &coarse_ie,
                &map_n,
                fine_to_coarse);


        if(err)
        {
            fprintf(
                stderr,
                "MPI aggregate query failed on rank %d rc=%d\n",
                rank,
                (int)err);

            local_failed = 1;
        }
        else
        {
            expected_local_nc =
                A->n/2;

            expected_global_nc =
                TEST_N/2;

            expected_coarse_is =
                A->is/2;

            expected_coarse_ie =
                expected_coarse_is
                +
                expected_local_nc;


            if(
                local_nc!=expected_local_nc
                ||
                global_nc!=expected_global_nc
                ||
                coarse_is!=expected_coarse_is
                ||
                coarse_ie!=expected_coarse_ie
                ||
                map_n!=A->np)
            {
                fprintf(
                    stderr,
                    "MPI aggregate ownership mismatch on rank %d: "
                    "local_nc=%d global_nc=%d coarse=[%d,%d) map_n=%d np=%d\n",
                    rank,
                    (int)local_nc,
                    (int)global_nc,
                    (int)coarse_is,
                    (int)coarse_ie,
                    (int)map_n,
                    (int)A->np);

                local_failed = 1;
            }


            for(i=0;i<A->n;i++)
            {
                expected_gid =
                    coarse_is
                    +
                    i/2;


                if(
                    fine_to_coarse[i]
                    !=
                    expected_gid)
                {
                    fprintf(
                        stderr,
                        "MPI owned aggregate ID mismatch on rank %d "
                        "local=%d got=%d expected=%d\n",
                        rank,
                        (int)i,
                        (int)fine_to_coarse[i],
                        (int)expected_gid);

                    local_failed = 1;
                }
            }


            for(i=A->n;i<A->np;i++)
            {
                global_fine =
                    A->l2g_map[
                        i
                        -
                        A->n];

                expected_gid =
                    global_fine/2;


                if(
                    fine_to_coarse[i]
                    !=
                    expected_gid)
                {
                    fprintf(
                        stderr,
                        "MPI halo aggregate ID mismatch on rank %d "
                        "fine=%d got=%d expected=%d\n",
                        rank,
                        (int)global_fine,
                        (int)fine_to_coarse[i],
                        (int)expected_gid);

                    local_failed = 1;
                }
            }
        }
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_failed =
        local_failed;

#endif


    if(rank==0)
    {
        printf(
            "MPI_AGGREGATE_OWNERSHIP=%s ranks=%d global_nc=%d\n",
            global_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs,
            (int)global_nc);
    }


    return global_failed;
}

static int
check_mfamg_strength_aggregation(
    int rank,
    int nprocs)
{
    LIS_MATRIX A;

    LIS_INT i;
    LIS_INT local_i;
    LIS_INT err;

    LIS_INT local_nc;
    LIS_INT global_nc;
    LIS_INT coarse_is;
    LIS_INT coarse_ie;
    LIS_INT map_n;

    LIS_INT fine_to_coarse[STRENGTH_N];

    LIS_INT expected;

    int local_failed;
    int global_failed;


    A = NULL;

    local_nc = 0;
    global_nc = 0;
    coarse_is = 0;
    coarse_ie = 0;
    map_n = 0;

    local_failed = 0;
    global_failed = 0;


    if(nprocs!=2)
        return 0;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_matrix_set_size(
            A,
            0,
            STRENGTH_N);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    if(
        A->gn!=STRENGTH_N
        ||
        A->n!=4)
    {
        local_failed = 1;
        goto cleanup;
    }


    for(i=A->is;i<A->ie;i++)
    {
        local_i =
            i
            -
            A->is;


        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)2.0,
                A);

        if(err)
        {
            local_failed = 1;
            goto cleanup;
        }


        if(local_i==0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+1,
                    (LIS_SCALAR)-1.0,
                    A);
        }
        else if(local_i==1)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-1.0,
                    A);

            if(!err)
            {
                err =
                    lis_matrix_set_value(
                        LIS_INS_VALUE,
                        i,
                        i+1,
                        (LIS_SCALAR)-1.0,
                        A);
            }
        }
        else if(local_i==2)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-1.0,
                    A);

            if(!err)
            {
                err =
                    lis_matrix_set_value(
                        LIS_INS_VALUE,
                        i,
                        i+1,
                        (LIS_SCALAR)-1.0e-6,
                        A);
            }
        }
        else if(local_i==3)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-1.0e-6,
                    A);
        }
        else
        {
            err = LIS_ERR_ILL_ARG;
        }


        if(err)
        {
            local_failed = 1;
            goto cleanup;
        }
    }


    err =
        lis_matrix_assemble(
            A);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_precon_mfamg_get_distributed_aggregate_info(
            A,
            STRENGTH_N,
            &local_nc,
            &global_nc,
            &coarse_is,
            &coarse_ie,
            &map_n,
            fine_to_coarse);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    if(
        local_nc!=3
        ||
        global_nc!=6
        ||
        coarse_is!=3*rank
        ||
        coarse_ie!=3*(rank+1)
        ||
        map_n!=4)
    {
        local_failed = 1;
        goto cleanup;
    }


    for(i=0;i<4;i++)
    {
        if(i<2)
        {
            expected =
                coarse_is;
        }
        else if(i==2)
        {
            expected =
                coarse_is
                +
                1;
        }
        else
        {
            expected =
                coarse_is
                +
                2;
        }


        if(fine_to_coarse[i]!=expected)
        {
            local_failed = 1;
            goto cleanup;
        }
    }


cleanup:

    if(A)
    {
        lis_matrix_destroy(
            A);
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_failed =
        local_failed;

#endif


    if(rank==0)
    {
        printf(
            "MPI_STRENGTH_AGGREGATION=%s "
            "ranks=%d global_nc=%d theta=%.2e\n",
            global_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs,
            (int)global_nc,
            (double)0.05);
    }


    return global_failed;
}



static int
check_mfamg_distributed_prolongation(
    LIS_MATRIX A,
    int rank,
    int nprocs)
{
    LIS_VECTOR z;

    LIS_SCALAR fine_prolong[TEST_N];
    LIS_SCALAR coarse_z[TEST_N];

    LIS_SCALAR expected;

    LIS_REAL z0;
    LIS_REAL z1;
    LIS_REAL denom;
    LIS_REAL diff;

    LIS_INT map_n;
    LIS_INT local_nc;
    LIS_INT global_nc;

    LIS_INT i;
    LIS_INT a;
    LIS_INT global_fine;
    LIS_INT pair_first;
    LIS_INT err;

    int local_failed;
    int global_failed;


    z = NULL;

    map_n = -1;
    local_nc = -1;
    global_nc = -1;

    local_failed = 0;
    global_failed = 0;


    if(
        A==NULL
        ||
        nprocs!=2
        ||
        TEST_N!=24
        ||
        A->n!=12)
    {
        local_failed = 1;

        goto reduce_result;
    }


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
    {
        local_failed = 1;

        goto reduce_result;
    }


    for(i=0;i<A->n;i++)
    {
        z->value[i] =
            (LIS_SCALAR)(
                A->is
                +
                i
                +
                1);
    }


    err =
        lis_precon_mfamg_get_distributed_prolongation_info(
            A,
            z,
            TEST_N,
            TEST_N,
            &map_n,
            &local_nc,
            &global_nc,
            fine_prolong,
            coarse_z);


    if(err)
    {
        fprintf(
            stderr,
            "MPI prolongation query failed on rank %d rc=%d\n",
            rank,
            (int)err);

        local_failed = 1;

        goto reduce_result;
    }


    if(
        map_n!=A->np
        ||
        local_nc!=6
        ||
        global_nc!=12)
    {
        fprintf(
            stderr,
            "MPI prolongation metadata mismatch on rank %d: "
            "map_n=%d np=%d local_nc=%d global_nc=%d\n",
            rank,
            (int)map_n,
            (int)A->np,
            (int)local_nc,
            (int)global_nc);

        local_failed = 1;
    }


    for(i=0;i<A->np;i++)
    {
        if(i<A->n)
        {
            global_fine =
                A->is
                +
                i;
        }
        else
        {
            global_fine =
                A->l2g_map[
                    i
                    -
                    A->n];
        }


        pair_first =
            2
            *
            (
                global_fine
                /
                2
            );


        z0 =
            (LIS_REAL)(
                pair_first
                +
                1);

        z1 =
            (LIS_REAL)(
                pair_first
                +
                2);


        denom =
            sqrt(
                z0*z0
                +
                z1*z1);


        expected =
            (LIS_SCALAR)(
                (
                    (LIS_REAL)(
                        global_fine
                        +
                        1)
                )
                /
                denom);


        diff =
            fabs(
                fine_prolong[i]
                -
                expected);


        if(diff>MATCH_TOL)
        {
            fprintf(
                stderr,
                "MPI prolongation mismatch on rank %d "
                "fine=%d local=%d got=%.16e expected=%.16e diff=%.3e\n",
                rank,
                (int)global_fine,
                (int)i,
                (double)fine_prolong[i],
                (double)expected,
                (double)diff);

            local_failed = 1;
        }
    }


    for(a=0;a<local_nc;a++)
    {
        pair_first =
            A->is
            +
            2*a;


        z0 =
            (LIS_REAL)(
                pair_first
                +
                1);

        z1 =
            (LIS_REAL)(
                pair_first
                +
                2);


        denom =
            sqrt(
                z0*z0
                +
                z1*z1);


        expected =
            (LIS_SCALAR)denom;


        diff =
            fabs(
                coarse_z[a]
                -
                expected);


        if(diff>MATCH_TOL)
        {
            fprintf(
                stderr,
                "MPI coarse nullspace mismatch on rank %d "
                "aggregate=%d got=%.16e expected=%.16e diff=%.3e\n",
                rank,
                (int)a,
                (double)coarse_z[a],
                (double)expected,
                (double)diff);

            local_failed = 1;
        }
    }


reduce_result:

#ifdef USE_MPI

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_failed =
        local_failed;

#endif


    if(z)
    {
        lis_vector_destroy(
            z);
    }


    if(rank==0)
    {
        printf(
            "MPI_DISTRIBUTED_PROLONGATION=%s ranks=%d global_nc=%d\n",
            global_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs,
            (int)global_nc);
    }


    return global_failed;
}


static int
check_mfamg_distributed_hierarchy(
    LIS_MATRIX A,
    int rank,
    int nprocs)
{
    LIS_VECTOR z;

    LIS_INT local_n[8];
    LIS_INT global_n[8];
    LIS_INT owns_A[8];

    LIS_INT row[TEST_N+1];
    LIS_INT global_cols[3*TEST_N];
    LIS_SCALAR values[3*TEST_N];

    LIS_SCALAR terminal_dense[64];
    LIS_SCALAR terminal_solution[8];

    LIS_INT nlevels;
    LIS_INT terminal_np;
    LIS_INT terminal_nnz;

    LIS_INT global_terminal_nnz;

    LIS_INT i;
    LIS_INT p;
    LIS_INT global_row;
    LIS_INT col;

    LIS_INT expected_nnz;

    LIS_INT seen_left;
    LIS_INT seen_diag;
    LIS_INT seen_right;

    LIS_INT boundary_local;
    LIS_INT boundary_hits;

    LIS_SCALAR expected;
    LIS_REAL diff;

    LIS_INT err;

    int local_failed;
    int global_failed;


    z = NULL;

    nlevels = 0;
    terminal_np = 0;
    terminal_nnz = 0;

    global_terminal_nnz = 0;

    boundary_local = 0;
    boundary_hits = 0;

    local_failed = 0;
    global_failed = 0;


    memset(
        local_n,
        0,
        sizeof(local_n));

    memset(
        global_n,
        0,
        sizeof(global_n));

    memset(
        owns_A,
        0,
        sizeof(owns_A));


    if(
        A==NULL
        ||
        nprocs!=2
        ||
        TEST_N!=24
        ||
        A->n!=12
        ||
        A->gn!=24)
    {
        local_failed = 1;
        goto reduce_result;
    }


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
    {
        local_failed = 1;
        goto reduce_result;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    err =
        lis_precon_mfamg_get_distributed_hierarchy_info(
            A,
            z,
            8,
            TEST_N+1,
            3*TEST_N,
            &nlevels,
            local_n,
            global_n,
            owns_A,
            &terminal_np,
            &terminal_nnz,
            row,
            global_cols,
            values,
            64,
            8,
            terminal_dense,
            terminal_solution);

    if(err)
    {
        fprintf(
            stderr,
            "MPI hierarchy query failed on rank %d rc=%d\n",
            rank,
            (int)err);

        local_failed = 1;
        terminal_nnz = 0;

        goto reduce_result;
    }


    if(
        nlevels!=3
        ||
        local_n[0]!=12
        ||
        local_n[1]!=6
        ||
        local_n[2]!=3
        ||
        global_n[0]!=24
        ||
        global_n[1]!=12
        ||
        global_n[2]!=6
        ||
        owns_A[0]!=0
        ||
        owns_A[1]!=1
        ||
        owns_A[2]!=1)
    {
        fprintf(
            stderr,
            "MPI hierarchy metadata mismatch on rank %d: "
            "nlevels=%d local=%d,%d,%d global=%d,%d,%d "
            "owns=%d,%d,%d\n",
            rank,
            (int)nlevels,
            (int)local_n[0],
            (int)local_n[1],
            (int)local_n[2],
            (int)global_n[0],
            (int)global_n[1],
            (int)global_n[2],
            (int)owns_A[0],
            (int)owns_A[1],
            (int)owns_A[2]);

        local_failed = 1;
    }


    if(
        terminal_np!=4
        ||
        terminal_nnz!=8
        ||
        row[0]!=0
        ||
        row[3]!=8)
    {
        fprintf(
            stderr,
            "MPI terminal metadata mismatch on rank %d: "
            "np=%d nnz=%d row0=%d row3=%d\n",
            rank,
            (int)terminal_np,
            (int)terminal_nnz,
            (int)row[0],
            (int)row[3]);

        local_failed = 1;
    }


    if(
        nlevels==3
        &&
        local_n[2]==3
        &&
        global_n[2]==6)
    {
        for(i=0;i<3;i++)
        {
            global_row =
                3*rank
                +
                i;

            expected_nnz =
                (
                    global_row==0
                    ||
                    global_row==5
                )
                ?
                2
                :
                3;


            if(
                row[i+1]-row[i]
                !=
                expected_nnz)
            {
                fprintf(
                    stderr,
                    "MPI terminal row nnz mismatch rank=%d "
                    "row=%d got=%d expected=%d\n",
                    rank,
                    (int)global_row,
                    (int)(
                        row[i+1]
                        -
                        row[i]),
                    (int)expected_nnz);

                local_failed = 1;
            }


            seen_left = 0;
            seen_diag = 0;
            seen_right = 0;


            for(
                p=row[i];
                p<row[i+1];
                p++)
            {
                col =
                    global_cols[p];


                if(col==global_row)
                {
                    expected =
                        (LIS_SCALAR)0.5;

                    seen_diag++;
                }
                else if(
                    global_row>0
                    &&
                    col==global_row-1)
                {
                    expected =
                        (LIS_SCALAR)-0.25;

                    seen_left++;


                    if(
                        global_row==3
                        &&
                        col==2)
                    {
                        boundary_local = 1;
                    }
                }
                else if(
                    global_row<5
                    &&
                    col==global_row+1)
                {
                    expected =
                        (LIS_SCALAR)-0.25;

                    seen_right++;


                    if(
                        global_row==2
                        &&
                        col==3)
                    {
                        boundary_local = 1;
                    }
                }
                else
                {
                    fprintf(
                        stderr,
                        "MPI terminal unexpected column rank=%d "
                        "row=%d col=%d value=%.16e\n",
                        rank,
                        (int)global_row,
                        (int)col,
                        (double)values[p]);

                    local_failed = 1;
                    continue;
                }


                diff =
                    fabs(
                        values[p]
                        -
                        expected);


                if(diff>MATCH_TOL)
                {
                    fprintf(
                        stderr,
                        "MPI terminal value mismatch rank=%d "
                        "row=%d col=%d got=%.16e expected=%.16e "
                        "diff=%.3e\n",
                        rank,
                        (int)global_row,
                        (int)col,
                        (double)values[p],
                        (double)expected,
                        (double)diff);

                    local_failed = 1;
                }
            }


            if(seen_diag!=1)
                local_failed = 1;

            if(
                global_row>0
                &&
                seen_left!=1)
            {
                local_failed = 1;
            }

            if(
                global_row==0
                &&
                seen_left!=0)
            {
                local_failed = 1;
            }

            if(
                global_row<5
                &&
                seen_right!=1)
            {
                local_failed = 1;
            }

            if(
                global_row==5
                &&
                seen_right!=0)
            {
                local_failed = 1;
            }
        }
    }


    if(boundary_local!=1)
    {
        fprintf(
            stderr,
            "MPI terminal boundary coupling missing on rank %d\n",
            rank);

        local_failed = 1;
    }


    if(
        nlevels==3
        &&
        global_n[2]==6)
    {
        for(i=0;i<6;i++)
        {
            for(p=0;p<6;p++)
            {
                if(i==p)
                {
                    expected =
                        (LIS_SCALAR)0.5;
                }
                else if(
                    p==i-1
                    ||
                    p==i+1)
                {
                    expected =
                        (LIS_SCALAR)-0.25;
                }
                else
                {
                    expected =
                        (LIS_SCALAR)0.0;
                }


                diff =
                    fabs(
                        terminal_dense[
                            i*6+p]
                        -
                        expected);


                if(diff>MATCH_TOL)
                {
                    fprintf(
                        stderr,
                        "MPI replicated terminal matrix mismatch "
                        "rank=%d row=%d col=%d got=%.16e "
                        "expected=%.16e diff=%.3e\n",
                        rank,
                        (int)i,
                        (int)p,
                        (double)terminal_dense[
                            i*6+p],
                        (double)expected,
                        (double)diff);

                    local_failed = 1;
                }
            }


            expected =
                (LIS_SCALAR)(
                    i+1);


            diff =
                fabs(
                    terminal_solution[i]
                    -
                    expected);


            if(diff>MATCH_TOL)
            {
                fprintf(
                    stderr,
                    "MPI replicated terminal solve mismatch "
                    "rank=%d i=%d got=%.16e expected=%.16e "
                    "diff=%.3e\n",
                    rank,
                    (int)i,
                    (double)terminal_solution[i],
                    (double)expected,
                    (double)diff);

                local_failed = 1;
            }
        }
    }
    else
    {
        local_failed = 1;
    }


reduce_result:

#ifdef USE_MPI

    MPI_Allreduce(
        &terminal_nnz,
        &global_terminal_nnz,
        1,
        LIS_MPI_INT,
        MPI_SUM,
        LIS_COMM_WORLD);

    MPI_Allreduce(
        &boundary_local,
        &boundary_hits,
        1,
        LIS_MPI_INT,
        MPI_SUM,
        LIS_COMM_WORLD);

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_terminal_nnz =
        terminal_nnz;

    boundary_hits =
        boundary_local;

    global_failed =
        local_failed;

#endif


    if(
        global_terminal_nnz!=16
        ||
        boundary_hits!=2)
    {
        global_failed = 1;
    }


    if(rank==0)
    {
        printf(
            "MPI_DISTRIBUTED_HIERARCHY=%s "
            "ranks=%d levels=%d global=24,12,6 "
            "local=12,6,3 terminal_nnz=%d boundary_edges=%d\n",
            global_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs,
            (int)nlevels,
            (int)global_terminal_nnz,
            (int)boundary_hits);
    }


    if(rank==0)
    {
        printf(
            "MPI_REPLICATED_TERMINAL=%s "
            "ranks=%d global_n=6 dense_entries=36\n",
            global_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs);
    }


    if(z)
        lis_vector_destroy(z);


    return global_failed;
}


static int
check_mfamg_distributed_galerkin(
    LIS_MATRIX A,
    int rank,
    int nprocs)
{
    LIS_VECTOR z;

    LIS_INT row[TEST_N+1];
    LIS_INT global_cols[3*TEST_N];
    LIS_SCALAR values[3*TEST_N];

    LIS_INT local_nc;
    LIS_INT global_nc;
    LIS_INT coarse_is;
    LIS_INT coarse_ie;
    LIS_INT coarse_np;
    LIS_INT local_nnz;

    LIS_INT global_nnz;

    LIS_INT i;
    LIS_INT p;
    LIS_INT global_row;
    LIS_INT col;

    LIS_INT expected_nnz;

    LIS_INT seen_left;
    LIS_INT seen_diag;
    LIS_INT seen_right;

    LIS_INT boundary_local;
    LIS_INT boundary_hits;

    LIS_SCALAR expected;

    LIS_REAL diff;

    LIS_INT err;

    int local_failed;
    int global_failed;


    z = NULL;

    local_nc = -1;
    global_nc = -1;
    coarse_is = -1;
    coarse_ie = -1;
    coarse_np = -1;
    local_nnz = 0;

    global_nnz = 0;

    boundary_local = 0;
    boundary_hits = 0;

    local_failed = 0;
    global_failed = 0;


    if(
        A==NULL
        ||
        nprocs!=2
        ||
        TEST_N!=24
        ||
        A->n!=12)
    {
        local_failed = 1;

        goto reduce_result;
    }


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
    {
        local_failed = 1;

        goto reduce_result;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    err =
        lis_precon_mfamg_get_distributed_coarse_info(
            A,
            z,
            TEST_N+1,
            3*TEST_N,
            &local_nc,
            &global_nc,
            &coarse_is,
            &coarse_ie,
            &coarse_np,
            &local_nnz,
            row,
            global_cols,
            values);


    if(err)
    {
        fprintf(
            stderr,
            "MPI Galerkin query failed on rank %d rc=%d\n",
            rank,
            (int)err);

        local_failed = 1;

        local_nnz = 0;

        goto reduce_result;
    }


    if(
        local_nc!=6
        ||
        global_nc!=12
        ||
        coarse_is!=6*rank
        ||
        coarse_ie!=6*(rank+1)
        ||
        coarse_np!=7
        ||
        local_nnz!=17)
    {
        fprintf(
            stderr,
            "MPI Galerkin metadata mismatch on rank %d: "
            "local_nc=%d global_nc=%d coarse=[%d,%d) "
            "np=%d nnz=%d\n",
            rank,
            (int)local_nc,
            (int)global_nc,
            (int)coarse_is,
            (int)coarse_ie,
            (int)coarse_np,
            (int)local_nnz);

        local_failed = 1;
    }


    if(
        local_nc>0
        &&
        (
            row[0]!=0
            ||
            row[local_nc]!=local_nnz
        ))
    {
        fprintf(
            stderr,
            "MPI Galerkin row pointer mismatch on rank %d\n",
            rank);

        local_failed = 1;
    }


    if(local_nc==6)
    {
        for(i=0;i<local_nc;i++)
        {
            global_row =
                coarse_is
                +
                i;


            expected_nnz =
                (
                    global_row==0
                    ||
                    global_row==global_nc-1
                )
                ?
                2
                :
                3;


            if(
                row[i+1]-row[i]
                !=
                expected_nnz)
            {
                fprintf(
                    stderr,
                    "MPI Galerkin row nnz mismatch on rank %d "
                    "row=%d got=%d expected=%d\n",
                    rank,
                    (int)global_row,
                    (int)(
                        row[i+1]
                        -
                        row[i]),
                    (int)expected_nnz);

                local_failed = 1;
            }


            seen_left = 0;
            seen_diag = 0;
            seen_right = 0;


            for(
                p=row[i];
                p<row[i+1];
                p++)
            {
                col =
                    global_cols[p];


                if(col==global_row)
                {
                    expected =
                        (LIS_SCALAR)1.0;

                    seen_diag++;
                }
                else if(
                    global_row>0
                    &&
                    col==global_row-1)
                {
                    expected =
                        (LIS_SCALAR)-0.5;

                    seen_left++;


                    if(
                        global_row==6
                        &&
                        col==5)
                    {
                        boundary_local = 1;
                    }
                }
                else if(
                    global_row+1<global_nc
                    &&
                    col==global_row+1)
                {
                    expected =
                        (LIS_SCALAR)-0.5;

                    seen_right++;


                    if(
                        global_row==5
                        &&
                        col==6)
                    {
                        boundary_local = 1;
                    }
                }
                else
                {
                    fprintf(
                        stderr,
                        "MPI Galerkin unexpected column on rank %d "
                        "row=%d col=%d value=%.16e\n",
                        rank,
                        (int)global_row,
                        (int)col,
                        (double)values[p]);

                    local_failed = 1;

                    continue;
                }


                diff =
                    fabs(
                        values[p]
                        -
                        expected);


                if(diff>MATCH_TOL)
                {
                    fprintf(
                        stderr,
                        "MPI Galerkin value mismatch on rank %d "
                        "row=%d col=%d got=%.16e expected=%.16e "
                        "diff=%.3e\n",
                        rank,
                        (int)global_row,
                        (int)col,
                        (double)values[p],
                        (double)expected,
                        (double)diff);

                    local_failed = 1;
                }
            }


            if(seen_diag!=1)
            {
                fprintf(
                    stderr,
                    "MPI Galerkin diagonal count mismatch on rank %d "
                    "row=%d count=%d\n",
                    rank,
                    (int)global_row,
                    (int)seen_diag);

                local_failed = 1;
            }


            if(
                global_row>0
                &&
                seen_left!=1)
            {
                fprintf(
                    stderr,
                    "MPI Galerkin left coupling mismatch on rank %d "
                    "row=%d count=%d\n",
                    rank,
                    (int)global_row,
                    (int)seen_left);

                local_failed = 1;
            }


            if(
                global_row==0
                &&
                seen_left!=0)
            {
                local_failed = 1;
            }


            if(
                global_row+1<global_nc
                &&
                seen_right!=1)
            {
                fprintf(
                    stderr,
                    "MPI Galerkin right coupling mismatch on rank %d "
                    "row=%d count=%d\n",
                    rank,
                    (int)global_row,
                    (int)seen_right);

                local_failed = 1;
            }


            if(
                global_row==global_nc-1
                &&
                seen_right!=0)
            {
                local_failed = 1;
            }
        }
    }


    if(boundary_local!=1)
    {
        fprintf(
            stderr,
            "MPI Galerkin boundary coupling missing on rank %d\n",
            rank);

        local_failed = 1;
    }


reduce_result:

#ifdef USE_MPI

    MPI_Allreduce(
        &local_nnz,
        &global_nnz,
        1,
        LIS_MPI_INT,
        MPI_SUM,
        LIS_COMM_WORLD);

    MPI_Allreduce(
        &boundary_local,
        &boundary_hits,
        1,
        LIS_MPI_INT,
        MPI_SUM,
        LIS_COMM_WORLD);

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_nnz =
        local_nnz;

    boundary_hits =
        boundary_local;

    global_failed =
        local_failed;

#endif


    if(
        global_nnz!=34
        ||
        boundary_hits!=2)
    {
        global_failed = 1;
    }


    if(rank==0)
    {
        printf(
            "MPI_DISTRIBUTED_GALERKIN=%s "
            "ranks=%d global_nc=%d global_nnz=%d "
            "boundary_edges=%d\n",
            global_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs,
            (int)global_nc,
            (int)global_nnz,
            (int)boundary_hits);
    }


    if(z)
    {
        lis_vector_destroy(
            z);
    }


    return global_failed;
}


static int
check_mfamg_distributed_psolve(
    LIS_MATRIX A,
    int rank,
    int nprocs)
{
    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR ax;
    LIS_VECTOR z;

    LIS_VECTOR basis[1];

    LIS_SOLVER solver;
    LIS_PRECON precon;

    LIS_REAL local_rhs2;
    LIS_REAL global_rhs2;

    LIS_REAL local_residual2;
    LIS_REAL global_residual2;

    LIS_REAL local_correction2;
    LIS_REAL global_correction2;

    LIS_REAL local_remote2;
    LIS_REAL global_remote2;

    LIS_REAL residual_ratio;
    LIS_REAL correction_norm;
    LIS_REAL remote_norm;

    LIS_REAL mag;

    LIS_INT global_i;
    LIS_INT remote_start;
    LIS_INT i;
    LIS_INT err;
    LIS_INT destroy_err;

    int local_failed;
    int global_failed;


    b = NULL;
    x = NULL;
    ax = NULL;
    z = NULL;

    solver = NULL;
    precon = NULL;

    local_rhs2 = 0.0;
    global_rhs2 = 0.0;

    local_residual2 = 0.0;
    global_residual2 = 0.0;

    local_correction2 = 0.0;
    global_correction2 = 0.0;

    local_remote2 = 0.0;
    global_remote2 = 0.0;

    residual_ratio = -1.0;
    correction_norm = -1.0;
    remote_norm = -1.0;

    local_failed = 0;
    global_failed = 0;


    if(
        A==NULL
        ||
        nprocs<2
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs!=nprocs
        ||
        A->n<=0
        ||
        A->gn<=0
        ||
        A->ranges==NULL
        ||
        A->ranges[0]!=0
        ||
        A->ranges[1]<=0
        ||
        A->ranges[1]>=A->gn)
    {
        local_failed = 1;
        goto cleanup;
    }


    remote_start =
        A->ranges[1];


    err =
        lis_vector_duplicate(
            A,
            &b);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            A,
            &x);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            A,
            &ax);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        ax);

    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    /*
     * Unit impulse on the first global row.
     *
     * The right half of the global domain starts with a zero
     * right-hand side.  A nonzero correction there therefore
     * demonstrates that the applied V-cycle crosses the rank
     * boundary.
     */
    if(
        A->is<=0
        &&
        A->ie>0)
    {
        b->value[
            -A->is] =
            (LIS_SCALAR)1.0;
    }


    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_solver_set_option(
            "-i gmres "
            "-p mfamg "
            "-scale none "
            "-print none",
            solver);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_solver_set_matrix(
            A,
            solver);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_precon_create(
            solver,
            &precon);

    if(err)
    {
        precon = NULL;
        local_failed = 1;
        goto cleanup;
    }


    solver->precon =
        precon;


    err =
        lis_psolve(
            solver,
            b,
            x);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_matvec(
            A,
            x,
            ax);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    for(i=0;i<A->n;i++)
    {
        global_i =
            A->is
            +
            i;


        mag =
            fabs(
                b->value[i]);

        local_rhs2 +=
            mag
            *
            mag;


        mag =
            fabs(
                b->value[i]
                -
                ax->value[i]);

        local_residual2 +=
            mag
            *
            mag;


        mag =
            fabs(
                x->value[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            local_failed = 1;
        }


        local_correction2 +=
            mag
            *
            mag;


        if(
            global_i
            >=
            remote_start)
        {
            local_remote2 +=
                mag
                *
                mag;
        }
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_rhs2,
        &global_rhs2,
        1,
        LIS_MPI_REAL,
        MPI_SUM,
        A->comm);

    MPI_Allreduce(
        &local_residual2,
        &global_residual2,
        1,
        LIS_MPI_REAL,
        MPI_SUM,
        A->comm);

    MPI_Allreduce(
        &local_correction2,
        &global_correction2,
        1,
        LIS_MPI_REAL,
        MPI_SUM,
        A->comm);

    MPI_Allreduce(
        &local_remote2,
        &global_remote2,
        1,
        LIS_MPI_REAL,
        MPI_SUM,
        A->comm);

#else

    global_rhs2 =
        local_rhs2;

    global_residual2 =
        local_residual2;

    global_correction2 =
        local_correction2;

    global_remote2 =
        local_remote2;

#endif


    if(!(global_rhs2>0.0))
    {
        local_failed = 1;
    }
    else
    {
        residual_ratio =
            sqrt(
                global_residual2
                /
                global_rhs2);
    }


    correction_norm =
        sqrt(
            global_correction2);

    remote_norm =
        sqrt(
            global_remote2);


    if(
        residual_ratio!=residual_ratio
        ||
        !(residual_ratio<1.0)
        ||
        correction_norm!=correction_norm
        ||
        !(correction_norm>0.0)
        ||
        remote_norm!=remote_norm
        ||
        !(remote_norm>1.0e-12))
    {
        local_failed = 1;
    }


cleanup:

    if(precon)
    {
        if(solver)
        {
            solver->precon =
                NULL;
        }


        destroy_err =
            lis_precon_destroy(
                precon);

        if(destroy_err)
        {
            local_failed = 1;
        }


        precon = NULL;
    }


    if(solver)
    {
        lis_solver_destroy(
            solver);
    }


    if(z)
    {
        lis_vector_destroy(
            z);
    }

    if(ax)
    {
        lis_vector_destroy(
            ax);
    }

    if(x)
    {
        lis_vector_destroy(
            x);
    }

    if(b)
    {
        lis_vector_destroy(
            b);
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_failed =
        local_failed;

#endif


    if(rank==0)
    {
        if(
            A!=NULL
            &&
            A->gn==TEST_N)
        {
            printf(
                "MPI_DISTRIBUTED_PSOLVE=%s "
                "ranks=%d relres=%.12e "
                "correction=%.12e remote=%.12e\n",
                global_failed
                ?
                "FAILED"
                :
                "PASSED",
                nprocs,
                (double)residual_ratio,
                (double)correction_norm,
                (double)remote_norm);
        }
        else
        {
            printf(
                "MPI_HARDENING_PSOLVE=%s "
                "ranks=%d n=%d relres=%.12e "
                "correction=%.12e remote=%.12e\n",
                global_failed
                ?
                "FAILED"
                :
                "PASSED",
                nprocs,
                A
                ?
                (int)A->gn
                :
                -1,
                (double)residual_ratio,
                (double)correction_norm,
                (double)remote_norm);
        }
    }


    return global_failed;
}

static int
check_mfamg_distributed_multivector_psolve(
    LIS_MATRIX A,
    int rank,
    int nprocs)
{
    enum
    {
        MODE_COUNT = 2,
        LOCAL_N = 12,
        LOCAL_AGGREGATES = 6,
        EXPECTED_LOCAL_COARSE = 9,
        EXPECTED_GLOBAL_COARSE = 18,
        TERMINAL_LIMIT = 16
    };

    LIS_VECTOR b;
    LIS_VECTOR x;
    LIS_VECTOR ax;
    LIS_VECTOR z0;
    LIS_VECTOR z1;

    LIS_VECTOR basis[MODE_COUNT];

    LIS_INT aggregate[LOCAL_N];

    LIS_INT *aggregate_rank;
    LIS_INT *aggregate_offset;
    LIS_INT *row_ptr;

    LIS_SCALAR *prolong;
    LIS_SCALAR *coarse_basis;

    LIS_INT local_coarse_n;
    LIS_INT global_first_coarse;

    LIS_SOLVER solver;
    LIS_PRECON precon;

    LIS_REAL local_rhs2;
    LIS_REAL global_rhs2;

    LIS_REAL local_residual2;
    LIS_REAL global_residual2;

    LIS_REAL local_correction2;
    LIS_REAL global_correction2;

    LIS_REAL local_remote2;
    LIS_REAL global_remote2;

    LIS_REAL residual_ratio;
    LIS_REAL correction_norm;
    LIS_REAL remote_norm;

    LIS_REAL mag;

    LIS_INT global_i;
    LIS_INT remote_start;
    LIS_INT i;
    LIS_INT a;
    LIS_INT err;
    LIS_INT destroy_err;

    int local_failed;
    int global_failed;


    b = NULL;
    x = NULL;
    ax = NULL;
    z0 = NULL;
    z1 = NULL;

    aggregate_rank = NULL;
    aggregate_offset = NULL;
    row_ptr = NULL;

    prolong = NULL;
    coarse_basis = NULL;

    local_coarse_n = -1;
    global_first_coarse = -1;

    solver = NULL;
    precon = NULL;

    local_rhs2 = 0.0;
    global_rhs2 = 0.0;

    local_residual2 = 0.0;
    global_residual2 = 0.0;

    local_correction2 = 0.0;
    global_correction2 = 0.0;

    local_remote2 = 0.0;
    global_remote2 = 0.0;

    residual_ratio = -1.0;
    correction_norm = -1.0;
    remote_norm = -1.0;

    local_failed = 0;
    global_failed = 0;


    if(
        A==NULL
        ||
        nprocs!=2
        ||
        A->matrix_type!=LIS_MATRIX_CSR
        ||
        A->nprocs!=nprocs
        ||
        A->n!=LOCAL_N
        ||
        A->gn!=TEST_N
        ||
        A->ranges==NULL
        ||
        A->ranges[0]!=0
        ||
        A->ranges[1]!=LOCAL_N
        ||
        EXPECTED_GLOBAL_COARSE<=TERMINAL_LIMIT)
    {
        local_failed = 1;
        goto cleanup;
    }


    remote_start =
        A->ranges[1];


    err =
        lis_vector_duplicate(
            A,
            &b);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            A,
            &x);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            A,
            &ax);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            A,
            &z0);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_vector_duplicate(
            A,
            &z1);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        ax);

    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z0);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        z1);


    for(i=0;i<LOCAL_N;i++)
    {
        aggregate[i] =
            i/2;

        if(
            aggregate[i]<3
            &&
            (i%2)!=0)
        {
            z1->value[i] =
                (LIS_SCALAR)1.0;
        }
    }


    basis[0] =
        z0;

    basis[1] =
        z1;


    err =
        lis_precon_mfamg_build_multivector_basis(
            LOCAL_N,
            LOCAL_AGGREGATES,
            aggregate,
            MODE_COUNT,
            basis,
            &aggregate_rank,
            &aggregate_offset,
            &row_ptr,
            &prolong,
            &coarse_basis,
            &local_coarse_n);

    if(err)
    {
        local_failed = 1;
    }


    if(
        !err
        &&
        (
            aggregate_rank==NULL
            ||
            aggregate_offset==NULL
            ||
            row_ptr==NULL
            ||
            prolong==NULL
            ||
            coarse_basis==NULL
            ||
            local_coarse_n!=EXPECTED_LOCAL_COARSE
            ||
            aggregate_offset[LOCAL_AGGREGATES]
                !=
                EXPECTED_LOCAL_COARSE
            ||
            row_ptr[LOCAL_N]!=18
        ))
    {
        local_failed = 1;
    }


    if(
        !err
        &&
        aggregate_rank!=NULL)
    {
        for(a=0;a<LOCAL_AGGREGATES;a++)
        {
            if(
                aggregate_rank[a]
                !=
                (
                    a<3
                    ?
                    2
                    :
                    1
                ))
            {
                local_failed = 1;
            }
        }
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        A->comm);

#else

    global_failed =
        local_failed;

#endif


    if(global_failed)
        goto cleanup;


#ifdef USE_MPI

    MPI_Allreduce(
        &local_coarse_n,
        &global_first_coarse,
        1,
        LIS_MPI_INT,
        MPI_SUM,
        A->comm);

#else

    global_first_coarse =
        local_coarse_n;

#endif


    if(
        global_first_coarse
        !=
        EXPECTED_GLOBAL_COARSE)
    {
        local_failed = 1;
        goto cleanup;
    }


    if(coarse_basis)
    {
        lis_free(coarse_basis);
        coarse_basis = NULL;
    }

    if(prolong)
    {
        lis_free(prolong);
        prolong = NULL;
    }

    if(row_ptr)
    {
        lis_free(row_ptr);
        row_ptr = NULL;
    }

    if(aggregate_offset)
    {
        lis_free(aggregate_offset);
        aggregate_offset = NULL;
    }

    if(aggregate_rank)
    {
        lis_free(aggregate_rank);
        aggregate_rank = NULL;
    }


    /*
     * Unit impulse on the first global row.
     *
     * The right half of the global domain starts with a zero
     * right-hand side.  A nonzero correction there therefore
     * demonstrates that the applied V-cycle crosses the rank
     * boundary.
     */
    if(
        A->is<=0
        &&
        A->ie>0)
    {
        b->value[
            -A->is] =
            (LIS_SCALAR)1.0;
    }


    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_solver_set_option(
            "-i gmres "
            "-p mfamg "
            "-scale none "
            "-print none",
            solver);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    basis[0] =
        z0;

    basis[1] =
        z1;


    err =
        lis_solver_set_near_nullspace(
            solver,
            MODE_COUNT,
            basis);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_solver_set_matrix(
            A,
            solver);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_precon_create(
            solver,
            &precon);

    if(rank==0)
    {
        printf(
            "MPI_MULTIVECTOR_CREATE=%s "
            "ranks=%d first_global_coarse=%d terminal_limit=%d\n",
            err
            ?
            "FAILED"
            :
            "PASSED",
            nprocs,
            (int)global_first_coarse,
            (int)TERMINAL_LIMIT);
    }


    if(err)
    {
        precon = NULL;
        local_failed = 1;
        goto cleanup;
    }


    solver->precon =
        precon;


    err =
        lis_psolve(
            solver,
            b,
            x);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_matvec(
            A,
            x,
            ax);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    for(i=0;i<A->n;i++)
    {
        global_i =
            A->is
            +
            i;


        mag =
            fabs(
                b->value[i]);

        local_rhs2 +=
            mag
            *
            mag;


        mag =
            fabs(
                b->value[i]
                -
                ax->value[i]);

        local_residual2 +=
            mag
            *
            mag;


        mag =
            fabs(
                x->value[i]);

        if(
            mag!=mag
            ||
            mag>LIS_SCALAR_MAX)
        {
            local_failed = 1;
        }


        local_correction2 +=
            mag
            *
            mag;


        if(
            global_i
            >=
            remote_start)
        {
            local_remote2 +=
                mag
                *
                mag;
        }
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_rhs2,
        &global_rhs2,
        1,
        LIS_MPI_REAL,
        MPI_SUM,
        A->comm);

    MPI_Allreduce(
        &local_residual2,
        &global_residual2,
        1,
        LIS_MPI_REAL,
        MPI_SUM,
        A->comm);

    MPI_Allreduce(
        &local_correction2,
        &global_correction2,
        1,
        LIS_MPI_REAL,
        MPI_SUM,
        A->comm);

    MPI_Allreduce(
        &local_remote2,
        &global_remote2,
        1,
        LIS_MPI_REAL,
        MPI_SUM,
        A->comm);

#else

    global_rhs2 =
        local_rhs2;

    global_residual2 =
        local_residual2;

    global_correction2 =
        local_correction2;

    global_remote2 =
        local_remote2;

#endif


    if(!(global_rhs2>0.0))
    {
        local_failed = 1;
    }
    else
    {
        residual_ratio =
            sqrt(
                global_residual2
                /
                global_rhs2);
    }


    correction_norm =
        sqrt(
            global_correction2);

    remote_norm =
        sqrt(
            global_remote2);


    if(
        residual_ratio!=residual_ratio
        ||
        !(residual_ratio<1.0)
        ||
        correction_norm!=correction_norm
        ||
        !(correction_norm>0.0)
        ||
        remote_norm!=remote_norm
        ||
        !(remote_norm>1.0e-12))
    {
        local_failed = 1;
    }


cleanup:

    if(precon)
    {
        if(solver)
        {
            solver->precon =
                NULL;
        }


        destroy_err =
            lis_precon_destroy(
                precon);

        if(destroy_err)
        {
            local_failed = 1;
        }


        precon = NULL;
    }


    if(solver)
    {
        lis_solver_destroy(
            solver);
    }


    if(coarse_basis)
        lis_free(coarse_basis);

    if(prolong)
        lis_free(prolong);

    if(row_ptr)
        lis_free(row_ptr);

    if(aggregate_offset)
        lis_free(aggregate_offset);

    if(aggregate_rank)
        lis_free(aggregate_rank);


    if(z1)
    {
        lis_vector_destroy(
            z1);
    }

    if(z0)
    {
        lis_vector_destroy(
            z0);
    }

    if(ax)
    {
        lis_vector_destroy(
            ax);
    }

    if(x)
    {
        lis_vector_destroy(
            x);
    }

    if(b)
    {
        lis_vector_destroy(
            b);
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_failed =
        local_failed;

#endif


    if(rank==0)
    {
        printf(
            "MPI_MULTIVECTOR_PSOLVE=%s "
            "ranks=%d relres=%.12e "
            "correction=%.12e remote=%.12e\n",
            global_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs,
            (double)residual_ratio,
            (double)correction_norm,
            (double)remote_norm);
    }


    return global_failed;
}


static int
check_mfamg_hardening_case(
    int rank,
    int nprocs)
{
    LIS_MATRIX A;

    LIS_INT i;
    LIS_INT err;

    LIS_INT local_n;
    LIS_INT min_local_n;
    LIS_INT max_local_n;

    int local_failed;
    int global_layout_failed;
    int psolve_failed;
    int global_failed;


    A = NULL;

    local_n = 0;
    min_local_n = 0;
    max_local_n = 0;

    local_failed = 0;
    global_layout_failed = 0;
    psolve_failed = 0;
    global_failed = 0;


    if(nprocs<2)
    {
        return 0;
    }


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    err =
        lis_matrix_set_size(
            A,
            0,
            HARD_N);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    for(i=A->is;i<A->ie;i++)
    {
        if(i>0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-1.0,
                    A);

            if(err)
            {
                local_failed = 1;
                goto cleanup;
            }
        }


        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)2.0,
                A);

        if(err)
        {
            local_failed = 1;
            goto cleanup;
        }


        if(i+1<HARD_N)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+1,
                    (LIS_SCALAR)-1.0,
                    A);

            if(err)
            {
                local_failed = 1;
                goto cleanup;
            }
        }
    }


    err =
        lis_matrix_assemble(
            A);

    if(err)
    {
        local_failed = 1;
        goto cleanup;
    }


    local_n =
        A->n;


#ifdef USE_MPI

    MPI_Allreduce(
        &local_n,
        &min_local_n,
        1,
        LIS_MPI_INT,
        MPI_MIN,
        A->comm);

    MPI_Allreduce(
        &local_n,
        &max_local_n,
        1,
        LIS_MPI_INT,
        MPI_MAX,
        A->comm);

#else

    min_local_n =
        local_n;

    max_local_n =
        local_n;

#endif


    if(
        A->gn!=HARD_N
        ||
        A->n<=0
        ||
        A->ranges==NULL
        ||
        A->is<0
        ||
        A->ie>A->gn
        ||
        A->ie-A->is!=A->n
        ||
        min_local_n<=0
        ||
        max_local_n<min_local_n
        ||
        max_local_n-min_local_n>1)
    {
        local_failed = 1;
    }


    if(
        HARD_N % nprocs != 0
        &&
        min_local_n==max_local_n)
    {
        local_failed = 1;
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_failed,
        &global_layout_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_layout_failed =
        local_failed;

#endif


    if(rank==0)
    {
        printf(
            "MPI_HARDENING_LAYOUT=%s "
            "ranks=%d n=%d local_min=%d local_max=%d\n",
            global_layout_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs,
            HARD_N,
            (int)min_local_n,
            (int)max_local_n);
    }


    if(global_layout_failed)
    {
        local_failed = 1;
        goto cleanup;
    }


    psolve_failed =
        check_mfamg_distributed_psolve(
            A,
            rank,
            nprocs);

    if(psolve_failed)
    {
        local_failed = 1;
    }


cleanup:

    if(A)
    {
        lis_matrix_destroy(
            A);
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &local_failed,
        &global_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_failed =
        local_failed;

#endif


    return global_failed;
}



static int
run_mpi_rejection(
    int rank,
    int nprocs)
{
    LIS_MATRIX A;
    LIS_MATRIX U;

    USER_CTX user_ctx;

    LIS_VECTOR z;
    LIS_VECTOR basis[1];

    LIS_SOLVER solver;
    LIS_PRECON precon;

    LIS_INT i;
    LIS_INT err;

    int row_failed;
    int global_row_failed;
    int local_failed;
    int failed;


    A = NULL;
    U = NULL;
    z = NULL;
    solver = NULL;
    precon = NULL;

    row_failed = 0;
    global_row_failed = 0;
    failed = 0;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    err =
        lis_matrix_set_size(
            A,
            0,
            TEST_N);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    for(i=A->is;i<A->ie;i++)
    {
        if(i>0)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-1.0,
                    A);

            if(err)
            {
                failed = 1;

                goto cleanup;
            }
        }


        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)2.0,
                A);

        if(err)
        {
            failed = 1;

            goto cleanup;
        }


        if(i+1<TEST_N)
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+1,
                    (LIS_SCALAR)-1.0,
                    A);

            if(err)
            {
                failed = 1;

                goto cleanup;
            }
        }
    }


    err =
        lis_matrix_assemble(
            A);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    if(
        A->np<=A->n
        ||
        A->l2g_map==NULL)
    {
        fprintf(
            stderr,
            "CSR MPI halo was not created on rank %d: n=%d np=%d\n",
            rank,
            (int)A->n,
            (int)A->np);

        row_failed = 1;
    }


    err =
        create_user(
            &user_ctx,
            &U);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    row_failed |=
        check_mfamg_global_rows(
            "CSR",
            A,
            rank);


    row_failed |=
        check_mfamg_global_rows(
            "USER",
            U,
            rank);


    if(nprocs==2)
    {
        row_failed |=
            check_mfamg_distributed_aggregates(
                A,
                rank,
                nprocs);

        row_failed |=
            check_mfamg_strength_aggregation(
                rank,
                nprocs);


        row_failed |=
            check_mfamg_distributed_prolongation(
                A,
                rank,
                nprocs);


        row_failed |=
            check_mfamg_distributed_hierarchy(
                A,
                rank,
                nprocs);


        row_failed |=
            check_mfamg_distributed_galerkin(
                A,
                rank,
                nprocs);
    }


#ifdef USE_MPI

    MPI_Allreduce(
        &row_failed,
        &global_row_failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    global_row_failed =
        row_failed;

#endif


    if(rank==0)
    {
        printf(
            "MPI_GLOBAL_ROW_ACCESS=%s ranks=%d\n",
            global_row_failed
            ?
            "FAILED"
            :
            "PASSED",
            nprocs);
    }


    if(global_row_failed)
    {
        failed = 1;
    }


    failed |=
        check_mfamg_distributed_psolve(
            A,
            rank,
            nprocs);


    if(nprocs==2)
    {
        failed |=
            check_mfamg_distributed_multivector_psolve(
                A,
                rank,
                nprocs);
    }


    failed |=
        check_mfamg_hardening_case(
            rank,
            nprocs);


    err =
        lis_vector_duplicate(
            A,
            &z);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        z);


    err =
        lis_solver_create(
            &solver);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    err =
        lis_solver_set_option(
            "-i gmres "
            "-p mfamg "
            "-scale none "
            "-print none",
            solver);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    basis[0] =
        z;


    err =
        lis_solver_set_near_nullspace(
            solver,
            1,
            basis);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    err =
        lis_solver_set_matrix(
            U,
            solver);

    if(err)
    {
        failed = 1;

        goto cleanup;
    }


    err =
        lis_precon_create(
            solver,
            &precon);


    if(err)
    {
        precon = NULL;
    }


    if(rank==0)
    {
        printf(
            "MPI_USER_REJECT rc=%d expected=%d ranks=%d\n",
            (int)err,
            (int)LIS_ERR_NOT_IMPLEMENTED,
            nprocs);
    }


    if(
        err
        !=
        LIS_ERR_NOT_IMPLEMENTED)
    {
        failed = 1;
    }


#ifdef USE_MPI

    local_failed =
        failed;

    MPI_Allreduce(
        &local_failed,
        &failed,
        1,
        MPI_INT,
        MPI_MAX,
        LIS_COMM_WORLD);

#else

    local_failed =
        failed;

    (void)local_failed;

#endif


cleanup:

    if(precon)
    {
        lis_precon_destroy(
            precon);
    }

    if(solver)
    {
        lis_solver_destroy(
            solver);
    }

    if(z)
    {
        lis_vector_destroy(
            z);
    }

    if(U)
    {
        lis_matrix_destroy(
            U);
    }

    if(A)
    {
        lis_matrix_destroy(
            A);
    }


    if(rank==0)
    {
        if(failed)
        {
            fprintf(
                stderr,
                "MFAMG MPI DISTRIBUTED TEST FAILED\n");
        }
        else
        {
            printf(
                "MFAMG MPI DISTRIBUTED TEST PASSED\n");
        }
    }


    return
        failed
        ?
        1
        :
        0;
}

/* ============================================================
 * Main
 * ============================================================ */

int
main(
    int argc,
    char *argv[])
{
    LIS_INT err;

    int rank;
    int nprocs;
    int failed;


    setvbuf(
        stdout,
        NULL,
        _IONBF,
        0);

    setvbuf(
        stderr,
        NULL,
        _IONBF,
        0);


    err =
        lis_initialize(
            &argc,
            &argv);

    if(err)
    {
        return 1;
    }


#ifdef USE_MPI

    MPI_Comm_rank(
        LIS_COMM_WORLD,
        &rank);

    MPI_Comm_size(
        LIS_COMM_WORLD,
        &nprocs);

#else

    rank = 0;
    nprocs = 1;

#endif


    if(nprocs==1)
    {
        failed =
            run_serial();
    }
    else
    {
        failed =
            run_mpi_rejection(
                rank,
                nprocs);
    }


    lis_finalize();


    return
        failed
        ?
        1
        :
        0;
}
