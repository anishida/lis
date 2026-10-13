#include <math.h>
#include <stdio.h>

#include "lislib.h"


extern LIS_INT
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


#define TEST_N 24
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


/* ============================================================
 * Serial contract
 * ============================================================ */

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


    row_failed |=
        check_mfamg_distributed_aggregates(
            A,
            rank,
            nprocs);


    row_failed |=
        check_mfamg_distributed_prolongation(
            A,
            rank,
            nprocs);


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
        precon = NULL;
    }


    if(rank==0)
    {
        printf(
            "MPI_REJECT rc=%d expected=%d ranks=%d\n",
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
                "MFAMG MPI ROW/REJECTION TEST FAILED\n");
        }
        else
        {
            printf(
                "MFAMG MPI ROW/REJECTION TEST PASSED\n");
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
