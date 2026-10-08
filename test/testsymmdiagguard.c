#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include "lis.h"

#include <math.h>
#include <stdio.h>

#define TEST_N 8


typedef struct
{
    LIS_INT n;
    LIS_INT is;
}
USER_ZERO;


/* ---------------------------------------------------------- */

static LIS_SCALAR
user_diag_value(
    USER_ZERO *ctx,
    LIS_INT local_i)
{
    LIS_INT global_i =
        ctx->is + local_i;

    if(global_i==0)
        return (LIS_SCALAR)0.0;

    return (LIS_SCALAR)4.0;
}


/* ---------------------------------------------------------- */

static LIS_INT
user_matvec(
    void *data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
    USER_ZERO *ctx =
        (USER_ZERO *)data;

    LIS_INT i;


    for(i=0;i<ctx->n;i++)
    {
        y[i] =
            user_diag_value(
                ctx,
                i)
            * x[i];
    }


    return LIS_SUCCESS;
}


/* ---------------------------------------------------------- */

static LIS_INT
user_matvech(
    void *data,
    const LIS_SCALAR *x,
    LIS_SCALAR *y)
{
#ifdef _COMPLEX

    USER_ZERO *ctx =
        (USER_ZERO *)data;

    LIS_INT i;


    for(i=0;i<ctx->n;i++)
    {
        y[i] =
            conj(
                user_diag_value(
                    ctx,
                    i))
            * x[i];
    }


    return LIS_SUCCESS;

#else

    return user_matvec(
        data,
        x,
        y);

#endif
}


/* ---------------------------------------------------------- */

static LIS_INT
user_get_diagonal(
    void *data,
    LIS_SCALAR *d)
{
    USER_ZERO *ctx =
        (USER_ZERO *)data;

    LIS_INT i;


    for(i=0;i<ctx->n;i++)
    {
        d[i] =
            user_diag_value(
                ctx,
                i);
    }


    return LIS_SUCCESS;
}


/* ---------------------------------------------------------- */

static LIS_INT
create_stored_diagonal(
    LIS_SCALAR first_diag,
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;

    LIS_INT is;
    LIS_INT ie;
    LIS_INT i;
    LIS_INT err;


    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if(err)
        return err;


    err = lis_matrix_set_size(
        A,
        0,
        TEST_N);

    if(err)
        goto fail;


    err = lis_matrix_set_type(
        A,
        LIS_MATRIX_CSR);

    if(err)
        goto fail;


    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if(err)
        goto fail;


    for(i=is;i<ie;i++)
    {
        LIS_SCALAR value =
            i==0
            ?
            first_diag
            :
            (LIS_SCALAR)4.0;


        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            value,
            A);

        if(err)
            goto fail;
    }


    err = lis_matrix_assemble(
        A);

    if(err)
        goto fail;


    *Aout = A;

    return LIS_SUCCESS;


fail:

    if(A)
        lis_matrix_destroy(A);

    return err;
}


/* ---------------------------------------------------------- */
/* Public lis_matrix_scale() must reject an exact zero.       */
/* ---------------------------------------------------------- */

static int
test_direct_zero(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR d = NULL;

    LIS_INT err;
    LIS_INT rc;

    int failed = 0;


    err = create_stored_diagonal(
        (LIS_SCALAR)0.0,
        &A);

    if(err)
    {
        printf(
            "DIRECT_ZERO setup=%d\n",
            (int)err);

        return 1;
    }


    err = lis_vector_duplicate(
        A,
        &b);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(
        A,
        &d);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        b);


    rc = lis_matrix_scale(
        A,
        b,
        d,
        LIS_SCALE_SYMM_DIAG);


    printf(
        "DIRECT_ZERO rc=%d expected=%d\n",
        (int)rc,
        (int)LIS_BREAKDOWN);


    if(rc!=LIS_BREAKDOWN)
    {
        fprintf(
            stderr,
            "FAIL direct zero diagonal was not rejected\n");

        failed = 1;
    }


cleanup:

    if(d)
        lis_vector_destroy(d);

    if(b)
        lis_vector_destroy(b);

    if(A)
        lis_matrix_destroy(A);


    return failed;
}


/* ---------------------------------------------------------- */
/* Very small but nonzero diagonal must remain legal.         */
/* ---------------------------------------------------------- */

static int
test_direct_tiny(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR d = NULL;

    LIS_INT is;
    LIS_INT ie;
    LIS_INT i;
    LIS_INT err;
    LIS_INT rc;

    int finite = 1;
    int saw_first = 0;
    int failed = 0;


    err = create_stored_diagonal(
        (LIS_SCALAR)1.0e-300,
        &A);

    if(err)
        return 1;


    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(
        A,
        &b);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_vector_duplicate(
        A,
        &d);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        b);


    rc = lis_matrix_scale(
        A,
        b,
        d,
        LIS_SCALE_SYMM_DIAG);


    if(rc==LIS_SUCCESS)
    {
        for(i=0;i<ie-is;i++)
        {
            LIS_REAL mag =
                fabs(d->value[i]);


            if(
                mag==0.0
                ||
                mag!=mag
                ||
                mag>LIS_SCALAR_MAX)
            {
                finite = 0;
            }


            if(is+i==0)
            {
                saw_first = 1;

                if(
                    mag<(LIS_REAL)1.0e149
                    ||
                    mag>(LIS_REAL)1.0e151)
                {
                    finite = 0;
                }
            }
        }
    }


    printf(
        "DIRECT_TINY rc=%d finite=%d saw_first=%d\n",
        (int)rc,
        finite,
        saw_first);


    if(
        rc!=LIS_SUCCESS
        ||
        !finite)
    {
        fprintf(
            stderr,
            "FAIL tiny nonzero diagonal was rejected or became non-finite\n");

        failed = 1;
    }


cleanup:

    if(d)
        lis_vector_destroy(d);

    if(b)
        lis_vector_destroy(b);

    if(A)
        lis_matrix_destroy(A);


    return failed;
}


/* ---------------------------------------------------------- */
/* Solver path for an explicitly stored matrix.               */
/* ---------------------------------------------------------- */

static int
test_solver_stored_zero(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT api;

    int failed = 0;


    err = create_stored_diagonal(
        (LIS_SCALAR)0.0,
        &A);

    if(err)
        return 1;


    if(
        lis_vector_duplicate(A,&b)
        ||
        lis_vector_duplicate(A,&x))
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    err = lis_solver_create(
        &solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_solver_set_option(
        "-i cg "
        "-p none "
        "-scale symm_diag "
        "-print none "
        "-maxiter 20",
        solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    api = lis_solve(
        A,
        b,
        x,
        solver);


    printf(
        "STORED_ZERO api=%d retcode=%d expected=%d\n",
        (int)api,
        (int)solver->retcode,
        (int)LIS_BREAKDOWN);


    if(
        api!=LIS_BREAKDOWN
        ||
        solver->retcode!=LIS_BREAKDOWN)
    {
        fprintf(
            stderr,
            "FAIL stored solver did not propagate scaling breakdown\n");

        failed = 1;
    }


cleanup:

    if(solver)
        lis_solver_destroy(solver);

    if(x)
        lis_vector_destroy(x);

    if(b)
        lis_vector_destroy(b);

    if(A)
        lis_matrix_destroy(A);


    return failed;
}


/* ---------------------------------------------------------- */
/* Matrix-free USER path must use the same guard.             */
/* ---------------------------------------------------------- */

static int
test_solver_user_zero(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;

    LIS_SOLVER solver = NULL;

    USER_ZERO ctx;

    LIS_INT is;
    LIS_INT ie;
    LIS_INT err;
    LIS_INT api;

    int failed = 0;


    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if(err)
        return 1;


    err = lis_matrix_set_size(
        A,
        0,
        TEST_N);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_matrix_get_range(
        A,
        &is,
        &ie);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    ctx.n  = ie-is;
    ctx.is = is;


    err = lis_matrix_set_user(
        A,
        &ctx,
        user_matvec,
        user_matvech);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_matrix_set_user_diagonal(
        A,
        user_get_diagonal);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    if(
        lis_vector_duplicate(A,&b)
        ||
        lis_vector_duplicate(A,&x))
    {
        failed = 1;
        goto cleanup;
    }


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    err = lis_solver_create(
        &solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    err = lis_solver_set_option(
        "-i cg "
        "-p none "
        "-scale symm_diag "
        "-print none "
        "-maxiter 20",
        solver);

    if(err)
    {
        failed = 1;
        goto cleanup;
    }


    api = lis_solve(
        A,
        b,
        x,
        solver);


    printf(
        "USER_ZERO api=%d retcode=%d expected=%d\n",
        (int)api,
        (int)solver->retcode,
        (int)LIS_BREAKDOWN);


    if(
        api!=LIS_BREAKDOWN
        ||
        solver->retcode!=LIS_BREAKDOWN)
    {
        fprintf(
            stderr,
            "FAIL USER solver did not propagate scaling breakdown\n");

        failed = 1;
    }


cleanup:

    if(solver)
        lis_solver_destroy(solver);

    if(x)
        lis_vector_destroy(x);

    if(b)
        lis_vector_destroy(b);

    if(A)
        lis_matrix_destroy(A);


    return failed;
}


/* ========================================================== */

int
main(
    int argc,
    char **argv)
{
    LIS_INT err;

    int failed = 0;


    err = lis_initialize(
        &argc,
        &argv);

    if(err)
        return 1;


    failed |=
        test_direct_zero();

    failed |=
        test_direct_tiny();

    failed |=
        test_solver_stored_zero();

    failed |=
        test_solver_user_zero();


    printf(
        "LIS_SYMM_DIAG_GUARD_STAGE6A %s\n",
        failed ? "FAILED" : "PASSED");


    lis_finalize();


    return failed ? 1 : 0;
}
