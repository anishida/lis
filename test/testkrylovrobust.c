#include "lis.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#define N 64
#define TEST_TOL 1.0e-10


typedef struct
{
    LIS_INT api;
    LIS_INT status;
    LIS_INT iter;

    LIS_REAL reported;
    LIS_REAL true_res;
    LIS_REAL xnorm;
}
RESULT;


/* ---------------------------------------------------------- */
/* Shifted 1-D Neumann operator                               */
/* ---------------------------------------------------------- */

static LIS_INT
build_neumann(
    LIS_REAL eps,
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;

    LIS_INT i;
    LIS_INT err;

    LIS_SCALAR diag;


    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if(err)
        return err;


    err = lis_matrix_set_size(
        A,
        0,
        N);

    if(err)
    {
        lis_matrix_destroy(A);
        return err;
    }


    for(i=0; i<N; i++)
    {
        diag =
            ((i==0 || i==N-1) ? 1.0 : 2.0)
            + eps;


        err = lis_matrix_set_value(
            LIS_INS_VALUE,
            i,
            i,
            diag,
            A);

        if(err)
        {
            lis_matrix_destroy(A);
            return err;
        }


        if(i>0)
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i-1,
                -1.0,
                A);

            if(err)
            {
                lis_matrix_destroy(A);
                return err;
            }
        }


        if(i<N-1)
        {
            err = lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i+1,
                -1.0,
                A);

            if(err)
            {
                lis_matrix_destroy(A);
                return err;
            }
        }
    }


    err = lis_matrix_assemble(A);

    if(err)
    {
        lis_matrix_destroy(A);
        return err;
    }


    *Aout = A;

    return LIS_SUCCESS;
}


/* ---------------------------------------------------------- */

static void
fill_rigid(
    LIS_VECTOR b)
{
    LIS_REAL value;

    value =
        1.0
        / sqrt((LIS_REAL)N);

    lis_vector_set_all(
        (LIS_SCALAR)value,
        b);
}


/* ---------------------------------------------------------- */

static LIS_REAL
true_relative_residual(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_VECTOR work)
{
    LIS_REAL rn = 0.0;
    LIS_REAL bn = 0.0;


    lis_matvec(
        A,
        x,
        work);

    lis_vector_axpy(
        (LIS_SCALAR)-1.0,
        b,
        work);

    lis_vector_nrm2(
        work,
        &rn);

    lis_vector_nrm2(
        b,
        &bn);


    if(bn!=0.0)
        return rn / bn;


    return rn;
}


/* ---------------------------------------------------------- */

static int
result_is_finite(
    const RESULT *r)
{
    return
        isfinite((double)r->reported)
        &&
        isfinite((double)r->true_res)
        &&
        isfinite((double)r->xnorm);
}


/* ---------------------------------------------------------- */
/* Generic Neumann solve                                      */
/* ---------------------------------------------------------- */

static LIS_INT
run_neumann(
    const char *solver_name,
    LIS_REAL eps,
    LIS_INT maxiter,
    LIS_INT restart,
    RESULT *result)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR work = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;

    char option[512];


    memset(
        result,
        0,
        sizeof(*result));

    result->status = -999;
    result->iter = -999;
    result->reported = -1.0;
    result->true_res = -1.0;
    result->xnorm = -1.0;


    err = build_neumann(
        eps,
        &A);

    if(err)
        goto cleanup;


    err = lis_vector_duplicate(
        A,
        &b);

    if(err)
        goto cleanup;


    err = lis_vector_duplicate(
        A,
        &x);

    if(err)
        goto cleanup;


    err = lis_vector_duplicate(
        A,
        &work);

    if(err)
        goto cleanup;


    fill_rigid(
        b);

    err = lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);

    if(err)
        goto cleanup;


    err = lis_solver_create(
        &solver);

    if(err)
        goto cleanup;


    snprintf(
        option,
        sizeof(option),
        "-i %s "
        "-p none "
        "-scale none "
        "-print none "
        "-tol 1.0e-10 "
        "-maxiter %d "
        "-maxiter_noimp 0 "
        "-restart %d",
        solver_name,
        (int)maxiter,
        (int)restart);


    err = lis_solver_set_option(
        option,
        solver);

    if(err)
        goto cleanup;


    result->api =
        lis_solve(
            A,
            b,
            x,
            solver);


    err = lis_solver_get_status(
        solver,
        &result->status);

    if(err)
        goto cleanup;


    err = lis_solver_get_iter(
        solver,
        &result->iter);

    if(err)
        goto cleanup;


    err = lis_solver_get_residualnorm(
        solver,
        &result->reported);

    if(err)
        goto cleanup;


    result->true_res =
        true_relative_residual(
            A,
            b,
            x,
            work);


    err = lis_vector_nrm2(
        x,
        &result->xnorm);

    if(err)
        goto cleanup;


    err = LIS_SUCCESS;


cleanup:

    if(solver)
        lis_solver_destroy(solver);

    if(work)
        lis_vector_destroy(work);

    if(x)
        lis_vector_destroy(x);

    if(b)
        lis_vector_destroy(b);

    if(A)
        lis_matrix_destroy(A);


    return err;
}


/* ---------------------------------------------------------- */
/* 1x1 identity lucky-breakdown case                          */
/* ---------------------------------------------------------- */

static LIS_INT
run_lucky(
    const char *solver_name,
    RESULT *result)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR work = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;

    char option[256];


    memset(
        result,
        0,
        sizeof(*result));

    result->status = -999;
    result->iter = -999;
    result->reported = -1.0;
    result->true_res = -1.0;
    result->xnorm = -1.0;


    err = lis_matrix_create(
        LIS_COMM_WORLD,
        &A);

    if(err)
        goto cleanup;


    err = lis_matrix_set_size(
        A,
        0,
        1);

    if(err)
        goto cleanup;


    err = lis_matrix_set_value(
        LIS_INS_VALUE,
        0,
        0,
        (LIS_SCALAR)1.0,
        A);

    if(err)
        goto cleanup;


    err = lis_matrix_assemble(
        A);

    if(err)
        goto cleanup;


    err = lis_vector_duplicate(
        A,
        &b);

    if(err)
        goto cleanup;


    err = lis_vector_duplicate(
        A,
        &x);

    if(err)
        goto cleanup;


    err = lis_vector_duplicate(
        A,
        &work);

    if(err)
        goto cleanup;


    lis_vector_set_all(
        (LIS_SCALAR)1.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    err = lis_solver_create(
        &solver);

    if(err)
        goto cleanup;


    snprintf(
        option,
        sizeof(option),
        "-i %s "
        "-p none "
        "-scale none "
        "-print none "
        "-tol 1.0e-12 "
        "-maxiter 10 "
        "-maxiter_noimp 0 "
        "-restart 5",
        solver_name);


    err = lis_solver_set_option(
        option,
        solver);

    if(err)
        goto cleanup;


    result->api =
        lis_solve(
            A,
            b,
            x,
            solver);


    lis_solver_get_status(
        solver,
        &result->status);

    lis_solver_get_iter(
        solver,
        &result->iter);

    lis_solver_get_residualnorm(
        solver,
        &result->reported);


    result->true_res =
        true_relative_residual(
            A,
            b,
            x,
            work);


    lis_vector_nrm2(
        x,
        &result->xnorm);


    err = LIS_SUCCESS;


cleanup:

    if(solver)
        lis_solver_destroy(solver);

    if(work)
        lis_vector_destroy(work);

    if(x)
        lis_vector_destroy(x);

    if(b)
        lis_vector_destroy(b);

    if(A)
        lis_matrix_destroy(A);


    return err;
}


/* ---------------------------------------------------------- */

static void
print_result(
    const char *name,
    const RESULT *r)
{
    printf(
        "CASE %-24s "
        "api=%d "
        "status=%d "
        "iter=%d "
        "reported=% .12e "
        "true=% .12e "
        "xnorm=% .12e "
        "finite=%d\n",
        name,
        (int)r->api,
        (int)r->status,
        (int)r->iter,
        (double)r->reported,
        (double)r->true_res,
        (double)r->xnorm,
        result_is_finite(r));
}


/* ========================================================== */

int
main(
    int argc,
    char **argv)
{
    RESULT r;

    LIS_INT err;

    int failed = 0;

    LIS_REAL gap;


    err = lis_initialize(
        &argc,
        &argv);

    if(err)
        return 1;


    /*
     * Exact singular and incompatible system.
     *
     * These Krylov methods must stop finitely rather than
     * propagating NaN through the iterate.
     */

    {
        const char *names[] =
        {
            "minres",
            "bicgstab",
            "gmres",
            "fgmres"
        };

        int i;


        for(i=0; i<4; i++)
        {
            err = run_neumann(
                names[i],
                0.0,
                500,
                20,
                &r);

            print_result(
                names[i],
                &r);


            if(
                err!=LIS_SUCCESS
                ||
                r.api!=LIS_SUCCESS
                ||
                r.status!=LIS_BREAKDOWN
                ||
                !result_is_finite(&r))
            {
                fprintf(
                    stderr,
                    "FAIL exact singular %s\n",
                    names[i]);

                failed = 1;
            }
        }
    }


    /*
     * Near-singular rigid direction.
     *
     * CG and MINRES previously reported successful recursive
     * convergence even though the independent residual was far
     * above the requested tolerance.
     */

    err = run_neumann(
        "cg",
        1.0e-14,
        500,
        20,
        &r);

    print_result(
        "cg near-singular",
        &r);


    if(
        err!=LIS_SUCCESS
        ||
        r.api!=LIS_SUCCESS
        ||
        r.status==LIS_SUCCESS
        ||
        !result_is_finite(&r)
        ||
        r.true_res<=TEST_TOL)
    {
        fprintf(
            stderr,
            "FAIL CG false-convergence guard\n");

        failed = 1;
    }


    err = run_neumann(
        "minres",
        1.0e-14,
        500,
        20,
        &r);

    print_result(
        "minres near-singular",
        &r);


    if(
        err!=LIS_SUCCESS
        ||
        r.api!=LIS_SUCCESS
        ||
        r.status==LIS_SUCCESS
        ||
        !result_is_finite(&r)
        ||
        r.true_res<=TEST_TOL)
    {
        fprintf(
            stderr,
            "FAIL MINRES false-convergence guard\n");

        failed = 1;
    }


    /*
     * GMRES restart residual replacement.
     *
     * maxiter=500 is an exact multiple of restart=20, so the
     * final recursive state has just crossed a restart boundary.
     * The reported and independently computed residuals must
     * therefore remain on the same scale.
     */

    err = run_neumann(
        "gmres",
        1.0e-14,
        500,
        20,
        &r);

    print_result(
        "gmres restart residual",
        &r);


    if(
        err!=LIS_SUCCESS
        ||
        r.api!=LIS_SUCCESS
        ||
        r.status!=LIS_MAXITER
        ||
        !result_is_finite(&r)
        ||
        r.reported<=0.0
        ||
        r.true_res<=0.0)
    {
        fprintf(
            stderr,
            "FAIL GMRES restart residual state\n");

        failed = 1;
    }
    else
    {
        gap =
            r.reported > r.true_res
            ?
            r.reported / r.true_res
            :
            r.true_res / r.reported;


        printf(
            "GMRES_RESTART_GAP %.12e\n",
            (double)gap);


        if(gap>2.0)
        {
            fprintf(
                stderr,
                "FAIL GMRES restart residual gap %.12e\n",
                (double)gap);

            failed = 1;
        }
    }


    /*
     * A zero Arnoldi norm can be a valid lucky breakdown.
     * The 1x1 identity problem must still converge exactly.
     */

    err = run_lucky(
        "gmres",
        &r);

    print_result(
        "gmres lucky",
        &r);


    if(
        err!=LIS_SUCCESS
        ||
        r.api!=LIS_SUCCESS
        ||
        r.status!=LIS_SUCCESS
        ||
        !result_is_finite(&r)
        ||
        r.true_res>1.0e-12)
    {
        fprintf(
            stderr,
            "FAIL GMRES lucky breakdown\n");

        failed = 1;
    }


    err = run_lucky(
        "fgmres",
        &r);

    print_result(
        "fgmres lucky",
        &r);


    if(
        err!=LIS_SUCCESS
        ||
        r.api!=LIS_SUCCESS
        ||
        r.status!=LIS_SUCCESS
        ||
        !result_is_finite(&r)
        ||
        r.true_res>1.0e-12)
    {
        fprintf(
            stderr,
            "FAIL FGMRES lucky breakdown\n");

        failed = 1;
    }


    printf(
        "LIS_KRYLOV_ROBUSTNESS_STAGE6F %s\n",
        failed ? "FAILED" : "PASSED");


    lis_finalize();


    return failed ? 1 : 0;
}
