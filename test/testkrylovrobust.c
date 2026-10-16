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

#ifndef _COMPLEX

static LIS_INT
check_fgmres_noimp_progress(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err;
    LIS_INT api;

    LIS_REAL true_residual;
    LIS_REAL xnorm;

    LIS_INT failed;


    failed =
        0;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        goto api_fail;


    err =
        lis_matrix_set_size(
            A,
            0,
            4);

    if( err )
        goto api_fail;


#define SET_A(row,col,val)                                      \
    do                                                          \
    {                                                           \
        err =                                                   \
            lis_matrix_set_value(                               \
                LIS_INS_VALUE,                                  \
                (row),                                          \
                (col),                                          \
                (LIS_SCALAR)(val),                              \
                A);                                             \
                                                                \
        if( err )                                               \
            goto api_fail;                                      \
    }                                                           \
    while(0)


    /*
     * Deterministic no-improvement system:
     *
     *     A =
     *
     *     [ 1  1  1  0 ]
     *     [ 1  1  1  1 ]
     *     [ 0  1  1  0 ]
     *     [ 0  0  1  0 ]
     *
     *     b = e0
     *
     * The least-squares residual is unchanged after the first
     * three Arnoldi iterations:
     *
     *     ||r_ls|| = 1/sqrt(2).
     *
     * With maxiter_noimp=1, FGMRES must terminate with MAXITER
     * at iteration 3, but must preserve the current Krylov
     * approximation instead of returning the initial x.
     */

    SET_A(0,0,1.0);
    SET_A(1,0,1.0);

    SET_A(0,1,1.0);
    SET_A(1,1,1.0);
    SET_A(2,1,1.0);

    SET_A(0,2,1.0);
    SET_A(1,2,1.0);
    SET_A(2,2,1.0);
    SET_A(3,2,1.0);

    SET_A(1,3,1.0);


#undef SET_A


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto api_fail;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &b);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if( err )
        goto api_fail;


    err =
        lis_vector_duplicate(
            A,
            &r);

    if( err )
        goto api_fail;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    b->value[0] =
        (LIS_SCALAR)1.0;


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;


    err =
        lis_solver_set_option(
            "-i fgmres "
            "-p none "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 3 "
            "-maxiter_noimp 1 "
            "-tol 1e-14",
            solver);

    if( err )
        goto api_fail;


    api =
        lis_solve(
            A,
            b,
            x,
            solver);


    if( api!=LIS_SUCCESS )
    {
        err =
            api;

        goto api_fail;
    }


    err =
        lis_matvec(
            A,
            x,
            r);

    if( err )
        goto api_fail;


    err =
        lis_vector_xpay(
            b,
            (LIS_SCALAR)-1.0,
            r);

    if( err )
        goto api_fail;


    err =
        lis_vector_nrm2(
            r,
            &true_residual);

    if( err )
        goto api_fail;


    err =
        lis_vector_nrm2(
            x,
            &xnorm);

    if( err )
        goto api_fail;


    printf(
        "FGMRES noimp progress "
        "api=%d status=%d iter=%d "
        "reported=% .17e true=% .17e xnorm=% .17e "
        "x=[% .17e,% .17e,% .17e,% .17e]\n",
        (int)api,
        (int)solver->retcode,
        (int)solver->iter,
        (double)solver->resid,
        (double)true_residual,
        (double)xnorm,
        (double)x->value[0],
        (double)x->value[1],
        (double)x->value[2],
        (double)x->value[3]);


    if(
        solver->retcode!=LIS_MAXITER ||
        solver->iter!=3 ||
        fabs(
            solver->resid
            -0.7071067811865476)
            >1.0e-12 ||
        fabs(
            true_residual
            -0.7071067811865476)
            >1.0e-12 ||
        fabs(
            xnorm
            -0.5)
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[0]
            -0.5)
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[1])
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[2])
            >1.0e-12 ||
        fabs(
            (LIS_REAL)x->value[3])
            >1.0e-12 )
    {
        fprintf(
            stderr,
            "FGMRES NO-IMPROVEMENT PROGRESS FAILED\n");

        failed =
            1;
    }
    else
    {
        printf(
            "FGMRES NO-IMPROVEMENT PROGRESS PASSED\n");
    }


    err =
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;

    goto cleanup;


api_fail:

    fprintf(
        stderr,
        "FGMRES no-improvement API failure: %d\n",
        (int)err);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( r )
        lis_vector_destroy(
            r);

    if( x )
        lis_vector_destroy(
            x);

    if( b )
        lis_vector_destroy(
            b);

    if( A )
        lis_matrix_destroy(
            A);


    return err;
}

#endif


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


#ifndef _COMPLEX
    if(
        check_fgmres_noimp_progress()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }
#endif


    lis_finalize();


    return failed ? 1 : 0;
}
