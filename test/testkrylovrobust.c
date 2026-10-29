#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#else
#ifdef HAVE_CONFIG_WIN_H
#include "lis_config_win.h"
#endif
#endif

#include <math.h>
#include "lis.h"

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

/*
 * Verify that an FGMRES projected convergence candidate is
 * certified against the explicit residual of the updated iterate.
 *
 * The projected residual can fall below tolerance before the
 * explicit residual does.  FGMRES must reject that candidate,
 * restart from the explicit residual, and continue the solve.
 */
static LIS_INT
run_reliable_case(
    LIS_INT n,
    const LIS_SCALAR *values,
    const LIS_SCALAR *rhs,
    char *options,
    RESULT *result)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR work = NULL;

    LIS_SOLVER solver = NULL;

    LIS_INT err = LIS_SUCCESS;
    LIS_INT i;
    LIS_INT j;


    memset(
        result,
        0,
        sizeof(*result));

    result->api = -999;
    result->status = -999;
    result->iter = -999;
    result->reported = -1.0;
    result->true_res = -1.0;
    result->xnorm = -1.0;


    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        goto cleanup;


    err =
        lis_matrix_set_size(
            A,
            0,
            n);

    if( err )
        goto cleanup;


    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            if(
                values[i*n+j]
                !=(LIS_SCALAR)0.0 )
            {
                err =
                    lis_matrix_set_value(
                        LIS_INS_VALUE,
                        i,
                        j,
                        values[i*n+j],
                        A);

                if( err )
                    goto cleanup;
            }
        }
    }


    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto cleanup;


    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto cleanup;


    err =
        lis_vector_duplicate(
            A,
            &b);

    if( err )
        goto cleanup;


    err =
        lis_vector_duplicate(
            A,
            &x);

    if( err )
        goto cleanup;


    err =
        lis_vector_duplicate(
            A,
            &work);

    if( err )
        goto cleanup;


    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        b);

    lis_vector_set_all(
        (LIS_SCALAR)0.0,
        x);


    for(i=0;i<n;i++)
    {
        b->value[i] =
            rhs[i];
    }


    err =
        lis_solver_create(
            &solver);

    if( err )
        goto cleanup;


    err =
        lis_solver_set_option(
            options,
            solver);

    if( err )
        goto cleanup;


    result->api =
        lis_solve(
            A,
            b,
            x,
            solver);


    err =
        lis_solver_get_status(
            solver,
            &result->status);

    if( err )
        goto cleanup;


    err =
        lis_solver_get_iter(
            solver,
            &result->iter);

    if( err )
        goto cleanup;


    err =
        lis_solver_get_residualnorm(
            solver,
            &result->reported);

    if( err )
        goto cleanup;


    result->true_res =
        true_relative_residual(
            A,
            b,
            x,
            work);


    err =
        lis_vector_nrm2(
            x,
            &result->xnorm);


cleanup:

    if( solver )
        lis_solver_destroy(
            solver);

    if( work )
        lis_vector_destroy(
            work);

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


/* ---------------------------------------------------------- */


static LIS_INT
run_fgmres_reliable_candidate(
    RESULT *result)
{
    static const LIS_SCALAR values[16] =
    {
         (LIS_SCALAR)9.6635805078063965e+11,
         (LIS_SCALAR)0.0,
         (LIS_SCALAR)1.0655399464330062e-07,
        -(LIS_SCALAR)2.9511571135475290e-12,

         (LIS_SCALAR)0.0,
         (LIS_SCALAR)8.0844035225373991e+03,
         (LIS_SCALAR)0.0,
         (LIS_SCALAR)0.0,

         (LIS_SCALAR)0.0,
         (LIS_SCALAR)0.0,
         (LIS_SCALAR)2.0245180553697102e-05,
         (LIS_SCALAR)0.0,

         (LIS_SCALAR)0.0,
         (LIS_SCALAR)1.1098830243238601e-03,
        -(LIS_SCALAR)1.9650356386922263e+02,
         (LIS_SCALAR)1.3661903277464329e+10
    };

    static const LIS_SCALAR rhs[4] =
    {
         (LIS_SCALAR)1.0,
        -(LIS_SCALAR)1.0,
        -(LIS_SCALAR)1.0,
        -(LIS_SCALAR)1.0
    };


    return
        run_reliable_case(
            4,
            values,
            rhs,
            "-i fgmres "
            "-p jacobi "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 20 "
            "-maxiter_noimp 0 "
            "-tol 1e-10",
            result);
}


/* ---------------------------------------------------------- */


/*
 * Verify that restarted right-preconditioned GMRES uses the
 * physical residual b-A*x at a restart boundary.
 *
 * The first restart leaves a residual component that Jacobi
 * preconditioning makes artificially small.  That must not be
 * accepted as convergence; GMRES must continue with the physical
 * residual and converge in a later cycle.
 */
static LIS_INT
check_gmres_reliable_restart(void)
{
    static const LIS_SCALAR values[4] =
    {
         (LIS_SCALAR)1.0e12,
        -(LIS_SCALAR)1.0,
         (LIS_SCALAR)0.0,
         (LIS_SCALAR)1.0
    };

    static const LIS_SCALAR rhs[2] =
    {
        (LIS_SCALAR)1.0,
        (LIS_SCALAR)1.0
    };

    RESULT r;

    LIS_INT err;
    LIS_INT failed = 0;


    err =
        run_reliable_case(
            2,
            values,
            rhs,
            "-i gmres "
            "-p jacobi "
            "-scale none "
            "-print none "
            "-restart 1 "
            "-maxiter 4 "
            "-maxiter_noimp 0 "
            "-tol 1e-10",
            &r);


    print_result(
        "gmres reliable restart",
        &r);


    if(
        err!=LIS_SUCCESS ||
        r.api!=LIS_SUCCESS ||
        r.status!=LIS_SUCCESS ||
        r.iter<=1 ||
        r.iter>4 ||
        !result_is_finite(&r) ||
        r.reported>1.0e-10 ||
        r.true_res>1.0e-10 )
    {
        fprintf(
            stderr,
            "GMRES RELIABLE RESTART FAILED\n");

        failed = 1;
    }
    else
    {
        printf(
            "GMRES RELIABLE RESTART PASSED\n");
    }


    return
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;
}


#endif


#if defined(USE_QUAD_PRECISION) && !defined(_COMPLEX)

/*
 * Run a small deterministic quad-precision GMRES-family problem
 * and return both the solver-reported and explicit residuals.
 */
static LIS_INT
run_quad_reliable_case(
    const LIS_SCALAR *values,
    const LIS_SCALAR *rhs,
    char *options,
    RESULT *result)
{
    return
        run_reliable_case(
            4,
            values,
            rhs,
            options,
            result);
}


/* ---------------------------------------------------------- */


/*
 * Quad GMRES restart regression.
 *
 * This matrix requires multiple restart cycles.  Reconstructing
 * the restart state algebraically can prevent convergence, while
 * restarting from the explicit physical residual b-A*x converges
 * to a certified solution.
 */
static LIS_INT
check_quad_gmres_reliable_restart(void)
{
    static const LIS_SCALAR values[16] =
    {
        (LIS_SCALAR) 1.22881149240986682e-07,
        (LIS_SCALAR) 0.00000000000000000e+00,
        (LIS_SCALAR) 1.54690992494933987e+10,
        (LIS_SCALAR) 1.52936057689083088e+06,

        (LIS_SCALAR) 0.00000000000000000e+00,
        (LIS_SCALAR)-1.51694359489564889e-01,
        (LIS_SCALAR) 1.92218355811301611e-11,
        (LIS_SCALAR) 1.63647935759591744e-14,

        (LIS_SCALAR) 0.00000000000000000e+00,
        (LIS_SCALAR) 0.00000000000000000e+00,
        (LIS_SCALAR) 1.09105007710154179e-07,
        (LIS_SCALAR) 1.37073376594539923e+02,

        (LIS_SCALAR) 1.72972035260890791e-10,
        (LIS_SCALAR) 0.00000000000000000e+00,
        (LIS_SCALAR) 0.00000000000000000e+00,
        (LIS_SCALAR)-6.56948468476498570e-01
    };

    static const LIS_SCALAR rhs[4] =
    {
        (LIS_SCALAR)-1.22930663336894730e+00,
        (LIS_SCALAR)-1.47749296734188285e+00,
        (LIS_SCALAR)-1.66024699245515350e+00,
        (LIS_SCALAR)-1.22799687471917718e+00
    };

    RESULT r;

    LIS_INT err;
    LIS_INT failed = 0;


    err =
        run_quad_reliable_case(
            values,
            rhs,
            "-i gmres "
            "-p jacobi "
            "-scale none "
            "-print none "
            "-restart 3 "
            "-maxiter 24 "
            "-maxiter_noimp 0 "
            "-tol 1e-10 "
            "-f quad",
            &r);


    print_result(
        "quad gmres reliable restart",
        &r);


    if(
        err!=LIS_SUCCESS
        ||
        r.api!=LIS_SUCCESS
        ||
        r.status!=LIS_SUCCESS
        ||
        r.iter<=3
        ||
        r.iter>24
        ||
        !result_is_finite(&r)
        ||
        r.reported>1.0e-10
        ||
        r.true_res>1.0e-10)
    {
        fprintf(
            stderr,
            "QUAD GMRES RELIABLE RESTART FAILED\n");

        failed = 1;
    }
    else
    {
        printf(
            "QUAD GMRES RELIABLE RESTART PASSED\n");
    }


    return
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;
}


/* ---------------------------------------------------------- */


/*
 * Quad FGMRES convergence-candidate regression.
 *
 * The projected residual can cross the requested tolerance while
 * the explicit residual remains large.  The kernel must reject
 * that candidate and continue instead of returning a false
 * success that is later rejected by the generic solve guard.
 */
static LIS_INT
check_quad_fgmres_reliable_candidate(void)
{
    static const LIS_SCALAR values[16] =
    {
        (LIS_SCALAR)-5.65082465591543848e+00,
        (LIS_SCALAR) 3.13406698132732377e-06,
        (LIS_SCALAR)-1.11737967063303498e-02,
        (LIS_SCALAR) 1.78564749677652944e+04,

        (LIS_SCALAR)-5.65082465597194705e+00,
        (LIS_SCALAR) 3.13406698129598313e-06,
        (LIS_SCALAR)-1.11737967064420886e-02,
        (LIS_SCALAR) 1.78564749675867279e+04,

        (LIS_SCALAR) 1.32952970409194200e-02,
        (LIS_SCALAR) 0.00000000000000000e+00,
        (LIS_SCALAR)-3.04880885448264641e-08,
        (LIS_SCALAR) 2.81881937980230432e-07,

        (LIS_SCALAR) 0.00000000000000000e+00,
        (LIS_SCALAR)-8.55093146947321366e+05,
        (LIS_SCALAR) 1.13729646288437489e+07,
        (LIS_SCALAR)-4.88471822752117291e+00
    };

    static const LIS_SCALAR rhs[4] =
    {
        (LIS_SCALAR)1.06421285574715152e+00,
        (LIS_SCALAR)8.57579685677863690e-01,
        (LIS_SCALAR)1.97352136137249756e+00,
        (LIS_SCALAR)5.95225387848174248e-01
    };

    RESULT r;

    LIS_INT err;
    LIS_INT failed = 0;
    LIS_INT terminal_ok = 0;


    err =
        run_quad_reliable_case(
            values,
            rhs,
            "-i fgmres "
            "-p jacobi "
            "-scale none "
            "-print none "
            "-restart 4 "
            "-maxiter 24 "
            "-maxiter_noimp 0 "
            "-tol 1e-14 "
            "-f quad",
            &r);


    print_result(
        "quad fgmres reliable candidate",
        &r);


    /*
     * The important invariant is that the projected candidate at
     * iteration 8 must not terminate the solve.  Later certified
     * convergence is valid; otherwise reaching MAXITER is also a
     * valid outcome for this deliberately difficult system.
     */
    if( r.status==LIS_SUCCESS )
    {
        terminal_ok =
            r.true_res<=1.0e-14;
    }
    else if( r.status==LIS_MAXITER )
    {
        terminal_ok =
            r.true_res>1.0e-14;
    }


    if(
        err!=LIS_SUCCESS
        ||
        r.api!=LIS_SUCCESS
        ||
        r.iter<=8
        ||
        !result_is_finite(&r)
        ||
        !terminal_ok)
    {
        fprintf(
            stderr,
            "QUAD FGMRES RELIABLE CANDIDATE FAILED\n");

        failed = 1;
    }
    else
    {
        printf(
            "QUAD FGMRES RELIABLE CANDIDATE PASSED\n");
    }


    return
        failed
        ? LIS_ERR_ILL_ARG
        : LIS_SUCCESS;
}


#endif


/* ---------------------------------------------------------- */


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
        !result_is_finite(&r)
        ||
        (
            r.status==LIS_SUCCESS
            &&
            r.true_res>TEST_TOL)
        ||
        (
            r.status!=LIS_SUCCESS
            &&
            r.true_res<=TEST_TOL))
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
        !result_is_finite(&r)
        ||
        (
            r.status==LIS_SUCCESS
            &&
            r.true_res>TEST_TOL)
        ||
        (
            r.status!=LIS_SUCCESS
            &&
            r.true_res<=TEST_TOL))
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
        !result_is_finite(&r)
        ||
        (
            r.status==LIS_SUCCESS
            &&
            (
                r.reported>TEST_TOL
                ||
                r.true_res>TEST_TOL))
        ||
        (
            r.status==LIS_MAXITER
            &&
            (
                r.reported<=0.0
                ||
                r.true_res<=0.0))
        ||
        (
            r.status!=LIS_SUCCESS
            &&
            r.status!=LIS_MAXITER))
    {
        fprintf(
            stderr,
            "FAIL GMRES restart residual state\n");

        failed = 1;
    }
    else if( r.status==LIS_MAXITER )
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
    else
    {
        printf(
            "GMRES restart case converged with certified residual\n");
    }



#ifndef _COMPLEX
    if(
        check_gmres_reliable_restart()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }
#endif


#ifndef _COMPLEX
    /*
     * A projected FGMRES convergence candidate must be certified
     * against the explicit residual of the updated iterate.
     *
     * Correctness is determined by the explicit residual of the
     * returned iterate rather than by a fixed iteration count.  The
     * latter can vary with floating-point reduction order, including
     * across OpenMP thread counts.
     */
    err =
        run_fgmres_reliable_candidate(
            &r);

    print_result(
        "fgmres reliable candidate",
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
        r.reported>1.0e-10
        ||
        r.true_res>1.0e-10)
    {
        fprintf(
            stderr,
            "FGMRES RELIABLE CANDIDATE FAILED\n");

        failed = 1;
    }
    else
    {
        printf(
            "FGMRES RELIABLE CANDIDATE PASSED\n");
    }
#endif


#if defined(USE_QUAD_PRECISION) && !defined(_COMPLEX)

    if(
        check_quad_gmres_reliable_restart()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }


    if(
        check_quad_fgmres_reliable_candidate()
        !=LIS_SUCCESS )
    {
        failed = 1;
    }

#endif


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
