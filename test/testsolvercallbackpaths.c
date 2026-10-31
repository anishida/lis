#include <stdio.h>
#include <math.h>
#include <string.h>

#include "lis_config.h"
#include "lis.h"

typedef struct
{
    const char *name;

    LIS_INT monitor_projected;
    LIS_INT monitor_true;

    LIS_INT convergence_projected;
    LIS_INT convergence_true;

    LIS_INT accepted_projected;
    LIS_INT accepted_true;
} CALLBACK_STATE;


static void
monitor_callback(
    LIS_SOLVER solver,
    LIS_INT iter,
    LIS_REAL residual,
    LIS_REAL tolerance,
    LIS_INT residual_kind,
    void *user_data)
{
    CALLBACK_STATE *state =
        (CALLBACK_STATE *)user_data;

    (void)solver;
    (void)iter;
    (void)residual;
    (void)tolerance;

    if( residual_kind==LIS_SOLVER_RESIDUAL_PROJECTED )
    {
        state->monitor_projected++;
    }
    else if( residual_kind==LIS_SOLVER_RESIDUAL_TRUE )
    {
        state->monitor_true++;
    }
}


static LIS_INT
convergence_callback(
    LIS_SOLVER solver,
    LIS_INT iter,
    LIS_REAL residual,
    LIS_REAL tolerance,
    LIS_INT residual_kind,
    void *user_data)
{
    CALLBACK_STATE *state =
        (CALLBACK_STATE *)user_data;

    LIS_INT accept;

    (void)solver;
    (void)iter;

    accept =
        residual<=tolerance
        ? LIS_TRUE
        : LIS_FALSE;

    if( residual_kind==LIS_SOLVER_RESIDUAL_PROJECTED )
    {
        state->convergence_projected++;

        if( accept )
            state->accepted_projected++;
    }
    else if( residual_kind==LIS_SOLVER_RESIDUAL_TRUE )
    {
        state->convergence_true++;

        if( accept )
            state->accepted_true++;
    }

    return accept;
}


static LIS_INT
build_matrix(
    LIS_MATRIX *Aout)
{
    LIS_MATRIX A = NULL;

    LIS_INT err;
    LIS_INT i;

    const LIS_INT n = 96;

    err =
        lis_matrix_create(
            LIS_COMM_WORLD,
            &A);

    if( err )
        return err;

    err =
        lis_matrix_set_size(
            A,
            0,
            n);

    if( err )
        goto fail;

    for(i=0;i<n;i++)
    {
        err =
            lis_matrix_set_value(
                LIS_INS_VALUE,
                i,
                i,
                (LIS_SCALAR)4.0,
                A);

        if( err )
            goto fail;

        if( i>0 )
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i-1,
                    (LIS_SCALAR)-1.0,
                    A);

            if( err )
                goto fail;
        }

        if( i<n-1 )
        {
            err =
                lis_matrix_set_value(
                    LIS_INS_VALUE,
                    i,
                    i+1,
                    (LIS_SCALAR)-0.65,
                    A);

            if( err )
                goto fail;
        }
    }

    err =
        lis_matrix_set_type(
            A,
            LIS_MATRIX_CSR);

    if( err )
        goto fail;

    err =
        lis_matrix_assemble(
            A);

    if( err )
        goto fail;

    *Aout = A;

    return LIS_SUCCESS;

fail:

    if( A )
        lis_matrix_destroy(A);

    return err;
}


static void
set_exact_solution(
    LIS_VECTOR u)
{
    LIS_INT i;

    for(i=0;i<u->n;i++)
    {
        LIS_REAL g =
            (LIS_REAL)(u->is+i+1);

        u->value[i] =
            (LIS_SCALAR)(
                1.0
                +0.07*sin(0.13*g)
                +0.03*cos(0.07*g));
    }
}


static LIS_REAL
true_relative_residual(
    LIS_MATRIX A,
    LIS_VECTOR b,
    LIS_VECTOR x,
    LIS_VECTOR r)
{
    LIS_REAL rnrm;
    LIS_REAL bnrm;

    lis_matvec(
        A,
        x,
        r);

    lis_vector_xpay(
        b,
        -1.0,
        r);

    lis_vector_nrm2(
        r,
        &rnrm);

    lis_vector_nrm2(
        b,
        &bnrm);

    if( bnrm>0.0 )
        return rnrm/bnrm;

    return rnrm;
}


static int
run_case(
    const char *name,
    const char *options)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR u = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    CALLBACK_STATE state;

    LIS_INT err;
    LIS_INT status;

    LIS_REAL tr;

    int failed = 0;

    memset(
        &state,
        0,
        sizeof(state));

    state.name = name;

    err =
        build_matrix(
            &A);

    if( err )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&u)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&b)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&x)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&r)) )
        goto api_fail;

    set_exact_solution(
        u);

    lis_matvec(
        A,
        u,
        b);

    lis_vector_set_all(
        0.0,
        x);

    err =
        lis_solver_create(
            &solver);

    if( err )
        goto api_fail;

    err =
        lis_solver_set_option(
            (char *)options,
            solver);

    if( err )
        goto api_fail;

    err =
        lis_solver_set_monitor(
            solver,
            monitor_callback,
            &state);

    if( err )
        goto api_fail;

    err =
        lis_solver_set_convergence(
            solver,
            convergence_callback,
            &state);

    if( err )
        goto api_fail;

    err =
        lis_solve(
            A,
            b,
            x,
            solver);

    lis_solver_get_status(
        solver,
        &status);

    tr =
        true_relative_residual(
            A,
            b,
            x,
            r);

    printf(
        "%-8s api=%d status=%d iter=%d "
        "mp=%d mt=%d cp=%d ct=%d "
        "ap=%d at=%d true=%.3e\n",
        name,
        (int)err,
        (int)status,
        (int)solver->iter,
        (int)state.monitor_projected,
        (int)state.monitor_true,
        (int)state.convergence_projected,
        (int)state.convergence_true,
        (int)state.accepted_projected,
        (int)state.accepted_true,
        (double)tr);

    if( err!=LIS_SUCCESS ||
        status!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "%s did not converge successfully\n",
            name);

        failed = 1;
    }

    if( state.monitor_projected<=0 ||
        state.convergence_projected<=0 )
    {
        fprintf(
            stderr,
            "%s did not report PROJECTED residuals\n",
            name);

        failed = 1;
    }

    if( state.monitor_true<=0 ||
        state.convergence_true<=0 )
    {
        fprintf(
            stderr,
            "%s did not report TRUE residuals\n",
            name);

        failed = 1;
    }

    /*
     * A successful solve must include an accepted TRUE residual.
     * PROJECTED acceptance alone is never sufficient.
     */
    if( state.accepted_true<=0 )
    {
        fprintf(
            stderr,
            "%s succeeded without accepted TRUE residual\n",
            name);

        failed = 1;
    }

    if( tr>1.0e-7 )
    {
        fprintf(
            stderr,
            "%s final true residual too large: %.6e\n",
            name,
            (double)tr);

        failed = 1;
    }

    goto cleanup;

api_fail:

    fprintf(
        stderr,
        "%s API failure rc=%d\n",
        name,
        (int)err);

    failed = 1;

cleanup:

    if( solver )
        lis_solver_destroy(solver);

    if( r )
        lis_vector_destroy(r);

    if( x )
        lis_vector_destroy(x);

    if( b )
        lis_vector_destroy(b);

    if( u )
        lis_vector_destroy(u);

    if( A )
        lis_matrix_destroy(A);

    return failed;
}



#if defined(USE_QUAD_PRECISION) && defined(_COMPLEX)
/*
 * Complex double-double is not implemented.  Verify that requesting it
 * is rejected before arithmetic runs, rather than producing an invalid
 * result or spending maxiter iterations on a broken recurrence.
 */
static int
check_unsupported_complex_quad(void)
{
    const char *options[] = {
        "-i gmres -p none -f quad -maxiter 40 -print none",
        "-i gmres -p none -f switch -maxiter 40 -print none",
        "-i cg -p none -f quad -maxiter 40 -print none"
    };
    LIS_MATRIX A = NULL;
    LIS_VECTOR b = NULL, x = NULL;
    LIS_SOLVER solver = NULL;
    LIS_INT err, status, i;
    int failed = 0;

    err = build_matrix(&A);
    if(err) goto fail;
    err = lis_vector_duplicate(A,&b);
    if(err) goto fail;
    err = lis_vector_duplicate(A,&x);
    if(err) goto fail;
    lis_vector_set_all((LIS_SCALAR)1.0,b);

    for(i=0;i<3;i++)
    {
        err = lis_solver_create(&solver);
        if(err) goto fail;
        err = lis_solver_set_option((char *)options[i],solver);
        if(err) goto fail;
        lis_vector_set_all((LIS_SCALAR)0.0,x);
        err = lis_solve(A,b,x,solver);
        status = -1;
        lis_solver_get_status(solver,&status);
        if(err!=LIS_ERR_NOT_IMPLEMENTED ||
           status!=LIS_ERR_NOT_IMPLEMENTED)
        {
            fprintf(stderr,
                    "Complex double-double was not rejected (%s): api=%d status=%d\n",
                    options[i],(int)err,(int)status);
            failed = 1;
        }
        lis_solver_destroy(solver);
        solver = NULL;
    }
    goto cleanup;

fail:
    fprintf(stderr,"Complex precision rejection test setup failed: %d\n",
            (int)err);
    failed = 1;

cleanup:
    if(solver) lis_solver_destroy(solver);
    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    return failed;
}
#endif

int
main(
    int argc,
    char **argv)
{
    int failed = 0;

    LIS_INT err;

    err =
        lis_initialize(
            &argc,
            &argv);

    if( err )
    {
        fprintf(
            stderr,
            "lis_initialize failed rc=%d\n",
            (int)err);

        return 1;
    }

    failed |=
        run_case(
            "GMRES",
            "-i gmres "
            "-p jacobi "
            "-restart 8 "
            "-maxiter 2000 "
            "-tol 1e-8 "
            "-print none");

    failed |=
        run_case(
            "FGMRES",
            "-i fgmres "
            "-p jacobi "
            "-restart 8 "
            "-maxiter 2000 "
            "-tol 1e-8 "
            "-print none");


#ifdef USE_QUAD_PRECISION

#ifdef _COMPLEX
    failed |= check_unsupported_complex_quad();
#else
    failed |=
        run_case(
            "GMRES_SW",
            "-i gmres "
            "-p jacobi "
            "-f switch "
            "-restart 8 "
            "-switch_maxiter 1 "
            "-switch_tol 9e-1 "
            "-maxiter 2000 "
            "-tol 1e-8 "
            "-print none");
#endif /* !_COMPLEX */

#endif

#ifndef _COMPLEX

    failed |=
        run_case(
            "GCRODR",
            "-i gcrodr "
            "-p jacobi "
            "-restart 8 "
            "-maxiter 2000 "
            "-tol 1e-8 "
            "-gcrodr_keep 2 "
            "-print none");

    failed |=
        run_case(
            "FGCRODR",
            "-i fgcrodr "
            "-p jacobi "
            "-restart 8 "
            "-maxiter 2000 "
            "-tol 1e-8 "
            "-fgcrodr_k 2 "
            "-print none");

#endif

    lis_finalize();

    printf(
        "LIS_SOLVER_CALLBACK_PATHS_STAGE6 %s\n",
        failed ? "FAILED" : "PASSED");

    return failed ? 1 : 0;
}
