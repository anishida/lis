#include <stdio.h>
#include <math.h>
#include <string.h>

#include "lis.h"

#define POLICY_ACCEPT_ALL  1
#define POLICY_REJECT_ALL  2


typedef struct
{
    LIS_INT mode;

    LIS_INT projected_calls;
    LIS_INT true_calls;
} POLICY_STATE;


static LIS_INT
policy_callback(
    LIS_SOLVER solver,
    LIS_INT iter,
    LIS_REAL residual,
    LIS_REAL tolerance,
    LIS_INT residual_kind,
    void *user_data)
{
    POLICY_STATE *state =
        (POLICY_STATE *)user_data;

    (void)solver;
    (void)iter;

    if( residual_kind==LIS_SOLVER_RESIDUAL_PROJECTED )
    {
        state->projected_calls++;
    }
    else if( residual_kind==LIS_SOLVER_RESIDUAL_TRUE )
    {
        state->true_calls++;
    }

    if( state->mode==POLICY_ACCEPT_ALL )
        return LIS_TRUE;

    if( state->mode==POLICY_REJECT_ALL )
        return LIS_FALSE;

    return residual<=tolerance
        ? LIS_TRUE
        : LIS_FALSE;
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
run_final_hard_certification(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR u = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    POLICY_STATE state;

    LIS_INT err = LIS_SUCCESS;
    LIS_INT status = LIS_SUCCESS;

    LIS_REAL tr = 0.0;

    int failed = 0;

    memset(
        &state,
        0,
        sizeof(state));

    state.mode =
        POLICY_ACCEPT_ALL;

    if( (err=build_matrix(&A)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&u)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&b)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&x)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&r)) )
        goto api_fail;

    set_exact_solution(u);

    lis_matvec(
        A,
        u,
        b);

    lis_vector_set_all(
        0.0,
        x);

    if( (err=lis_solver_create(&solver)) )
        goto api_fail;

    err =
        lis_solver_set_option(
            "-i gmres "
            "-p jacobi "
            "-restart 8 "
            "-maxiter 20 "
            "-tol 1e-14 "
            "-print none",
            solver);

    if( err )
        goto api_fail;

    err =
        lis_solver_set_convergence(
            solver,
            policy_callback,
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
        "FINAL_CERT "
        "api=%d status=%d "
        "projected=%d true=%d "
        "residual=%.3e\n",
        (int)err,
        (int)status,
        (int)state.projected_calls,
        (int)state.true_calls,
        (double)tr);

    /*
     * The callback accepts the first convergence candidate.
     * The outer physical certification must still reject
     * a final true residual above the solver tolerance.
     */
    if( err!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "final certification API failure rc=%d\n",
            (int)err);

        failed = 1;
    }

    if( status!=LIS_BREAKDOWN )
    {
        fprintf(
            stderr,
            "hard final certification failed status=%d\n",
            (int)status);

        failed = 1;
    }

    if( state.projected_calls<=0 ||
        state.true_calls<=0 )
    {
        fprintf(
            stderr,
            "accept-all callback did not exercise "
            "PROJECTED and TRUE paths\n");

        failed = 1;
    }

    if( tr<=1.0e-14 )
    {
        fprintf(
            stderr,
            "hard-cert residual unexpectedly small: %.6e\n",
            (double)tr);

        failed = 1;
    }

    goto cleanup;

api_fail:

    fprintf(
        stderr,
        "final certification setup failed rc=%d\n",
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


static int
run_sequence_veto(void)
{
    LIS_MATRIX A = NULL;

    LIS_VECTOR u = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR r = NULL;

    LIS_SOLVER solver = NULL;

    POLICY_STATE state;

    LIS_INT err = LIS_SUCCESS;
    LIS_INT status = LIS_SUCCESS;

    LIS_REAL tr = 0.0;

    int failed = 0;

    memset(
        &state,
        0,
        sizeof(state));

    state.mode =
        POLICY_REJECT_ALL;

    if( (err=build_matrix(&A)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&u)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&b)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&x)) )
        goto api_fail;

    if( (err=lis_vector_duplicate(A,&r)) )
        goto api_fail;

    set_exact_solution(u);

    lis_matvec(
        A,
        u,
        b);

    lis_vector_set_all(
        0.0,
        x);

    if( (err=lis_solver_create(&solver)) )
        goto api_fail;

    err =
        lis_solver_set_option(
            "-i gmres "
            "-p jacobi "
            "-restart 8 "
            "-maxiter 1 "
            "-tol 9.9e-1 "
            "-print none",
            solver);

    if( err )
        goto api_fail;

    err =
        lis_solver_set_convergence(
            solver,
            policy_callback,
            &state);

    if( err )
        goto api_fail;

    err =
        lis_solve_sequence(
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
        "SEQUENCE_VETO "
        "api=%d status=%d "
        "projected=%d true=%d "
        "residual=%.3e\n",
        (int)err,
        (int)status,
        (int)state.projected_calls,
        (int)state.true_calls,
        (double)tr);

    if( err!=LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "sequence API failure rc=%d\n",
            (int)err);

        failed = 1;
    }

    if( tr>9.9e-1 )
    {
        fprintf(
            stderr,
            "sequence test did not reach physical tolerance: %.6e\n",
            (double)tr);

        failed = 1;
    }

    if( status==LIS_SUCCESS )
    {
        fprintf(
            stderr,
            "sequence promotion ignored callback veto\n");

        failed = 1;
    }

    if( state.projected_calls<=0 ||
        state.true_calls<=0 )
    {
        fprintf(
            stderr,
            "sequence callback paths were not exercised\n");

        failed = 1;
    }

    goto cleanup;

api_fail:

    fprintf(
        stderr,
        "sequence veto setup failed rc=%d\n",
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


int
main(
    int argc,
    char **argv)
{
    LIS_INT err;

    int failed = 0;

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
        run_final_hard_certification();

    failed |=
        run_sequence_veto();

    lis_finalize();

    printf(
        "LIS_SOLVER_CALLBACK_CERT_STAGE7 %s\n",
        failed ? "FAILED" : "PASSED");

    return failed ? 1 : 0;
}
