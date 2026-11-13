#include <stdio.h>
#include <math.h>
#include "lislib.h"

static int marker_a = 11;
static int marker_b = 22;

static LIS_INT monitor_calls = 0;
static LIS_INT monitor_iter = -1;
static LIS_INT monitor_kind = -1;
static LIS_REAL monitor_residual = 0.0;
static LIS_REAL monitor_tolerance = 0.0;
static void *monitor_seen_data = NULL;

static LIS_INT convergence_calls = 0;
static LIS_INT convergence_iter = -1;
static LIS_INT convergence_kind = -1;
static LIS_REAL convergence_residual = 0.0;
static LIS_REAL convergence_tolerance = 0.0;
static void *convergence_seen_data = NULL;

static void
monitor_callback(
    LIS_SOLVER solver,
    LIS_INT iter,
    LIS_REAL residual,
    LIS_REAL tolerance,
    LIS_INT residual_kind,
    void *user_data)
{
    (void)solver;

    monitor_calls++;
    monitor_iter = iter;
    monitor_residual = residual;
    monitor_tolerance = tolerance;
    monitor_kind = residual_kind;
    monitor_seen_data = user_data;
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
    (void)solver;

    convergence_calls++;
    convergence_iter = iter;
    convergence_residual = residual;
    convergence_tolerance = tolerance;
    convergence_kind = residual_kind;
    convergence_seen_data = user_data;

    /*
     * Deliberately independent of default_converged.
     * This verifies that an installed callback owns the
     * convergence candidate decision.
     */
    return iter>=3 ? LIS_TRUE : LIS_FALSE;
}

static int
real_equal(LIS_REAL a, LIS_REAL b)
{
    return fabs(a-b) <= 1.0e-15;
}

int
main(int argc, char **argv)
{
    LIS_SOLVER solver = NULL;
    LIS_INT err;
    LIS_INT converged;
    int failed = 0;

    err = lis_initialize(&argc,&argv);

    if( err )
    {
        fprintf(
            stderr,
            "initialize failed rc=%d\n",
            (int)err);
        return 1;
    }

    if( LIS_SOLVER_RESIDUAL_RECURSIVE!=0 ||
        LIS_SOLVER_RESIDUAL_PROJECTED!=1 ||
        LIS_SOLVER_RESIDUAL_TRUE!=2 )
    {
        fprintf(
            stderr,
            "residual kind constants mismatch\n");
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_create(&solver);

    if( err || solver==NULL )
    {
        fprintf(
            stderr,
            "solver create failed rc=%d\n",
            (int)err);
        failed = 1;
        goto cleanup;
    }

    if( solver->monitor!=NULL ||
        solver->monitor_data!=NULL ||
        solver->convergence!=NULL ||
        solver->convergence_data!=NULL )
    {
        fprintf(
            stderr,
            "default callback state mismatch\n");
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_monitor(
        solver,
        monitor_callback,
        &marker_a);

    if( err )
    {
        fprintf(
            stderr,
            "monitor install failed rc=%d\n",
            (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_convergence(
        solver,
        convergence_callback,
        &marker_b);

    if( err )
    {
        fprintf(
            stderr,
            "convergence install failed rc=%d\n",
            (int)err);
        failed = 1;
        goto cleanup;
    }

    lis_solver_monitor_residual(
        solver,
        7,
        0.125,
        0.010,
        LIS_SOLVER_RESIDUAL_PROJECTED);

    if( monitor_calls!=1 ||
        monitor_iter!=7 ||
        monitor_kind!=LIS_SOLVER_RESIDUAL_PROJECTED ||
        !real_equal(monitor_residual,0.125) ||
        !real_equal(monitor_tolerance,0.010) ||
        monitor_seen_data!=&marker_a )
    {
        fprintf(
            stderr,
            "monitor helper forwarding mismatch\n");
        failed = 1;
        goto cleanup;
    }

    /*
     * Installed convergence callback overrides the caller's
     * default decision.
     */
    converged = lis_solver_test_convergence(
        solver,
        2,
        0.005,
        0.010,
        LIS_SOLVER_RESIDUAL_RECURSIVE,
        LIS_TRUE);

    if( converged!=LIS_FALSE )
    {
        fprintf(
            stderr,
            "custom convergence failed to reject default success\n");
        failed = 1;
        goto cleanup;
    }

    converged = lis_solver_test_convergence(
        solver,
        3,
        0.500,
        0.010,
        LIS_SOLVER_RESIDUAL_TRUE,
        LIS_FALSE);

    if( converged!=LIS_TRUE )
    {
        fprintf(
            stderr,
            "custom convergence failed to accept candidate\n");
        failed = 1;
        goto cleanup;
    }

    if( convergence_calls!=2 ||
        convergence_iter!=3 ||
        convergence_kind!=LIS_SOLVER_RESIDUAL_TRUE ||
        !real_equal(convergence_residual,0.500) ||
        !real_equal(convergence_tolerance,0.010) ||
        convergence_seen_data!=&marker_b )
    {
        fprintf(
            stderr,
            "convergence helper forwarding mismatch\n");
        failed = 1;
        goto cleanup;
    }

    /*
     * Sequence reset must preserve both callbacks.
     */
    err = lis_solver_reset_sequence(solver);

    if( err ||
        solver->monitor!=monitor_callback ||
        solver->monitor_data!=&marker_a ||
        solver->convergence!=convergence_callback ||
        solver->convergence_data!=&marker_b )
    {
        fprintf(
            stderr,
            "sequence reset did not preserve callbacks rc=%d\n",
            (int)err);
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_convergence(
        solver,
        NULL,
        &marker_b);

    if( err ||
        solver->convergence!=NULL ||
        solver->convergence_data!=NULL )
    {
        fprintf(
            stderr,
            "convergence clear failed rc=%d\n",
            (int)err);
        failed = 1;
        goto cleanup;
    }

    /*
     * Without a custom callback, preserve the caller's exact
     * existing convergence decision.
     */
    converged = lis_solver_test_convergence(
        solver,
        9,
        1.0,
        2.0,
        LIS_SOLVER_RESIDUAL_RECURSIVE,
        LIS_TRUE);

    if( converged!=LIS_TRUE )
    {
        fprintf(
            stderr,
            "default true convergence was not preserved\n");
        failed = 1;
        goto cleanup;
    }

    converged = lis_solver_test_convergence(
        solver,
        9,
        1.0,
        2.0,
        LIS_SOLVER_RESIDUAL_RECURSIVE,
        LIS_FALSE);

    if( converged!=LIS_FALSE )
    {
        fprintf(
            stderr,
            "default false convergence was not preserved\n");
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_monitor(
        solver,
        NULL,
        &marker_a);

    if( err ||
        solver->monitor!=NULL ||
        solver->monitor_data!=NULL )
    {
        fprintf(
            stderr,
            "monitor clear failed rc=%d\n",
            (int)err);
        failed = 1;
        goto cleanup;
    }

    /*
     * NULL monitor helper is intentionally a no-op.
     */
    lis_solver_monitor_residual(
        solver,
        10,
        1.0,
        0.1,
        LIS_SOLVER_RESIDUAL_TRUE);

    if( monitor_calls!=1 )
    {
        fprintf(
            stderr,
            "cleared monitor was unexpectedly invoked\n");
        failed = 1;
        goto cleanup;
    }

    err = lis_solver_set_convergence(
        NULL,
        convergence_callback,
        &marker_b);

    if( err!=LIS_ERR_ILL_ARG )
    {
        fprintf(
            stderr,
            "NULL convergence solver rc=%d expected=%d\n",
            (int)err,
            (int)LIS_ERR_ILL_ARG);
        failed = 1;
        goto cleanup;
    }

    printf(
        "CALLBACK_HELPERS "
        "monitor=pass "
        "custom_override=pass "
        "default_preserve=pass "
        "reset=pass "
        "clear=pass\n");

cleanup:

    if( solver )
    {
        lis_solver_destroy(solver);
    }

    lis_finalize();

    printf(
        "LIS_SOLVER_CALLBACK_STAGE5 %s\n",
        failed ? "FAILED" : "PASSED");

    return failed ? 1 : 0;
}
