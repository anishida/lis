#include "lis.h"
#include "lis_solver_options.h"
#include <stdio.h>
#include <math.h>

static void test_monitor(
    struct LIS_SOLVER_STRUCT *solver,
    LIS_INT iter,
    LIS_REAL residual,
    LIS_REAL tolerance,
    LIS_INT residual_kind,
    void *user_data)
{
    (void)solver; (void)iter; (void)residual; (void)tolerance;
    (void)residual_kind; (void)user_data;
}

static LIS_INT test_convergence(
    struct LIS_SOLVER_STRUCT *solver,
    LIS_INT iter,
    LIS_REAL residual,
    LIS_REAL tolerance,
    LIS_INT residual_kind,
    void *user_data)
{
    (void)solver; (void)iter; (void)residual; (void)tolerance;
    (void)residual_kind; (void)user_data;
    return LIS_FALSE;
}

static int build_diagonal_problem(
    LIS_MATRIX *A,
    LIS_VECTOR *b,
    LIS_VECTOR *x)
{
    LIS_INT err,i;

    *A = NULL;
    *b = NULL;
    *x = NULL;

    err = lis_matrix_create(LIS_COMM_WORLD,A);
    if(err) return 1;
    err = lis_matrix_set_size(*A,0,4);
    if(err) return 1;

    for(i=0;i<4;i++)
    {
        err = lis_matrix_set_value(
            LIS_INS_VALUE,i,i,(LIS_SCALAR)(2.0+i),*A);
        if(err) return 1;
    }

    err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
    if(err) return 1;
    err = lis_matrix_assemble(*A);
    if(err) return 1;

    err = lis_vector_duplicate(*A,b);
    if(err) return 1;
    err = lis_vector_duplicate(*A,x);
    if(err) return 1;
    err = lis_vector_set_all((LIS_SCALAR)1.0,*b);
    if(err) return 1;
    err = lis_vector_set_all((LIS_SCALAR)0.0,*x);
    if(err) return 1;

    return 0;
}

int main(int argc, char **argv)
{
    LIS_SOLVER solver = NULL;
    LIS_MATRIX A = NULL;
    LIS_VECTOR b = NULL;
    LIS_VECTOR x = NULL;
    LIS_VECTOR q = NULL;
    LIS_VECTOR basis[1];
    LIS_PRECON saved_precon = NULL;
    LIS_SOLVE_REPORT report;
    LIS_REAL tolerance,residual;
    LIS_INT err,method,precon,policy,age,threads,half_mode,status,changed;
    LIS_INT maxiter,iter,iter_double,iter_quad;
    double time,itime,ptime,p_c_time,p_i_time;
    int monitor_cookie = 11;
    int convergence_cookie = 22;
    int failed = 0;

    err = lis_initialize(&argc,&argv);
    if(err) return 1;

    err = lis_solver_create(&solver);
    if(err) { lis_finalize(); return 1; }

    err = lis_solver_set_sequence_policy(solver,LIS_SEQUENCE_POLICY_AUTO);
    if(err) failed = 1;
    err = lis_solver_set_sequence_precon_max_age(solver,3);
    if(err) failed = 1;
    err = lis_solver_set_num_threads(solver,1);
    if(err) failed = 1;
    err = lis_solver_set_precon_half(solver,LIS_PRECON_HALF_ON);
    if(err) failed = 1;

    if(!failed)
        err = lis_solver_set_option_replace(
            "-i fgmres -p ilu -maxiter 77 -restart 12 -tol 1e-7",solver);
    if(err) failed = 1;

    lis_solver_get_solver(solver,&method);
    lis_solver_get_precon(solver,&precon);
    if(method!=LIS_SOLVER_FGMRES || precon!=LIS_PRECON_TYPE_ILU ||
       solver->options[LIS_OPTIONS_MAXITER]!=77 ||
       solver->options[LIS_OPTIONS_RESTART]!=12)
    {
        fprintf(stderr,"first replacement mismatch\n");
        failed = 1;
    }

    /*
     * Build real persistent sequence/preconditioner state before installing
     * application-owned callback/nullspace objects.  The basis used below is
     * intentionally only a preservation marker; it is not a nullspace of this
     * diagonal regression matrix and therefore must not alter this solve.
     */
    if(!failed) failed |= build_diagonal_problem(&A,&b,&x);
    if(!failed)
    {
        err = lis_solve_sequence(A,b,x,solver);
        if(err) failed = 1;
    }
    if(!failed)
    {
        err = lis_solver_get_status(solver,&status);
        if(err || status!=LIS_SUCCESS || solver->sequence_precon==NULL)
        {
            fprintf(stderr,
                "failed to create persistent sequence state: err=%d status=%d precon=%p\n",
                (int)err,(int)status,(void *)solver->sequence_precon);
            failed = 1;
        }
    }

    /* The compact report must match the existing individual diagnostics. */
    if(!failed)
    {
        err = lis_solver_get_report(solver,&report);
        if(err) failed = 1;
    }
    if(!failed)
    {
        lis_solver_get_status(solver,&status);
        lis_solver_get_iterex(solver,&iter,&iter_double,&iter_quad);
        lis_solver_get_residualnorm(solver,&residual);
        lis_solver_get_timeex(
            solver,&time,&itime,&ptime,&p_c_time,&p_i_time);

        if(report.status!=status ||
           report.iterations!=iter ||
           report.iterations_double!=iter_double ||
           report.iterations_quad!=iter_quad ||
           report.residual_norm!=residual ||
           report.time!=time ||
           report.iteration_time!=itime ||
           report.preconditioner_time!=ptime ||
           report.preconditioner_create_time!=p_c_time ||
           report.preconditioner_iteration_time!=p_i_time)
        {
            fprintf(stderr,"compact solve report mismatch\n");
            failed = 1;
        }
    }

    /*
     * Reapplying an equivalent complete profile must be a no-op.
     * Use a different textual ordering to verify comparison of the
     * parsed effective profile rather than the original string.
     */
    if(!failed)
    {
        saved_precon = solver->sequence_precon;
        solver->recycle_valid = LIS_TRUE;
        changed = -1;

        err = lis_solver_set_option_replace_if_changed(
            "-p ilu -restart 12 -tol 1e-7 -maxiter 77 -i fgmres",
            solver,
            &changed);

        if(err || changed!=LIS_FALSE ||
           solver->sequence_precon!=saved_precon ||
           solver->recycle_valid!=LIS_TRUE)
        {
            fprintf(stderr,
                "equivalent replacement invalidated reusable state\n");
            failed = 1;
        }
    }

    /* Typed solve-time controls must preserve reusable solver state. */
    if(!failed)
    {
        err = lis_solver_get_tolerance(solver,&tolerance);
        if(err || tolerance!=(LIS_REAL)1.0e-7)
            failed = 1;
        err = lis_solver_get_maxiter(solver,&maxiter);
        if(err || maxiter!=77)
            failed = 1;

        err = lis_solver_set_tolerance(solver,(LIS_REAL)2.0e-7);
        if(err) failed = 1;
        err = lis_solver_set_maxiter(solver,88);
        if(err) failed = 1;

        lis_solver_get_tolerance(solver,&tolerance);
        lis_solver_get_maxiter(solver,&maxiter);
        if(tolerance!=(LIS_REAL)2.0e-7 || maxiter!=88 ||
           solver->sequence_precon!=saved_precon ||
           solver->recycle_valid!=LIS_TRUE)
        {
            fprintf(stderr,
                "typed solve-time update mismatch or reusable state lost\n");
            failed = 1;
        }
    }

    /* Invalid typed updates must leave the active values unchanged. */
    if(!failed)
    {
        err = lis_solver_set_tolerance(solver,(LIS_REAL)-1.0);
        if(err==LIS_SUCCESS) failed = 1;
        err = lis_solver_set_tolerance(solver,(LIS_REAL)NAN);
        if(err==LIS_SUCCESS) failed = 1;
        err = lis_solver_set_tolerance(solver,(LIS_REAL)INFINITY);
        if(err==LIS_SUCCESS) failed = 1;
        err = lis_solver_set_tolerance(solver,(LIS_REAL)-INFINITY);
        if(err==LIS_SUCCESS) failed = 1;
        err = lis_solver_set_maxiter(solver,-1);
        if(err==LIS_SUCCESS) failed = 1;

        lis_solver_get_tolerance(solver,&tolerance);
        lis_solver_get_maxiter(solver,&maxiter);
        if(tolerance!=(LIS_REAL)2.0e-7 || maxiter!=88)
        {
            fprintf(stderr,"failed typed update changed active values\n");
            failed = 1;
        }
    }

    /* Install application-owned state only after the sequence solve. */
    if(!failed) err = lis_solver_set_monitor(solver,test_monitor,&monitor_cookie);
    if(err) failed = 1;
    if(!failed) err = lis_solver_set_convergence(solver,test_convergence,&convergence_cookie);
    if(err) failed = 1;

    if(!failed) err = lis_vector_create(LIS_COMM_WORLD,&q);
    if(err) failed = 1;
    if(!failed) err = lis_vector_set_size(q,0,4);
    if(err) failed = 1;
    if(!failed) err = lis_vector_set_all((LIS_SCALAR)0.5,q);
    if(err) failed = 1;
    basis[0] = q;

    if(!failed) err = lis_solver_set_nullspace(solver,1,basis);
    if(err) failed = 1;
    if(!failed) err = lis_solver_set_left_nullspace(solver,1,basis);
    if(err) failed = 1;
    if(!failed) err = lis_solver_set_near_nullspace(solver,1,basis);
    if(err) failed = 1;

    /* Also mark generic GCRO-DR recycle validity stale. */
    solver->recycle_valid = LIS_TRUE;

    if(!failed)
        err = lis_solver_set_option_replace("-i cg -p none",solver);
    if(err) failed = 1;

    lis_solver_get_solver(solver,&method);
    lis_solver_get_precon(solver,&precon);
    lis_solver_get_sequence_policy(solver,&policy);
    lis_solver_get_sequence_precon_max_age(solver,&age);
    lis_solver_get_num_threads(solver,&threads);
    lis_solver_get_precon_half(solver,&half_mode);

    if(method!=LIS_SOLVER_CG || precon!=LIS_PRECON_TYPE_NONE ||
       solver->options[LIS_OPTIONS_MAXITER]!=1000 ||
       solver->options[LIS_OPTIONS_RESTART]!=40 ||
       policy!=LIS_SEQUENCE_POLICY_AUTO || age!=3 || threads!=1 ||
       half_mode!=LIS_PRECON_HALF_ON)
    {
        fprintf(stderr,
            "default restoration/application policy mismatch: method=%d precon=%d maxiter=%d restart=%d policy=%d age=%d threads=%d half=%d\n",
            (int)method,(int)precon,
            (int)solver->options[LIS_OPTIONS_MAXITER],
            (int)solver->options[LIS_OPTIONS_RESTART],
            (int)policy,(int)age,(int)threads,(int)half_mode);
        failed = 1;
    }

    if(solver->monitor!=test_monitor || solver->monitor_data!=&monitor_cookie ||
       solver->convergence!=test_convergence ||
       solver->convergence_data!=&convergence_cookie ||
       solver->nullspace_dim!=1 || solver->nullspace==NULL ||
       solver->left_nullspace_dim!=1 || solver->left_nullspace==NULL ||
       solver->near_nullspace_dim!=1 || solver->near_nullspace==NULL ||
       solver->sequence_precon!=NULL || solver->precon!=NULL ||
       solver->sequence_precon_age!=0 || solver->recycle_valid!=LIS_FALSE)
    {
        fprintf(stderr,
            "application-state preservation or sequence-state invalidation mismatch\n");
        failed = 1;
    }

    /* Invalid replacement must leave the installed profile and app state unchanged. */
    err = lis_solver_set_option_replace(
        "-i this_solver_does_not_exist",solver);
    if(err==LIS_SUCCESS)
    {
        fprintf(stderr,"invalid replacement unexpectedly succeeded\n");
        failed = 1;
    }

    lis_solver_get_solver(solver,&method);
    lis_solver_get_precon(solver,&precon);
    lis_solver_get_sequence_policy(solver,&policy);
    lis_solver_get_num_threads(solver,&threads);
    lis_solver_get_precon_half(solver,&half_mode);
    if(method!=LIS_SOLVER_CG || precon!=LIS_PRECON_TYPE_NONE ||
       solver->options[LIS_OPTIONS_MAXITER]!=1000 ||
       solver->monitor!=test_monitor || solver->nullspace_dim!=1 ||
       solver->left_nullspace_dim!=1 || solver->near_nullspace_dim!=1 ||
       policy!=LIS_SEQUENCE_POLICY_AUTO || threads!=1 ||
       half_mode!=LIS_PRECON_HALF_ON)
    {
        fprintf(stderr,"failed replacement modified active solver/application state\n");
        failed = 1;
    }

    /* A genuinely different effective profile must report a change. */
    if(!failed)
    {
        changed = -1;
        err = lis_solver_set_option_replace_if_changed(
            "-i cg -p none -maxiter 123",
            solver,
            &changed);

        if(err || changed!=LIS_TRUE ||
           solver->options[LIS_OPTIONS_MAXITER]!=123)
        {
            fprintf(stderr,
                "changed replacement was not reported correctly\n");
            failed = 1;
        }
    }

    /* Invalid input must report no completed change. */
    if(!failed)
    {
        changed = -1;
        err = lis_solver_set_option_replace_if_changed(
            "-i this_solver_still_does_not_exist",
            solver,
            &changed);

        if(err==LIS_SUCCESS || changed!=LIS_FALSE ||
           solver->options[LIS_OPTIONS_MAXITER]!=123)
        {
            fprintf(stderr,
                "invalid conditional replacement changed active profile\n");
            failed = 1;
        }
    }

    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    if(q) lis_vector_destroy(q);
    lis_solver_destroy(solver);
    lis_finalize();

    if(failed)
    {
        fprintf(stderr,"SOLVER_OPTION_REPLACE FAILED\n");
        return 1;
    }

    printf("SOLVER_OPTION_REPLACE PASSED\n");
    return 0;
}
