#include "lis.h"
#include "lis_solver_options.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

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

static int check_report_snapshot(
    LIS_SOLVER solver,LIS_INT expected_status,const char *stage)
{
    LIS_SOLVE_REPORT snapshot;
    LIS_REAL residual;
    LIS_INT status,iter,iter_double,iter_quad;
    double time,itime,ptime,p_c_time,p_i_time;
    LIS_INT err;

    err=lis_solver_get_report(solver,&snapshot);
    if(err) return 1;
    if(lis_solver_get_status(solver,&status) ||
       lis_solver_get_iterex(solver,&iter,&iter_double,&iter_quad) ||
       lis_solver_get_residualnorm(solver,&residual) ||
       lis_solver_get_timeex(solver,&time,&itime,&ptime,&p_c_time,&p_i_time))
        return 1;
    if(snapshot.status!=expected_status || snapshot.status!=status ||
       snapshot.iterations!=iter ||
       snapshot.iterations_double!=iter_double ||
       snapshot.iterations_quad!=iter_quad ||
       snapshot.residual_norm!=residual ||
       snapshot.time!=time || snapshot.iteration_time!=itime ||
       snapshot.preconditioner_time!=ptime ||
       snapshot.preconditioner_create_time!=p_c_time ||
       snapshot.preconditioner_iteration_time!=p_i_time)
    {
        fprintf(stderr,"report snapshot mismatch: %s\n",stage);
        return 1;
    }
    return 0;
}

/* Isolate the typed API from legacy text-option parsing for leak checks.
 * The full regression is still run separately without suppressing failures.
 */
static int test_typed_api_isolated(int argc,char **argv)
{
    LIS_SOLVER solver=NULL;
    LIS_MATRIX A=NULL;
    LIS_VECTOR b=NULL,x=NULL;
    LIS_REAL tol;
    LIS_INT maxiter,err;
    int failed=0;

    argc=1; /* --typed-only is a test harness flag, not a LIS option. */
    err=lis_initialize(&argc,&argv);
    if(err) return 1;
    err=lis_solver_create(&solver);
    if(err){failed=1;goto done;}
    if(check_report_snapshot(solver,LIS_SUCCESS,"typed-only-before-solve"))
        failed=1;
    if(lis_solver_set_tolerance(solver,(LIS_REAL)1.0e-8) ||
       lis_solver_set_maxiter(solver,15))
        failed=1;
    if(lis_solver_set_tolerance(solver,(LIS_REAL)0.0) ||
       lis_solver_set_maxiter(solver,0) ||
       lis_solver_get_tolerance(solver,&tol) ||
       lis_solver_get_maxiter(solver,&maxiter) ||
       tol!=(LIS_REAL)0.0 || maxiter!=0)
        failed=1;
    if(lis_solver_set_tolerance(solver,(LIS_REAL)NAN)==LIS_SUCCESS ||
       lis_solver_set_tolerance(solver,(LIS_REAL)INFINITY)==LIS_SUCCESS ||
       lis_solver_set_maxiter(solver,-1)==LIS_SUCCESS ||
       lis_solver_get_tolerance(solver,&tol) ||
       lis_solver_get_maxiter(solver,&maxiter) ||
       tol!=(LIS_REAL)0.0 || maxiter!=0)
        failed=1;
    if(!failed && build_diagonal_problem(&A,&b,&x))
        failed=1;
    if(!failed)
    {
        solver->options[LIS_OPTIONS_SOLVER]=0;
        err=lis_solve(A,b,x,solver);
        if(err==LIS_SUCCESS ||
           check_report_snapshot(solver,err,"typed-only-after-failure"))
            failed=1;
    }
done:
    if(x)lis_vector_destroy(x);
    if(b)lis_vector_destroy(b);
    if(A)lis_matrix_destroy(A);
    if(solver)lis_solver_destroy(solver);
    lis_finalize();
    printf("SOLVER_TYPED_API_ISOLATED %s\n",failed?"FAILED":"PASSED");
    return failed;
}


/*
 * GCRO-DR options have dedicated solver fields, outside the normal option
 * arrays.  Replacing a complete profile must copy all of them, detect
 * changes to any one, restore omitted defaults and avoid invalidating
 * reusable state for an equivalent effective profile.
 */
static int check_gcrodr_profile_settings(void)
{
    LIS_SOLVER s=NULL;
    LIS_INT err,changed;
    int failed=0;

    err=lis_solver_create(&s);
    if(err) return 1;

    err=lis_solver_set_option_replace(
        "-i gcrodr -p none -gcrodr_keep 5 "
        "-gcrodr_keep_quality 0.82 -gcrodr_reject_quality 0.22",s);
    if(err || s->recycle_target_dim!=5 ||
       s->recycle_keep_threshold!=(LIS_REAL)0.82 ||
       s->recycle_reject_threshold!=(LIS_REAL)0.22)
    {
        fprintf(stderr,"GCRO-DR profile fields not copied\n");
        failed=1; goto done;
    }

    /* Effective profile equality must not clear an existing recycle. */
    s->recycle_valid=LIS_TRUE;
    changed=-1;
    err=lis_solver_set_option_replace_if_changed(
        "-gcrodr_reject_quality 0.22 -p none -gcrodr_keep 5 "
        "-gcrodr_keep_quality 0.82 -i gcrodr",s,&changed);
    if(err || changed!=LIS_FALSE || s->recycle_valid!=LIS_TRUE)
    {
        fprintf(stderr,"equivalent GCRO-DR profile invalidated recycle\n");
        failed=1; goto done;
    }

    /* Change only the number of vectors. */
    changed=-1;
    err=lis_solver_set_option_replace_if_changed(
        "-i gcrodr -p none -gcrodr_keep 6 "
        "-gcrodr_keep_quality 0.82 -gcrodr_reject_quality 0.22",s,&changed);
    if(err || changed!=LIS_TRUE || s->recycle_valid!=LIS_FALSE ||
       s->recycle_target_dim!=6)
    {
        fprintf(stderr,"GCRO-DR keep-only profile change not detected\n");
        failed=1; goto done;
    }

    /* Change only the keep-quality policy. */
    s->recycle_valid=LIS_TRUE;
    changed=-1;
    err=lis_solver_set_option_replace_if_changed(
        "-i gcrodr -p none -gcrodr_keep 6 "
        "-gcrodr_keep_quality 0.75 -gcrodr_reject_quality 0.22",s,&changed);
    if(err || changed!=LIS_TRUE || s->recycle_valid!=LIS_FALSE ||
       s->recycle_keep_threshold!=(LIS_REAL)0.75)
    {
        fprintf(stderr,"GCRO-DR keep-quality-only profile change not detected\n");
        failed=1; goto done;
    }

    /* Change only the reject-quality policy. */
    s->recycle_valid=LIS_TRUE;
    changed=-1;
    err=lis_solver_set_option_replace_if_changed(
        "-i gcrodr -p none -gcrodr_keep 6 "
        "-gcrodr_keep_quality 0.75 -gcrodr_reject_quality 0.10",s,&changed);
    if(err || changed!=LIS_TRUE || s->recycle_valid!=LIS_FALSE ||
       s->recycle_reject_threshold!=(LIS_REAL)0.10)
    {
        fprintf(stderr,"GCRO-DR reject-quality-only profile change not detected\n");
        failed=1; goto done;
    }

    /* Parse failure must leave the active profile and recycle intact. */
    s->recycle_valid=LIS_TRUE;
    changed=-1;
    err=lis_solver_set_option_replace_if_changed(
        "-i gcrodr -p none -gcrodr_keep 9 "
        "-gcrodr_keep_quality 1.7 -gcrodr_reject_quality 0.10",s,&changed);
    if(err==LIS_SUCCESS || changed!=LIS_FALSE ||
       s->recycle_valid!=LIS_TRUE || s->recycle_target_dim!=6 ||
       s->recycle_keep_threshold!=(LIS_REAL)0.75 ||
       s->recycle_reject_threshold!=(LIS_REAL)0.10)
    {
        fprintf(stderr,"invalid GCRO-DR profile changed active settings\n");
        failed=1; goto done;
    }

    /* Omitted options in a full profile restore documented defaults. */
    changed=-1;
    err=lis_solver_set_option_replace_if_changed("-i gcrodr -p none",s,&changed);
    if(err || changed!=LIS_TRUE || s->recycle_valid!=LIS_FALSE ||
       s->recycle_target_dim!=8 ||
       s->recycle_keep_threshold!=(LIS_REAL)0.90 ||
       s->recycle_reject_threshold!=(LIS_REAL)0.35)
    {
        fprintf(stderr,"GCRO-DR omitted defaults not restored\n");
        failed=1; goto done;
    }

    s->recycle_valid=LIS_TRUE;
    changed=-1;
    err=lis_solver_set_option_replace_if_changed("-p none -i gcrodr",s,&changed);
    if(err || changed!=LIS_FALSE || s->recycle_valid!=LIS_TRUE)
    {
        fprintf(stderr,"GCRO-DR default profile unexpectedly invalidated\n");
        failed=1;
    }

done:
    lis_solver_destroy(s);
    printf("GCRODR_OPTION_PROFILE %s\n",failed?"FAILED":"PASSED");
    return failed;
}

int main(int argc, char **argv)
{
    if(argc>1 && strcmp(argv[1],"--typed-only")==0)
        return test_typed_api_isolated(argc,argv);

    LIS_SOLVER solver = NULL;
    LIS_SOLVER failing_solver = NULL;
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

    /* Report is a getter snapshot, including before the first solve. */
    if(check_report_snapshot(solver,LIS_SUCCESS,"before-first-solve"))
        failed=1;
    if(lis_solver_get_report(NULL,&report)==LIS_SUCCESS ||
       lis_solver_get_report(solver,NULL)==LIS_SUCCESS)
        failed=1;

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
        err=lis_solver_create(&failing_solver);
        if(err) failed=1;
    }
    if(!failed)
    {
        /* Force a rejected solver type through the normal solve path.
         * The report must describe the failed attempt, not imply success. */
        failing_solver->options[LIS_OPTIONS_SOLVER]=0;
        err=lis_solve(A,b,x,failing_solver);
        if(err==LIS_SUCCESS ||
           check_report_snapshot(failing_solver,err,"after-failed-solve"))
        {
            fprintf(stderr,"failed-solve report mismatch\n");
            failed=1;
        }
        lis_solver_destroy(failing_solver);
        failing_solver=NULL;
    }
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

    /* Zero retains LIS option semantics; neither setter destroys state. */
    if(!failed)
    {
        err=lis_solver_set_tolerance(solver,(LIS_REAL)0.0);
        if(err) failed=1;
        err=lis_solver_set_maxiter(solver,0);
        if(err) failed=1;
        lis_solver_get_tolerance(solver,&tolerance);
        lis_solver_get_maxiter(solver,&maxiter);
        if(tolerance!=(LIS_REAL)0.0 || maxiter!=0 ||
           solver->sequence_precon!=saved_precon ||
           solver->recycle_valid!=LIS_TRUE)
            failed=1;
        err=lis_solver_set_tolerance(solver,(LIS_REAL)2.0e-7);
        if(err) failed=1;
        err=lis_solver_set_maxiter(solver,88);
        if(err) failed=1;
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

    /* Isolated GCRO-DR profile regression, independent of the solve above. */
    failed |= check_gcrodr_profile_settings();

    if(x) lis_vector_destroy(x);
    if(b) lis_vector_destroy(b);
    if(A) lis_matrix_destroy(A);
    if(q) lis_vector_destroy(q);
    if(failing_solver) lis_solver_destroy(failing_solver);
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
