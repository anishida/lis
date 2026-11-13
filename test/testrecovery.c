#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"
#include "lis_solver.h"

#include <stdio.h>

static int create_count = 0;
static int destroy_count = 0;
static int psolve_failures = 0;

static int recovery_mode = 0;
static int fail_once = 0;

static LIS_INT
recovery_create(LIS_SOLVER solver, LIS_PRECON precon)
{
	(void)solver;
	(void)precon;

	create_count++;
	return LIS_SUCCESS;
}

static LIS_INT
recovery_psolve(
	LIS_SOLVER solver,
	LIS_VECTOR b,
	LIS_VECTOR x)
{
	(void)solver;

	if( recovery_mode==1 && fail_once )
	{
		fail_once = 0;
		psolve_failures++;
		return LIS_BREAKDOWN;
	}

	if( recovery_mode==2 )
	{
		psolve_failures++;
		return LIS_BREAKDOWN;
	}

	return lis_vector_copy(b,x);
}

static LIS_INT
recovery_destroy(LIS_PRECON precon)
{
	(void)precon;

	destroy_count++;
	return LIS_SUCCESS;
}

static LIS_INT
create_matrix(LIS_MATRIX *A)
{
	LIS_INT err;
	LIS_INT i,is,ie;

	err = lis_matrix_create(LIS_COMM_WORLD,A);
	if( err ) return err;

	err = lis_matrix_set_size(*A,0,2);
	if( err ) return err;

	err = lis_matrix_get_range(*A,&is,&ie);
	if( err ) return err;

	for(i=is;i<ie;i++)
	{
		err = lis_matrix_set_value(
			LIS_INS_VALUE,
			i,i,
			(LIS_SCALAR)(i==0 ? 2.0 : 3.0),
			*A);
		if( err ) return err;
	}

	err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
	if( err ) return err;

	return lis_matrix_assemble(*A);
}

static LIS_INT
create_solver(
	LIS_SOLVER *solver,
	LIS_INT bad_option)
{
	LIS_INT err;

	err = lis_solver_create(solver);
	if( err ) return err;

	if( bad_option )
	{
		return lis_solver_set_option(
			"-i gmres -p recs7 "
			"-scale none -print none "
			"-restart 2 -maxiter 20 "
			"-maxiter_noimp -1 -tol 1e-12",
			*solver);
	}

	return lis_solver_set_option(
		"-i gmres -p recs7 "
		"-scale none -print none "
		"-restart 2 -maxiter 20 "
		"-maxiter_noimp 0 -tol 1e-12",
		*solver);
}

static void
reset_counts(void)
{
	create_count = 0;
	destroy_count = 0;
	psolve_failures = 0;
	recovery_mode = 0;
	fail_once = 0;
}

static int
check_sequence_diag(
    const char *name,
    LIS_SOLVER solver,
    LIS_INT attempts,
    LIS_INT retries,
    LIS_INT resets,
    LIS_INT rebuilds,
    LIS_INT first_a,
    LIS_INT first_b,
    LIS_INT final_status,
    LIS_INT initial_warm,
    LIS_INT last_warm,
    LIS_INT initial_valid,
    LIS_INT true_valid,
    LIS_UNSIGNED_INT total_solves,
    LIS_UNSIGNED_INT total_attempts,
    LIS_UNSIGNED_INT total_retries,
    LIS_UNSIGNED_INT total_resets,
    LIS_UNSIGNED_INT total_rebuilds)
{
    LIS_INT err;
    LIS_SEQUENCE_DIAGNOSTICS d;

    err = lis_solver_get_sequence_diagnostics(
        solver,&d);

    if( err )
    {
        fprintf(
            stderr,
            "%s diagnostics getter failed rc=%d\n",
            name,
            (int)err);
        return 1;
    }

    if( d.attempts!=attempts ||
        d.retries!=retries ||
        d.recycle_resets!=resets ||
        d.precon_rebuilds!=rebuilds )
    {
        fprintf(
            stderr,
            "%s per-call mismatch "
            "attempts=%d/%d retries=%d/%d "
            "resets=%d/%d rebuilds=%d/%d\n",
            name,
            (int)d.attempts,
            (int)attempts,
            (int)d.retries,
            (int)retries,
            (int)d.recycle_resets,
            (int)resets,
            (int)d.precon_rebuilds,
            (int)rebuilds);
        return 1;
    }

    if( d.first_retcode!=first_a &&
        d.first_retcode!=first_b )
    {
        fprintf(
            stderr,
            "%s first status=%d expected=%d/%d\n",
            name,
            (int)d.first_retcode,
            (int)first_a,
            (int)first_b);
        return 1;
    }

    if( d.final_retcode!=final_status )
    {
        fprintf(
            stderr,
            "%s final status=%d expected=%d\n",
            name,
            (int)d.final_retcode,
            (int)final_status);
        return 1;
    }

    if( initial_warm>=0 &&
        d.initial_warm_start_kept!=initial_warm )
    {
        fprintf(
            stderr,
            "%s initial warm=%d expected=%d\n",
            name,
            (int)d.initial_warm_start_kept,
            (int)initial_warm);
        return 1;
    }

    if( last_warm>=0 &&
        d.last_warm_start_kept!=last_warm )
    {
        fprintf(
            stderr,
            "%s last warm=%d expected=%d\n",
            name,
            (int)d.last_warm_start_kept,
            (int)last_warm);
        return 1;
    }

    if( (initial_valid>=0 &&
         d.initial_residual_valid!=initial_valid) ||
        (true_valid>=0 &&
         d.true_residual_valid!=true_valid) )
    {
        fprintf(
            stderr,
            "%s residual validity=%d/%d expected=%d/%d\n",
            name,
            (int)d.initial_residual_valid,
            (int)d.true_residual_valid,
            (int)initial_valid,
            (int)true_valid);
        return 1;
    }

    if( d.total_solves!=total_solves ||
        d.total_attempts!=total_attempts ||
        d.total_retries!=total_retries ||
        d.total_recycle_resets!=total_resets ||
        d.total_precon_rebuilds!=total_rebuilds )
    {
        fprintf(
            stderr,
            "%s lifetime mismatch "
            "solves=%llu/%llu "
            "attempts=%llu/%llu "
            "retries=%llu/%llu "
            "resets=%llu/%llu "
            "rebuilds=%llu/%llu\n",
            name,
            (unsigned long long)d.total_solves,
            (unsigned long long)total_solves,
            (unsigned long long)d.total_attempts,
            (unsigned long long)total_attempts,
            (unsigned long long)d.total_retries,
            (unsigned long long)total_retries,
            (unsigned long long)d.total_recycle_resets,
            (unsigned long long)total_resets,
            (unsigned long long)d.total_precon_rebuilds,
            (unsigned long long)total_rebuilds);
        return 1;
    }

    printf(
        "%s DIAG attempts=%d retries=%d "
        "resets=%d rebuilds=%d "
        "first=%d final=%d warm=%d/%d valid=%d/%d "
        "initial=%.6e final_resid=%.6e "
        "totals=%llu/%llu/%llu/%llu/%llu\n",
        name,
        (int)d.attempts,
        (int)d.retries,
        (int)d.recycle_resets,
        (int)d.precon_rebuilds,
        (int)d.first_retcode,
        (int)d.final_retcode,
        (int)d.initial_warm_start_kept,
        (int)d.last_warm_start_kept,
        (int)d.initial_residual_valid,
        (int)d.true_residual_valid,
        (double)d.initial_residual,
        (double)d.true_residual,
        (unsigned long long)d.total_solves,
        (unsigned long long)d.total_attempts,
        (unsigned long long)d.total_retries,
        (unsigned long long)d.total_recycle_resets,
        (unsigned long long)d.total_precon_rebuilds);

    return 0;
}

int main(int argc, char **argv)
{
	LIS_MATRIX A = NULL;
	LIS_VECTOR b = NULL;
	LIS_VECTOR x = NULL;
	LIS_SOLVER solver = NULL;
	LIS_INT err;
	int failed = 0;

        LIS_INT method_index;
        char method_name[2][8] = {"GMRES","FGMRES"};
        char method_option[2][16] = {"-i gmres","-i fgmres"};

	err = lis_initialize(&argc,&argv);
	if( err ) return 1;

	err = lis_precon_register_ex(
		"recs7",
		recovery_create,
		recovery_psolve,
		recovery_psolve,
		recovery_destroy);
	if( err )
	{
		fprintf(stderr,"preconditioner registration failed\n");
		failed = 1;
		goto cleanup;
	}

	err = create_matrix(&A);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(A,&b);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(A,&x);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_set_all((LIS_SCALAR)2.0,b);


        /* ======================================================== */
        /* Cases 1-2: ordinary GMRES-family recovery policy.    */
        /* ======================================================== */

        for(method_index=0;method_index<2;method_index++)
        {
                char label[64];

                /* Fresh preconditioner: retryable failure must not retry. */

                reset_counts();

                err = create_solver(&solver,LIS_FALSE);
                if( err )
                {
                        failed = 1;
                        goto cleanup;
                }

                err = lis_solver_set_option(
                        method_option[method_index],
                        solver);
                if( err )
                {
                        failed = 1;
                        goto cleanup;
                }

                recovery_mode = 2;
                fail_once = 0;

                lis_vector_set_all((LIS_SCALAR)0.0,x);

                err = lis_solve_sequence(A,b,x,solver);

                snprintf(
                        label,
                        sizeof(label),
                        "FRESH_%s",
                        method_name[method_index]);

                if( err!=LIS_SUCCESS ||
                    (solver->retcode!=LIS_BREAKDOWN &&
                     solver->retcode!=LIS_MAXITER) )
                {
                        fprintf(
                                stderr,
                                "%s unexpected api=%d status=%d\n",
                                label,
                                (int)err,
                                (int)solver->retcode);
                        failed = 1;
                        goto cleanup;
                }

                if( create_count!=1 ||
                    destroy_count!=0 ||
                    psolve_failures<1 )
                {
                        fprintf(
                                stderr,
                                "%s lifecycle create=%d destroy=%d failures=%d\n",
                                label,
                                create_count,
                                destroy_count,
                                psolve_failures);
                        failed = 1;
                        goto cleanup;
                }

                if( check_sequence_diag(
                        label,
                        solver,
                        1,0,0,0,
                        LIS_BREAKDOWN,LIS_MAXITER,
                        solver->retcode,
                        -1,-1,
                        -1,-1,
                        1,1,0,0,0) )
                {
                        failed = 1;
                        goto cleanup;
                }

                printf("%s PASSED\n",label);

                err = lis_solver_reset_sequence(solver);
                if( err || destroy_count!=1 )
                {
                        fprintf(stderr,"%s reset failed\n",label);
                        failed = 1;
                        goto cleanup;
                }

                lis_solver_destroy(solver);
                solver = NULL;


                /* Inherited preconditioner: rebuild exactly once. */

                reset_counts();

                err = create_solver(&solver,LIS_FALSE);
                if( err )
                {
                        failed = 1;
                        goto cleanup;
                }

                err = lis_solver_set_option(
                        method_option[method_index],
                        solver);
                if( err )
                {
                        failed = 1;
                        goto cleanup;
                }

                recovery_mode = 0;
                fail_once = 0;

                lis_vector_set_all((LIS_SCALAR)0.0,x);

                err = lis_solve_sequence(A,b,x,solver);

                if( err!=LIS_SUCCESS ||
                    solver->retcode!=LIS_SUCCESS ||
                    solver->sequence_precon==NULL )
                {
                        fprintf(
                                stderr,
                                "INHERITED_%s setup failed "
                                "api=%d status=%d precon=%p\n",
                                method_name[method_index],
                                (int)err,
                                (int)solver->retcode,
                                (void *)solver->sequence_precon);
                        failed = 1;
                        goto cleanup;
                }

                recovery_mode = 2;

                lis_vector_set_all((LIS_SCALAR)0.0,x);

                err = lis_solve_sequence(A,b,x,solver);

                snprintf(
                        label,
                        sizeof(label),
                        "INHERITED_%s",
                        method_name[method_index]);

                if( err!=LIS_SUCCESS ||
                    (solver->retcode!=LIS_BREAKDOWN &&
                     solver->retcode!=LIS_MAXITER) )
                {
                        fprintf(
                                stderr,
                                "%s unexpected api=%d status=%d\n",
                                label,
                                (int)err,
                                (int)solver->retcode);
                        failed = 1;
                        goto cleanup;
                }

                if( create_count!=2 ||
                    destroy_count!=1 ||
                    psolve_failures<2 )
                {
                        fprintf(
                                stderr,
                                "%s lifecycle create=%d destroy=%d failures=%d\n",
                                label,
                                create_count,
                                destroy_count,
                                psolve_failures);
                        failed = 1;
                        goto cleanup;
                }

                if( check_sequence_diag(
                        label,
                        solver,
                        2,1,0,1,
                        LIS_BREAKDOWN,LIS_MAXITER,
                        solver->retcode,
                        -1,-1,
                        -1,-1,
                        2,3,1,0,1) )
                {
                        failed = 1;
                        goto cleanup;
                }

                printf("%s PASSED\n",label);

                err = lis_solver_reset_sequence(solver);
                if( err || destroy_count!=2 )
                {
                        fprintf(stderr,"%s reset failed\n",label);
                        failed = 1;
                        goto cleanup;
                }

                lis_solver_destroy(solver);
                solver = NULL;
        }


        /* ======================================================== */
        /* Case 3: GCRO-DR owns reusable Krylov recycle state.     */
        /* ======================================================== */

        reset_counts();

        err = create_solver(&solver,LIS_FALSE);
        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        {
                char gcrodr_option[] =
                        "-i gcrodr -gcrodr_keep 1";

                err = lis_solver_set_option(
                        gcrodr_option,
                        solver);
        }

        if( err )
        {
                failed = 1;
                goto cleanup;
        }

        solver->options[LIS_OPTIONS_MAXITER] = 1;

        recovery_mode = 0;
        fail_once = 0;

        lis_vector_set_all((LIS_SCALAR)0.0,x);

        err = lis_solve_sequence(A,b,x,solver);

        if( err!=LIS_SUCCESS )
        {
                fprintf(
                        stderr,
                        "GCRODR_RECOVERY api=%d status=%d\n",
                        (int)err,
                        (int)solver->retcode);
                failed = 1;
                goto cleanup;
        }

        if( create_count!=1 ||
            destroy_count!=0 ||
            psolve_failures!=0 )
        {
                fprintf(
                        stderr,
                        "GCRODR_RECOVERY lifecycle "
                        "create=%d destroy=%d failures=%d\n",
                        create_count,
                        destroy_count,
                        psolve_failures);
                failed = 1;
                goto cleanup;
        }

        if( check_sequence_diag(
                "GCRODR_RECOVERY",
                solver,
                2,1,1,0,
                LIS_BREAKDOWN,LIS_MAXITER,
                solver->retcode,
                -1,-1,
                -1,-1,
                1,2,1,1,0) )
        {
                failed = 1;
                goto cleanup;
        }

        printf(
                "GCRODR_RECOVERY PASSED status=%d\n",
                (int)solver->retcode);

        err = lis_solver_reset_sequence(solver);
        if( err || destroy_count!=1 )
        {
                fprintf(stderr,"GCRODR_RECOVERY reset failed\n");
                failed = 1;
                goto cleanup;
        }

        lis_solver_destroy(solver);
        solver = NULL;


	/* ======================================================== */
	/* Case 3: configuration error is not retryable.            */
	/* ======================================================== */

	reset_counts();

	err = create_solver(&solver,LIS_TRUE);
	if( err )
	{
		/*
		 * Current option parser accepts this value and the
		 * kernel rejects it.  If parser validation becomes
		 * stricter later, that is also a valid no-retry result.
		 */
		printf(
			"NONRETRYABLE rejected during option setup rc=%d\n",
			(int)err);

		lis_solver_destroy(solver);
		solver = NULL;
	}
	else
	{
		recovery_mode = 0;

		lis_vector_set_all((LIS_SCALAR)0.0,x);

		err = lis_solve_sequence(A,b,x,solver);

		if( err!=LIS_ERR_ILL_ARG ||
		    solver->retcode!=LIS_ERR_ILL_ARG )
		{
			fprintf(
				stderr,
				"NONRETRYABLE expected ILL_ARG "
				"api=%d status=%d\n",
				(int)err,
				(int)solver->retcode);
			failed = 1;
			goto cleanup;
		}

        if( check_sequence_diag(
                "NONRETRYABLE",
                solver,
                1,0,0,0,
                LIS_ERR_ILL_ARG,LIS_ERR_ILL_ARG,
                LIS_ERR_ILL_ARG,
                LIS_FALSE,LIS_FALSE,
                LIS_TRUE,LIS_FALSE,
                1,1,0,0,0) )
        {
            failed = 1;
            goto cleanup;
        }

		if( create_count!=1 ||
		    destroy_count!=0 ||
		    psolve_failures!=0 )
		{
			fprintf(
				stderr,
				"NONRETRYABLE retried unexpectedly "
				"create=%d destroy=%d failures=%d\n",
				create_count,
				destroy_count,
				psolve_failures);
			failed = 1;
			goto cleanup;
		}

		printf(
			"NONRETRYABLE PASSED "
			"create=%d destroy=%d\n",
			create_count,
			destroy_count);

		err = lis_solver_reset_sequence(solver);
		if( err || destroy_count!=1 )
		{
			fprintf(stderr,"NONRETRYABLE reset failed\n");
			failed = 1;
			goto cleanup;
		}

		lis_solver_destroy(solver);
		solver = NULL;
	}


	printf("LIS_SEQUENCE_RECOVERY_S7 PASSED\n");


cleanup:

	if( solver ) lis_solver_destroy(solver);
	if( x ) lis_vector_destroy(x);
	if( b ) lis_vector_destroy(b);
	if( A ) lis_matrix_destroy(A);

	lis_finalize();

	return failed ? 1 : 0;
}
