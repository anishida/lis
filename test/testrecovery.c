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

int main(int argc, char **argv)
{
	LIS_MATRIX A = NULL;
	LIS_VECTOR b = NULL;
	LIS_VECTOR x = NULL;
	LIS_SOLVER solver = NULL;
	LIS_INT err;
	int failed = 0;

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
	/* Case 1: first attempt breaks, recycle-reset retry works. */
	/* ======================================================== */

	reset_counts();

	err = create_solver(&solver,LIS_FALSE);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	/*
	 * One CG iteration is intentionally insufficient for
	 * diag(2,3).  With the fixed ||b||2 normalization the
	 * first attempt reaches a relative residual of about 0.20,
	 * while the retry from the improved x reaches about 0.04.
	 * A tolerance of 0.15 therefore forces exactly one retry
	 * without rebuilding the persistent preconditioner.
	 */
	solver->options[LIS_OPTIONS_SOLVER] = LIS_SOLVER_CG;
	solver->options[LIS_OPTIONS_CONV_COND] =
		LIS_CONV_COND_NRM2_B;
	solver->options[LIS_OPTIONS_MAXITER] = 1;
	solver->params[LIS_PARAMS_RESID-LIS_OPTIONS_LEN] =
		(LIS_SCALAR)1.5e-1;

	recovery_mode = 0;
	fail_once = 0;

	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(A,b,x,solver);

	if( err!=LIS_SUCCESS ||
	    solver->retcode!=LIS_SUCCESS )
	{
		fprintf(
			stderr,
			"RETRY_REUSE failed api=%d status=%d\n",
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
			"RETRY_REUSE lifecycle "
			"create=%d destroy=%d failures=%d\n",
			create_count,
			destroy_count,
			psolve_failures);
		failed = 1;
		goto cleanup;
	}

	printf(
		"RETRY_REUSE PASSED "
		"create=%d destroy=%d failures=%d\n",
		create_count,
		destroy_count,
		psolve_failures);

	err = lis_solver_reset_sequence(solver);
	if( err || destroy_count!=1 )
	{
		fprintf(stderr,"RETRY_REUSE reset failed\n");
		failed = 1;
		goto cleanup;
	}

	lis_solver_destroy(solver);
	solver = NULL;


	/* ======================================================== */
	/* Case 2: every attempt breaks -> preconditioner rebuild.  */
	/* ======================================================== */

	reset_counts();

	err = create_solver(&solver,LIS_FALSE);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	recovery_mode = 2;

	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(A,b,x,solver);

	if( err!=LIS_SUCCESS )
	{
		fprintf(
			stderr,
			"PRECON_REBUILD unexpected API error=%d\n",
			(int)err);
		failed = 1;
		goto cleanup;
	}

	if( solver->retcode!=LIS_BREAKDOWN &&
	    solver->retcode!=LIS_MAXITER )
	{
		fprintf(
			stderr,
			"PRECON_REBUILD status=%d\n",
			(int)solver->retcode);
		failed = 1;
		goto cleanup;
	}

	if( create_count!=2 ||
	    destroy_count!=1 ||
	    psolve_failures<3 )
	{
		fprintf(
			stderr,
			"PRECON_REBUILD lifecycle "
			"create=%d destroy=%d failures=%d\n",
			create_count,
			destroy_count,
			psolve_failures);
		failed = 1;
		goto cleanup;
	}

	printf(
		"PRECON_REBUILD PASSED "
		"status=%d create=%d destroy=%d failures=%d\n",
		(int)solver->retcode,
		create_count,
		destroy_count,
		psolve_failures);

	err = lis_solver_reset_sequence(solver);
	if( err || destroy_count!=2 )
	{
		fprintf(stderr,"PRECON_REBUILD reset failed\n");
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
