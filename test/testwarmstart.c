#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_solver.h"

#include <math.h>
#include <stdio.h>

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
			(LIS_SCALAR)2.0,
			*A);
		if( err ) return err;
	}

	err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
	if( err ) return err;

	return lis_matrix_assemble(*A);
}

static int
check_small(LIS_REAL value, LIS_REAL limit, const char *label)
{
	if( value!=value || fabs(value)>limit )
	{
		fprintf(
			stderr,
			"%s value=%e limit=%e\n",
			label,
			(double)value,
			(double)limit);
		return 1;
	}

	return 0;
}

int main(int argc, char **argv)
{
	LIS_MATRIX A = NULL;
	LIS_VECTOR b = NULL;
	LIS_VECTOR x = NULL;
	LIS_SOLVER solver = NULL;
	LIS_REAL rnorm,bnorm;
	LIS_INT err;
	LIS_INT saved_initguess_option;
	int failed = 0;

	err = lis_initialize(&argc,&argv);
	if( err ) return 1;

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

	err = lis_solver_create(&solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_set_option(
		"-i gmres -p none "
		"-scale none -print none "
		"-restart 2 -maxiter 20 "
		"-maxiter_noimp 0 -tol 1e-12",
		solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	saved_initguess_option =
		solver->options[LIS_OPTIONS_INITGUESS_ZEROS];

	/*
	 * ZERO:
	 * x=0 gives ||b-Ax|| == ||b||, so AUTO must reject
	 * the candidate and use the zero initial guess.
	 */
	lis_vector_set_all((LIS_SCALAR)2.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(A,b,x,solver);
	if( err || solver->retcode!=LIS_SUCCESS )
	{
		fprintf(stderr,"AUTO_ZERO solve failed\n");
		failed = 1;
		goto cleanup;
	}

	if( solver->sequence_warm_start_kept!=LIS_FALSE )
	{
		fprintf(stderr,"AUTO_ZERO unexpectedly kept x\n");
		failed = 1;
		goto cleanup;
	}

	if( fabs(
		solver->sequence_initial_residual -
		solver->sequence_rhs_norm) > 1.0e-12 )
	{
		fprintf(stderr,"AUTO_ZERO residual/rhs mismatch\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_small(
		solver->sequence_true_relative_residual,
		1.0e-12,
		"AUTO_ZERO true residual");
	if( failed ) goto cleanup;

	printf(
		"AUTO_ZERO initial=%e rhs=%e final=%e\n",
		(double)solver->sequence_initial_residual,
		(double)solver->sequence_rhs_norm,
		(double)solver->sequence_true_relative_residual);

	/*
	 * KEEP:
	 * exact x=1 for A=2I and b=2 must be retained.
	 */
	lis_vector_set_all((LIS_SCALAR)2.0,b);
	lis_vector_set_all((LIS_SCALAR)1.0,x);

	err = lis_solve_sequence(A,b,x,solver);
	if( err || solver->retcode!=LIS_SUCCESS )
	{
		fprintf(stderr,"AUTO_KEEP solve failed\n");
		failed = 1;
		goto cleanup;
	}

	if( solver->sequence_warm_start_kept!=LIS_TRUE )
	{
		fprintf(stderr,"AUTO_KEEP rejected exact x\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_small(
		solver->sequence_initial_residual,
		1.0e-14,
		"AUTO_KEEP initial residual");

	failed |= check_small(
		solver->sequence_true_relative_residual,
		1.0e-12,
		"AUTO_KEEP final residual");

	if( failed ) goto cleanup;

	printf(
		"AUTO_KEEP initial=%e rhs=%e final=%e\n",
		(double)solver->sequence_initial_residual,
		(double)solver->sequence_rhs_norm,
		(double)solver->sequence_true_relative_residual);

	/*
	 * PARTIAL KEEP:
	 * A=2I, b=[2,2], x=[0.75,0.75].
	 *
	 * ||b-Ax||2 = sqrt(0.5^2+0.5^2)
	 *            < ||b||2,
	 * so AUTO must retain this useful but non-exact guess.
	 */
	lis_vector_set_all((LIS_SCALAR)2.0,b);
	lis_vector_set_all((LIS_SCALAR)0.75,x);

	err = lis_solve_sequence(A,b,x,solver);
	if( err || solver->retcode!=LIS_SUCCESS )
	{
		fprintf(stderr,"AUTO_PARTIAL solve failed\n");
		failed = 1;
		goto cleanup;
	}

	if( solver->sequence_warm_start_kept!=LIS_TRUE )
	{
		fprintf(stderr,"AUTO_PARTIAL rejected useful x\n");
		failed = 1;
		goto cleanup;
	}

	if( solver->sequence_initial_residual<=0.0 ||
	    solver->sequence_initial_residual>=
	        solver->sequence_rhs_norm )
	{
		fprintf(
			stderr,
			"AUTO_PARTIAL invalid initial/rhs ratio %e/%e\n",
			(double)solver->sequence_initial_residual,
			(double)solver->sequence_rhs_norm);
		failed = 1;
		goto cleanup;
	}

	failed |= check_small(
		solver->sequence_true_relative_residual,
		1.0e-12,
		"AUTO_PARTIAL final residual");
	if( failed ) goto cleanup;

	if( solver->options[LIS_OPTIONS_INITGUESS_ZEROS]!=
	    saved_initguess_option )
	{
		fprintf(
			stderr,
			"AUTO_PARTIAL did not restore initguess option\n");
		failed = 1;
		goto cleanup;
	}

	printf(
		"AUTO_PARTIAL initial=%e rhs=%e ratio=%e final=%e\n",
		(double)solver->sequence_initial_residual,
		(double)solver->sequence_rhs_norm,
		(double)(solver->sequence_initial_residual /
		         solver->sequence_rhs_norm),
		(double)solver->sequence_true_relative_residual);

	/*
	 * REJECT:
	 * x=10 is worse than the zero guess.
	 */
	lis_vector_set_all((LIS_SCALAR)2.0,b);
	lis_vector_set_all((LIS_SCALAR)10.0,x);

	err = lis_solve_sequence(A,b,x,solver);
	if( err || solver->retcode!=LIS_SUCCESS )
	{
		fprintf(stderr,"AUTO_REJECT solve failed\n");
		failed = 1;
		goto cleanup;
	}

	if( solver->sequence_warm_start_kept!=LIS_FALSE )
	{
		fprintf(stderr,"AUTO_REJECT kept bad x\n");
		failed = 1;
		goto cleanup;
	}

	if( solver->sequence_initial_residual <=
	    solver->sequence_rhs_norm )
	{
		fprintf(stderr,"AUTO_REJECT test vector is not worse\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_small(
		solver->sequence_true_relative_residual,
		1.0e-12,
		"AUTO_REJECT final residual");
	if( failed ) goto cleanup;

	printf(
		"AUTO_REJECT initial=%e rhs=%e final=%e\n",
		(double)solver->sequence_initial_residual,
		(double)solver->sequence_rhs_norm,
		(double)solver->sequence_true_relative_residual);

	/*
	 * ZERO RHS + exact x=0:
	 * special case residual==0 must KEEP rather than fail the
	 * strict residual<rhs comparison.
	 */
	lis_vector_set_all((LIS_SCALAR)0.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(A,b,x,solver);
	if( err || solver->retcode!=LIS_SUCCESS )
	{
		fprintf(stderr,"ZERO_RHS_EXACT solve failed\n");
		failed = 1;
		goto cleanup;
	}

	if( solver->sequence_warm_start_kept!=LIS_TRUE )
	{
		fprintf(stderr,"ZERO_RHS_EXACT rejected exact x\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_small(
		solver->sequence_true_residual,
		1.0e-14,
		"ZERO_RHS_EXACT final residual");
	if( failed ) goto cleanup;

	/*
	 * Direct helper check.
	 */
	err = lis_solver_compute_true_residual(
		A,b,x,&rnorm,&bnorm);
	if( err )
	{
		fprintf(stderr,"TRUE_RESIDUAL helper failed\n");
		failed = 1;
		goto cleanup;
	}

	if( rnorm!=0.0 || bnorm!=0.0 )
	{
		fprintf(
			stderr,
			"TRUE_RESIDUAL expected 0/0 got %e/%e\n",
			(double)rnorm,
			(double)bnorm);
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_reset_sequence(solver);
	if( err )
	{
		fprintf(stderr,"RESET failed\n");
		failed = 1;
		goto cleanup;
	}

	if( solver->sequence_initial_residual!=0.0 ||
	    solver->sequence_rhs_norm!=0.0 ||
	    solver->sequence_true_residual!=0.0 ||
	    solver->sequence_true_relative_residual!=0.0 ||
	    solver->sequence_warm_start_kept!=LIS_FALSE )
	{
		fprintf(stderr,"RESET left S6 state\n");
		failed = 1;
		goto cleanup;
	}

	printf("LIS_AUTO_WARM_START_S6 PASSED\n");

cleanup:

	if( solver ) lis_solver_destroy(solver);
	if( x ) lis_vector_destroy(x);
	if( b ) lis_vector_destroy(b);
	if( A ) lis_matrix_destroy(A);

	lis_finalize();

	return failed ? 1 : 0;
}
