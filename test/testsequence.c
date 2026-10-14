#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#endif

#include "lis.h"
#include "lis_precon.h"
#include "lis_solver.h"

#include <stdio.h>

static int create_count = 0;
static int destroy_count = 0;

static LIS_INT
sequence_create(LIS_SOLVER solver, LIS_PRECON precon)
{
	(void)solver;
	(void)precon;

	create_count++;
	return LIS_SUCCESS;
}

static LIS_INT
sequence_psolve(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	(void)solver;
	return lis_vector_copy(b,x);
}

static LIS_INT
sequence_destroy(LIS_PRECON precon)
{
	(void)precon;

	destroy_count++;
	return LIS_SUCCESS;
}

static LIS_INT
create_matrix(LIS_SCALAR diagonal, LIS_MATRIX *A)
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
			diagonal,
			*A);
		if( err ) return err;
	}

	err = lis_matrix_set_type(*A,LIS_MATRIX_CSR);
	if( err ) return err;

	return lis_matrix_assemble(*A);
}

static int
check_counts(
	int creates,
	int destroys,
	const char *label)
{
	if( create_count!=creates ||
	    destroy_count!=destroys )
	{
		fprintf(
			stderr,
			"%s: create/destroy=%d/%d expected=%d/%d\n",
			label,
			create_count,
			destroy_count,
			creates,
			destroys);
		return 1;
	}

	return 0;
}

int main(int argc, char **argv)
{
	LIS_MATRIX A0 = NULL;
	LIS_MATRIX A1 = NULL;
	LIS_SOLVER solver = NULL;
	LIS_SOLVER legacy_solver = NULL;
	LIS_SOLVER sequence_solver = NULL;
	LIS_PRECON sequence_precon0 = NULL;
	void *sequence_work0 = NULL;
	LIS_VECTOR b = NULL;
	LIS_VECTOR x = NULL;
	LIS_PRECON p0 = NULL;
	LIS_PRECON p1 = NULL;
	LIS_INT err;
	int failed = 0;

	err = lis_initialize(&argc,&argv);
	if( err ) return 1;

	err = lis_precon_register_ex(
		"seqs4",
		sequence_create,
		sequence_psolve,
		sequence_psolve,
		sequence_destroy);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = create_matrix((LIS_SCALAR)2.0,&A0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = create_matrix((LIS_SCALAR)3.0,&A1);
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
		"-i gmres -p seqs4 "
		"-scale none -print none "
		"-restart 2 -maxiter 20 "
		"-maxiter_noimp 0 -tol 1e-12",
		solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	/* BUILD */
	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p0);
	if( err || p0==NULL )
	{
		fprintf(stderr,"BUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(1,0,"BUILD");
	if( failed ) goto cleanup;

	/* REUSE */
	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p1);
	if( err || p1!=p0 )
	{
		fprintf(stderr,"REUSE failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(1,0,"REUSE");
	if( failed ) goto cleanup;

	/*
	 * Configuration change must invalidate the persistent
	 * preconditioner even when matrix generations are unchanged.
	 */
	err = lis_solver_set_option(
		"-ilu_fill 1",
		solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"CONFIG_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(2,1,"CONFIG_REBUILD");
	if( failed ) goto cleanup;

	/* Numerical generation change -> REBUILD. */
	err = lis_matrix_mark_values_changed(A0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"VALUE_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(3,2,"VALUE_REBUILD");
	if( failed ) goto cleanup;

	/* Structural generation change -> REBUILD. */
	err = lis_matrix_mark_structure_changed(A0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_sequence_prepare_precon(
		A0,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"STRUCTURE_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(
		4,3,"STRUCTURE_REBUILD");
	if( failed ) goto cleanup;

	/*
	 * Different matrix object with compatible generations
	 * must also force a rebuild.
	 */
	err = lis_solver_sequence_prepare_precon(
		A1,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"MATRIX_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(5,4,"MATRIX_REBUILD");
	if( failed ) goto cleanup;

	/* Explicit reset. */
	err = lis_solver_sequence_precon_destroy(solver);
	if( err )
	{
		fprintf(stderr,"RESET failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(5,5,"RESET");
	if( failed ) goto cleanup;

	if( solver->sequence_precon!=NULL ||
	    solver->precon!=NULL )
	{
		fprintf(
			stderr,
			"RESET left live preconditioner state\n");
		failed = 1;
		goto cleanup;
	}

	/*
	 * Final BUILD is intentionally left alive.
	 * lis_solver_destroy() must release it.
	 */
	err = lis_solver_sequence_prepare_precon(
		A1,solver,&p1);
	if( err || p1==NULL )
	{
		fprintf(stderr,"FINAL_BUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(6,5,"FINAL_BUILD");
	if( failed ) goto cleanup;

	lis_solver_destroy(solver);
	solver = NULL;

	failed |= check_counts(
		6,6,"SOLVER_DESTROY");

	/*
	 * Legacy lis_solve() still owns a temporary preconditioner.
	 * It must destroy it and leave solver->precon clear.
	 */
	err = lis_solver_create(&legacy_solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_set_option(
		"-i gmres -p seqs4 "
		"-scale none -print none "
		"-restart 2 -maxiter 20 "
		"-maxiter_noimp 0 -tol 1e-12",
		legacy_solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(A0,&b);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(A0,&x);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_set_all((LIS_SCALAR)1.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve(A0,b,x,legacy_solver);
	if( err )
	{
		fprintf(stderr,"LEGACY_SOLVE failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(7,7,"LEGACY_SOLVE");
	if( failed ) goto cleanup;

	if( legacy_solver->precon!=NULL )
	{
		fprintf(stderr,
			"legacy lis_solve left dangling preconditioner\n");
		failed = 1;
		goto cleanup;
	}

	printf("LEGACY_SOLVE PASSED\n");

	/* ---------------------------------------------------- */
	/* S5 public sequence-solve API                         */
	/* ---------------------------------------------------- */

	err = lis_solver_create(&sequence_solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err = lis_solver_set_option(
		"-i gmres -p seqs4 "
		"-scale none -print none "
		"-restart 2 -maxiter 20 "
		"-maxiter_noimp 0 -tol 1e-12",
		sequence_solver);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	/* First solve -> BUILD. */
	lis_vector_set_all((LIS_SCALAR)1.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(
		A0,b,x,sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_0 failed\n");
		failed = 1;
		goto cleanup;
	}

	if( sequence_solver->retcode!=LIS_SUCCESS )
	{
		fprintf(
			stderr,
			"SEQUENCE_SOLVE_0 iterative status=%d\n",
			(int)sequence_solver->retcode);
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(8,7,"SEQUENCE_SOLVE_0");
	if( failed ) goto cleanup;

	if( sequence_solver->sequence_precon==NULL ||
	    sequence_solver->work==NULL )
	{
		fprintf(
			stderr,
			"SEQUENCE_SOLVE_0 did not retain reusable state\n");
		failed = 1;
		goto cleanup;
	}

	sequence_precon0 = sequence_solver->sequence_precon;
	sequence_work0 = (void *)sequence_solver->work;

	/* Second solve -> REUSE preconditioner and workspace. */
	lis_vector_set_all((LIS_SCALAR)2.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(
		A0,b,x,sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_1 failed\n");
		failed = 1;
		goto cleanup;
	}

	if( sequence_solver->retcode!=LIS_SUCCESS )
	{
		fprintf(
			stderr,
			"SEQUENCE_SOLVE_1 iterative status=%d\n",
			(int)sequence_solver->retcode);
		failed = 1;
		goto cleanup;
	}

	if( sequence_solver->sequence_precon!=sequence_precon0 )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_1 rebuilt preconditioner\n");
		failed = 1;
		goto cleanup;
	}

	if( (void *)sequence_solver->work!=sequence_work0 )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_1 rebuilt workspace\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(8,7,"SEQUENCE_SOLVE_1");
	if( failed ) goto cleanup;

	/* Matrix value generation change -> REBUILD. */
	err = lis_matrix_mark_values_changed(A0);
	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	lis_vector_set_all((LIS_SCALAR)3.0,b);
	lis_vector_set_all((LIS_SCALAR)0.0,x);

	err = lis_solve_sequence(
		A0,b,x,sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_REBUILD failed\n");
		failed = 1;
		goto cleanup;
	}

	if( sequence_solver->retcode!=LIS_SUCCESS )
	{
		fprintf(
			stderr,
			"SEQUENCE_SOLVE_REBUILD iterative status=%d\n",
			(int)sequence_solver->retcode);
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(
		9,8,"SEQUENCE_SOLVE_REBUILD");
	if( failed ) goto cleanup;

	if( sequence_solver->sequence_precon==NULL )
	{
		fprintf(stderr,"SEQUENCE_SOLVE_REBUILD lost preconditioner\n");
		failed = 1;
		goto cleanup;
	}

	/* Explicit sequence reset. */
	err = lis_solver_reset_sequence(sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_RESET failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(9,9,"SEQUENCE_RESET");
	if( failed ) goto cleanup;

	if( sequence_solver->sequence_precon!=NULL ||
	    sequence_solver->precon!=NULL ||
	    sequence_solver->work!=NULL ||
	    sequence_solver->rhistory!=NULL ||
	    sequence_solver->A!=NULL ||
	    sequence_solver->b!=NULL ||
	    sequence_solver->x!=NULL ||
	    sequence_solver->xx!=NULL )
	{
		fprintf(
			stderr,
			"SEQUENCE_RESET left persistent/transient state\n");
		failed = 1;
		goto cleanup;
	}

	/* Reset must be idempotent. */
	err = lis_solver_reset_sequence(sequence_solver);
	if( err )
	{
		fprintf(stderr,"SEQUENCE_RESET_2 failed\n");
		failed = 1;
		goto cleanup;
	}

	failed |= check_counts(9,9,"SEQUENCE_RESET_2");
	if( failed ) goto cleanup;

	printf("LIS_SEQUENCE_SOLVE_S5 PASSED\n");

cleanup:

	if( sequence_solver )
	{
		lis_solver_destroy(sequence_solver);
		sequence_solver = NULL;
	}

	if( x )
	{
		lis_vector_destroy(x);
		x = NULL;
	}

	if( b )
	{
		lis_vector_destroy(b);
		b = NULL;
	}

	if( legacy_solver )
	{
		lis_solver_destroy(legacy_solver);
		legacy_solver = NULL;
	}

	if( solver )
	{
		lis_solver_destroy(solver);
		solver = NULL;
	}

	if( A1 ) lis_matrix_destroy(A1);
	if( A0 ) lis_matrix_destroy(A0);

	printf(
		"LIS_SEQUENCE_PRECON_S4_S5 %s "
		"create=%d destroy=%d\n",
		failed ? "FAILED" : "PASSED",
		create_count,
		destroy_count);

	lis_finalize();

	return failed ? 1 : 0;
}
