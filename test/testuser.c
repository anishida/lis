#include "lis.h"
#include <stdio.h>

typedef struct
{
	LIS_INT n;
	LIS_SCALAR diag;
	LIS_SCALAR offdiag;
} USER_TRIDIAG;

static LIS_INT shell_matvec(void *ctx, const LIS_SCALAR *x, LIS_SCALAR *y)
{
	USER_TRIDIAG *A = (USER_TRIDIAG *)ctx;
	LIS_INT i;

	if( A->n==1 )
	{
		y[0] = A->diag*x[0];
		return LIS_SUCCESS;
	}

	y[0] = A->diag*x[0] + A->offdiag*x[1];
	for(i=1;i<A->n-1;i++)
	{
		y[i] = A->offdiag*x[i-1] + A->diag*x[i] + A->offdiag*x[i+1];
	}
	y[A->n-1] = A->offdiag*x[A->n-2] + A->diag*x[A->n-1];
	return LIS_SUCCESS;
}

/* Real symmetric test operator: A^H == A. */
static LIS_INT shell_matvech(void *ctx, const LIS_SCALAR *x, LIS_SCALAR *y)
{
	return shell_matvec(ctx,x,y);
}

static LIS_INT userdiag_create(LIS_SOLVER solver, LIS_PRECON precon)
{
	(void)solver;
	(void)precon;
	return LIS_SUCCESS;
}

static LIS_INT userdiag_psolve(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	USER_TRIDIAG *A;
	void *ctx = NULL;
	LIS_INT i, err;

	err = lis_matrix_get_user_data(solver->A,&ctx);
	if( err ) return err;
	A = (USER_TRIDIAG *)ctx;

	for(i=0;i<b->n;i++) x->value[i] = b->value[i] / A->diag;
	return LIS_SUCCESS;
}

static LIS_INT userdiag_psolveh(LIS_SOLVER solver, LIS_VECTOR b, LIS_VECTOR x)
{
	return userdiag_psolve(solver,b,x);
}

int main(int argc, char *argv[])
{
	LIS_MATRIX A;
	LIS_VECTOR b, x, exact;
	LIS_SOLVER solver;
	USER_TRIDIAG shell;
	LIS_INT i, err;
	LIS_REAL nrm2;
	const LIS_INT n = 100;

	lis_initialize(&argc,&argv);

	shell.n = n;
	shell.diag = 2.0;
	shell.offdiag = -1.0;

	err = lis_precon_register("userdiag",
	                          userdiag_create,
	                          userdiag_psolve,
	                          userdiag_psolveh);
	if( err )
	{
		fprintf(stderr,"lis_precon_register failed: %d\n",(int)err);
		lis_finalize();
		return 1;
	}

	lis_matrix_create(LIS_COMM_WORLD,&A);
	lis_matrix_set_size(A,0,n);
	err = lis_matrix_set_user(A,&shell,shell_matvec,shell_matvech);
	if( err )
	{
		fprintf(stderr,"lis_matrix_set_user failed: %d\n",(int)err);
		lis_matrix_destroy(A);
		lis_finalize();
		return 1;
	}

	lis_vector_duplicate(A,&b);
	lis_vector_duplicate(A,&x);
	lis_vector_duplicate(A,&exact);
	for(i=0;i<n;i++) exact->value[i] = 1.0;

	err = lis_matvec(A,exact,b);
	if( err )
	{
		fprintf(stderr,"lis_matvec failed: %d\n",(int)err);
		return 1;
	}
	lis_vector_set_all(0.0,x);

	lis_solver_create(&solver);
	err = lis_solver_set_option(
		"-i bicg -p userdiag -tol 1.0e-12 -maxiter 1000 -print none",
		solver);
	if( err )
	{
		fprintf(stderr,"lis_solver_set_option failed: %d\n",(int)err);
		return 1;
	}

	err = lis_solve(A,b,x,solver);
	if( err )
	{
		fprintf(stderr,"lis_solve failed: %d\n",(int)err);
		return 1;
	}

	lis_vector_axpy(-1.0,exact,x);
	lis_vector_nrm2(x,&nrm2);
	printf("LIS_MATRIX_USER solution error = %.15e\n",(double)nrm2);

	lis_solver_destroy(solver);
	lis_vector_destroy(exact);
	lis_vector_destroy(x);
	lis_vector_destroy(b);
	lis_matrix_destroy(A);
	lis_precon_register_free();
	lis_finalize();

	if( nrm2>1.0e-8 )
	{
		fprintf(stderr,"FAILED\n");
		return 1;
	}
	printf("PASSED\n");
	return 0;
}
