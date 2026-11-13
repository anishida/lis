#ifdef HAVE_CONFIG_H
#include "lis_config.h"
#elif defined(HAVE_CONFIG_WIN_H)
#include "lis_config_win.h"
#endif

#include <math.h>
#include <stdio.h>

#include "lislib.h"


static LIS_INT
create_matrix(
	LIS_MATRIX *A)
{
	const LIS_INT n = 64;

	LIS_INT i;
	LIS_INT err;

	err =
		lis_matrix_create(
			LIS_COMM_WORLD,
			A);

	if( err ) return err;

	err =
		lis_matrix_set_size(
			*A,
			0,
			n);

	if( err ) return err;

	for(i=0;i<n;i++)
	{
		if( i>1 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i-2,
					(LIS_SCALAR)-0.05,
					*A);

			if( err ) return err;
		}

		if( i>0 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i-1,
					(LIS_SCALAR)-0.72,
					*A);

			if( err ) return err;
		}

		err =
			lis_matrix_set_value(
				LIS_INS_VALUE,
				i,
				i,
				(LIS_SCALAR)
				(4.75 +
				 0.011*(double)i),
				*A);

		if( err ) return err;

		if( i<n-1 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i+1,
					(LIS_SCALAR)-0.61,
					*A);

			if( err ) return err;
		}

		if( i<n-2 )
		{
			err =
				lis_matrix_set_value(
					LIS_INS_VALUE,
					i,
					i+2,
					(LIS_SCALAR)-0.04,
					*A);

			if( err ) return err;
		}
	}

	err =
		lis_matrix_set_type(
			*A,
			LIS_MATRIX_CSR);

	if( err ) return err;

	return lis_matrix_assemble(*A);
}


static double
vector_norm(
	LIS_VECTOR X)
{
	LIS_INT i;
	double n;

	n = 0.0;

	for(i=0;i<X->n;i++)
	{
		n =
			hypot(
				n,
				(double)X->value[i]);
	}

	return n;
}


static double
vector_checksum(
	LIS_VECTOR X)
{
	LIS_INT i;
	double s;

	s = 0.0;

	for(i=0;i<X->n;i++)
	{
		s +=
			(double)(i+1) *
			(double)X->value[i];
	}

	return s;
}


int
main(
	int argc,
	char **argv)
{
	LIS_MATRIX A;
	LIS_SOLVER solver;
	LIS_PRECON precon;

	LIS_VECTOR B;
	LIS_VECTOR X;
	LIS_VECTOR XH;

	LIS_INT i;
	LIS_INT err;

	int failed;

	A = NULL;
	solver = NULL;
	precon = NULL;

	B = NULL;
	X = NULL;
	XH = NULL;

	failed = 0;

	err =
		lis_initialize(
			&argc,
			&argv);

	if( err )
	{
		return 1;
	}

	err = create_matrix(&A);

	if( err )
	{
		fprintf(
			stderr,
			"matrix rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	err =
		lis_solver_create(
			&solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err =
		lis_solver_set_option(
			"-i gmres "
			"-p ilu "
			"-ilu_fill 0 "
			"-scale none "
			"-print none "
			"-precon_half 2",
			solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err =
		lis_solver_set_matrix(
			A,
			solver);

	if( err )
	{
		failed = 1;
		goto cleanup;
	}

	err =
		lis_precon_create(
			solver,
			&precon);

	if( err )
	{
		fprintf(
			stderr,
			"precon create rc=%d\n",
			(int)err);

		precon = NULL;
		failed = 1;
		goto cleanup;
	}

	solver->precon = precon;

	if( precon->half==NULL )
	{
		fprintf(
			stderr,
			"half ILUK not accepted\n");

		failed = 1;
		goto cleanup;
	}

	err = lis_vector_duplicate(A,&B);
	if( err ) goto vector_fail;

	err = lis_vector_duplicate(A,&X);
	if( err ) goto vector_fail;

	err = lis_vector_duplicate(A,&XH);
	if( err ) goto vector_fail;

	for(i=0;i<A->n;i++)
	{
		B->value[i] =
			(LIS_SCALAR)
			(0.25 +
			 0.019*(double)i);
	}

	err =
		lis_psolve_iluk_csr(
			solver,
			B,
			X);

	if( err )
	{
		fprintf(
			stderr,
			"psolve rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	err =
		lis_psolveh_iluk_csr(
			solver,
			B,
			XH);

	if( err )
	{
		fprintf(
			stderr,
			"psolveh rc=%d\n",
			(int)err);

		failed = 1;
		goto cleanup;
	}

	printf(
		"HALF_ILUK_SIMD "
		"norm=%.17g checksum=%.17g "
		"normh=%.17g checksumh=%.17g "
		"half=1\n",
		vector_norm(X),
		vector_checksum(X),
		vector_norm(XH),
		vector_checksum(XH));

	goto cleanup;


vector_fail:

	fprintf(
		stderr,
		"vector allocation rc=%d\n",
		(int)err);

	failed = 1;


cleanup:

	if( XH!=NULL )
		lis_vector_destroy(XH);

	if( X!=NULL )
		lis_vector_destroy(X);

	if( B!=NULL )
		lis_vector_destroy(B);

	if( precon!=NULL )
	{
		solver->precon = NULL;
		lis_precon_destroy(precon);
	}

	if( solver!=NULL )
		lis_solver_destroy(solver);

	if( A!=NULL )
		lis_matrix_destroy(A);

	lis_finalize();

	return failed ? 1 : 0;
}
